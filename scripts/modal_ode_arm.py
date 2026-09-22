"""Run ODE arms on Modal: shards, then replicates, then finalize, for each network.

    modal run --detach -m scripts.modal_ode_arm::main --run-name primary \
        --networks lambda_phage,davidich_yeast --construction hillcube_normalized --hill-n 10 \
        --solver experiments/ode_v2/solver_primary.json \
        --classifier experiments/ode_v2/classifier_primary.json --n-init 32 --seed 42 --keep-states

    modal run -m scripts.modal_ode_arm::smoke      # the pipeline end to end on lambda_phage

A registered launch (`main`) first runs `scripts.verify_freeze.require` and refuses unless the
checkout is the tagged freeze commit with every registered file unchanged; the tag and commits it
returns go into the launch record and every record's launch manifest. The smoke test runs
before the freeze and is not gated.

The pipeline is idempotent. A shard resumes from its checkpoint and a finished shard returns at
once; a replicate chunk is reused only if its file holds the same coalitions under the same
fingerprint. If the orchestrator reaches its 24-hour limit, or anything is interrupted, running
the same command again continues where the volume left off.

Every file the run imports or reads is copied into the image. Its SHA-256 is computed here, at
launch, and again inside every container; a container whose copy differs refuses to run, so a
record cannot come from code other than the code on disk when the run was launched.

Each container that works on a shard or a replicate chunk leaves a sidecar file naming its
Modal task id, so the record can show that replicates ran in containers other than the shards'.

State lives on the volume `knockout-ode-v2`, under /results/<run name>/. To check a run, look
at the app state first (`modal app list`, `modal app logs <app id>`), then the files
(`modal run -m scripts.modal_ode_arm::status --run-name primary`). A stopped app leaves
behind the files it had already written, which look like slow progress.
"""

import dataclasses
import importlib
import json
import os
import shutil
from datetime import datetime, timezone
from pathlib import Path

import modal

# First among the project's imports, so that in the shadow image the stand-in data_utils is loaded
# before scripts/ode_coalition_sweep.py puts the repository root first on the import path.
import data_utils  # noqa: F401
from scripts.modal_provenance import ODE_ARM_FILES as FILES
from scripts.modal_provenance import ProvenanceError, file_hashes, verify
from scripts.paths import PROJECT_ROOT, RESULTS
from scripts.verify_freeze import require as require_freeze

REMOTE_ROOT = Path("/root/repo")
VOLUME_ROOT = Path("/results")
COALITIONS_PER_SHARD = 4096
CHECKPOINT_EVERY = 64
REPLICATE_COALITIONS_PER_CHUNK = 256
DAY = 24 * 60 * 60

image = modal.Image.debian_slim(python_version="3.13").pip_install(
    "numpy==2.2.6", "scipy==1.15.3", "tqdm==4.67.1").env({"PYTHONPATH": str(REMOTE_ROOT)})
for relative in FILES:
    image = image.add_local_file(PROJECT_ROOT / relative, str(REMOTE_ROOT / relative), copy=True)

# include_source=False: Modal would otherwise mount the local `scripts` package over the copies
# above, and the container would import files the hash check never saw.
app = modal.App("knockout-ode-v2", image=image, include_source=False)
# The same image with a stand-in data_utils placed ahead of the hashed copies on the import path.
# The smoke test runs one function in it, which must refuse to proceed.
SHADOW = "tests/fixtures/shadow/data_utils.py"
SHADOW_REMOTE = "/root/shadow/data_utils.py"
shadow_image = image.add_local_file(PROJECT_ROOT / SHADOW, SHADOW_REMOTE, copy=True).env(
    {"PYTHONPATH": f"/root/shadow:{REMOTE_ROOT}", "SHADOWED_ROOT": str(REMOTE_ROOT)})
volume = modal.Volume.from_name("knockout-ode-v2", create_if_missing=True)

with image.imports():
    import numpy as np

    from scripts.run_ode_arm import (
        REPLICATE_FRACTION,
        ArmSpec,
        ShardError,
        finalize,
        plan_replicates,
        plan_shards,
        replicate_file,
        replicate_is_current,
        run_replicate,
        run_shard,
        shard_file,
    )


class SimulatedPreemption(Exception):
    """Raised once, after a checkpoint, by a shard asked to test its own resumption."""


def now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def log(message: str) -> None:
    print(f"[{now()}] {message}", flush=True)


def remote_spec(fields: dict, expected_hashes: dict[str, str]) -> "ArmSpec":
    verify(REMOTE_ROOT, FILES, expected_hashes)
    classifier = fields["classifier_path"]
    return ArmSpec(**{**fields, "solver_path": REMOTE_ROOT / fields["solver_path"],
                      "classifier_path": REMOTE_ROOT / classifier if classifier else None})


def leave_sidecar(path: Path, **fields) -> None:
    """Record which container is working on a shard or chunk: one file per container, never overwritten."""
    task = os.environ.get("MODAL_TASK_ID", "unknown")
    sidecar = path.with_name(f"{path.stem}.task_{task}.json")
    previous = json.loads(sidecar.read_text()) if sidecar.exists() else {}
    sidecar.write_text(json.dumps({**previous, "task_id": task, **fields}, indent=2))
    volume.commit()


@app.function(cpu=1.0, memory=1024, timeout=DAY, volumes={str(VOLUME_ROOT): volume},
              retries=modal.Retries(max_retries=3, initial_delay=10.0))
def shard(fields: dict, hashes: dict[str, str], run_name: str, start: int, end: int,
          checkpoint_every: int = CHECKPOINT_EVERY, preempt_after_checkpoints: int = 0) -> list[int]:
    spec = remote_spec(fields, hashes)
    out = VOLUME_ROOT / run_name
    out.mkdir(parents=True, exist_ok=True)
    volume.reload()  # a retried shard resumes from its last committed checkpoint
    path = shard_file(out, spec, start, end)
    marker = path.with_name(path.stem + ".preempted")
    leave_sidecar(path, started_at=now())
    log(f"{spec.stem} shard [{start}, {end}) starting")
    checkpoints = [0]

    def on_checkpoint() -> None:
        volume.commit()
        checkpoints[0] += 1
        if preempt_after_checkpoints and checkpoints[0] == preempt_after_checkpoints and not marker.exists():
            marker.write_text(os.environ.get("MODAL_TASK_ID", "unknown"))
            volume.commit()
            raise SimulatedPreemption(f"shard [{start}, {end}) stopped after checkpoint {checkpoints[0]}")

    run_shard(spec, start, end, out, checkpoint_every=checkpoint_every, on_checkpoint=on_checkpoint)
    volume.commit()
    leave_sidecar(path, finished_at=now())
    log(f"{spec.stem} shard [{start}, {end}) complete")
    return [start, end]


@app.function(cpu=1.0, memory=1024, timeout=DAY, volumes={str(VOLUME_ROOT): volume},
              retries=modal.Retries(max_retries=3, initial_delay=10.0))
def replicate(fields: dict, hashes: dict[str, str], run_name: str, chunk: int, coalitions: list[int]) -> int:
    spec = remote_spec(fields, hashes)
    path = replicate_file(VOLUME_ROOT / run_name, spec, chunk)
    volume.reload()
    coalitions = np.array(coalitions, dtype=np.int64)
    if replicate_is_current(spec, chunk, coalitions, path):
        log(f"{spec.stem} replicate chunk {chunk} already on the volume under this run's fingerprint")
        return chunk
    leave_sidecar(path, started_at=now())
    log(f"{spec.stem} replicate chunk {chunk}: {len(coalitions)} coalitions starting")
    run_replicate(spec, chunk, coalitions, path)
    volume.commit()
    leave_sidecar(path, finished_at=now())
    log(f"{spec.stem} replicate chunk {chunk} complete")
    return chunk


@app.function(cpu=1.0, memory=8192, timeout=DAY, volumes={str(VOLUME_ROOT): volume})
def plan(fields: dict, hashes: dict[str, str], run_name: str, n_chunks: int) -> list[list[int]]:
    volume.reload()
    return [c.tolist() for c in plan_replicates(remote_spec(fields, hashes), VOLUME_ROOT / run_name, n_chunks)]


def container_record(out: Path, stem: str) -> dict:
    """Which containers worked on the shards and on the replicate chunks, from the sidecars."""
    def tasks(kind: str) -> dict[str, list[str]]:
        found: dict[str, list[str]] = {}
        for f in sorted(out.glob(f"{stem}.{kind}_*.task_*.json")):
            found.setdefault(f.name.split(".task_")[0], []).append(json.loads(f.read_text())["task_id"])
        return found
    shards, replicates = tasks("shard"), tasks("replicate")
    shard_tasks = {t for ts in shards.values() for t in ts}
    replicate_tasks = {t for ts in replicates.values() for t in ts}
    return {"shards": shards, "replicates": replicates,
            "replicates_ran_in_containers_distinct_from_shards": not (shard_tasks & replicate_tasks),
            "shards_resumed_in_a_second_container": sorted(k for k, v in shards.items() if len(set(v)) > 1)}


@app.function(cpu=1.0, memory=16384, timeout=DAY, volumes={str(VOLUME_ROOT): volume})
def finish(fields: dict, hashes: dict[str, str], run_name: str, freeze: dict | None = None) -> dict:
    volume.reload()
    spec = remote_spec(fields, hashes)
    out = VOLUME_ROOT / run_name
    record = finalize(spec, out, {"run_name": run_name, "spec": fields, "file_sha256": hashes, "freeze": freeze})
    containers = container_record(out, spec.stem)
    (out / f"{spec.stem}.containers.json").write_text(json.dumps(containers, indent=2))
    volume.commit()
    log(f"{spec.stem} finalized: scored {record['scored']}")
    return {"stem": spec.stem, "scored": record["scored"], "summary": record["summary"],
            "replicate": record["replicate"], "configuration_sha256": record["configuration_sha256"],
            "containers": containers}


@app.function(cpu=1.0, memory=2048, timeout=DAY, volumes={str(VOLUME_ROOT): volume})
def orchestrate(specs: list[dict], hashes: dict[str, str], run_name: str,
                coalitions_per_shard: int = COALITIONS_PER_SHARD, checkpoint_every: int = CHECKPOINT_EVERY,
                preempt_after_checkpoints: int = 0, freeze: dict | None = None) -> list[dict]:
    """Shards for every network at once, then each network's replicates, then its record.
    With preempt_after_checkpoints set, the first shard of each network stops once after that
    many checkpoints, and its retry must resume it."""
    shard_calls = {}
    for fields in specs:
        spec = remote_spec(fields, hashes)
        ranges = plan_shards(len(spec.net().node_names), coalitions_per_shard)
        shard_calls[spec.stem] = [shard.spawn(fields, hashes, run_name, s, e, checkpoint_every,
                                              preempt_after_checkpoints if i == 0 else 0)
                                  for i, (s, e) in enumerate(ranges)]
        log(f"{spec.stem}: {len(ranges)} shards spawned")
    records = []
    for fields in specs:
        spec = remote_spec(fields, hashes)
        for call in shard_calls[spec.stem]:
            call.get()
        n_nodes = len(spec.net().node_names)
        n_chunks = max(1, round(REPLICATE_FRACTION * 2**n_nodes / REPLICATE_COALITIONS_PER_CHUNK))
        chunks = plan.remote(fields, hashes, run_name, n_chunks)
        list(replicate.starmap([(fields, hashes, run_name, i, c) for i, c in enumerate(chunks)]))
        records.append(finish.remote(fields, hashes, run_name, freeze))
        log(f"{spec.stem}: done")
    return records


@app.function(cpu=1.0, memory=4096, timeout=DAY, volumes={str(VOLUME_ROOT): volume})
def refuses_stale_files(fields: dict, hashes: dict[str, str], run_name: str) -> dict:
    """Copy a finished run, plant a shard and then a replicate produced under another seed in
    each copy, and check that finalization refuses both."""
    volume.reload()
    spec = remote_spec(fields, hashes)
    other = dataclasses.replace(spec, seed=spec.seed + 1)
    results = {}
    for kind in ("shard", "replicate"):
        copy = VOLUME_ROOT / f"{run_name}-stale-{kind}"
        if copy.exists():
            shutil.rmtree(copy)
        shutil.copytree(VOLUME_ROOT / run_name, copy)
        if kind == "shard":
            first = sorted(copy.glob(f"{spec.stem}.shard_*-*.npz"))[0]
            with np.load(first) as z:
                start, end = json.loads(str(z["fingerprint"]))["coalition_range"]
            first.unlink()
            run_shard(other, start, end, copy)
        else:
            path = replicate_file(copy, spec, 0)
            run_replicate(other, 0, np.load(path)["coalitions"], path)
        try:
            finalize(spec, copy)
            results[kind] = {"refused": False}
        except ShardError as err:
            results[kind] = {"refused": True, "message": str(err)}
        shutil.rmtree(copy)
    volume.commit()
    return results


@app.function(image=shadow_image, cpu=1.0, memory=2048, timeout=60 * 60)
def shadowed_module_is_refused(fields: dict, hashes: dict[str, str]) -> dict:
    """Negative control: in the shadow image every file still hashes correctly, and the check must
    refuse because data_utils is imported from the stand-in. Reports where data_utils was loaded
    from, so that a pass cannot come from a stand-in that was never imported."""
    loaded_from = importlib.import_module("data_utils").__file__
    try:
        remote_spec(fields, hashes)
    except ProvenanceError as err:
        return {"refused": True, "message": str(err), "data_utils_loaded_from": loaded_from}
    return {"refused": False, "data_utils_loaded_from": loaded_from}


@app.function(cpu=1.0, memory=2048, timeout=60 * 60, volumes={str(VOLUME_ROOT): volume})
def inspect_volume(run_name: str) -> dict:
    volume.reload()
    out = VOLUME_ROOT / run_name
    report = {}
    for f in sorted(out.glob("*.npz")):
        if ".shard_" in f.name and not f.name.endswith(".partial.npz"):
            with np.load(f) as z:
                start, end = json.loads(str(z["fingerprint"]))["coalition_range"]
                report[f.name] = {"completed": int(z["completed"]), "of": end - start}
        elif ".replicate_" in f.name and not f.name.endswith(".partial.npz"):
            report[f.name] = "complete"
    for f in sorted(out.glob("*.json")):
        if ".task_" not in f.name and not f.name.endswith(".partial.json"):
            report[f.name] = "written"
    return report


def launch_record(run_name: str, specs: list[dict], hashes: dict[str, str], **extra) -> dict:
    launch = {"run_name": run_name, "launched_at": now(), "specs": specs, "file_sha256": hashes, **extra}
    launches = RESULTS / "ode_v2" / "launches"
    launches.mkdir(parents=True, exist_ok=True)
    (launches / f"{run_name}__{launch['launched_at'].replace(':', '')}.json").write_text(json.dumps(launch, indent=2))
    return launch


@app.local_entrypoint()
def main(run_name: str, networks: str, construction: str, hill_n: float, solver: str,
         n_init: int, seed: int, classifier: str = "", clamp_value: float = 0.0, keep_states: bool = False,
         hill_k: float = 0.5):
    freeze = require_freeze()
    hashes = file_hashes(PROJECT_ROOT, FILES)
    for path in (solver, classifier):
        if path and path not in FILES:
            raise SystemExit(f"{path} is not among the files copied into the image")
    specs = [{"network": n, "construction": construction, "hill_n": hill_n, "hill_k": hill_k, "solver_path": solver,
              "classifier_path": classifier or None, "n_init": n_init, "seed": seed,
              "clamp_value": clamp_value, "keep_states": keep_states} for n in networks.split(",")]
    launch_record(run_name, specs, hashes, freeze=freeze)
    log(f"{run_name}: launching {len(specs)} networks from {freeze['tag']} at {freeze['commit'][:12]}")
    for record in orchestrate.remote(specs, hashes, run_name, freeze=freeze):
        rep = record["replicate"]
        log(f"{record['stem']}: scored {record['scored']}, trajectories {record['summary']['trajectories']}, "
            f"replicate mismatches {rep['status_mismatches']}/{rep['class_mismatches']}/{rep['value_mismatches']}")


@app.local_entrypoint()
def smoke():
    """lambda_phage, 2 initial states, 4 shards of 32 coalitions checkpointed every 8, the first
    shard preempted once after its first checkpoint; then the stale-file checks. Writes
    results/ode_v2/smoke/<run name>.json locally."""
    hashes = file_hashes(PROJECT_ROOT, FILES)
    run_name = f"smoke-{now().replace(':', '')}"
    fields = {"network": "lambda_phage", "construction": "hillcube_normalized", "hill_n": 10.0, "hill_k": 0.5,
              "solver_path": "experiments/ode_v2/solver_primary.json",
              "classifier_path": "experiments/ode_v2/classifier_primary.json",
              "n_init": 2, "seed": 42, "clamp_value": 0.0, "keep_states": True}
    launch_record(run_name, [fields], hashes, purpose="smoke test")
    record = orchestrate.remote([fields], hashes, run_name, 32, 8, 1)[0]
    stale = refuses_stale_files.remote(fields, hashes, run_name)
    shadow = shadowed_module_is_refused.remote(fields, hashes)
    containers = record["containers"]
    checks = {
        "a_shard_was_preempted_and_resumed_in_a_second_container": bool(containers["shards_resumed_in_a_second_container"]),
        "replicates_ran_in_containers_distinct_from_shards": containers["replicates_ran_in_containers_distinct_from_shards"],
        "replicates_agree": (record["replicate"]["status_mismatches"], record["replicate"]["class_mismatches"],
                             record["replicate"]["value_mismatches"]) == (0, 0, 0),
        "record_scored": record["scored"],
        "stale_shard_refused": stale["shard"]["refused"],
        "stale_replicate_refused": stale["replicate"]["refused"],
        "shadowed_module_refused": (shadow["refused"] and "data_utils" in shadow.get("message", "")
                                    and shadow["data_utils_loaded_from"] == SHADOW_REMOTE),
    }
    report = {"run_name": run_name, "checks": checks, "passed": all(checks.values()),
              "record": record, "stale_files": stale, "shadowed_module": shadow}
    out = RESULTS / "ode_v2" / "smoke"
    out.mkdir(parents=True, exist_ok=True)
    (out / f"{run_name}.json").write_text(json.dumps(report, indent=2))
    for name, ok in checks.items():
        log(f"{'PASS' if ok else 'FAIL'}  {name}")
    log(f"wrote {out / (run_name + '.json')}")


@app.local_entrypoint()
def status(run_name: str):
    for name, state in inspect_volume.remote(run_name).items():
        print(f"{name:<90} {state}")
