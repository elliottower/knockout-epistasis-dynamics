"""Run chunks of the Boolean cycle-handling audit on Modal, one container per chunk.

    uv run --with modal==1.4.3 modal run -m scripts.modal_boolean_audit \
        --networks calzone_cell_fate,grieco_bladder

(Run through the project environment, which the local comparison needs for NumPy.)

Each container computes one chunk with `audit_boolean_cycle_handling.compute_chunk`, the
function the local audit uses, and writes it to the volume `knockout-boolean-audit`; a chunk
already on the volume is not recomputed. The local entrypoint then downloads every chunk. Where
the local audit already wrote the same chunk, the two must agree exactly, and the comparison is
recorded; where it has not, the Modal chunk is written into the local shard directory, where
the local audit picks it up. Boolean dynamics on int8 states are platform-independent, so any
disagreement is a defect, not rounding.

Every file the chunk imports or reads is copied into the image and hashed at launch; a container
whose copy differs refuses to run. Writes results/audit/boolean_cycle_handling/modal_crosscheck.json.
After this, `uv run python -m scripts.audit_boolean_cycle_handling --networks ...` assembles the
records from the shards without recomputing any.
"""

import io
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import modal
import numpy as np

from scripts.audit_boolean_cycle_handling import CHUNK, compute_chunk
from scripts.modal_provenance import BOOLEAN_AUDIT_FILES as FILES
from scripts.modal_provenance import file_hashes, verify
from scripts.ode_coalition_sweep import ALL_MODELS
from scripts.paths import AUDIT, PROJECT_ROOT, RESULTS

REMOTE_ROOT = Path("/root/repo")
VOLUME_ROOT = Path("/audit")
LOCAL_SHARDS = AUDIT / "boolean_cycle_handling" / "shards"
CROSSCHECK = AUDIT / "boolean_cycle_handling" / "modal_crosscheck.json"

image = modal.Image.debian_slim(python_version="3.13").pip_install(
    "numpy==2.2.6", "scipy==1.15.3", "tqdm==4.67.1").env({"PYTHONPATH": str(REMOTE_ROOT)})
for relative in FILES:
    image = image.add_local_file(PROJECT_ROOT / relative, str(REMOTE_ROOT / relative), copy=True)

# include_source=False: Modal would otherwise mount the local `scripts` package over the copies
# above, and the container would import files the hash check never saw.
app = modal.App("knockout-boolean-audit", image=image, include_source=False)
volume = modal.Volume.from_name("knockout-boolean-audit", create_if_missing=True)


def now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


@app.function(cpu=1.0, memory=4096, timeout=24 * 60 * 60, volumes={str(VOLUME_ROOT): volume},
              retries=modal.Retries(max_retries=3, initial_delay=10.0))
def chunk(name: str, start: int, end: int, hashes: dict[str, str]) -> dict:
    verify(REMOTE_ROOT, FILES, hashes)
    path = VOLUME_ROOT / name / f"{start:07d}.npz"
    volume.reload()
    if not path.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
        impl, reg, totals = compute_chunk(name, start, end)
        tmp = path.with_name(path.stem + ".partial.npz")
        np.savez_compressed(tmp, impl=impl, reg=reg, start=start, end=end, counters=json.dumps(totals))
        tmp.replace(path)
        volume.commit()
    print(f"[{now()}] {name} [{start}, {end}) on the volume", flush=True)
    return {"name": name, "start": start, "task_id": os.environ.get("MODAL_TASK_ID")}


def chunk_ranges(names: list[str]) -> list[tuple[str, int, int]]:
    jobs = []
    for name in names:
        n = len(ALL_MODELS[name]["rules"])
        jobs.extend((name, s, min(s + CHUNK, 2**n)) for s in range(0, 2**n, CHUNK))
    return jobs


@app.local_entrypoint()
def main(networks: str):
    names = networks.split(",")
    table = json.loads((RESULTS / "paper_number_reconciliation.json").read_text())["per_network_table"]
    unknown = [m for m in names if m not in table]
    if unknown:
        raise SystemExit(f"not paper networks: {unknown}")
    hashes = file_hashes(PROJECT_ROOT, FILES)
    jobs = chunk_ranges(names)
    print(f"[{now()}] {len(jobs)} chunks across {names}", flush=True)
    tasks = [r for r in chunk.starmap([(n, s, e, hashes) for n, s, e in jobs])]

    report = json.loads(CROSSCHECK.read_text()) if CROSSCHECK.exists() else {}
    for name in names:
        compared, written, disagreeing = 0, 0, []
        for n, start, end in (j for j in jobs if j[0] == name):
            data = b"".join(volume.read_file(f"{name}/{start:07d}.npz"))
            local = LOCAL_SHARDS / name / f"{start:07d}.npz"
            if local.exists():
                with np.load(io.BytesIO(data)) as remote, np.load(local) as mine:
                    same = (all(np.array_equal(remote[k], mine[k]) for k in ("impl", "reg", "start", "end"))
                            and json.loads(str(remote["counters"])) == json.loads(str(mine["counters"])))
                compared += 1
                if not same:
                    disagreeing.append(start)
            else:
                local.parent.mkdir(parents=True, exist_ok=True)
                tmp = local.with_name(local.stem + ".partial.npz")
                tmp.write_bytes(data)
                tmp.replace(local)
                written += 1
        report[name] = {"checked_at": now(), "n_chunks": sum(1 for j in jobs if j[0] == name),
                        "compared_with_local": compared, "disagreeing_chunks": disagreeing,
                        "written_from_modal": written,
                        "task_ids": sorted({t["task_id"] for t in tasks if t["name"] == name}),
                        "file_sha256": hashes}
        print(f"[{now()}] {name}: compared {compared}, disagreeing {len(disagreeing)}, written {written}", flush=True)
    CROSSCHECK.write_text(json.dumps(report, indent=2))
    print(f"wrote {CROSSCHECK}")
