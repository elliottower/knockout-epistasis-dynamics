"""Write the identities the ODE rerun registrations pin: the trajectory-code digest, the SHA-256
of each frozen configuration file, and every paper network's model digest.

    uv run python -m scripts.record_registered_identities

Run it at the freeze, and commit its output with the registrations. Every sweep record carries
the same digests in its configuration, so a record can be checked against this file.

Writes experiments/ode_v2/registered_identities.json.
"""

import hashlib
import json

from scripts.ode_coalition_sweep import ALL_MODELS
from scripts.ode_engine import TRAJECTORY_CODE, code_digest, compile_model, model_digest
from scripts.paths import PROJECT_ROOT, RESULTS

OUT = PROJECT_ROOT / "experiments" / "ode_v2" / "registered_identities.json"
CONFIGS = ("experiments/ode_v2/solver_primary.json", "experiments/ode_v2/classifier_primary.json")


def build() -> dict:
    """The identities, computed from the files as they are now."""
    table = json.loads((RESULTS / "paper_number_reconciliation.json").read_text())["per_network_table"]
    record = {
        "generated_by": "uv run python -m scripts.record_registered_identities",
        "code_sha256": code_digest(),
        "code_files": [p.relative_to(PROJECT_ROOT).as_posix() for p in TRAJECTORY_CODE],
        "config_sha256": {f: hashlib.sha256((PROJECT_ROOT / f).read_bytes()).hexdigest() for f in CONFIGS},
        "networks": {},
    }
    for name in sorted(table, key=lambda m: (table[m]["n"], m)):
        info = ALL_MODELS[name]
        net = compile_model(info["rules"], info["output_nodes"])
        record["networks"][name] = {"n_nodes": len(net.node_names), "output_nodes": list(info["output_nodes"]),
                                    "model_sha256": model_digest(net),
                                    "boolean_delta_3plus_pp_published": table[name]["delta_3plus_pp"]}
    return record


def main() -> None:
    record = build()
    OUT.write_text(json.dumps(record, indent=2) + "\n")
    print(f"wrote {OUT.relative_to(PROJECT_ROOT)}: code {record['code_sha256'][:12]}, {len(record['networks'])} networks")


if __name__ == "__main__":
    main()
