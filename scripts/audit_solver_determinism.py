"""Does the legacy ODE engine reproduce itself on identical inputs?

Runs scripts/ode_coalition_sweep.simulate_ode_output, unchanged and with its wall-clock
timeout disabled, several times on the same coalition and the same initial states, and
reports how many distinct results come back. Any spread here is upstream of the new
engine: it comes from NumPy/SciPy, not from either construction.

    uv run python -m scripts.audit_solver_determinism

Writes results/audit/solver_determinism.json. The engine tests compare at a relative
tolerance justified by the spread this measures.
"""

import json
import os

import numpy as np

from scripts.ode_coalition_sweep import ALL_MODELS, simulate_ode_output
from scripts.ode_engine import clamp_mask_for, initial_states
from scripts.paths import AUDIT, PROJECT_ROOT

OUT = AUDIT / "solver_determinism.json"
NETWORKS = ("arellano_rootstem", "davidich_yeast", "lambda_phage")
REPEATS = 6
N_INIT = 4
SEED = 42


def probe(name: str) -> list[dict]:
    info = ALL_MODELS[name]
    rules = info["rules"]
    names = list(rules)
    init = initial_states(N_INIT, len(rules), SEED)
    output_indices = [names.index(o) for o in info["output_nodes"]]
    rows = []
    for coalition in (2 ** len(rules) - 1, 2 ** (len(rules) - 1) + 3):
        mask = clamp_mask_for(coalition, len(rules))
        reps = np.stack([simulate_ode_output(rules, names, mask, 0, output_indices, init,
                                             t_max=30.0, t_tail=10.0, hill_n=10.0, hill_k=0.5,
                                             per_solve_timeout=1e9) for _ in range(REPEATS)])
        scale = np.maximum(np.abs(reps[0]), np.finfo(float).tiny)
        rows.append({
            "coalition": coalition,
            "distinct_results": len({tuple(r.tolist()) for r in reps}),
            "max_relative_spread": float(np.max(np.ptp(reps, axis=0) / scale)),
            "max_absolute_spread": float(np.max(np.ptp(reps, axis=0))),
        })
    return rows


def main() -> None:
    blas = np.show_config(mode="dicts")["Build Dependencies"]["blas"]
    result = {
        "question": "Does the legacy engine reproduce itself on identical inputs?",
        "platform": {"numpy": np.__version__, "blas": blas.get("name"),
                     "thread_env": {k: os.environ.get(k) for k in
                                    ("OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "OMP_NUM_THREADS")}},
        "inputs": {"networks": NETWORKS, "repeats": REPEATS, "n_init": N_INIT, "seed": SEED,
                   "hill_n": 10.0, "t_max": 30.0, "t_tail": 10.0, "per_solve_timeout": "disabled"},
        "networks": {name: probe(name) for name in NETWORKS},
    }
    spreads = [r["max_relative_spread"] for rows in result["networks"].values() for r in rows]
    result["max_relative_spread"] = max(spreads)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(result, indent=2))
    print(f"wrote {OUT.relative_to(PROJECT_ROOT)}  blas={blas.get('name')}  "
          f"max relative spread {result['max_relative_spread']:.3g}")


if __name__ == "__main__":
    main()
