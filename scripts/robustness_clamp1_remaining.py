"""Run clamp-to-1 for the 3 models that didn't complete, append to summary."""
import json
import sys
from datetime import datetime, timezone

sys.path.insert(0, "/Users/elliottower/Documents/GitHub/epistasis-bench")
from scripts.robustness_clamp_to_1 import run_model, RESULTS_DIR

REMAINING = ["hematopoiesis_aging", "irons_cardiac", "calzone_cell_fate"]


def timestamp():
    return datetime.now(timezone.utc).isoformat()


if __name__ == "__main__":
    summary_path = f"{RESULTS_DIR}/robustness_clamp1_summary.json"
    with open(summary_path) as f:
        all_results = json.load(f)

    done = {r["model"] for r in all_results}
    todo = [m for m in REMAINING if m not in done]
    print(f"[{timestamp()}] {len(todo)} models remaining: {todo}")

    for model_name in todo:
        print(f"\n{'='*60}")
        print(f"[{timestamp()}] === {model_name} ===")
        print(f"{'='*60}")
        result = run_model(model_name, n_workers=10)
        all_results.append(result)

        with open(summary_path, "w") as f:
            json.dump(all_results, f, indent=2)
        print(f"[{timestamp()}] Saved ({len(all_results)}/27)")

    print(f"\n[{timestamp()}] Done. {len(all_results)} total models.")
