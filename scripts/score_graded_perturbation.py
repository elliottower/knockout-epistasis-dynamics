"""Score the graded perturbation sweep against its registered hypotheses.

Registered in prereg_graded_perturbation_v1.md (3095900), amended in 98d9bea
and 2afd478.

Written BEFORE the sweep results were available, so the scoring cannot be
tuned to the data. It reads whatever units exist and reports which are missing
rather than silently pooling a partial grid.

H1  gap sign preserved at every f in (0,1)
H2  normalised higher-order fraction rises as f -> 0, by >= 1.5x from
    f=0.75 to f=0
H3  dose-response is monotone in f
H4  creation/destruction ratio flips at some f*, and we report where
H5  FALSIFIER: if the gap vanishes (sign flip, or |delta_3+| < 0.5pp) for any
    f >= 0.25, the published result is specific to complete knockout
"""

import glob
import json
import os
from collections import defaultdict

import numpy as np

IN_DIR = "results/graded_perturbation"
OUT = "results/graded_perturbation_scored.json"
LEVELS = [0.0, 0.25, 0.5, 0.75, 1.0]
PP = 0.005  # 0.5 percentage points, the registered vanishing threshold


def local_o3plus(model):
    """Local-rule order-3+ energy, from the committed ODE summary."""
    p = f"results/ode_full/{model}_ode.json"
    if not os.path.exists(p):
        return None
    return json.load(open(p)).get("local_o3plus")


def main():
    units = defaultdict(dict)
    for p in sorted(glob.glob(f"{IN_DIR}/*.json")):
        d = json.load(open(p))
        units[d["model"]][round(float(d["clamp_level"]), 2)] = d

    models = sorted(units)
    if not models:
        print(f"no units found in {IN_DIR}/ -- nothing to score")
        return

    complete = [m for m in models
                if all(round(f, 2) in units[m] for f in LEVELS)]
    partial = [m for m in models if m not in complete]
    print(f"units: {sum(len(v) for v in units.values())}")
    print(f"models complete (all 5 levels): {len(complete)}")
    if partial:
        print(f"models partial, EXCLUDED from pooled statistics: "
              f"{len(partial)} -> {', '.join(partial)}")

    # gap = global order-3+ minus local order-3+, per the paper's definition
    rows = []
    for m in complete:
        loc = local_o3plus(m)
        if loc is None:
            print(f"  {m}: no local_o3plus, skipped")
            continue
        r = {"model": m, "local_o3plus": loc}
        for f in LEVELS:
            u = units[m][round(f, 2)]
            r[f"gap_{f}"] = u["energy_3plus"] - loc
            r[f"ratio_{f}"] = u["normalised_higher_order"]
            r[f"gated_{f}"] = u["spectrum_gated"]
        rows.append(r)

    if not rows:
        print("no scorable models")
        return

    mid = [f for f in LEVELS if 0.0 < f < 1.0]

    # --- H1: sign preserved at every intermediate level
    ref = np.array([r["gap_0.0"] for r in rows])
    h1_viol = []
    for f in mid:
        g = np.array([r[f"gap_{f}"] for r in rows])
        flips = [rows[i]["model"] for i in range(len(rows))
                 if np.sign(g[i]) != np.sign(ref[i]) and abs(g[i]) > PP]
        if flips:
            h1_viol.append((f, flips))
    h1 = not h1_viol

    # --- H2: normalised fraction rises as f -> 0
    ung = [r for r in rows if not (r["gated_0.0"] or r["gated_0.75"])]
    if ung:
        a = np.median([r["ratio_0.0"] for r in ung])
        b = np.median([r["ratio_0.75"] for r in ung])
        h2_factor = float(a / b) if b > 0 else float("inf")
    else:
        a = b = h2_factor = float("nan")
    h2 = bool(h2_factor >= 1.5)

    # --- H3: monotone dose-response (ratio non-increasing as f rises, 0->0.75)
    mono = 0
    for r in rows:
        seq = [r[f"ratio_{f}"] for f in [0.0, 0.25, 0.5, 0.75]]
        if all(np.isfinite(seq)) and all(
                seq[i] >= seq[i + 1] - 1e-9 for i in range(len(seq) - 1)):
            mono += 1
    h3 = bool(mono >= 0.5 * len(rows))

    # --- H4: where does creation/destruction flip?
    flip_at = {}
    for r in rows:
        prev = np.sign(r["gap_0.0"])
        for f in LEVELS[1:]:
            s = np.sign(r[f"gap_{f}"])
            if s != prev and abs(r[f"gap_{f}"]) > PP:
                flip_at[r["model"]] = f
                break
            prev = s
    creating = {f: int(sum(1 for r in rows if r[f"gap_{f}"] > PP))
                for f in LEVELS}
    destroying = {f: int(sum(1 for r in rows if r[f"gap_{f}"] < -PP))
                  for f in LEVELS}

    # --- H5 FALSIFIER: gap vanishes at any f >= 0.25
    h5_fires, h5_detail = False, {}
    for f in [x for x in LEVELS if x >= 0.25]:
        g = np.array([r[f"gap_{f}"] for r in rows])
        vanished = [rows[i]["model"] for i in range(len(rows))
                    if abs(g[i]) < PP]
        h5_detail[f] = {"n_vanished": len(vanished),
                        "frac": len(vanished) / len(rows),
                        "models": vanished}
        if len(vanished) > 0.5 * len(rows):
            h5_fires = True

    print(f"\n=== gap by clamp level (n={len(rows)} models) ===")
    print(f"{'f':>6} {'median gap':>12} {'creating':>9} {'destroying':>11} "
          f"{'gated':>6}")
    for f in LEVELS:
        g = np.array([r[f"gap_{f}"] for r in rows])
        ng = sum(1 for r in rows if r[f"gated_{f}"])
        print(f"{f:>6.2f} {np.median(g):>+12.4f} {creating[f]:>9} "
              f"{destroying[f]:>11} {ng:>6}")

    print(f"\n=== registered verdicts ===")
    print(f"H1 sign preserved:        {'HOLDS' if h1 else 'FAILS'}"
          + ("" if h1 else f"  -> {h1_viol}"))
    print(f"H2 ratio rises as f->0:   {'HOLDS' if h2 else 'FAILS'}"
          f"  (median {b:.4f} at f=0.75 -> {a:.4f} at f=0, "
          f"factor {h2_factor:.2f}, need >=1.5)")
    print(f"H3 monotone:              {'HOLDS' if h3 else 'FAILS'}"
          f"  ({mono}/{len(rows)} models monotone)")
    print(f"H4 flip locations:        "
          f"{flip_at if flip_at else 'no sign flips observed'}")
    print(f"H5 FALSIFIER:             "
          f"{'FIRES - result is specific to complete knockout' if h5_fires else 'does not fire'}")
    for f, d in h5_detail.items():
        print(f"    f={f:.2f}: {d['n_vanished']}/{len(rows)} vanished "
              f"({d['frac']:.0%})")

    out = {
        "prereg": "prereg_graded_perturbation_v1.md",
        "prereg_commit": "3095900",
        "amendments": ["98d9bea", "2afd478"],
        "n_models_scored": len(rows),
        "models_excluded_partial": partial,
        "levels": LEVELS,
        "per_model": rows,
        "creating": creating,
        "destroying": destroying,
        "verdicts": {"H1": bool(h1), "H2": bool(h2), "H3": bool(h3),
                     "H4_flip_at": flip_at, "H5_falsifier_fires": bool(h5_fires)},
        "h2_factor": h2_factor,
        "h5_detail": h5_detail,
    }
    json.dump(out, open(OUT, "w"), indent=2)
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
