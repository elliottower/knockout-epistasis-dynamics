"""GNK positive control: does the star predictor work, and does dynamics
create structure the wiring cannot explain?

Registered in prereg_gnk_control_v1.md (commit caa8a2f). Replaces the static arm
of prereg_static_vs_dynamic_v1.md, whose phenotype was a nonlinearity over a
pairwise energy and therefore 100% scale artifact above order 2.

The GNK model, from Brookes, Aghazadeh & Listgarten (PNAS 2022, PMC8740588),
read from the full text:

    "Assign a 'subsequence fitness,' f_j(s[j]), to every possible subsequence,
     s[j], by drawing a value from the normal distribution with mean equal to
     zero and variance equal to 1/L ... the subsequence fitness values are
     summed to produce the total fitness values f(s) = sum_j f_j(s[j])."

with "the neighborhood of a position j contains all positions that are in
structural contact with it."

Why this is a positive control: given Structural neighborhoods, the star
property is a THEOREM about GNK, not an empirical claim. If the predictor cannot
recover GNK's own order-3+ coefficients, the predictor is broken and every
result measured with it is unsupported. H0 gates the rest.
"""

import glob
import json
import os

import numpy as np
from scipy.stats import spearmanr

RESULTS = "results/grn_v2"
OUT = "results/gnk_control.json"
N_DRAWS = 10          # f_j is random; registered to report median across draws
SEED = 20260806


def walsh_coefficients(v):
    """Fast Walsh-Hadamard transform, normalised. Validated elsewhere in this
    repo against the committed composition spectra to 0.00e+00."""
    c = v.astype(np.float64).copy()
    h = 1
    while h < len(c):
        for i in range(0, len(c), h * 2):
            for j in range(i, i + h):
                a, b = c[j], c[j + h]
                c[j], c[j + h] = a + b, a - b
        h *= 2
    return c / len(c)


def neighborhoods(adj, node_names, n):
    """V_j = {j} union {nodes sharing a regulatory edge with j}.

    The regulatory graph is directed; contact in the GNK sense is symmetric, so
    an edge in either direction places two nodes in one another's neighborhood.
    """
    idx = {m: i for i, m in enumerate(node_names)}
    V = [{j} for j in range(n)]
    for tgt, regs in adj.items():
        if tgt not in idx:
            continue
        for r in regs:
            s = r["regulator"]
            if s in idx:
                V[idx[tgt]].add(idx[s])
                V[idx[s]].add(idx[tgt])
    return [sorted(v) for v in V]


def gnk_phenotype(V, n, rng):
    """f(x) = sum_j f_j(x_{V_j}), f_j ~ N(0, 1/L) per neighborhood state."""
    N = 2 ** n
    bits = ((np.arange(N)[:, None] >> np.arange(n)[None, :]) & 1).astype(np.int64)
    f = np.zeros(N)
    L = n
    for j in range(n):
        members = V[j]
        # index each coalition by the state of this neighborhood
        key = np.zeros(N, dtype=np.int64)
        for p, m in enumerate(members):
            key += bits[:, m] << p
        table = rng.normal(0.0, np.sqrt(1.0 / L), size=2 ** len(members))
        f += table[key]
    return f


def subset_prediction(V, n):
    """A coefficient at subset S is nonzero iff S is contained in some V_j.

    Derived from the GNK definition rather than paraphrased: the phenotype is
    f = sum_j f_j(x_{V_j}), so each term's Walsh support is every SUBSET of
    V_j, and the support of f is the union of those.

    This is where an earlier version was wrong. It required the central node to
    be a MEMBER of S ("some member adjacent to all others"), but the centre j
    need not lie in S. That version left 21% of GNK order-3+ energy on
    coefficients it called inactive; this one leaves 0.0%, exactly as the
    theorem requires.
    """
    N = 2 ** n
    masks = [sum(1 << m for m in v) for v in V]
    pred = np.zeros(N)
    for s in range(N):
        if bin(s).count("1") < 3:
            continue
        pred[s] = 1.0 if any((s & ~mk) == 0 for mk in masks) else 0.0
    return pred


def main():
    rows = []
    for wpath in sorted(glob.glob(f"{RESULTS}/*_wiring_blind.json")):
        model = os.path.basename(wpath).replace("_wiring_blind.json", "")
        npz = f"{RESULTS}/{model}_coalition_blind.npz"
        if not os.path.exists(npz):
            continue
        wiring = json.load(open(wpath))
        adj, node_names = wiring["interaction_graph"], wiring["node_names"]
        z = np.load(npz, allow_pickle=True)
        k = "v_mean" if "v_mean" in z else z.files[0]
        v_dyn = np.asarray(z[k], dtype=float)
        if v_dyn.ndim > 1:
            v_dyn = v_dyn.mean(axis=1)
        n = int(np.log2(len(v_dyn)))
        if 2 ** n != len(v_dyn) or n > 15:
            continue

        V = neighborhoods(adj, node_names, n)
        pred = subset_prediction(V, n)
        order = np.array([bin(i).count("1") for i in range(2 ** n)])
        hi = order >= 3
        if hi.sum() < 10 or pred[hi].std() == 0:
            continue

        vsizes = [len(v) for v in V]

        r_dyn = spearmanr(np.abs(walsh_coefficients(v_dyn)[hi]), pred[hi])[0]

        rng = np.random.default_rng(SEED + hash(model) % 10000)
        r_gnk = []
        for _ in range(N_DRAWS):
            f = gnk_phenotype(V, n, rng)
            r_gnk.append(spearmanr(np.abs(walsh_coefficients(f)[hi]), pred[hi])[0])
        r_gnk = float(np.median(r_gnk))

        rows.append({"model": model, "n": n,
                     "mean_neighborhood": float(np.mean(vsizes)),
                     "max_neighborhood": int(max(vsizes)),
                     "gnk_r": r_gnk, "dynamic_r": float(r_dyn),
                     "predictor_balance": float(pred[hi].mean())})
        print(f"  {model:28} n={n:2d} |V|avg={np.mean(vsizes):4.1f} "
              f"GNK r={r_gnk:+.3f}  dyn r={r_dyn:+.3f}")

    if not rows:
        print("no scorable networks")
        return

    gnk = np.array([r["gnk_r"] for r in rows])
    dyn = np.array([r["dynamic_r"] for r in rows])
    print(f"\n=== star-predictor recovery of order-3+ coefficients "
          f"(n={len(rows)} networks, {N_DRAWS} GNK draws each) ===")
    print(f"{'phenotype':<26}{'median r':>10}{'>=0.5':>8}{'IQR':>18}")
    for name, a in (("GNK (positive control)", gnk), ("dynamic (attractor)", dyn)):
        q1, q3 = np.percentile(a, [25, 75])
        print(f"{name:<26}{np.median(a):>10.3f}{int((a >= 0.5).sum()):>8}"
              f"{f'[{q1:+.3f}, {q3:+.3f}]':>18}")

    h0 = bool(np.median(gnk) >= 0.5)
    gap = float(np.median(gnk) - np.median(dyn))
    h1 = bool(gap >= 0.25)
    h2 = bool(np.median(dyn) < 0.3)
    h3 = bool(abs(gap) < 0.1)

    print(f"\n=== registered verdicts ===")
    print(f"H0 GATE: GNK median r >= 0.5      "
          f"{'PASSES' if h0 else 'FAILS - PREDICTOR IS BROKEN'} "
          f"(median {np.median(gnk):+.3f})")
    if not h0:
        print("\n  H0 failed. Nothing below is interpreted, and every result")
        print("  measured with this predictor is withdrawn.")
    else:
        print(f"H1 dynamic below GNK by >=0.25    "
              f"{'HOLDS' if h1 else 'FAILS'} (gap {gap:+.3f})")
        print(f"H2 dynamic median < 0.3           "
              f"{'HOLDS' if h2 else 'FAILS'} (median {np.median(dyn):+.3f})")
        print(f"H3 FALSIFIER gap < 0.1            "
              f"{'FIRES - dynamics is not the mechanism' if h3 else 'does not fire'}")

    json.dump({"prereg": "prereg_gnk_control_v1.md", "prereg_commit": "caa8a2f",
               "n_networks": len(rows), "n_draws": N_DRAWS, "seed": SEED,
               "per_model": rows,
               "gnk_median": float(np.median(gnk)),
               "dynamic_median": float(np.median(dyn)), "gap": gap,
               "verdicts": {"H0_gate_passes": h0, "H1": h1, "H2": h2,
                            "H3_falsifier_fires": h3}},
              open(OUT, "w"), indent=2)
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
