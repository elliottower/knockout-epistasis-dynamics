"""Does the composition gap require dynamics?

Registered in prereg_static_vs_dynamic_v1.md (commit 62d7a59).

The paper's proposed mechanism is that composition through TIME generates
higher-order structure absent from the wiring. The evidence offered is a
comparison with proteins (Brookes et al., PNAS 2022), where higher-order
epistasis is localised on the pairwise contact graph. That is a
cross-literature comparison, and proteins differ from gene regulatory networks
in many ways besides being static.

This holds the wiring fixed and builds a STATIC phenotype on the same graph, so
the dynamic/static contrast is measured inside one system.

The measured quantity is not whether higher-order energy exists -- both systems
have some -- but whether it is PREDICTABLE FROM LOCAL STRUCTURE, which is the
Brookes quantity.
"""

import glob
import json
import os

import numpy as np
from scipy.stats import spearmanr

RESULTS = "results/grn_v2"
OUT = "results/static_vs_dynamic.json"

# Registered: report at least two saturating nonlinearities and show the
# conclusion is invariant to the choice.
NONLINEARITIES = {
    "tanh": np.tanh,
    "logistic": lambda z: 1.0 / (1.0 + np.exp(-z)),
}


def walsh_spectrum(v, n):
    """Energy by interaction order from a value vector over 2^n coalitions."""
    c = v.astype(np.float64).copy()
    h = 1
    while h < len(c):
        for i in range(0, len(c), h * 2):
            for j in range(i, i + h):
                a, b = c[j], c[j + h]
                c[j], c[j + h] = a + b, a - b
        h *= 2
    c /= len(c)
    order = np.array([bin(i).count("1") for i in range(len(c))])
    spec = np.zeros(n + 1)
    for k in range(n + 1):
        spec[k] = float((c[order == k] ** 2).sum())
    tot = spec.sum()
    return (spec / tot if tot > 0 else spec), c


def static_phenotype(adj, node_names, n, phi):
    """phi(E(x)) over all 2^n knockout coalitions.

    E(x) = sum over edges of J_ij x_i x_j, with J signed by edge character.
    A knocked-out node has x = 0; a retained node has x = 1.
    """
    idx = {m: i for i, m in enumerate(node_names)}
    edges = []
    for tgt, regs in adj.items():
        for r in regs:
            src = r["regulator"]
            if src in idx and tgt in idx:
                edges.append((idx[src], idx[tgt], float(r["sign"])))
    if not edges:
        return None

    N = 2 ** n
    bits = ((np.arange(N)[:, None] >> np.arange(n)[None, :]) & 1).astype(float)
    E = np.zeros(N)
    for i, j, s in edges:
        E += s * bits[:, i] * bits[:, j]
    # Scale into the nonlinearity's responsive range; without this the
    # saturating function is effectively linear or fully saturated and the
    # test degenerates.
    sd = E.std()
    if sd > 0:
        E = E / sd
    return phi(E)


def wiring_prediction(adj, node_names, n, mode="star"):
    """Order-3+ coefficients predicted from the wiring alone.

    Brookes et al. (PNAS 2022) state the construction precisely: an rth-order
    interaction is active when r-1 positions are in contact with a CENTRAL
    position. That is a star, not general connectivity, and the distinction
    matters from order 4 upward.

    mode="star"      faithful to Brookes; a subset is active when some member is
                     adjacent to every other member.
    mode="connected" weaker variant, retained for comparison. It reports
                     correlations roughly four times lower, so using it would
                     understate how well wiring predicts higher-order terms.
    """
    idx = {m: i for i, m in enumerate(node_names)}
    nb = {i: set() for i in range(n)}
    for tgt, regs in adj.items():
        for r in regs:
            s = r["regulator"]
            if s in idx and tgt in idx:
                nb[idx[s]].add(idx[tgt])
                nb[idx[tgt]].add(idx[s])
    N = 2 ** n
    pred = np.zeros(N)
    for m in range(N):
        members = [i for i in range(n) if (m >> i) & 1]
        if len(members) < 3:
            continue
        ms = set(members)
        if mode == "star":
            pred[m] = 1.0 if any((ms - {c}) <= nb[c] for c in members) else 0.0
        else:
            seen, stack = {members[0]}, [members[0]]
            while stack:
                u = stack.pop()
                for w in nb[u] & ms:
                    if w not in seen:
                        seen.add(w)
                        stack.append(w)
            pred[m] = 1.0 if len(seen) == len(members) else 0.0
    return pred


def main():
    rows = []
    for cpath in sorted(glob.glob(f"{RESULTS}/*_composition_blind.json")):
        model = os.path.basename(cpath).replace("_composition_blind.json", "")
        wpath = f"{RESULTS}/{model}_wiring_blind.json"
        npz = f"{RESULTS}/{model}_coalition_blind.npz"
        if not (os.path.exists(wpath) and os.path.exists(npz)):
            continue

        wiring = json.load(open(wpath))
        adj = wiring["interaction_graph"]
        node_names = wiring["node_names"]
        d = np.load(npz, allow_pickle=True)
        key = "v_mean" if "v_mean" in d else d.files[0]
        v_dyn = np.asarray(d[key], dtype=float)
        if v_dyn.ndim > 1:
            v_dyn = v_dyn.mean(axis=1)
        n = int(np.log2(len(v_dyn)))
        if 2 ** n != len(v_dyn) or n > 18:
            continue

        pred = wiring_prediction(adj, node_names, n)
        order = np.array([bin(i).count("1") for i in range(2 ** n)])
        hi = order >= 3
        if hi.sum() < 10 or pred[hi].std() == 0:
            continue

        _, coef_dyn = walsh_spectrum(v_dyn, n)
        r_dyn, _ = spearmanr(np.abs(coef_dyn[hi]), pred[hi])

        row = {"model": model, "n": n, "dynamic_r": float(r_dyn)}
        for name, phi in NONLINEARITIES.items():
            v_st = static_phenotype(adj, node_names, n, phi)
            if v_st is None:
                continue
            _, coef_st = walsh_spectrum(v_st, n)
            r_st, _ = spearmanr(np.abs(coef_st[hi]), pred[hi])
            row[f"static_{name}_r"] = float(r_st)
        rows.append(row)
        print(f"  {model:28} n={n:2d} dyn r={r_dyn:+.3f} "
              + " ".join(f"{k.replace('static_','').replace('_r','')} "
                         f"r={row[k]:+.3f}" for k in row if k.startswith("static_")))

    if not rows:
        print("no scorable networks")
        return

    dyn = np.array([r["dynamic_r"] for r in rows])
    print(f"\n=== predictability of order-3+ coefficients from wiring "
          f"(n={len(rows)} networks) ===")
    print(f"{'system':<22} {'median r':>10} {'>=0.5':>8} {'<0.3':>8}")
    print(f"{'dynamic (attractor)':<22} {np.median(dyn):>10.3f} "
          f"{int((dyn >= 0.5).sum()):>8} {int((dyn < 0.3).sum()):>8}")
    stats = {}
    for name in NONLINEARITIES:
        k = f"static_{name}_r"
        st = np.array([r[k] for r in rows if k in r])
        stats[name] = st
        print(f"{'static (' + name + ')':<22} {np.median(st):>10.3f} "
              f"{int((st >= 0.5).sum()):>8} {int((st < 0.3).sum()):>8}")

    print(f"\n=== registered verdicts ===")
    verdicts = {}
    for name, st in stats.items():
        h1 = int((st >= 0.5).sum()) >= 20
        h2 = bool(np.median(dyn) < 0.3)
        paired = st - dyn[: len(st)]
        h3 = bool(np.median(paired) > 0) and int((paired > 0).sum()) > len(paired) / 2
        h4 = not h1
        print(f"  [{name}]")
        print(f"    H1 static predictable (>=20/28 with r>=0.5): "
              f"{'HOLDS' if h1 else 'FAILS'} ({int((st>=0.5).sum())}/{len(st)})")
        print(f"    H2 dynamic NOT predictable (median r<0.3):   "
              f"{'HOLDS' if h2 else 'FAILS'} (median {np.median(dyn):+.3f})")
        print(f"    H3 paired difference favours static:         "
              f"{'HOLDS' if h3 else 'FAILS'} "
              f"({int((paired>0).sum())}/{len(paired)} networks)")
        print(f"    H4 FALSIFIER (static also unpredictable):    "
              f"{'FIRES - dynamics is not the mechanism' if h4 else 'does not fire'}")
        verdicts[name] = {"H1": bool(h1), "H2": bool(h2), "H3": bool(h3),
                          "H4_falsifier_fires": bool(h4)}

    json.dump({"prereg": "prereg_static_vs_dynamic_v1.md",
               "prereg_commit": "62d7a59",
               "n_networks": len(rows), "per_model": rows,
               "verdicts": verdicts}, open(OUT, "w"), indent=2)
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
