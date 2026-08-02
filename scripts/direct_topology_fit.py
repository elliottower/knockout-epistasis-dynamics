"""Direct fit of topological features to the composition gap.

Registered in prereg_direct_topology_fit_v1.md (commit 05b33a9).

The paper establishes "the gap is not predictable from standard topological
features" from the failure of pre-registered predictions (10/21 = 48%). That
inference cannot separate "the features carry no information" from "the
predictor was poor". This fits the features directly, so no predictor is
involved.

Feature families are taken from the paper's own text -- feedback loop counts,
rule complexity, in-degree distribution -- not chosen here.
"""

import glob
import itertools
import json
import os

import numpy as np
from scipy.stats import spearmanr
from sklearn.linear_model import RidgeCV
from sklearn.model_selection import LeaveOneOut
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

RESULTS = "results/grn_v2"
OUT = "results/direct_topology_fit.json"
MAX_CYCLE_LEN = 4


def order3plus(spectrum):
    """Fraction of spectral energy at order 3 and above."""
    return float(sum(spectrum[2:]))


def cycles_up_to(adj, max_len):
    """Count simple directed cycles up to max_len, by sign.

    adj maps node -> list of (regulator, sign). Edge direction is
    regulator -> node.
    """
    nodes = list(adj)
    idx = {n: i for i, n in enumerate(nodes)}
    out = {i: [] for i in range(len(nodes))}
    for tgt, regs in adj.items():
        for r in regs:
            src = r["regulator"]
            if src in idx:
                out[idx[src]].append((idx[tgt], r["sign"]))

    pos = neg = total = 0
    for start in range(len(nodes)):
        # paths from start, only visiting nodes >= start to count each cycle once
        stack = [(start, [start], 1)]
        while stack:
            node, path, sign = stack.pop()
            if len(path) > max_len:
                continue
            for nxt, s in out[node]:
                if nxt == start and len(path) >= 1:
                    total += 1
                    if sign * s > 0:
                        pos += 1
                    else:
                        neg += 1
                elif nxt > start and nxt not in path:
                    stack.append((nxt, path + [nxt], sign * s))
    return pos, neg, total


def features_for(wiring):
    """Feature vector from a wiring JSON, grouped by the paper's families."""
    graph = wiring["interaction_graph"]
    fourier = wiring.get("rule_fourier", {})
    per_gene = fourier.get("per_gene", {}) or {}

    indeg = np.array([len(v) for v in graph.values()], dtype=float)
    n = len(graph)
    n_edges = float(indeg.sum())

    pos, neg, total = cycles_up_to(graph, MAX_CYCLE_LEN)

    # rule complexity: how many interaction terms each gene's rule carries.
    # per_gene entries are {"regulators": [...], "interactions": [...]}.
    inter_counts, reg_counts = [], []
    for g in graph:
        entry = per_gene.get(g, {}) or {}
        inter_counts.append(len(entry.get("interactions", []) or []))
        reg_counts.append(len(entry.get("regulators", []) or []))
    inter_counts = np.array(inter_counts, dtype=float)
    reg_counts = np.array(reg_counts, dtype=float)
    nonlinear = inter_counts > 0
    n_pairwise = float(fourier.get("n_pairwise", 0) or 0)
    n_triples = float(fourier.get("n_triples", 0) or 0)

    return {
        # in-degree distribution
        "indeg_mean": float(indeg.mean()),
        "indeg_max": float(indeg.max()),
        "indeg_var": float(indeg.var()),
        "edge_density": float(n_edges / (n * n)) if n else 0.0,
        # feedback loops
        "cycles_pos": float(pos),
        "cycles_neg": float(neg),
        "cycles_total": float(total),
        # rule complexity
        "rule_inter_mean": float(inter_counts.mean()) if len(inter_counts) else 0.0,
        "rule_inter_max": float(inter_counts.max()) if len(inter_counts) else 0.0,
        "rule_reg_mean": float(reg_counts.mean()) if len(reg_counts) else 0.0,
        "n_pairwise_terms": n_pairwise,
        "n_triple_terms": n_triples,
        "frac_nonlinear": float(nonlinear.mean()) if len(nonlinear) else 0.0,
        # size, as a stated nuisance covariate
        "n_nodes": float(n),
    }


def holm(pvals):
    order = np.argsort(pvals)
    m = len(pvals)
    adj = np.empty(m)
    running = 0.0
    for rank, i in enumerate(order):
        running = max(running, (m - rank) * pvals[i])
        adj[i] = min(1.0, running)
    return adj


def main():
    rows, names = [], []
    for cpath in sorted(glob.glob(f"{RESULTS}/*_composition_blind.json")):
        model = os.path.basename(cpath).replace("_composition_blind.json", "")
        wpath = f"{RESULTS}/{model}_wiring_blind.json"
        if not os.path.exists(wpath):
            print(f"skip {model}: no wiring file")
            continue
        comp = json.load(open(cpath))
        wiring = json.load(open(wpath))
        spec = comp["energy_spectrum"]
        gap = order3plus(spec["global"]) - order3plus(spec["local_rules"])
        feats = features_for(wiring)
        rows.append(feats)
        names.append((model, gap))

    feat_names = sorted(rows[0])
    X = np.array([[r[k] for k in feat_names] for r in rows])
    y = np.array([g for _, g in names])
    print(f"networks: {len(y)}   features: {len(feat_names)}")
    print(f"gap range: {y.min():+.4f} to {y.max():+.4f}   mean {y.mean():+.4f}")

    # H1: leave-one-out cross-validated R^2 of the joint fit
    model = make_pipeline(
        StandardScaler(), RidgeCV(alphas=np.logspace(-3, 3, 25))
    )
    preds = np.empty_like(y)
    for tr, te in LeaveOneOut().split(X):
        model.fit(X[tr], y[tr])
        preds[te] = model.predict(X[te])
    ss_res = float(((y - preds) ** 2).sum())
    ss_tot = float(((y - y.mean()) ** 2).sum())
    cv_r2 = 1.0 - ss_res / ss_tot

    # H2: per-feature Spearman, Holm corrected
    rhos, praw = [], []
    for j in range(X.shape[1]):
        r, p = spearmanr(X[:, j], y)
        rhos.append(float(r))
        praw.append(float(p))
    padj = holm(np.array(praw))

    print(f"\nH1  LOO cross-validated R^2 = {cv_r2:+.4f}   (predicted < 0.20)")
    print(f"    in-sample R^2 (reference) = "
          f"{model.fit(X, y).score(X, y):+.4f}")
    print(f"\nH2  per-feature Spearman vs gap (Holm-adjusted p):")
    order = np.argsort(-np.abs(rhos))
    for j in order:
        star = "  <-- |rho| >= 0.5" if abs(rhos[j]) >= 0.5 else ""
        print(f"    {feat_names[j]:<18} rho={rhos[j]:+.3f}  "
              f"p_adj={padj[j]:.3f}{star}")

    h1 = cv_r2 < 0.20
    h3 = cv_r2 >= 0.40
    sig = [feat_names[j] for j in range(len(feat_names))
           if abs(rhos[j]) >= 0.5 and padj[j] < 0.05]
    h2 = len(sig) == 0
    h4 = abs(rhos[feat_names.index("n_nodes")]) < 0.3

    print(f"\n=== registered verdicts ===")
    print(f"H1 (cv R^2 < 0.20):            {'HOLDS' if h1 else 'FAILS'}")
    print(f"H2 (no feature |rho|>=0.5):    {'HOLDS' if h2 else 'FAILS'}"
          + (f"  -> {sig}" if sig else ""))
    print(f"H3 FALSIFIER (cv R^2 >= 0.40): {'FIRES - published claim is wrong' if h3 else 'does not fire'}")
    print(f"H4 (n_nodes |rho| < 0.3):      {'HOLDS' if h4 else 'FAILS'}")

    out = {
        "prereg": "prereg_direct_topology_fit_v1.md",
        "prereg_commit": "05b33a9",
        "n_networks": int(len(y)),
        "features": feat_names,
        "cv_r2_loo": cv_r2,
        "spearman": {feat_names[j]: {"rho": rhos[j], "p_raw": praw[j],
                                     "p_holm": float(padj[j])}
                     for j in range(len(feat_names))},
        "gap_by_model": {m: float(g) for m, g in names},
        "verdicts": {"H1": bool(h1), "H2": bool(h2),
                     "H3_falsifier_fires": bool(h3), "H4": bool(h4)},
    }
    os.makedirs("results", exist_ok=True)
    json.dump(out, open(OUT, "w"), indent=2)
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
