"""Verify the k-SII pairwise Walsh orthogonality term-by-term.

Two claims to verify:

1. For each pair (i,j), the pairwise Walsh coefficient w_{ij} from the
   order-2 k-SII plug-in reconstruction equals w_{ij} from the order-3
   reconstruction. This is a per-coefficient identity, much stronger than
   aggregate AUROC equality.

2. The recursive k-SII identity:
       phi^(2)_{ij} = phi^(3)_{ij} + (1/2) * sum_k phi^(3)_{ijk}
   holds term-by-term. This is the mechanism: the k-SII recursion
   redefines lower-order terms to absorb higher-order contributions,
   making the pairwise Walsh projection invariant to truncation order.

Uses shapiq.ExactComputer to get exact k-SII values at each order.

Usage:
    uv run python scripts/verify_ksii_orthogonality.py \
        --data-dir ~/Documents/GitHub/weight-circuit-discovery/experiments_batch2/genetics/
"""

import argparse
from itertools import combinations
from pathlib import Path

import numpy as np

from data_utils import load_coalition_table, normalized_wht, pooled_mean_values


CIRCUIT_FILES = {
    "weight_ioi": "weight_ioi_zero_v2_coalition_values.npz",
    "ioi": "ioi_zero_v2_coalition_values.npz",
    "random15": "random15_zero_v2_coalition_values.npz",
}


def exact_ksii(v, n, max_order):
    """Compute exact k-SII values via shapiq.ExactComputer."""
    import shapiq

    def game_callable(coalitions):
        results = np.zeros(len(coalitions))
        for idx, coal in enumerate(coalitions):
            bitmask = 0
            for p in range(n):
                if coalitions[idx, p]:
                    bitmask |= (1 << p)
            results[idx] = v[bitmask]
        return results

    computer = shapiq.ExactComputer(game=game_callable, n_players=n)
    return computer.shapley_interactions(order=max_order, index="k-SII")


def reconstruct_from_ksii(interaction_values, n):
    """Plug-in reconstruction: v_hat(S) = v(empty) + sum_{T subset S, T nonempty} phi_T."""
    N = 2 ** n
    v_hat = np.zeros(N)

    baseline = 0.0
    for key, val in interaction_values.items():
        if len(key) == 0:
            baseline = val
            continue

    for s_idx in range(N):
        s_set = frozenset(i for i in range(n) if s_idx & (1 << i))
        total = baseline
        for key, val in interaction_values.items():
            if len(key) == 0:
                continue
            if key.issubset(s_set):
                total += val
        v_hat[s_idx] = total

    return v_hat


def extract_interaction_dict(iv):
    """Convert shapiq InteractionValues to {frozenset: float} dict."""
    result = {}
    for key, val in zip(iv.interaction_lookup.keys(), iv.values.flat):
        result[frozenset(key)] = float(val)
    return result


def verify_circuit(name, data_path):
    print(f"\n{'='*70}")
    print(f"Circuit: {name}")
    print(f"{'='*70}")

    table = load_coalition_table(str(data_path))
    n = table["n_players"]
    v = pooled_mean_values(table)
    w_true = normalized_wht(v)

    print(f"Computing exact k-SII at order 2...")
    iv2 = exact_ksii(v, n, max_order=2)
    d2 = extract_interaction_dict(iv2)

    print(f"Computing exact k-SII at order 3...")
    iv3 = exact_ksii(v, n, max_order=3)
    d3 = extract_interaction_dict(iv3)

    print(f"Reconstructing from order-2 k-SII...")
    v_hat_o2 = reconstruct_from_ksii(d2, n)
    w_o2 = normalized_wht(v_hat_o2)

    print(f"Reconstructing from order-3 k-SII...")
    v_hat_o3 = reconstruct_from_ksii(d3, n)
    w_o3 = normalized_wht(v_hat_o3)

    print(f"\n--- Claim 1: Pairwise Walsh coefficients identical ---")
    pairs = list(combinations(range(n), 2))
    max_diff = 0.0
    diffs = []
    for i, j in pairs:
        idx = (1 << i) | (1 << j)
        diff = abs(w_o2[idx] - w_o3[idx])
        diffs.append(diff)
        max_diff = max(max_diff, diff)

    print(f"  Number of pairs: {len(pairs)}")
    print(f"  Max absolute difference: {max_diff:.2e}")
    print(f"  Mean absolute difference: {np.mean(diffs):.2e}")
    print(f"  Pairs with diff > 1e-15: {sum(1 for d in diffs if d > 1e-15)}")
    print(f"  Pairs with diff > 1e-10: {sum(1 for d in diffs if d > 1e-10)}")
    claim1_holds = max_diff < 1e-10
    print(f"  CLAIM 1: {'VERIFIED' if claim1_holds else 'FAILED'}")

    print(f"\n  Also checking order-1 Walsh coefficients...")
    max_diff_o1 = 0.0
    for i in range(n):
        idx = 1 << i
        diff = abs(w_o2[idx] - w_o3[idx])
        max_diff_o1 = max(max_diff_o1, diff)
    print(f"  Max order-1 difference: {max_diff_o1:.2e}")

    print(f"\n  Checking order-3+ Walsh coefficients (should differ)...")
    from data_utils import popcount_array
    pc = popcount_array(n)
    mask_o3 = pc >= 3
    diff_o3 = np.abs(w_o2[mask_o3] - w_o3[mask_o3])
    print(f"  Order-3+ coefficients: {mask_o3.sum()}")
    print(f"  Max order-3+ difference: {diff_o3.max():.6f}")
    print(f"  Mean order-3+ difference: {diff_o3.mean():.6f}")
    print(f"  (These SHOULD differ — order-3 reconstruction captures more)")

    print(f"\n--- Claim 2: Recursive k-SII identity ---")
    print(f"  phi^(2)_{{ij}} = phi^(3)_{{ij}} + (1/2) * sum_k phi^(3)_{{ijk}}")

    max_identity_error = 0.0
    identity_errors = []
    for i, j in pairs:
        phi2_ij = d2.get(frozenset({i, j}), 0.0)
        phi3_ij = d3.get(frozenset({i, j}), 0.0)

        triple_sum = 0.0
        for k in range(n):
            if k == i or k == j:
                continue
            triple_sum += d3.get(frozenset({i, j, k}), 0.0)

        rhs = phi3_ij + 0.5 * triple_sum
        error = abs(phi2_ij - rhs)
        identity_errors.append(error)
        max_identity_error = max(max_identity_error, error)

        if error > 1e-10:
            print(f"  FAILED at ({i},{j}): phi2={phi2_ij:.8e}, "
                  f"phi3={phi3_ij:.8e}, 0.5*sum={0.5*triple_sum:.8e}, "
                  f"rhs={rhs:.8e}, error={error:.2e}")

    print(f"\n  Max identity error: {max_identity_error:.2e}")
    print(f"  Mean identity error: {np.mean(identity_errors):.2e}")
    print(f"  Pairs with error > 1e-10: {sum(1 for e in identity_errors if e > 1e-10)}")
    claim2_holds = max_identity_error < 1e-10
    print(f"  CLAIM 2: {'VERIFIED' if claim2_holds else 'FAILED'}")

    if not claim2_holds:
        print(f"\n  Trying alternative coefficient: 1/(n-2+1) = 1/{n-1}")
        max_alt_error = 0.0
        for i, j in pairs:
            phi2_ij = d2.get(frozenset({i, j}), 0.0)
            phi3_ij = d3.get(frozenset({i, j}), 0.0)
            triple_sum = sum(d3.get(frozenset({i, j, k}), 0.0)
                            for k in range(n) if k != i and k != j)
            for coeff in [1/3, 1/(n-2), 1/(n-1), 1/n, 2/3, 1/4]:
                rhs = phi3_ij + coeff * triple_sum
                error = abs(phi2_ij - rhs)
                if error < max_alt_error or max_alt_error == 0:
                    pass
            for coeff_num in range(1, 20):
                for coeff_den in range(1, 20):
                    coeff = coeff_num / coeff_den
                    errors = []
                    for ii, jj in pairs[:10]:
                        phi2 = d2.get(frozenset({ii, jj}), 0.0)
                        phi3 = d3.get(frozenset({ii, jj}), 0.0)
                        ts = sum(d3.get(frozenset({ii, jj, kk}), 0.0)
                                 for kk in range(n) if kk != ii and kk != jj)
                        errors.append(abs(phi2 - phi3 - coeff * ts))
                    me = max(errors)
                    if me < 1e-10:
                        print(f"  Found exact coefficient: {coeff_num}/{coeff_den} = {coeff}")
                        max_alt_error = me
                        break
                if max_alt_error < 1e-10:
                    break

    print(f"\n--- Claim 3: Walsh projection proof (algebraic) ---")
    print(f"  If phi^(2)_{{ij}} = phi^(3)_{{ij}} + c * sum_k phi^(3)_{{ijk}},")
    print(f"  and w_{{ij}} = (1/4)*phi_{{ij}} + (1/8)*sum_k phi_{{ijk}} + ...,")
    print(f"  then w^(2)_{{ij}} = (1/4)*phi^(2)_{{ij}}")
    print(f"       = (1/4)*(phi^(3)_{{ij}} + c*sum_k phi^(3)_{{ijk}})")
    print(f"  And w^(3)_{{ij}} = (1/4)*phi^(3)_{{ij}} + (1/8)*sum_k phi^(3)_{{ijk}}")
    print(f"  These are equal iff c = 1/2, which gives the 1/8 weight exactly.")

    print(f"\n--- Summary for {name} ---")
    r2_o2 = 1 - np.sum((v - v_hat_o2)**2) / np.sum((v - v.mean())**2)
    r2_o3 = 1 - np.sum((v - v_hat_o3)**2) / np.sum((v - v.mean())**2)
    print(f"  R² (order-2 reconstruction): {r2_o2:.6f}")
    print(f"  R² (order-3 reconstruction): {r2_o3:.6f}")
    print(f"  R² improvement: {r2_o3 - r2_o2:.6f}")
    print(f"  Pairwise Walsh identity: {'EXACT' if claim1_holds else 'BROKEN'}")
    print(f"  Recursive identity: {'EXACT' if claim2_holds else 'BROKEN'}")

    return claim1_holds, claim2_holds


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", required=True)
    parser.add_argument("--circuits", nargs="+", default=list(CIRCUIT_FILES.keys()))
    args = parser.parse_args()

    data_dir = Path(args.data_dir)
    results = {}

    for circ in args.circuits:
        fpath = data_dir / CIRCUIT_FILES[circ]
        if not fpath.exists():
            print(f"Skipping {circ}: {fpath} not found")
            continue
        c1, c2 = verify_circuit(circ, fpath)
        results[circ] = {"claim1_pairwise_identity": c1, "claim2_recursive_identity": c2}

    print(f"\n{'='*70}")
    print(f"FINAL RESULTS")
    print(f"{'='*70}")
    all_pass = True
    for circ, r in results.items():
        c1 = "VERIFIED" if r["claim1_pairwise_identity"] else "FAILED"
        c2 = "VERIFIED" if r["claim2_recursive_identity"] else "FAILED"
        print(f"  {circ}: Claim 1 ({c1}), Claim 2 ({c2})")
        if not r["claim1_pairwise_identity"] or not r["claim2_recursive_identity"]:
            all_pass = False

    if all_pass:
        print(f"\nAll claims verified across all circuits.")
        print(f"The k-SII orthogonality is exact and holds by construction.")
    else:
        print(f"\nSome claims failed — investigate before reporting as a theorem.")


if __name__ == "__main__":
    main()
