"""Batch 2b: Additional web-sourced models from PyBoolNet/Biodivine repositories.

These supplement the batch 2 sweep with smaller (fast) models found via
web search of curated Boolean model repositories.
"""
import hashlib
import json
import multiprocessing as mp
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from tqdm import tqdm

sys.path.insert(0, str(Path(__file__).parent.parent))

from grn_coalition_sweep import (
    _find_regulators,
    compile_network,
    extract_interaction_graph,
    extract_rule_fourier,
    identify_input_nodes,
    simulate_sync_output,
)
from composition_scorer import score_composition


RESULTS_DIR = Path(__file__).parent.parent / "results" / "grn_v2"


def timestamp():
    return datetime.now(timezone.utc).isoformat()


# ── Additional web-sourced models ───────────────────────────────────

EXTRA_MODELS = {
    "lambda_phage": {
        "citation": "Thieffry & Thomas, Bull Math Biol 57.2 (1995): 277-297",
        "description": "Lambda phage lysis/lysogeny decision (7 nodes, Booleanized multi-valued)",
        "output_nodes": ["CI_b1"],
        "rules": {
            "CII": "!CI_b1 & !Cro_b1 & N | !CI_b1 & Cro_b1 & !Cro_b2 & N | !CI_b1 & Cro_b1 & Cro_b2 & !Cro_b3 & N | CI_b1 & !CI_b2 & !Cro_b1 & N | CI_b1 & !CI_b2 & Cro_b1 & !Cro_b2 & N | CI_b1 & !CI_b2 & Cro_b1 & Cro_b2 & !Cro_b3 & N",
            "CI_b1": "!CI_b1 & !Cro_b1 | !CI_b1 & Cro_b1 & CII | CI_b1 & !CI_b2 & !Cro_b1 | CI_b1 & !CI_b2 & Cro_b1 & CII | CI_b1 & CI_b2",
            "CI_b2": "CI_b1 & !Cro_b1 | CI_b1 & Cro_b1 & CII",
            "Cro_b1": "!CI_b1 | CI_b1 & !CI_b2 | CI_b1 & CI_b2 & Cro_b1 & Cro_b2",
            "Cro_b2": "!CI_b1 & Cro_b1 | CI_b1 & !CI_b2 & Cro_b1 | CI_b1 & CI_b2 & Cro_b1 & Cro_b2 & Cro_b3",
            "Cro_b3": "!CI_b1 & Cro_b1 & Cro_b2 & !Cro_b3 | CI_b1 & !CI_b2 & Cro_b1 & Cro_b2 & !Cro_b3",
            "N": "!CI_b1 & !Cro_b1 | !CI_b1 & Cro_b1 & !Cro_b2",
        },
    },
    "cell_cycle_transcription": {
        "citation": "Orlando et al., Nature 453.7197 (2008): 944-947",
        "description": "Yeast cell cycle transcription oscillator (9 nodes, no inputs)",
        "output_nodes": ["SFF"],
        "rules": {
            "ACE2": "SFF",
            "CLN3": "SWI5 & ACE2 & !YOX1 & !YHP1",
            "HCM1": "MBF & SBF",
            "MBF": "CLN3",
            "SBF": "MBF & !YHP1 & !YOX1 | CLN3 & !YHP1 & !YOX1",
            "SFF": "SBF & HCM1",
            "SWI5": "SFF",
            "YHP1": "SBF | MBF",
            "YOX1": "MBF & SBF",
        },
    },
    "asymmetric_cell_division": {
        "citation": "Sanchez-Osorio et al., Comput Methods Syst Biol, Springer (2017): 3-19",
        "description": "Caulobacter crescentus asymmetric cell division (9 nodes, no inputs)",
        "output_nodes": ["CtrAb"],
        "rules": {
            "CckA": "DivL",
            "ChpT": "CckA",
            "ClpXP_RcdA": "!CpdR",
            "CpdR": "ChpT",
            "CtrAb": "ChpT & !ClpXP_RcdA",
            "DivJ": "DivK & !PleC",
            "DivK": "!PleC & DivJ",
            "DivL": "!DivK",
            "PleC": "!DivK",
        },
    },
    "arellano_rootstem": {
        "citation": "Azpeitia et al., BMC Syst Biol 4 (2010): 134",
        "description": "Arabidopsis root stem cell niche patterning (9 nodes, 2 inputs)",
        "output_nodes": ["WOX"],
        "rules": {
            "AUXINS": "AUXINS",
            "SHR": "SHR",
            "ARF": "!IAA",
            "IAA": "!AUXINS",
            "JKD": "SHR & SCR",
            "MGP": "!WOX & SHR & SCR",
            "SCR": "SHR & SCR & !MGP | SHR & SCR & JKD",
            "WOX": "WOX & SHR & SCR & ARF | SHR & SCR & !MGP & ARF",
            "PLT": "ARF",
        },
    },
    "pair_rule_module": {
        "citation": "Sanchez & Thieffry, J Theor Biol 224.4 (2003): 517-537",
        "description": "Drosophila pair-rule gene module for segmentation (11 nodes, no inputs)",
        "output_nodes": ["Eve_b1"],
        "rules": {
            "Eve_b1": "!Eve_b1 & Prd_b1 & !Slp & !Odd | Eve_b1 & !Eve_b2 & !Run & !Slp | Eve_b1 & !Eve_b2 & Run & !Slp & !Odd | Eve_b1 & Eve_b2",
            "Eve_b2": "Eve_b1 & !Eve_b2 & Prd_b1 & !Run & !Slp & !Odd | Eve_b1 & Eve_b2 & !Eve_b3 & Prd_b1 & !Run & !Slp & !Odd | Eve_b1 & Eve_b2 & Eve_b3",
            "Eve_b3": "Eve_b1 & Eve_b2 & Prd_b1 & !Run & !Slp & !Odd",
            "Ftz_b1": "!Eve_b1 & !Slp & !Ftz_b1 & !Odd | !Eve_b1 & !Slp & Ftz_b1 & !Ftz_b2 & !Odd | !Eve_b1 & !Slp & Ftz_b1 & Ftz_b2 | !Eve_b1 & Slp & Ftz_b1 & Ftz_b2 | Eve_b1 & !Eve_b2 & !Slp & !Ftz_b1 & !Odd | Eve_b1 & !Eve_b2 & !Slp & Ftz_b1 & !Ftz_b2 & !Odd | Eve_b1 & !Eve_b2 & !Slp & Ftz_b1 & Ftz_b2 | Eve_b1 & !Eve_b2 & Slp & Ftz_b1 & Ftz_b2 | Eve_b1 & Eve_b2 & Ftz_b1 & Ftz_b2",
            "Ftz_b2": "!Eve_b1 & !Slp & Ftz_b1 & !Odd | Eve_b1 & !Eve_b2 & !Slp & Ftz_b1 & !Odd",
            "Odd": "!Eve_b1 & !Prd_b1",
            "Ppa": "!Eve_b1 | Eve_b1 & !Eve_b2",
            "Prd_b1": "!Prd_b1 & !Odd | Prd_b1 & !Prd_b2 & !Odd | Prd_b1 & Prd_b2",
            "Prd_b2": "Prd_b1 & !Ppa & !Odd",
            "Run": "!Eve_b1 & Prd_b1 & !Odd | Eve_b1 & !Eve_b2 & Prd_b1 & !Odd",
            "Slp": "!Eve_b1 & !Ftz_b1 & !Odd | Eve_b1 & !Eve_b2 & !Ftz_b1 & !Odd",
        },
    },
    "calzone_cellfate_reduced": {
        "citation": "Calzone et al., PLoS Comput Biol 6.3 (2010): e1000702",
        "description": "Reduced cell fate decision: apoptosis/necrosis/survival (11 nodes, 2 inputs)",
        "output_nodes": ["C3"],
        "rules": {
            "TNF": "TNF",
            "FAS": "FAS",
            "RIP1": "!C8 & TNF | !C8 & FAS",
            "NFkB": "cIAP & RIP1 & !C3",
            "C8": "TNF & !NFkB | FAS & !NFkB | C3 & !NFkB",
            "cIAP": "NFkB & !MOMP | cIAP & !MOMP",
            "ATP": "!MPT",
            "C3": "ATP & MOMP & !NFkB",
            "ROS": "!NFkB & RIP1 | !NFkB & MPT",
            "MOMP": "MPT | C8 & !NFkB",
            "MPT": "ROS & !NFkB",
        },
    },
    "saadatpour_guardcell": {
        "citation": "Li et al., PLoS Biol 4.10 (2006): e312; Saadatpour et al., SIAM 2013",
        "description": "Guard cell ABA signal transduction for stomatal closure (13 nodes, no inputs)",
        "output_nodes": ["KEV"],
        "rules": {
            "ADPRc": "NO",
            "CIS": "cGMP & cADPR | InsP3",
            "Ca2": "!Ca2ATP & CIS",
            "Ca2ATP": "Ca2",
            "GC": "NO",
            "InsP3": "PLC",
            "NO": "NOS",
            "NOS": "Ca2",
            "PLC": "Ca2",
            "cADPR": "ADPRc",
            "cGMP": "GC",
            "KAP": "!Ca2",
            "KEV": "Ca2",
        },
    },
    "morphogenetic_checkpoint": {
        "citation": "Faure et al., Mol Biosyst 5.12 (2009): 1787-1796",
        "description": "Yeast morphogenesis checkpoint (12 nodes, 1 input)",
        "output_nodes": ["Clb2_b2"],
        "rules": {
            "MASS_b1": "MASS_b1",
            "BUD": "MASS_b1",
            "Clb2_b1": "!Clb2_b1 & MASS_b1 | Clb2_b1 & !Clb2_b2 & MASS_b1 | Clb2_b1 & Clb2_b2",
            "Clb2_b2": "!Swe1_b1 & Clb2_b1 & MASS_b1 | Swe1_b1 & !Swe1_b2 & !Mih1_b1 & Clb2_b1 & MASS_b1 & MASS_b2 | Swe1_b1 & !Swe1_b2 & Mih1_b1 & !Mih1_b2 & Clb2_b1 & MASS_b1 & MASS_b2 | Swe1_b1 & !Swe1_b2 & Mih1_b1 & Mih1_b2 & Clb2_b1 & MASS_b1 | Swe1_b1 & Swe1_b2 & Mih1_b1 & !Mih1_b2 & Clb2_b1 & MASS_b1 & MASS_b2 | Swe1_b1 & Swe1_b2 & Mih1_b1 & Mih1_b2 & Clb2_b1 & MASS_b1",
            "Hsl1": "BUD",
            "MASS_b2": "MASS_b1 & MASS_b2",
            "Mih1_b1": "!Mih1_b1 & !Clb2_b1 & !Mpk1 | !Mih1_b1 & Clb2_b1 | Mih1_b1 & !Mih1_b2 & !Clb2_b1 & !Mpk1 | Mih1_b1 & !Mih1_b2 & Clb2_b1 | Mih1_b1 & Mih1_b2",
            "Mih1_b2": "Mih1_b1 & Clb2_b1 & !Mpk1",
            "Mpk1": "!BUD",
            "SBF": "!Clb2_b1 & MASS_b1 | Clb2_b1 & !Clb2_b2 & MASS_b1",
            "Swe1_b1": "!SBF & Swe1_b1 & Swe1_b2 | SBF & !Swe1_b1 & !Clb2_b1 | SBF & !Swe1_b1 & Clb2_b1 & !Clb2_b2 | SBF & !Swe1_b1 & Clb2_b1 & Clb2_b2 & !Hsl1 | SBF & Swe1_b1 & !Swe1_b2 & !Clb2_b1 | SBF & Swe1_b1 & !Swe1_b2 & Clb2_b1 & !Clb2_b2 | SBF & Swe1_b1 & !Swe1_b2 & Clb2_b1 & Clb2_b2 & !Hsl1 | SBF & Swe1_b1 & Swe1_b2",
            "Swe1_b2": "SBF & Swe1_b1 & !Clb2_b1 & !Hsl1 | SBF & Swe1_b1 & Clb2_b1 & !Clb2_b2 & !Hsl1",
        },
    },
    "lac_operon": {
        "citation": "Veliz-Cuba & Stigler, J Comput Biol 18.6 (2011): 783-794",
        "description": "Lac operon bistability: lactose utilization (13 nodes, 3 inputs)",
        "output_nodes": ["P"],
        "rules": {
            "Ge": "Ge",
            "Le": "Le",
            "Lem": "Lem",
            "A": "L & B",
            "Am": "Lm | L",
            "B": "M",
            "C": "!Ge",
            "L": "Le & P & !Ge",
            "Lm": "Lem & P & !Ge | Le & !Ge",
            "M": "C & !R & !Rm",
            "P": "M",
            "R": "!Am & !A",
            "Rm": "R | !R & !Am & !A",
        },
    },
    "hematopoiesis_aging": {
        "citation": "Herault et al., Comput Struct Biotechnol J 21 (2023): 21-33",
        "description": "Early hematopoiesis aging network (15 nodes, no inputs)",
        "output_nodes": ["Gata1"],
        "rules": {
            "Bclaf1": "Myc",
            "CDK46CycD": "Bclaf1 | Myc",
            "CIPKIP": "Junb",
            "Cebpa": "Gata2 & !Ikzf1 | Spi1 & !Ikzf1",
            "Egr1": "Gata2 & Junb",
            "Fli1": "Junb | Gata1 & !Klf1",
            "Gata1": "Fli1 | Gata2 & !Spi1 | Gata1 & !Ikzf1 & !Spi1",
            "Gata2": "Gata2 & !Gata1 & !Zfpm1 | Egr1 & !Gata1 & !Zfpm1 & !Spi1",
            "Ikzf1": "Gata2",
            "Junb": "Egr1 | Myc",
            "Klf1": "Gata1 & !Fli1",
            "Myc": "Cebpa & Bclaf1",
            "Spi1": "Spi1 & !Gata1 | Cebpa & !Gata1 & !Gata2",
            "Tal1": "Gata1 & !Spi1",
            "Zfpm1": "Gata1",
        },
    },
}


# ── Reuse analysis/sweep/scoring from batch 2 ──────────────────────

from scripts.run_batch2_blind_sweep import (
    analyze_model_generic,
    make_blind_prediction,
    sweep_coalitions_parallel,
)


def run_model(model_name, rules, output_nodes, n_init=512, max_steps=200,
              seed=42, n_workers=10, clamp_value=0):
    """Run coalition sweep + composition scoring for one model."""
    rf = extract_rule_fourier(rules)
    wiring = extract_interaction_graph(rules)
    node_names = list(rules.keys())
    n = len(node_names)

    print(f"  [{timestamp()}] Local Fourier: {rf['n_pairwise']} pairwise, {rf['n_triples']} triples")
    print(f"  [{timestamp()}] Starting coalition sweep (2^{n} = {2**n} coalitions)...")

    result = sweep_coalitions_parallel(
        rules, output_nodes, n_init=n_init, max_steps=max_steps,
        seed=seed, clamp_value=clamp_value, n_workers=n_workers,
    )

    v_mean = result["values"].mean(axis=1)
    n_unique = len(np.unique(np.round(v_mean, 6)))
    total_var = float(np.var(v_mean))
    conv = result["convergence"]

    print(f"  Value function: unique={n_unique}, var={total_var:.4e}")
    print(f"  Fixed: {conv['total_fixed_point']}, Cycling: {conv['total_cycling']} ({conv['cycling_fraction']:.1%})")

    # Save coalition NPZ
    npz_path = RESULTS_DIR / f"{model_name}_coalition_blind.npz"
    np.savez_compressed(
        npz_path,
        target_logits=result["values"],
        foil_logits=np.zeros_like(result["values"]),
        coalition_indices=np.arange(2**n, dtype=np.int64),
        circuit_heads=np.array(node_names, dtype=object),
        n_players=np.int64(n),
        n_prompts=np.int64(n_init),
        model_name=model_name,
        output_nodes=np.array(output_nodes, dtype=object),
    )
    print(f"  Saved: {npz_path}")

    # Save wiring JSON
    wiring_path = RESULTS_DIR / f"{model_name}_wiring_blind.json"
    wiring_output = {
        "model": model_name,
        "n_nodes": n,
        "node_names": node_names,
        "output_nodes": output_nodes,
        "interaction_graph": wiring,
        "rule_fourier": rf,
        "convergence": conv,
        "n_unique_values": n_unique,
        "total_variance": total_var,
        "timestamp": timestamp(),
    }
    with open(wiring_path, "w") as f:
        json.dump(wiring_output, f, indent=2)
    print(f"  Saved: {wiring_path}")

    # Score composition
    print(f"  [{timestamp()}] Scoring composition...")
    scores = score_composition(v_mean, n, node_names, rf)
    pw = scores["pairwise"]
    spec = scores["energy_spectrum"]
    g3p = sum(spec["global"][3:])
    l3p = sum(spec["local_rules"][3:])
    tri = scores.get("triples", {})

    print(f"  Spearman rho={pw['spearman_rho']:.4f} (p={pw['spearman_pvalue']:.2e})")
    print(f"  CI: [{pw['spearman_ci_95'][0]:.4f}, {pw['spearman_ci_95'][1]:.4f}]")
    print(f"  Global o3+={g3p:.1%}, Local o3+={l3p:.1%}, Delta={g3p-l3p:+.1%}")

    # Save composition JSON
    comp_path = RESULTS_DIR / f"{model_name}_composition_blind.json"
    comp_output = {
        "model": model_name,
        "n_players": n,
        "n_init": n_init,
        "node_names": node_names,
        "timestamp": timestamp(),
        **scores,
    }
    with open(comp_path, "w") as f:
        json.dump(comp_output, f, indent=2)
    print(f"  Saved: {comp_path}")

    tri_rho = tri.get("spearman_rho", None)
    return {
        "model": model_name,
        "n": n,
        "spearman_rho": pw["spearman_rho"],
        "spearman_ci": pw["spearman_ci_95"],
        "spearman_p": pw["spearman_pvalue"],
        "global_o3plus": g3p,
        "local_o3plus": l3p,
        "delta_o3plus": g3p - l3p,
        "creation_or_destruction": "creation" if g3p > l3p else "destruction",
        "cycling_fraction": conv["cycling_fraction"],
        "triple_rho": tri_rho,
        "global_spectrum": spec["global"],
        "local_spectrum": spec["local_rules"],
    }


def main():
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    N_INIT = 512
    MAX_STEPS = 200
    SEED = 42
    N_WORKERS = 10
    CLAMP_VALUE = 0

    print(f"[{timestamp()}] Batch 2b: {len(EXTRA_MODELS)} additional models")

    # ── Phase 1: Structural analysis + blind predictions ──
    print(f"\n{'='*70}")
    print("PHASE 1: Structural analysis and blind predictions")
    print(f"{'='*70}")

    predictions = {}
    for model_name, model_info in EXTRA_MODELS.items():
        rules = model_info["rules"]
        output_nodes = model_info["output_nodes"]
        n = len(rules)

        print(f"\n--- {model_name} (n={n}) ---")
        analysis = analyze_model_generic(
            model_name, rules, output_nodes, model_info["description"]
        )
        prediction = make_blind_prediction(analysis)
        predictions[model_name] = {
            "model": model_name,
            "source": "web_extra",
            "citation": model_info["citation"],
            "n": n,
            "structural": {
                "n_edges": analysis["n_edges"],
                "edge_density": analysis["edge_density"],
                "max_and_arity": analysis["max_and_arity"],
                "pathway_depth": analysis["pathway_depth"],
                "n_positive_loops": analysis["n_positive_loops"],
                "n_negative_loops": analysis["n_negative_loops"],
                "n_mixed_loops": analysis["n_mixed_loops"],
                "n_input_nodes": len(analysis["input_nodes"]),
                "input_nodes": analysis["input_nodes"],
                "output_regulators": analysis["output_regulators"],
                "n_output_regulators": analysis["n_output_regulators"],
                "local_energy_spectrum": analysis["local_energy_spectrum"],
                "local_o3plus": sum(analysis["local_energy_spectrum"][3:]) if len(analysis["local_energy_spectrum"]) > 3 else 0.0,
                "n_local_pairwise": analysis["n_local_pairwise"],
                "n_local_triples": analysis["n_local_triples"],
            },
            "prediction": prediction,
        }

        p = prediction
        print(f"  Prediction: {p['creation_or_destruction']}, "
              f"rho=[{p['spearman_rho_range'][0]:.2f}, {p['spearman_rho_range'][1]:.2f}], "
              f"global_o3+={p['global_o3plus_estimate']:.4f}")

    # Save predictions with SHA hash
    pred_path = RESULTS_DIR / "blind_predictions_batch2b.json"
    pred_output = {
        "experiment": "batch2b_blind_predictions",
        "timestamp": timestamp(),
        "n_models": len(predictions),
        "models": predictions,
    }
    pred_json = json.dumps(pred_output, indent=2, sort_keys=True)
    pred_sha = hashlib.sha256(pred_json.encode()).hexdigest()
    pred_output["sha256"] = pred_sha
    with open(pred_path, "w") as f:
        json.dump(pred_output, f, indent=2, sort_keys=True)
    print(f"\n[{timestamp()}] Predictions frozen: {pred_path}")
    print(f"  SHA-256: {pred_sha}")

    # ── Phase 2: Coalition sweeps (sorted by size, smallest first) ──
    print(f"\n{'='*70}")
    print("PHASE 2: Coalition sweeps")
    print(f"{'='*70}")

    sorted_models = sorted(EXTRA_MODELS.items(), key=lambda x: len(x[1]["rules"]))

    all_results = []
    for model_name, model_info in sorted_models:
        rules = model_info["rules"]
        output_nodes = model_info["output_nodes"]
        n = len(rules)

        print(f"\n{'='*60}")
        print(f"[{timestamp()}] === {model_name} (n={n}, 2^{n}={2**n}) ===")
        print(f"{'='*60}")

        try:
            result = run_model(
                model_name, rules, output_nodes,
                n_init=N_INIT, max_steps=MAX_STEPS, seed=SEED,
                n_workers=N_WORKERS, clamp_value=CLAMP_VALUE,
            )
            all_results.append(result)
        except Exception as e:
            print(f"  ERROR: {e}")
            import traceback
            traceback.print_exc()
            all_results.append({"model": model_name, "n": n, "error": str(e)})

    # ── Phase 3: Summary ──
    summary_path = RESULTS_DIR / "blind_batch2b_summary.json"
    summary = {
        "experiment": "batch2b_extra_models",
        "timestamp": timestamp(),
        "prediction_sha256": pred_sha,
        "n_models": len(all_results),
        "parameters": {
            "n_init": N_INIT,
            "max_steps": MAX_STEPS,
            "seed": SEED,
            "n_workers": N_WORKERS,
            "clamp_value": CLAMP_VALUE,
        },
        "results": all_results,
    }
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"\n[{timestamp()}] Summary saved: {summary_path}")

    # Print results table
    print(f"\n{'='*70}")
    print("RESULTS")
    print(f"{'='*70}")
    for r in sorted(all_results, key=lambda x: x.get("n", 0)):
        if "error" in r:
            print(f"  {r['model']:30s} n={r['n']:2d}  ERROR")
        else:
            print(f"  {r['model']:30s} n={r['n']:2d}  rho={r['spearman_rho']:+.3f}  "
                  f"o3+={r['global_o3plus']:.4f}  {r['creation_or_destruction']}")

    # ── Phase 4: Scorecard ──
    scorecard_lines = [
        "# Blind Batch 2b Scorecard (Extra Models)",
        "",
        f"Prediction SHA-256: `{pred_sha}`",
        "",
        "| Model | n | Pred Dir | Actual Dir | Pred rho | Actual rho | Pred o3+ | Actual o3+ | rho OK | Dir OK |",
        "|-------|---|---------|-----------|---------|-----------|---------|-----------|--------|--------|",
    ]

    n_dir_hits = 0
    n_rho_hits = 0
    n_scored = 0

    for result in all_results:
        model_name = result["model"]
        if "error" in result:
            continue

        n_scored += 1
        pred = predictions[model_name]["prediction"]
        actual_dir = result["creation_or_destruction"]
        pred_dir = pred["creation_or_destruction"]
        dir_hit = actual_dir == pred_dir
        if dir_hit:
            n_dir_hits += 1

        actual_rho = result["spearman_rho"]
        rho_lo, rho_hi = pred["spearman_rho_range"]
        rho_hit = rho_lo <= actual_rho <= rho_hi
        if rho_hit:
            n_rho_hits += 1

        scorecard_lines.append(
            f"| {model_name} | {result['n']} | {pred_dir} | {actual_dir} | "
            f"[{rho_lo:.2f},{rho_hi:.2f}] | {actual_rho:.3f} | "
            f"{pred['global_o3plus_estimate']:.4f} | {result['global_o3plus']:.4f} | "
            f"{'YES' if rho_hit else 'NO'} | {'YES' if dir_hit else 'NO'} |"
        )

    scorecard_lines.extend([
        "",
        f"Direction: {n_dir_hits}/{n_scored}, Rho range: {n_rho_hits}/{n_scored}",
    ])

    sc_path = RESULTS_DIR / "blind_batch2b_scorecard.md"
    with open(sc_path, "w") as f:
        f.write("\n".join(scorecard_lines) + "\n")
    print(f"\n[{timestamp()}] Scorecard: {sc_path}")
    print(f"  Direction: {n_dir_hits}/{n_scored} ({n_dir_hits/max(1,n_scored):.0%})")
    print(f"  Rho range: {n_rho_hits}/{n_scored} ({n_rho_hits/max(1,n_scored):.0%})")


if __name__ == "__main__":
    main()
