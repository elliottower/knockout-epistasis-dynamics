"""Batch 2 blind experiment: sweep remaining BUILTIN models + web-sourced models.

Runs coalition sweeps, scores composition, and compiles blind results.
All predictions are frozen to file BEFORE any dynamics are run.
"""
import hashlib
import json
import multiprocessing as mp
import sys
from collections import defaultdict
from datetime import datetime, timezone
from itertools import combinations
from pathlib import Path

import numpy as np
from tqdm import tqdm

sys.path.insert(0, str(Path(__file__).parent.parent))

from grn_coalition_sweep import (
    BUILTIN_MODELS,
    _find_regulators,
    compile_network,
    extract_interaction_graph,
    extract_rule_fourier,
    identify_input_nodes,
    simulate_sync_output,
)
from composition_scorer import score_composition
from scripts.structural_analysis import (
    analyze_model,
    find_feedback_loops,
    max_and_arity,
    compute_pathway_depth,
)


RESULTS_DIR = Path(__file__).parent.parent / "results" / "grn_v2"


def timestamp():
    return datetime.now(timezone.utc).isoformat()


# ── Web-sourced models ───────────────────────────────────────────────
# These are published Boolean GRN models sourced from the literature.
# Each has n <= 18 nodes, synchronous update, and a biologically
# meaningful output node.

WEB_MODELS = {
    "li_budding_yeast": {
        "citation": "Li et al., PNAS 101.14 (2004): 4781-4786",
        "description": "Budding yeast cell cycle (11 nodes, 1 input)",
        "output_nodes": ["Clb1_2"],
        "rules": {
            "Cln3": "Cln3",
            "MBF": "!Clb1_2 & (Cln3 | MBF & !Clb1_2 & !Cln3) & !MBF",
            "SBF": "!Clb1_2 & (Cln3 | SBF & !Clb1_2 & !Cln3) & !MBF",
            "Cln1_2": "SBF",
            "Cdh1": "!Cln1_2 & !Clb1_2 & !Clb5_6 | Cdc20",
            "Swi5": "!Clb1_2 & (Cdc20 | Mcm1)",
            "Cdc20": "Clb1_2 & Mcm1",
            "Clb5_6": "MBF & !(Cdc20 | Cdh1)",
            "Sic1": "!Cln1_2 & !Clb1_2 & !Clb5_6 | Cdc20 | Swi5",
            "Clb1_2": "!(Sic1 | Cdc20 | Cdh1) & (Clb5_6 | Clb1_2 | Mcm1)",
            "Mcm1": "Clb1_2",
        },
    },
    "mendoza_thelper": {
        "citation": "Mendoza & Xenarios, Bioinformatics 22.11 (2006): 1375-1382",
        "description": "T-helper cell Th1/Th2 differentiation (13 nodes, 3 inputs)",
        "output_nodes": ["GATA3"],
        "rules": {
            "APC": "APC",
            "IFNg_e": "IFNg_e",
            "IL12_e": "IL12_e",
            "IFNg": "T_bet & !GATA3 | IFNg_e & !GATA3 | STAT1 & !GATA3",
            "IFNgR": "IFNg & !SOCS1",
            "IL12R": "IL12_e & !GATA3 | T_bet & !GATA3",
            "IL4": "GATA3 & !T_bet",
            "IL4R": "IL4 & !SOCS1",
            "STAT1": "IFNgR",
            "STAT6": "IL4R",
            "SOCS1": "STAT1 | T_bet",
            "T_bet": "STAT1 & !GATA3 | T_bet & !GATA3",
            "GATA3": "STAT6 & !T_bet | GATA3 & !T_bet",
        },
    },
    "albert_segment_polarity": {
        "citation": "Albert & Othmer, J Theor Biol 223.1 (2003): 1-18",
        "description": "Drosophila segment polarity single-cell (11 nodes, 2 inputs)",
        "output_nodes": ["wg"],
        "rules": {
            "SLP": "SLP",
            "WG": "WG",
            "wg": "CIA & !CIR & SLP | wg & CIA & !CIR",
            "WG_ext": "WG",
            "en": "WG_ext & !SLP",
            "EN": "en",
            "hh": "EN & !CIR",
            "HH_ext": "hh",
            "PTC": "!EN & !HH_ext | PTC & !EN & !HH_ext",
            "CIA": "!PTC | HH_ext",
            "CIR": "PTC & !HH_ext",
        },
    },
    "calzone_cell_fate": {
        "citation": "Calzone et al., PLoS Comput Biol 6.9 (2010): e1000912",
        "description": "Cell fate decision: survival vs apoptosis (12 nodes, 2 inputs)",
        "output_nodes": ["Apoptosis"],
        "rules": {
            "TNF": "TNF",
            "FADD": "FADD",
            "DISC_TNF": "TNF & FADD & !FLIP",
            "DISC_FAS": "TNF & FADD & !FLIP",
            "Caspase8": "DISC_TNF | DISC_FAS & !cIAP | DISC_FAS & !BAR",
            "cIAP": "NFkB & !SMAC",
            "BAR": "BAR",
            "RIP1": "DISC_TNF & !cIAP",
            "NFkB": "RIP1 | IKK",
            "IKK": "RIP1 & !Caspase3",
            "FLIP": "NFkB & !DISC_TNF",
            "Caspase3": "Caspase8 & !cIAP | Caspase8 & SMAC",
            "SMAC": "MOMP",
            "MOMP": "Caspase8 & !BCL2 | BAX & !BCL2",
            "BAX": "Caspase8 & !BCL2",
            "BCL2": "NFkB & !BAX & !Caspase3",
            "Apoptosis": "Caspase3",
        },
    },
    "grieco_bladder": {
        "citation": "Grieco et al., Mol BioSyst 9.7 (2013): 1646-1659",
        "description": "Bladder cancer signaling: growth vs apoptosis (18 nodes, 3 inputs)",
        "output_nodes": ["Apoptosis"],
        "rules": {
            "GF": "GF",
            "EGFR_stimulus": "EGFR_stimulus",
            "FGFR3_stimulus": "FGFR3_stimulus",
            "EGFR": "EGFR_stimulus & !FGFR3 & !GRB2",
            "FGFR3": "FGFR3_stimulus & !GRB2",
            "GRB2": "EGFR | FGFR3",
            "RAS": "GRB2 | GF",
            "PI3K": "GRB2 | RAS",
            "AKT": "PI3K",
            "PTEN": "TP53",
            "RB1": "!CyclinD1 & !p16",
            "p16": "!RB1",
            "CyclinD1": "RAS & !p21",
            "p21": "TP53 & !AKT",
            "TP53": "!MDM2",
            "MDM2": "AKT & !TP53 | TP53 & !AKT",
            "Growth_arrest": "RB1 | p21",
            "Apoptosis": "TP53 & !AKT",
        },
    },
    "irons_cardiac": {
        "citation": "Irons & Monk, J Theor Biol 259.4 (2009): 760-769",
        "description": "Cardiac gene regulatory network (15 nodes, no inputs)",
        "output_nodes": ["Nkx2_5"],
        "rules": {
            "BMP2": "!Noggin & Smad",
            "Dkk1": "BMP2 & !Wnt",
            "FGF8": "!BMP2 & Nkx2_5 | Wnt & Nkx2_5",
            "Gata4": "BMP2 & Nkx2_5 | Gata4 & Nkx2_5",
            "HAS2": "BMP2 & !FGF8",
            "Hex": "BMP2 | Wnt",
            "Mesp1": "Wnt & !Hex | BMP2 & !Hex",
            "Nkx2_5": "Gata4 | BMP2 & Smad | Mesp1",
            "Noggin": "BMP2 & Smad",
            "Smad": "BMP2 & !Noggin",
            "Snail": "BMP2 & !FGF8 | Wnt & !FGF8",
            "Tbx5": "Nkx2_5 & Gata4",
            "Wnt": "!Dkk1 & Mesp1 | !Dkk1 & Wnt",
            "bCatenin": "Wnt",
            "eHand": "Nkx2_5 & !Hex & Gata4",
        },
    },
    "zanudo_tlgl": {
        "citation": "Zanudo & Albert, PLoS Comput Biol 11.4 (2015): e1004193",
        "description": "T-LGL leukemia survival reduced (12 nodes, 1 input)",
        "output_nodes": ["Apoptosis"],
        "rules": {
            "Stimuli": "Stimuli",
            "PDGFR": "Stimuli | S1P",
            "S1P": "PDGFR & !Ceramide",
            "Ceramide": "Fas & !S1P",
            "Fas": "!S1P & Ceramide | !sFas",
            "sFas": "NFkB",
            "DISC": "Fas & !FLIP",
            "FLIP": "NFkB",
            "Caspase": "DISC & !IAP",
            "IAP": "NFkB",
            "NFkB": "Stimuli",
            "Apoptosis": "Caspase",
        },
    },
    "remy_p53_mdm2": {
        "citation": "Abou-Jaoude et al., BMC Syst Biol 3 (2009): 100",
        "description": "p53-Mdm2 DNA damage response network (10 nodes, 1 input)",
        "output_nodes": ["p53"],
        "rules": {
            "DNAdam": "DNAdam",
            "ATM": "DNAdam & !p53",
            "p53": "!Mdm2cyt & (ATM | p53 | p21)",
            "Mdm2cyt": "p53 & !ATM",
            "Mdm2nuc": "Mdm2cyt & !ATM",
            "p21": "p53 & !Mdm2nuc",
            "CycE": "!p21 & !Rb",
            "Rb": "!CycE & !CycE",
            "E2F1": "!Rb & !CycE",
            "Caspase3": "p53 & !Mdm2nuc & CycE",
        },
    },
    "fumia_cellcycle": {
        "citation": "Fumia & Martins, PLoS ONE 8.7 (2013): e69008",
        "description": "Mammalian cell cycle regulation simplified (14 nodes, 1 input)",
        "output_nodes": ["CycB"],
        "rules": {
            "GF": "GF",
            "CycD": "GF & !p21 & !p27",
            "CycE": "E2F & !p21 & !p27 & !Rb",
            "CycA": "E2F & !Cdc20 & !Rb & !cdh1 & !p21",
            "CycB": "!Cdc20 & !cdh1 & !p21 & CycA",
            "Rb": "!CycD & !CycE & !CycA & !CycB | p27 & !CycD & !CycB",
            "E2F": "!Rb & !CycB & !CycA | p27 & !Rb & !CycB",
            "p21": "!CycE & !CycA & !CycB & p53",
            "p27": "!CycE & !CycA & !CycB & !CycD",
            "p53": "!Mdm2 | !AKT",
            "AKT": "GF",
            "Mdm2": "p53 & AKT | p53 & !p21",
            "Cdc20": "CycB",
            "cdh1": "!CycA & !CycB | Cdc20 | p21",
        },
    },
}


# ── Structural analysis helper ──────────────────────────────────────

def analyze_model_generic(model_name, rules, output_nodes, description=""):
    """Structural analysis for any model (builtin or web-sourced)."""
    node_names = list(rules.keys())
    n = len(node_names)

    inputs = identify_input_nodes(rules)
    graph = extract_interaction_graph(rules)
    n_edges = sum(len(edges) for edges in graph.values())
    rf = extract_rule_fourier(rules)
    loops = find_feedback_loops(rules, max_length=4)
    pos_loops = [l for l in loops if l[1] > 0]
    neg_loops = [l for l in loops if l[1] < 0]
    mixed_loops = [l for l in loops if l[1] == 0]
    arity = max_and_arity(rules)
    depth = compute_pathway_depth(rules, output_nodes)

    output_rule = rules[output_nodes[0]]
    output_regs = _find_regulators(output_rule, node_names)

    reg_counts = {}
    for node in node_names:
        regs = _find_regulators(rules[node], node_names)
        reg_counts[node] = len(regs)

    return {
        "model": model_name,
        "description": description,
        "n_nodes": n,
        "n_coalitions": 2**n,
        "output_nodes": output_nodes,
        "output_rule": output_rule,
        "output_regulators": output_regs,
        "n_output_regulators": len(output_regs),
        "input_nodes": inputs,
        "n_edges": n_edges,
        "edge_density": n_edges / n,
        "max_and_arity": arity,
        "pathway_depth": depth,
        "n_positive_loops": len(pos_loops),
        "n_negative_loops": len(neg_loops),
        "n_mixed_loops": len(mixed_loops),
        "local_energy_spectrum": rf["local_energy_spectrum"],
        "n_local_pairwise": rf["n_pairwise"],
        "n_local_triples": rf["n_triples"],
        "per_node_regulator_count": reg_counts,
    }


# ── Blind prediction engine ─────────────────────────────────────────

def make_blind_prediction(analysis):
    """Generate blind prediction from structural analysis alone.

    Heuristics:
    - Positive feedback dominated → bistable → likely DESTROYS higher-order
    - Balanced or negative feedback → complex attractors → likely CREATES
    - High edge density + deep pathways → more composition → more creation
    - High local o3+ with low feedback → local structure preserved (high rho)
    - Output with many regulators → more sensitivity to knockouts
    """
    n = analysis["n_nodes"]
    n_pos = analysis["n_positive_loops"]
    n_neg = analysis["n_negative_loops"]
    n_loops = n_pos + n_neg + analysis["n_mixed_loops"]
    density = analysis["edge_density"]
    depth = analysis["pathway_depth"]
    arity = analysis["max_and_arity"]
    n_inputs = len(analysis["input_nodes"])

    spec = analysis["local_energy_spectrum"]
    local_o3plus = sum(spec[3:]) if len(spec) > 3 else 0.0
    local_o2 = spec[2] if len(spec) > 2 else 0.0

    # Feedback balance ratio: positive/(positive+negative)
    if n_pos + n_neg > 0:
        pos_ratio = n_pos / (n_pos + n_neg)
    else:
        pos_ratio = 0.5  # no loops = neutral

    # Prediction 1: creation or destruction
    # Positive-feedback-dominated → bistable → destroys
    # Balanced/negative → complex → creates
    # Deep pathways + high density favor creation
    creation_score = (
        (1 - pos_ratio) * 0.4  # negative loops favor creation
        + min(depth / 10, 1.0) * 0.3  # deep paths favor creation
        + min(density / 4, 1.0) * 0.2  # high density favors creation
        + (1 if n_neg > n_pos else 0) * 0.1  # explicit negative dominance
    )
    predicts_creation = creation_score > 0.45

    # Prediction 2: Spearman rho range
    # High local o3+ → more structure to correlate with → potentially higher rho
    # But high feedback → rearranges which pairs matter → lower rho
    # Input nodes reduce player count, often leaving cleaner structure
    base_rho = 0.4
    rho_adj = 0.0
    rho_adj += min(local_o3plus * 5, 0.15)  # more local structure → higher rho
    rho_adj -= min(n_loops / 30, 0.2)  # many loops → lower rho
    rho_adj += n_inputs * 0.03  # inputs simplify dynamics
    rho_adj -= min(depth / 15, 0.15)  # deep paths → lower rho
    rho_center = max(-0.1, min(0.9, base_rho + rho_adj))
    rho_range = [round(max(-1.0, rho_center - 0.2), 2),
                 round(min(1.0, rho_center + 0.2), 2)]

    # Prediction 3: global o3+ energy fraction
    # Start from local o3+, adjust for dynamics
    if predicts_creation:
        global_o3_est = local_o3plus * 1.5 + 0.02
    else:
        global_o3_est = local_o3plus * 0.6
    global_o3_est = round(max(0.001, min(0.3, global_o3_est)), 4)

    return {
        "predicts_creation": predicts_creation,
        "creation_or_destruction": "creation" if predicts_creation else "destruction",
        "spearman_rho_range": rho_range,
        "global_o3plus_estimate": global_o3_est,
        "reasoning": {
            "pos_ratio": round(pos_ratio, 3),
            "creation_score": round(creation_score, 3),
            "local_o3plus": round(local_o3plus, 4),
            "n_loops": n_loops,
            "depth": depth,
            "density": round(density, 2),
        },
    }


# ── Parallel coalition sweep ────────────────────────────────────────

def _run_one_coalition(args):
    (coalition, player_indices, init_states, compiled,
     clamp_value, output_indices, max_steps, n_total) = args
    clamp_mask = np.zeros(n_total, dtype=bool)
    for bit_pos, node_idx in enumerate(player_indices):
        if not (coalition & (1 << bit_pos)):
            clamp_mask[node_idx] = True
    output, info = simulate_sync_output(
        init_states, compiled, clamp_mask, clamp_value,
        output_indices, max_steps=max_steps,
    )
    return coalition, output, info["n_fixed_point"], info["n_cycling"]


def sweep_coalitions_parallel(rules, output_nodes, n_init=512, max_steps=200,
                               seed=42, clamp_value=0, n_workers=10):
    compiled, node_names = compile_network(rules)
    n_total = len(node_names)
    name_to_idx = {name: i for i, name in enumerate(node_names)}
    output_indices = [name_to_idx[o] for o in output_nodes]

    player_indices = list(range(n_total))
    player_names = node_names
    n_players = len(player_indices)
    N = 2**n_players

    rng = np.random.default_rng(seed)
    init_states = rng.integers(0, 2, size=(n_init, n_total), dtype=np.int8)

    values = np.zeros((N, n_init), dtype=np.float64)
    total_cycling = 0
    total_fixed = 0

    task_args = [
        (c, player_indices, init_states, compiled, clamp_value,
         output_indices, max_steps, n_total)
        for c in range(N)
    ]

    with mp.Pool(n_workers) as pool:
        for coalition, output, n_fixed, n_cyc in tqdm(
            pool.imap_unordered(_run_one_coalition, task_args, chunksize=64),
            total=N,
            desc=f"Coalitions ({n_players}n, {n_workers}w)",
        ):
            values[coalition] = output
            total_fixed += n_fixed
            total_cycling += n_cyc

    convergence = {
        "update_scheme": "sync",
        "clamp_value": clamp_value,
        "total_simulations": N * n_init,
        "total_fixed_point": total_fixed,
        "total_cycling": total_cycling,
        "cycling_fraction": total_cycling / (N * n_init),
    }

    return {
        "values": values,
        "node_names": node_names,
        "player_names": player_names,
        "output_nodes": output_nodes,
        "n_players": n_players,
        "n_total_nodes": n_total,
        "n_init_states": n_init,
        "convergence": convergence,
    }


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


# ── Main ─────────────────────────────────────────────────────────────

def main():
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    N_INIT = 512
    MAX_STEPS = 200
    SEED = 42
    N_WORKERS = 10
    CLAMP_VALUE = 0

    # Collect all models to sweep
    already_swept = {
        "faure_cellcycle", "tournier_apoptosis", "davidich_yeast",
        "drosophila_cellcycle", "fanconi_anemia", "arabidopsis_cellcycle",
    }

    all_models = {}
    # Remaining builtins
    for name, info in BUILTIN_MODELS.items():
        if name not in already_swept:
            all_models[name] = {
                "rules": info["rules"],
                "output_nodes": info["output_nodes"],
                "description": info["description"],
                "citation": info["citation"],
                "source": "builtin",
            }
    # Web-sourced models
    for name, info in WEB_MODELS.items():
        all_models[name] = {
            "rules": info["rules"],
            "output_nodes": info["output_nodes"],
            "description": info["description"],
            "citation": info["citation"],
            "source": "web",
        }

    print(f"[{timestamp()}] Batch 2 blind experiment")
    print(f"  Total models to sweep: {len(all_models)}")
    print(f"  Builtin remaining: {[k for k, v in all_models.items() if v['source'] == 'builtin']}")
    print(f"  Web-sourced: {[k for k, v in all_models.items() if v['source'] == 'web']}")

    # ── Phase 1: Structural analysis + blind predictions ──
    print(f"\n{'='*70}")
    print(f"PHASE 1: Structural analysis and blind predictions")
    print(f"{'='*70}")

    analyses = {}
    predictions = {}
    for model_name, model_info in all_models.items():
        rules = model_info["rules"]
        output_nodes = model_info["output_nodes"]
        n = len(rules)

        print(f"\n--- {model_name} (n={n}) ---")
        analysis = analyze_model_generic(
            model_name, rules, output_nodes, model_info["description"]
        )
        analyses[model_name] = analysis

        prediction = make_blind_prediction(analysis)
        predictions[model_name] = {
            "model": model_name,
            "source": model_info["source"],
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

    # Save predictions with SHA hash BEFORE any sweeps
    pred_path = RESULTS_DIR / "blind_predictions_batch2.json"
    pred_output = {
        "experiment": "batch2_blind_predictions",
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

    # ── Phase 2: Coalition sweeps ──
    print(f"\n{'='*70}")
    print(f"PHASE 2: Coalition sweeps")
    print(f"{'='*70}")

    # Sort by size (smallest first for quick wins)
    sorted_models = sorted(all_models.items(), key=lambda x: len(x[1]["rules"]))

    all_results = []
    for model_name, model_info in sorted_models:
        rules = model_info["rules"]
        output_nodes = model_info["output_nodes"]
        n = len(rules)

        print(f"\n{'='*60}")
        print(f"[{timestamp()}] === {model_name} (n={n}, 2^{n}={2**n} coalitions) ===")
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
            all_results.append({
                "model": model_name,
                "n": n,
                "error": str(e),
            })

    # ── Phase 3: Compile summary ──
    print(f"\n{'='*70}")
    print(f"PHASE 3: Compile summary")
    print(f"{'='*70}")

    summary_path = RESULTS_DIR / "blind_batch2_summary.json"
    summary = {
        "experiment": "batch2_blind_sweep",
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

    # ── Phase 4: Scorecard ──
    scorecard_lines = [
        "# Blind Batch 2 Scorecard",
        "",
        f"Prediction SHA-256: `{pred_sha}`",
        "",
        "## Predictions vs Actuals",
        "",
        "| Model | n | Pred Direction | Actual Direction | Pred rho Range | Actual rho | Pred o3+ | Actual o3+ | rho Hit? | Dir Hit? |",
        "|-------|---|---------------|-----------------|---------------|------------|----------|-----------|---------|---------|",
    ]

    n_dir_hits = 0
    n_rho_hits = 0
    n_scored = 0

    for result in all_results:
        model_name = result["model"]
        if "error" in result:
            scorecard_lines.append(
                f"| {model_name} | {result['n']} | - | ERROR | - | - | - | - | - | - |"
            )
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

        actual_o3 = result["global_o3plus"]
        pred_o3 = pred["global_o3plus_estimate"]

        scorecard_lines.append(
            f"| {model_name} | {result['n']} | {pred_dir} | {actual_dir} | "
            f"[{rho_lo:.2f}, {rho_hi:.2f}] | {actual_rho:.3f} | "
            f"{pred_o3:.4f} | {actual_o3:.4f} | "
            f"{'YES' if rho_hit else 'NO'} | {'YES' if dir_hit else 'NO'} |"
        )

    scorecard_lines.extend([
        "",
        "## Summary",
        "",
        f"- Models scored: {n_scored}",
        f"- Direction hits: {n_dir_hits}/{n_scored} ({n_dir_hits/max(1,n_scored):.0%})",
        f"- Rho in range: {n_rho_hits}/{n_scored} ({n_rho_hits/max(1,n_scored):.0%})",
    ])

    scorecard_path = RESULTS_DIR / "blind_batch2_scorecard.md"
    with open(scorecard_path, "w") as f:
        f.write("\n".join(scorecard_lines) + "\n")
    print(f"[{timestamp()}] Scorecard saved: {scorecard_path}")

    # Print summary table
    print(f"\n{'='*70}")
    print(f"RESULTS SUMMARY")
    print(f"{'='*70}")
    for r in sorted(all_results, key=lambda x: x.get("n", 0)):
        if "error" in r:
            print(f"  {r['model']:30s} n={r['n']:2d}  ERROR: {r['error']}")
        else:
            print(f"  {r['model']:30s} n={r['n']:2d}  rho={r['spearman_rho']:+.3f}  "
                  f"o3+={r['global_o3plus']:.4f}  {r['creation_or_destruction']}")

    print(f"\n[{timestamp()}] Batch 2 complete.")


if __name__ == "__main__":
    main()
