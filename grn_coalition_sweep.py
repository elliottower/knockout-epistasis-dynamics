"""Generate exhaustive coalition tables from Boolean gene regulatory networks.

GRN models provide two objects to compare:
1. Attractor-level interaction structure (Walsh decomposition of the
   coalition value function, measured via knockout sweep)
2. Local rule interaction structure (Walsh decomposition of each
   node's truth table, computed exactly — representation-independent)

The scientific question is: how much of the local rule Fourier
structure survives dynamical composition into attractor-level epistasis?
This is NOT validation (local != global), but it IS a question the GRN
regime is uniquely positioned to answer.

The wiring diagram (edge list) is a lossy projection of the rules,
reported as a secondary comparison.

Knockout = ablation. A gene knockout fixes a node's update rule to a
constant (0 for knockout, 1 for constitutive activation), just like
zero/mean ablation fixes a head's output. Coalition S = set of nodes
whose update rules are intact (not clamped).

Value function: for each coalition x initial state, run Boolean
dynamics until attractor convergence. If the trajectory enters a
limit cycle, average the output over one full cycle period (not an
arbitrary phase). Averaging over initial states gives a continuous
value in [0,1].

Supports both synchronous and asynchronous update schemes. Models
published as asynchronous (tournier_apoptosis)
should be run under both schemes; disagreement is flagged.

Usage:
    uv run python grn_coalition_sweep.py \
        --model faure_cellcycle --output results/faure_cellcycle_coalition.npz

    uv run python grn_coalition_sweep.py \
        --model tournier_apoptosis --update-scheme async \
        --output results/tournier_apoptosis_async_coalition.npz

    uv run python grn_coalition_sweep.py \
        --model faure_cellcycle --clamp-value 1 \
        --output results/faure_cellcycle_clamp1_coalition.npz
"""

import argparse
import json
import re
from datetime import datetime, timezone
from itertools import combinations
from pathlib import Path

import numpy as np
from tqdm import tqdm

from data_utils import wht


def timestamp():
    return datetime.now(timezone.utc).isoformat()


# ── Built-in models (embedded .bnet rules) ──────────────────────────

BUILTIN_MODELS = {
    "faure_cellcycle": {
        "citation": "Faure et al., Bioinformatics 22.14 (2006): e124-e131",
        "description": "Mammalian cell cycle control (10 nodes)",
        "output_nodes": ["CycB"],
        "rules": {
            "CycD": "CycD",
            "Cdc20": "CycB",
            "CycA": "!cdh1&!Rb&E2F&!Cdc20 | !UbcH10&!Rb&E2F&!Cdc20 | !cdh1&!Rb&CycA&!Cdc20 | !UbcH10&!Rb&CycA&!Cdc20",
            "CycB": "!cdh1&!Cdc20",
            "CycE": "!Rb&E2F",
            "E2F": "p27&!Rb&!CycB | !Rb&!CycB&!CycA",
            "Rb": "!CycE&!CycD&!CycB&!CycA | p27&!CycD&!CycB",
            "UbcH10": "UbcH10&CycB | UbcH10&CycA | UbcH10&Cdc20 | !cdh1",
            "cdh1": "p27&!CycB | !CycB&!CycA | Cdc20",
            "p27": "p27&!CycE&!CycD&!CycB | p27&!CycD&!CycB&!CycA | !CycE&!CycD&!CycB&!CycA",
        },
    },
    "tournier_apoptosis": {
        "citation": "Tournier & Chaves, J Theor Biol 260.2 (2009): 196-209",
        "description": "Apoptosis signaling network (12 nodes)",
        "output_nodes": ["C3a"],
        "rules": {
            "TNF": "TNF",
            "A20a": "TNF&NFkBnuc",
            "C3a": "!IAP&C8a",
            "C8a": "T2&!CARP | !CARP&C3a",
            "CARP": "!TNF&NFkB | !TNF&!C3a | NFkBnuc&!C3a",
            "FLIP": "NFkBnuc",
            "IAP": "!TNF&NFkBnuc | !TNF&!C3a | NFkBnuc&!C3a",
            "IKKa": "TNF&!C3a&!A20a",
            "IkB": "!TNF&NFkBnuc | !TNF&!IKKa | NFkBnuc&!IKKa",
            "NFkB": "!IkB",
            "NFkBnuc": "NFkB&!IkB",
            "T2": "TNF&!FLIP",
        },
    },
    "davidich_yeast": {
        "citation": "Davidich & Bornholdt, PLoS ONE 3.2 (2008): e1672",
        "description": "Fission yeast cell cycle (10 nodes)",
        "output_nodes": ["Cdc2_Cdc13"],
        "rules": {
            "Start": "Start",
            "SK": "Start",
            "Cdc2_Cdc13": "!Ste9&!Rum1&!Slp1 | !Ste9&!Rum1&Cdc2_Cdc13_A",
            "Ste9": "!SK&!Cdc2_Cdc13&!Cdc2_Cdc13_A | !SK&!Cdc2_Cdc13&PP | !Cdc2_Cdc13&!Cdc2_Cdc13_A&PP",
            "Rum1": "!SK&!Cdc2_Cdc13&!Cdc2_Cdc13_A | !SK&!Cdc2_Cdc13&PP | !Cdc2_Cdc13&!Cdc2_Cdc13_A&PP",
            "PP": "Slp1",
            "Cdc25": "Cdc2_Cdc13 | Cdc2_Cdc13_A",
            "Slp1": "Cdc2_Cdc13_A",
            "Wee1_Mik1": "!Cdc2_Cdc13 & !Cdc2_Cdc13_A",
            "Cdc2_Cdc13_A": "Cdc2_Cdc13&Cdc25&!Wee1_Mik1",
        },
    },
    "drosophila_cellcycle": {
        "citation": "Faure et al., J Theor Biol 250.2 (2008): 332-343",
        "description": "Drosophila cell cycle (11 dynamic + 3 inputs)",
        "output_nodes": ["CycB"],
        "rules": {
            "Ago": "Ago",
            "CycD": "CycD",
            "Notch": "Notch",
            "CycA": "!Rb&E2F&!Fzy&!Fzr",
            "CycB": "!Fzy&!Fzr&!Wee1 | !Fzy&!Fzr&Wee1&Stg",
            "CycE": "!Rb&E2F&!Dap | !Rb&E2F&Dap&!Ago",
            "Dap": "!CycE&!Notch | CycE",
            "E2F": "!Rb&!CycA&!Rux&!CycB | !Rb&!CycA&Rux | !Rb&CycA&Rux",
            "Fzr": "!CycE&!CycA&!Rux&!CycB | !CycE&!CycA&Rux | !CycE&CycA&Rux | CycE&!CycA&!Rux&!CycB&Notch | CycE&!CycA&Rux&Notch | CycE&CycA&Rux&Notch",
            "Fzy": "!Rux&CycB",
            "Rb": "!CycD&!CycE&!CycA&!Rux&!CycB | !CycD&!CycE&!CycA&Rux | !CycD&!CycE&CycA&Rux",
            "Rux": "!CycD&!CycA&!Rux&!CycB | !CycD&!CycA&Rux | !CycD&CycA&Rux&!CycB | CycD&!CycE&!CycA&!Rux&!CycB | CycD&!CycE&!CycA&Rux | CycD&!CycE&CycA&Rux&!CycB",
            "Stg": "!Rb&!E2F&!Rux&CycB&!Notch | !Rb&E2F&!Notch | Rb&!Rux&CycB&!Notch",
            "Wee1": "!Rux&!CycB | Rux",
        },
    },
    "myeloid_progenitors": {
        "citation": "Krumsiek et al., PLoS ONE 6.7 (2011): e22649",
        "description": "Myeloid progenitor cell fate (11 nodes, no inputs)",
        "output_nodes": ["GATA1"],
        "rules": {
            "CEBPA": "CEBPA & !(GATA1 & FOG1 & SCL)",
            "EGR_NAB": "PU1 & JUN & !GFI1",
            "EKLF": "GATA1 & !FLI1",
            "FLI1": "GATA1 & !EKLF",
            "FOG1": "GATA1",
            "GATA1": "GATA1 | GATA2 | FLI1 & !PU1",
            "GATA2": "GATA2 & !(GATA1 & FOG1) & !PU1",
            "GFI1": "CEBPA & !EGR_NAB",
            "JUN": "PU1 & !GFI1",
            "PU1": "CEBPA | PU1 & !(GATA1 | GATA2)",
            "SCL": "GATA1 & !PU1",
        },
    },
    "blood_stem_cell": {
        "citation": "Thalheim et al., Cells 10.12 (2021): 3533",
        "description": "Blood stem cell heterogeneity (11 nodes, no inputs)",
        "output_nodes": ["GATA1"],
        "rules": {
            "ERG": "(FLI1 | ERG | GATA2 | RUNX1) & !(SCL & ETO2)",
            "ETO2": "SCL & GATA2",
            "FLI1": "(FLI1 | GATA2 | ERG | GATA2&SCL) & !GATA1",
            "GATA1": "SCL & GATA1 & !(PU1 & GATA1)",
            "GATA2": "(FLI1 | ERG | SCL&GATA2) & !(GATA1&ZFPM1) & !(GATA2&HHEX)",
            "HHEX": "FLI1 | ERG | SCL&GATA2",
            "PU1": "(FLI1 | RUNX1 | PU1) & !(PU1 & GATA1)",
            "RUNX1": "(PU1 | RUNX1 | FLI1 | ERG | GATA2&SCL) & !(RUNX1 & SMAD6)",
            "SCL": "FLI1 | ERG | FLI1&GATA2 | SCL&GATA1 | GATA2&SCL",
            "SMAD6": "ERG | FLI1 | SCL&GATA2",
            "ZFPM1": "GATA2 & SCL",
        },
    },
    "emt_switch": {
        "citation": "Steinway et al., PLoS Comput Biol 10.8 (2014): e1003762",
        "description": "Epithelial-mesenchymal transition switch (12 nodes, no inputs)",
        "output_nodes": ["Ecadherin_mRNA"],
        "rules": {
            "Ecadherin_mRNA": "!(ZEB1_H & ZEB1 & SNAI1 & SNAI2 & Twist)",
            "LEF1": "ZEB1 & !miR_200 | ZEB1_H",
            "N_bcatenin_H": "!miR_34 & !miR_200",
            "SNAI1": "(ZEB1_H | ZEB1) & !miR_34 | ZEB1_H & ZEB1",
            "SNAI2": "Twist & (SNAI2 | N_bcatenin_H & LEF1)",
            "TGFb_secr": "b_catenin_TCF4 & !miR_200",
            "Twist": "SNAI1 & !miR_34",
            "ZEB1": "SNAI2 | b_catenin_TCF4 & !miR_200",
            "ZEB1_H": "ZEB1 & N_bcatenin_H & LEF1 & (SNAI2 | !miR_200)",
            "b_catenin_TCF4": "N_bcatenin_H & SNAI1 & SNAI2",
            "miR_200": "!(Twist & ZEB1_H & SNAI1) & !(ZEB1 & !miR_200)",
            "miR_34": "!SNAI1 | !(ZEB1 | ZEB1_H)",
        },
    },
    "fanconi_anemia": {
        "citation": "Rodriguez et al., BMC Syst Biol 9 (2015): 79",
        "description": "Fanconi anemia and checkpoint recovery (15 nodes, no inputs)",
        "output_nodes": ["CHKREC"],
        "rules": {
            "ADD": "(NUC1 | NUC2 | NUC1&NUC2) & !TLS",
            "ATM": "(ATR | DSB) & !(FAcore | CHKREC)",
            "ATR": "(ICL | ATM) & !CHKREC",
            "CHKREC": "(TLS&!DSB | HRR2&!DSB | FAHRR&!DSB | NHEJ&!DSB) | !(ICL|NHEJ|FAHRR|ADD|CHKREC|DSB|HRR2|TLS)",
            "DSB": "(NUC1 | NUC2) & !(FAHRR | HRR2 | NHEJ)",
            "FAHRR": "FANCD2I & DSB & !(NHEJ & CHKREC)",
            "FANCD2I": "FAcore & (DSB&(ATM|ATR) | ATM | ATR) & !CHKREC",
            "FAcore": "ICL & (ATM | ATR) & !CHKREC",
            "HRR2": "TLS&NUC2&DSB & !(NHEJ|CHKREC|FAHRR) | NUC2&ICL&NHEJ&DSB & !(CHKREC|FAHRR)",
            "ICL": "ICL & !DSB",
            "NHEJ": "NUC2 & DSB & !(HRR2 | FAHRR | CHKREC)",
            "NUC1": "FANCD2I & ICL",
            "NUC2": "ICL&(ATM|ATR) & !(FAcore&FANCD2I) | NUC1&ICL&p53 & !(FAcore&FANCD2I)",
            "TLS": "(FAcore&ADD | ADD) & !CHKREC",
            "p53": "(NHEJ | ATM | ATR) & !CHKREC",
        },
    },
    "arabidopsis_cellcycle": {
        "citation": "Ortiz-Gutierrez et al., PLoS Comput Biol 11.9 (2015): e1004486",
        "description": "Arabidopsis thaliana cell cycle (14 nodes, no inputs)",
        "output_nodes": ["CYCB1_1"],
        "rules": {
            "APC_C": "(E2Fa & !(E2Fe | RBR)) | (MYB3R1_4 & !E2Fe) | (MYB77 & !E2Fe)",
            "CDKB1_1": "MYB77 | (E2Fb & (!RBR | CYCD3_1 & !KRP1) & !E2Fc) | MYB3R1_4",
            "CYCA2_3": "(MYB77 | MYB3R1_4) & !APC_C",
            "CYCB1_1": "(MYB77 & !APC_C) | (E2Fb & (!RBR | CYCD3_1 & !KRP1) & !(E2Fc | APC_C)) | (MYB3R1_4 & !APC_C)",
            "CYCD3_1": "!SCF",
            "E2Fa": "(CYCA2_3 & !E2Fc & !CDKB1_1 & !E2Fa) | (E2Fa & !(CDKB1_1 & CYCA2_3)) | (CDKB1_1 & !E2Fc & !CYCA2_3 & !E2Fa) | !(E2Fc | CYCA2_3 | CDKB1_1 | E2Fa)",
            "E2Fb": "E2Fa & !RBR",
            "E2Fc": "(MYB3R1_4 & !(SCF & CYCD3_1 & !KRP1)) | (E2Fa & !(SCF & CYCD3_1 & !KRP1 | RBR))",
            "E2Fe": "((((((((((((((((((CYCD3_1 & ((((!E2Fc & !MYB77) & !E2Fb) & !RBR) & !KRP1)) | ((CYCD3_1 & (E2Fb & RBR)) & ((!E2Fc & !MYB77) & !KRP1))) | ((CYCD3_1 & ((E2Fb & RBR) & KRP1)) & (!E2Fc & !MYB77))) | ((E2Fb & (E2Fc & KRP1)) & ((!MYB77 & !CYCD3_1) & !RBR))) | ((E2Fb & KRP1) & (((!E2Fc & !MYB77) & !CYCD3_1) & !RBR))) | ((E2Fb & RBR) & (((!E2Fc & !MYB77) & !CYCD3_1) & !KRP1))) | (RBR & ((((!E2Fc & !MYB77) & !E2Fb) & !CYCD3_1) & !KRP1))) | (E2Fb & ((((!E2Fc & !MYB77) & !CYCD3_1) & !RBR) & !KRP1))) | MYB77) | ((CYCD3_1 & (RBR & KRP1)) & ((!E2Fc & !MYB77) & !E2Fb))) | (E2Fb & (!RBR | (!KRP1 & CYCD3_1)))) | ((E2Fb & (RBR & KRP1)) & ((!E2Fc & !MYB77) & !CYCD3_1))) | ((CYCD3_1 & KRP1) & (((!E2Fc & !MYB77) & !E2Fb) & !RBR))) | ((CYCD3_1 & (E2Fb & KRP1)) & ((!E2Fc & !MYB77) & !RBR))) | ((RBR & KRP1) & (((!E2Fc & !MYB77) & !E2Fb) & !CYCD3_1))) | ((CYCD3_1 & RBR) & (((!E2Fc & !MYB77) & !E2Fb) & !KRP1))) | (KRP1 & ((((!E2Fc & !MYB77) & !E2Fb) & !CYCD3_1) & !RBR))) | !(((((E2Fc | MYB77) | E2Fb) | CYCD3_1) | RBR) | KRP1))",
            "KRP1": "(MYB77 | MYB3R1_4) & !(CDKB1_1 & CYCA2_3 & SCF)",
            "MYB3R1_4": "MYB77 | (MYB3R1_4 & CYCB1_1 & !KRP1)",
            "MYB77": "E2Fb & (CYCD3_1 & !KRP1 | !RBR)",
            "RBR": "(MYB3R1_4 & (KRP1 | !CYCD3_1)) | (E2Fa & (KRP1 | !CYCD3_1) & !(RBR & (!CYCD3_1 | KRP1)))",
            "SCF": "(MYB3R1_4 & !APC_C) | (E2Fb & (!RBR | CYCD3_1 & !KRP1) & !APC_C)",
        },
    },
}


# ── Input node detection ───────────────────────────────────────────

def identify_input_nodes(rules):
    """Detect nodes whose update rule is self-referential (f(x) = x).

    These persist indefinitely once set — external inputs, not dynamically
    regulated genes. Examples: CycD (G1), TNF (G2), Start (G4).
    """
    return [node for node, expr in rules.items() if expr.strip() == node]


# ── .bnet parser ────────────────────────────────────────────────────

def parse_bnet(text):
    rules = {}
    for line in text.strip().split("\n"):
        line = line.strip()
        if not line or line.startswith("#") or line.lower().startswith("targets"):
            continue
        parts = line.split(",", 1)
        if len(parts) != 2:
            continue
        target = parts[0].strip()
        expression = parts[1].strip()
        rules[target] = expression

    all_refs = set()
    for expr in rules.values():
        all_refs |= set(re.findall(r"[A-Za-z_]\w*", expr)) - {"not", "and", "or"}
    external = all_refs - set(rules.keys())
    for ext in sorted(external):
        rules[ext] = ext

    return rules


# ── Boolean expression compiler ─────────────────────────────────────

def _find_regulators(expression, node_names):
    regs = []
    for name in node_names:
        if re.search(r"\b" + re.escape(name) + r"\b", expression):
            regs.append(name)
    return regs


def _eval_bool_expr(expression, state_dict):
    expr = expression
    for node in sorted(state_dict.keys(), key=len, reverse=True):
        expr = re.sub(
            r"\b" + re.escape(node) + r"\b",
            str(int(state_dict[node])),
            expr,
        )
    expr = expr.replace("!", " not ").replace("&", " and ").replace("|", " or ")
    return int(eval(expr))


def compile_network(rules):
    """Compile Boolean rules into truth tables for vectorized simulation.

    Returns list of (node_name, regulator_indices, truth_table) tuples,
    plus the ordered list of node names.
    """
    node_names = list(rules.keys())
    name_to_idx = {name: i for i, name in enumerate(node_names)}
    n = len(node_names)

    compiled = []
    for node in node_names:
        expression = rules[node]
        regs = _find_regulators(expression, node_names)
        reg_indices = [name_to_idx[r] for r in regs]
        k = len(regs)

        truth_table = np.zeros(2**k, dtype=np.int8)
        for bits in range(2**k):
            state = {}
            for j, reg_name in enumerate(regs):
                state[reg_name] = (bits >> j) & 1
            truth_table[bits] = _eval_bool_expr(expression, state)

        compiled.append((node, reg_indices, truth_table))

    return compiled, node_names


def extract_interaction_graph(rules):
    """Extract the wiring diagram: which nodes regulate which, and how.

    Returns adjacency dict: {target: [(regulator, sign), ...]}
    where sign is +1 (activator), -1 (repressor), or 0 (ambiguous).
    """
    node_names = list(rules.keys())
    graph = {}

    for target in node_names:
        expression = rules[target]
        regs = _find_regulators(expression, node_names)
        edges = []
        for reg in regs:
            if reg == target:
                continue
            state_base = {n: 0 for n in regs}
            effects = []
            for bits in range(2 ** (len(regs) - 1)):
                other_regs = [r for r in regs if r != reg]
                state_off = dict(state_base)
                state_on = dict(state_base)
                for j, other in enumerate(other_regs):
                    val = (bits >> j) & 1
                    state_off[other] = val
                    state_on[other] = val
                state_off[reg] = 0
                state_on[reg] = 1
                val_off = _eval_bool_expr(expression, state_off)
                val_on = _eval_bool_expr(expression, state_on)
                effects.append(val_on - val_off)

            pos = sum(1 for e in effects if e > 0)
            neg = sum(1 for e in effects if e < 0)
            if pos > 0 and neg == 0:
                sign = 1
            elif neg > 0 and pos == 0:
                sign = -1
            else:
                sign = 0
            edges.append({"regulator": reg, "sign": sign})
        graph[target] = edges

    return graph


def extract_rule_fourier(rules):
    """Extract exact interaction structure from Boolean rules via local WHT.

    For each node's update rule (a Boolean function of k regulators),
    takes the Walsh-Hadamard transform of the 2^k truth table. Nonzero
    coefficients at order >= 2 indicate genuine multi-way interactions
    among regulators, independent of how the rule is written.

    This replaces DNF clause-splitting, which is representation-dependent
    (e.g. a&b | a&!b reduces to a, but clause-splitting reports a
    spurious (a,b) interaction).

    Returns dict with:
      pairwise: list of [node_i, node_j] pairs with nonzero local w_{i,j}
      triples: list of [node_i, node_j, node_k] triples with nonzero w_{i,j,k}
      per_gene: per-node details (regulators, nonzero coefficients, magnitudes)
      pairwise_strengths: {(i,j): sum of |local w| across all genes}
    """
    compiled, node_names = compile_network(rules)

    pairwise = set()
    triples = set()
    pairwise_strength = {}
    per_gene = {}

    max_order = max(len(_find_regulators(rules[n], list(rules.keys()))) for n in rules)
    local_energy = np.zeros(max_order + 1, dtype=np.float64)

    for node, reg_indices, truth_table in compiled:
        regs = [node_names[ri] for ri in reg_indices]
        k = len(regs)

        if k == 0:
            per_gene[node] = {"regulators": regs, "interactions": []}
            continue

        w = wht(truth_table.astype(np.float64)) / (2 ** k)

        for T in range(2 ** k):
            order = bin(T).count("1")
            if order <= max_order:
                local_energy[order] += w[T] ** 2

        interactions = []
        for T in range(1, 2 ** k):
            if abs(w[T]) < 1e-12:
                continue
            bits = [j for j in range(k) if T & (1 << j)]
            order = len(bits)
            if order < 2:
                continue

            reg_names = tuple(sorted(regs[j] for j in bits))
            interactions.append({
                "regulators": list(reg_names),
                "order": order,
                "coefficient": float(w[T]),
            })

            key = frozenset(reg_names)
            if order == 2:
                pairwise.add(key)
                pairwise_strength[key] = pairwise_strength.get(key, 0.0) + abs(w[T])
            elif order == 3:
                triples.add(key)

        per_gene[node] = {"regulators": regs, "interactions": interactions}

    total_energy = local_energy.sum()
    local_spectrum = (local_energy / total_energy).tolist() if total_energy > 0 else local_energy.tolist()

    return {
        "pairwise": [sorted(p) for p in sorted(pairwise)],
        "triples": [sorted(t) for t in sorted(triples)],
        "per_gene": per_gene,
        "n_pairwise": len(pairwise),
        "n_triples": len(triples),
        "pairwise_strengths": {
            str(sorted(k)): v for k, v in sorted(pairwise_strength.items())
        },
        "local_energy_spectrum": local_spectrum,
    }


# ── Vectorized synchronous simulation ───────────────────────────────

def update_all(states, compiled):
    """One synchronous step across all states.

    states: (N, n_nodes) int8 array
    Returns: (N, n_nodes) updated states
    """
    N = states.shape[0]
    new_states = np.empty_like(states)

    for col, (_, reg_indices, truth_table) in enumerate(compiled):
        if len(reg_indices) == 0:
            new_states[:, col] = truth_table[0]
            continue
        idx = np.zeros(N, dtype=np.int64)
        for j, ri in enumerate(reg_indices):
            idx += states[:, ri].astype(np.int64) << j
        new_states[:, col] = truth_table[idx]

    return new_states


def simulate_sync_output(states, compiled, clamp_mask, clamp_value,
                         output_indices, max_steps=200):
    """Run synchronous dynamics and extract output, handling limit cycles.

    For fixed points, returns the output node value directly.
    For limit cycles, averages the output over one full cycle period.

    Returns: (N_init,) float64 output values, dict with convergence stats
    """
    n_init, n_nodes = states.shape
    current = states.copy()
    current[:, clamp_mask] = clamp_value

    converged = np.zeros(n_init, dtype=bool)

    for step in range(max_steps):
        new = update_all(current, compiled)
        new[:, clamp_mask] = clamp_value
        newly_converged = ~converged & np.all(new == current, axis=1)
        converged |= newly_converged
        current = new
        if converged.all():
            break

    output = np.zeros(n_init, dtype=np.float64)
    for oi in output_indices:
        output += current[:, oi].astype(np.float64)
    output /= len(output_indices)

    n_cycling = int((~converged).sum())
    if n_cycling > 0:
        cycling_idx = np.where(~converged)[0]
        for i in cycling_idx:
            anchor = current[i].copy()
            cycle_output = []
            test_state = anchor.copy()
            for _ in range(max_steps):
                val = sum(float(test_state[oi]) for oi in output_indices) / len(output_indices)
                cycle_output.append(val)
                test_state = update_all(test_state.reshape(1, -1), compiled)[0]
                test_state[clamp_mask] = clamp_value
                if np.array_equal(test_state, anchor):
                    output[i] = np.mean(cycle_output)
                    break

    info = {
        "update_scheme": "sync",
        "n_fixed_point": int(converged.sum()),
        "n_cycling": n_cycling,
        "convergence_rate": float(converged.mean()),
    }
    return output, info


def update_async(states, compiled, clamp_mask, clamp_value, rng):
    """One asynchronous step: update all non-clamped nodes in random order.

    Each step applies a random permutation of the active nodes, updating
    each node once using the current state (which includes prior updates
    within the same permutation).
    """
    N = states.shape[0]
    active_nodes = np.where(~clamp_mask)[0]
    if len(active_nodes) == 0:
        return states.copy()

    current = states.copy()
    order = rng.permutation(active_nodes)

    for col in order:
        _, reg_indices, truth_table = compiled[col]
        if len(reg_indices) == 0:
            current[:, col] = truth_table[0]
            continue
        idx = np.zeros(N, dtype=np.int64)
        for j, ri in enumerate(reg_indices):
            idx += current[:, ri].astype(np.int64) << j
        current[:, col] = truth_table[idx]

    current[:, clamp_mask] = clamp_value
    return current


def simulate_async_output(states, compiled, clamp_mask, clamp_value,
                          output_indices, max_steps=200, n_replicates=10,
                          seed=42):
    """Run async dynamics: n_replicates independent trajectories, averaged.

    Each replicate uses a different random permutation sequence. For
    each, runs until fixed point or max_steps. Non-converged trajectories
    average the output over the last 50 steps.
    """
    n_init = states.shape[0]
    output_sum = np.zeros(n_init, dtype=np.float64)
    tail_window = min(50, max_steps)

    for rep in range(n_replicates):
        rng = np.random.default_rng(seed + rep * 997)
        current = states.copy()
        current[:, clamp_mask] = clamp_value

        converged = np.zeros(n_init, dtype=bool)
        maybe_converged = np.zeros(n_init, dtype=bool)
        tail_sum = np.zeros(n_init, dtype=np.float64)
        tail_count = 0

        for step in range(max_steps):
            new = update_async(current, compiled, clamp_mask, clamp_value, rng)
            no_change = ~converged & np.all(new == current, axis=1)
            newly_converged = no_change & maybe_converged
            converged |= newly_converged
            maybe_converged = no_change
            current = new

            if step >= max_steps - tail_window:
                for oi in output_indices:
                    tail_sum += current[:, oi].astype(np.float64)
                tail_count += len(output_indices)

            if converged.all():
                break

        rep_output = np.zeros(n_init, dtype=np.float64)
        for oi in output_indices:
            rep_output += current[:, oi].astype(np.float64)
        rep_output /= len(output_indices)

        if (~converged).any() and tail_count > 0:
            rep_output[~converged] = tail_sum[~converged] / tail_count

        output_sum += rep_output

    output = output_sum / n_replicates

    info = {
        "update_scheme": "async",
        "max_steps": max_steps,
        "n_replicates": n_replicates,
    }
    return output, info


# ── Coalition sweep ─────────────────────────────────────────────────

def sweep_coalitions(rules, output_nodes, n_init=512, max_steps=200, seed=42,
                     update_scheme="sync", clamp_value=0,
                     async_replicates=10, exclude_inputs=False,
                     input_config=None):
    """Generate exhaustive coalition table for a Boolean network.

    For each of the 2^n_players coalitions (subsets of active player nodes):
    - Player nodes NOT in the coalition are clamped to clamp_value
    - Input nodes (f(x)=x) are optionally excluded from the player set
      and fixed to input_config (0 or 1) for the entire sweep
    - Run from n_init random initial states under the chosen update scheme
    - Record output node values at the attractor (sync) or time-average (async)

    Returns dict with coalition table arrays matching the transformer format.
    """
    compiled, node_names = compile_network(rules)
    n_total = len(node_names)
    name_to_idx = {name: i for i, name in enumerate(node_names)}
    output_indices = [name_to_idx[o] for o in output_nodes]

    input_nodes = identify_input_nodes(rules) if exclude_inputs else []
    input_indices = [name_to_idx[inp] for inp in input_nodes]
    player_indices = [i for i in range(n_total) if i not in input_indices]
    player_names = [node_names[i] for i in player_indices]
    n_players = len(player_indices)
    N = 2**n_players

    rng = np.random.default_rng(seed)
    init_states = rng.integers(0, 2, size=(n_init, n_total), dtype=np.int8)

    if input_indices and input_config is not None:
        init_states[:, input_indices] = input_config

    values = np.zeros((N, n_init), dtype=np.float64)
    total_cycling = 0
    total_fixed = 0

    desc = f"Coalitions ({update_scheme}, clamp={clamp_value})"
    if input_nodes:
        desc += f", inputs={input_nodes}={input_config}"
    for coalition in tqdm(range(N), desc=desc):
        clamp_mask = np.zeros(n_total, dtype=bool)
        for bit_pos, node_idx in enumerate(player_indices):
            if not (coalition & (1 << bit_pos)):
                clamp_mask[node_idx] = True

        if update_scheme == "sync":
            output, info = simulate_sync_output(
                init_states, compiled, clamp_mask, clamp_value,
                output_indices, max_steps=max_steps,
            )
            total_fixed += info["n_fixed_point"]
            total_cycling += info["n_cycling"]
        else:
            output, info = simulate_async_output(
                init_states, compiled, clamp_mask, clamp_value,
                output_indices, max_steps=max_steps,
                n_replicates=async_replicates, seed=seed + coalition,
            )

        values[coalition] = output

    convergence = {
        "update_scheme": update_scheme,
        "clamp_value": clamp_value,
        "total_simulations": N * n_init,
    }
    if update_scheme == "sync":
        convergence["total_fixed_point"] = total_fixed
        convergence["total_cycling"] = total_cycling
        convergence["cycling_fraction"] = total_cycling / (N * n_init)

    return {
        "values": values,
        "node_names": node_names,
        "player_names": player_names,
        "output_nodes": output_nodes,
        "n_players": n_players,
        "n_total_nodes": n_total,
        "n_init_states": n_init,
        "input_nodes": input_nodes,
        "input_config": input_config,
        "convergence": convergence,
    }


def main():
    parser = argparse.ArgumentParser(
        description="Generate exhaustive coalition tables from Boolean GRN models"
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--model", choices=list(BUILTIN_MODELS.keys()),
                       help="Built-in model name")
    group.add_argument("--bnet-file", help="Path to .bnet file")
    parser.add_argument("--output", required=True, help="Output .npz path")
    parser.add_argument("--output-node", nargs="+", default=None,
                        help="Output node name(s) for readout (required for --bnet-file)")
    parser.add_argument("--n-init", type=int, default=512,
                        help="Number of random initial states (default: 512)")
    parser.add_argument("--max-steps", type=int, default=200)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--update-scheme", choices=["sync", "async"], default="sync",
                        help="Update scheme: sync (default) or async")
    parser.add_argument("--clamp-value", type=int, choices=[0, 1], default=0,
                        help="Value to clamp knocked-out nodes to: 0=knockout, 1=constitutive")
    parser.add_argument("--async-replicates", type=int, default=10,
                        help="Independent async trajectories per initial state (default: 10)")
    parser.add_argument("--exclude-inputs", action="store_true",
                        help="Exclude input nodes (f(x)=x) from the player set")
    parser.add_argument("--input-config", type=int, choices=[0, 1], default=None,
                        help="Clamp value for excluded input nodes (required with --exclude-inputs)")
    parser.add_argument("--wiring-output", default=None,
                        help="Path to save wiring diagram JSON (second ground truth)")
    args = parser.parse_args()

    if args.exclude_inputs and args.input_config is None:
        parser.error("--input-config required when using --exclude-inputs")

    if args.model:
        model_info = BUILTIN_MODELS[args.model]
        rules = model_info["rules"]
        output_nodes = args.output_node or model_info["output_nodes"]
        model_name = args.model
        citation = model_info["citation"]
        description = model_info["description"]
        print(f"[{timestamp()}] Model: {args.model} — {description}")
        print(f"  Citation: {citation}")
    else:
        bnet_text = Path(args.bnet_file).read_text()
        rules = parse_bnet(bnet_text)
        output_nodes = args.output_node
        if not output_nodes:
            parser.error("--output-node required when using --bnet-file")
        model_name = Path(args.bnet_file).stem
        citation = "custom"
        description = f"Custom model from {args.bnet_file}"
        print(f"[{timestamp()}] Loading custom model: {args.bnet_file}")

    node_names = list(rules.keys())
    n = len(node_names)
    input_nodes = identify_input_nodes(rules)

    print(f"  Nodes ({n}): {node_names}")
    print(f"  Output: {output_nodes}")
    if input_nodes:
        print(f"  Input nodes (f(x)=x): {input_nodes}")
        if args.exclude_inputs:
            n_players = n - len(input_nodes)
            print(f"  Players (excluding inputs): {n_players}")
            print(f"  Input config: {args.input_config}")
            print(f"  Coalitions: {2**n_players}")
        else:
            print(f"  Coalitions: {2**n} (inputs included as players)")
    else:
        print(f"  Coalitions: {2**n}")
    print(f"  Initial states: {args.n_init}")

    for out in output_nodes:
        if out not in rules:
            parser.error(f"Output node '{out}' not in model nodes: {list(rules.keys())}")

    print(f"\n[{timestamp()}] Extracting interaction graph (ground truth wiring)...")
    wiring = extract_interaction_graph(rules)
    n_edges = sum(len(edges) for edges in wiring.values())
    print(f"  {n_edges} regulatory edges found")

    print(f"[{timestamp()}] Extracting local rule Fourier structure...")
    rule_fourier = extract_rule_fourier(rules)
    print(f"  {rule_fourier['n_pairwise']} pairwise interactions (exact local WHT)")
    print(f"  {rule_fourier['n_triples']} triple interactions (exact local WHT)")

    print(f"\n[{timestamp()}] Running exhaustive coalition sweep...")
    print(f"  Update scheme: {args.update_scheme}")
    print(f"  Clamp value: {args.clamp_value} ({'knockout' if args.clamp_value == 0 else 'constitutive activation'})")
    result = sweep_coalitions(
        rules, output_nodes,
        n_init=args.n_init, max_steps=args.max_steps, seed=args.seed,
        update_scheme=args.update_scheme, clamp_value=args.clamp_value,
        async_replicates=args.async_replicates,
        exclude_inputs=args.exclude_inputs,
        input_config=args.input_config,
    )

    v_mean = result["values"].mean(axis=1)
    n_unique = len(np.unique(np.round(v_mean, 6)))
    total_var = float(np.var(v_mean))
    n_players = result["n_players"]
    player_names = result["player_names"]

    print(f"\n[{timestamp()}] Value function stats:")
    print(f"  v(empty):    {v_mean[0]:.4f}")
    print(f"  v(all):      {v_mean[-1]:.4f}")
    print(f"  mean:        {v_mean.mean():.4f}")
    print(f"  std:         {v_mean.std():.4f}")
    print(f"  unique vals: {n_unique}")
    print(f"  total var:   {total_var:.6e}")

    conv = result["convergence"]
    if conv["update_scheme"] == "sync":
        print(f"  fixed point: {conv['total_fixed_point']}/{conv['total_simulations']}")
        print(f"  cycling:     {conv['total_cycling']}/{conv['total_simulations']} ({conv['cycling_fraction']:.1%})")

    variance_floor_warning = False
    if n_unique < 20:
        print(f"  WARNING: Only {n_unique} unique values — energy spectrum fractions are numerically unstable")
        variance_floor_warning = True
    if total_var < 1e-6:
        print(f"  WARNING: Total variance {total_var:.2e} below threshold — model may be degenerate")
        variance_floor_warning = True

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        output_path,
        target_logits=result["values"],
        foil_logits=np.zeros_like(result["values"]),
        coalition_indices=np.arange(2**n_players, dtype=np.int64),
        circuit_heads=np.array(player_names, dtype=object),
        n_players=np.int64(n_players),
        n_prompts=np.int64(args.n_init),
        model_name=model_name,
        output_nodes=np.array(output_nodes, dtype=object),
        update_scheme=args.update_scheme,
        clamp_value=np.int64(args.clamp_value),
        variance_floor_warning=variance_floor_warning,
        n_unique_values=np.int64(n_unique),
        total_variance=np.float64(total_var),
        input_nodes=np.array(result["input_nodes"], dtype=object),
        input_config=np.int64(args.input_config) if args.input_config is not None else np.int64(-1),
        exclude_inputs=args.exclude_inputs,
    )
    print(f"  Saved coalition table: {output_path}")

    wiring_path = args.wiring_output
    if wiring_path is None:
        wiring_path = str(output_path).replace(".npz", "_wiring.json")
    wiring_output = {
        "model": model_name,
        "citation": citation,
        "description": description,
        "n_nodes": n,
        "node_names": node_names,
        "player_names": player_names,
        "output_nodes": output_nodes,
        "input_nodes": result["input_nodes"],
        "input_config": args.input_config,
        "exclude_inputs": args.exclude_inputs,
        "interaction_graph": wiring,
        "n_edges": n_edges,
        "rule_fourier": rule_fourier,
        "update_scheme": args.update_scheme,
        "clamp_value": args.clamp_value,
        "convergence": conv,
        "variance_floor_warning": variance_floor_warning,
        "n_unique_values": n_unique,
        "total_variance": total_var,
        "timestamp": timestamp(),
    }
    with open(wiring_path, "w") as f:
        json.dump(wiring_output, f, indent=2)
    print(f"  Saved wiring diagram: {wiring_path}")

    print(f"\n[{timestamp()}] Done.")


if __name__ == "__main__":
    main()
