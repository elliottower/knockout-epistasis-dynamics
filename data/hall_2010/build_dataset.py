"""
Build the Hall, Agan & Pope (2010) fitness epistasis dataset.

Source: Hall DW, Agan M, Pope SC. 2010. Fitness epistasis among 6 biosynthetic
loci in the budding yeast Saccharomyces cerevisiae. J Heredity 101(suppl_1):S75-S84.
DOI: 10.1093/jhered/esq007

Data extracted from the harmslab/notebooks-nonlinear-high-order-epistasis GitHub
repository (https://github.com/harmslab/notebooks-nonlinear-high-order-epistasis),
which was used in:
  Sailer ZR, Harms MJ. 2017. Detecting high-order epistasis in nonlinear
  genotype-phenotype maps. Genetics 205(3):1079-1088.

The Sailer & Harms paper acknowledges David Hall for providing the complete data
(datasets VI and VII in their paper).

GENE NAME ASSIGNMENT:
The harmslab data uses positional indices 0-5 without gene names. The original
paper (behind paywall) identifies these as 6 of the 7 standard S. cerevisiae
auxotrophic biosynthetic markers from the BY strain background (Brachmann et al.
1998): ADE2, HIS3, LEU2, LYS2, MET15, TRP1, URA3.

The EXACT mapping of positions 0-5 to gene names requires the original paper.
We label them as locus_0 through locus_5 in the output. These are gene DELETION
knockouts -- "1" = wildtype (functional gene), "0" = knockout (deleted gene).

GENOTYPE ENCODING:
- Binary string of length 6, e.g., "111111" = wildtype (all genes functional)
- "0" at position i means locus i is knocked out
- "000000" = all 6 genes knocked out

FITNESS COMPONENTS (4 measured):
1. haploid_growth_rate - growth rate of haploid strains
2. diploid_growth_rate - growth rate of homozygous diploid strains
3. mating_efficiency - haploid mating efficiency
4. sporulation_efficiency - diploid sporulation efficiency

All values are relative to wildtype. The data includes log-transformed values
(the harmslab repository notes log_transform: true for all datasets).

For sporulation efficiency, log(1+fitness) was used since some raw values were zero.
"""

import json
import os

# Raw data from harmslab repository
# Each dataset has: genotypes, phenotypes (mean), stdeviations, n_replicates

HAPLOID_GROWTH = {
    "genotypes": [
        "111111", "111110", "111101", "111011", "110111", "101111", "011111",
        "111100", "111010", "111001", "110110", "110101", "110011", "101110",
        "101101", "101011", "100111", "011110", "011101", "011011", "010111",
        "001111", "111000", "110100", "110010", "110001", "101100", "101010",
        "101001", "100110", "100101", "100011", "011100", "011010", "011001",
        "010110", "010101", "010011", "001110", "001101", "001011", "000111",
        "110000", "101000", "100100", "100010", "100001", "011000", "010100",
        "010010", "010001", "001100", "001010", "001001", "000110", "000101",
        "000011", "100000", "010000", "001000", "000100", "000010", "000001",
        "000000"
    ],
    "phenotypes": [
        1.017502524, 0.962880306, 0.970186617, 1.020017150, 1.033741380, 0.986513217,
        0.818847212, 0.965518781, 0.998638388, 0.982267866, 0.989507080, 0.973392049,
        0.975981305, 0.957856199, 1.001578445, 1.034003512, 1.023197519, 0.809682791,
        0.830424627, 0.900205964, 0.926179370, 0.876069972, 0.935346601, 0.949474916,
        0.991830634, 0.985172468, 1.015666444, 0.986199414, 0.961311461, 0.997582101,
        0.995037749, 0.996480637, 0.841439517, 0.861259728, 0.868121387, 0.844085387,
        0.889491684, 0.785200634, 0.849394513, 0.827782925, 0.811695376, 0.773902264,
        0.985981885, 0.995262062, 0.974855411, 0.957577220, 0.976160612, 0.818688552,
        0.836375885, 0.712625109, 0.794632577, 0.875824798, 0.728634447, 0.866675480,
        0.751475727, 0.868890078, 0.822795289, 0.968103209, 0.864458991, 0.802302365,
        0.806531009, 0.845905946, 0.877606558, 0.872302955
    ],
    "stdeviations": [
        0.057827954, 0.024590893, 0.034205426, 0.067228420, 0.030038917, 0.060975988,
        0.051228136, 0.051466669, 0.073701591, 0.045159856, 0.049438525, 0.027109622,
        0.046498315, 0.076502073, 0.043169932, 0.080646569, 0.039700019, 0.026366135,
        0.056572857, 0.052456292, 0.063515762, 0.036943085, 0.028873284, 0.042015152,
        0.043756423, 0.049377825, 0.053574826, 0.050427576, 0.039690270, 0.022239580,
        0.036321880, 0.063499760, 0.045702627, 0.036302426, 0.055972192, 0.062243356,
        0.040877577, 0.037813782, 0.070959912, 0.036940479, 0.053528075, 0.071521361,
        0.055842617, 0.051686000, 0.037666693, 0.026053979, 0.044582577, 0.040262736,
        0.026681951, 0.048175883, 0.033742095, 0.031341599, 0.050793352, 0.029547307,
        0.028995455, 0.048111283, 0.041115279, 0.057158890, 0.065647552, 0.042202336,
        0.036552504, 0.051182774, 0.046797738, 0.027150360
    ],
    "n_replicates": 10
}

DIPLOID_GROWTH = {
    "genotypes": [
        "111111", "111110", "111101", "111011", "110111", "101111", "011111",
        "111100", "111010", "111001", "110110", "110101", "110011", "101110",
        "101101", "101011", "100111", "011110", "011101", "011011", "010111",
        "001111", "111000", "110100", "110010", "110001", "101100", "101010",
        "101001", "100110", "100101", "100011", "011100", "011010", "011001",
        "010110", "010101", "010011", "001110", "001101", "001011", "000111",
        "110000", "101000", "100100", "100010", "100001", "011000", "010100",
        "010010", "010001", "001100", "001010", "001001", "000110", "000101",
        "000011", "100000", "010000", "001000", "000100", "000010", "000001",
        "000000"
    ],
    "phenotypes": [
        1.034522633, 1.026586726, 1.017473395, 0.952735149, 1.017330988, 1.034462098,
        0.794367560, 0.931867793, 0.994957239, 0.888565504, 0.965518060, 0.895727944,
        1.029944551, 0.891500855, 0.907931706, 1.001330941, 1.007757387, 0.843431775,
        0.838071280, 0.896882946, 0.950558531, 0.873678837, 0.883903372, 0.921787287,
        0.987984102, 0.941149756, 0.983195375, 0.959300932, 0.991816054, 0.974051400,
        0.932158878, 0.990386772, 0.740364694, 0.838400703, 0.841337910, 0.784381959,
        0.825578903, 0.765667155, 0.803731520, 0.778770294, 0.869890142, 0.757682550,
        0.935532860, 0.952064995, 0.923726275, 0.967007263, 0.978219665, 0.724512595,
        0.738600747, 0.866740715, 0.731135545, 0.819333757, 0.813574410, 0.833000877,
        0.869529723, 0.813272126, 0.850612978, 0.983919038, 0.838290476, 0.705999564,
        0.713362526, 0.831060544, 0.806422544, 0.783539183
    ],
    "stdeviations": [
        0.051865178, 0.040237655, 0.055212208, 0.053673517, 0.034272913, 0.027847893,
        0.055225415, 0.045876191, 0.050039530, 0.043146360, 0.033338594, 0.028695337,
        0.032208399, 0.034154449, 0.036434929, 0.044085773, 0.050788030, 0.050958562,
        0.070048367, 0.046860584, 0.056977126, 0.055280001, 0.032476252, 0.033524394,
        0.070728949, 0.050102102, 0.063085730, 0.052738071, 0.027351978, 0.069666726,
        0.062623241, 0.042127645, 0.070723738, 0.053795584, 0.050264800, 0.034009618,
        0.041070500, 0.046790577, 0.063985009, 0.034756024, 0.034632929, 0.064412905,
        0.028353955, 0.064876228, 0.033236211, 0.056241242, 0.045316388, 0.055204933,
        0.046864458, 0.048109144, 0.061078735, 0.052287195, 0.048032564, 0.053751130,
        0.051611572, 0.047638141, 0.067915480, 0.053238777, 0.064681821, 0.076262697,
        0.064151197, 0.038252029, 0.055268567, 0.036580740
    ],
    "n_replicates": 10
}

MATING_EFFICIENCY = {
    "genotypes": [
        "111111", "111110", "111101", "111011", "110111", "101111", "011111",
        "111100", "111010", "111001", "110110", "110101", "110011", "101110",
        "101101", "101011", "100111", "011110", "011101", "011011", "010111",
        "001111", "111000", "110100", "110010", "110001", "101100", "101010",
        "101001", "100110", "100101", "100011", "011100", "011010", "011001",
        "010110", "010101", "010011", "001110", "001101", "001011", "000111",
        "110000", "101000", "100100", "100010", "100001", "011000", "010100",
        "010010", "010001", "001100", "001010", "001001", "000110", "000101",
        "000011", "100000", "010000", "001000", "000100", "000010", "000001",
        "000000"
    ],
    "phenotypes": [
        1.000000000, 0.447655147, 0.638674801, 0.993473526, 0.970738073, 0.959106451,
        0.335432787, 0.393634459, 0.715797579, 0.737100887, 0.723673467, 1.063473448,
        1.002409301, 0.448414750, 0.961183380, 1.136994845, 1.263448612, 0.308878115,
        1.072838812, 0.708428410, 0.472041940, 0.511019569, 0.243856800, 0.477331157,
        0.833250151, 1.162134198, 0.836345270, 0.915738403, 0.958130037, 0.684394111,
        1.118730670, 1.096764315, 0.358163915, 0.397224914, 0.184133542, 0.465707999,
        0.077658742, 0.424548645, 0.434484235, 0.498208613, 0.997341717, 0.328936648,
        0.559761005, 0.540917620, 0.658664819, 0.829613671, 1.303497443, 0.698902457,
        0.770729136, 0.788671075, 0.793359178, 0.443309316, 0.733926168, 0.674489159,
        0.647006345, 0.549899998, 0.559565947, 0.407353539, 0.451299729, 0.728326516,
        0.512251455, 0.465794866, 0.700034652, 0.557833421
    ],
    "stdeviations": [
        0.126347391, 0.018439154, 0.068037611, 0.093897532, 0.092905423, 0.053463816,
        0.167079081, 0.054819181, 0.089089714, 0.183339281, 0.054579813, 0.069366873,
        0.240632726, 0.156530397, 0.139102054, 0.034955794, 0.045621374, 0.057288321,
        0.096736034, 0.103445949, 0.062911740, 0.102527703, 0.076093843, 0.026765174,
        0.057040314, 0.077113360, 0.070107251, 0.216525833, 0.201077356, 0.107373544,
        0.057237732, 0.087586127, 0.037356391, 0.160029542, 0.034880196, 0.036907995,
        0.036907995, 0.091558610, 0.079277335, 0.070017578, 0.047268016, 0.125600788,
        0.058463759, 0.021500282, 0.097804380, 0.176002630, 0.055287589, 0.022282024,
        0.234893777, 0.085154637, 0.075191117, 0.091678208, 0.030641631, 0.089848191,
        0.157566328, 0.123034494, 0.086914067, 0.029697236, 0.066261258, 0.238772446,
        0.037288204, 0.082077777, 0.247781302, 0.054229687
    ],
    "n_replicates": 3
}

SPORULATION_EFFICIENCY = {
    "genotypes": [
        "111111", "111110", "111101", "111011", "110111", "101111", "011111",
        "111100", "111010", "111001", "110110", "110101", "110011", "101110",
        "101101", "101011", "100111", "011110", "011101", "011011", "010111",
        "001111", "111000", "110100", "110010", "110001", "101100", "101010",
        "101001", "100110", "100101", "100011", "011100", "011010", "011001",
        "010110", "010101", "010011", "001110", "001101", "001011", "000111",
        "110000", "101000", "100100", "100010", "100001", "011000", "010100",
        "010010", "010001", "001100", "001010", "001001", "000110", "000101",
        "000011", "100000", "010000", "001000", "000100", "000010", "000001",
        "000000"
    ],
    "phenotypes": [
        2.000000000, 2.007077667, 1.800506288, 1.773466253, 1.880976223, 2.032417118,
        1.230768880, 1.537633442, 1.448885249, 1.675422577, 1.460695783, 1.529213892,
        1.017445900, 1.826269696, 1.791341960, 1.921067632, 1.396131240, 1.216476066,
        1.000000000, 1.307046808, 1.380026207, 1.892345013, 1.862577899, 1.000000000,
        1.884800660, 1.523439415, 1.683481890, 1.685946935, 1.462201026, 1.655843137,
        1.000000000, 1.700791408, 1.000000000, 1.063564278, 1.005760611, 1.147003462,
        1.000000000, 1.000000000, 1.110826526, 1.000000000, 1.127290175, 1.193310018,
        1.600617460, 1.596130850, 1.678934122, 1.478311472, 1.467340839, 1.000000000,
        1.000000000, 1.018203529, 1.000000000, 1.006106513, 1.006177851, 1.000000000,
        1.126892288, 1.012712117, 1.384300886, 1.492251297, 1.000000000, 1.012343318,
        1.005075087, 1.205119932, 1.000000000, 1.000000000
    ],
    "stdeviations": [
        0.017343001, 0.032948128, 0.003172019, 0.021948841, 0.024347720, 0.053828283,
        0.013071988, 0.010857200, 0.004156192, 0.006962362, 0.013621163, 0.005051591,
        0.002119305, 0.016167953, 0.033691198, 0.015743313, 0.009917713, 0.032694884,
        0.000000000, 0.014414236, 0.024242913, 0.013673816, 0.018531621, 0.000000000,
        0.007316327, 0.003241261, 0.028018582, 0.025184307, 0.013772544, 0.003468407,
        0.000000000, 0.015134972, 0.000000000, 0.004681732, 0.001257844, 0.008487365,
        0.000000000, 0.000000000, 0.010451162, 0.000000000, 0.005624836, 0.004927489,
        0.007120167, 0.005298101, 0.005264799, 0.010144931, 0.009393359, 0.000000000,
        0.000000000, 0.002355092, 0.000000000, 0.001333372, 0.001348949, 0.000000000,
        0.002964160, 0.002775722, 0.005184531, 0.018761875, 0.000000000, 0.001347930,
        0.001108158, 0.022394275, 0.000000000, 0.000000000
    ],
    "n_replicates": 3
}

# The 6 loci are biosynthetic genes from the S. cerevisiae BY strain background.
# The 7 standard markers available are: ADE2, HIS3, LEU2, LYS2, MET15, TRP1, URA3.
# The paper uses 6 of these 7, but the positional mapping requires the original paper.
LOCUS_NAMES = ["locus_0", "locus_1", "locus_2", "locus_3", "locus_4", "locus_5"]

# Candidate genes (6 of these 7 were used; exact mapping unknown without paper)
CANDIDATE_GENES = ["ADE2", "HIS3", "LEU2", "LYS2", "MET15", "TRP1", "URA3"]


def build_fitness_data():
    """Build the combined fitness data for all 64 genotypes."""
    records = []

    # Use haploid growth genotype list as canonical ordering
    for i, genotype in enumerate(HAPLOID_GROWTH["genotypes"]):
        # Determine which loci are knocked out
        knocked_out = []
        for j, bit in enumerate(genotype):
            if bit == "0":
                knocked_out.append(LOCUS_NAMES[j])

        n_mutations = genotype.count("0")

        record = {
            "genotype": genotype,
            "n_mutations": n_mutations,
            "loci_knocked_out": knocked_out,
            "haploid_growth_rate": {
                "mean": HAPLOID_GROWTH["phenotypes"][i],
                "std": HAPLOID_GROWTH["stdeviations"][i],
                "n_replicates": HAPLOID_GROWTH["n_replicates"],
            },
            "diploid_growth_rate": {
                "mean": DIPLOID_GROWTH["phenotypes"][i],
                "std": DIPLOID_GROWTH["stdeviations"][i],
                "n_replicates": DIPLOID_GROWTH["n_replicates"],
            },
            "mating_efficiency": {
                "mean": MATING_EFFICIENCY["phenotypes"][i],
                "std": MATING_EFFICIENCY["stdeviations"][i],
                "n_replicates": MATING_EFFICIENCY["n_replicates"],
            },
            "sporulation_efficiency": {
                "mean": SPORULATION_EFFICIENCY["phenotypes"][i],
                "std": SPORULATION_EFFICIENCY["stdeviations"][i],
                "n_replicates": SPORULATION_EFFICIENCY["n_replicates"],
            },
        }
        records.append(record)

    return {
        "metadata": {
            "source": "Hall DW, Agan M, Pope SC. 2010. Fitness epistasis among 6 biosynthetic loci in the budding yeast Saccharomyces cerevisiae. J Heredity 101(suppl_1):S75-S84.",
            "doi": "10.1093/jhered/esq007",
            "data_source": "https://github.com/harmslab/notebooks-nonlinear-high-order-epistasis/tree/master/datasets",
            "organism": "Saccharomyces cerevisiae",
            "strain_background": "BY series (Brachmann et al. 1998, derived from S288C)",
            "n_loci": 6,
            "n_genotypes": 64,
            "locus_names": LOCUS_NAMES,
            "locus_name_note": (
                "Positions 0-5 correspond to 6 of the 7 standard yeast auxotrophic "
                "biosynthetic markers (ADE2, HIS3, LEU2, LYS2, MET15, TRP1, URA3). "
                "The exact positional mapping requires the original paper (paywalled). "
                "The harmslab repository data uses only positional indices."
            ),
            "candidate_genes": CANDIDATE_GENES,
            "genotype_encoding": "Binary string: '1' = wildtype (functional gene), '0' = knockout (gene deleted)",
            "wildtype_genotype": "111111",
            "fitness_components": {
                "haploid_growth_rate": {
                    "description": "Growth rate of haploid strains",
                    "n_replicates": 10,
                    "log_transformed": True,
                },
                "diploid_growth_rate": {
                    "description": "Growth rate of homozygous diploid strains",
                    "n_replicates": 10,
                    "log_transformed": True,
                },
                "mating_efficiency": {
                    "description": "Haploid mating efficiency",
                    "n_replicates": 3,
                    "log_transformed": True,
                },
                "sporulation_efficiency": {
                    "description": "Diploid sporulation efficiency",
                    "n_replicates": 3,
                    "log_transformed": True,
                    "note": "log(1+fitness) was used since some raw values were zero",
                },
            },
        },
        "genotypes": records,
    }


def build_pathway_interactions():
    """Build the pathway interaction data based on yeast biochemistry."""
    return {
        "metadata": {
            "source": "Standard yeast biochemistry (SGD, KEGG)",
            "note": (
                "The 6 loci studied are biosynthetic genes for amino acids, "
                "nucleotides, or vitamins. Each catalyzes a step in a distinct "
                "biosynthetic pathway, but the pathways share common precursors "
                "from central carbon metabolism. The specific pathway interactions "
                "below are based on the 7 candidate genes; the actual 6 used in "
                "the study are a subset of these."
            ),
        },
        "candidate_genes": {
            "ADE2": {
                "full_name": "phosphoribosylaminoimidazole carboxylase",
                "pathway": "de novo purine biosynthesis",
                "product": "adenine (AMP)",
                "precursor_from_central_metabolism": "PRPP (from pentose phosphate pathway)",
                "sgd_id": "S000005654",
            },
            "HIS3": {
                "full_name": "imidazoleglycerol-phosphate dehydratase",
                "pathway": "histidine biosynthesis",
                "product": "histidine",
                "precursor_from_central_metabolism": "PRPP and ATP",
                "sgd_id": "S000005728",
            },
            "LEU2": {
                "full_name": "beta-isopropylmalate dehydrogenase",
                "pathway": "leucine biosynthesis (branched-chain amino acid)",
                "product": "leucine",
                "precursor_from_central_metabolism": "pyruvate (from glycolysis)",
                "sgd_id": "S000000523",
            },
            "LYS2": {
                "full_name": "alpha-aminoadipate reductase",
                "pathway": "lysine biosynthesis (aminoadipate pathway)",
                "product": "lysine",
                "precursor_from_central_metabolism": "alpha-ketoglutarate (from TCA cycle) and acetyl-CoA",
                "sgd_id": "S000000714",
            },
            "MET15": {
                "full_name": "O-acetylhomoserine O-acetylserine sulfhydrylase",
                "pathway": "methionine/cysteine biosynthesis (sulfur amino acid)",
                "product": "methionine",
                "precursor_from_central_metabolism": "oxaloacetate (from TCA cycle, via aspartate)",
                "sgd_id": "S000004294",
            },
            "TRP1": {
                "full_name": "phosphoribosylanthranilate isomerase",
                "pathway": "tryptophan biosynthesis (aromatic amino acid)",
                "product": "tryptophan",
                "precursor_from_central_metabolism": "chorismate (from shikimate pathway, uses PEP + E4P)",
                "sgd_id": "S000002414",
            },
            "URA3": {
                "full_name": "orotidine-5'-phosphate decarboxylase",
                "pathway": "de novo pyrimidine biosynthesis",
                "product": "uracil (UMP)",
                "precursor_from_central_metabolism": "aspartate (from TCA cycle) and carbamoyl phosphate",
                "sgd_id": "S000000747",
            },
        },
        "pathway_connections": {
            "shared_precursor_PRPP": {
                "genes": ["ADE2", "HIS3"],
                "description": (
                    "Both purine (ADE2) and histidine (HIS3) biosynthesis use "
                    "PRPP (5-phosphoribosyl-1-pyrophosphate) as a key precursor. "
                    "PRPP is produced from ribose-5-phosphate by the pentose "
                    "phosphate pathway."
                ),
                "interaction_type": "shared_precursor",
            },
            "shared_precursor_aspartate_oxaloacetate": {
                "genes": ["MET15", "URA3"],
                "description": (
                    "Methionine (MET15) and pyrimidine (URA3) biosynthesis both "
                    "draw on aspartate/oxaloacetate from the TCA cycle."
                ),
                "interaction_type": "shared_precursor",
            },
            "shared_precursor_TCA_cycle": {
                "genes": ["LYS2", "MET15", "URA3"],
                "description": (
                    "Lysine (LYS2, via alpha-ketoglutarate), methionine (MET15, "
                    "via oxaloacetate/aspartate), and pyrimidine (URA3, via "
                    "aspartate) biosynthesis all draw precursors from the TCA cycle."
                ),
                "interaction_type": "shared_precursor",
            },
            "purine_histidine_link": {
                "genes": ["ADE2", "HIS3"],
                "description": (
                    "Histidine biosynthesis produces AICAR "
                    "(5-aminoimidazole-4-carboxamide ribotide), which feeds into "
                    "the purine biosynthetic pathway. This creates a direct "
                    "metabolic link between HIS3 and ADE2 pathways."
                ),
                "interaction_type": "metabolite_feed",
            },
            "all_pathways_share_central_metabolism": {
                "genes": ["ADE2", "HIS3", "LEU2", "LYS2", "MET15", "TRP1", "URA3"],
                "description": (
                    "All biosynthetic pathways ultimately draw precursors from "
                    "central carbon metabolism (glycolysis, pentose phosphate "
                    "pathway, TCA cycle). Competition for these shared metabolic "
                    "resources creates potential for indirect epistatic interactions "
                    "among all loci, though the strength depends on metabolic flux."
                ),
                "interaction_type": "indirect_competition",
            },
        },
        "pairwise_pathway_relationships": [
            {"gene1": "ADE2", "gene2": "HIS3", "same_pathway": False, "shared_precursor": True, "metabolite_link": True, "notes": "Both use PRPP; HIS produces AICAR for purine pathway"},
            {"gene1": "ADE2", "gene2": "LEU2", "same_pathway": False, "shared_precursor": False, "metabolite_link": False, "notes": "Different biosynthetic branches"},
            {"gene1": "ADE2", "gene2": "LYS2", "same_pathway": False, "shared_precursor": False, "metabolite_link": False, "notes": "Different biosynthetic branches"},
            {"gene1": "ADE2", "gene2": "MET15", "same_pathway": False, "shared_precursor": False, "metabolite_link": False, "notes": "Different biosynthetic branches"},
            {"gene1": "ADE2", "gene2": "TRP1", "same_pathway": False, "shared_precursor": False, "metabolite_link": False, "notes": "Different biosynthetic branches"},
            {"gene1": "ADE2", "gene2": "URA3", "same_pathway": False, "shared_precursor": False, "metabolite_link": True, "notes": "Purine and pyrimidine biosynthesis share some regulatory connections"},
            {"gene1": "HIS3", "gene2": "LEU2", "same_pathway": False, "shared_precursor": False, "metabolite_link": False, "notes": "Different biosynthetic branches"},
            {"gene1": "HIS3", "gene2": "LYS2", "same_pathway": False, "shared_precursor": False, "metabolite_link": False, "notes": "Different biosynthetic branches"},
            {"gene1": "HIS3", "gene2": "MET15", "same_pathway": False, "shared_precursor": False, "metabolite_link": False, "notes": "Different biosynthetic branches"},
            {"gene1": "HIS3", "gene2": "TRP1", "same_pathway": False, "shared_precursor": False, "metabolite_link": False, "notes": "Different biosynthetic branches"},
            {"gene1": "HIS3", "gene2": "URA3", "same_pathway": False, "shared_precursor": False, "metabolite_link": False, "notes": "Different biosynthetic branches"},
            {"gene1": "LEU2", "gene2": "LYS2", "same_pathway": False, "shared_precursor": False, "metabolite_link": False, "notes": "Different biosynthetic branches"},
            {"gene1": "LEU2", "gene2": "MET15", "same_pathway": False, "shared_precursor": False, "metabolite_link": False, "notes": "Different biosynthetic branches"},
            {"gene1": "LEU2", "gene2": "TRP1", "same_pathway": False, "shared_precursor": False, "metabolite_link": False, "notes": "Different biosynthetic branches"},
            {"gene1": "LEU2", "gene2": "URA3", "same_pathway": False, "shared_precursor": False, "metabolite_link": False, "notes": "Different biosynthetic branches"},
            {"gene1": "LYS2", "gene2": "MET15", "same_pathway": False, "shared_precursor": True, "metabolite_link": False, "notes": "Both draw from TCA cycle intermediates"},
            {"gene1": "LYS2", "gene2": "TRP1", "same_pathway": False, "shared_precursor": False, "metabolite_link": False, "notes": "Different biosynthetic branches"},
            {"gene1": "LYS2", "gene2": "URA3", "same_pathway": False, "shared_precursor": True, "metabolite_link": False, "notes": "Both draw from TCA cycle intermediates"},
            {"gene1": "MET15", "gene2": "TRP1", "same_pathway": False, "shared_precursor": False, "metabolite_link": False, "notes": "Different biosynthetic branches"},
            {"gene1": "MET15", "gene2": "URA3", "same_pathway": False, "shared_precursor": True, "metabolite_link": False, "notes": "Both use aspartate/oxaloacetate from TCA cycle"},
            {"gene1": "TRP1", "gene2": "URA3", "same_pathway": False, "shared_precursor": False, "metabolite_link": False, "notes": "Different biosynthetic branches"},
        ],
    }


def main():
    outdir = os.path.dirname(os.path.abspath(__file__))

    fitness_data = build_fitness_data()
    fitness_path = os.path.join(outdir, "fitness_data.json")
    with open(fitness_path, "w") as f:
        json.dump(fitness_data, f, indent=2)
    print(f"Wrote {fitness_path} ({len(fitness_data['genotypes'])} genotypes)")

    pathway_data = build_pathway_interactions()
    pathway_path = os.path.join(outdir, "pathway_interactions.json")
    with open(pathway_path, "w") as f:
        json.dump(pathway_data, f, indent=2)
    print(f"Wrote {pathway_path}")


if __name__ == "__main__":
    main()
