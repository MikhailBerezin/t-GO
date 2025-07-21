# -----------------------------------------------
# 0 | Imports, constants, output folders
# -----------------------------------------------
import os, numpy as np, pandas as pd, matplotlib.pyplot as plt, seaborn as sns
from gprofiler import GProfiler

gp = GProfiler(return_dataframe=True)

GENE_FILES = [
              "1790 genes heart FC 1.5.txt"         ]

#TOP_N_TERMS = 700
THEMES = {
    "Stress & cytokine response": [
        "stress", "interferon", "cytokine", "inflammatory", "defense"
    ],
    "Inflammation & immune signaling": [
        "inflammation", "inflammatory", "tnf", "il-1", "il-6", "nf-kb", "toll-like",
        "interleukin", "chemokine", "ccl", "cxcl", "immune response",
        "inflammasome", "pattern recognition", "pathogen response"
    ],
    "Oxidative stress & redox regulation": [
        "oxidative", "redox", "reactive oxygen", "ros", "nitrosative", "nrf2",
        "antioxidant", "glutathione", "superoxide", "peroxidase", "peroxiredoxin",
        "sod", "catalase", "thioredoxin", "oxidoreductase"
    ],
    "Extracellular matrix & adhesion": [
        "extracellular", "matrix", "adhesion", "integrin", "collagen",
        "remodeling", "fibronectin", "laminin", "basement membrane",
        "mmp", "matrix metalloproteinase", "tenascin", "focal adhesion",
        "ecm", "tissue remodeling", "stromal", "scaffold", "matrisome",
        "cell junction", "cell adhesion", "cell-matrix", "desmosome"
    ],
    "Metabolic re-wiring": [
        "metabolic", "oxidoreductase", "catabolic", "fatty",
        "one-carbon", "biosynthetic", "glycolysis", "glyco"
    ],
    "Hematopoietic & immune commitment": [
        "hematopoiet", "myeloid", "lymphoid", "leukocyte", "granulocyte",
        "erythro", "megakary", "erythropoiet", "myelopoiet", "thrombopoiet",
        "lymphocyte", "monocyte", "neutrophil", "eosinophil", "basophil",
        "platelet", "erythrocyte", "anemia", "cytopenia", "pancytopenia",
        "thrombocytopenia", "leukopenia", "neutropenia", "immune cell",
        "blood cell", "hematologic", "hematopoiesis", "stem cell", "hsc"
    ],
    "Cell-cycle & Apoptosis": [
        "cell cycle", "mitotic", "chromosome", "checkpoint",
        "dna replication", "nuclear division", "apoptosis",
        "programmed cell death", "caspase"
    ],
    "Neuronal Excitability & Synapse": [
        "axon", "dendrite", "synapse", "neurotransmitter", "vesicle",
        "action potential", "ion channel", "potassium", "sodium", "calcium",
        "glutamate", "gaba", "synaptic", "neurogenesis", "axonogenesis"
    ],
    "Neurotrophic Signaling & Growth Factors": [
        "neurotrophin", "ngf", "bdnf", "ntf", "trk", "trka", "trkb", "gdnf",
        "growth factor", "igf", "egf", "fgf", "receptor tyrosine kinase"
    ],
    "Immune-Neuronal Crosstalk": [
        "microglia", "macrophage", "satellite glia", "neuroimmune", "neuroinflammation",
        "cd11b", "cd68", "csf1", "tslp", "complement", "ccr", "cxcr"
    ],
    "Pain & Nociception": [
        "pain", "nociception", "nociceptor", "hyperalgesia", "allodynia",
        "trpv1", "trpa1", "scn9a", "piezo", "itch", "sensory perception", "neuropeptide"
    ],
    "Oxidative Phosphorylation & Mitochondria": [
        "mitochondrial", "oxidative phosphorylation", "electron transport chain",
        "atp synthase", "complex i", "respiratory chain", "mitophagy"
    ],
    "Autophagy & Proteostasis": [
        "autophagy", "lysosome", "proteasome", "ubiquitin", "protein folding", "chaperone"
    ],
    "Myelination & Schwann Cell Biology": [
        "myelin", "schwann cell", "mbp", "mpz", "prx", "pmp22", "node of ranvier"
    ],
"Cardiac & Muscle Function": [
    "heart", "cardiac", "cardiomyocyte", "myocardium", "ventricle",
    "atrium", "electrophysiology", "contraction", "contractile", "sarcomere",
    "troponin", "myosin", "actin", "titin", "desmin", "ryr2", "calmodulin",
    "cardiovascular", "vasculature", "blood pressure", "stroke volume",
    "ejection fraction", "arrhythmia", "bradycardia", "tachycardia",
    "heart failure", "hypertrophy", "cardiomyopathy"
]
}



os.makedirs("go_theme_outputs", exist_ok=True)
os.makedirs("theme_correlations", exist_ok=True)


# -----------------------------------------------
# 1 | Helper functions (enrichment & theming)
# -----------------------------------------------
def load_genes(txt):
    genes = [g.strip() for g in open(txt) if g.strip()]
    print("Loaded", txt, "(", len(genes), "genes )")
    return genes


def enrich(genes, p_thresh=1e-2):
    df = gp.profile(organism="mmusculus", query=genes)
    df = df[df["p_value"] < p_thresh].sort_values("p_value").copy()
    df["Score"] = -np.log10(df["p_value"])
    print(f"Significant enriched terms (p < {p_thresh}): {len(df)} rows")
    return df
#Only works when TOP_N_TERMS is active and specified.
#def enrich(genes):
#    df = gp.profile(organism="mmusculus", query=genes)
#   df = df.sort_values("p_value").head(TOP_N_TERMS).copy()
#   df["Score"] = -np.log10(df["p_value"])
#    return df

def assign_theme(name):
    low = name.lower()
    for th, kws in THEMES.items():
        if any(kw in low for kw in kws):
            return th
    return None


def aggregate(df):
    df["Theme"] = df["name"].apply(assign_theme)
    themed = (df.dropna(subset=["Theme"])
              .groupby("Theme")
              .agg(Score=("Score", "sum"),
                   Terms=("Theme", "count"))
              .sort_values("Score", ascending=False))
    return themed


# -----------------------------------------------
# 2 | Enrichment plots + TSV tables (unchanged)
# -----------------------------------------------
for path in GENE_FILES:
    genes = load_genes(path)
    enr = enrich(genes)
    themed = aggregate(enr)

    # save TSV
    tsv_out = os.path.join("go_theme_outputs",
                           path.replace(' ', '_').replace('.txt', '_themes.tsv'))
    themed.to_csv(tsv_out, sep='\t')
    print("table saved →", tsv_out)

    # bar-plot
    plt.figure(figsize=(8, 4.5))
    plt.barh(themed.index, themed.Score, color=plt.cm.Set2.colors[:len(themed)])
    plt.gca().invert_yaxis()
    plt.xlabel("Cumulative –log$_{10}$(p$_{value}$)", size=11)
    plt.title("Thematic processes  (" + path + ")", loc="left", weight="bold")
    plt.tight_layout()
    png_out = os.path.join("go_theme_outputs",
                           path.replace(' ', '_').replace('.txt', '_themes.png'))
    plt.savefig(png_out, dpi=300)
    plt.show()
    print("plot saved  →", png_out, "\n")

# -----------------------------------------------
#