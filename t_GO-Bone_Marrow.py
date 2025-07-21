# -----------------------------------------------
# 0 | Imports, constants, output folders
# -----------------------------------------------
import os, numpy as np, pandas as pd, matplotlib.pyplot as plt, seaborn as sns
from gprofiler import GProfiler

gp = GProfiler(return_dataframe=True)

GENE_FILES = [
              "3690 genes.txt"         ]

#TOP_N_TERMS = 700
THEMES = {
    "Stress & cytokine response":
        ["stress", "interferon", "cytokine", "inflammatory", "defense"],
    "Extracellular matrix & adhesion":
        ["extracellular", "matrix", "adhesion", "integrin", "collagen"],
    "Metabolic re-wiring":
        ["metabolic", "oxidoreductase", "catabolic", "fatty",
         "one-carbon", "biosynthetic"],
    "Hematopoietic & immune commitment":
        ["hematopoiet", "myeloid", "lymphoid", "leukocyte", "granulocyte",
         "erythro", "megakary", "erythropoiet", "myelopoiet", "thrombopoiet",
         "lymphocyte", "monocyte", "neutrophil", "eosinophil", "basophil",
         "platelet", "erythrocyte", "anemia", "cytopenia", "pancytopenia",
         "thrombocytopenia", "leukopenia", "neutropenia", "immune cell",
         "blood cell", "hematologic", "hematopoiesis", "stem cell", "hsc"],
    "Cell-cycle & Apoptosis":
        ["cell cycle", "mitotic", "chromosome", "checkpoint",
         "dna replication", "nuclear division", "apoptosis",
         "programmed cell death", "caspase"]
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