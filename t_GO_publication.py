from __future__ import annotations

from pathlib import Path
from typing import Dict, Iterable, List, Optional

import urllib.request

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

try:
    from gprofiler import GProfiler
except ImportError as e:
    raise ImportError(
        "Could not import GProfiler.\n\n"
        "Install the official g:Profiler client:\n"
        "  pip install gprofiler-official\n\n"
        "If you previously installed the wrong package 'gprofiler', remove it:\n"
        "  pip uninstall gprofiler\n"
    ) from e

GENE_FILES: List[str] = [
    "DRG_oxa_271_genes_FC1.5.txt",
]

OUT_DIR = Path("go_theme_outputs")
CORR_DIR = Path("theme_correlations")
SUBTERM_DIR = OUT_DIR / "subterm_barplots"

OUT_DIR.mkdir(parents=True, exist_ok=True)
CORR_DIR.mkdir(parents=True, exist_ok=True)
SUBTERM_DIR.mkdir(parents=True, exist_ok=True)

ORGANISM = "mmusculus"
# ORGANISM = "hsapiens"

PVAL_COL_PREFERRED = "p_value_adjusted"
P_THRESH = 5e-2

GO_ASPECT = "BP"   # "BP", "MF", "CC", or "ALL"

THEMES: Dict[str, List[str]] = {
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
        "sod", "catalase", "thioredoxin", "oxidoreductase",  "superoxide",
    "hydrogen peroxide",    "peroxide",    "nitric oxide",    "peroxynitrite",     "NADPH oxidase",
    "mitochondrial ROS",    "electron transport chain",    "mitochondrial dysfunction",
    # Redox damage & consequences
    "oxidative damage",
    "protein oxidation",
    "lipid peroxidation",
    "DNA oxidation",
    "redox imbalance"
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
        "one-carbon", "biosynthetic"
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
     "myelin", "schwann cell", "mbp", "mpz", "prx", "pmp22", "node of ranvier", "myelination",   "myelin sheath",   "myelin assembly",   "myelin maintenance",    "axon ensheathment",
    "axonal insulation",  "Schwann cell differentiation",  "Schwann cell proliferation",  "Schwann cell migration",
    "glial cell", "glial cell differentiation",   "peripheral glial cell",   "oligodendrocyte",
    "axon-glia interaction", "neurofilament organization","lipid biosynthesis", "cholesterol biosynthesis",
    "sphingolipid metabolism", "fatty acid metabolism","nerve development", "peripheral nervous system development",
    "axon development",    "axon guidance",    "nerve regeneration",    "remyelination",    "demyelination"
    ],
"Fibrosis": [
    "fibrosis","fibrotic","extracellular matrix", "matrix organization","matrix remodeling",
    "collagen", "collagen fibril","collagen biosynthesis","collagen organization",
    "fibronectin","laminin","proteoglycan", "elastin","fibroblast activation",
    "fibroblast proliferation", "myofibroblast", "myofibroblast differentiation",
    "tissue remodeling","wound healing","scar formation","transforming growth factor beta","TGF beta",
    "SMAD signaling","profibrotic signaling", "epithelial to mesenchymal transition",
    "EMT","endothelial to mesenchymal transition","EndMT", "lysyl oxidase", "matrix crosslinking",
    "tissue stiffness", "focal adhesion", "integrin signaling"
],
"Adipose Tissue Development": [
    "adipose tissue", "adipogenesis","adipocyte","adipocyte differentiation", "adipocyte development",
    "preadipocyte", "preadipocyte differentiation","fat cell differentiation","lipid droplet","lipid storage",
    "triglyceride metabolism","fatty acid uptake","fatty acid storage","lipogenesis","lipid biosynthetic process",
    "PPAR gamma", "C/EBP", "insulin signaling", "glucose uptake", "brown adipose tissue", "white adipose tissue",
    "beige adipocyte", "thermogenesis", "energy homeostasis", "metabolic regulation"
],
"Allergy": [
    "allergy",  "allergic",    "allergic response",    "hypersensitivity",
    "type I hypersensitivity",    "IgE",    "IgE-mediated",    "Fc epsilon receptor",    "FcεRI",
    "mast cell",    "mast cell activation",    "mast cell degranulation",    "basophil",
    "basophil activation",    "histamine",    "histamine release",    "eosinophil",    "eosinophil activation",
    "type 2 immune response",    "Th2",    "IL-4",    "IL-5",    "IL-13",    "cytokine-mediated signaling",
    "leukotriene",    "prostaglandin",    "inflammatory mediator release",    "immune hypersensitivity"
]
}

gp = GProfiler(return_dataframe=True)

def load_genes(path: str | Path) -> List[str]:
    """Load gene symbols from a text file (one per line)."""
    path = Path(path)
    genes = [line.strip() for line in path.read_text().splitlines() if line.strip()]
    print(f"Loaded {path.name}: {len(genes)} genes")
    return genes

def _pick_pval_column(df: pd.DataFrame) -> str:
    """Pick adjusted p-value column if present; otherwise fall back to raw p-value."""
    if PVAL_COL_PREFERRED in df.columns:
        return PVAL_COL_PREFERRED
    if "p_value" in df.columns:
        return "p_value"
    raise KeyError(
        "No p-value column found in g:Profiler results. "
        "Expected 'p_value' and/or 'p_value_adjusted'."
    )

def filter_go_aspect(enr_df: pd.DataFrame, go_aspect: str) -> pd.DataFrame:
    go_aspect = go_aspect.upper().strip()
    if go_aspect in ("ALL", ""):
        return enr_df

    source_map = {"BP": "GO:BP", "MF": "GO:MF", "CC": "GO:CC"}
    if go_aspect not in source_map:
        raise ValueError("go_aspect must be one of: 'BP', 'MF', 'CC', 'ALL'")

    return enr_df[enr_df["source"] == source_map[go_aspect]].copy()

def enrich(genes: Iterable[str], p_thresh: float = P_THRESH) -> pd.DataFrame:
    """Run g:Profiler enrichment and compute a Score = -log10(p).

    Uses adjusted p-values if available; otherwise uses raw p-values.
    """
    df = gp.profile(organism=ORGANISM, query=list(genes))
    if df is None or df.empty:
        return pd.DataFrame()

    pcol = _pick_pval_column(df)
    df = df.sort_values(pcol).copy()
    df = df[df[pcol] < p_thresh].copy()
    df["Score"] = -np.log10(df[pcol].astype(float))
    df["p_col_used"] = pcol
    print(f"Significant terms ({pcol} < {p_thresh}): {len(df)}")
    return df


def assign_themes(term_name: str) -> list[str]:
    low = term_name.lower()
    matched = [
        theme
        for theme, keywords in THEMES.items()
        if any(kw in low for kw in keywords)
    ]
    return matched

def aggregate_themes(enr_df: pd.DataFrame) -> pd.DataFrame:
    """
    Aggregate enrichment results by Theme.
    Includes ALL predefined themes, even if they have zero significant terms.
    """
    # Initialize all themes with zeros
    out = pd.DataFrame(
        {
            "Theme": list(THEMES.keys()),
            "Score": 0.0,
            "Terms": 0,
        }
    ).set_index("Theme")

    if enr_df.empty:
        return out

    df = enr_df.copy()
#    df["Themes"] = df["name"].apply(assign_themes)
#    df = df[df["Themes"].map(len) > 0]
#    df = df.explode("Themes").rename(columns={"Themes": "Theme"})

    df = enr_df.copy()
    df = df.dropna(subset=["Theme"])


    agg = (
        df.groupby("Theme", sort=False)
        .agg(
            Score=("Score", "sum"),
            Terms=("Theme", "count"),
        )
    )

    # Update initialized table
    out.update(agg)

    # Sort by score (optional)
    out = out.sort_values("Score", ascending=False)

    return out


def save_theme_table(themed: pd.DataFrame, prefix: str) -> Path:
    """Save theme summary table as TSV."""
    out = OUT_DIR / f"{prefix}_themes.tsv"
    themed.to_csv(out, sep="\t")
    return out

def plot_theme_bar(themed: pd.DataFrame, title: str, outfile: Path) -> None:
    """Publication-style horizontal bar plot for theme scores."""
    if themed.empty:
        return

    sns.set_context("paper", font_scale=1.0)
    plt.rcParams.update(
        {
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "font.size": 10,
            "axes.titlesize": 11,
            "axes.labelsize": 10,
        }
    )

    fig, ax = plt.subplots(figsize=(7.5, 4.2))
    ax.barh(themed.index[::-1], themed["Score"][::-1])
    ax.set_xlabel(r"Cumulative $-\log_{10}(p)$")
    ax.set_title(title, loc="left", weight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(outfile, dpi=600)
    plt.close(fig)

def plot_subterms_bar(enr_df: pd.DataFrame, theme_name: str, prefix: str) -> Optional[Path]:
    """Save a per-theme bar plot of subterms ranked by Score."""
  #  sub = enr_df.copy()
   # sub["Themes"] = sub["name"].apply(assign_themes)
 #   sub = sub[sub["Themes"].apply(lambda x: theme_name in x)]
#  sub = sub.explode("Themes")
 #   sub = sub[sub["Themes"] == theme_name]
  #  sub = sub.sort_values("Score", ascending=True)

    sub = enr_df[enr_df["Theme"] == theme_name].sort_values("Score", ascending=True)

    if sub.empty:
        return None

    sns.set_context("paper", font_scale=0.9)
    plt.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42, "font.size": 9})

    height = max(2.2, 0.28 * len(sub))
    fig, ax = plt.subplots(figsize=(7.5, height))
    ax.barh(sub["name"], sub["Score"])
    ax.set_xlabel(r"$-\log_{10}(p)$")
    ax.set_title(f"{theme_name}: enriched terms", loc="left", weight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()

    theme_safe = theme_name.replace(" ", "_").replace("&", "and").replace("/", "_")
    out = SUBTERM_DIR / f"{prefix}_{theme_safe}_subterms.png"
    fig.savefig(out, dpi=600)
    plt.close(fig)
    return out

def run_one_gene_file(path: str | Path) -> None:
    path = Path(path)
    prefix = path.stem.replace(" ", "_")

    genes = load_genes(path)
    enr = enrich(genes, p_thresh=P_THRESH)
    enr = filter_go_aspect(enr, GO_ASPECT)

    enr["All_Themes"] = enr["name"].apply(
        lambda x: "; ".join(assign_themes(x))
    )

    if enr.empty:
        print(f"No significant enrichment for {path.name}.")
        return

    # Export the supplementary table
    # Keep only GO terms (NCBI gene2go is GO-only; KEGG/Reactome cannot be expanded with gene2go)
    enr_go = enr[enr["native"].astype(str).str.startswith("GO:")].copy()

    if enr_go.empty:
        print("No GO terms (native starts with 'GO:') found in enrichment results; "
              "cannot build GO→genes intersections via NCBI gene2go.")
    else:

     export_supplementary_go_table_with_intersections(
        enr_df=enr,
        gene_file=path,
        prefix=f"{prefix}_{GO_ASPECT}",  # optional: keeps filenames distinct
        out_dir=OUT_DIR,
        categories=ncbi_categories_for_aspect(GO_ASPECT),
    )
    enr["Themes"] = enr["name"].apply(assign_themes)
    enr = enr.explode("Themes").rename(columns={"Themes": "Theme"})

    themed = aggregate_themes(enr)

    tsv_path = save_theme_table(themed, prefix)
    print(f"Theme table saved → {tsv_path}")

    png_path = OUT_DIR / f"{prefix}_themes.png"

    themed_nonzero = themed[themed["Score"] > 0]

    plot_theme_bar(
        themed_nonzero,
        title=f"Thematic processes ({path.name})",
        outfile=png_path
    )

    for theme in themed.index:
        plot_subterms_bar(enr, theme, prefix)


def main() -> None:
    for f in GENE_FILES:
        run_one_gene_file(f)

from pathlib import Path
import gzip
import urllib.request
import pandas as pd

NCBI_GENE2GO_URL = "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz"
NCBI_GENE_INFO_URL = "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz"
MOUSE_TAXID = 10090
# HUMAN_TAXID = 9606

def _download_if_missing(url: str, out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    if out_path.exists() and out_path.stat().st_size > 0:
        return
    print(f"Downloading → {out_path.name}")
    urllib.request.urlretrieve(url, out_path)

def _load_mouse_entrez_to_symbol(gene_info_gz: Path, taxid: int = MOUSE_TAXID) -> dict[int, str]:
    """
    Load mapping: Entrez GeneID -> Symbol (mouse only) from NCBI gene_info.gz
    """
    mapping: dict[int, str] = {}
    with gzip.open(gene_info_gz, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            if int(parts[0]) != taxid:
                continue
            gene_id = int(parts[1])
            symbol = parts[2]
            mapping[gene_id] = symbol
    return mapping

def ncbi_categories_for_aspect(go_aspect: str) -> set[str] | None:
    go_aspect = go_aspect.upper().strip()
    if go_aspect in ("ALL", ""):
        return None
    return {"BP": {"Process"}, "MF": {"Function"}, "CC": {"Component"}}[go_aspect]

def _build_goid_to_all_genes(
    gene2go_gz: Path,
    entrez_to_symbol: dict[int, str],
    taxid: int = MOUSE_TAXID,
    categories: set[str] | None = None,
) -> dict[str, set[str]]:
    """
    Build mapping: GO_ID -> set(Symbols) for all genes annotated to that GO term
    using NCBI gene2go.gz (filtered to the specified taxid and optional categories).

    Caches the final mapping to disk to avoid rebuilding on every run.
    """
    import gzip
    import pickle

    # Cache file name depends on taxid + categories to avoid mixing caches
    cat_tag = "ALL" if categories is None else "_".join(sorted(categories))
    cache_file = Path("ncbi_cache") / f"go_to_genes_tax{taxid}_{cat_tag}.pkl"
    cache_file.parent.mkdir(parents=True, exist_ok=True)

    if cache_file.exists() and cache_file.stat().st_size > 0:
        print(f"Loading cached GO→gene mapping → {cache_file}")
        with open(cache_file, "rb") as f:
            return pickle.load(f)

    go_to_genes: dict[str, set[str]] = {}

    with gzip.open(gene2go_gz, "rt", encoding="utf-8", errors="replace") as f:
        header = f.readline().rstrip("\n").split("\t")

        # Normalize header names (NCBI uses '#tax_id')
        header_norm = [h.lstrip("#") for h in header]
        col = {name: i for i, name in enumerate(header_norm)}

        required = ["tax_id", "GeneID", "GO_ID", "Category"]
        missing = [c for c in required if c not in col]
        if missing:
            raise KeyError(
                f"gene2go header missing columns: {missing}. Found (normalized): {header_norm}"
            )

        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= max(col.values()):
                continue

            if int(parts[col["tax_id"]]) != taxid:
                continue

            cat = parts[col["Category"]]
            if categories is not None and cat not in categories:
                continue

            entrez = int(parts[col["GeneID"]])
            symbol = entrez_to_symbol.get(entrez)
            if not symbol:
                continue

            go_id = parts[col["GO_ID"]]
            go_to_genes.setdefault(go_id, set()).add(symbol)

    # Cache once, after fully built
    with open(cache_file, "wb") as f:
        pickle.dump(go_to_genes, f)

    print(f"Cached GO→gene mapping → {cache_file}")
    return go_to_genes

def export_supplementary_go_table_with_intersections(
    enr_df: pd.DataFrame,
    gene_file: str | Path,
    prefix: str,
    out_dir: str | Path,
    cache_dir: str | Path = "ncbi_cache",
    categories: set[str] | None = None,   # e.g., {"Process"} to restrict to BP
) -> Path | None:
    """
    Export a GO-term supplementary table including:
      - full gene set size from NCBI annotation
      - intersection gene list with your input genes (case-insensitive intersection)

    Requires 'native' column in enr_df (GO IDs).
    """
    if enr_df is None or enr_df.empty:
        print("Enrichment dataframe is empty; skipping supplementary export.")
        return None

    if "native" not in enr_df.columns:
        raise KeyError(
            "enr_df must include a 'native' column containing GO IDs (e.g., GO:0006954)."
        )

    gene_file = Path(gene_file)

    # Read input genes (case-insensitive)
    in_genes_raw = [g.strip() for g in gene_file.read_text().splitlines() if g.strip()]
    in_genes_upper = {g.upper() for g in in_genes_raw}

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    cache_dir = Path(cache_dir)
    gene2go_gz = cache_dir / "gene2go.gz"
    gene_info_gz = cache_dir / "gene_info.gz"

    # Download annotation files once (cached)
    _download_if_missing(NCBI_GENE2GO_URL, gene2go_gz)
    _download_if_missing(NCBI_GENE_INFO_URL, gene_info_gz)

    # Build annotation maps
    print("Loading mouse GeneID→Symbol mapping...")
    entrez_to_symbol = _load_mouse_entrez_to_symbol(gene_info_gz, taxid=MOUSE_TAXID)

    print("Building GO_ID→all genes mapping (mouse)...")
    go_to_all_genes = _build_goid_to_all_genes(
        gene2go_gz,
        entrez_to_symbol,
        taxid=MOUSE_TAXID,
        categories=categories,
    )

    # Make sure Theme exists
    if "Theme" not in enr_df.columns:
        enr_df = enr_df.copy()
        enr_df["Themes"] = enr_df["name"].apply(assign_themes)
        enr_df = enr_df.explode("Themes").rename(columns={"Themes": "Theme"})

        if enr_df.empty:
            print("No enriched terms matched any theme keywords; skipping themed supplementary export.")
            return None

    # Choose p-value column used
    pcol = _pick_pval_column(enr_df)

    rows = []
    for _, r in enr_df.iterrows():
        go_id = r["native"]
        all_genes = go_to_all_genes.get(go_id, set())

        if not all_genes:
            inter = []
        else:
            # Case-insensitive intersection, but keep canonical symbols from NCBI
            all_upper_to_canonical: dict[str, str] = {}
            for g in all_genes:
                gu = g.upper()
                # keep the first seen canonical symbol (stable enough for output)
                if gu not in all_upper_to_canonical:
                    all_upper_to_canonical[gu] = g

            inter_upper = sorted(set(all_upper_to_canonical.keys()).intersection(in_genes_upper))
            inter = [all_upper_to_canonical[u] for u in inter_upper]

        rows.append({
            "GO_ID": go_id,
            "GO_term": r["name"],
            "source": r.get("source", np.nan),

            "Theme": r.get("Theme", None),
            "p_value_used": r[pcol],
            "Score": r["Score"],

            # Sizes
            "term_size_gProfiler": r.get("term_size", np.nan),
            "Term_size_annotation": len(all_genes),

            # Intersection
            "Intersection_size": len(inter),
            "Intersection_genes": ",".join(inter),
        })

    sup = pd.DataFrame(rows).sort_values(["Theme", "Score"], ascending=[True, False])

    out_path = out_dir / f"{prefix}_GO_terms_with_intersection_genes.tsv"
    sup.to_csv(out_path, sep="\t", index=False)

    print(f"Supplementary table saved → {out_path}")
    #return out_path

    sup = pd.DataFrame(rows).sort_values(["Theme", "Score"], ascending=[True, False])

    out_tsv = out_dir / f"{prefix}_GO_terms_with_intersection_genes.tsv"
    out_csv = out_dir / f"{prefix}_GO_terms_with_intersection_genes.csv"
    out_xlsx = out_dir / f"{prefix}_GO_terms_with_intersection_genes.xlsx"

    print("Saved files:")
    print(out_tsv.resolve())
    print(out_csv.resolve())
    print(out_xlsx.resolve())

    sup.to_csv(out_tsv, sep="\t", index=False)
    sup.to_csv(out_csv, index=False)

    # Excel export requires: pip install openpyxl
    sup.to_excel(out_xlsx, index=False)

    print(f"Supplementary tables saved → {out_tsv} | {out_csv} | {out_xlsx}")
    return out_tsv


if __name__ == "__main__":
    main()