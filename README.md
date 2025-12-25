# Thematic  Gene Ontology (ThematicGO)



ThematicGO is a transparent, user-customizable framework for reorganizing Gene Ontology (GO) enrichment results into biologically meaningful themes. It is designed to reduce redundancy in standard GO outputs while preserving statistical rigor, gene-level traceability, and interpretability. The method is particularly useful for transcriptomic studies where biologically relevant signals are distributed across many overlapping GO terms.



### Key Features:



Top-down thematic aggregation of GO terms using user-defined keyword sets

Transparent gene-level traceability via independent reconstruction of GO annotations from NCBI

Quantitative theme scoring based on cumulative enrichment significance

Full auditability through automatically generated summary tables and subterm plots

Flexible design supporting Biological Process (BP), Molecular Function (MF), Cellular Component (CC), or ALL



### Workflow Overview:



Input gene list

One gene symbol per line (case-insensitive)

GO enrichment analysis

Performed using g:Profiler

Adjusted p-values preferred when available

User-defined significance threshold



### Theme assignment:



Enriched GO terms are assigned to predefined biological themes using keyword matching

A GO term may map to multiple themes



### Theme aggregation:



Theme score = sum of −log₁₀(p) across all contributing GO terms

All predefined themes are evaluated

Themes without significant terms receive Score = 0 and Terms = 0



### Annotation reconstruction (transparency layer):



GO → gene mappings reconstructed from NCBI gene2go and gene\_info

Input genes intersected case-insensitively with annotated genes



### Visualization \& outputs:



Theme bar plots (only themes with non-zero scores)

Per-theme subterm plots (GO terms contributing to each theme)

Summary tables documenting GO terms, genes, and intersections



##### &nbsp;	Output Files:



For an input gene list named:

use this file as an example: DRG\_oxa\_271\_genes\_FC1.5.txt



##### &nbsp;	Main Outputs:



\*\_themes.tsv	Theme summary table (includes all predefined themes, zero scores allowed)

\*\_themes.png	Bar plot of themes with non-zero scores

\*\_GO\_terms\_with\_intersection\_genes.tsv	GO-term-level table with gene intersections

subterm\_barplots/\*\_subterms.png	Per-theme GO subterm plots

Tables include all themes, even those with no significant signal

Plots exclude zero-score themes for clarity





## Required packages:



Python ≥ 3.9

pip install gprofiler-official pandas numpy matplotlib seaborn openpyxl



Internet access required on first run to download NCBI annotation files:



gene2go.gz

gene\_info.gz

These files are cached locally after download.



### Configuration:



Key parameters are defined at the top of the script:

ORGANISM = "mmusculus"   # or "hsapiens"

GO\_ASPECT = "BP"        # BP, MF, CC, or ALL

P\_THRESH = 0.05

Themes and their keyword definitions are fully editable via the THEMES dictionary.





### Citation:



If you use ThematicGO, please cite the associated manuscript: Zhimu Wang, Leland C. Sudlow, Junwei Du, Mikhail Y. Berezin; thematicGO: A Keyword-Based Framework for Interpreting Gene Ontology Enrichment via Biological Themes (in preparation).



### Contact:



For questions, suggestions, or contributions, please contact the Mikhail Berezin (berezinm@wustl.edu) at Washington Unoversity in St. Louis, School of Medicine or open an issue in the repository.

