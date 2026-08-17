# results/

Analysis outputs backing the numbers reported in the manuscript. All files are Seurat
`FindMarkers()` or DESeq2 exports produced by the scripts in `../scripts/`.

## Nucleus-level differential expression, primary cohort (GSE157783)

Columns: gene symbol (row name), `p_val`, `avg_log2FC`, `pct.1` (fraction of PD nuclei
expressing), `pct.2` (fraction of control nuclei expressing), `p_val_adj` (Bonferroni).
Note that `pct.1` and `pct.2` are **fractions, not percentages** — 0.484 means 48.4%.

| File | Cell type | Backs |
|---|---|---|
| `DE_Microglia.csv` | Microglia | Table 1, Figure 3A |
| `DE_Astrocytes.csv` | Astrocytes | Tables 1–2 (NEAT1), Figures 2C, 3C |
| `DE_Oligodendrocytes.csv` | Oligodendrocytes | Tables 1–2 (MOBP), Figures 2D, 3D |
| `DE_Dopaminergic_Neurons.csv` | Dopaminergic neurons | Tables 1–2 (PBX1, EN2, LINGO1, TH), Figures 2A–2B, 3B |
| `DE_OPCs.csv` | OPCs | not reported in the manuscript |

`DE_OPCs.csv` is included for completeness. It was previously written to a filename implying
donor-level pseudobulk results; it is a nucleus-level `FindMarkers()` export and should not be
read as pseudobulk.

## Donor-level pseudobulk, primary cohort (GSE157783)

`Astrocytes_Pseudobulk_DESeq2.csv` — DESeq2 output, columns `baseMean`, `log2FoldChange`,
`lfcSE`, `stat`, `pvalue`, `padj`, `Gene_Symbol`. Backs Section 3.4. No gene reaches
significance after Benjamini–Hochberg correction across 25,503 tested genes; the minimum
adjusted p-value in the file is 1.0. This is the reported result, not a truncated export.

Equivalent DESeq2 exports for microglia and OPCs are not yet available; see open item 7 in
the top-level README.

## Provenance note on the dopaminergic export

Two vintages of `DE_Dopaminergic_Neurons.csv` exist. The file committed here has roughly
49 PD and 196 control dopaminergic nuclei and is the authoritative export for this
manuscript: every dopaminergic statistic in Tables 1–3 comes from it. A later export
produced after a QC change contains only 5 PD nuclei and is not used; see open item 6 in
the top-level README for why.

## Not committed

Seurat objects (`.rds`), raw matrices, and GEO downloads are excluded by `.gitignore` on
account of size. Both source datasets are public: GSE157783 and GSE184950.
