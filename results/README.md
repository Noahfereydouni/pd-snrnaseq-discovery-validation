# pd-snrnaseq-discovery-validation

Analysis code for **"Cell-Type-Resolved Transcriptomic Signatures of Parkinson's Disease:
Discovery in a Primary Cohort and Independent Cross-Cohort Validation."**

Noah Fereydouni, Nika Abravan

## Overview

A discovery-and-validation single-nucleus RNA-seq study across two independent, publicly
available human postmortem cohorts:

| | Accession | Source | Cohort |
|---|---|---|---|
| Discovery | GSE157783 | Smajić et al., 2022 (Brain) | Midbrain, 11 donors (5 PD, 6 control) |
| Validation | GSE184950 | Wang et al., 2024 (Sci Adv) | Substantia nigra, ~84,000 nuclei after subsetting |

A six-gene candidate panel (CHI3L1, NEAT1, MOBP, PBX1, EN2, LINGO1) was nominated in the
discovery cohort and tested for replication in the validation cohort, alongside TH.

## Repository structure

```
scripts/
  01_primary_cohort_QC_clustering_DE.Rmd
      Loading, doublet removal (scDblFinder), QC filtering, normalization, PCA,
      Harmony integration, clustering and cell-type annotation, and per-cell-type
      differential expression. Produces the DE_*.csv exports and the violin plots.
      → manuscript Table 1, Table 2, Figure 1, Figure 2

  02_primary_cohort_pseudobulk_DESeq2.Rmd
      Donor-level pseudobulk aggregation (AggregateExpression) and DESeq2 differential
      expression for microglia, astrocytes, and OPCs within the discovery cohort.
      → manuscript Section 3.4

  03_primary_cohort_pathway_enrichment.Rmd
      Gene Ontology (clusterProfiler::enrichGO) and GSEA (gseGO, seed = 1234) by cell
      type, consuming the DE_*.csv exports from script 01.
      → manuscript Figure 3, Figure 4

  04_validation_cohort_cross_cohort_analysis.Rmd
      Loading and subsetting of the validation object, dopaminergic-neuron differential
      expression, volcano plot, enrichR GO enrichment, and extraction of the six-gene
      panel + TH from per-cell-type FindMarkers results.
      → manuscript Table 3, Figure 5, Figure 6, Figure 7

results/
      Committed analysis outputs: per-cell-type nucleus-level DE for the primary cohort
      (DE_Microglia, DE_Astrocytes, DE_Oligodendrocytes, DE_Dopaminergic_Neurons, DE_OPCs)
      and the donor-level astrocyte pseudobulk DESeq2 result. See results/README.md for
      column definitions, which table or figure each file backs, and provenance notes.
```

## Provenance

These are the original R Markdown sources that produced the manuscript's numbers, included
unmodified. Local paths (`~/Desktop/PD_Project_Fresh/...`) are left in place deliberately so
the provenance chain is auditable — replace them with your own paths to rerun.

## Data availability

Both datasets are public via NCBI GEO: [GSE157783](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE157783),
[GSE184950](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE184950).

## Software

R 4.5.1. Key packages: Seurat 5.5.0, SeuratObject 5.4.0, harmony 2.0.5, scDblFinder 1.24.10,
SingleCellExperiment 1.32.0, DESeq2 1.50.2, clusterProfiler 4.18.4, enrichplot 1.30.5,
DOSE 4.4.0, org.Hs.eg.db 3.22.0, GO.db 3.22.0, enrichR, GEOquery 2.78.0, ggplot2 4.0.3,
ggrepel 0.9.8, dplyr 1.2.1, data.table 1.18.4, Matrix 1.7-5, patchwork 1.3.2.

## Code ↔ manuscript reconciliation

The manuscript text was checked line by line against these scripts. Items resolved during
that check are recorded below so the provenance of each reported number is traceable.

**Resolved**

1. **Sex correction.** Earlier drafts claimed the primary-cohort DE used logistic regression
   with donor sex as a latent variable. It does not — the `FindMarkers()` calls in
   `01_...Rmd` run Seurat's default Wilcoxon test. The GSE157783 object as distributed
   carries no sex annotation (metadata holds `orig.ident`, `Patient_ID`, `Condition`,
   `Doublet_Class` and clustering resolutions only), so the covariate was unavailable rather
   than omitted by choice. The manuscript now states this and treats it as a limitation.
   `04_...Rmd` (validation cohort) does correctly apply `test.use = "LR", latent.vars = "Sex"`.

2. **Cluster identities.** Methods now match the annotation dictionary in `01_...Rmd`:
   22 communities, dopaminergic neurons at cluster 19, oligodendrocytes across clusters
   0/2/3/7. Vascular, stromal and low-quality clusters were annotated and then removed;
   Figure 1 shows the filtered object.

3. **CHI3L1.** The published Table 2 row for CHI3L1 in microglia (+2.89, 66.6% vs 29.1%,
   padj 2.7e-130) was traced to the OPC file, which contains CHI3L1 at +2.96, 67.5% vs 29.3%,
   padj 1.5e-131. In microglia CHI3L1 is +1.35 at 1.9% vs 0.8% with padj = 1. CHI3L1 has been
   removed from the candidate panel, which is now five genes.

4. **Endothelial-cell subsetting.** Removal is an explicit `subset(idents = c(...))` whitelist,
   unrelated to scDblFinder. Because every DE test runs within a single annotated cell type,
   the subsetting cannot change any reported statistic — confirmed by the July 2026 validation
   run against the unsubset object, which reproduced the Table 3 values exactly.

5. **Pseudobulk file collision.** `01_...Rmd` wrote a Seurat `FindMarkers()` result to
   `OPCs_Pseudobulk_Results.csv`, the same path used by the DESeq2 script in `02_...Rmd`,
   so the nucleus-level export overwrote the donor-level one. That file is now named
   `DE_OPCs.csv`. The astrocyte pseudobulk (authoritative August run) returns no gene
   significant after correction across 25,503 tested genes, and Section 3.4 reports this as
   a negative result.

**Open**

6. **Dopaminergic nucleus count — export selected.** Two vintages of
   `DE_Dopaminergic_Neurons.csv` exist. The version committed in `results/` has ~49 PD and
   ~196 control dopaminergic nuclei and is the authoritative export for this manuscript; all
   dopaminergic statistics in Tables 1–3 come from it. A later export, produced after a QC
   change, contains only 5 PD and ~71 control nuclei — `pct.1` there takes only the values
   0, 0.2, 0.4, 0.6, 0.8, 1.0. In that version PBX1 and LINGO1 saturate at 5/5 and reverse
   direction with uncorrected p-values of 0.43 and 0.55, TH falls below the `min.pct`
   detection pre-filter and is absent from the file entirely, and all 44 significant genes
   are upregulated with none downregulated. Since TH is one of the markers used to define the
   dopaminergic cluster in the first place (Methods 2.2), its loss is evidence about the
   surviving nuclei rather than about TH. The later run is treated as over-filtering of the
   rarest population in the dataset and is not used. The specific QC thresholds that differ
   between the two runs remain to be documented.

7. **Donor-level pseudobulk for microglia and OPCs.** DESeq2 results for these two cell types
   have not been exported to distinct filenames. Section 3.4 currently reports astrocytes only.

8. **Figure provenance for the validation UMAP.** The panel in the manuscript is the
   pre-subsetting UMAP (title "Cellular Ecosystem of the Human Brain GSE184950", endothelial
   population present), which is what its caption describes. `04_...Rmd` renders a different,
   post-subset version. The script that produced the published panel was overwritten and is
   not in this repository.

## License

Code in this repository is released under the MIT License (see `LICENSE`).
The accompanying preprint is distributed under CC BY.
