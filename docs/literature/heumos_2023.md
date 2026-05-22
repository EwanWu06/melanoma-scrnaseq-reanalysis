# Heumos et al. 2023 — Reading Notes

> Author collaboration note: the English summary was structured from the
> author's Chinese reading notes with AI assistance; the pipeline
> interpretation and project-specific choices are the author's own.
> Read status: Sections 1-2 read; full paper used as reference.

## Citation

Heumos, L., Schaar, A. C., Lance, C., et al. (2023). "Best practices for
single-cell analysis across modalities." *Nature Reviews Genetics*, 24,
550–572. DOI: 10.1038/s41576-023-00586-w.

## What This Paper Is
A comprehensive review of best practices for scRNA-seq analysis, covering
the full pipeline from raw data to downstream interpretation.

## Pipeline Steps (For Our Project)

| Step | What It Does | What We Do for Tirosh Data |
|------|--------------|---------------------------|
| 1. QC | Filter low-quality cells | Likely skip — data already QC'd |
| 2. Normalization + HVG | Correct depth, select genes | Skip normalize/log (already done); HVG with `flavor="seurat"` |
| 3. Dim reduction + Batch | PCA + Harmony/scVI | PCA(30 PCs) + Harmony for classical; scVI for deep learning |
| 4. Clustering + Annotation | Leiden + marker-based | Leiden; annotate with Tsoi 2018 markers |
| 5. Downstream | DEG, trajectory, communication | Skip for core project |

## Key Concepts Learned
- **Preprocessing & QC.** Filter cells with an excessive mitochondrial-gene
  fraction (dying cells) and remove doublets, so only single, healthy,
  living cells enter downstream analysis.
- **Normalization & feature selection.** Normalization corrects systematic
  sequencing-depth differences between cells; selecting highly variable
  genes (HVGs) concentrates downstream computation on genes that carry
  genuine biological variation, lowering noise.
- **Dimensionality reduction is mathematical preparation for clustering**,
  not just for plotting. It addresses the "curse of dimensionality" —
  distance metrics between cells break down in very high-dimensional space,
  so clustering on raw gene space fails. The 2-D UMAP/t-SNE plot is a
  *visualization byproduct* of dimensionality reduction, not the substrate
  the clustering runs on. Batch correction at this step prevents cells from
  grouping by experimental batch instead of by biology.
- **Clustering & annotation.** Graph-based clustering (e.g. Leiden) on the
  reduced space partitions the cells; each cluster is then given a
  biological identity via marker genes or a reference dataset.
- **Downstream analysis.** With labelled clusters in hand, further analyses
  (differential expression, trajectory inference, cell-cell communication)
  probe specific biological mechanisms or clinical questions.

## Relevance to Our Project
- This is our **methodological reference** for the full analysis pipeline.
- Whenever we make a design choice (QC threshold, normalization, etc.),
  we will cite Heumos 2023 as the source of best-practice guidance.
