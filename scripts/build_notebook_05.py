"""
Notebook generator for Q3.2 — Seurat RPCA on Tirosh malignant cells.

Builds notebooks/05_seurat_rpca.ipynb with R kernel (`ir-melanoma`)
mirroring Stage 2 Q2.4's structure for direct comparability. The notebook
is generated programmatically so that re-generation is reproducible and
the cell list is easy to review without diffing raw JSON.

Run:
  conda activate melanoma-scrnaseq
  python scripts/build_notebook_05.py
"""

from pathlib import Path

import nbformat as nbf

REPO = Path(__file__).resolve().parents[1]
OUT = REPO / "notebooks" / "05_seurat_rpca.ipynb"

cells = []


def md(text: str) -> None:
    cells.append(nbf.v4.new_markdown_cell(text))


def code(text: str) -> None:
    cells.append(nbf.v4.new_code_cell(text))


# ============================================================================
# Header
# ============================================================================
md(r"""# Notebook 05 — Seurat RPCA Integration on Tirosh Malignant Cells

**Stage 3 / Q3.2** · R kernel (`ir-melanoma`) · 2026-05-19

## Goal

Re-integrate the Tirosh 2016 malignant subset (1,257 cells, 15 patients)
using **Seurat v5's reciprocal-PCA (RPCA) anchor-based integration**, in
parallel with Stage 2 Q2.4's Harmony pipeline. The output is a 30-dim
RPCA embedding + Leiden clusters at 5 resolutions + the same diagnostic
figure family that Stage 2 produced — enabling direct method comparison.

## Design rationale

- **Why RPCA?** Per `docs/stage3_feasibility_audit.md` §2, RPCA is the
  linear, anchor-based comparator to Harmony's linear post-PCA correction
  and scGen's non-linear VAE — three method families, one log-input
  surface, one evaluation set. This is one of three Core methods in
  `CLAUDE.md` Stage 3 Working Spec.

- **Why log-TPM directly into the `data` slot?** The audit's recommended
  pattern (skip `NormalizeData()`, write log-TPM into `obj[["RNA"]]$data`)
  was validated by the smoke test (`scripts/smoke_test/02_seurat_data_slot_test.R`)
  on a 100×100 subset. The full pipeline reuses the same pattern at full
  scale.

- **Evaluation scope:** malignant-only (1,257 cells), matching Stage 2 +
  Balderson 2024 for direct comparability. See `CLAUDE.md` Stage 3
  Working Spec.

- **Embedding dimension:** 30 dims (`npcs=50` for the initial PCA, then
  `dims=1:30` for RPCA integration and downstream), aligning with
  Stage 2's `n_pcs=30`.

- **Annotation:** *No automatic Tsoi-state assignment.* The notebook
  produces diagnostic plots and a per-cluster marker × patient-dominance
  summary table; the project author (Ewan) performs the manual
  annotation as in Q2.4.

## Method family note

Seurat RPCA does its correction in **low-dimensional PCA space** (anchor-
based), not in expression space. Output is a corrected cell embedding —
*not* corrected counts. Downstream UMAP / Leiden therefore operate on the
`integrated.rpca` reduction.
""")

# ============================================================================
# Phase 1 — Load and validate
# ============================================================================
md("## Phase 1 — Load and validate")

code(r"""# CWD guard - ensure consistent path resolution regardless of launch dir.
# jupyter nbconvert starts the kernel in the notebook's directory (notebooks/);
# jumping up one level matches the project-root-relative path convention used
# in the smoke test and all data/results references below.
if (basename(getwd()) == "notebooks") {
  setwd("..")
}
cat("Working directory:", getwd(), "\n")

# Library imports + version sanity check
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(SingleCellExperiment)
  library(zellkonverter)
  library(Matrix)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
})

cat("R               :", R.version.string, "\n")
cat("Seurat          :", as.character(packageVersion("Seurat")), "\n")
cat("SeuratObject    :", as.character(packageVersion("SeuratObject")), "\n")
cat("zellkonverter   :", as.character(packageVersion("zellkonverter")), "\n")
cat("SingleCellExperiment:", as.character(packageVersion("SingleCellExperiment")), "\n")
cat("Matrix          :", as.character(packageVersion("Matrix")), "\n")
cat("ggplot2         :", as.character(packageVersion("ggplot2")), "\n")
cat("patchwork       :", as.character(packageVersion("patchwork")), "\n")

# Set RNG seed for reproducibility (UMAP / Leiden have randomness)
set.seed(42)

# Ensure figure output dir exists
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
""")

code(r"""# Load AnnData -> SCE via zellkonverter
in_path <- "data/processed/tirosh_malignant_annotated.h5ad"
stopifnot(file.exists(in_path))

sce <- zellkonverter::readH5AD(in_path, X_name = "X", verbose = FALSE)
cat("SCE shape (genes x cells):", paste(dim(sce), collapse = " x "), "\n")
cat("Assays:", paste(assayNames(sce), collapse = ", "), "\n")
cat("colData columns:", paste(colnames(colData(sce)), collapse = ", "), "\n")
cat("reducedDims:", paste(reducedDimNames(sce), collapse = ", "), "\n")
""")

code(r"""# SCE -> Seurat (manual construction so we control which matrix lands in which slot)
expr <- assay(sce, "X")       # gene x cell, log2(TPM/10+1)
meta <- as.data.frame(colData(sce))

# Build Seurat object. We pass the log-TPM matrix as `counts` to satisfy the
# CreateSeuratObject signature, then immediately overwrite the data slot
# (and leave counts as-is). This is the smoke-test-validated pattern.
obj <- CreateSeuratObject(
  counts    = expr,
  meta.data = meta,
  project   = "tirosh_malignant_rpca",
  min.cells = 0,
  min.features = 0
)

# Write log-TPM into the data slot. Downstream Seurat functions read from `data`.
obj[["RNA"]]$data <- expr

cat("Seurat object built.\n")
cat("  dim:", paste(dim(obj), collapse = " x "), "\n")
cat("  assays:", paste(Assays(obj), collapse = ", "), "\n")
cat("  data layer range:", paste(round(range(LayerData(obj, layer='data')), 3), collapse = " - "), "\n")
""")

code(r"""# Verify Stage 2 metadata columns survived the conversion
expected <- c("tumor_id", "patient", "malignant_status", "cell_type_code",
              "tsoi_state",
              "leiden_r0.1", "leiden_r0.2", "leiden_r0.3",
              "leiden_r0.5", "leiden_r0.8", "leiden_r1.0")
have <- colnames(obj@meta.data)
missing <- setdiff(expected, have)

cat("meta.data columns (", length(have), "):\n", sep = "")
cat(paste("  ", have, collapse = "\n"), "\n")

if (length(missing) > 0) {
  stop("Missing expected meta.data columns: ", paste(missing, collapse = ", "),
       "\n  -> AnnData -> Seurat conversion lost metadata. Stop and fix before proceeding.")
} else {
  cat("\n[OK] All expected Stage 2 columns present.\n")
}

# Quick numeric summary
cat("\nshape:", paste(dim(obj), collapse = " x "), "(genes x cells)\n")
cat("patients:", length(unique(obj$patient)), "(expected 15)\n")
cat("tsoi_state levels:", paste(levels(factor(obj$tsoi_state)), collapse = ' | '), "\n")
""")

# ============================================================================
# Phase 2 — RPCA workflow
# ============================================================================
md(r"""## Phase 2 — RPCA integration

The split-then-integrate pattern: split the `RNA` assay layer by patient,
run the standard per-layer preprocessing (HVG / scale / PCA), then
`IntegrateLayers(method = RPCAIntegration)` to produce a corrected
embedding in `integrated.rpca`.

**Watch-out** (per `docs/stage3_feasibility_audit.md` §2.6 and the Stage 2
report): patients with very few cells (e.g. `pt75 = 3`, `pt65 = 4`)
may cause `FindVariableFeatures` to fail per-layer. We print the
distribution first so we can spot the risk; if `FindVariableFeatures`
errors, we stop and report instead of silently merging patients.""")

code(r"""# Patient cell-count distribution (small patients are a known Stage 2 limitation)
pt_counts <- sort(table(obj$patient))
cat("Cells per patient (smallest first):\n")
print(as.data.frame(pt_counts))

small <- pt_counts[pt_counts < 10]
if (length(small) > 0) {
  cat("\nWARNING: ", length(small), " patient(s) with <10 cells: ",
      paste(names(small), small, sep="=", collapse=", "), "\n",
      "  These may trip up per-layer FindVariableFeatures.\n", sep="")
}
""")

code(r"""# Split RNA assay layers by patient
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$patient)
cat("Layers after split:\n")
print(Layers(obj))
""")

code(r"""# Per-layer preprocessing: HVG -> Scale -> PCA
# Skip NormalizeData() — our data slot is already log-TPM.

obj <- FindVariableFeatures(obj, selection.method = "vst",
                            nfeatures = 2000, verbose = FALSE)
cat("[1/3] FindVariableFeatures done.\n")
cat("    VariableFeatures count (joint, after split):", length(VariableFeatures(obj)), "\n")

obj <- ScaleData(obj, verbose = FALSE)
cat("[2/3] ScaleData done.\n")

obj <- RunPCA(obj, npcs = 50, verbose = FALSE)
cat("[3/3] RunPCA done.\n")
cat("    PCA embedding dim:", paste(dim(Embeddings(obj, 'pca')), collapse = ' x '), "\n")
""")

code(r"""# RPCA integration in PCA space
# k.anchor = 5 (default). Per audit doc §2.3, can be raised to 20 if alignment is weak;
# we start at default for parity with the Stage 2 Harmony default-parameter run.

t0 <- Sys.time()
obj <- IntegrateLayers(
  object         = obj,
  method         = RPCAIntegration,
  orig.reduction = "pca",
  new.reduction  = "integrated.rpca",
  dims           = 1:30,
  k.anchor       = 5,
  verbose        = TRUE
)
cat("\nIntegration runtime:", round(as.numeric(Sys.time() - t0, units = 'secs'), 1), "sec\n")

stopifnot("integrated.rpca" %in% Reductions(obj))
cat("integrated.rpca dim:",
    paste(dim(Embeddings(obj, 'integrated.rpca')), collapse = ' x '), "\n")
""")

# ============================================================================
# Phase 3 — Downstream (join, neighbors, clusters, UMAP)
# ============================================================================
md("""## Phase 3 — Downstream: join layers + neighbors + clusters + UMAP

Joining the split layers back into a single assay before clustering /
UMAP, then standard `FindNeighbors` → `FindClusters` (5 resolutions,
matching Stage 2) → `RunUMAP`, all on the `integrated.rpca` reduction.""")

code(r"""# Re-join the split layers
obj <- JoinLayers(obj)
cat("After JoinLayers - layers:", paste(Layers(obj), collapse = ", "), "\n")
""")

code(r"""# Neighbors + multi-resolution Leiden + UMAP, all on integrated.rpca
obj <- FindNeighbors(obj, reduction = "integrated.rpca", dims = 1:30,
                     verbose = FALSE)
cat("[1/3] FindNeighbors done.\n")

obj <- FindClusters(obj,
                    resolution = c(0.1, 0.2, 0.3, 0.5, 0.8),
                    algorithm  = 4,    # 4 = Leiden (requires reticulate + leidenalg)
                    method     = "igraph",
                    cluster.name = paste0("rpca_leiden_r", c(0.1, 0.2, 0.3, 0.5, 0.8)),
                    verbose = FALSE)
cat("[2/3] FindClusters done. Resolutions stored as:\n")
res_cols <- grep("^rpca_leiden_r", colnames(obj@meta.data), value = TRUE)
for (rc in res_cols) {
  cat(sprintf("    %s : %d clusters\n", rc, nlevels(factor(obj@meta.data[[rc]]))))
}

obj <- RunUMAP(obj, reduction = "integrated.rpca", dims = 1:30,
               reduction.name = "umap.rpca", verbose = FALSE)
cat("[3/3] RunUMAP done.  umap.rpca dim:",
    paste(dim(Embeddings(obj, 'umap.rpca')), collapse = ' x '), "\n")
""")

# ============================================================================
# Phase 4 — Diagnostics (6 figures)
# ============================================================================
md("""## Phase 4 — Diagnostic figures (6)

All figures saved under `results/figures/` with `05_rpca_*` prefix,
parallel to Stage 2's `04_*` naming. Inline display follows each save.""")

code(r"""# Helper: a plain UMAP scatter on the RPCA embedding
umap_df <- as.data.frame(Embeddings(obj, "umap.rpca"))
colnames(umap_df) <- c("UMAP_1", "UMAP_2")
umap_df$patient    <- obj$patient
umap_df$tsoi_state <- obj$tsoi_state
for (rc in res_cols) umap_df[[rc]] <- factor(obj@meta.data[[rc]])
cat("UMAP coord df:\n")
str(umap_df)
""")

code(r"""# Figure 1: UMAP colored by patient
p1 <- ggplot(umap_df, aes(UMAP_1, UMAP_2, color = patient)) +
  geom_point(size = 0.6, alpha = 0.8) +
  theme_minimal(base_size = 11) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  labs(title = "Q3.2 RPCA UMAP — colored by patient (n=15)",
       subtitle = "Visual check of batch integration") +
  theme(legend.position = "right")

ggsave("results/figures/05_rpca_umap_patient.png", plot = p1,
       width = 9, height = 7, dpi = 150)
p1
""")

code(r"""# Figure 2: UMAP colored by Stage 2 tsoi_state
p2 <- ggplot(umap_df, aes(UMAP_1, UMAP_2, color = tsoi_state)) +
  geom_point(size = 0.7, alpha = 0.85) +
  scale_color_manual(values = c(
    "Melanocytic"      = "#1f77b4",
    "Undifferentiated" = "#d62728",
    "Other (batch?)"   = "#7f7f7f",
    "Transitory"       = "#ff7f0e",
    "Neural crest-like"= "#2ca02c"
  ), na.value = "lightgrey") +
  theme_minimal(base_size = 11) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  labs(title = "Q3.2 RPCA UMAP — colored by Stage 2 tsoi_state",
       subtitle = "Biological-signal-preservation visual check")

ggsave("results/figures/05_rpca_umap_stage2_tsoi.png", plot = p2,
       width = 9, height = 7, dpi = 150)
p2
""")

code(r"""# Figure 3: UMAP grid by Leiden cluster (one panel per resolution)
panels <- list()
for (rc in res_cols) {
  res_label <- sub("^rpca_leiden_", "", rc)
  pi <- ggplot(umap_df, aes(UMAP_1, UMAP_2, color = !!sym(rc))) +
    geom_point(size = 0.5, alpha = 0.85) +
    theme_minimal(base_size = 9) +
    labs(title = sprintf("Leiden %s (%d clusters)",
                         res_label, nlevels(umap_df[[rc]]))) +
    theme(legend.position = "right",
          legend.key.size = unit(0.3, "cm"),
          legend.text = element_text(size = 7))
  panels[[rc]] <- pi
}
combined <- wrap_plots(panels, ncol = 3)
ggsave("results/figures/05_rpca_umap_leiden_grid.png", plot = combined,
       width = 16, height = 10, dpi = 150)
combined
""")

code(r"""# Figure 4: Dotplot of MITF/SOX10/NGFR/AXL by cluster, one row per resolution
tsoi_markers <- c("MITF", "SOX10", "NGFR", "AXL")
have <- intersect(tsoi_markers, rownames(obj))
if (length(have) < length(tsoi_markers)) {
  cat("[note] missing Tsoi markers in this matrix:",
      paste(setdiff(tsoi_markers, have), collapse=", "), "\n")
}

dot_panels <- list()
for (rc in res_cols) {
  res_label <- sub("^rpca_leiden_", "", rc)
  Idents(obj) <- rc
  pi <- DotPlot(obj, features = have) +
    RotatedAxis() +
    labs(title = sprintf("RPCA %s", res_label),
         x = NULL, y = NULL) +
    theme(plot.title = element_text(size = 10),
          axis.text = element_text(size = 8),
          legend.position = "right")
  dot_panels[[rc]] <- pi
}
dot_combined <- wrap_plots(dot_panels, ncol = 1)
ggsave("results/figures/05_rpca_dotplot_by_res.png", plot = dot_combined,
       width = 7, height = 3 * length(dot_panels), dpi = 150)
dot_combined
""")

code(r"""# Figure 5: Cluster x patient crosstab heatmap, row-normalized, one panel per resolution
heat_panels <- list()
for (rc in res_cols) {
  res_label <- sub("^rpca_leiden_", "", rc)
  tab <- table(cluster = obj@meta.data[[rc]], patient = obj$patient)
  df  <- as.data.frame(tab) %>%
    group_by(cluster) %>%
    mutate(prop = Freq / sum(Freq)) %>%
    ungroup()
  pi <- ggplot(df, aes(patient, cluster, fill = prop)) +
    geom_tile() +
    scale_fill_viridis_c(limits = c(0, 1)) +
    theme_minimal(base_size = 9) +
    labs(title = sprintf("RPCA %s", res_label),
         fill = "row prop") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
          axis.text.y = element_text(size = 7),
          plot.title  = element_text(size = 10))
  heat_panels[[rc]] <- pi
}
heat_combined <- wrap_plots(heat_panels, ncol = 2)
ggsave("results/figures/05_rpca_crosstab_by_res.png", plot = heat_combined,
       width = 14, height = 3.5 * ceiling(length(heat_panels)/2), dpi = 150)
heat_combined
""")

code(r"""# Figure 6: Stage 2 Harmony UMAP vs Q3.2 RPCA UMAP — side by side, colored by tsoi_state
# Stage 2 X_umap was loaded via zellkonverter into reducedDims(sce, 'X_umap')
stage2_emb <- reducedDim(sce, "X_umap")
stopifnot(nrow(stage2_emb) == ncol(obj))

stage2_df <- data.frame(
  UMAP_1 = stage2_emb[, 1],
  UMAP_2 = stage2_emb[, 2],
  tsoi_state = obj$tsoi_state,
  source = "Stage 2 (Harmony)"
)
rpca_df_min <- data.frame(
  UMAP_1 = umap_df$UMAP_1,
  UMAP_2 = umap_df$UMAP_2,
  tsoi_state = umap_df$tsoi_state,
  source = "Q3.2 (Seurat RPCA)"
)
combo_df <- rbind(stage2_df, rpca_df_min)

p6 <- ggplot(combo_df, aes(UMAP_1, UMAP_2, color = tsoi_state)) +
  geom_point(size = 0.6, alpha = 0.85) +
  scale_color_manual(values = c(
    "Melanocytic"      = "#1f77b4",
    "Undifferentiated" = "#d62728",
    "Other (batch?)"   = "#7f7f7f",
    "Transitory"       = "#ff7f0e",
    "Neural crest-like"= "#2ca02c"
  ), na.value = "lightgrey") +
  facet_wrap(~ source, nrow = 1, scales = "free") +
  theme_minimal(base_size = 11) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  labs(title = "Stage 2 Harmony vs. Q3.2 Seurat RPCA — same cells, same tsoi colors",
       subtitle = "Side-by-side visual: do the methods agree on cell-state geometry?")

ggsave("results/figures/05_rpca_vs_harmony_sidebyside.png", plot = p6,
       width = 14, height = 6, dpi = 150)
p6
""")

# ============================================================================
# Phase 5 — Annotation prep + export
# ============================================================================
md("""## Phase 5 — Annotation prep + export

Per the Stage 2 / Q2.4 convention:
1. Build a per-cluster summary (Tsoi marker mean expression + dominant
   patient + size) at each resolution — handed to the project author for
   manual Tsoi-state annotation. **No automatic state assignment.**
2. Export the Seurat object (RPCA embedding + UMAP + new cluster columns)
   back into AnnData via `zellkonverter::writeH5AD`, preserving all
   Stage 2 metadata.""")

code(r"""# Per-cluster summary table for manual annotation
make_summary <- function(obj, res_col, markers) {
  cl <- obj@meta.data[[res_col]]
  cells_in <- split(colnames(obj), cl)
  expr_mat <- LayerData(obj, layer = "data")

  rows <- lapply(names(cells_in), function(k) {
    cells <- cells_in[[k]]
    pat_tab <- sort(table(obj$patient[cl == k]), decreasing = TRUE)
    top_pat <- names(pat_tab)[1]
    top_share <- round(100 * pat_tab[1] / sum(pat_tab), 1)
    marker_means <- sapply(markers, function(g) {
      if (g %in% rownames(expr_mat)) round(mean(expr_mat[g, cells]), 3) else NA
    })
    data.frame(
      resolution = sub("^rpca_leiden_", "", res_col),
      cluster = k,
      n_cells = length(cells),
      top_patient = top_pat,
      top_patient_pct = top_share,
      t(marker_means),
      check.names = FALSE
    )
  })
  do.call(rbind, rows)
}

markers <- c("MITF", "SOX10", "NGFR", "AXL")
summary_all <- do.call(rbind, lapply(res_cols, function(rc) make_summary(obj, rc, markers)))
summary_all$my_annotation <- ""    # left empty — for manual fill

cat("Per-cluster summary (manual annotation column empty by design):\n")
print(summary_all, row.names = FALSE)

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
write.csv(summary_all, "results/tables/05_rpca_cluster_summary.csv", row.names = FALSE)
cat("\nWritten: results/tables/05_rpca_cluster_summary.csv\n")
""")

code(r"""# Export Seurat -> SCE -> h5ad, preserving Stage 2 metadata + new RPCA artifacts
# Pre-export housekeeping: rename reductions so they end up with intuitive obsm keys

# Re-name the RPCA reduction to 'X_rpca' so zellkonverter writes obsm['X_rpca']
obj[["X_rpca"]] <- CreateDimReducObject(
  embeddings = Embeddings(obj, "integrated.rpca"),
  key = "RPCA_", assay = "RNA"
)
# Re-name the RPCA UMAP to 'X_umap_rpca'
obj[["X_umap_rpca"]] <- CreateDimReducObject(
  embeddings = Embeddings(obj, "umap.rpca"),
  key = "UMAPRPCA_", assay = "RNA"
)
# Drop the duplicates so the output isn't cluttered
obj[["integrated.rpca"]] <- NULL
obj[["umap.rpca"]] <- NULL
# 'pca' is the per-cell PCA used as RPCA input — keep but don't shadow Stage 2 X_pca

# Convert and write
sce_out <- as.SingleCellExperiment(obj, assay = "RNA")
cat("Output SCE assays:", paste(assayNames(sce_out), collapse=", "), "\n")
cat("Output SCE reducedDims:", paste(reducedDimNames(sce_out), collapse=", "), "\n")
cat("Output SCE colData cols (subset):",
    paste(head(colnames(colData(sce_out)), 20), collapse=", "), " ...\n")

out_path <- "data/processed/tirosh_malignant_rpca.h5ad"
zellkonverter::writeH5AD(sce_out, file = out_path, verbose = FALSE)
cat("\nWritten:", out_path, "\n")
""")

md(r"""## Annotation decisions (TO FILL)

> Per Q2.4's workflow: do **not** auto-assign Tsoi states. Review the
> diagnostic figures + the per-cluster summary table
> (`results/tables/05_rpca_cluster_summary.csv`) and fill in the
> `my_annotation` column with one of `Melanocytic` / `Transitory` /
> `Neural crest-like` / `Undifferentiated` / `Other (batch?)` per
> cluster, justified against the MITF / SOX10 / NGFR / AXL means.
>
> Then compare against Stage 2's Harmony assignment:
> - Where do RPCA and Harmony agree?
> - Where do they disagree, and which is more defensible?
> - Did RPCA recover any of the missing states (Transitory, Neural
>   crest-like) that Harmony missed?

**Chosen resolution:** *to fill — typically the lowest resolution that
cleanly separates the known states without over-fragmenting batch
artifacts.*

**Per-cluster annotation (resolution =  ):**

| cluster | tsoi_state | rationale |
|---|---|---|
|   |   |   |

**Cross-method agreement summary:**

*To fill after annotation.*
""")

# ============================================================================
# Build + write
# ============================================================================
nb = nbf.v4.new_notebook()
nb["cells"] = cells
nb["metadata"] = {
    "kernelspec": {
        "display_name": "R (melanoma-r)",
        "language": "R",
        "name": "ir-melanoma",
    },
    "language_info": {
        "codemirror_mode": "r",
        "file_extension": ".r",
        "mimetype": "text/x-r-source",
        "name": "R",
        "pygments_lexer": "r",
    },
}

OUT.parent.mkdir(parents=True, exist_ok=True)
with OUT.open("w", encoding="utf-8") as f:
    nbf.write(nb, f)

print(f"Wrote {OUT.relative_to(REPO)}  ({len(cells)} cells)")
