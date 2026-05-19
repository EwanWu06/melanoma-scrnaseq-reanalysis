# Q3.2 prep — Seurat v5 smoke test
# -----------------------------------------------------------------------------
# Verifies the data-slot workflow proposed in the Q3.1 feasibility audit:
#
#   1. R + Seurat v5 are installed and importable
#   2. A log-TPM matrix can be written to obj[["RNA"]]$data directly,
#      bypassing NormalizeData()
#   3. The downstream Seurat pipeline (FindVariableFeatures -> ScaleData ->
#      RunPCA) runs without error on log-TPM input
#
# Run:
#   conda activate melanoma-r
#   Rscript scripts/smoke_test/02_seurat_data_slot_test.R
#
# This is a sanity check; not part of the Stage 3 main analysis.
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

# --- 0. Header ----------------------------------------------------------------
cat("====================================================================\n")
cat("Seurat v5 data-slot smoke test\n")
cat("--------------------------------------------------------------------\n")
cat("R version:      ", R.version.string, "\n")
cat("Seurat version: ", as.character(packageVersion("Seurat")), "\n")
cat("Matrix version: ", as.character(packageVersion("Matrix")), "\n")
cat("====================================================================\n\n")

results <- list()
record <- function(step, ok, detail = "") {
  results[[step]] <<- list(ok = ok, detail = detail)
  status <- if (ok) "PASS" else "FAIL"
  cat(sprintf("[%s] %s%s\n", status, step,
              if (nchar(detail) > 0) paste0(" -- ", detail) else ""))
}

# --- 1. Load the 100x100 subset -----------------------------------------------
expr_path <- "data/processed/smoke_test/subset_100x100.tsv"
meta_path <- "data/processed/smoke_test/subset_metadata.tsv"

if (!file.exists(expr_path)) {
  stop("Subset not found. Run scripts/smoke_test/01_prepare_subset.py first.")
}

expr_df <- read.table(expr_path, sep = "\t", header = TRUE,
                      row.names = 1, check.names = FALSE)
meta_df <- read.table(meta_path, sep = "\t", header = TRUE,
                      row.names = 1, check.names = FALSE)

expr_mat <- as.matrix(expr_df)
rng <- range(expr_mat)
record("load_subset",
       ok = identical(dim(expr_mat), c(100L, 100L)) && rng[1] >= 0 && rng[2] < 16,
       detail = sprintf("dim=%dx%d, value range=[%.3f, %.3f]",
                        nrow(expr_mat), ncol(expr_mat), rng[1], rng[2]))

# Convert to a sparse matrix (matches Seurat's expected storage).
expr_sparse <- as(expr_mat, "CsparseMatrix")

# --- 2. Build a Seurat object and write log-TPM into the data slot ------------
# Strategy: CreateSeuratObject populates the `counts` layer; we then
# overwrite the `data` layer with the same log-TPM matrix. Downstream
# Seurat functions (FindVariableFeatures, ScaleData, RunPCA) read from
# `data`, so this is the intended pre-normalized-input pattern.
obj <- tryCatch({
  CreateSeuratObject(counts = expr_sparse,
                     project = "smoke",
                     min.cells = 0, min.features = 0)
}, error = function(e) {
  record("create_object", FALSE, conditionMessage(e))
  NULL
})

if (is.null(obj)) {
  cat("\nFATAL: object creation failed; cannot continue.\n")
  quit(save = "no", status = 1)
}
record("create_object", TRUE,
       sprintf("nFeatures=%d, nCells=%d",
               nrow(obj), ncol(obj)))

# Attach the smoke-test metadata so the workflow has a batch identifier.
obj$tumor <- factor(meta_df[colnames(obj), "tumor"])

# Write log-TPM into the data slot.
write_data_ok <- tryCatch({
  obj[["RNA"]]$data <- expr_sparse
  TRUE
}, error = function(e) {
  record("write_data_slot", FALSE, conditionMessage(e))
  FALSE
})

if (write_data_ok) {
  data_layer <- LayerData(obj, layer = "data", assay = "RNA")
  data_rng <- range(data_layer)
  record("write_data_slot", TRUE,
         sprintf("data layer range=[%.3f, %.3f] (should match log-TPM scale)",
                 data_rng[1], data_rng[2]))
}

# --- 3. Downstream Seurat workflow on the data slot ---------------------------
# Skip NormalizeData() entirely — our values are already log-normalized.

ok_hvf <- tryCatch({
  obj <- FindVariableFeatures(obj, selection.method = "vst",
                              nfeatures = 50, verbose = FALSE)
  TRUE
}, error = function(e) {
  record("FindVariableFeatures", FALSE, conditionMessage(e))
  FALSE
})
if (ok_hvf) {
  n_hvg <- length(VariableFeatures(obj))
  record("FindVariableFeatures", TRUE,
         sprintf("identified %d variable features (requested 50)", n_hvg))
}

ok_scale <- tryCatch({
  obj <- ScaleData(obj, features = rownames(obj), verbose = FALSE)
  TRUE
}, error = function(e) {
  record("ScaleData", FALSE, conditionMessage(e))
  FALSE
})
if (ok_scale) {
  scale_layer <- LayerData(obj, layer = "scale.data", assay = "RNA")
  record("ScaleData", TRUE,
         sprintf("scaled matrix dim=%dx%d, mean~0 (got %.3f)",
                 nrow(scale_layer), ncol(scale_layer),
                 mean(scale_layer[1, ])))
}

ok_pca <- tryCatch({
  # npcs must be < min(features, cells); use 10 for a tiny 100x100 matrix.
  obj <- RunPCA(obj, npcs = 10, verbose = FALSE)
  TRUE
}, error = function(e) {
  record("RunPCA", FALSE, conditionMessage(e))
  FALSE
})
if (ok_pca) {
  pca_emb <- Embeddings(obj, reduction = "pca")
  record("RunPCA", TRUE,
         sprintf("PCA embedding dim=%dx%d",
                 nrow(pca_emb), ncol(pca_emb)))
}

# --- 4. Final verdict ---------------------------------------------------------
cat("\n====================================================================\n")
all_ok <- all(vapply(results, function(r) r$ok, logical(1)))
if (all_ok) {
  cat("OVERALL: PASS -- Q3.2 (Seurat RPCA on full Tirosh data) can proceed.\n")
  cat("\nNext step: notebooks/05_seurat_rpca.ipynb (R kernel) using the same\n")
  cat("data-slot pattern, but on the full malignant subset (1257 cells) with\n")
  cat("patient as the split-layer key.\n")
} else {
  cat("OVERALL: FAIL -- see failed steps above before proceeding to Q3.2.\n")
}
cat("====================================================================\n")

quit(save = "no", status = if (all_ok) 0 else 1)
