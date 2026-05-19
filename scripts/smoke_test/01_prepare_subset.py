"""
Q3.2 prep — Seurat smoke-test subset builder.

Goal: produce a tiny 100-cell x 100-gene slice of the Tirosh GSE72056 data
for the R/Seurat compatibility test in `02_seurat_data_slot_test.R`.

This is a *sanity check artifact*, not part of the main analysis. The full
Q3.2 RPCA notebook will read the full AnnData via anndata2ri or a fresh
load.

Outputs (gitignored):
  data/processed/smoke_test/subset_100x100.tsv  — gene-by-cell log-TPM matrix
  data/processed/smoke_test/subset_metadata.tsv — per-cell patient + class

Run:
  conda activate melanoma-scrnaseq
  python scripts/smoke_test/01_prepare_subset.py
"""

from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
RAW = REPO / "data" / "raw" / "GSE72056_melanoma_single_cell_revised_v2.txt.gz"
OUT_DIR = REPO / "data" / "processed" / "smoke_test"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def main() -> None:
    print(f"Reading {RAW.name} ...")
    # Per CLAUDE.md: header=0 makes the first row (Cell IDs) the column header.
    # Rows 0-2 of the DataFrame are then patient / malignant flag / cell-type
    # code, and gene rows start at DataFrame index 3.
    df = pd.read_csv(RAW, sep="\t", header=0)
    df = df.set_index(df.columns[0])  # first column is the row label ("tumor", "malignant(...)", "non-malignant(...)", then gene symbols)
    print(f"  Loaded shape (rows x cols): {df.shape}  (rows = 3 meta + genes; cols = cells)")

    # Split metadata rows from gene-expression rows.
    meta = df.iloc[:3]
    expr = df.iloc[3:]
    expr = expr[~expr.index.duplicated(keep="first")]  # MARCH1/MARCH2 dedup
    print(f"  Genes: {expr.shape[0]}  Cells: {expr.shape[1]}")

    # Deterministic subset: first 100 cells, top-100 most variable genes
    # across those cells. Variance-based selection just to avoid all-zero
    # genes; we are NOT doing real HVG selection here.
    sub_cells = expr.columns[:100]
    cell_expr = expr[sub_cells]
    top_genes = cell_expr.var(axis=1).sort_values(ascending=False).index[:100]
    sub = cell_expr.loc[top_genes]
    sub_meta = meta[sub_cells].T  # cells x 3 meta columns
    sub_meta.columns = ["tumor", "malignant_flag", "non_malignant_type"]

    expr_path = OUT_DIR / "subset_100x100.tsv"
    meta_path = OUT_DIR / "subset_metadata.tsv"

    sub.to_csv(expr_path, sep="\t")
    sub_meta.to_csv(meta_path, sep="\t")

    print(f"  Subset shape: {sub.shape}  (genes x cells)")
    print(f"  Value range: [{sub.values.min():.3f}, {sub.values.max():.3f}]"
          f"  (expected: log2(TPM/10+1), ~0 to ~14)")
    print(f"  Wrote: {expr_path.relative_to(REPO)}")
    print(f"  Wrote: {meta_path.relative_to(REPO)}")


if __name__ == "__main__":
    main()
