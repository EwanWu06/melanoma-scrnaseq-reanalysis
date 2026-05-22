"""Q3.3 — progressive smoke test for the GitHub-master build of scGen.

Tests scgen 2.1.1 (installed from git+https://github.com/theislab/scgen.git)
against the env's scvi-tools 1.4.2, on tiny fake data. Each step prints a
marker before it runs, so a crash shows exactly which step failed. No
fixes are attempted — the point is to learn how far the GitHub build gets.

Touches no project data. Safe to delete.
"""
import warnings

warnings.filterwarnings("ignore")

# --- Step 2a: imports -------------------------------------------------
print("--- Step 2a: import scgen + scvi ---")
import importlib.metadata as md

import scgen
import scvi

print("  scvi-tools:", scvi.__version__)
print("  scgen     :", md.version("scgen"))
print("PASSED 2a: import OK")

# --- Step 2b: build a SCGEN object on fake data -----------------------
print("--- Step 2b: build SCGEN object ---")
import anndata
import numpy as np
import scanpy as sc  # noqa: F401  (imported to confirm it still co-exists)

n_cells, n_genes = 100, 50
adata = anndata.AnnData(
    X=np.random.poisson(1, (n_cells, n_genes)).astype(np.float32),
    obs={"batch": ["A"] * 50 + ["B"] * 50,
         "celltype": ["x"] * 30 + ["y"] * 70},
    var={"gene_id": [f"g{i}" for i in range(n_genes)]},
)
adata.obs["batch"] = adata.obs["batch"].astype("category")
adata.obs["celltype"] = adata.obs["celltype"].astype("category")

scgen.SCGEN.setup_anndata(adata, batch_key="batch", labels_key="celltype")
model = scgen.SCGEN(adata)
print("PASSED 2b: SCGEN object created")

# --- Step 2c: train 1 epoch ------------------------------------------
print("--- Step 2c: train 1 epoch (CPU) ---")
model.train(max_epochs=1, accelerator="cpu")
print("PASSED 2c: 1 epoch training completed")

# --- Step 2d: get_latent_representation() ----------------------------
print("--- Step 2d: get_latent_representation() ---")
latent = model.get_latent_representation()
print(f"PASSED 2d: get_latent_representation() -> shape {latent.shape}")

# --- Step 2e: history loss keys --------------------------------------
print("--- Step 2e: model.history keys ---")
keys = list(model.history.keys())
print("  model.history keys:", keys)
print("PASSED 2e")

print("=== ALL STEPS PASSED ===")
