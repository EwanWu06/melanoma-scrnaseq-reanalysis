"""Q3.3 — final test: does scGen 2.1.1's batch_removal() work vs scvi-tools 1.4.2?

Rebuilds the same tiny fake-data model as q33_scgen_github_smoke.py
(100 cells x 50 genes, 1-epoch CPU train) — the previous smoke-test process
has exited, so the model object is reconstructed identically here — then
calls batch_removal(). Touches no project data. Safe to delete.
"""
import warnings

warnings.filterwarnings("ignore")

import anndata
import numpy as np
import scgen

# --- rebuild the identical model from the previous smoke test ---------
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
model.train(max_epochs=1, accelerator="cpu")
print("--- model rebuilt + trained; calling batch_removal() ---")

# --- the actual test --------------------------------------------------
try:
    corrected = model.batch_removal()
    print("PASSED: batch_removal() executed")
    print(f"Returned type: {type(corrected)}")
    if hasattr(corrected, "shape"):
        print(f"  shape       : {corrected.shape}")
    if hasattr(corrected, "obsm"):
        print(f"  obsm keys   : {list(corrected.obsm.keys())}")
    if hasattr(corrected, "X") and corrected.X is not None:
        print(f"  .X shape    : {corrected.X.shape}")
    if hasattr(corrected, "obs"):
        print(f"  obs columns : {list(corrected.obs.columns)}")
except Exception as e:  # noqa: BLE001
    print(f"FAILED: {type(e).__name__}: {e}")
    import traceback

    traceback.print_exc()
