"""Q3.3 Phase 0 smoke test — verify scGen / scvi-tools API on the melanoma env.

Runs a 2-epoch scGen training on the real Tirosh data (full 4645 x 2000 HVG)
purely to exercise the API surface the Q3.3 notebook (06_scgen.ipynb) depends on.
Modifies nothing on disk. Safe to delete.
"""
import sys
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")
REPO = Path(__file__).resolve().parents[2]

print("=" * 64)
print("[1] VERSIONS  (torch AFTER install)")
import importlib.metadata as md
import torch
print("  torch     :", torch.__version__)
import scvi
print("  scvi-tools:", scvi.__version__)
import scgen
try:
    scgen_ver = scgen.__version__
except AttributeError:
    scgen_ver = md.version("scgen")
print("  scgen     :", scgen_ver)
import scanpy as sc
print("  scanpy    :", sc.__version__)

print("=" * 64)
print("[2] LOAD + HVG  (mirrors the real Phase 2 step)")
adata = sc.read_h5ad(REPO / "data/processed/tirosh_anndata_raw.h5ad")
print("  loaded shape        :", adata.shape)
print("  tumor_id dtype      :", adata.obs["tumor_id"].dtype)
print("  cell_type_label dt  :", adata.obs["cell_type_label"].dtype)
print("  obsm keys pre-HVG   :", list(adata.obsm.keys()))
sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat", subset=True)
print("  shape after HVG     :", adata.shape, " (expect (4645, 2000))")

print("=" * 64)
print("[3] scGen setup + 2-epoch train (CPU)")
scvi.settings.seed = 42
scgen.SCGEN.setup_anndata(adata, batch_key="tumor_id", labels_key="cell_type_label")
model = scgen.SCGEN(adata, n_latent=30, n_layers=1, n_hidden=128)
model.train(max_epochs=2, batch_size=32, accelerator="cpu")
print("  train OK")

print("=" * 64)
print("[4] API VERIFICATION")

# 4a — get_latent_representation()
try:
    latent = model.get_latent_representation()
    print("  get_latent_representation() : EXISTS")
    print("    -> type :", type(latent).__name__)
    print("    -> shape:", latent.shape, " (expect (4645, 30))")
except Exception as e:  # noqa: BLE001
    print("  get_latent_representation() : FAILED ->", repr(e))

# 4b — model.history loss keys
try:
    keys = list(model.history.keys())
    print("  model.history keys          :", keys)
    loss_like = [k for k in keys if "loss" in k.lower() or "elbo" in k.lower()]
    print("    -> loss-like keys         :", loss_like)
except Exception as e:  # noqa: BLE001
    print("  model.history               : FAILED ->", repr(e))

# 4c — obsm keys (does scGen auto-write any obsm?)
print("  adata.obsm keys post-train  :", list(adata.obsm.keys()))

print("=" * 64)
print("SMOKE TEST COMPLETE")
