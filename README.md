# Melanoma scRNA-seq Reanalysis

Reanalysis of single-cell RNA sequencing (scRNA-seq) data from melanoma to
characterize the tumor microenvironment, immune cell states, and
transcriptional programs associated with the disease.

## Project status

- [x] Stage 0: Setup & Foundation
- [x] Stage 1: Foundation Building
- [x] Stage 2: Classical Pipeline Replication
- [ ] Stage 3: Deep Learning Methods Application
- [ ] Stage 4: Synthesis & Deliverables

**Stage 2 complete:** the classical pipeline (QC → HVG → PCA → Harmony → Leiden)
re-analyzed the Tirosh malignant cells; only 2 of Tsoi's 4 states (Melanocytic,
Undifferentiated) were robustly recovered — full write-up in
[`docs/stage2_report.md`](docs/stage2_report.md).

## Environment setup

This project uses a reproducible [conda](https://conda-forge.org/) environment.

```bash
# 1. Install Miniforge if you don't have conda (https://github.com/conda-forge/miniforge)
# 2. Create the environment from the spec file
conda env create -f environment.yml

# 3. Activate it
conda activate melanoma-scrnaseq

# 4. Sanity check
python -c "import scanpy as sc; print(sc.__version__)"
```

## Analysis stack

Built on the standard Python single-cell ecosystem, pinned in
[`environment.yml`](environment.yml):

- **scanpy** / **anndata** — single-cell data structures, QC, clustering, DE
- **scvi-tools** — probabilistic deep-learning models (integration, label transfer)
- **scikit-learn**, **numpy**, **pandas**, **scipy**, **statsmodels** — numerics
- **matplotlib**, **seaborn** — plotting
- **jupyter** / **jupyterlab** — interactive analysis

## Directory structure

```
melanoma-scrnaseq-reanalysis/
├── data/
│   ├── raw/         # Original, immutable downloaded data (git-ignored)
│   └── processed/   # Cleaned / intermediate data (git-ignored)
├── notebooks/       # Jupyter notebooks for exploration & analysis
├── src/             # Reusable Python modules and pipeline code
├── results/
│   ├── figures/     # Generated plots
│   └── tables/      # Generated tables / summary stats
├── docs/            # Notes, methods, references
├── environment.yml  # Conda environment specification
└── README.md
```

Large data files are intentionally **not** tracked in git (see `.gitignore`).
Each data/results subfolder keeps a `.gitkeep` so the structure is preserved.

## License

Released under the [MIT License](LICENSE).
