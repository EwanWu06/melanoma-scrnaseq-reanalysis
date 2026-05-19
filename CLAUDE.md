# CLAUDE.md — Project Context for Claude Code

> This file provides project-specific context for Claude Code.
> It is automatically loaded when Claude Code operates within this repository.
> Updated: 2026-05-17

---

## Project Overview

**Project name:** Melanoma scRNA-seq Re-analysis

**Working title (academic):** *Learning to Re-analyze: A Methodological Study of
Deep Learning Approaches on the Tirosh et al. (2016) Melanoma Single-Cell Dataset*

**One-line description:**
A learning-oriented re-analysis of a foundational 2016 melanoma scRNA-seq dataset
(Tirosh et al., GSE72056), applying modern deep learning methods (scVI, scANVI,
and potentially foundation models) and comparing results against a recent classical
re-analysis (Balderson et al., 2024).

**Purpose:**
This project is part of an undergraduate transfer application portfolio.
Target programs include UCLA Computational Biology, UCSD Public Health
(Biostatistics concentration), UNC Biostatistics, and UMich Informatics.
The project prioritizes genuine learning, methodological rigor, and
reproducibility over novel discovery.

---

## Current Stage

**Active stage:** Stage 2 — Classical Pipeline Replication (Q2.1–Q2.2 done; Q2.3 Harmony next)

**Stage progression:**
- [x] Stage 0: Setup & Foundation
  - [x] Repository structure created
  - [x] Tirosh 2016 data downloaded (GSE72056)
  - [x] Initial correspondence with Dr. Tirosh sent
  - [x] Core literature read (Tirosh 2016, Tsoi 2018, Balderson 2024, Heumos 2023)
- [x] Stage 1: Foundation Building
- [ ] Stage 2: Classical Pipeline Replication
  - [x] Q2.1: Data loading + AnnData + basic QC (notebook 01)
  - [x] Q2.2: HVG selection + PCA + batch-effect diagnostic (notebook 02)
  - [ ] Q2.3: Harmony batch correction
  - [ ] Q2.4+: Leiden clustering + Tsoi 4-state annotation
- [ ] Stage 3: Deep Learning Methods Application
- [ ] Stage 4: Synthesis & Deliverables
- [ ] Stretch: Cross-dataset Robustness (Jerby-Arnon GSE115978)

---

## Critical Data Constraints

**Primary dataset:** GSE72056 (Tirosh et al., 2016)
- File: `data/raw/GSE72056_melanoma_single_cell_revised_v2.txt.gz`
  (**gzip-compressed**; read directly with `pd.read_csv(..., sep="\t")` —
  pandas auto-detects `.gz`. Do not expect a plain `.txt`.)
- Shape: 4,645 cells × 23,686 genes
- Technology: Smart-seq2 (full-length, plate-based, non-UMI)
- 19 patients
- Mixed cell composition: malignant + immune + stromal/endothelial

**Data format (IMPORTANT):**
- Values are **log2(TPM/10 + 1)** — NOT raw counts
- Already normalized and log-transformed
- Min ≈ 0, Max ≈ 14 across the full matrix (≈11.11 in a 100×100 corner sample),
  all float64
- ~81% of values are zero (sparse)

**Implications for analysis:**
1. Do NOT apply `sc.pp.normalize_total()` or `sc.pp.log1p()` — data is already
   normalized and log-transformed.
2. scVI, scANVI, and most foundation models **cannot directly use this data** —
   they require raw integer counts (negative binomial-based likelihood models).
3. When selecting highly variable genes, use `sc.pp.highly_variable_genes(flavor="seurat")`
   — NOT `flavor="seurat_v3"` which requires raw counts.
4. Plan B (if raw counts cannot be obtained from Dr. Tirosh): use methods that
   accept log-normalized data (PCA, Harmony, classical Seurat workflow).

**Metadata structure (rows 1-4 of the file):**
- Row 1: Cell ID (encodes tumor + CD45 FACS sort + plate well)
- Row 2: tumor ID (patient number, integer)
- Row 3: malignant status (0=unresolved, 1=non-malignant, 2=malignant)
- Row 4: non-malignant cell type code (1=T, 2=B, 3=Macrophage, 4=Endothelial,
  5=CAF, 6=NK; set to 0 if malignant or unresolved)

> Loading note: with `pd.read_csv(header=0)` the `Cell` row becomes the column
> header, so only 3 metadata rows remain in the DataFrame (rows 0–2) and the
> first gene starts at DataFrame row 3. Also call `var_names_make_unique()` —
> gene symbols `MARCH1` and `MARCH2` are duplicated.

---

## Project Scope (Three Tiers)

**Core (must complete):**
Re-analyze Tirosh 2016 with 3-4 methods (classical baseline + scVI + scANVI),
compare against Balderson 2024's four-state result. Deliverables: GitHub repo,
technical report, dashboard.

**Stretch (if time allows):**
Extend analysis to Jerby-Arnon 2018 (GSE115978) for cross-dataset robustness
assessment.

**Reach (probably skip, documented as future work):**
Foundation model evaluation (scGPT, Geneformer, UCE), Smart-seq2 vs 10x
distribution shift analysis.

---

## Reference Framework

The project's biological "ground truth" for cell states is the **Tsoi et al. 2018
four-state model** (*Cancer Cell* 33:890-904):

| State | Key markers | Differentiation level |
|-------|-------------|----------------------|
| Melanocytic | MITF↑ SOX10↑ NGFR↓ AXL↓ | Most differentiated |
| Transitory | MITF↑ SOX10↑ NGFR↑ AXL↓ | Intermediate–high `[reconstructed — verify vs Tsoi 2018]` |
| Neural crest-like | MITF↓ SOX10↑ NGFR↑ AXL↑ | Low |
| Undifferentiated | MITF↓ SOX10↓ NGFR↓ AXL↑ | Least differentiated |

**Balderson et al. 2024** (BFG, elad055) already re-analyzed Tirosh 2016 with
classical methods (PCA + Monocle) and recovered these four states. Our project's
novelty is **applying deep learning methods** to the same dataset and benchmarking
against this classical result.

---

## Directory Structure

```
melanoma-scrnaseq-reanalysis/
├── CLAUDE.md                  # This file — project context for Claude Code
├── README.md                  # Public project overview
├── LICENSE                    # MIT
├── .gitignore
├── environment.yml            # Conda environment spec (scanpy/scvi-tools stack)
├── data/
│   ├── raw/                   # GSE72056_*.txt.gz — git-ignored, immutable
│   └── processed/             # Derived AnnData / intermediates — git-ignored
├── notebooks/
│   ├── 00_data_exploration.ipynb     # Public English notebooks (committed)
│   └── .local_zh/                    # Personal Chinese working notebooks
│       └── 00_data_exploration_zh.ipynb
├── src/                       # Reusable Python modules (loaders, pipeline)
├── results/
│   ├── figures/               # Generated plots — git-ignored
│   └── tables/                # Generated tables — git-ignored
└── docs/
    └── correspondence/        # e.g. tirosh_email_draft.md
```

---

## Workflow Conventions

**Code generation principles:**
1. The student designs the analysis; Claude Code implements.
2. Every analytical decision must be documented in the relevant notebook
   or design doc before code is written.
3. Always explain code logic when generating — the student must understand
   what they commit.
4. Prefer clear, well-commented code over clever one-liners.

**Notebook conventions:**
- English notebooks (`notebooks/NN_title.ipynb`) are public and committed.
- Chinese working notebooks (`notebooks/.local_zh/`) are personal aids.
- Number notebooks sequentially (00_, 01_, 02_, …).

**Commit message style:**
- Concise English, present tense ("Add QC pipeline" not "Added QC pipeline")
- Reference Stage/Quest when relevant ("Stage 0: Add data exploration notebook")

**What NOT to do:**
- Do not re-normalize or re-log-transform the Tirosh data.
- Do not attempt to run scVI / scANVI / foundation models directly on the
  log-TPM data without addressing the count requirement.
- Do not commit large data files to git (`data/raw/` is gitignored).
- Do not invent biological interpretations — defer to Tsoi 2018 markers and
  established melanoma literature.

---

## External References

- Tirosh et al. 2016 (*Science*): DOI 10.1126/science.aad0501
- Tsoi et al. 2018 (*Cancer Cell*): DOI 10.1016/j.ccell.2018.03.017
- Balderson et al. 2024 (*BFG*): DOI 10.1093/bfgp/elad055
- Heumos et al. 2023 (*Nat Rev Genet*): DOI 10.1038/s41576-023-00586-w
- Jerby-Arnon et al. 2018 (*Cell*) — stretch dataset: DOI 10.1016/j.cell.2018.09.006

---

## Open Issues / TODOs

- [ ] Confirm whether Dr. Tirosh can provide raw count matrix (email sent 2026-05-17)
- [ ] Update Balderson 2024 reading notes Q3/Q4 (current draft has generic limitations)
- [ ] Q2.3: Harmony batch correction — 15 patient batches (not 19); decide
      ultra-small-patient handling (malignant cells: 75=3, 65=4, 60=9, 94=10)
- [ ] Decide handling of immune contamination in Tirosh "malignant" cells
      (top HVGs are immune genes) — kept as documented limitation for now to
      preserve Balderson 2024 comparability; revisit after Q2.3

---

## Development

This project uses Claude Code with project-specific context (see `CLAUDE.md`)
for AI-assisted implementation. All analytical decisions and interpretations
are made by the author.
