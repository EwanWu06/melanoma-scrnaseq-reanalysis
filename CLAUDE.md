# CLAUDE.md — Project Context for Claude Code

> This file provides project-specific context for Claude Code.
> It is automatically loaded when Claude Code operates within this repository.
> Updated: 2026-05-19 (post Q3.1 feasibility audit)

---

## Project Overview

**Project name:** Melanoma scRNA-seq Re-analysis

**Working title (academic):** *Learning to Re-analyze: A Methodological Study of
Deep Learning Approaches on the Tirosh et al. (2016) Melanoma Single-Cell Dataset*

**One-line description:**
A learning-oriented re-analysis of a foundational 2016 melanoma scRNA-seq dataset
(Tirosh et al., GSE72056), comparing three integration methods that accept
log-normalized data (Harmony, Seurat RPCA, scGen) for their ability to recover
the Tsoi four-state model, benchmarked against a recent classical re-analysis
(Balderson et al., 2024). Count-based methods (scVI / scANVI / UCE /
count-tokenized foundation models) were ruled out — raw counts are DUOS
controlled-access (see `docs/decision_log.md`).

**Purpose:**
This project is part of an undergraduate transfer application portfolio.
Target programs include UCLA Computational Biology, UCSD Public Health
(Biostatistics concentration), UNC Biostatistics, and UMich Informatics.
The project prioritizes genuine learning, methodological rigor, and
reproducibility over novel discovery.

---

## Current Stage

**Active stage:** Stage 3 preparation — Stage 2 (Classical Pipeline Replication) complete; see `docs/stage2_report.md`

**Stage progression:**
- [x] Stage 0: Setup & Foundation
  - [x] Repository structure created
  - [x] Tirosh 2016 data downloaded (GSE72056)
  - [x] Initial correspondence with Dr. Tirosh sent
  - [x] Core literature read (Tirosh 2016, Tsoi 2018, Balderson 2024, Heumos 2023)
- [x] Stage 1: Foundation Building
- [x] Stage 2: Classical Pipeline Replication
  - [x] Q2.1: Data loading + AnnData + basic QC (notebook 01)
  - [x] Q2.2: HVG selection + PCA + batch-effect diagnostic (notebook 02)
  - [x] Q2.3: Harmony batch correction + UMAP (notebook 03)
  - [x] Q2.4: Leiden clustering + Tsoi-state annotation (notebook 04)
  - [x] Q2.5: Stage 2 mini-report (`docs/stage2_report.md`)
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
2. scVI, scANVI, UCE, and most count-tokenized foundation models
   (scGPT / Geneformer) **cannot use this data** — they require raw integer
   counts (negative-binomial likelihoods or count-derived gene sentences),
   which are DUOS-gated for GSE72056. Stage 3 therefore does **not** use
   them. UCE specifically was ruled out by Q3.1 feasibility audit (see
   `docs/decision_log.md` and `docs/stage3_feasibility_audit.md`).
3. When selecting highly variable genes, use `sc.pp.highly_variable_genes(flavor="seurat")`
   — NOT `flavor="seurat_v3"` which requires raw counts.
4. Raw counts **exist but are DUOS controlled-access** (confirmed by Dr.
   Tirosh, 2026-05-19); not obtainable on this project's timeline. Stage 3
   therefore uses **integration methods that accept log-normalized data**
   (Harmony, Seurat RPCA, scGen) — **not** count-based scVI / scANVI / UCE.
   Authoritative rationale: `docs/decision_log.md`.

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
Re-analyze Tirosh 2016 with three log-input integration methods spanning
three method families: **Harmony** (linear post-PCA correction — Stage 2
baseline), **Seurat RPCA** (linear anchor-based), and **scGen** (non-linear
VAE). Compare their ability to recover the Tsoi four states against
Balderson 2024's result. (scVI / scANVI / UCE ruled out — counts are
DUOS-gated; see `docs/decision_log.md`.) Deliverables: GitHub repo,
technical report, dashboard.

**Stretch (if time allows):**
- Add **BBKNN** as a fourth method (graph-based batch correction, accepts
  log-input, no GPU) — decision deferred until Q3.3 finishes.
- Extend analysis to Jerby-Arnon 2018 (GSE115978) for cross-dataset
  robustness assessment.

**Reach (probably skip, documented as future work):**
Count- and count-tokenization-based foundation models (**scGPT, Geneformer,
UCE**) — same DUOS raw-count constraint. UCE specifically was investigated
in Q3.1 and ruled out (also requires A100-class GPU). Smart-seq2 vs 10x
distribution shift analysis.

---

## Stage 3 Working Spec

> Frozen 2026-05-19 (post Q3.1 feasibility audit). Used by Claude Code as
> the authoritative spec for Stage 3 implementation work.

**Headline:** Method comparison on log-normalized data.

**Core methods (three families, three implementations):**

| # | Method | Family | Status |
|---|---|---|---|
| 1 | Harmony | Linear post-PCA correction | Done (Stage 2 baseline) |
| 2 | Seurat RPCA | Linear anchor-based | Q3.2 next |
| 3 | scGen | Non-linear VAE | Q3.3 after |

**Stretch:**
- BBKNN (4th method, graph-based) — decide after Q3.3 completes.

**Removed from Core:**
- UCE — raw-count blocker (see `docs/decision_log.md` 2026-05-19 (b) and
  `docs/stage3_feasibility_audit.md`).

**Evaluation scope:** malignant-only subset (1,257 cells) for **all**
cross-method comparisons (state recovery *and* batch quality), consistent
with Stage 2 and Balderson 2024. Running batch-quality metrics on the
**full** dataset (4,645 cells) is deferred to Stretch — it would address a
different question (general-purpose batch correction across heterogeneous
cell types) than the Tsoi-recovery question driving Stage 3.

**Embedding dimensions:** 30 dimensions across all methods (parity with
Stage 2 PCA / Harmony output, and with `n_latent = 30` for scGen).

**Evaluation metrics (Stage 3 comparison notebook):**

- **A. Batch Integration Quality**
  - A1: Patient silhouette score on integrated embedding (lower = better mixing)
  - A2: kNN patient purity (lower = better mixing; theoretical baseline
    depends on patient size distribution)
- **B. Biological Signal Preservation**
  - B1: Cluster-level expression of canonical markers (MITF / SOX10 / NGFR / AXL)
    visualized as dotplot per resolution
  - B2: Tsoi state silhouette score on integrated embedding (higher = states
    more separable)
- **C. Cell State Recovery**
  - C1: Manual annotation — per-method Leiden clustering + Tsoi state
    assignment, following Stage 2 Q2.4 procedure (no algorithmic state
    assignment)
  - C2: Number of biologically interpretable Tsoi states recovered + number
    of patient-driven (batch artifact) clusters identified

**scGen-specific procedure (not a separate metric category):**

- Trained on the full dataset (all 4,645 cells) using coarse Tirosh-authored
  labels as `labels_key` (malignant flag + 6 non-malignant cell type codes).
- Evaluated on the malignant subset (1,257 cells) — embedding subsetted
  downstream for fair comparison with Harmony / RPCA.
- `n_latent = 30`, matching Harmony / RPCA dimensionality for direct
  silhouette comparability.

**Implementation order (notebook layout):**
- `notebooks/05_seurat_rpca.ipynb` — Q3.2 (R kernel; Seurat v5; AnnData ↔ Seurat conversion documented in-notebook)
- `notebooks/06_scgen.ipynb` — Q3.3 (Python; scGen with coarse-label setup)
- `notebooks/07_bbknn.ipynb` — Stretch (optional)
- `notebooks/08_method_comparison.ipynb` — Cross-method evaluation:
  metrics A1/A2/B1/B2/C1/C2, plus cross-method cluster agreement metrics
  (specific choice — e.g. ARI, NMI, or alternative — to be finalized in
  notebook 08)

**Environment:** Stage 3 introduces a second conda env `melanoma-r` for
the R/Seurat workflow. The Python env `melanoma-scrnaseq` remains the
analysis driver; data exchange between R and Python happens via files
(CSV or AnnData h5ad).

> *Note: the `melanoma-r` conda environment was set up during the Q3.1
> Seurat smoke test as the cleanest path for R/Seurat work; this spec
> documents it retroactively.*

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
novelty is **applying modern log-input integration methods across three
different families** (linear post-PCA correction / linear anchor-based /
non-linear VAE) to the same dataset and benchmarking against this classical
result.

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
- Do not pursue count-based scVI / scANVI / count foundation models — raw
  counts are DUOS-gated and out of scope (see `docs/decision_log.md`);
  Stage 3 uses log-input integration methods only.
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

- [x] Dr. Tirosh raw-count request — **resolved 2026-05-19**: Dr. Tirosh
      replied promptly; GSE72056 raw counts **exist but are DUOS
      controlled-access** (Broad portal), requiring institutional sponsorship
      unavailable on this timeline. Stage 3 therefore **pivots to log-input
      integration methods** — authoritative rationale in
      `docs/decision_log.md` (2026-05-19 entry).
- [ ] Update Balderson 2024 reading notes Q3/Q4 (current draft has generic limitations)
- [x] Stage 3 Q3.1 — feasibility audit complete (`docs/stage3_feasibility_audit.md`);
      UCE removed from Core (raw-count blocker, A100 GPU requirement)
- [ ] Stage 3 Q3.2 — Seurat RPCA on Tirosh malignant cells (after R-env
      smoke test passes)
- [ ] Stage 3 Q3.3 — scGen with Tirosh coarse labels (`labels_key`),
      `n_latent = 30` for parity
- [ ] Stage 3 carryover from Stage 2: revisit Harmony theta sensitivity,
      ultra-small patients (malignant 75=3, 65=4, 60=9, 94=10), and immune
      contamination — all documented as Stage 2 limitations in
      `docs/stage2_report.md`

---

## Development

This project uses Claude Code with project-specific context (see `CLAUDE.md`)
for AI-assisted implementation. All analytical decisions and interpretations
are made by the author.
