# Stage 3 Feasibility Audit — Integration Methods for Log-Normalized Tirosh 2016 Data

*Project: Learning to Re-analyze — A Methodological Study of Deep Learning
Approaches on the Tirosh et al. (2016) Melanoma Single-Cell Dataset*
Stage 3 preparation · Q3.1 · 2026-05-19

> **Purpose.** Before any code is written, verify that each candidate Stage 3
> integration method (Seurat RPCA, scGen, UCE) is actually compatible with
> our $\log_2(\text{TPM}/10+1)$ Tirosh data, our hardware (no NVIDIA GPU,
> Apple Silicon Mac), and our annotation status (no ground-truth cell types).
> This is a research deliverable — no models are trained at this stage.
>
> **Method.** Live documentation review for each candidate (latest READMEs,
> source code, and official tutorials), explicitly bypassing potentially-stale
> training-data knowledge. Each method is graded along the six axes specified
> in the Stage 3 design document.
>
> **Author collaboration note:** Findings and verdicts are Ewan Wu's own;
> the document is structured with AI assistance (Claude).

---

## 1. Executive Summary

| Method | Accepts log-input | GPU required | Needs cell-type labels | Verdict |
|---|---|---|---|---|
| **Seurat RPCA** | Yes (with one minor adaptation) | No | No | **CORE** — proceed |
| **scGen** | Yes (explicitly recommended) | Optional, helpful | **Yes** (workaround required) | **CORE** — proceed with caveats |
| **UCE** | **No** — internal pipeline assumes raw counts | Yes (≥1× A100-class) | No | **DROP from Core** — defer / document as ruled out |

**Headline finding.** The previously-assumed inclusion of UCE in the
log-normalized Core group was based on outdated information. Inspection of
the live UCE preprocessing source confirms that the model's internal
pipeline applies its own count-style normalization
(`log1p(x / sum(x) * 1000)`) and that the developers explicitly check the
input scale ("a nice check to make sure it's counts"). Feeding
$\log_2(\text{TPM}/10+1)$ values into this pipeline would produce a
mathematically meaningless double-log transformation. Combined with an
A100-class GPU requirement, UCE belongs in the same bucket as scVI / scANVI
(also raw-count-dependent) and should be **moved out of Core**, leaving the
project with a cleaner three-method comparison (Harmony baseline + Seurat
RPCA + scGen) rather than four.

**Decision-log impact.** The 2026-05-19 entry of `docs/decision_log.md`
currently lists UCE in the log-normalized Core group. Once this audit is
accepted, that entry should be amended (or a follow-up entry added) to
reflect that UCE was investigated and ruled out for the same reason
class as scVI / scANVI — *raw-count dependence*.

---

## 2. Method 1 — Seurat RPCA (Reciprocal PCA)

### 2.1 Input Requirements

Seurat's standard workflow assumes the user starts from raw counts and
applies `NormalizeData()` (log + size-factor normalization) before
`FindVariableFeatures()` → `ScaleData()` → `RunPCA()`. Our Tirosh data is
already in $\log_2(\text{TPM}/10+1)$ form, which is structurally equivalent
to what `NormalizeData()` produces (a log-transformed, depth-corrected
matrix) — just with a different log base and a different scaling factor.

**Compatibility verdict:** Compatible *with one minor adaptation*. Two
options:

1. **Skip `NormalizeData()` and write the values directly to the `data`
   slot.** Seurat v5 stores normalized values in the `data` layer (separate
   from the `counts` layer). We can place our log-TPM matrix in
   `obj[["RNA"]]$data` and leave `counts` empty or NA-filled. RPCA's
   downstream steps (`FindVariableFeatures`, `ScaleData`, `RunPCA`) all
   read from `data`, so this should "just work."
2. **Treat log-TPM as if it were raw counts and skip `NormalizeData()`
   explicitly** by calling `obj <- FindVariableFeatures(obj)` directly
   after object construction. Cosmetically uglier but functionally
   equivalent.

The cleaner option is (1). The Stage 2 Harmony pipeline already established
this principle — no extra normalization on log-TPM — so we are not opening
new methodological ground here.

### 2.2 Computational Cost

R-based, anchor-based integration in low-dimensional PCA space. The
[Seurat v5 RPCA documentation][s5-int] notes RPCA is "faster and more
conservative" than CCA. For our scale (~1,257 malignant cells × 15
patients, or ~4,645 full cells × 19 patients), expected runtime is on the
order of minutes on a laptop. Memory pressure is negligible. **No GPU
required.**

### 2.3 Batch Conditioning

Patient identity is encoded by splitting layers:

```r
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$patient)
```

After splitting, each patient becomes a separate layer; the integration
function re-anchors across these layers. The strength of alignment is
controlled by `k.anchor` (default 5); the [Seurat RPCA tutorial][s5-rpca]
specifically notes that increasing `k.anchor` (e.g. to 20) helps when
populations are weakly aligned — a parameter we may want to vary for the
Tirosh dataset's small patients.

### 2.4 Output Format

`IntegrateLayers(method = RPCAIntegration, ...)` produces a new dimensional
reduction named `integrated.rpca` — a cell-embedding (not corrected
expression). Default 30 dimensions, matching our Stage 2 Harmony choice.
Downstream UMAP / Leiden clustering operates on this reduction:

```r
obj <- FindNeighbors(obj, reduction = "integrated.rpca", dims = 1:30)
obj <- FindClusters(obj)
obj <- RunUMAP(obj, reduction = "integrated.rpca", dims = 1:30)
```

### 2.5 Implementation Path

The friction point is **R / Python interop**. Our existing pipeline is
Python-native (scanpy + AnnData). Three integration paths:

1. **Use Seurat in R via a Jupyter R-kernel notebook.** Convert AnnData
   ↔ Seurat with the `anndata2ri` (Python→R) or `SeuratDisk` /
   `convert2anndata` (R→Python) packages. Highest fidelity.
2. **Wrap the R call from Python via `rpy2`.** Less common in 2026; risk
   of subtle version mismatches.
3. **Use the `scanpy-external` wrapper `sc.external.pp.bbknn` or `mnn`** —
   *not* RPCA equivalents. **Rejected** — would change the method.

Recommended: option (1). Add Seurat (R) to `environment.yml` via
`r-essentials` + `r-seurat` (CRAN). Document the Python→R conversion step
clearly in the Stage 3 notebook so the comparison with Harmony stays
methodologically clean.

### 2.6 Risk Assessment

| Risk | Severity | Mitigation |
|---|---|---|
| R/Python interop (AnnData ↔ Seurat conversion) | Medium | Use `anndata2ri`; verify cell × gene counts match before and after conversion. |
| Forgetting to skip `NormalizeData()` and double-log-transforming | High if missed, easy to spot | Add a printout of `obj[["RNA"]]$data` summary stats; confirm range matches our log-TPM scale (~0–14). |
| Seurat v5 layer API breaking changes | Low (API stabilized) | Pin Seurat version in `environment.yml`. |
| Ultra-small patients (pt75 = 3 cells, etc.) breaking anchor finding | Medium (already a Stage 2 limitation) | Mirror Stage 2 — keep all patients, document the same limitation. |

**Overall risk:** **Low-medium.** Workflow well-documented, deterministic,
no GPU dependency. The R interop is the largest practical friction but is
a known and surmountable engineering problem.

---

## 3. Method 2 — scGen

### 3.1 Input Requirements

The [scGen README][scgen-readme] is explicit: *"We recommend to use
normalized data for the training."* Example code:

```python
import scanpy as sc
adata = sc.read(data)
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
```

This matches the **shape and scale** of our Tirosh data — both are
depth-normalized + log-transformed. We will skip the
`normalize_total + log1p` calls since our data already satisfies the
recommendation.

**Compatibility verdict:** Compatible directly, with the standard "do not
re-normalize already-normalized data" caveat.

### 3.2 Computational Cost

scGen is built on top of [scvi-tools][scvi-tools] (an autoencoder
framework, PyTorch backend). The default training schedule is
`max_epochs=100`, `batch_size=32`, with early stopping. For our ~1,257
malignant cells (or ~4,645 full), training is small by deep-learning
standards. GPU is *helpful* but not required — CPU training should
complete in tens of minutes.

**Apple Silicon (MPS) support:** scvi-tools supports MPS as of mid-2024,
but support is occasionally flaky. If MPS errors arise, fall back to CPU
explicitly (`accelerator="cpu"` in `model.train()`). Document the choice
in the Stage 3 notebook.

### 3.3 Batch Conditioning — and the Critical Caveat

```python
scgen.SCGEN.setup_anndata(adata, batch_key="batch", labels_key="cell_type")
```

**`scGen.batch_removal()` requires `labels_key` with real cell-type
annotations.** Per the official [scGen batch-removal tutorial][scgen-batch],
the labels must be meaningful categorical assignments (a warning is
triggered for categories with fewer than 3 cells, confirming the function
expects real cluster structure, not placeholders).

**This is a genuine methodological issue for us.** We do not have
ground-truth cell-type labels for the Tirosh data — the whole point of the
re-analysis is *to discover* the Tsoi state structure. Available label
sources, in order of decreasing circularity:

| Source | Provenance | Circularity risk |
|---|---|---|
| Original Tirosh malignant flag (0/1/2) | Authors' authoritative annotation | None — coarse but uncontested |
| Original Tirosh non-malignant type code (T/B/Mac/Endo/CAF/NK) | Authors' authoritative annotation | None — coarse but uncontested |
| Stage 2 Harmony-derived Tsoi states (Melanocytic / Transitory / etc.) | Our own marker-based annotation | **High** — using these as scGen input then evaluating scGen's recovery of them is circular |

**Recommended workaround:** Use the **coarse Tirosh-authored labels**
(malignant vs. each non-malignant type) as `labels_key`. This is honest
("we provide scGen with the same labels the Tirosh authors published"),
non-circular (we are not feeding it our Stage 2 marker assignments), and
sufficient for scGen's batch correction mechanism (which only needs labels
to *anchor* the latent space, not to define the states we evaluate).

This means scGen will be run *on the full Tirosh dataset*, not just the
malignant subset — different from the Stage 2 Harmony pipeline. The
malignant subset can still be analyzed downstream by subsetting the
corrected embedding.

### 3.4 Output Format

`model.batch_removal()` returns a new AnnData object containing:
- `adata.obsm['latent']` — the raw latent representation
- `adata.obsm['corrected_latent']` — batch-corrected low-dimensional
  embedding (this is what we feed to UMAP / Leiden)

The default `n_latent` dimensionality is **not stated in the tutorial** —
needs to be confirmed by inspecting `scgen.SCGEN.__init__` defaults or by
explicitly setting `n_latent=10` (the conventional scvi-tools default) to
match our existing 30-dim PCA / Harmony embedding for fair comparison.

> **TODO before coding:** Confirm `n_latent` default and decide on a
> deliberate value for the Stage 3 comparison.

### 3.5 Implementation Path

Pure Python, integrates directly with scanpy. Install via
`pip install scgen` (or conda). Add to `environment.yml`. Verify
compatibility with the current scvi-tools / PyTorch versions.

```python
import scgen

scgen.SCGEN.setup_anndata(
    adata,
    batch_key="patient",
    labels_key="tirosh_celltype_coarse",  # malignant + 6 immune/stromal types
)
model = scgen.SCGEN(adata)
model.train(max_epochs=100, batch_size=32, early_stopping=True,
            early_stopping_patience=25)
corrected = model.batch_removal()
# corrected.obsm['corrected_latent'] -> UMAP / Leiden
```

### 3.6 Risk Assessment

| Risk | Severity | Mitigation |
|---|---|---|
| `labels_key` requirement → potential circularity | High if misused | Use coarse Tirosh-authored labels, never Stage 2 Tsoi assignments. Document. |
| Apple Silicon MPS instability | Low | Fall back to CPU explicitly. |
| Stochasticity (different seeds → different latent spaces) | Medium | Set `scvi.settings.seed`; run 3 seeds, report consensus. |
| n_latent default unknown | Low | Set explicitly. |
| scvi-tools / PyTorch version mismatch in env | Medium | Pin versions; test environment build before notebook work. |

**Overall risk:** **Medium.** The methodological pitfall (label
circularity) is real but resolvable by a deliberate choice of labels. The
engineering risks are standard for deep-learning Python tooling.

---

## 4. Method 3 — UCE (Universal Cell Embedding)

### 4.1 Input Requirements

The [UCE README][uce-readme] states: *"The `.X` slot of the file should be
scRNA-seq counts."* The relevant preprocessing code in `data_proc/`
applies an internal normalization of the form
`log1p(expression / sum(1) * 1000)` and includes a developer comment
(`print(arr.max()); # a nice check to make sure it's counts`) that
betrays the underlying assumption — UCE's pipeline is built around
**integer count input**.

**What happens if we feed log-TPM anyway?** The internal pipeline
computes:

$$
\text{UCE input} \;=\; \log_1p\!\left(\frac{\log_2(\text{TPM}/10+1)}{\sum_g \log_2(\text{TPM}_g/10+1)} \times 1000\right)
$$

This is a **double-log transformation** applied to row-normalized log
values. Genes will still be ranked roughly correctly (monotonicity is
preserved at the cell level), but the absolute magnitudes the model was
trained on no longer correspond to anything biologically meaningful, and
the model's gene-embedding "sentence" representation depends on absolute
expression values (not just rank — confirmed by the source-code finding
that "data remains in continuous normalized space rather than being
thresholded"). The result would be a latent embedding of unclear meaning,
defensibility, and reviewability.

**Compatibility verdict:** **Not compatible.** Same blocker class as scVI
/ scANVI: trained on count-based likelihoods (or count-derived
representations) and not robust to non-count input substitution.

### 4.2 Computational Cost

| Model variant | GPU memory | Realistic hardware |
|---|---|---|
| 4-layer model | batch_size ≤ 100 on 80 GB GPU | 1× consumer GPU (16–24 GB) probably works at reduced batch_size; CPU not supported |
| 33-layer model | batch_size ≤ 25 on 80 GB GPU | 1× A100 80 GB (or H100); **not viable on consumer hardware** |

The user-facing 33-layer model (the one used in the UCE preprint's
benchmarks) is **inaccessible on this project's hardware** (Apple Silicon
Mac, no NVIDIA GPU). The 4-layer model is smaller but the README warns
its embeddings are **not compatible** with the 33-layer model's
embeddings, so we could not directly compare against the published UCE
benchmarks even if it ran.

### 4.3 Batch Conditioning

UCE is a **zero-shot foundation model**: there is no explicit `batch_key`.
Patient effects are expected to be handled implicitly by the model's
broad pretraining distribution. This is an interesting research question
on its own (does pretraining mitigate batch effects without explicit
conditioning?), but only investigable if the input format problem is
first solved.

### 4.4 Output Format

`adata.obsm["X_uce"]` — a 1280-dimensional cell embedding (default
`--output_dim 1280`).

### 4.5 Implementation Path

Even if input incompatibility were resolved:
- Install `snap-stanford/UCE` (PyTorch + HuggingFace Accelerator stack)
- Download 33-layer checkpoint (manual Figshare download, sizable)
- Run `eval_single_anndata.py` from CLI
- Requires CUDA — Apple Silicon MPS not supported in the published scripts

### 4.6 Risk Assessment

| Risk | Severity | Mitigation |
|---|---|---|
| Input format mismatch (raw counts required) | **Blocker** | None viable without DUOS-gated counts |
| GPU hardware unavailable | **Blocker** | Cloud GPU rental possible but expensive and out of scope |
| Defensibility — even if it ran, the results would be hard to defend in a methodological writeup | High | n/a |

**Overall risk:** **Blocked.** Two simultaneous blockers (data format,
hardware) and a defensibility problem on top. UCE belongs in the same
methodological bucket as scVI / scANVI / count-tokenized foundation
models (Geneformer, scGPT) that the project already ruled out.

---

## 5. Cross-Method Comparison Matrix

| Axis | Harmony (Stage 2) | Seurat RPCA | scGen | UCE |
|---|---|---|---|---|
| Family | Linear, post-PCA correction | Anchor-based, linear (RPCA) | Non-linear (VAE) | Foundation model (zero-shot transformer) |
| Accepts log-TPM | Yes (used in Stage 2) | Yes (with `data`-slot adaptation) | Yes (explicitly recommended) | **No** |
| Needs cell-type labels | No | No | **Yes** (workaround: Tirosh coarse labels) | No |
| Explicit batch key | Yes (`theta` controls strength) | Yes (`f` for split layers; `k.anchor`) | Yes (`batch_key`) | No (implicit) |
| Default output dim | 30 (matches PCA) | 30 | TBD (10? 100? — confirm) | 1280 |
| Compute footprint | Trivial (seconds) | Minutes (CPU, R) | Minutes-hours (CPU; faster on GPU) | Hours on A100; infeasible without it |
| Stochasticity | Low (deterministic per seed) | Low | **Medium-High** (random init) | Medium |
| Interpretability of embedding | Good (PCA axes) | Good (PCA axes) | Poor (latent dims) | Poor (1280 black-box dims) |
| Methodological fit for Stage 3 goal | Done | **Strong** — direct linear comparison | **Strong** — adds non-linear axis | Blocked |

---

## 6. Recommendation for Stage 3 Scope

### 6.1 Revised Core Plan

| Method | Status | Notes |
|---|---|---|
| Harmony | Done (Stage 2) | Baseline |
| Seurat RPCA | **Core — implement next** | Classical anchor-based comparator |
| scGen | **Core — implement after RPCA** | Non-linear VAE comparator |
| UCE | **Removed from Core — document as ruled out** | Same raw-count blocker class as scVI / scANVI |

The resulting comparison still provides a clean methodological story:
**linear-PCA-space methods (Harmony, RPCA)** versus **non-linear
latent-space methods (scGen)**, all on the same log-input data. Three
methods, three families. This is a sharper and more honest story than
four methods where one is hand-wavy.

### 6.2 Optional Fourth Method (Stretch)

If a fourth log-input integration method is desired to round out the
comparison, two readily-available candidates accept log-normalized input
without raw counts:

- **BBKNN** (Polański et al. 2020) — graph-based batch correction; very
  different mechanism from Harmony/RPCA; available via
  `sc.external.pp.bbknn`. Fast, deterministic, no GPU.
- **MNN-correct** / `mnnpy` — mutual-nearest-neighbor correction; older
  but still a recognized comparator.

Either could occupy the UCE slot if four methods are wanted, *without
introducing the count-dependency blocker*. Recommendation: prefer BBKNN
if added — it is mechanistically more distinct from PCA-anchor methods,
and currently more actively maintained than mnnpy.

### 6.3 Suggested Implementation Order

1. **Notebook `05`** — Seurat RPCA. R-kernel notebook; document the
   AnnData ↔ Seurat conversion explicitly. Output: RPCA embedding +
   UMAP + Leiden clustering + Tsoi-marker dotplot, parallel to Stage 2
   notebook `04`.
2. **Notebook `06`** — scGen. Python; use coarse Tirosh labels for
   `labels_key`. Output: corrected latent + UMAP + Leiden + dotplot.
3. **Notebook `07`** (optional) — BBKNN if pursuing the 4th-method
   stretch.
4. **Notebook `08`** — Cross-method comparison: Tsoi state recovery,
   silhouette by patient, kNN patient purity, agreement between methods
   (ARI / NMI of cluster assignments).

This staged order ensures each method is fully validated before the next
is added, and the cross-method comparison can begin even if scGen runs
into engineering trouble.

---

## 7. Critical Open Questions (resolve before coding)

| # | Question | Where to answer it |
|---|---|---|
| Q1 | Does writing log-TPM directly to `obj[["RNA"]]$data` and skipping `NormalizeData()` actually produce valid Seurat behavior in v5? | Quick smoke test in an R session before notebook `05`. |
| Q2 | What is the default `n_latent` for `scgen.SCGEN()`? Should we set it explicitly to 10 or 30 for comparability with our PCA / Harmony embedding? | Inspect `scgen.SCGEN.__init__` or run a quick `help(scgen.SCGEN)`. |
| Q3 | Do we use coarse Tirosh labels for the **entire dataset** (all 4,645 cells), and then subset to malignant downstream? Or do we run scGen on the malignant subset only (in which case `labels_key` becomes uninformative)? | Methodological decision — recommend the former; document in notebook `06`. |
| Q4 | Should the Stage 3 comparison be conducted on the malignant-only subset (matching Stage 2 / Balderson) or the full dataset? | Methodological decision — recommend malignant-only for the Tsoi recovery analysis; full dataset only for batch-quality metrics. |
| Q5 | Decision-log amendment: update `docs/decision_log.md` to reflect UCE's ruling out, or add a fresh dated entry? | Recommend a fresh dated entry referencing this audit, preserving the original record. |

---

## 8. References

- Seurat v5 integrative analysis: <https://satijalab.org/seurat/articles/seurat5_integration> [s5-int]
- Seurat RPCA tutorial (original): <https://satijalab.org/seurat/articles/integration_rpca> [s5-rpca]
- scGen GitHub: <https://github.com/theislab/scgen> [scgen-readme]
- scGen batch-removal tutorial: <https://scgen.readthedocs.io/en/stable/tutorials/scgen_batch_removal.html> [scgen-batch]
- scvi-tools: <https://scvi-tools.org> [scvi-tools]
- UCE GitHub: <https://github.com/snap-stanford/UCE> [uce-readme]
- UCE preprint (Rosen et al. 2023, bioRxiv 2023.11.28.568918): <https://www.biorxiv.org/content/10.1101/2023.11.28.568918>
- BBKNN (Polański et al. 2020): <https://doi.org/10.1093/bioinformatics/btz625>

[s5-int]: https://satijalab.org/seurat/articles/seurat5_integration
[s5-rpca]: https://satijalab.org/seurat/articles/integration_rpca
[scgen-readme]: https://github.com/theislab/scgen
[scgen-batch]: https://scgen.readthedocs.io/en/stable/tutorials/scgen_batch_removal.html
[scvi-tools]: https://scvi-tools.org
[uce-readme]: https://github.com/snap-stanford/UCE
