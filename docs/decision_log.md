# Decision Log

> Chronological record of the project's key methodological and strategic
> decisions. Each entry states the trigger, the decision, the reasoning, and
> its strategic significance. Intended to be referenced in application
> materials (PIQ / essays).
>
> **Author collaboration note:** Decisions and judgments are Ewan Wu's own;
> entries are structured with AI assistance (Claude).

---

## 2026-05-19: Stage 3 Methodological Pivot

> *Supersedes the earlier "data not usable" wording in commit `151704f`; that
> phrasing was imprecise — the raw counts exist but are DUOS controlled-access.
> CLAUDE.md and the Stage 2 report were reconciled to this entry.*

### Trigger
Email from Dr. Itay Tirosh (2026-05-19) confirmed that raw counts for
GSE72056 are accessible only through DUOS (Broad Institute's controlled-
access portal at duos.broadinstitute.org).

### Decision
Abandon scVI / scANVI / foundation models that require raw counts as input.
Pivot Stage 3 toward integration and embedding methods that accept
log-normalized data.

### Reason
- DUOS access typically requires institutional sponsorship (PI affiliation,
  Data Use Agreement) which is unavailable as a community college student.
- Estimated approval timeline (2-6 weeks) is incompatible with the 6-month
  project timeline.
- Cost-benefit analysis: pursuing DUOS access risks losing 4-6 weeks
  with non-trivial probability of rejection. The opportunity cost is
  higher than the methodological gain.

### New Stage 3 Framing
"Comparison of integration methods that accept log-normalized expression
data — Harmony, Seurat RPCA, scGen, UCE — in their ability to recover the
Tsoi 4-state model from Tirosh 2016 data. A real-world methodological
assessment of what is achievable under data access constraints common to
historical scRNA-seq datasets."

### Strategic Significance
This constraint is not a project weakness but a real-world research
condition. Many practitioners working with historical scRNA-seq datasets
face the same limitation. The constraint becomes part of the project's
honest framing — what can be learned about cell state structure when
only normalized expression is publicly available.

### Documentation
Original Tirosh email saved at docs/.private/correspondence/
(gitignored — contains personal contact information).

> **Note:** Revised by 2026-05-19 (b) entry below — UCE was later ruled out
> after the Q3.1 feasibility audit; Stage 3 Core narrowed from four methods
> to three (Harmony / Seurat RPCA / scGen).

---

## 2026-05-19 (b): UCE removed from Stage 3 Core

> *This entry refines the 2026-05-19 (a) pivot decision. Following a
> formal feasibility audit (docs/stage3_feasibility_audit.md), UCE was
> investigated and found to share the same raw-count dependency as
> scVI/scANVI, contrary to my initial expectation. This entry records
> that finding.*

### Trigger
Stage 3 feasibility audit (Q3.1) examined the live source code and
preprocessing pipeline of UCE (snap-stanford/UCE). The README explicitly
states `.X` should contain scRNA-seq counts, and the internal pipeline
applies `log1p(x / sum(x) * 1000)` — a transformation that assumes
integer counts as input.

### Decision
Remove UCE from the Stage 3 Core method comparison. Document UCE in the
same "investigated but ruled out due to raw-count dependence" category
as scVI / scANVI / count-tokenized foundation models (Geneformer, scGPT).

### Reason
1. Feeding log2(TPM/10+1) values into UCE's internal log1p normalization
   would produce a double-log-transformed input — biologically meaningless
   relative to the model's pretraining distribution.
2. The 33-layer UCE checkpoint requires an A100-class GPU (80 GB), which
   is not available on this project's hardware.
3. Even if both issues were resolved, the resulting embedding would be
   difficult to defend methodologically in a final write-up.

### Revised Stage 3 Core (3 methods, 3 families)
- Harmony (Stage 2, completed) — linear post-PCA correction
- Seurat RPCA — linear anchor-based integration
- scGen — non-linear variational autoencoder

### Strategic Significance
Q3.1 feasibility audit served exactly its design purpose: prevent
multi-week implementation cost on a method that would have failed at
the input-format gate. This kind of pre-implementation audit is itself
a transferable engineering skill — verifying assumptions against
ground-truth source code before committing.

---

## 2026-05-19 (c): Switch Stage 3 Q3.2 from Seurat RPCA (R) to Scanorama (Python)

> *This entry refines the Stage 3 Core set by 2026-05-19 (b). Following an
> implementation attempt of the Seurat RPCA pipeline that the Q3.1
> feasibility audit had cleared, the R-Python bridge proved fragile in a
> way the paper audit did not anticipate. The method-family slot is
> preserved by switching to Scanorama — a Python-native RPCA-equivalent
> that lives in the same algorithm family but in the project's already-
> validated Python stack.*

### Trigger
First execution of `notebooks/05_seurat_rpca.ipynb` (R kernel
`ir-melanoma`, conda env `melanoma-r`) via `jupyter nbconvert` triggered
R's `reticulate` package — pulled in transitively by `zellkonverter` and
the R `anndata` wrapper — to auto-provision a private Python interpreter
via `pyenv`, downloading and source-building Python 3.14.0 along with
openssl-3.6.0 and readline-8.3. The provisioning exceeded nbconvert's
300-second cell timeout (Cell 3, `zellkonverter::readH5AD`) and left an
orphaned `~/.pyenv/versions/3.14.0/` directory. Diagnostic inspection
also surfaced reticulate's modern provisioning mechanism, which creates
a parallel uv-managed Python under `~/Library/Caches/.../reticulate/uv/`.
The R-side pipeline therefore requires multi-step Python-provisioning
intervention before the Tirosh data even loads.

### Decision
Switch the Stage 3 Q3.2 implementation from **Seurat v5 RPCA** (R) to
**Scanorama** (Hie et al. 2019, *Nature Biotechnology* —
[github.com/brianhie/scanorama](https://github.com/brianhie/scanorama)).
Scanorama is the same method family as Seurat RPCA — mutual-nearest-
neighbors anchoring in batch-split PCA space — implemented natively for
the scanpy stack. It runs entirely in the project's existing
`melanoma-scrnaseq` conda env with a single `pip install scanorama`.

Naming convention going forward:
- Notebook: `notebooks/05_rpca.ipynb` (replaces the historical
  `05_seurat_rpca.ipynb`, which is preserved as a record of the
  attempted R-side path)
- `adata.obsm` key: `X_rpca` (preserves downstream notebook 08
  cross-method comparison contract)
- Method label in CLAUDE.md spec, audit doc updates, and the final
  write-up: **"RPCA-family integration (Scanorama)"** — implementation
  is named explicitly so the framework choice is honest in reporting.

### Reason
1. The R-Python bridge fragility — discovered at runtime, not at audit
   time — represents an open-ended environment-setup cost (manual Python
   provisioning, reticulate reconfiguration, possible repeat across
   sessions) with no methodological return. The same RPCA-family
   integration is available natively in Python.
2. The `melanoma-scrnaseq` env already hosts the full required stack
   (scanpy, anndata, numpy, scikit-learn) and is the same env that
   Q3.3 scGen will use. Single-env Q3.2 + Q3.3 is operationally simpler.
3. Scanorama is widely cited in the scanpy / scvi-tools integration-
   benchmark literature, so the cross-method comparison narrative remains
   defensible.
4. The Q3.1 audit's core conclusions still stand — the audit was a
   *paper-feasibility* check, and Scanorama satisfies the same gates
   (accepts log-normalized input, no GPU required, supports explicit
   batch key, runs on the project's hardware).

### Method-family note (for the final write-up)
Stage 3 still spans three method families and three implementations:

| Family | Implementation (Stage 3) |
|---|---|
| Linear post-PCA correction | Harmony (Stage 2 baseline) |
| Linear anchor-based / MNN in PCA space | **Scanorama** (Q3.2; *was* Seurat RPCA) |
| Non-linear VAE | scGen (Q3.3) |

The implementation switch within the "anchor-based / MNN" family is a
documented tooling choice, not a methodological change. The final
write-up will state the implementation tool by name and reference this
decision log entry.

### Strategic Significance
Real-research outcome: when one implementation path proves
environmentally fragile, the right move is to switch tooling while
preserving the methodological role. Documenting the switch — rather than
silently swapping notebooks and pretending the R attempt never happened —
is the academic-honesty version of this. The previously-committed
Seurat-side artifacts (commits `f32d296` through `d605024`) are
preserved as historical record of the attempted path.

### Housekeeping (deliberately *not* done)
- `~/.pyenv/versions/3.14.0/` (empty orphan dir): left in place. Not
  cleaned because (a) it is harmless and (b) the failure mode would
  recur if any future R-side reticulate work is attempted without a
  `RETICULATE_PYTHON` guard.
- `~/Library/Caches/org.R-project.R/R/reticulate/uv/cache/...` (~150 MB):
  left in place for the same reason.
- `melanoma-r` conda env: left installed. Available if a future R-only
  method (Stage 3 stretch or beyond) is needed; not used for Q3.2.
- `notebooks/05_seurat_rpca.ipynb` + `scripts/build_notebook_05.py`:
  left committed as historical record. The new Q3.2 work uses a new
  notebook file (`05_rpca.ipynb`) without a build script per project
  lead's instruction.

### Documentation
- New (Python) notebook: `notebooks/05_rpca.ipynb`
- New dependency: `scanorama` added to `environment.yml` (pip section)
- CLAUDE.md Stage 3 Working Spec updated in the same commit as this entry.

---

## 2026-05-19 (d): Q3.2 algorithmic pivot — Scanorama to BBKNN

> *Further refines Q3.2 implementation. Scanorama was investigated and
> found to fail silently on Tirosh data; BBKNN is selected as the
> RPCA-family replacement.*

### Trigger

Scanorama integration was attempted with both default `dimred=30` and
reduced `dimred=5`. In both configurations:
- Verbose log printed only `Processing datasets (0, 1)` — Scanorama
  processed only the first pair of patients (pt75=3 cells, pt65=4 cells,
  sorted by ascending cell count) and terminated without processing
  the remaining 13 patients.
- Patient silhouette score on integrated embedding: 0.34 (`dimred=30`)
  and 0.46 (`dimred=5`) — both higher than the raw PCA baseline (0.33),
  meaning Scanorama did not integrate the data.
- The integrated embedding was structurally equivalent to per-batch
  PCA projection, with no cross-batch alignment.

### Decision

Abandon Scanorama for Stage 3. Reassign the second slot of the method
comparison to BBKNN (Polański et al. 2020, batch-balanced k-nearest
neighbors), which operates at the kNN graph level and does not require
per-batch SVD.

### Reason

1. Scanorama's algorithm requires per-batch principal component analysis
   to compute reciprocal projections. With 6 of 15 patients having fewer
   than 30 cells (pt75=3, pt65=4, pt60=9, pt94=10, pt84=14, pt53=16),
   the small-batch SVD enters degenerate regimes. The silent failure
   to process pairs beyond (0, 1) is consistent with such a degeneracy
   in pair-selection or termination logic.
2. BBKNN's algorithm is structurally tolerant of small batches: it
   operates on the post-PCA neighbor graph and treats small batches
   as a few additional graph nodes, not as data subsets requiring
   independent dimensionality reduction.
3. BBKNN remains in a different algorithmic family from Harmony
   (post-PCA cluster shift) and scGen (non-linear VAE), preserving the
   3-family comparison structure.

### Method comparison structure after this pivot

- Method 1: Harmony — linear post-PCA cluster shift correction
- Method 2: BBKNN — graph-based batch-balanced neighbor correction
  (previously planned: Seurat RPCA → Scanorama → BBKNN)
- Method 3: scGen — non-linear conditional VAE

### What changes in CLAUDE.md

- Notebook 05 renamed from `05_rpca.ipynb` to `05_bbknn.ipynb`
  (commit will reflect the rename)
- Stage 3 Working Spec "Implementation order" updated to reflect
  BBKNN as method 2

### Strategic significance

The Scanorama investigation produced two empirical findings worth
documenting in the Stage 3 mini-report:
1. Scanorama silently fails on datasets with very small batches in a
   way that does not raise errors or warnings. Users should validate
   integration with patient silhouette diagnostics.
2. The algorithmic assumption of per-batch SVD becomes a hard constraint
   on small-sample scRNA-seq integration. BBKNN's graph-level approach
   sidesteps this constraint.

Both findings are research observations about the algorithms themselves,
not project failures.

---

## 2026-05-19 (e): BBKNN approximation vs exact nearest neighbors

> Documents a reproducibility finding during Q3.2 implementation. Adds
> approx=False as a project standard for BBKNN.

### Trigger

Initial Q3.2 BBKNN runs used the default approx=True (annoy-based approximate
nearest neighbors). Multiple executions of the same code produced different
cluster partitions — res=0.8 fluctuated between 6 and 7 clusters across runs,
and the canonical run produced a 4-cluster partition at res=0.5 with two
batch-pure clusters (pt79=98%, pt59=98%).

### Decision

For all BBKNN runs in this project, set approx=False (exact nearest neighbors).

### Reason

1. Reproducibility: For a portfolio project intended to be re-runnable by
   reviewers or by the author at a future date, stochastic clustering output
   is unacceptable. With 1257 cells, exact kNN is computationally negligible.

2. Genuine algorithmic difference, not just noise: approx=False at res=0.5
   produced a 3-cluster partition with zero batch-pure clusters; approx=True
   produced a 4-cluster partition with two batch-pure clusters (pt79, pt59).
   The "extra" approx=True clusters appear to be artifacts of the
   approximate kNN graph rather than biological structure. This is itself
   a methodological finding worth recording.

3. The corrected (approx=False) partition is also the cleaner method-comparison
   result: BBKNN integrates pt79 and pt59 into biological clusters rather than
   leaving them as patient-specific artifacts. This strengthens the Q3.2
   finding that graph-based batch correction outperforms post-PCA cluster shift
   on this dataset for patient integration.

### Strategic significance

This finding will inform the Stage 3 mini-report's discussion of reproducibility
in graph-based clustering — a common but often-undocumented source of variability
in scRNA-seq pipelines.

---

## 2026-05-22 (f): Q3.3 scGen ruled out — API-contract drift and library abandonment

> *Stage 3 method comparison changes from 3-method to 2-method following
> runtime verification that scGen is unusable with modern scvi-tools.*

### Trigger

Q3.3 implementation began with installation of scGen into the
melanoma-scrnaseq Python environment. Two installation paths were tested:

**Path 1 — PyPI release (scgen 2.1.0)**: Installation succeeded.
`import scgen` failed immediately with `ModuleNotFoundError: No module
named 'scvi._compat'`. The module had been removed from scvi-tools 1.0+
(2023), and scgen 2.1.0 still imports it.

**Path 2 — GitHub master (scgen 2.1.1, one patch ahead of PyPI)**: GitHub
issue #97 reports that the master branch installs past the PyPI import
failure. Runtime testing on minimal fake data confirmed this and went
further:

- ✅ `import scgen` succeeds (`_compat` reference patched)
- ✅ `scgen.SCGEN.setup_anndata()` succeeds
- ✅ `scgen.SCGEN()` object construction succeeds (`LossRecorder → LossOutput` rename patched)
- ✅ `model.train(max_epochs=1)` completes
- ❌ `model.get_latent_representation()` fails:

```
File scvi/model/base/_vaemixin.py:348, in get_latent_representation
    qz: Distribution = Normal(qzm, qzv.sqrt())
                              ^^^^^^^^
AttributeError: 'NoneType' object has no attribute 'sqrt'
```

- ❌ `model.batch_removal()` fails with the same error (`scgen/_scgen.py:220`
  calls `get_latent_representation()` internally)

### Root cause — API-contract drift, not architecture

scGen's module **is a genuine variational autoencoder** — class `SCGENVAE`,
a probabilistic encoder (`scvi.nn.Encoder`), a computed posterior variance
`qz_v`, and a KL term in its loss. Its `inference()` returns the posterior
as separate tensors keyed `qz_m` / `qz_v` / `z` — the scvi-tools 0.x
convention.

scvi-tools 1.x changed the contract: `VAEMixin.get_latent_representation()`
expects either a single `qz` `Distribution` object, or the renamed keys
`qzm` / `qzv` (no underscore). scGen still emits `qz_m` / `qz_v`, satisfies
neither, so scvi-tools 1.4.2 resolves both to `None` and crashes at
`Normal(qzm, qzv.sqrt())`. `batch_removal()` fails identically — it calls
`get_latent_representation()` internally (`scgen/_scgen.py:220`).

The gap is small in principle (scGen's `inference()` would need to emit a
`qz` distribution — a few lines). It is not fixable for this project:
patching means forking/vendoring a dependency; upstream is unmaintained
(#90 closed without fix, #97 open with no maintainer response); a
self-patched integration method cannot be honestly defended in a
methodological write-up. The cause is **library drift + abandonment**, not
architectural impossibility.

### Open GitHub issues confirming abandonment

- **Issue #90** (opened 2023-10-24, closed) — "tutorial does not run,
  version issues with scvi-tools". Closed without resolution.
- **Issue #97** (opened 2024-08-01, still open) — "dependencies update
  required". A user reports that GitHub master installation passes
  `_compat`, but no maintainer has responded for ~21 months. Our
  verification confirms the user's report and also demonstrates that the
  remaining latent-extraction failure is not solvable at the user level.

### Decision

scGen is removed from the Stage 3 Core method comparison. Stage 3 final
method set:

- **Method 1**: Harmony — linear post-PCA correction (completed in Stage 2)
- **Method 2**: BBKNN — graph-based batch correction (completed in Q3.2)

scGen joins UCE (entry b) and Scanorama (entry d) as "investigated but
ruled out". Each has a distinct failure category:

| Method | Failure category | Reference |
|--------|------------------|-----------|
| scVI / scANVI | raw-count dependency | entry a |
| UCE | raw-count dependency + GPU requirement | entry b |
| Seurat RPCA | R/Python interop instability | entry c |
| Scanorama | algorithm-data mismatch on small batches | entry d |
| scGen | API-contract drift + library abandonment (scvi-tools 0.x → 1.x convention change unaddressed upstream) | this entry |

### Ecosystem-level observation

The "log-input + non-linear (VAE-style) + actively maintained + compatible
with current scvi-tools" intersection is **empty** in the 2026 Python
single-cell ecosystem:

- scVI / scANVI: actively maintained, but require raw counts
- scGen: accepts log data (and is a legitimate VAE), but emits an outdated
  API contract (`qz_m` / `qz_v` naming) incompatible with current
  scvi-tools (`qzm` / `qzv`); upstream unmaintained
- DESC and other scvi-derived methods: likely share the same
  incompatibility (not tested)

For log-input scRNA-seq integration in 2026, the actively maintained
methods are dominated by **linear** (Harmony) and **graph-based** (BBKNN)
approaches. Non-linear VAE-based integration on log-input data requires
either deprecated tooling (scGen) or raw counts (scVI / scANVI). This is
itself a methodological finding worth documenting.

### Lessons learned (for future feasibility audits)

The Q3.1 feasibility audit evaluated scGen by reading documentation and
source code. It did NOT include a runtime install + smoke test
verification. As a result, the incompatibility was discovered only at Q3.3
implementation time (~1.5 hours into setup).

**Audit procedure for any future method consideration in this project (or
others)** must include:

1. Documentation / API review (what Q3.1 already does)
2. Runtime install verification — does the package install cleanly into the
   project environment?
3. Smoke test — can the package's primary API be called on minimal fake
   data without errors?

Steps 2 and 3 catch incompatibilities that documentation review alone
cannot. This is added as a project methodology note.

### Verification evidence retained in repo

- `scripts/smoke_test/q33_scgen_smoke.py` — original PyPI smoke test (failed
  at import)
- `scripts/smoke_test/q33_scgen_github_smoke.py` — GitHub master smoke test
  (failed at latent extraction)
- `scripts/smoke_test/q33_scgen_batch_removal_smoke.py` — batch_removal
  alternative-path test (failed at the same internal call)
- `docs/environment_backup_before_scgen.yml` — environment state before the
  scGen install was attempted

These artifacts allow any future reviewer to independently verify the
ruled-out decision.

### Pollution status

The melanoma-scrnaseq environment retains ~33 packages added during the
scGen install attempt (scgen, scvi-tools, lightning, pyro, etc.). These are
inert. Decision: leave them. Will be cleaned at Stage 4 project
finalization.

---

## 2026-05-26 (g): Stage 4 scVelo ruled out (audit-stage)

> Documented at Stage 4 feasibility audit. scVelo is ruled out before 
> any implementation work, following the Q3.1 audit pattern (entry b 
> for UCE).

### Trigger

Stage 4 audit (`docs/stage4_feasibility_audit.md` §3.3) evaluated 
scVelo as a trajectory inference method candidate for addressing 
P1b's continuum framing.

### Primary blocker

scVelo computes RNA velocity from spliced/unspliced read count ratios, 
which require splice-aware alignment of raw FASTQ reads (via velocyto, 
kallisto-bustools, or STAR-solo). Both Stage 4 candidate datasets 
(GSE72056 Tirosh, GSE115978 Jerby-Arnon) place raw FASTQ behind 
controlled-access portals (DUOS and dbGaP respectively). Neither 
public deposit contains the spliced/unspliced separation scVelo 
requires.

### Secondary methodological concern

Even with FASTQ access, RNA velocity on Smart-seq2 (full-length 
non-UMI) data is methodologically off-label — scVelo's published 
validations are predominantly on 10x Chromium UMI data, where 
spliced/unspliced distinction is cleaner. Smart-seq2 velocity is 
possible in principle but adds layered methodological caveats.

### Decision

scVelo ruled out at Stage 4 audit. Same blocker class as Stage 3 
entries a (scVI/scANVI), b (UCE) — raw reads required but 
controlled-access.

### Stage 4 trajectory inference path

PAGA (scanpy core) selected as the trajectory inference method for 
Stage 4 C2 (verified runtime-importable; no new dependencies; uses 
existing kNN graph). See `docs/stage4_feasibility_audit.md` §3.2.

### Strategic significance

Stage 4 now has 7 ruled-out methods documented across entries (a–g), 
spanning multiple distinct failure categories: raw-count dependency 
(a, b, g), R/Python interop (c), small-batch failure (d), API drift 
(f). Pre-implementation audit-stage rulings (b, g) demonstrate the 
Q3.1 forward fix (runtime verification before implementation) in 
sustained practice.

---

## 2026-05-26 (h): Stage 4 C1 framing revised — Jerby-Arnon cohort composition

> Discovered during Boundary B (Phase 2 schema inspection of GSE115978). 
> Audit framing assumed Jerby-Arnon was "an independent dataset"; cohort 
> composition reveals 58% Tirosh-overlap. Stage 4 C1 redesigned as 
> 3-layer framework.

### Trigger

Stage 4 Boundary B Phase 2 (schema inspection of 
`GSE115978_cell.annotations.csv.gz`) revealed that the `Cohort` 
column splits cells into two distinct subsets:

- `Cohort == 'Tirosh'` (4199 cells, 58%) — Tirosh 2016 cells re-annotated by Jerby-Arnon lab
- `Cohort == 'New'` (2987 cells, 42%) — Jerby-Arnon's newly acquired cells

Evidence: Tirosh-cohort cell IDs follow Tirosh naming convention 
(`cy78_*`, `cy79_*`, `CY88_*`); `samples` values match Tirosh patient 
identifiers (Mel79 ↔ cy79). Confirmed not coincidental.

### Audit framing error

Stage 4 feasibility audit (commit 9cf2c69) §2.1 and §2.6 framed 
Jerby-Arnon as "the natural complement: a larger melanoma scRNA-seq 
cohort" for "cross-dataset validation". This framing was incomplete: 
58% of cells are not an independent dataset but a re-annotation 
of Tirosh.

### Revised framing — three layers

See `docs/stage4_dataset_schema.md` §5 for detail. Summary:

1. **Layer 1 (Tirosh-subset, 4199 cells)**: Same cells, annotation-framework meta-comparison (our Tsoi vs Jerby-Arnon cell.types)
2. **Layer 2 (New cohort, 2987 cells)**: True independent validation of P1's dataset-level NGFR claim
3. **Layer 3 (7186 cells)**: Cross-cohort heterogeneity within one publication

This is arguably *richer* than the original binary framing — three 
scientific questions instead of one — but the original "independent 
validation" framing must not be carried forward.

### Forward fix — audit procedure extension

Stage 3 Q3.1 audit forward fix (entry f) added "runtime install + 
smoke test" to method-level audit. This entry adds:

**Audit procedures must include dataset composition verification** 
when a dataset is framed as independent / complementary. Specifically:
- For any "cross-dataset" framing, inspect `Cohort` / source columns 
  in annotation metadata before declaring independence.
- Cell-ID naming convention cross-referencing with potential overlap 
  datasets.
- Sample-/patient-identifier overlap analysis.

The Stage 4 audit (entry g context) verified method runtime (PAGA 
import) but did not include dataset composition verification. 
Forward Stage 5+ audits include this step.

### Stage 4 audit doc handling

Stage 4 feasibility audit (`docs/stage4_feasibility_audit.md`, commit 
9cf2c69) is left unmodified to preserve the audit's original framing 
as part of the project's transparent research trajectory. This 
schema document and decision log entry serve as the documented 
revision; future readers should consult `docs/stage4_dataset_schema.md` 
§5 for the operating Stage 4 C1 framework.

---

## 2026-05-26 (i): Stage 4 Phase 2.6 statistical naïveté — mentor χ² narrative reach

> Self-documented mentor error. Phase 2.6 Layer 2 treatment hypothesis 
> chi-square returned p < 0.0001 with 100% post-treatment in NC cluster. 
> Mentor immediately reached for "treatment-induced biology" narrative 
> framing for Stage 3 P1 revision. Phase 2.7 sanity checks revealed the 
> result reflects patient/treatment confounding, not biology.

### Trigger

Phase 2.6 (`notebooks/08_jerby_arnon_load.ipynb` Cell 31-32): NC + 
Transitory clusters showed extreme treatment-group composition (NC 
72/72 = 100% post, Trans 21/28 = 75% post; combined NGFR-high 93/100 = 
93% post). Chi-square test: χ² = 36.58, p < 0.0001.

Mentor's immediate framing was: "NGFR-high states are treatment-induced 
dedifferentiation, consistent with Tsoi 2018 + BRAF resistance 
literature. P1 revised to 'Tirosh failed to recover NGFR states because 
cohort lacks treatment-exposed cells.'"

### What Phase 2.7 sanity checks revealed

3 checks invalidated the inference:

1. **Patient breakdown**: NC + Transitory cells 74% from single patient 
   Mel110, 93-97% from 2 patients (Mel110 + Mel106). Not "treatment-driven 
   biology" — "Mel110-driven biology".

2. **Treatment-patient confound**: All 11 Layer 2 patients are 100% 
   one-way (no within-patient paired naive/post samples). Treatment 
   effect mathematically inseparable from patient identity at this 
   design.

3. **Tirosh cohort treatment status**: Jerby-Arnon's Tirosh cohort 
   contains 1803 post-treatment cells (43% of 4199). The proposed 
   "Tirosh = naive cohort" hypothesis underlying the treatment-biology 
   narrative is empirically false.

### Mentor's error

I (mentor) reached for biological narrative based on χ² magnitude 
without analyzing experimental design. Specifically:

- I did not request patient breakdown before interpreting χ² (Phase 2.6 
  → 2.7 progression was prompted by Ewan + Claude Code's sanity 
  surfacing, not mentor-initiated)
- I treated pseudo-replicated cell-level data as if independent — 
  860 cells from 11 patients are not 860 independent observations
- I did not verify the Tirosh-cohort treatment status assumption before 
  building narrative on it
- I was prior-anchored: Phase 1.2 verdict was "P1 strengthened"; when 
  Phase 2 data appeared to contradict, I over-corrected and latched 
  onto the first available alternative (treatment biology) rather than 
  continuing to test

### Forward fix

For all future statistical claims with χ², ANOVA, regression, etc., on 
multi-patient single-cell data:

1. **Verify independence assumption first**. Pseudo-replicated cells 
   from N patients should not be treated as N_cells independent 
   observations. Report patient-level summaries; use 
   patient-as-random-effect models when inferential statistics matter.

2. **Pre-commit to design analysis before result interpretation**. When 
   potential confounds (patient ↔ treatment, batch ↔ condition) exist, 
   analyze the confound structure first. Do not interpret the result 
   until confound is ruled out.

3. **Statistical effect sizes that look "too clean" (e.g., 100% one-way 
   in 72 cells) are suspicious, not evidence**. 100% post-treatment in 
   72 cells from 2 patients is much weaker evidence than 70% post in 
   72 cells from 10 patients.

4. **Mentor role is to slow down narrative when verdict is too clean, 
   not accelerate**. χ² p < 0.0001 should trigger "verify this isn't 
   pseudo-replication" instinct, not "build narrative" instinct.

### Stage 4 narrative consequence

After Phase 2.7-2.8 sanity checks, Stage 4 Layer 2's final finding is:

- Stage 3 P1 (Tirosh NGFR weakness) survives patient-level scrutiny
- Layer 2 NGFR-cluster recovery is single-tumor-driven (Mel110, 91%)
- P1 refined to cohort-composition mechanism (not uniform NGFR 
  suppression, but absence of NGFR-rich tumors in Tirosh)
- Treatment-biology hypothesis ruled non-attributable at Layer 2 design

### Reference

`notebooks/08_jerby_arnon_load.ipynb` Cell 30-38 (Phase 2.6 + 2.7 + 2.8 
sequence).
