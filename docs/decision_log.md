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

## 2026-05-19 (b): UCE Removed from Stage 3 Core

### Trigger
Q3.1 feasibility audit (`docs/stage3_feasibility_audit.md`) reviewed live
documentation and source code for each candidate Stage 3 integration method
(Seurat RPCA, scGen, UCE) before any implementation began. The audit
surfaced a compatibility problem with UCE that contradicts the assumption
made in the 2026-05-19 (a) decision.

### Decision
Remove UCE from Stage 3 Core. Stage 3 proceeds with three integration
methods: Harmony (Stage 2 baseline), Seurat RPCA, and scGen. UCE is
documented as "investigated but ruled out" under the same blocker class as
scVI / scANVI.

### Reason
- The UCE README specifies `.X` must contain scRNA-seq counts. The
  `data_proc/` preprocessing source applies an internal
  `log1p(x / sum(x) * 1000)` normalization that assumes count-scale input.
  A developer comment in the source (`print(arr.max()); # a nice check to
  make sure it's counts`) confirms the intended input domain.
- Feeding our $\log_2(\text{TPM}/10+1)$ data into this pipeline would
  produce a mathematically meaningless double-log transformation, with no
  defensibility in the project writeup.
- The published UCE benchmarks use the 33-layer model, which requires an
  80 GB GPU (A100-class). This is unavailable on the project's hardware
  (Apple Silicon Mac, no NVIDIA GPU). The 4-layer model has a lower
  footprint but produces embeddings explicitly noted as incompatible with
  the 33-layer benchmarks.
- UCE thus shares the same blocker class as scVI / scANVI / count-tokenized
  foundation models (raw-count dependence), already ruled out in
  2026-05-19 (a).

### Strategic Significance
This is the value of running Q3.1 as a feasibility audit *before* any
implementation: a one-day documentation review prevented an estimated
1–2 weeks of likely-fruitless effort on UCE. The remaining three-method
comparison is also methodologically cleaner — three distinct families
(linear post-PCA correction / linear anchor-based / non-linear VAE) on
the same log-input data — than four methods where one is hand-wavy.

### Revised Stage 3 Scope
| Method | Family | Status |
|---|---|---|
| Harmony | Linear post-PCA correction | Done (Stage 2 baseline) |
| Seurat RPCA | Linear anchor-based | Q3.2 (next) |
| scGen | Non-linear VAE | Q3.3 |
| BBKNN | Graph-based (optional 4th) | Stretch — decide after Q3.3 |
| UCE | Foundation model | **Removed from Core** — see above |
| scVI / scANVI / scGPT / Geneformer | Count-based | Already ruled out (2026-05-19 a) |

### Documentation
Full audit: `docs/stage3_feasibility_audit.md`.
