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
