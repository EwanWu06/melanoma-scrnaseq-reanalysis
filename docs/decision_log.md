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
