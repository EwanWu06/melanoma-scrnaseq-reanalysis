# Stage 4 Feasibility Audit — Candidate Directions for Future Work

*Project: Learning to Re-analyze — A Methodological Study of Deep Learning
Approaches on the Tirosh et al. (2016) Melanoma Single-Cell Dataset*
Stage 4 preparation · 2026-05-26

> **Purpose.** Stage 3 mini-report (`docs/stage3_report.md`) names three Stage 4
> candidate directions. Before any implementation is committed, verify each
> candidate's feasibility against this project's data access constraints
> (DUOS-gated Tirosh raw FASTQ, see decision log entry a), environment
> footprint (Python-only `melanoma-scrnaseq` env after Q3.2's R-bridge
> ruling), and methodological fit with the Tsoi state recovery framing.
>
> **Method.** Live documentation review (GEO record for the candidate
> dataset, Bioconductor / PyPI / GitHub for the trajectory methods),
> minimal-cost runtime sanity check on the one trajectory method that lives
> in our existing env (PAGA via scanpy), and a paper-level evaluation of
> the third candidate which requires no external dependency.
>
> **Author collaboration note:** Findings and verdicts are Ewan Wu's own;
> the document is structured with AI assistance (Claude).

---

## 1. Executive Summary

| Candidate | Data access | New env / deps | Verdict |
|---|---|---|---|
| **1. Cross-dataset validation — Jerby-Arnon 2018 (GSE115978)** | Counts + TPM CSV **publicly distributed**; raw FASTQ in dbGaP (not needed for Stage 2-style replication) | None — existing `melanoma-scrnaseq` env sufficient | **VIABLE — strongest scientific lever for P1 dataset-level claim** |
| **2a. Trajectory — Slingshot (R)** | Works on existing reduced embedding (no raw counts required) | R env (currently dormant; Q3.2 RPCA experience suggests bridge cost) | **VIABLE w/ caveats — R/Python bridge friction must be re-priced** |
| **2b. Trajectory — PAGA (Python, scanpy)** | Works on existing kNN graph + cluster labels (no raw counts required) | **None — already importable in `melanoma-scrnaseq` env** (verified) | **VIABLE — lowest implementation cost** |
| **2c. Trajectory — scVelo (Python)** | **Requires spliced + unspliced count matrices from velocyto on raw FASTQ** — both candidate datasets gate raw reads | scvelo (PyPI 0.3.4, actively maintained as of 2026-02) | **BLOCKED — same class blocker as Stage 3 scVI (raw-read controlled access)** |
| **3. Multi-resolution sweep — Harmony + BBKNN** | Existing processed h5ad — no external data | None — uses notebook 07 pipeline directly | **VIABLE — lowest scope, directly addresses Stage 3 L4** |

**Headline findings.**

1. **Jerby-Arnon's public artifact profile is materially different from
   Tirosh's.** Tirosh distributes only `log2(TPM/10+1)`; Jerby-Arnon
   distributes both a raw-count CSV and a TPM CSV (plus cell annotations)
   as supplementary files on GEO. Only the FASTQ reads are dbGaP-gated, and
   FASTQ is not required for any of the Stage 4 questions on the table.
   Cross-dataset validation is **not** re-blocked by the Stage 3
   raw-count constraint class.
2. **scVelo on either melanoma dataset is blocked at the data layer.**
   scVelo's input contract requires spliced and unspliced count matrices,
   which are produced by `velocyto` (or equivalent) from aligned BAM files
   derived from raw FASTQ. Tirosh's FASTQ is DUOS-gated (entry a);
   Jerby-Arnon's is dbGaP-gated. The blocker class is identical to the
   Stage 3 scVI / scANVI / UCE ruling; scVelo lands in the same bucket.
3. **PAGA is the lowest-cost trajectory route** — it is already
   importable in the project's existing conda env (runtime verified) and
   consumes the kNN graph and cluster labels already produced by notebooks
   03 / 04 / 05. No new dependency, no R bridge.
4. **Multi-resolution sweep (candidate 3) is the lowest-scope, highest-
   confidence extension.** It reuses Stage 3's pipeline end-to-end on
   existing inputs and directly addresses L4 of the Stage 3 mini-report.
   The scientific question is narrower than candidates 1 or 2 but the
   feasibility risk is essentially zero.

---

## 2. Candidate 1 — Cross-Dataset Validation: Jerby-Arnon 2018 (GSE115978)

### 2.1 Scientific question

Stage 3 mini-report [P1] established the Transitory / Neural-crest-like
recovery failure as a **dataset-level** NGFR signal limitation, based on
two structurally different integration methods (Harmony + BBKNN) both
failing on the same Tirosh data. This is a strong claim that is properly
testable only against a second dataset. Jerby-Arnon 2018 is the natural
complement: a larger melanoma scRNA-seq cohort acquired with the same
core technology (Smart-seq2). If NGFR is similarly weak across all
clusters there, the dataset-level framing of P1 is empirically reinforced;
if NGFR is robustly expressed in some Jerby-Arnon clusters and intermediate
Tsoi states emerge, the framing must be revised — P1 would become
"Tirosh-specific" rather than "dataset-level".

### 2.2 Data access / format

GEO record GSE115978 distributes three processed supplementary files:
- **Counts CSV** — gene × cell raw count matrix (publicly downloadable)
- **TPM CSV** — TPM-normalized matrix (publicly downloadable)
- **Cell annotations CSV** — author-curated cell type labels

Per the GEO record, "raw reads will be made available through dbGaP" —
raw FASTQ is controlled-access, but **the count matrix itself is public**.
This is a critical asymmetry vs. Tirosh GSE72056, where only
`log2(TPM/10+1)` is distributed publicly.

| Item | GSE72056 (Tirosh) | GSE115978 (Jerby-Arnon) |
|---|---|---|
| Public matrix | `log2(TPM/10+1)` only | **counts + TPM** |
| Raw FASTQ | DUOS-gated | dbGaP-gated |
| Implication for Stage 4 Q's | Count-based methods blocked | **Stage 2-style replication unblocked** |

Cohort size: 7,186 cells across "31 melanoma tumors" per the GEO record
title (the published paper reports 33 tumors — a 2-tumor discrepancy that
should be reconciled during implementation; likely 2 tumors were added
after GEO submission or excluded from the public deposit).

Platform: "Modified full length SMART-Seq2 protocol" on Illumina NextSeq
500 — same core technology as Tirosh, so distributional shift relative to
Tirosh is minimized.

### 2.3 Computational cost

Stage 2 replicated on Tirosh's 4,645 cells × 23,686 genes in trivial time
(Harmony converges in seconds on this scale). Jerby-Arnon's 7,186 cells is
~1.5× larger; gene count unknown from the GEO page but expected to be the
same order of magnitude. Total Stage 2-equivalent pipeline runtime
estimate: **minutes**, fully on a laptop, no GPU.

### 2.4 Implementation effort

Mirror Stage 2 notebooks (00–04) end-to-end on Jerby-Arnon, then re-run
Q3.4-style cross-method comparison if BBKNN is also applied. Estimated
implementation effort:

| Step | Reuse from | New work |
|---|---|---|
| Loading + QC + AnnData | notebook 01 | Adapt CSV path; verify cell-type-annotation column matches |
| HVG + PCA | notebook 02 | Re-run; cohort size differs but parameters likely transferable |
| Harmony + UMAP | notebook 03 | Re-run with `batch_key="sample"` (terminology in this dataset is `samples` not `patient`; verify column name) |
| Leiden + Tsoi annotation | notebook 04 | Re-run; the Tsoi marker dotplot is the key Stage 4 deliverable |
| (Optional) BBKNN | notebook 05 | Re-run; same `approx=False` convention |
| (Optional) Cross-method comparison | notebook 07 | Re-run; same A1/A2/B1/B2/C1/C2 metrics |

Pipeline ports are mostly path / column-name changes. No new dependencies.

### 2.5 Risk assessment

| Risk | Severity | Mitigation |
|---|---|---|
| Cell-annotation column name / schema differs from Tirosh | Low — easy to spot at load time | Inspect annotations CSV before pipeline run |
| Counts file size > Tirosh; memory pressure on laptop | Low (still small-data regime) | Sparse-load if needed |
| Sample-size imbalance differs (Jerby-Arnon's 31-tumor distribution is unknown) | Medium — if some tumors have very few cells, will replicate Tirosh's small-batch issues for non-Harmony methods | Inspect per-tumor cell counts before integrator choice; same as Stage 2's malignant-only `pt75=3 pt65=4` situation |
| 31 vs 33 tumor discrepancy between GEO and paper | Low | Document the chosen count source and reason in the Stage 4 notebook |
| Spliced/unspliced separation absent in the public deposit | Not relevant for candidate 1's questions (only relevant if scVelo is attempted) | n/a |

### 2.6 Strategic significance

This is the strongest single lever for the Stage 3 mini-report's headline
claim. P1's "dataset-level" framing was made cautiously precisely because
it rested on a single dataset; an independent replication on Jerby-Arnon
either confirms it as a general finding or scopes it down to
Tirosh-specific. Either outcome is methodologically informative and
publishable as a Stage 4 deliverable.

Cross-dataset validation is also the candidate most directly aligned with
the existing CLAUDE.md Stretch-tier scope (CLAUDE.md already names
Jerby-Arnon as the cross-dataset stretch dataset).

---

## 3. Candidate 2 — Trajectory Inference

The Stage 3 mini-report [P1b] and [P5] surface a continuum-vs-discrete
framing tension. Trajectory inference is the methodological response:
characterize the dedifferentiation axis as a continuous pseudotime instead
of approximating it with discrete cluster labels. Three candidates were
named; each is evaluated separately, followed by a per-method comparison.

### 3.1 Slingshot (R, Bioconductor)

#### 3.1.1 Maintenance & versioning
- Current version 2.20.0 on Bioconductor 3.23 (R 4.6)
- In Bioconductor "since BioC 3.8 (R-3.5)" — 7.5 years of continuous maintenance
- Maintainer: Kelly Street (active)
- Verdict: **maintained, stable**

#### 3.1.2 Input requirements
Slingshot operates on a **reduced-dimensional embedding plus cluster
labels** — explicitly designed for "after dimensionality reduction and
clustering". It does **not** require raw counts. Our existing
`X_pca_harmony[:, :30]` (or 2D UMAP) and `leiden_r0.3` outputs are
direct-compatible inputs.

#### 3.1.3 Implementation effort
Slingshot is R-only. Two paths:
1. **R-kernel notebook**, converting AnnData ↔ SingleCellExperiment via
   `anndata2ri` or `zellkonverter`. This is exactly the R-Python bridge
   that bit Q3.2 (decision log entry c) when `reticulate` auto-provisioned
   pyenv Python 3.14 mid-notebook.
2. **Slim bridge** — instead of converting the whole AnnData, export just
   the embedding matrix and cluster labels to disk (CSV / RDS), run
   Slingshot in a standalone R script, write pseudotime back to disk, and
   re-import into the AnnData. This sidesteps reticulate entirely because
   no AnnData object crosses the language boundary.

Path 2 is much lighter than the Q3.2 Seurat RPCA attempt. The bridge cost
of Q3.2 was largely the AnnData ↔ Seurat conversion via `zellkonverter`
that triggered `reticulate`; a Slingshot-only call needs no AnnData
conversion.

#### 3.1.4 Risk assessment
| Risk | Severity | Mitigation |
|---|---|---|
| Re-encountering the Q3.2 R env / reticulate issue | Medium (only if path 1 is chosen) | Use path 2 (disk-bridge) to avoid reticulate |
| `melanoma-r` conda env was left dormant after Q3.2 (decision log entry c "Housekeeping") | Low | Verify env still activates; rebuild if not |
| Smaller chance of methodological wrong-fit since Slingshot is mainstream | Low | Slingshot is widely used for melanoma dedifferentiation trajectories — direct fit |

**Implementation friction is the only material risk, and it is significantly
lower than the Q3.2 RPCA attempt suggested in absolute terms.**

### 3.2 PAGA (Python, in scanpy)

#### 3.2.1 Maintenance & versioning
PAGA is a first-class scanpy tool (`sc.tl.paga`). The project's
`melanoma-scrnaseq` env has scanpy 1.11.5 (runtime verified). PAGA's
canonical reference is Wolf et al. 2019 (*Genome Biology*) and it is
included in scanpy's core release cycle.

#### 3.2.2 Input requirements
PAGA operates on the **kNN connectivity graph and cluster labels** already
produced by `sc.pp.neighbors` and Leiden clustering. Our existing Harmony
pipeline (notebook 03) writes `adata.obsp['connectivities']` and notebook
04 writes `adata.obs['leiden_r0.3']` — both are direct PAGA inputs. No
raw counts required.

#### 3.2.3 Runtime sanity check (this audit)
```
$ conda run -n melanoma-scrnaseq python -c \
    "import scanpy as sc; print(sc.__version__); \
     print(callable(sc.tl.paga))"
scanpy 1.11.5
paga callable: True
paga signature: (adata, groups=None, *, use_rna_velocity=False,
                 model='v1.2', neighbors_key=None, copy=False)
```

The default `use_rna_velocity=False` confirms PAGA can run **without**
spliced/unspliced velocity input — the standard mode operates on
connectivities alone. The optional `use_rna_velocity=True` would require
scVelo upstream, but that mode is opt-in.

#### 3.2.4 Implementation effort
~1 notebook of work. Pipeline: existing h5ad → `sc.tl.paga` → connectivity
graph + per-cluster transitions → optional pseudotime via `sc.tl.dpt`
(diffusion pseudotime, also in scanpy). All on existing env, existing
data. Adds zero new dependencies.

#### 3.2.5 Risk assessment
| Risk | Severity | Mitigation |
|---|---|---|
| Continuum framing may not yield clean pseudotime if data is genuinely multi-modal | Medium — but this would itself be a finding (P1b's "continuum" framing might over-state) | Run PAGA + DPT, evaluate by Tsoi marker plotting along pseudotime |
| PAGA's discrete-cluster abstraction (still a graph of clusters) may not satisfy "continuous trajectory" framing as fully as Slingshot's principal curves | Low | Pair PAGA cluster-graph with `sc.tl.dpt` cell-level pseudotime for the continuous view |
| Stochasticity in UMAP / Leiden affects downstream PAGA structure | Low (use existing seed-locked pipeline) | Reuse the `random_state=0` convention from Stage 2 |

#### 3.2.6 Strategic significance
PAGA is the **lowest-cost path** to a trajectory result. It directly tests
[P1b]'s continuum framing using existing data, existing env, and a
2-decade-cited scanpy tool. The deliverable is a pseudotime-colored UMAP
plus a Tsoi-marker-along-pseudotime line plot — both new visual evidence
beyond Stage 3's discrete-state framing.

### 3.3 scVelo (Python)

#### 3.3.1 Maintenance & versioning
scVelo 0.3.4 (PyPI, released 2026-02-24); recent prior releases 0.3.3
(2024-12) and 0.3.2 (2024-03). Maintained by Theis Lab. Active.

#### 3.3.2 Input requirements
scVelo computes RNA velocity from the ratio of **spliced (exonic) to
unspliced (intronic)** read counts. Producing these two matrices requires:

1. Raw FASTQ reads
2. Splice-aware alignment to a reference with intron annotation
3. `velocyto` (or equivalent: `kallisto-bustools`, `STAR-solo`) to count
   reads as spliced or unspliced per cell

For Tirosh GSE72056: raw FASTQ is DUOS-gated (decision log entry a).
For Jerby-Arnon GSE115978: raw FASTQ is dbGaP-gated (this audit, §2.2).

**Neither candidate dataset's public deposit contains the spliced /
unspliced separation** scVelo needs. Both datasets place the upstream
input (raw reads) behind controlled-access portals. The blocker is
identical in class to the Stage 3 ruling on scVI / scANVI / UCE.

#### 3.3.3 Secondary concern: Smart-seq2 compatibility
Even if FASTQ access were obtained, RNA velocity on Smart-seq2 (full-length
non-UMI) is methodologically off-label — scVelo's published validations
are predominantly on 10x Chromium UMI data, where spliced / unspliced
distinction is cleaner. Smart-seq2 velocity is possible in principle but
adds a layer of methodological caveats that the Stage 3 honest-framing
rubric would have to address.

#### 3.3.4 Verdict
**BLOCKED.** Same blocker class as Stage 3 scVI / scANVI / UCE — raw
reads required, raw reads controlled-access. Document as ruled-out for
Stage 4. Goes to the same bucket as those Stage 3 methods.

### 3.4 Trajectory comparison matrix

| Axis | Slingshot | PAGA | scVelo |
|---|---|---|---|
| Algorithm family | Principal curves in reduced space | Cluster graph + connectivities | Splicing kinetics (RNA velocity) |
| Continuous pseudotime | Yes (per-cell) | Yes (via `sc.tl.dpt`) | Yes (per-cell, directed) |
| Cluster-graph output | No | Yes (PAGA's primary output) | No |
| Raw counts required | No | No | **Yes (spliced + unspliced)** |
| Existing env compatible | No (R) | **Yes (verified)** | Yes (would need pip install) |
| New deps | R env reactivation; `slingshot` BioC package | **None** | scvelo |
| Bridge friction | Disk-bridge feasible (avoids reticulate) | None | None — but data is blocked |
| Maintenance | Active | Active (scanpy core) | Active |
| Recommended path | Path 2 (disk-bridge) only | **Yes (no friction)** | Ruled out |

---

## 4. Candidate 3 — Multi-Resolution Sweep

### 4.1 Scientific question

Stage 3 [L4] flags single-resolution evaluation as a limitation: Harmony
was evaluated at `res=0.3` (6 clusters) and BBKNN at `res=0.5`
(3 clusters), each chosen by independent annotation workflows in their
respective stages. A multi-resolution sweep tests whether [P5]'s
"Undifferentiated" label divergence (33 vs 635 cells, 2-cell overlap) is
robust to resolution choice — i.e., does the label-divergence finding
generalize, or is it a function of the specific resolution selected?

### 4.2 Data / input requirements

None beyond what notebooks 03 / 05 / 07 already produce: the saved
`tirosh_malignant_annotated.h5ad` (Harmony) and
`tirosh_malignant_bbknn.h5ad` (BBKNN) AnnData files already contain
multiple leiden columns (`leiden_r0.1` … `leiden_r1.0` for Harmony; the
BBKNN h5ad has its own resolution columns).

### 4.3 Implementation effort

Re-execute notebook 07's A1/A2/B1/B2/C1/C2 metrics for each of (Harmony
res ∈ {0.1, 0.2, 0.3, 0.5, 0.8, 1.0}, BBKNN res ∈ corresponding sweep).
The C2 ARI/NMI cross-method agreement specifically becomes a function of
(harmony_res, bbknn_res) — a 2D table or heatmap. The cell-level
crosstab of `tsoi_state × tsoi_state_bbknn` similarly extends to a
per-resolution-pair analysis.

Tsoi state annotation per cluster requires manual judgment (per the Stage
2 convention — no algorithmic state assignment) at each new resolution.
This is the most labor-intensive part: ~6 resolutions × 2 methods = 12
new annotation tables.

| Step | Existing | New work |
|---|---|---|
| Multi-resolution Leiden output | Already in h5ad | None |
| Per-resolution dotplot + crosstab | Notebook 04 template | Wrap in a loop |
| Per-resolution Tsoi annotation | Notebook 04 template | Manual review × 12 |
| C2 ARI/NMI heatmap | Notebook 07 single-pair version | Generalize to 2D sweep |

### 4.4 Risk assessment

| Risk | Severity | Mitigation |
|---|---|---|
| Annotation burden at 12+ resolution combinations | Medium (this is the main cost) | Limit to 4 × 4 = 16 combinations or fewer; document selection |
| Resolution-dependent results may not generalize cleanly into the report | Low | Even if results are messy, the messiness is itself the answer to "is P5 resolution-dependent" |
| No new dependencies, no new env | n/a | n/a |

### 4.5 Strategic significance

Narrower in scope than candidates 1 or 2, but **highest feasibility
confidence**. It directly resolves L4 of the Stage 3 mini-report without
introducing any new external dependency or data. The scientific deliverable
is a paragraph in a Stage 4 supplementary report: "P5's label divergence
is / is not robust to resolution choice across the tested range." A
useful addition to the Stage 3 narrative, even if the headline novelty is
modest compared to candidate 1.

---

## 5. Cross-Candidate Comparison

| Axis | C1 (Jerby-Arnon) | C2-PAGA | C2-Slingshot | C2-scVelo | C3 (multi-res sweep) |
|---|---|---|---|---|---|
| Addresses Stage 3 finding | P1 dataset-level claim | P1b continuum framing, P5 label divergence | P1b continuum framing | (blocked) | L4, indirectly P5 |
| Scientific novelty | High — second-dataset validation | Medium — adds continuous-time view | Medium — same as PAGA | n/a | Low — robustness check |
| Data acquisition cost | Download GSE115978 supp files (public) | None | None | Blocked at FASTQ access | None |
| New dependency cost | None | None (verified) | R env reactivation | scvelo + raw-FASTQ pipeline | None |
| Implementation cost | ~Stage 2 redo (5-6 notebooks) | ~1 notebook | ~1 notebook + R bridge | n/a | ~1 notebook of loops + 12 manual annotations |
| Risk of being blocked | Low | Very low | Low-medium (R bridge re-priced) | **High — same as Stage 3 raw-count class** | Very low |
| Strategic fit | CLAUDE.md Stretch — explicitly named | Direct continuation of Q3.4 mini-report's tension | Same as PAGA | n/a | Cleanup item for the Stage 3 report's own L4 |

---

## 6. Recommendation

This audit does not select a winner; the choice depends on Stage 4 scope
intent (single deep extension vs. multiple shallow ones), available time
budget, and portfolio framing preferences. Three honest framings are
offered for decision support, not as a default.

### Framing A — Single deep extension: **Cross-dataset validation only**
Choose C1 only. Strongest scientific lever for the Stage 3 P1 headline
claim. Highest implementation cost of the three viable candidates (~Stage
2 redo), but no new dependencies and no data-access friction. The
Stage 4 deliverable is a second mini-report or a unified Stage 3+4 report
comparing Tirosh vs Jerby-Arnon four-state recovery side by side.

### Framing B — Two complementary extensions: **C1 + C2-PAGA**
Choose C1 plus PAGA. C1 validates P1; PAGA reframes P1b's continuum
observation into an actual continuous-time deliverable. PAGA is
incremental (~1 notebook, no new deps, already-verified import), so the
marginal cost on top of C1 is small. The combined deliverable tells a
more complete story: cross-dataset robustness + continuum-vs-discrete
framing resolution.

### Framing C — Three shallow extensions: **C1 + C2-PAGA + C3**
Choose all three viable candidates. C3 is cheap and cleanly closes L4 of
the Stage 3 mini-report. Total Stage 4 work expands by ~25% over Framing
B for the addition of C3.

### Methods to *not* attempt in Stage 4
- **C2-scVelo** — same blocker class as Stage 3 scVI / scANVI / UCE.
  Document in `decision_log.md` as Stage 4 ruled-out, following the same
  pattern as Stage 3 entries a / b / f. The audit verification record is
  this document.
- **C2-Slingshot** — viable but adds R-bridge re-pricing risk that PAGA
  avoids entirely. Recommend only if a methodological reason demands
  Slingshot's principal-curve formalism specifically over PAGA's
  cluster-graph formalism. The Q3.2 R-bridge experience (decision log
  entry c) suggests choosing PAGA when both are scientifically equivalent
  for the question at hand.

---

## 7. Critical Open Questions (resolve before coding)

| # | Question | Where to answer |
|---|---|---|
| Q1 | Jerby-Arnon GSE115978 — what is the exact column-name schema in the public annotations CSV (cell type column, sample/tumor identifier, malignant flag equivalent)? | Inspect the downloaded CSV before notebook 01-equivalent is written. |
| Q2 | The 31 vs 33 tumor discrepancy between GEO record and published paper — is the public count matrix scoped to all 33 or only 31? Does this matter for our purposes? | Verify at download time and document. |
| Q3 | If PAGA is adopted, does the "Other (batch?)" / "Ambiguous" catch-all category from Stage 3 need explicit handling in the cluster-graph (drop them, label them, or keep them as a node)? | Decide before running PAGA. |
| Q4 | Multi-resolution sweep — how many resolutions are enough to make the robustness claim convincingly? 4 × 4 = 16 pairs? 6 × 6 = 36? Resolution sweep design decision. | Decide before running. |
| Q5 | Should scVelo's ruling-out be recorded as a fresh decision log entry (e.g., 2026-05-26 (g)) or rolled into a "Stage 4 audit" entry that cites this document? | Methodological-record decision. |
| Q6 | Is Stage 4 a separate report or a continuation of Stage 3's mini-report? Affects framing of deliverables. | Author decision. |

---

## 8. References

- Tirosh et al. 2016 (*Science*): DOI 10.1126/science.aad0501
- Tsoi et al. 2018 (*Cancer Cell*): DOI 10.1016/j.ccell.2018.03.017
- Jerby-Arnon et al. 2018 (*Cell*): DOI 10.1016/j.cell.2018.09.006
- Balderson et al. 2024 (*BFG*): DOI 10.1093/bfgp/elad055
- Slingshot (Street et al. 2018, *BMC Genomics*): <https://bioconductor.org/packages/release/bioc/html/slingshot.html>
- PAGA (Wolf et al. 2019, *Genome Biology*): <https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.paga.html>
- scVelo (Bergen et al. 2020, *Nature Biotechnology*): <https://pypi.org/project/scvelo/>
- GSE115978 GEO record: <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE115978>

### Internal cross-references
- `docs/stage3_report.md` (P1, P1b, P5, L4)
- `docs/decision_log.md` (entries a–f, especially a for raw-count DUOS, c for R-bridge experience, f for runtime-test mandate)
- `docs/stage3_feasibility_audit.md` (Stage 3 audit precedent; this document follows the same structure)
- `environment.yml` (current `melanoma-scrnaseq` env composition)
- `notebooks/07_method_comparison.ipynb` (Q3.4 evaluation pipeline that C3 would extend)
