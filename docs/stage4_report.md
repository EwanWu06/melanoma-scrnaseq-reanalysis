# Cross-Dataset Test of Stage 3 P1: From "Dataset-Level NGFR Limitation" to a Cohort-Composition Mechanism

## 1. Abstract

Stage 3 reached conclusion P1: NGFR signal is broadly weak across Tirosh
2016 malignant cells, leaving the two NGFR-dependent Tsoi states
(Transitory, Neural crest-like) unrecoverable, and this limitation was
dataset-level rather than method-level (both Harmony and BBKNN failed
identically). Stage 4 tests P1 directly against Jerby-Arnon 2018
(GSE115978, 7,186 cells). A key fact surfaced during testing: the
dataset's Cohort column splits cells into Tirosh (4,199, 58%) and New
(2,987, 42%), where the Tirosh portion is a re-annotation of Tirosh
2016 rather than independent data. The analysis was accordingly
restructured into a three-layer framework, of which this report
completes Layers 1 and 2. Layer 2 (New cohort, 860 malignant cells)
does recover Neural crest-like (72) and Transitory (28), but 90.7% of
its NGFR-positive cells come from a single tumor, Mel110. Layer 1 (the
same Tirosh tumors, re-derived from raw counts, 1,158 cells) again
yields only Melanocytic and Undifferentiated, no NC/Transitory. P1 is
therefore refined, not refuted: the limitation is not uniform NGFR
suppression but the absence, in the Tirosh cohort, of a dedifferentiated
tumor in which NGFR is elevated and MITF/SOX10 collapse together.

## 2. Background

Tsoi et al. (2018) defined four melanoma malignant states along a
dedifferentiation axis using four markers (MITF / SOX10 / NGFR / AXL):
Melanocytic (MITF↑ SOX10↑), Transitory (adds NGFR↑), Neural crest-like
(MITF↓, NGFR↑ AXL↑), and Undifferentiated (MITF↓ SOX10↓ AXL↑).
Identification of the two intermediate states — Transitory and Neural
crest-like — both depend on elevated NGFR.

Stage 2 (Harmony + Leiden) recovered only Melanocytic and
Undifferentiated on Tirosh's 1,257 malignant cells. Stage 3 compared
Harmony with BBKNN, two structurally different integration methods, on
the same subset; both recovered only the same two states and both
failed at NC/Transitory. Because the failure was identical across two
mechanisms, Stage 3 elevated the limitation from a Harmony-specific
ceiling to a dataset-level NGFR signal limitation (P1) — a conclusion a
single-method analysis could not reach.

Stage 4's motivation is to test P1's "dataset-level" framing: apply the
analysis to a second melanoma scRNA-seq dataset potentially carrying
stronger NGFR signal, and check whether NC/Transitory become
recoverable. Jerby-Arnon 2018 (GSE115978) is the natural candidate.
Methods requiring raw counts or splicing information (scVI / scANVI /
UCE / scVelo) are excluded because raw counts or FASTQ are
controlled-access; the corresponding rulings are recorded in
`docs/decision_log.md` (Stage 3 entries a–b and Stage 4 entry g).

## 3. Methods

### 3a. Cohort discovery and the three-layer framework

GSE115978 contains 7,186 cells and 23,686 genes. The feasibility audit
originally treated it wholesale as an "independent dataset" for
cross-dataset validation. Schema inspection revealed that the Cohort
column splits cells into Tirosh (4,199, 58%) and New (2,987, 42%): the
Tirosh portion's cell IDs follow Tirosh naming conventions (cy79_*,
CY88_*) and its samples values map one-to-one to Tirosh patient
identifiers (Mel79 ↔ cy79), confirming this portion is a re-annotation
of Tirosh 2016, not an independent acquisition. This invalidated the
original "independent validation" framing (see
`docs/stage4_dataset_schema.md` §4 and decision_log entry h) and
restructured Stage 4 C1 into three layers: Layer 1 (Tirosh-subset,
annotation/processing control), Layer 2 (New cohort, the genuinely
independent test of P1), and Layer 3 (joint 7,186 cells, cross-cohort
heterogeneity). This report completes Layers 1 and 2; Layer 3 and C2
(PAGA trajectory inference) are future work.

### 3b. Layer 2 (New cohort)

Filtering to Cohort = New and cell.types = Mal gives 860 malignant
cells across 11 patients. The pipeline mirrors Stage 2:
normalize_total(1e6) + log1p, HVG (flavor='seurat', 2,000 genes),
scale, PCA (50), Harmony (batch_key='samples', θ=2), neighbors (15,
n_pcs=30), Leiden multi-resolution sweep; working resolution r=0.1 (6
clusters). Because this dataset's normalization scheme differs from
Tirosh's log2(TPM/10+1), marker thresholds were recalibrated per
dataset, and cross-dataset comparison relies on "% expressing" and
distribution shape rather than direct mean comparison.

### 3c. Layer 1 (Tirosh-subset, cluster-level comparison)

Filtering to Cohort = Tirosh and cell.types = Mal gives 1,158 malignant
cells across 12 patients — the same Tirosh tumors corresponding to
Stage 2's 1,257 cells, but here reconstructed and normalized from
GSE115978's raw counts. The pipeline matches Layer 2; to avoid the
z-score comparability problem from Layer 2 (where AXL entered the HVG
set and was scaled), all four marker statistics were computed on
un-scaled normalized data. Comparison is cluster-level rather than
cell-level: the GSE72056 and GSE115978 barcode naming schemes are
incompatible, with direct cell-ID overlap only 366/1,158 (~32%), too
low for cell-by-cell mapping; patient-level correspondence (Mel79 ↔ 79)
is exact and reliable, so Tsoi state distributions are compared at the
cluster/cohort level.

### 3d. Normalization control (dual-source)

To rule out normalization scheme or source matrix as the cause of
NC/Transitory failure, NGFR-positive rate was compared across two
independent sources of Tirosh malignant cells: GSE72056
(log2(TPM/10+1), n=1,257) gives 10.0%, distributed across 9 patients
with top patient 27.0%; GSE115978 (raw counts + normalize+log1p,
n=1,158) gives 8.5%, also across 9 patients with top patient 31.6%. The
leading shared NGFR-contributing patients are Mel79/89/81 (with Mel53
positive in both sources). This agreement indicates NGFR signal
characteristics are robust to source and normalization.

## 4. Results

### [R1] Layer 2 recovers NC and Transitory, but they are Mel110-driven

At r=0.1 (6 clusters), Layer 2's Tsoi distribution is: Melanocytic 376
(cl1+cl2), Ambiguous-Mel 342 (cl0), Neural crest-like 72 (cl3),
Ambiguous 42 (cl4), Transitory 28 (cl5), totaling 860. Unlike Stage
2/Stage 3, which recovered only two states, both NGFR-dependent
intermediate states appear here. The NC cluster (cl3, n=72) has NGFR
positivity 61.1% (44/72), MITF only 26.4%, SOX10 only 2.8% — a marker
profile matching Tsoi's NC definition; the Transitory cluster (cl5,
n=28) has NGFR 60.7% (17/28), MITF 53.6%, SOX10 10.7%.

Patient-level inspection, however, shows this recovery is single-tumor
dominated: the NC cluster is 75.0% (54/72) from Mel110, and Transitory
is 71.4% (20/28) from Mel110. Combining both states into an NGFR-high
group (n=100), Mel110 contributes 74%, and adding Mel106 covers nearly
all of it. Taking all NGFR-positive cells (n=86) as the denominator,
Mel110 accounts for 90.7%. In other words, Layer 2 can recover NGFR
states not because the New cohort is broadly NGFR-strong, but because
it happened to sample Mel110, a dedifferentiated tumor.

### [R2] Layer 1 still fails to recover NC/Transitory on the same Tirosh tumors

Re-running the identical pipeline from raw counts on the same Tirosh
tumors (1,158 cells), r=0.1 yields 7 clusters with Tsoi distribution
Melanocytic 1,069 and Undifferentiated 89 (cl5) — no Neural crest-like
or Transitory, consistent with Stage 2 on GSE72056. The highest-NGFR
cluster (cl6, n=86) has NGFR positivity only 27.9% (24/86), but retains
MITF 76.7% and SOX10 81.4%, and is therefore annotated Melanocytic
rather than NC (true NC requires MITF/SOX10 collapse). cl6 is dominated
by Mel81 (81/86 = 94.2%).

### [R3] Mel81 versus Mel110: the distinction is dedifferentiation, not NGFR itself

Placing Layer 1's Mel81 beside Layer 2's Mel110, both are
"single-tumor-dominated NGFR-elevated clusters" (Mel110 75–91% of its
cluster/NGFR+ pool, Mel81 83–94%), but their biological meaning differs
sharply:

| | Layer 2 Mel110 (NC cl3) | Layer 1 Mel81 (cl6) |
|---|---|---|
| Single-patient share | 75.0% (cl3) / 90.7% (all NGFR+) | 94.2% (cl6) / 83.3% (NGFR+ within cl6) |
| NGFR % | 61.1% | 27.9% |
| MITF % | 26.4% | 76.7% |
| SOX10 % | 2.8% | 81.4% |
| Dedifferentiated? | Yes → true NC/Transitory | No → Melanocytic |

This refines the mechanism. Tirosh is not simply "missing an NGFR-rich
tumor" (it has Mel81, an NGFR-elevated tumor); it is missing a tumor in
which NGFR is elevated and MITF/SOX10 simultaneously collapse — a
dedifferentiated tumor. NGFR elevation alone (Mel81) stops at
Melanocytic; dedifferentiation (Mel110) produces NC/Transitory.

### [R4] Normalization and source matrix ruled out

The dual-source NGFR-positive rates (GSE72056 10.0% vs GSE115978 8.5%)
and patient distributions (both 9 patients, top only ~27–32%, matching
contributor identities) are highly consistent, indicating the
non-recoverability of NC/Transitory is independent of normalization
scheme and source matrix, and further confirming P1's reading that NGFR
is spread across multiple patients with no single outlier driver.

## 5. Discussion

### 5a. Refining P1: from dataset-level to cohort composition

P1 holds under patient-level scrutiny: Tirosh's NGFR-positive cells are
spread across 9 patients with the top patient at only ~32%, not driven
by 1–2 outliers, so the cohort-level conclusion that "NGFR is
insufficient to form NC/Transitory clusters" is correct. Layer 2
contributes a mechanistic refinement: the limitation is not uniform
NGFR suppression across all cells but cohort composition — Tirosh lacks
a dedifferentiated tumor capable of forming an NGFR-defined cluster. By
contrast, Layer 2's Mel110 provides exactly such a tumor; Tirosh's
Mel53, though extremely NGFR-positive, has only 14 cells, too few to
form a stable cluster. Whether NGFR states can be recovered at the
cohort level depends on whether sufficient NGFR-rich dedifferentiated
tumors are sampled.

### 5b. NGFR states are tumor-specific

Taking Layers 1 and 2 together, NGFR-defined Tsoi states in melanoma
scRNA-seq appear tumor-specific (concentrated in individual tumors)
rather than evenly distributed across a cohort. This explains why
different studies report inconsistent recovery of these intermediate
states: recoverability depends largely on which tumors were sampled.
This hypothesis is presently supported by only two tumors (Mel53,
Mel110) and is hypothesis-generating (see §6 L1).

### 5c. Methodological note

Layer 2's conclusion went through several revisions (Phase 1.2 → 2 →
2.6 → 2.8). At Phase 2.6 a result appeared — the NC cluster was 100%
post-treatment, χ² = 36.58, p < 0.0001 — and was briefly read as
"treatment-induced dedifferentiation." Phase 2.7 sanity checks
(patient breakdown, the fully confounded structure in which every
patient is 100% one treatment arm, and the fact that the Tirosh cohort
actually contains 43% post-treatment cells) showed this was
patient/treatment confounding rather than a biological signal. This
trajectory is recorded as transparent research practice in decision_log
entry (i); it is not developed as a main thread, noting only why the
conclusion ultimately rests on patient/cohort composition rather than a
treatment effect.

## 6. Limitations

**L1 — Only two informative tumors.** Essentially all NGFR-state signal
derives from two tumors: Mel110 (dedifferentiated → NC/Transitory) and
Mel81 (NGFR-elevated but still Melanocytic). Layer 2's NGFR-high group
(n=100) depends on Mel110 + Mel106. n=2 is hypothesis-generating, not
confirmatory.

**L2 — Single-patient / treatment confounding.** All 11 Layer 2
patients are 100% one treatment arm (no within-patient naive/post
pairing), so treatment effect and patient identity are mathematically
inseparable; the 12 Tirosh-subset patients are likewise 100% one-way.
Whether treatment participates in NGFR-state formation can be neither
established nor excluded.

**L3 — cluster 6 batch caveat.** Layer 1's cl6 is 94.2% from a single
patient (Mel81), and may carry patient-specific residual batch signal
(analogous to Stage 2's "Other (batch?)" catch-all). But whether cl6 is
genuine Mel81 biology or residual batch, the conclusion is unchanged:
it has no dedifferentiation signature (MITF/SOX10 retained) and does not
constitute recovered NC/Transitory.

**L4 — Cell-ID mapping not viable.** The GSE72056 and GSE115978
cell-barcode naming schemes are incompatible, with direct overlap only
366/1,158 (~32%), precluding cell-level verification; Layer 1 compares
only at the cluster level. Patient-level correspondence is reliable,
but cell-level identity alignment is unavailable.

**L5 — Deferred work.** Layer 3 (joint cross-cohort analysis), C2
(PAGA trajectory inference), and a definitive test on a larger
multi-patient NGFR-rich cohort are all deferred; this report is scoped
to Layers 1 and 2.

## 7. Conclusions and Future Directions

Stage 4 C1 refines Stage 3 P1: the failure to recover NC/Transitory in
Tirosh is confirmed at the cohort level and survives patient-level
scrutiny, but its mechanism is cohort composition — specifically, the
absence of a dedifferentiated tumor (NGFR-high with MITF/SOX10
collapse) — rather than uniform NGFR suppression, method ceiling,
normalization scheme, or source matrix. NGFR-defined Tsoi states appear
to be tumor-specific dedifferentiation phenomena whose cohort-level
recoverability depends on tumor sampling.

Future directions: Layer 3 (joint cross-cohort analysis to test whether
the Cohort axis is batch or biology), C2-PAGA (continuous-pseudotime
characterization of the dedifferentiation axis, addressing Stage 3
P1b), and a larger multi-patient cohort enriched for NGFR-high
dedifferentiated tumors for a definitive test of the tumor-specificity
hypothesis.

---

*Author: Ewan Wu (Youyou Wu)*  
*Source artifacts: `notebooks/08_jerby_arnon_load.ipynb` (Layer 2), `notebooks/09_layer1_tirosh_subset.ipynb` (Layer 1), `docs/decision_log.md` (entries g, h, i), `docs/stage4_dataset_schema.md`, `docs/stage3_report.md` (P1)*
