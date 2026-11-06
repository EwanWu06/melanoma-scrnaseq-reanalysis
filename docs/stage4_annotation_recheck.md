# Tsoi-State Annotation Recheck (2026-11-05)

Decision log: entry (j). Script: `scripts/recheck_tsoi_annotation.py`
(read-only on `data/processed/`; no notebook or processed file was changed).

## 1. Why

Two inconsistencies were flagged on 2026-11-03:

1. The project's own Tsoi table (`docs/literature/tsoi_2018.md`) defines
   both intermediate states with **SOX10↑**: Transitory = MITF↑ SOX10↑ NGFR↑
   AXL↓; Neural crest-like = MITF↓ SOX10↑ NGFR↑ AXL↑. Tsoi et al. (2018)
   distinguish the Undifferentiated subtype from the Neural crest-like
   subtype precisely by lower NGFR and SOX10. Yet the Stage 4 Layer 2
   clusters labeled Neural crest-like (cl3) and Transitory (cl5) had SOX10
   detection of only 2.8% and 10.7%, and the Stage 4 report and notebook
   09 stated the opposite rule ("true NC requires MITF/SOX10 collapse").
2. Single-patient handling was inconsistent. Stage 2 set aside its
   highest-NGFR cluster (cl3, 78% patient 80) as a batch artifact, while
   Stage 4 accepted an equally single-patient cluster (Layer 2 cl3, 75%
   Mel110) as a biological state.

## 2. Rules (written before the new computations; see caveat in §6)

The original clusterings are reused, not recomputed: Stage 2 Harmony
`leiden_r0.3`, Stage 3 BBKNN `leiden_bbknn_r0.5`, Stage 4 Layer 1 and
Layer 2 `leiden_r0.1`. In each case this is the resolution at which every
cluster carries a single final label.

**Cluster level.** For each cluster, compute the detection rate (fraction
of cells with value > 0) of MITF, SOX10, NGFR and AXL, then call each
marker's direction two ways:

- **Rule 1 (relative, primary):** ↑ if ≥ 1.25 × the analysis-wide
  malignant detection rate; ↓ if ≤ 0.75 ×; otherwise "~". A "~" can
  match either direction.
- **Rule 2 (absolute):** ↑ if ≥ 50%, else ↓.

A cluster matches a state only if all four markers agree with the table.
AXL is detected in most cells in every dataset and barely discriminates,
so a three-marker match (ignoring AXL) is also reported. A cluster is
flagged **single-patient** if one patient contributes ≥ 75% of its cells.
The same threshold is applied to every dataset.

**Cell level (descriptive).** Each cell is placed in one of three marker
patterns using detection only:
- Transitory pattern: NGFR+ SOX10+ MITF+
- NC pattern: NGFR+ SOX10+ MITF−
- NGFR+ SOX10−

## 3. Cluster-level results

Detection rates are in %. "Rule 1" is the four-marker match; "Rule 2" is
the three-marker match (Rule 2 with all four markers matches nothing in
the Tirosh data, because AXL is ≥ 50% in every cluster). ● marks a
single-patient cluster. \*Melanocytic/Undifferentiated ties are an
artifact of Rule 1: when base detection is above 80%, ×1.25 exceeds
100%, so MITF and SOX10 can never be called ↑ (see §6).

**Stage 2, Harmony (GSE72056, 1,257 cells).** Analysis-wide detection:
MITF 84.6, SOX10 83.1, NGFR 10.0, AXL 72.2.

| cl | n | original label | MITF | SOX10 | NGFR | AXL | Rule 1 | Rule 2 | top patient |
|---|---|---|---|---|---|---|---|---|---|
| 0 | 451 | Melanocytic | 93.6 | 86.9 | 6.0 | 77.4 | Mel/Undiff\* | Mel | pt88 25% |
| 1 | 69 | Melanocytic | 89.9 | 84.1 | 7.2 | 65.2 | Mel/Undiff\* | Mel | pt80 36% |
| 2 | 179 | Other (batch?) | 54.2 | 69.8 | 24.6 | 69.3 | NC | Mel | pt81 73% |
| 3 | 107 | Other (batch?) | 84.1 | 86.0 | 19.6 | 60.7 | T/NC | Mel | pt80 78% ● |
| 4 | 33 | Undifferentiated | 60.6 | 81.8 | 21.2 | 81.8 | NC | Mel | pt79 21% |
| 5 | 418 | Other (batch?) | 89.0 | 84.0 | 5.3 | 71.3 | Mel/Undiff\* | Mel | pt79 100% ● |

**Stage 3, BBKNN (GSE72056, 1,257 cells).** Same analysis-wide rates.

| cl | n | original label | MITF | SOX10 | NGFR | AXL | Rule 1 | Rule 2 | top patient |
|---|---|---|---|---|---|---|---|---|---|
| 0 | 404 | Melanocytic | 93.6 | 87.1 | 5.4 | 75.0 | Mel/Undiff\* | Mel | pt78 24% |
| 1 | 635 | Undifferentiated | 78.1 | 78.7 | 11.7 | 72.3 | any (all "~") | Mel | pt79 68% |
| 2 | 218 | Ambiguous | 86.7 | 88.5 | 13.8 | 67.0 | T/NC | Mel | pt80 44% |

**Stage 4 Layer 1 (same Tirosh tumors, GSE115978 counts, 1,158 cells).**
Analysis-wide: MITF 85.0, SOX10 82.9, NGFR 8.5, AXL 72.7.

| cl | n | original label | MITF | SOX10 | NGFR | AXL | Rule 1 | Rule 2 | top patient |
|---|---|---|---|---|---|---|---|---|---|
| 0 | 348 | Melanocytic | 90.8 | 84.2 | 6.3 | 65.2 | Mel/Undiff\* | Mel | Mel79 100% ● |
| 1 | 283 | Melanocytic | 91.2 | 86.2 | 1.8 | 82.3 | Mel/Undiff\* | Mel | Mel88 34% |
| 2 | 159 | Melanocytic | 85.5 | 86.2 | 17.6 | 79.9 | T/NC | Mel | Mel80 36% |
| 3 | 100 | Melanocytic | 67.0 | 71.0 | 6.0 | 46.0 | Mel | Mel | Mel79 37% |
| 4 | 93 | Melanocytic | 94.6 | 89.2 | 5.4 | 76.3 | Mel/Undiff\* | Mel | Mel78 43% |
| 5 | 89 | Undifferentiated | 59.6 | 69.7 | 9.0 | 80.9 | NC/Undiff | Mel | Mel79 44% |
| 6 | 86 | Melanocytic | 76.7 | 81.4 | 27.9 | 76.7 | T/NC | Mel | Mel81 94% ● |

**Stage 4 Layer 2 (independent New cohort, 860 cells).** Analysis-wide:
MITF 64.1, SOX10 49.7, NGFR 10.0, AXL 94.8.

| cl | n | original label | MITF | SOX10 | NGFR | AXL | Rule 1 | Rule 2 | top patient |
|---|---|---|---|---|---|---|---|---|---|
| 0 | 342 | Ambiguous-Mel | 49.1 | 34.5 | 2.3 | 88.3 | Undiff | Undiff | Mel106 22% |
| 1 | 271 | Melanocytic | 84.5 | 76.8 | 5.2 | 100.0 | Mel | Mel | Mel103 32% |
| 2 | 105 | Melanocytic | 94.3 | 84.8 | 1.0 | 100.0 | Mel | Mel | Mel102 89% ● |
| 3 | 72 | Neural crest-like | 26.4 | **2.8** | 61.1 | 98.6 | **none** | **none** | Mel110 75% ● |
| 4 | 42 | Ambiguous | 50.0 | 16.7 | 4.8 | 100.0 | Undiff | none | Mel102 40% |
| 5 | 28 | Transitory | 53.6 | **10.7** | 60.7 | 85.7 | **none** | **none** | Mel110 71% |

## 4. Cell-level results

| dataset | Transitory pattern | NC pattern | NGFR+ SOX10− |
|---|---|---|---|
| GSE72056 (Stage 2/3, 1,257) | 71 (5.6%), 5 patients ≥ 3 cells, top pt89 30% | 35 (2.8%), 4 patients, top pt79 37% | 20 (1.6%), 3 patients, top pt53 35% |
| Layer 1 (1,158) | 58 (5.0%), 3 patients, top Mel89 34% | 25 (2.2%), 3 patients, top Mel79 48% | 15 (1.3%), 3 patients, top Mel53 40% |
| Layer 2 (860) | 3 (0.3%), no patient ≥ 3 | 3 (0.3%), no patient ≥ 3 | **80 (9.3%), 77 from Mel110 (96%)** |

## 5. Findings

1. **The Layer 2 "Neural crest-like" and "Transitory" labels are
   withdrawn.** Under both rules, with three or four markers, neither
   cluster matches any Tsoi state. Both have lost SOX10 (2.8% and 10.7%,
   against 77–85% in the same cohort's Melanocytic clusters and 49.7%
   cohort-wide). Only 6 of 860 cells (0.7%) show an intermediate-state
   marker pattern. What the independent cohort contains is an
   **NGFR+ SOX10− population, 96% from one tumor (Mel110)**, which matches
   none of the four states.
2. **Cells with the intermediate-state marker pattern do exist in Tirosh
   2016, but they are sparse and scattered.** About 8% of GSE72056
   malignant cells (106) show the Transitory or NC pattern. They are spread
   across 4–5 patients, and no patient contributes more than 37%. This
   explains why they do not form discrete clusters, and it is consistent
   with Balderson et al. (2024) recovering all four states as positions
   along a continuous trajectory.
3. **Whether the Tirosh data contain intermediate-state *clusters* depends
   on an NGFR threshold that the original analysis never fixed.** Under
   Rule 2, no cluster qualifies, because NGFR detection is ≤ 28% in every
   cluster. Under Rule 1, several clusters qualify as Transitory/NC
   candidates. Some are single-patient (Stage 2 cl3, Layer 1 cl6) or nearly
   so (Stage 2 cl2, 73%). Others are cross-patient: Stage 2 cl4, BBKNN
   cl2, and Layer 1 cl2, each with NGFR detection of only 14–21%.
4. **The Stage 2/3 "Undifferentiated" labels do not meet the full
   definition.** Stage 2 cl4 has SOX10 detection of 81.8% (about the
   cohort average) and NGFR of 21.2% (twice the average). The table
   requires both to be low. The original label rested on MITF-low and
   AXL-high mean expression only. BBKNN's 635-cell "Undifferentiated"
   cluster is indistinguishable from the cohort average on all four
   detection rates. In the independent cohort, by contrast, a
   cross-patient cluster (Layer 2 cl0) does match Undifferentiated under
   both rules.
5. **Single-patient handling, made consistent.** Under the same
   ≥ 75% rule, Stage 2 cl3 (patient 80) and Layer 2 cl3 (Mel110) are both
   single-patient clusters and should be treated alike.

**Revised reading of P1.** The Tirosh data contain no cluster with strong
NGFR signal under either integration method; this still holds. The
mechanism proposed in Stage 4 is withdrawn: Tirosh does not "lack a
dedifferentiated tumor whose recovery appears in the independent cohort."
The independent cohort's NGFR-high cells are not intermediate states.
Tirosh's intermediate-pattern cells exist but are too sparse and too
dispersed across patients to form discrete clusters. The four-state
annotation in all stages depended on marker thresholds and single-patient
rules that were never written down. Re-annotating every stage under one
written rule set is the open item (entry j).

## 6. Limitations of this recheck

- **Detection only.** The recheck uses detection (> 0) and ignores
  expression level. Stage 2 annotated on mean expression, and the two
  views can disagree.
- **A flaw in Rule 1.** For markers already detected in more than 80% of
  cells (MITF and SOX10 in the Tirosh data), Rule 1 can never call ↑.
  These markers therefore never discriminate Melanocytic from
  Undifferentiated there. This flaw was not anticipated when the rules
  were written.
- **Dropout.** Cell-level patterns are affected by dropout and are
  descriptive only. The Layer 2 SOX10 loss, however, is a cluster-level
  contrast (2.8% vs 77–85%) and does not depend on single-cell dropout.
- **Not blind.** The rules were written after the Layer 1 and Layer 2
  per-cluster detection rates had been seen (notebooks 08 and 09). The
  Stage 2 and Stage 3 per-cluster rates and all cell-level counts were
  computed for the first time here.
