# Stage 4 Dataset Schema — Jerby-Arnon 2018 (GSE115978)

*Stage 4 preparation · 2026-05-26*
*Documents discovery during Boundary B (Phase 1+2 of Stage 4 C1 implementation prep)*

## 1. GEO Accession & Download

Source: NCBI GEO accession GSE115978
URL: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE115978
FTP: https://ftp.ncbi.nlm.nih.gov/geo/series/GSE115nnn/GSE115978/suppl/

Files downloaded (per Stage 4 audit §2.2 + this preparation):

| File | Compressed | MD5 |
|---|---|---|
| GSE115978_counts.csv.gz | 49.1 MB | cd4d37fdb7e9c51a64f51222027a7008 |
| GSE115978_cell.annotations.csv.gz | 98.9 KB | 066fa72772ae0b2cf065f03aa9c89a05 |

TPM CSV intentionally skipped: we re-normalize from counts using 
sc.pp.normalize_total + sc.pp.log1p (standard scanpy convention). 
TPM CSV remains available at the GEO accession for users requiring 
byte-level comparison with the original publication.

Checksums tracked at `data/jerby_arnon_CHECKSUMS.md5`. Actual data 
files at `data/raw/jerby_arnon/*.csv.gz` (gitignored).

## 2. Annotation Schema (Audit Q1 Resolved)

Cell annotations file (`GSE115978_cell.annotations.csv.gz`): 7186 rows × 7 columns.

| Column | Type | Notes |
|---|---|---|
| `cells` | object | Cell ID (e.g., `cy78_CD45_neg_1_B04_S496_comb`) |
| `samples` | object | Sample identifier (e.g., `Mel79`, `Mel129pa`) |
| `cell.types` | object | Author-curated cell type (10 unique values) |
| `treatment.group` | object | naive / post-treatment (drug exposure context) |
| `Cohort` | object | `Tirosh` or `New` — source dataset |
| `no.of.genes` | int64 | Per-cell gene count (QC) |
| `no.of.reads` | int64 | Per-cell read count (QC) |

`cell.types` distribution (10 unique values):

```
Mal           2018   (malignant)
T.CD8         1759
T.CD4          856
B.cell         818
T.cell         706   (likely unsubtyped T)
Macrophage     420
?              307   (Jerby-Arnon's uncertain bucket, 4.3%)
CAF            106
Endo.          104
NK              92
```

**Malignant filter**: `cell.types == 'Mal'` (2018 cells). 
No separate boolean — single-column convention differs from 
Tirosh's `malignant_status` (0/1/2) scheme.

## 3. Tumor/Sample/Patient Counts (Audit Q2 Resolved)

| Granularity | Count | Source |
|---|---|---|
| Distinct `samples` values | 32 | This dataset |
| Distinct patients (unique donors) | 31 | Mel129pa + Mel129pb collapsed |
| Published paper headline | 33 | Jerby-Arnon 2018 (*Cell*) |

`Mel129pa` and `Mel129pb` are almost certainly two timepoints from 
patient `Mel129` (`pa`/`pb` = pre/post-treatment, consistent with 
the `treatment.group` column having both `naive` and 
`post-treatment` values). The 1-tumor discrepancy between this 
public deposit (32 samples / 31 patients) and the paper's 33 figure 
likely reflects 1 patient not included in the GEO supplementary 
files.

Recommendation for Stage 4 notebook: report both granularities 
(31 patients, 32 samples) and note the 33 paper discrepancy.

## 4. Cohort Composition — Critical Finding

**Audit assumption (incorrect)**: Stage 4 audit §2.1 and §2.6 framed 
Jerby-Arnon as "an independent dataset" for cross-dataset validation 
of Stage 3 P1.

**Actual cohort structure**:

The `Tirosh` cohort consists of cells from Tirosh 2016 (GSE72056) 
re-annotated by the Jerby-Arnon lab. Evidence: cell IDs follow 
Tirosh naming convention (e.g., `cy78_*`, `cy79_*`, `CY88_*`), 
and `samples` values in this cohort match Tirosh patient identifiers 
(`Mel79` ↔ `cy79`, `Mel75` ↔ `cy75`, etc.).

This invalidates the audit's "independent cross-dataset validation" 
framing for Jerby-Arnon as a whole. See `docs/decision_log.md` 
entry (h) for the corrected framing.

## 5. Revised Stage 4 C1 Framework — Three Layers

Stage 4 C1 (Jerby-Arnon analysis) is redesigned as 3 layered 
comparisons rather than a single "Tirosh vs Jerby-Arnon" binary:

### Layer 1 — Jerby-Arnon-Tirosh-subset (4199 cells)
**Same cells, two annotation frameworks.** Compare Stage 2's 
Tsoi-state annotation (our work) with Jerby-Arnon's `cell.types` 
annotation on the same cell population. Specific questions:
- What does our "Other (batch?)" / "Ambiguous" map to in 
  Jerby-Arnon's labels?
- Does our Mel/Undiff Tsoi partition match Jerby-Arnon's 
  Mal cell subdivision?
- Annotation-framework meta-comparison.

### Layer 2 — Jerby-Arnon-New (2987 cells)
**True independent validation of Stage 3 P1.** Apply Stage 2-style 
pipeline (Harmony → Leiden → Tsoi annotation) to the New cohort 
only. Test:
- Is NGFR signal stronger in this independent acquisition?
- Do Transitory / Neural crest-like states emerge with sufficient 
  NGFR expression?
- Direct test of P1's "dataset-level NGFR limitation" claim.

### Layer 3 — Joint (7186 cells)
**Cross-cohort heterogeneity within one publication.** Apply 
integration to the full Jerby-Arnon dataset and compare how 
`Cohort` partitions the integrated space. Tests:
- Does Jerby-Arnon-New cluster separately from Jerby-Arnon-Tirosh 
  in their own embedding?
- Is the `Cohort` axis (within one lab) a meaningful biological 
  signal or a technical batch?

## 6. Implementation Notes for Stage 4 Notebook(s)

- For malignant subset analyses: `adata = adata[adata.obs['cell.types'] == 'Mal']` (2018 cells)
- For the 307 `?` cells (uncertain bucket): handling deferred to per-notebook judgment (recommend: filter out for Tsoi-state work, retain for joint analyses)
- For `Mel129pa`/`Mel129pb`: treat as separate samples (preserve `treatment.group` signal); collapse only when patient-level statistics are required.
- Tirosh-subset cells map to Stage 2's 1,257 malignant subset via cell-ID intersection on the Tirosh cohort filter.
