# Tirosh et al. 2016 — Reading Notes

> **Author collaboration note:** Primary reading and conceptual understanding
> by Ewan Wu. Notes structured with AI assistance (Claude); all biological
> and methodological interpretations are the author's own.
>
> **Last updated:** 2026-05-17
> **Read status:** ✅ Complete (Abstract, Intro, Methods, Discussion)

---

## Citation

Tirosh I, Izar B, Prakadan SM, et al. **Dissecting the multicellular ecosystem
of metastatic melanoma by single-cell RNA-seq.** *Science* 352(6282):189-196
(April 8, 2016). DOI: [10.1126/science.aad0501](https://doi.org/10.1126/science.aad0501)

---

## 1. Experimental Design and Technology

The authors profiled tumor tissues from **19 metastatic melanoma patients**
using **fluorescence-activated cell sorting (FACS)** to isolate single cells,
followed by **Smart-seq2** single-cell RNA sequencing. The final dataset
contains **4,645 single cells × 23,686 genes**, with a mixed composition of
malignant cells, immune cells, stromal cells, and endothelial cells.

Key design choice: cells were FACS-sorted by CD45 status before sequencing,
which strongly affected the immune cell representation in the final dataset.

---

## 2. Main Findings

**Malignant cell heterogeneity:**
- Confirmed within-tumor co-existence of two transcriptional programs:
  - **MITF-high state** (proliferative)
  - **MITF-low / AXL-high state** (associated with targeted therapy resistance)

**Tumor microenvironment:**
- Characterized non-malignant cell types in the TME
- Detailed analysis of **T cell exhaustion programs** and their relationship
  to activation
- Identified communication networks between cancer-associated fibroblasts
  (CAFs) and immune cells

---

## 3. Methodology — Critical Detail for Our Project

The authors used a **mixed approach** to cell state analysis:

- **Unsupervised methods** (hierarchical clustering, t-SNE) for:
  - Distinguishing malignant vs. non-malignant cells
  - Clustering immune cell subpopulations

- **Supervised methods** for malignant cell states:
  - They did **not** perform unsupervised clustering on malignant cells.
  - Instead, they applied **pre-defined gene signatures** from Hoek 2008 and
    Konieczkowski 2014 to compute MITF-program and AXL-program scores for
    each cell.
  - This is the **methodological gap** that Balderson 2024 later exploited.

---

## 4. Why This Paper Matters

This was among the first applications of scRNA-seq to a **complex human solid
tumor**, enabling simultaneous resolution of malignant clonal heterogeneity
and microenvironmental cell states at single-cell resolution. It established
the dataset (GSE72056) that has since become a foundational benchmark in
melanoma single-cell biology.

---

## 5. Relevance to Our Project

This is the **primary dataset** for our re-analysis. The key implications for
our project design are:

1. **Data format issue:** The GEO-deposited data is `log2(TPM/10+1)`, not
   raw counts — this constrains which deep learning methods we can apply.
2. **Patient-origin batch effect:** 19 patients implies strong batch
   structure that must be addressed (Harmony, scVI batch conditioning).
3. **The supervised cell state assignment is the gap:** Our project asks
   whether unsupervised deep learning methods recover the same (or different)
   cell state structure as the original supervised analysis and Balderson's
   classical re-analysis.
