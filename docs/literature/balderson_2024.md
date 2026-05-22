# Balderson et al. 2024 — Reading Notes

> Author collaboration note: Q3 and Q4 were initially drafted with AI
> assistance but contained generic limitations that did not apply to this
> paper. Sections were rewritten by the author after a focused re-reading
> of the Methods and Discussion.

## Citation

Balderson et al. (2024). "Systematic analysis of the transcriptional
landscape of melanoma reveals drug-target expression plasticity."
*Briefings in Functional Genomics*, elad055. DOI: 10.1093/bfgp/elad055.

## 1. What They Did
- Used **Tirosh 2016 scRNA-seq data** (malignant cells, ~1,252 cells from
  19 patients)
- Applied **unsupervised methods**: PCA for dimensionality reduction +
  Monocle for trajectory inference — unsupervised clustering and trajectory
  analysis, with no pre-defined signature scores imposed on the data.

## 2. Main Findings
- Re-analysing the Tirosh data, they directly recovered **four distinct
  malignant cell states** at single-cell resolution.
- These four states **map onto the Tsoi 2018 four-state model** — the
  melanoma differentiation axis from Melanocytic to Undifferentiated.
- Demonstrated that the Tirosh "binary" finding was a methodological
  artifact of supervised signature scoring.

## 3. Limitations of Their Approach (Our Entry Point)

Although Balderson 2024 successfully recovered four cell states from the
Tirosh 2016 dataset, the analysis relies heavily on classical
dimensionality reduction and trajectory inference (PCA + Monocle).
Classical methods of this kind — linear PCA in particular — struggle with
scRNA-seq data, which is extremely high-dimensional, highly sparse, and
governed by complex non-linear relationships. Such methods can lose deep
biological features during dimensionality reduction, yielding a less
precise representation.

Key technical limitations:
1. **PCA is linear** — assumes major variation can be captured by linear
   projection, but scRNA-seq biology is highly non-linear.
2. **No built-in batch correction** — Monocle does not inherently handle
   batch effects, and Tirosh's 19 patients are a strong batch source.
3. **Single-dataset analysis** — no cross-cohort validation of robustness.

## 4. How Our Project Differs

Our project is a methodological upgrade of — and a friendly challenge to —
Balderson 2024. We use exactly the same underlying data (Tirosh 2016
malignant cells), but rather than the traditional route we re-analyse it
with state-of-the-art deep-learning / deep generative models (scVI,
scANVI). The comparison we want to make: relative to the classical
PCA/Monocle pipeline, does a neural-network-based latent space recover the
four cell states more robustly and precisely — and can it surface data
features the classical methods failed to capture?

We build on Balderson's work by asking: **Is their four-state recovery
method-independent?** Specifically:
- Apply **variational autoencoder-based methods** (scVI, scANVI) that
  model non-linear structure and condition on batch variables
- If we recover the same four states → strengthens the four-state model
  as robust biology
- If we recover different structure → suggests the four-state model may
  have method-dependent components

## 5. Relevance to Our Project
- This is the **direct predecessor** of our work — we explicitly
  position our project as a methodological extension.
- Their result becomes our "classical baseline" for comparison.
