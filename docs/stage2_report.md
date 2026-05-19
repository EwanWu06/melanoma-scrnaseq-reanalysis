# Stage 2 Mini-Report — Classical Pipeline Replication

*Project: Learning to Re-analyze — A Methodological Study of Deep Learning
Approaches on the Tirosh et al. (2016) Melanoma Single-Cell Dataset*
Stage 2 (Q2.1–Q2.5) · Notebooks `01`–`04`

## 1. Background

Malignant melanoma cells exhibit transcriptional state heterogeneity along a
dedifferentiation axis. Tsoi et al. 2018 (*Cancer Cell*) defined four distinct
states using four marker genes ($MITF$ / $SOX10$ / $NGFR$ / $AXL$):
Melanocytic ($MITF^{high}SOX10^{high}NGFR^{low}AXL^{low}$), Transitory
($+NGFR^{high}$), Neural crest-like ($MITF^{low}SOX10^{high}NGFR^{high}AXL^{high}$),
and Undifferentiated ($MITF^{low}SOX10^{low}NGFR^{low}AXL^{high}$). Balderson
et al. 2024 recapitulated these four states in malignant cells from the Tirosh
2016 dataset (GSE72056) using a conventional workflow (PCA + Monocle trajectory
inference). However, their analysis was conducted without any batch correction.
Given that the dataset is derived from 19 patients and features highly
imbalanced sample sizes per patient, this omission represents a glaring
methodological weakness. The research question for the current stage is: does
the Tirosh dataset still support a clean, four-state Tsoi structure under a
classic workflow that incorporates batch correction (Harmony) and an upgraded
clustering algorithm (Leiden)?

## 2. Methods

**Data Source.** GSE72056 (Tirosh et al. 2016, Smart-seq2, full-length non-UMI
protocol), comprising 4,645 cells × 23,686 genes. Expression values are
recorded as $\log_2(\text{TPM}/10+1)$, i.e. already normalized and
log-transformed. Consequently, no additional `normalize_total`, `log1p`, or
`scale` operations were applied across the entire workflow (an extra scaling
step on log-TPM data distorts inter-genic relationships). The dataset includes
19 patients with highly skewed sample sizes, ranging from a maximum of 896
cells (pt79) to a minimum of 63 cells (pt65)—an approximately 14-fold
disparity.

**Malignant Cell Subset.** To maintain strict comparability with Balderson
2024, only the 1,257 cells explicitly labeled as malignant were retained,
directly adopting the original Tirosh malignant annotations. Notably, these
1,257 cells span only 15 of the 19 patients (4 patients had zero malignant
cells captured).

**Highly Variable Genes (HVGs).** HVGs were identified using
`sc.pp.highly_variable_genes(flavor="seurat", n_top_genes=2000)`. The
`"seurat"` flavor is specifically adapted for log-transformed inputs, unlike
`"seurat_v3"`, which strictly requires raw integer counts. Selected HVGs were
only flagged, not filtered out, ensuring that the Tsoi marker genes could be
retrieved by name downstream.

**Dimensionality Reduction (PCA).** PCA was computed on the HVG subset with 50
principal components, applying mean-centering only (no `scale`). The first 30
components were used for all downstream graph construction, chosen from the
scree plot (a clear elbow followed by a noise floor; the first 30 PCs
cumulatively explain ≈ 42.5 % of variance).

**Harmony Batch Correction.** Harmony was applied with `patient` as the batch
variable (15 levels), `theta = 2` (default), and `random_state = 0`, producing
the corrected embedding `X_pca_harmony`. `harmonypy.run_harmony` was called
directly rather than through `scanpy.external.pp.harmony_integrate`, because
the scanpy wrapper is incompatible with the PyTorch rewrite in
`harmonypy 0.2.0` (the `Z_corr` orientation changed); calling `run_harmony`
directly is the version-safe path.

**Leiden Clustering and Annotation.** Clustering was performed on the kNN graph
generated from the first 30 dimensions of `X_pca_harmony` using
`flavor="igraph"` across a multi-resolution grid
$[0.1, 0.2, 0.3, 0.5, 0.8, 1.0]$, with UMAP `min_dist = 0.3`. For each
resolution, two diagnostics were deployed simultaneously: a dotplot of the 4
marker genes × clusters, and a cluster × patient crosstabulation (any cluster
with a single-patient proportion $>0.90$ was flagged as a suspected residual
batch artifact). Tsoi state labels were assigned manually by the author based
on biological judgment rather than automated classification. Clusters with
single-patient dominance $>70\%$ were classified as Other (batch?), while
clusters failing to match any known Tsoi state signature were labeled
Ambiguous.

## 3. Results

**HVG Diagnostics.** The 2,000 HVGs were cleanly selected along the upper
boundary of normalized dispersion. $MITF$ and $SOX10$ were not selected;
because the `"seurat"` flavor relies on normalized dispersion, their variance
was not exceptionally high relative to their high baseline expression
levels—an expected technical behavior rather than an error. $AXL$ and $NGFR$
were successfully captured. The top HVGs were heavily dominated by
immune-related genes ($HLA\text{-}DRA$, $CCL4$, $CXCL13$, $CD14$, $HBB$, etc.),
indicating that a fraction of the "malignant" cells harbor residual immune
contamination (noted as a limitation, but left uncleaned to preserve strict
comparability with Balderson).

**Before vs. After Harmony Integration.** Prior to correction, the PCA and
UMAP embeddings fractured into more than 10 isolated patient-specific islands,
demonstrating that the primary axis of variation was patient identity rather
than biological state. Post-correction, the majority of patients integrated
into a single connected manifold. Quantitatively, the patient silhouette score
decreased from 0.33 to 0.08, and the fraction of intra-patient kNN neighbors
dropped from 0.965 to 0.725 (theoretical perfect-mixing baseline = 0.188).
This indicates successful batch correction without over-correction, though
pt79 and pt81/82 remained noticeably unmixed.

**Resolution Comparison.** The number of clusters scaled with resolution as
follows: $\text{res } 0.1/0.2/0.3/0.5/0.8/1.0 = 3/5/6/9/11/15$ clusters.
Over-clustering occurred at $\text{res} \ge 0.5$, where newly emerging clusters
predominantly represented single-patient artifacts (at $\text{res}=1.0$, 5
clusters triggered the $>0.90$ single-patient dominance flag).
$\text{res}=0.3$ (6 clusters) provided the optimal balance, capturing 3 robust
cross-patient biological clusters and 3 batch-driven pseudo-clusters.

**Cluster × Patient Crosstabulation Insights ($\text{res}=0.3$).** Clusters
cl0 (pt88 25 %), cl1 (pt80 36 %), and cl4 (pt79 21 %) were supported by
multiple patients. Conversely, cl2 (pt81 73 %), cl3 (pt80 78 %), and cl5
(pt79 100 %) were dominated by single patients. This underscores that
clustering validity must be assessed by synthesizing dotplots and crosstabs
together: even if cl5 exhibits highly uniform internal marker expression, it
cannot be ascribed a biological interpretation due to complete patient
confounding.

**Final Tsoi Annotation.**

- **cl0, cl1 → Melanocytic ($n=520$):** Characterized by $MITF^{high}$ (mean
  expression 4.32 / 4.08), $SOX10^{\sim}$, $NGFR^{low}$ (mean expression
  0.11 / 0.26, expressing-cell fraction ~7 %), and $AXL^{low}$. These clusters
  are robustly shared across patients. cl1 was initially considered for a
  Transitory designation; however, because the Transitory state is strictly
  defined by $NGFR^{high}$ rather than $AXL^{high}$, and the $NGFR$ signal is
  nearly absent here, it was ultimately classified as a Melanocytic variant.
- **cl4 → Undifferentiated ($n=33$):** Characterized by $MITF^{low}$ (2.35)
  and $AXL^{high\uparrow}$ (1.62, the highest across the entire dataset), with
  cross-patient representation. Notably, $SOX10$ expression here is moderate
  rather than strictly low, and the small sample size ($n=33$) is recorded as
  a caveat.
- **cl2, cl3, cl5 → Other (batch?) ($n=704$):** All exclusively dominated by
  single patients. Among them, cl3 (78 % pt80 cells) exhibited the highest
  $NGFR$ expression (0.83) in the dataset—the precise signal required to
  distinguish the Transitory and Neural crest-like states. However, under
  Harmony's default of $\theta=2$, this signal remains completely inseparable
  from the pt80 patient-specific effect.

Side-by-side visualization of UMAP by `tsoi_state` and UMAP by patient reveals
that the Other (batch?) domain perfectly overlaps with the isolated spatial
pockets of pt79 (as well as pt81/82). This serves as direct visual evidence of
batch-induced pseudo-clustering.

## 4. Discussion

**Conclusion: only two of the four Tsoi states were robustly recapitulated**
(Melanocytic and Undifferentiated). Even after incorporating Harmony batch
correction and upgrading to Leiden clustering, the Tirosh malignant cells
failed to exhibit the clean four-state structure reported by Balderson 2024.
This honest, partial replication stands as the core finding of this research
phase.

**Why intermediate states remain unresolvable.** The absolute expression of
$NGFR$ is exceedingly weak across all clusters (mean $< 0.85$, with the
percentage of expressing cells typically $< 25\%$). Consequently, the
Transitory and Neural crest-like states—which fundamentally depend on
$NGFR^{high}$—lack the requisite marker-gene validation. The only cluster
displaying a relative elevation in $NGFR$ (cl3) is heavily confounded by the
pt80 patient effect and cannot be isolated. For marker-based annotation on the
Tirosh dataset, these intermediate states are inherently unidentifiable,
imposing a hard methodological ceiling on this dataset.

**Why certain clusters are patient-dominated.** Harmony running at its default
parameter ($\theta=2$) fails to sufficiently integrate patients with high cell
counts (pt79/80/81). More fundamentally, malignant cells harbor
patient-specific copy-number variations (CNVs) and distinct clonal genetic
backgrounds; thus the "batch effect" observed here largely represents true
patient biology rather than pure technical variance. Linear integration
frameworks such as Harmony struggle to decouple the two (and aggressively
increasing $\theta$ risks erasing authentic biological states). This specific
challenge provides the motivation for Stage 3, where we will transition to
scVI/scANVI: explicitly modeling batch as a conditional variable within a
non-linear latent space offers a promising avenue to recover states such as
cl3 (pt80, highest $NGFR$) that are currently submerged by patient-specific
signatures.

## 5. Limitations & Future Work

The limitations of the current analysis include: (1) the $NGFR$ signal is weak
across the entire dataset, imposing a hard upper limit on marker-based
annotation of the intermediate states; (2) no parameter sweep was performed
for Harmony's $\theta$ (kept at default to maintain direct comparability with
Balderson); (3) residual immune contamination within the "malignant" subset
was left uncleaned, again for comparability; (4) the analysis relies on a
single dataset (Tirosh), with cross-dataset robustness validation using the
Jerby-Arnon dataset (GSE115978) remaining a stretch goal. Stage 3 plans to
re-implement integration and state identification using scVI/scANVI; however,
a technical constraint must be noted: this dataset consists of
$\log_2(\text{TPM}/10+1)$ values rather than raw integer counts, and because
scVI's negative-binomial likelihood requires raw counts, the prerequisite for
Stage 3 is to obtain the original count matrix from Dr. Tirosh (email request
sent) or, alternatively, to implement a Plan B using count approximation or a
method that accepts log-transformed input.
