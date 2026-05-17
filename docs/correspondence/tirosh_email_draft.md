# Correspondence — Request for GSE72056 raw counts

**Task:** Q0.3
**Purpose:** Contact the original data author to request raw read counts.
GEO only provides the processed `log2(TPM/10+1)` matrix, while several
benchmarked methods (scVI/scANVI) are designed for raw counts.
**Status:** Sent 2026-05-17 (see Sending Log below).

---

## Recipient

| Field | Value |
|---|---|
| Recipient | Prof. Itay Tirosh, Weizmann Institute of Science |
| Email | Verified against the lab website prior to sending |
| CC | None |

## Subject

```
Request for raw count matrix of GSE72056 — undergraduate research project
```

## Body

```
Dear Dr. Tirosh,

I am an undergraduate student working on a methodological re-analysis of
your 2016 Science paper "Dissecting the multicellular ecosystem of
metastatic melanoma by single-cell RNA-seq," as part of my preparation
for transfer to a four-year university.

I have downloaded the processed expression matrix from GEO (GSE72056),
which is provided as log2(TPM/10+1). My project benchmarks modern deep
learning methods (scVI, scANVI, and potentially single-cell foundation
models) against the classical approaches used in the original paper and
in Balderson et al. (2024). Several of these methods are designed for
raw read counts rather than log-normalized values.

Would it be possible to access the raw read count matrix (e.g.,
unnormalized RSEM/featureCounts estimates) for this dataset? I fully
understand it may not be readily available; in that case, any guidance
on the most appropriate way to approximate counts from the provided TPM
values would be very helpful.

For transparency, the analysis code and documentation for this project
are openly available at:
https://github.com/EwanWu06/melanoma-scrnaseq-reanalysis

Thank you for your time, and for making this foundational dataset
publicly available.

Best regards,
Ewan Wu
UCLA Extension
https://github.com/EwanWu06/melanoma-scrnaseq-reanalysis
```

---

## Sending Log

| Field | Value |
|---|---|
| Date sent | **2026-05-17** |
| Method | Email (personal account) |
| Recipient address | _(verified on lab website prior to sending)_ |
| CC | None |
| Reply received | Pending |
| Reply summary | Pending |
| Next action | If no reply by **2026-05-31**, send one polite follow-up |

> No reply is itself a valid outcome. If no response is received after a
> follow-up, this attempt is documented as part of the methodology, and
> the analysis proceeds with a clearly stated count-approximation strategy.
> Documenting the outreach is itself a marker of methodological rigor.
