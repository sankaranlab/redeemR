# redeemR 2.0.0

## redeemR2.0 framework

- Defines redeemR2.0 as the downstream framework for processing REDEEM-V mitochondrial consensus variant calls before lineage reconstruction.
- Introduces filter2 as the default post-consensus filtering workflow, replacing the original filter1 strategy for most downstream analyses.
- Adds edge trimming of raw genotype calls to reduce residual fragment-end artifacts.
- Adds per-variant binomial goodness-of-fit filtering to remove residual technical noise after consensus calling.
- Adds annotation and filtering utilities for RSRS50, blacklist regions, population features, homopolymer context, hypermutable sites, amino-acid consequences, mitochondrial disease associations, and transition/transversion class.
- Adds matched depth-matrix construction from `QualifiedTotalCts` for downstream lineage and model workflows.
- Adds median-depth and UMI-support metrics for quality control.
- Adds filter2 QC summaries, including mutation spectra, transition/transversion ratios, depth distributions, variant-support plots, and matrix-dimension logs.

## Documentation

- Reworked the README around the five core redeemR2.0 tasks: parsing REDEEM-V calls, post-consensus filtering/QC, edge trimming, binomial noise filtering, and variant annotation/flagging.
- Added a filter2 quickstart with command-line and R examples.
- Added a canonical `vignettes/redeemr-workflow.Rmd` article for the REDEEM-V-to-filter2-object workflow.
- Added `vignettes/filter2-qc.Rmd` and `vignettes/variant-annotation.Rmd` articles.
- Added `_pkgdown.yml` for a pkgdown documentation site with workflow, QC, annotation, reference, and changelog navigation.
