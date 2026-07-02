![GitHub R package version](https://img.shields.io/github/r-package/v/chenweng1991/redeemR?label=ReDeeM)
![GitHub last commit](https://img.shields.io/github/last-commit/chenweng1991/redeemR)
![GitHub](https://img.shields.io/github/license/chenweng1991/redeemR)

# redeemR

`redeemR` is the downstream R framework for processing REDEEM-V mitochondrial
consensus variant calls into filtered, annotated, quality-controlled objects for
single-cell lineage analysis. The package supports post-consensus variant
filtering, annotation, quality control, cell-by-variant matrix construction, and
utilities for downstream lineage reconstruction and multiome integration.

ReDeeM stands for single-cell **Re**gulatory multi-omics with **Dee**p
**M**itochondrial mutation profiling. The broader ReDeeM analysis stack has two
parts:

- [`redeemV`](https://github.com/chenweng1991/redeemV): mapping and deep
  mitochondrial consensus mutation calling from raw sequencing data.
- `redeemR`: downstream processing of REDEEM-V output before lineage
  reconstruction and single-cell integrative analysis.

![ReDeeM overview](https://github.com/chenweng1991/redeemR/assets/43254272/da3c9a70-53c8-4861-b3ac-3f351a1b540f)

## What is redeemR2.0?

`redeemR2.0` is the updated downstream analysis framework used to process
REDEEM-V mitochondrial consensus variant calls before lineage reconstruction. It
converts consensus-filtered mtDNA calls into lineage-ready `redeemR` objects
through five major steps:

1. Parse consensus-filtered mutation calls from REDEEM-V.
2. Apply post-consensus variant filtering and quality control.
3. Remove residual fragment-end artifacts with edge trimming.
4. Remove residual technical artifacts with binomial noise filtering.
5. Flag inherited or pre-existing variants and annotate functional impact.

In this branch, `filter2` refers specifically to the default redeemR2.0
post-consensus filtering workflow. The term distinguishes this updated strategy
from the original redeemR filtering strategy, referred to as `filter1`. Filter2
is applied after UMI-based consensus variant calling and is designed to retain
high-fidelity somatic mtDNA variants while removing residual technical artifacts,
inherited or pre-existing polymorphisms, and variants less suitable for neutral
lineage tracing.

## Installation

Install the current development branch:

```r
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("chenweng1991/redeemR", ref = "redeemR2.0")
library(redeemR)
```

After `redeemR2.0` is merged into the default branch, installation will use:

```r
devtools::install_github("chenweng1991/redeemR")
```

## Quick start

Filter2 starts from a REDEEM-V `final/` output directory. Unless otherwise
specified, analyses use the stringent consensus output (`thr = "S"`).

A REDEEM-V `final/` directory should contain:

```text
QualifiedTotalCts
RawGenotypes.Sensitive.StrandBalance
RawGenotypes.Specific.StrandBalance
RawGenotypes.Total.StrandBalance
RawGenotypes.VerySensitive.StrandBalance
```

Run filter2 preprocessing from the command line:

```bash
Rscript scripts/redeemR2.0_preprocess.R \
  --name sample1 \
  --input /path/to/redeemV/final \
  --output preprocessing/filter2/sample1/sample1.S.redeemR_filter2_adddepth.rds \
  --thr S \
  --edge-trim 9 \
  --min-variant-depth 5 \
  --do-median-depth-filter \
  --do-qc
```

The command writes a cleaned `redeemR` object. With `--do-qc`, it also writes QC
plots, matrix-dimension logs, and summary metrics under
`preprocessing/filter2/sample1/plots/`.

## R workflow

The same workflow can be run directly in R:

```r
library(redeemR)

redeemv_final <- "/path/to/redeemV/final"

variants <- redeemR.read.trim(
  path = redeemv_final,
  thr = "S",
  edge_trim = 9
)

obj <- Create_redeemR_model(
  variants,
  qualifiedCellCut = 10,
  VAFcut = 1,
  Cellcut = 2
)

obj <- clean_redeem(obj, fdr = 0.05, min_cell_per_variant = 2)
obj <- add_annotation_redeem(obj)
obj <- clean_redeem_remove_blacklist_RSRS50(obj)
obj <- Add_DepthMatrix_filter2(obj)
obj <- add_median_depth_to_redeemR(obj)
obj <- clean_redeem_remove_low_median_depth(obj, min_median_depth = 5)
obj <- add_prop_2_3_to_redeemR(obj)

saveRDS(obj, "sample1.S.redeemR_filter2_adddepth.rds")
```

To generate filter2 QC summaries from an existing object:

```r
obj <- add_raw_fragment(obj, raw = "RawGenotypes.Sensitive.StrandBalance")
qc <- run_redeem_qc(obj, obj@HomoVariants)
```

## Filter2 workflow

Filter2 applies the following sequence:

1. Import REDEEM-V consensus variant calls from `final/`.
2. Annotate raw genotype calls by distance to the nearest fragment end.
3. Remove calls within the edge-trimmed window. The default threshold is 9 bp.
4. Reconstruct genotype summaries after trimming, including per-cell variant UMI
   counts, site depth, and heteroplasmy for each retained cell-variant record.
5. Create an initial `redeemR` object from trimmed genotype summaries.
6. Retain cells with mean mitochondrial coverage >= 10.
7. Initially retain candidate variants detected in at least two cells.
8. Annotate homoplasmic or near-homoplasmic variants using broad detection, high
   mean heteroplasmy, and low variability.
9. Build heteroplasmic cell-by-variant count and binary matrices.
10. Apply a per-variant binomial goodness-of-fit test and retain variants passing
    FDR <= 0.05.
11. Annotate retained variants for population, sequence-context, and functional
    features.
12. Remove heteroplasmic RSRS50 variants and variants in predefined mitochondrial
    blacklist regions for downstream lineage analysis.
13. Generate a matched depth matrix from `QualifiedTotalCts`.
14. Add median-depth and UMI-support metrics.
15. Optionally remove variants with median depth < 5 and depth-corrected possible
    homoplasmic variants.

## Default parameters

| Step | Default |
|---|---|
| Consensus input | Stringent consensus output, `thr = "S"` |
| Edge trimming | 9 bp from fragment ends |
| Cell coverage filter | Mean mtDNA coverage >= 10 |
| Initial variant cell filter | Detected in >= 2 cells |
| Supporting UMI threshold | At least 1 supporting UMI in an individual cell |
| Homoplasmic or near-homoplasmic annotation | `CellNPCT > 0.75`, `PositiveMean > 0.75`, `CV < 0.01` |
| Binomial filter | Per-variant goodness-of-fit test |
| FDR threshold | `q <= 0.05` |
| Optional median-depth filter | Remove variants with median depth < 5 |
| Downstream lineage exclusions | RSRS50 variants and mitochondrial blacklist regions |

## Variant annotations

`redeemR2.0` annotates retained variants with population, context, and functional
features, including:

- MITOMAP population frequencies
- RSRS50 ancestral-state status
- haplogroup marker counts
- mitochondrial blacklist regions
- homopolymer context
- hypermutable-site labels
- amino-acid consequence and predicted coding impact
- mitochondrial disease annotations
- transition/transversion class

## Outputs

The final `redeemR2.0` object contains:

| Slot | Description |
|---|---|
| `@V.fitered` | Filtered variant-level summary and annotations |
| `@GTsummary.filtered` | Filtered cell-variant genotype records |
| `@Cts.Mtx` | Cell-by-variant mutant allele count matrix |
| `@Cts.Mtx.bi` | Binarized cell-by-variant mutation matrix |
| `@Ctx.Mtx.depth` | Matched cell-by-variant depth matrix |
| `@HomoVariants` | Homoplasmic or near-homoplasmic variants identified during preprocessing |

QC outputs include mutation-spectrum summaries, transition/transversion ratios,
mtDNA depth distributions, variant-support plots, matrix-dimension logs across
filtering steps, and filter2 diagnostic plots.

## Main functions

| Task | Functions |
|---|---|
| Parse REDEEM-V output | `redeemR.read()`, `redeemR.read.trim()` |
| Create redeemR object | `Create_redeemR()`, `Create_redeemR_model()` |
| Filter variants | `clean_redeem()`, `clean_redeem_remove_blacklist_RSRS50()`, `clean_redeem_remove_low_median_depth()` |
| Add depth and QC metrics | `Add_DepthMatrix_filter2()`, `add_median_depth_to_redeemR()`, `add_prop_2_3_to_redeemR()` |
| Annotate variants | `add_annotation_redeem()`, `annotate_all_variants()` |
| Convert formats | `convert_redeem_matrix_long()` |
| Validate objects | `check_redeem()`, `print_redeemR_matrix_dims()` |

## Documentation

- [redeemR2.0 workflow](./vignettes/redeemr-workflow.Rmd)
- [Getting started](./vignettes/01_Get_Started.ipynb)
- [Use filter2](./vignettes/02_redeem_filter2.ipynb)
- [Filter2 QC report](./vignettes/filter2-qc.Rmd)
- [Variant annotation](./vignettes/variant-annotation.Rmd)
- [Legacy variant annotation notebook](./vignettes/03_variant_annotation.ipynb)
- [redeemR2.0 preprocess script](./scripts/)
- [ReDeeM paper analysis reproducibility](https://github.com/chenweng1991/redeem_reproducibility)
- [Extended ReDeeM robustness analysis](https://github.com/chenweng1991/redeem_robustness_reproducibility)
- [ReDeeM filtering principles and strategies](https://github.com/chenweng1991/redeemR/wiki/ReDeeM-filtering-strategies)

## Citation

Please cite the ReDeeM study: [Deciphering cell states and genealogies of human
hematopoiesis](https://doi.org/10.1038/s41586-024-07066-z).

## Contact

Questions and feedback are welcome. Contact Chen Weng at
cweng@broadinstitute.org.
