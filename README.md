![GitHub R package version](https://img.shields.io/github/r-package/v/chenweng1991/redeemR?label=ReDeeM)
![GitHub last commit](https://img.shields.io/github/last-commit/chenweng1991/redeemR)
![GitHub](https://img.shields.io/github/license/chenweng1991/redeemR)

# redeemR

<!-- badges: start -->
<!-- badges: end -->


## Introduction
ReDeeM: single-cell **Re**gulatory multi-omics with **Dee**p **M**itochondrial mutation profiling. ReDeeM is a single-cell multiomics platform featuring ultra-sensitive mitochondrial DNA (mtDNA) variant calling and joint RNA+ATAC profiling. ReDeeM enables fine-scale lineage tracing at single cell level, allowing for subclonal and phylogenetic analyses, with simultaneous integrative analyses of cell-state and gene regulatory circuits.</br> 


The analytical pipelines for ReDeeM analysis includes two parts:
- [redeemV](https://github.com/chenweng1991/redeemV) is set of in-house Bash pipeline and python scripts for mapping and deep mitochondrial mutation calling. (Input from rawdata)
- [redeemR](https://github.com/chenweng1991/redeemR) is an in-house R package for downstream lineage tracing and single cell integrative analysis. (**This page**, Input from redeemV)

**redeemR** is designed to refine the somatic mitochondiral mutations, build the single-cell phylogenetic tree and facilitate genealogical and clonal/subclonal-resolved single-cell multiomics analysis.
![Github figs](https://github.com/chenweng1991/redeemR/assets/43254272/da3c9a70-53c8-4861-b3ac-3f351a1b540f)


## Installation

You can install the development version of redeemR:

``` r
devtools::install_github("chenweng1991/redeemR")
library(redeemR)
```
## Script
- [redeemR2.0 preprocess](./scripts/)
## Update
- 2024-8-24 We provide an additional filtering option in the ReDeeM-R package. This approach effectively eliminates the position biases, leads to a mutational signature
indistinguishable from bona fide mitochondrial mutations, and removes excess low molecule high connectedness mutations with high sensitivity. Please see Document Use filter-2 below

## New Functions in redeemR 2.0

### Filtering / Object Manipulation
| Function | Description |
|---|---|
| `Add_DepthMatrix_filter2` | Builds the depth matrix slot (`@Ctx.Mtx.depth`) |
| `add_annotation_redeem` | Adds variant annotations to a redeemR object |
| `clean_redeem_remove_blacklist_RSRS50` | Removes RSRS50 blacklisted variants |
| `clean_redeem_remove_low_median_depth` | Filters variants by minimum median depth |
| `check_redeem` | Validates a redeemR object |
| `make_V.fitered_from_mtx` | Reconstructs `@V.fitered` from a matrix |
| `subset_redeem` | Subsets a redeemR object by a cell whitelist |
| `filter_redeemR_by_cells` | Filters by cell whitelist with a minimum cells-per-variant threshold |
| `print_redeemR_matrix_dims` / `safe_dim` / `append_dim_row` | Matrix dimension reporting utilities |

### Variant Annotation
| Function | Description |
|---|---|
| `Annotate_base_change` | Annotates base change type for each variant |
| `annotate_variants_aachange` | Adds amino acid change annotations |
| `annotate_variants_blacklist` | Flags blacklisted variant positions |
| `annotate_variants_homopolymer` | Flags homopolymer-associated variants |
| `annotate_variants_hypermutable` | Flags hypermutable sites |
| `annotate_variants_mito_disease` | Annotates known mitochondrial disease associations |
| `annotate_variants_population_stats` | Adds population-level frequency statistics |
| `annotate_all_variants` | Runs the full annotation suite in one call |

### Depth & Statistics
| Function | Description |
|---|---|
| `add_median_depth_to_redeemR` | Adds per-variant median depth to the object |
| `add_prop_2_3_to_redeemR` | Adds UMI proportion statistics |

### Format Conversion
| Function | Description |
|---|---|
| `convert_redeem_matrix_long` | Converts variant and depth matrices to long format |


## Documentation 
### Notebooks
- [Getting started](./vignettes/Get_Started.ipynb)
- [Use filter-2](./vignettes/redeem_filter2.ipynb)
- [ReDeeM paper analysis reproducibility](https://github.com/chenweng1991/redeem_reproducibility)
- [Extended ReDeeM robustness analysis](https://github.com/chenweng1991/redeem_robustness_reproducibility)
### WiKi
- [redeem filtering principles and strategies](https://github.com/chenweng1991/redeemR/wiki/ReDeeM-filtering-strategies)

## Citation
Please check out our study of human hematopoiesis using ReDeeM [Deciphering cell states and genealogies of human hematopoiesis](https://doi.org/10.1038/s41586-024-07066-z)

## Contact
If you have any questions or suggestions, please feel free to contact us. Feedbacks are very welcome! (Chen Weng, cweng@broadinstitute.org)
