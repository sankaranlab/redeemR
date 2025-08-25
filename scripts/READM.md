# redeemR2.0_preprocess.R
End-to-end R script to preprocess **redeemV** outputs into a cleaned **redeemR** object, optionally run **filter2 QC**, and produce a compact multi-page PDF with depth/variant/QC plots plus CSV summaries.

---

## Features

- Read & trim variants from a redeemV “final” folder
- Build a `redeemR` model and apply FDR-based variant filtering
- Annotate variants; remove blacklist & RSRS50 sites
- Add depth matrix (required for **MitoDrift**) and median-depth stats
- *(Optional, recommended)* Filter variants by **median depth**
- Add `prop_2_3` metric
- *(Optional)* Run **filter2 QC** and write plots & summary CSVs
- Log matrix dimensions at each step

---

## Requirements

- **R** ≥ 4.3  
- Packages: `argparse`, `devtools`, `stringr`, `glue`, `patchwork`, `gridExtra`, `ggplot2`
- **redeemR** available (installed or loadable via a local path)

> The script currently calls following for development purpose:
> 
> ```r
> devtools::load_all("/lab/solexa_weissman/cweng/Packages/redeemR/")
> ```
> 
> Replace with `library(redeemR)` if the package is installed, or point to your local checkout.

---

## Installation

Clone this repository and ensure the required R packages are installed. If `redeemR` is not installed system-wide, adjust the `devtools::load_all()` path in the script to your local checkout or change to `library(redeemR)`.

---

## Command-line Usage

```bash
Rscript preprocess_redeemr.R --name NAME --input DIR --output FILE [options]
```

### Required arguments
| Flag | Description |
|---|---|
| `-n, --name` | Sample name (used in file prefixes) |
| `-i, --input` | Path to **redeemV final** folder (directory) |
| `-o, --output` | Output `.rds` file |

### Optional arguments
| Flag | Default | Description |
|---|---:|---|
| `-t, --thr` | `S` | Threshold: one of `{T, LS, S, VS}` |
| `-e, --edge-trim` | `9` | Minimum edge distance to trim |
| `-d, --min-variant-depth` | `5` | **Median** depth cutoff for variant filtering |
| `--do-median-depth-filter` | *(off)* | **Enable** median-depth filtering *(recommended)* |
| `--do-qc` | *(off)* | Run filter2 QC and include QC report & plots |

> ⚠️ **Notes on flags**
> - `--do-median-depth-filter` **enables** filtering (the in-script help string says “Skip…”, but passing the flag sets filtering **on**).
> - The code default for `--min-variant-depth` is **5** (help text elsewhere may mention 10).

---

## Quick Start

**Recommended run (median-depth filtering + QC):**
```bash
Rscript preprocess_redeemr.R   --name DN12   --input /path/to/redeemV_final   --output ./outputs/   --thr S   --edge-trim 9   --min-variant-depth 5   --do-median-depth-filter   --do-qc
```

**Minimal run (no QC, no depth filtering):**
```bash
Rscript preprocess_redeemr.R   -n DN12   -i /data/redeemV/DN12_final   -o ./out/
```

---

## Programmatic Use (R)

```r
# devtools::load_all("/path/to/redeemR/")  # or library(redeemR)
source("preprocess_redeemr.R")

res <- preprocessed_redeemr(
  name                  = "DN12",
  input                 = "/data/redeemV/DN12_final",
  thr                   = "S",
  edge_trim             = 9,
  min_variant_depth     = 5,
  plots_dir             = "./out/plots",
  do_median_depth_filter= TRUE,
  do_qc                 = TRUE
)

saveRDS(res, "./out/DN12.redeemR.rds")
```

---

## Outputs

- **RDS** at `--output`
  - With `--do-qc`: a **list** containing `redeemR` and `report` (filter2 QC results)
  - Without `--do-qc`: the **redeemR** object
- **Plots** in `dirname(--output)/plots/`
  - `<name>_basic_qc.pdf` (multi-page):
    1. Mutation profile
    2. Depth summary (combined)
    3. Variant metrics (combined)
    4. Filter2 QC panel
- **CSV summaries** in the same `plots` directory:
  - `<name>_matrix_dims.csv` — matrix dimensions across pipeline steps
  - `<name>_transversion_rate.csv` — basic QC metric

---

## Processing Steps (Detail)

1. `redeemR.read.trim(input, thr, edge_trim)`  
2. `Create_redeemR_model(qualifiedCellCut=10, VAFcut=1, Cellcut=2)`  
3. `clean_redeem(fdr=0.05, min_cell_per_variant=2)`  
4. `add_annotation_redeem(...)`  
5. `clean_redeem_remove_blacklist_RSRS50(...)`  
6. `Add_DepthMatrix_filter2(...)` → adds `@Ctx.Mtx.depth`  
7. `add_median_depth_to_redeemR(...)`  
8. *(Optional)* `clean_redeem_remove_low_median_depth(min_median_depth = <—min-variant-depth>)`  
9. `add_prop_2_3_to_redeemR(...)`  
10. *(Optional)* QC:  
    - `add_raw_fragment(redeemR, raw = "RawGenotypes.Sensitive.StrandBalance")`  
    - `run_redeem_qc(redeemR, redeemR@HomoVariants)`  
    - Plot helpers: `MutationProfile.bulk`, `plot_depth`, `plot_variant`  

Matrix dimensions are appended via `append_dim_row()` and saved as CSV.

---

## Known Quirks / TODO

- Make the `redeemR` load path configurable or switch to `library(redeemR)`.


---
