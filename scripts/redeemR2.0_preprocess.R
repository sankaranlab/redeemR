##
#' @param input Input folder
#' @param thr Threshold
#' @param edge_trim Edge trim
#' @param min_variant_depth Minimum variant depth
#' @param do_median_depth_filter Whether to filter variants by median depth
#' @param do_filter2_qc Whether to run filter2 QC
#' @return A list with redeemR object and QC report
#' @export
#' do_median_depth_filter is recommended to be TRUE, but if you want to skip it, set it to FALSE,
#' do_median_depth_filter also remove homoplasmy variants that are not detected before the depth correction
preprocessed_redeemr <- function(name,
                                 input,
                                 thr = c("T","LS","S","VS"),
                                 edge_trim        = 9,
                                 min_variant_depth = 5,
                                 plots_dir = NULL,
                                 do_median_depth_filter = T,
                                 do_qc = T){
    # 2) standardize thr
    thr <- match.arg(thr)
    dim_log <- list()
  
    # 3) load & trim variants
    message("[", name, "] reading & trimming genotypes (thr=", thr, ", edge_trim=", edge_trim, ")")
    VariantsGTSummary <- redeemR.read.trim(input, 
                                           thr = thr, 
                                           edge_trim = edge_trim)
    message("[redeemR.read.trim-->]", "VariantsGTSummary is read in. ")
    
    # 4) create redeemR model
    redeemR <- Create_redeemR_model(VariantsGTSummary,
                                    qualifiedCellCut=10,
                                    VAFcut=1,
                                    Cellcut=2) 
    message("[Create_redeemR_model-->]", "redeemR is created by Create_redeemR_model. Cts.Mtx and Cts.Mtx.bi are created")
    .tmp_dim <- append_dim_row(redeemR, "after_create_model"); if (is.data.frame(.tmp_dim)) dim_log[[length(dim_log)+1]] <- .tmp_dim
    
    # 5) filter variants
    redeemR <- clean_redeem(redeemR,
                            fdr=0.05, 
                            min_cell_per_variant=2)
    message("[clean_redeem-->]", "clean_redeem by binomial fdr 0.05. Homoplasmy mutations are not in V.fitered nor in GTsummary.filtered. \
            Cts.Mtx and Cts.Mtx.bi are updated: Now, the matries dimensions are:")
    print_redeemR_matrix_dims(redeemR)
    .tmp_dim <- append_dim_row(redeemR, "after_clean_redeem"); if (is.data.frame(.tmp_dim)) dim_log[[length(dim_log)+1]] <- .tmp_dim

    # 6) add annotation
    redeemR <- add_annotation_redeem(redeemR)
    message("[add_annotation_redeem-->]", "redeemR is updated with add_annotation_redeem for variant population genetics and functional impact-> @V.fitered")
   
    # 7) remove blacklist and RSR50
    redeemR <- clean_redeem_remove_blacklist_RSRS50(redeemR)
    message("[clean_redeem_remove_blacklist_RSRS50-->]", "redeemR is updated with removing blacklist and RSR50, see @V.fitered. \
             @V.fitered, @GTsummary.filtered, Cts.Mtx and Cts.Mtx.bi are all updated: Now, the matries dimensions are:")
    print_redeemR_matrix_dims(redeemR)
    .tmp_dim <- append_dim_row(redeemR, "after_blacklist_RSRS50"); if (is.data.frame(.tmp_dim)) dim_log[[length(dim_log)+1]] <- .tmp_dim

    # 8) Add depth matrix (required for mitodrift)
    redeemR <- Add_DepthMatrix_filter2(redeemR)
    message("[Add_DepthMatrix_filter2-->]", "@Ctx.Mtx.depth added. Required for Mitodrift")
    .tmp_dim <- append_dim_row(redeemR, "after_add_depth_matrix"); if (is.data.frame(.tmp_dim)) dim_log[[length(dim_log)+1]] <- .tmp_dim
    
    # 9) Add median depth information
    redeemR <- add_median_depth_to_redeemR(redeemR)
    message("[add_median_depth_to_redeemR-->]", "Median depth information added to @V.fitered")
    
    # 10) Filter variants by median depth (optional)
    if (do_median_depth_filter) {
        redeemR <- clean_redeem_remove_low_median_depth(redeemR, min_median_depth = min_variant_depth)
        message("[clean_redeem_remove_low_median_depth-->]", "Variants with median depth < ", min_variant_depth, " are filtered out. Now, the matrices dimensions are:")
        print_redeemR_matrix_dims(redeemR)
        .tmp_dim <- append_dim_row(redeemR, "after_median_depth_filter"); if (is.data.frame(.tmp_dim)) dim_log[[length(dim_log)+1]] <- .tmp_dim
    } else {
        message("[median_depth_filter-->]", "Skipping median depth filtering as requested")
    }
    
    # 11) Add prop_2_3 to redeemR
    redeemR <- add_prop_2_3_to_redeemR(redeemR)
    message("[add_prop_2_3_to_redeemR-->]", "prop_2_3 added to @V.fitered")
    .tmp_dim <- append_dim_row(redeemR, "after_add_prop_2_3"); if (is.data.frame(.tmp_dim)) dim_log[[length(dim_log)+1]] <- .tmp_dim
    
    # 9) (Optional) Add filter2 qc plots and metrics
    if (do_qc){
        redeemR <- add_raw_fragment(redeemR, raw = "RawGenotypes.Sensitive.StrandBalance")
        print("raw fragments added for QC")
        message("[add_raw_fragment and run_redeem_qc -->]", "filter2 QC is performed, return a list with redeemR and QC results")
        report <- run_redeem_qc(redeemR, redeemR@HomoVariants)
        # ensure plots_dir exists and save outputs
        if (!dir.exists(plots_dir)) {
            dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
        }
        # Generate basic redeemR QC plots
        # Mutation profile bulk plot (returns ggplot)
        mutation_profile_plot <- MutationProfile.bulk(redeemR@UniqueV)
        
        # Depth plot (now returns list with combined and individual plots)
        depth_plots <- plot_depth(redeemR)
        
        # Variant plot (now returns list with combined and individual plots)
        variant_plots <- plot_variant(redeemR)
        
        # Save all plots in one multi-page PDF file
        pdf_file <- file.path(plots_dir, paste0(name, "_basic_qc.pdf"))
            pdf(pdf_file, width = 14, height = 5)
            
            # Page 1: Mutation profile plot
            print(mutation_profile_plot)
            
            # Page 2: Depth summary plots (combined)
            print(depth_plots$combined)
            
            # Page 3: Variant metrics plots (combined)
            gridExtra::grid.arrange(variant_plots$p1,
                                    variant_plots$p2, 
                                    variant_plots$p3, 
                                    variant_plots$p4, 
                                    ncol = 4)
            
            # Page 4: Filter2 QC plots (from report$plots)
            filter2_qc <- (report$plots$p_pos | report$plots$pos_1mol) /
                        (report$plots$p_cell_maxcts | report$plots$p_cell_meancts) +
                patchwork::plot_annotation(
                    title = paste0(name, " Filter2 QC"),
                    theme = ggplot2::theme(plot.title = ggplot2::element_text(family = "sans"))
                )
            print(filter2_qc)
        
        dev.off()
        

        dim_df <- if (length(dim_log) > 0) do.call(rbind, dim_log) else data.frame()
        write.csv(dim_df,
                  file = file.path(plots_dir, paste0(name, "_matrix_dims.csv")),
                  row.names = FALSE)
        # save simple metrics
        write.csv(as.data.frame(t(report$transversion_rate)),
                         file = file.path(plots_dir, paste0(name, "_transversion_rate.csv")),
                         row.names = FALSE)
        return(list(redeemR = redeemR, report = report))
    }
    message(sprintf("[%s] Done. Returning redeemR object.", name))
    return(redeemR)
    
}

 

# === argparse setup ===
library(argparse)
parser <- ArgumentParser(description = "Preprocess redeemR data")
parser$add_argument("-n", "--name",  required=TRUE, help="*required* Sample name")
parser$add_argument("-i", "--input", required=TRUE, help="*required* redeemV final folder")
parser$add_argument("-o", "--output", required=TRUE, help="*required* Output RDS file")
parser$add_argument("-t", "--thr",   default="S",   choices=c("T","LS","S","VS"),
                    help="Threshold (T, LS, S, VS)")
parser$add_argument("-e", "--edge-trim", type="integer", default=9,
                    help="Minimum edge distance")
parser$add_argument("-d", "--min-variant-depth", type="integer", default=5,
                    help="Minimum median depth threshold for variant filtering (default: 10)")
parser$add_argument("--do-median-depth-filter", action="store_true",
                    help="Skip median depth filtering step")
parser$add_argument("--do-qc", action="store_true",
                    help="Run filter2 QC and return report")

opts <- parser$parse_args()

# === run ===
# library(redeemR)
library(stringr)
library(glue)
library(patchwork)
devtools::load_all("/lab/solexa_weissman/cweng/Packages/redeemR/")
## derive plots_dir: always use '<dirname(output)>/plots'
out_dir <- dirname(opts$output)
plots_dir <- file.path(out_dir, "plots")

result <- preprocessed_redeemr(
  name              = opts$name,
  input             = opts$input,
  thr               = opts$thr,
  edge_trim         = opts$edge_trim,
  min_variant_depth = opts$min_variant_depth,
  plots_dir         = plots_dir,
  do_median_depth_filter = opts$do_median_depth_filter,
  do_qc             = opts$do_qc
)

# === save ===
saveRDS(result, opts$output)
message("Saved to ", opts$output)
