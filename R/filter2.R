
## Internal function for Create_redeemR_model
goodness_of_fit_test <- function(trials,coverage = 30) {
  # Number of trials
  n_trials <- length(trials)
  # Handle possible coverage exceed
  if (max(trials)>coverage){
    coverage <- max(trials)
  }  
  # Observed frequencies
  observed_frequencies <- table(factor(trials, levels = 0:coverage))
  # Single coin model
  total_heads <- sum(trials)
  total_flips <- n_trials * coverage
  p_hat <- total_heads / total_flips
  # Expected frequencies under the single binomial model
  expected_frequencies <- n_trials * dbinom(0:coverage, size = coverage, prob = p_hat)
  # Handle the 0 in expected frequency
  expected_frequencies[expected_frequencies==0]<- 1e-318  
  # Chi-squared test statistic
  chi2_stat <- sum((observed_frequencies - expected_frequencies)^2 / expected_frequencies)
  # p-value from chi-squared distribution with coverage degrees of freedom (coverage + 1 possible outcomes - 1 estimated parameter)
  p_value <- pchisq(chi2_stat, df = coverage, lower.tail = FALSE)
  return(list(chi2_stat = chi2_stat, p_value = p_value))
}

## Internal function, executed by Create_redeemR_model
run_binomial_noise_removal <- function(redeem){
    require(qvalue)
    if ("chi" %in% colnames(redeem@V.fitered)){
        stop("goodness_of_fit_test has been run on this dataset")
    }else{    
        Mtx <- redeem@Cts.Mtx
        stats <- c()
        pvalues <- c()
        for (i in 1:ncol(Mtx)){
            v<- as.integer(redeem@Cts.Mtx[,i])
            pos <- as.numeric(gsub("Variants([0-9]+)[A-Za-z]{2}$", "\\1", colnames(Mtx)))
            cov<-as.integer(as.numeric(redeem@DepthSummary$Pos.MeanCov[pos[i],"meanCov"]))
            res<-goodness_of_fit_test(v,cov)
            stats<-c(stats,res$chi2_stat)
            pvalues<-c(pvalues,res$p_value)
        }
        qvalues <- qvalue(pvalues)$qvalues
        redeem@V.fitered<-merge(redeem@V.fitered,data.frame(Variants=convert_variant(colnames(Mtx)),pvalues=pvalues, qvalues=qvalues, chi=stats),all = T)
        # Also asign this back to redeem@V.ini
        redeem@V.ini <- redeem@V.fitered
        return(redeem)
        }    
}

#' Create a redeemR object from raw variant summaries
#'
#' This builds a redeemR S4 object by:
#' 1) filtering cells by median coverage (`qualifiedCellCut`),  
#' 2) filtering variants by cell-count and VAF (`Cellcut`, `VAFcut`),  
#' 3) initializing all slots (`GTsummary.ini`, `V.ini`, `GTsummary.filtered`, `V.fitered`, etc.),  
#' 4) running `Make_matrix()` to build count matrices,  
#' 5) and performing binomial noise removal.
#'
#' @param VariantsGTSummary A `data.frame` produced by `redeemR.read.trim()`, with attributes  
#'   - `thr` (threshold),  
#'   - `depth` (output of `DepthSummary()`),  
#'   - (optionally) `edge_trim`, `path`, `combined`, `suffix`.
#' @param qualifiedCellCut Numeric; minimum **median** mitochondrial coverage for a cell to be kept (default 10).
#' @param VAFcut Numeric in (0,1]; only variants with **variant allele frequency** <= `VAFcut` are considered (default 1).
#' @param Cellcut Integer >= 1; only variants seen in at least `Cellcut` cells are kept (default 2).
#'
#' @return An object of S4 class **redeemR**, with slots  
#'   - `GTsummary.ini`, `V.ini`: the unfiltered genotype & variant tables  
#'   - `GTsummary.filtered`, `V.fitered`: the post-filter tables  
#'   - `CellMeta`, `DepthSummary`, `HomoVariants`, `UniqueV`, `para`, `attr`, ...  
#'   and with count matrices in `@Cts.Mtx` / `@Cts.Mtx.bi` populated.
#'
#' @examples
#' \dontrun{
#'   vgtsum <- redeemR.read.trim("path/to/data", thr="S", edge_trim=9)
#'   rObj   <- Create_redeemR_model(vgtsum, qualifiedCellCut=10, VAFcut=0.5, Cellcut=3)
#'   summary(rObj)
#' }
#'
#' @export
#' @import Seurat ape phytools phangorn tidytree ggtreeExtra
#' @importFrom ggtree ggtree

Create_redeemR_model<-function(VariantsGTSummary=VariantsGTSummary,qualifiedCellCut=10,VAFcut=1,Cellcut=2){
 if ("edge_trim" %in% names(attributes(VariantsGTSummary))){
        edge_trim <- as.numeric(attr(VariantsGTSummary,"edge_trim"))
    }else{
        edge_trim <- 0
    }
CellMeta<-subset(attr(VariantsGTSummary,"depth")[["Cell.MeanCov"]],meanCov>=qualifiedCellCut)
names(CellMeta)[1]<-"Cell"
VariantsGTSummary.feature<-Vfilter_v4(VariantsGTSummary,Min_Cells = Cellcut, Max_Count_perCell = 1, QualifyCellCut = qualifiedCellCut)
GTsummary.filtered<-subset(VariantsGTSummary,Variants %in% VariantsGTSummary.feature$Variants & Cell %in% CellMeta$Cell)
ob<-new("redeemR")
ob@CellMeta<-CellMeta
## Populate GTsummary.ini and V.ini
ob@GTsummary.ini<-GTsummary.filtered
ob@V.ini=VariantsGTSummary.feature
## Populate GTsummary.filtered and V.fitered, which are the same as ini here.
ob@GTsummary.filtered<-GTsummary.filtered
ob@V.fitered <-VariantsGTSummary.feature
ob@HomoVariants<-attr(VariantsGTSummary.feature,"HomoVariants")
ob@UniqueV<-VariantsGTSummary.feature$Variants
ob@DepthSummary<-attr(VariantsGTSummary,"depth")
ob@para<-c(Threhold=attr(VariantsGTSummary,"thr"),qualifiedCellCut=qualifiedCellCut,VAFcut=VAFcut,Cellcut=Cellcut,edge_trim=edge_trim)
ob@attr<-list(path=attr(VariantsGTSummary,"path"))

if ("combined" %in% names(attributes(VariantsGTSummary))){
    ob@attr <- c(ob@attr,combine_sample=list(combine=attr(VariantsGTSummary,"combined"),suffix=attr(VariantsGTSummary,"suffix")))
}
ob<-Make_matrix(ob,onlyhetero=T)
ob<-run_binomial_noise_removal(ob)
return(ob)
}


#' clean_redeem
#'
#' This function is to clean redeem by filtering both V.fitered and GTsummary.filtered by qvalues
#' 
#' @param ob redeem object
#' @param fdr fdr cutoff, default is 0.05
#' @export
clean_redeem <-function(ob,fdr=0.05, min_cell_per_variant=2){
    require(dplyr)
    ob@V.fitered <- ob@V.ini %>% 
                    filter(HomoTag == "Hetero") %>% 
                    filter(qvalues<=fdr) %>%
                    filter(CellN >= min_cell_per_variant)
    ob@GTsummary.filtered<-subset(ob@GTsummary.ini, Variants %in% ob@V.fitered$Variants)
    message("Homoplasmic variants are removed from V.fitered and GTsummary.filtered")
    ob<-Make_matrix(ob,onlyhetero=T)  
    ob@UniqueV <- ob@V.fitered$Variants
    return(ob)
}

#' add_annotation_redeem
#'
#' This function is to simply add all annotation in 
#' 
#' @param ob redeem object
#' @param fdr fdr cutoff, default is 0.05
#' @export
add_annotation_redeem <-function(ob){
    ob@V.fitered <- annotate_all_variants(ob@V.fitered)
    return(ob)
}


#' clean_redeem_remove_blacklist_RSRS50
#'
#' This function is to clean redeem by filtering both V.fitered and GTsummary.filtered by qvalues
#' 
#' @param ob redeem object
#' @export
clean_redeem_remove_blacklist_RSRS50 <-function(ob){
    variant_to_remove<- ob@V.fitered %>% 
                        filter((HomoTag == "Hetero" & RSRS50 == "Yes") | blacklist_region == "blacklist_region") %>% 
                        pull(Variants)
    ob@V.fitered <- subset(ob@V.fitered,!Variants %in% variant_to_remove)
    ob@GTsummary.filtered<-subset(ob@GTsummary.filtered, !Variants %in% variant_to_remove)
    ob<-Make_matrix(ob,onlyhetero=T)
    ob@UniqueV <- ob@V.fitered$Variants
    return(ob)
}

#' clean_redeem_removehomo
#'
#' This function is to clean redeem by filtering both V.fitered and GTsummary.filtered by qvalues
#' 
#' @param ob redeem object
#' @param hotcall fdr cutoff, default is 0.05
#' @export
clean_redeem_removehomo <-function(ob){
  message("Removing homoplasmy variants from V.fitered and GTsummary.filtered")
    
    ob@V.fitered <- subset(ob@V.fitered,HomoTag != "Homo")
    ob@GTsummary.filtered<-subset(ob@GTsummary.filtered, Variants %in% ob@V.fitered$Variants)
    ob<-Make_matrix(ob,onlyhetero=T)
    ob@UniqueV <- ob@V.fitered$Variants
    return(ob)
}

#' clean_redeem_removehotcall
#'
#' This function is to clean redeem by filtering both V.fitered and GTsummary.filtered by qvalues
#' 
#' @param ob redeem object
#' @param hotcall fdr cutoff, default is 0.05
#' @export
clean_redeem_removehot <-function(ob,hotcall= c("310_T_C","3109_T_C")){
    ob@V.fitered <- subset(ob@V.fitered,!Variants %in% hotcall)
    ob@GTsummary.filtered<-subset(ob@GTsummary.filtered, !Variants %in% hotcall)
    ob<-Make_matrix(ob,onlyhetero=T)
    ob@UniqueV <- ob@V.fitered$Variants
    return(ob)
}




#' Add Filter-2-adjusted depth matrix to a redeemR object
#'
#' This function reads the “QualifiedTotalCts” matrix from disk (if not already supplied),
#' applies the Filter-2 adjustment—subtracting counts for UMIs removed by edge trimming
#' (zeros and unaffected entries remain unchanged)—then reshapes the result into a matrix
#' matching the dimensions of `object@Cts.Mtx.bi` and stores it in
#' `object@Ctx.Mtx.depth`.
#'
#' @param object A \code{redeemR} object
#' @param QualifiedTotalCts Optional data.frame of the same format as the
#'   “QualifiedTotalCts” file. If \code{NULL}, it will be read from
#'   \code{file.path(object@attr$path, "QualifiedTotalCts")}.
#' @return The input \code{redeemR} object, with
#'   \code{object@Ctx.Mtx.depth} populated.
#' @export
Add_DepthMatrix_filter2 <- function(object, QualifiedTotalCts = NULL) {
    # 1. load the table if needed
    if (is.null(QualifiedTotalCts)) {
        QualifiedTotalCts<-fread(paste(object@attr$path,"/QualifiedTotalCts",sep=""))
    }
    colnames(QualifiedTotalCts)<-c("Cell","Pos","T","LS","S","VS")
    message("[Step 1] QualifiedTotalCts loaded with columns: ", paste(colnames(QualifiedTotalCts), collapse = ", "))
    
    message("[Step 2] Build DepthMatrix")
    GTsummary.filtered = object@GTsummary.filtered %>% 
        select(Variants, Freq, depth, Cell)

    GTsummary.filtered.complete.adj = GTsummary.filtered %>%
        tidyr::separate(Variants, c('Pos', 'Ref', 'Alt'), sep = '_', remove = FALSE) %>%
        mutate(Pos = as.integer(Pos)) %>%
        tidyr::complete(
            Cell, 
            tidyr::nesting(Variants, Pos, Ref, Alt),
            fill = list('Freq' = 0)
        ) %>%
        left_join(QualifiedTotalCts, by = join_by(Cell, Pos)) %>%
        mutate(depth = ifelse(is.na(depth), S, depth))

    # convert to matrix format
    DepthMatrix<-reshape2::dcast(GTsummary.filtered.complete.adj,Cell~Variants, value.var = 'depth') %>% 
        tibble::column_to_rownames("Cell") %>% as.matrix
    
    # change back the names
    colnames(DepthMatrix) = colnames(DepthMatrix) %>% stringr::str_remove_all('_') %>% paste0('Variants', .)

    if (all(colnames(object@Cts.Mtx) %in% colnames(DepthMatrix)) & all(rownames(object@Cts.Mtx) %in% rownames(DepthMatrix))){
        message("[Step 3] Assigning DepthMatrix to object@Ctx.Mtx.depth")
        object@Ctx.Mtx.depth<-DepthMatrix[rownames(object@Cts.Mtx),colnames(object@Cts.Mtx)]
    }else{
        print(dim(object@Cts.Mtx))
        print(dim(DepthMatrix))
        print("Missing variants or cells in DepthMatrix")
    }
    
    return(object)
}

#' Remove variants with low median depth and depth-corrected homoplasmy
#'
#' @description
#' Removes variants whose per‑variant median depth (from \code{@V.fitered$median_depth})
#' is below a threshold, and also removes depth‑corrected possible homoplasmy
#' (\code{CellNPCT > 0.75}). After filtering, rebuilds matrices and the depth matrix.
#'
#' @param ob A \code{redeemR} object that already has \code{median_depth} and \code{CellNPCT}
#'   in \code{@V.fitered} (see \code{add_median_depth_to_redeemR}).
#' @param min_median_depth Numeric. Minimum per‑variant median depth to keep (default: 5).
#'
#' @return The updated \code{redeemR} object with variants filtered and matrices rebuilt.
#'
#' @examples
#' # ob <- add_median_depth_to_redeemR(ob)
#' # ob <- clean_redeem_remove_low_median_depth(ob, min_median_depth = 10)
#'
#' @export
clean_redeem_remove_low_median_depth <- function(ob, min_median_depth = 5) {
    if (is.null(ob@V.fitered$median_depth)) {
        stop("No 'median_depth' column found in ob@V.fitered. Please run add_median_depth_to_redeemR() first.")
    }
    
    # Identify variants to remove
    variant_to_remove <- ob@V.fitered %>%
        filter(median_depth < min_median_depth) %>%
        pull(Variants)
        
    # Also remove depth-corrected possible homoplasmy after fixing CellPCT considering the depth
    # Remove variants where CellNPCT (CellN/cellN_depth_gt0) > 0.75 (i.e., possible homoplasmy)
    depth_corrected_homo_variants <- ob@V.fitered %>%
        filter(CellNPCT > 0.75) %>%
        pull(Variants)
    if (length(depth_corrected_homo_variants) > 0) {
        message("Additionally removing ", length(depth_corrected_homo_variants), " depth-corrected possible homoplasmy variants (CellNPCT > 0.75)")
    }
    remove_variants <- union(variant_to_remove, depth_corrected_homo_variants)
    # Filter V.fitered and GTsummary.filtered
    ob@V.fitered <- subset(ob@V.fitered, !Variants %in% remove_variants)
    ob@GTsummary.filtered <- subset(ob@GTsummary.filtered, !Variants %in% remove_variants)
    
    # Update UniqueV
    ob@UniqueV <- ob@V.fitered$Variants
    
    # Rebuild matrices using Make_matrix
    ob <- Make_matrix(ob, onlyhetero = T)
    
    # Rebuild depth matrix using Add_DepthMatrix_filter2
    ob <- Add_DepthMatrix_filter2(ob)
    
    message("Removed ", length(variant_to_remove), " variants with median depth < ", min_median_depth)
    message("Matrices and depth matrix have been rebuilt")
    
    return(ob)
}

# ============================================================================
# Extra filtering functions
# ============================================================================

#' Update redeemR object from GTsummary.filtered (helper)
#'
#' @description
#' Synchronizes \code{@V.fitered} with \code{@GTsummary.filtered}: keeps only present variants,
#' computes per‑variant summaries (CellN, PositiveMean, PositiveMean_cts, maxcts, CV,
#' TotalVcount), merges them into \code{@V.fitered}, updates \code{@CellMeta} and \code{@UniqueV},
#' and rebuilds matrices and depth matrix.
#'
#' @param redeemR_obj A \code{redeemR} object with a populated \code{@GTsummary.filtered}.
#' @param QualifiedTotalCts Optional “QualifiedTotalCts” table passed to \code{Add_DepthMatrix_filter2}.
#'
#' @return The updated \code{redeemR} object.
#'
#' @examples
#' # ob@GTsummary.filtered <- dplyr::filter(ob@GTsummary.filtered, Freq >= 2)
#' # ob <- update_redeemR_from_GTsummary(ob)
#'
#' @export
update_redeemR_from_GTsummary <- function(redeemR_obj, QualifiedTotalCts = NULL, update_depth_matrix = TRUE) {
  ob <- redeemR_obj

  if (!("GTsummary.filtered" %in% slotNames(ob)) || !is.data.frame(ob@GTsummary.filtered)) {
    stop("redeemR@GTsummary.filtered is required to update V.fitered")
  }

  gts <- ob@GTsummary.filtered
  required_cols <- c("Variants", "hetero", "Freq")
  if (!all(required_cols %in% names(gts))) {
    stop("GTsummary.filtered must contain columns: ", paste(required_cols, collapse = ", "))
  }

  var_feat <- gts %>%
    dplyr::group_by(Variants) %>%
    dplyr::summarise(
      CellN = dplyr::n(),
      PositiveMean = mean(hetero, na.rm = TRUE),
      PositiveMean_cts = mean(Freq, na.rm = TRUE),
      maxcts = max(Freq, na.rm = TRUE),
      TotalVcount = sum(Freq, na.rm = TRUE),
      CV = {
        m <- mean(hetero, na.rm = TRUE)
        if (is.na(m) || m == 0) NA_real_ else stats::sd(hetero, na.rm = TRUE) / m
      },
      .groups = "drop"
    )
  feature_cols <- c("CellN", "PositiveMean", "PositiveMean_cts", "maxcts", "CV", "TotalVcount", "meanCov")
  orig_cols <- names(ob@V.fitered)
  ob@V.fitered <- ob@V.fitered %>%
    dplyr::semi_join(gts, by = "Variants") %>%
    dplyr::select(-dplyr::any_of(feature_cols)) %>%
    dplyr::left_join(var_feat, by = "Variants") %>%
    dplyr::select(dplyr::any_of(orig_cols), dplyr::everything())

  # Update @CellMeta using semi_join based on Cell in GTsummary
  if ("CellMeta" %in% slotNames(ob) && is.data.frame(ob@CellMeta) && "Cell" %in% names(ob@CellMeta)) {
    ob@CellMeta <- ob@CellMeta %>%
      dplyr::semi_join(gts, by = "Cell")
  }

  # Update @UniqueV to unique variants from GTsummary
  if ("UniqueV" %in% slotNames(ob)) {
    ob@UniqueV <- unique(gts$Variants)
  }

  # Convert tibble to data frame to avoid Matrix compatibility issues
  # The Matrix package doesn't handle tbl_df objects well, causing errors like:
  # "*(<dgCMatrix>, <tbl_df>) is not yet implemented"
  ob@GTsummary.filtered <- as.data.frame(ob@GTsummary.filtered)

  # Rebuild matrices from filtered GTsummary
  # onlyhetero = TRUE ensures only heteroplasmic mutations are used
  ob <- Make_matrix(ob, onlyhetero = TRUE)
  if (update_depth_matrix) {
    ob <- Add_DepthMatrix_filter2(ob, QualifiedTotalCts)
  }
  ob
}

#' Filter by LIS score and VAF threshold (Wang et al., 2025, Genome Biology)
#'
#' @description
#' Keeps variants with VAF (\code{hetero}) ≥ \code{hetero_cutoff} and LIS > \code{lis_cutoff},
#' enforces a minimum number of cells per variant, updates \code{@GTsummary.filtered} and
#' \code{@V.fitered}, then rebuilds matrices and depth matrix.
#'
#' @param redeemR_obj A \code{redeemR} object.
#' @param lis_cutoff Numeric LIS threshold (default 0.6). LIS ≈ mean(hetero) / var(hetero).
#' @param hetero_cutoff Numeric VAF (hetero) threshold (default 0.05).
#' @param min_cells_per_variant Integer; minimum rows per variant (default 2).
#' @param QualifiedTotalCts Optional “QualifiedTotalCts” table for depth rebuilding.
#' @param filter_name Optional label passed to \code{print_redeemR_matrix_dims}.
#'
#' @return The updated \code{redeemR} object.
#'
#' @export
filter_redeemR_by_LIS <- function(redeemR_obj, lis_cutoff = 0.6, hetero_cutoff = 0.05, min_cells_per_variant = 2, QualifiedTotalCts = NULL, filter_name = NULL) {
  ob <- redeemR_obj
  gts <- ob@GTsummary.filtered
  gts <- gts %>% dplyr::filter(hetero > hetero_cutoff)
  lis_tbl <- gts %>%
    dplyr::group_by(Variants) %>%
    dplyr::summarise(
      mean_hetero = mean(hetero, na.rm = TRUE),
      var_hetero  = stats::var(hetero, na.rm = TRUE),
      LIS         = ifelse(is.na(var_hetero) | var_hetero == 0, NA_real_, mean_hetero / (1+var_hetero)),
      .groups = "drop"
    )

  keep_vars <- lis_tbl %>%
    dplyr::filter(!is.na(LIS) & LIS >= lis_cutoff) %>%
    dplyr::pull(Variants)

  gts2 <- gts %>% dplyr::filter(Variants %in% keep_vars) %>%
    dplyr::group_by(Variants) %>%
    dplyr::filter(dplyr::n() >= min_cells_per_variant) %>%
    dplyr::ungroup()
  ob@GTsummary.filtered <- gts2
  
  ob <- update_redeemR_from_GTsummary(ob, QualifiedTotalCts)
  
  # Print matrix dimensions
  print_redeemR_matrix_dims(ob, filter_name)
  
  ob
}

#' Filter by VAF (hetero) threshold on GTsummary
#'
#' @description
#' Keeps rows with \code{hetero >= hetero_cutoff} in \code{@GTsummary.filtered}, enforces a
#' minimum number of rows per variant, updates \code{@V.fitered}, and rebuilds matrices and depth matrix.
#'
#' @param redeemR_obj A \code{redeemR} object.
#' @param hetero_cutoff Numeric VAF threshold.
#' @param min_cells_per_variant Integer; minimum rows per variant (default 2).
#' @param QualifiedTotalCts Optional “QualifiedTotalCts” for depth rebuilding.
#' @param filter_name Optional label for matrix-dimension printing.
#'
#' @return The updated \code{redeemR} object.
#'
#' @export
filter_redeemR_by_hetero <- function(redeemR_obj, hetero_cutoff, min_cells_per_variant = 2, QualifiedTotalCts = NULL, filter_name = NULL) {
  ob <- redeemR_obj
  gts <- ob@GTsummary.filtered

  ob@GTsummary.filtered <- gts %>%
    dplyr::filter(hetero >= hetero_cutoff) %>%
    dplyr::group_by(Variants) %>%
    dplyr::filter(dplyr::n() >= min_cells_per_variant) %>%
    dplyr::ungroup()

  ob <- update_redeemR_from_GTsummary(ob, QualifiedTotalCts)
  
  # Print matrix dimensions
  print_redeemR_matrix_dims(ob, filter_name)
  
  ob
}

#' Filter by UMI count (Freq) threshold on GTsummary
#'
#' @description
#' Keeps rows with \code{Freq >= umi_cutoff} in \code{@GTsummary.filtered}, enforces a minimum
#' number of rows per variant, updates \code{@V.fitered}, and rebuilds matrices and depth matrix.
#'
#' @param redeemR_obj A \code{redeemR} object.
#' @param umi_cutoff Integer UMI threshold (default 2).
#' @param min_cells_per_variant Integer; minimum rows per variant (default 2).
#' @param QualifiedTotalCts Optional “QualifiedTotalCts” for depth rebuilding.
#' @param filter_name Optional label for matrix-dimension printing.
#'
#' @return The updated \code{redeemR} object.
#'
#' @export
filter_redeemR_by_UMI <- function(redeemR_obj, umi_cutoff = 2, min_cells_per_variant = 2, QualifiedTotalCts = NULL, filter_name = NULL) {
  ob <- redeemR_obj
  gts <- ob@GTsummary.filtered

  ob@GTsummary.filtered <- gts %>%
    dplyr::filter(Freq >= umi_cutoff) %>%
    dplyr::group_by(Variants) %>%
    dplyr::filter(dplyr::n() >= min_cells_per_variant) %>%
    dplyr::ungroup()

  ob <- update_redeemR_from_GTsummary(ob, QualifiedTotalCts)
  
  # Print matrix dimensions
  print_redeemR_matrix_dims(ob, filter_name)
  
  ob
}

#' Filter by cell subset
#'
#' @description
#' Keeps only rows in \code{@GTsummary.filtered} where the cell is in \code{cells_to_keep},
#' enforces a minimum number of rows per variant, updates \code{@V.fitered}, and rebuilds
#' matrices and depth matrix.
#'
#' @param redeemR_obj A \code{redeemR} object.
#' @param cells_to_keep Character vector of cell names to retain.
#' @param min_cells_per_variant Integer; minimum rows per variant (default 2). This is needed because after filterout some cells, some mutation may only show in 1 cell
#' @param QualifiedTotalCts Optional "QualifiedTotalCts" for depth rebuilding.
#' @param filter_name Optional label for matrix-dimension printing.
#'
#' @return The updated \code{redeemR} object.
#'
#' @export
filter_redeemR_by_cells <- function(redeemR_obj, cells_to_keep, min_cells_per_variant = 2, QualifiedTotalCts = NULL, filter_name = NULL, update_depth_matrix = TRUE) {
  ob <- redeemR_obj
  gts <- ob@GTsummary.filtered

  ob@GTsummary.filtered <- gts %>%
    dplyr::filter(Cell %in% cells_to_keep) %>%
    dplyr::group_by(Variants) %>%
    dplyr::filter(dplyr::n() >= min_cells_per_variant) %>%
    dplyr::ungroup()

  ob <- update_redeemR_from_GTsummary(ob, QualifiedTotalCts, update_depth_matrix = update_depth_matrix)
  
  # Print matrix dimensions
  print_redeemR_matrix_dims(ob, filter_name)
  
  ob
}

#' Filter by mean UMI count threshold on GTsummary
#'
#' @description
#' Keeps variants with \code{mean(Freq) >= mean_count_cutoff} across rows in \code{@GTsummary.filtered},
#' enforces a minimum number of rows per variant, updates \code{@V.fitered}, and rebuilds
#' matrices and depth matrix.
#'
#' @param redeemR_obj A \code{redeemR} object.
#' @param mean_count_cutoff Numeric mean UMI threshold (default 2).
#' @param min_cells_per_variant Integer; minimum rows per variant (default 2).
#' @param QualifiedTotalCts Optional “QualifiedTotalCts” for depth rebuilding.
#' @param filter_name Optional label for matrix-dimension printing.
#'
#' @return The updated \code{redeemR} object.
#'
#' @export
filter_redeemR_by_meancount <- function(redeemR_obj, mean_count_cutoff = 2, min_cells_per_variant = 2, QualifiedTotalCts = NULL, filter_name = NULL) {
  ob <- redeemR_obj
  gts <- ob@GTsummary.filtered

  # Calculate mean UMI count per variant
  variant_means <- gts %>%
    dplyr::group_by(Variants) %>%
    dplyr::summarise(
      mean_freq = mean(Freq, na.rm = TRUE),
      .groups = "drop"
    )

  # Keep variants with mean UMI count >= threshold
  keep_vars <- variant_means %>%
    dplyr::filter(mean_freq >= mean_count_cutoff) %>%
    dplyr::pull(Variants)

  # Filter GTsummary and enforce minimum rows per variant
  ob@GTsummary.filtered <- gts %>%
    dplyr::filter(Variants %in% keep_vars) %>%
    dplyr::group_by(Variants) %>%
    dplyr::filter(dplyr::n() >= min_cells_per_variant) %>%
    dplyr::ungroup()

  ob <- update_redeemR_from_GTsummary(ob, QualifiedTotalCts)
  
  # Print matrix dimensions
  print_redeemR_matrix_dims(ob, filter_name)
  
  ob
}

#' Filter by predefined rules using GTsummary data
#'
#' @description
#' Applies rule‑based filtering on \code{@GTsummary.filtered} using UMI counts and depth:
#' keeps rows with \code{Freq >= 2}; for \code{Freq == 1}, keeps only rows from \code{good_variants}
#' (variants meeting UMI‑based criteria) and with \code{depth <= max_depth}. Enforces minimum
#' rows per variant, then updates \code{@V.fitered} and rebuilds matrices.
#'
#' @param redeemR_obj A \code{redeemR} object.
#' @param rule Character; one of \code{"ruleA"}, \code{"ruleB"}, \code{"ruleC"}.
#' @param min_cells_per_variant Integer; minimum rows per variant (default 2).
#' @param QualifiedTotalCts Optional “QualifiedTotalCts” for depth rebuilding.
#' @param filter_name Optional label for matrix-dimension printing.
#'
#' @return The updated \code{redeemR} object.
#'
#' @export
filter_redeemR_by_rules <- function(redeemR_obj, rule = "ruleA", min_cells_per_variant = 2, QualifiedTotalCts = NULL, filter_name = NULL) {
  ob <- redeemR_obj
  gts <- ob@GTsummary.filtered

  # Define rule parameters
  rule_params <- list(
    ruleA = list(max_mean_umi_gt1 = 4, min_prop_2_3 = 0.3, max_depth = 70),
    ruleB = list(max_mean_umi_gt1 = 5, min_prop_2_3 = 0.1, max_depth = 100),
    ruleC = list(max_mean_umi_gt1 = 8, min_prop_2_3 = 0.1, max_depth = 150)
  )
  
  if (!rule %in% names(rule_params)) {
    stop("Invalid rule. Must be one of: ", paste(names(rule_params), collapse = ", "))
  }
  
  params <- rule_params[[rule]]
  
  # Calculate variant-level features from GTsummary Freq column (UMI counts)
  variant_features <- gts %>%
    dplyr::group_by(Variants) %>%
    dplyr::summarise(
      mean_umi_gt1 = {
        umi_gt1 <- Freq[Freq > 1]
        if (length(umi_gt1) == 0) 1 else mean(umi_gt1, na.rm = TRUE)
      },
      prop_2_3 = sum(Freq %in% c(2, 3)) / dplyr::n(),
      .groups = "drop"
    )
  
  # Identify good variants based on UMI criteria
  good_variants <- variant_features %>%
    dplyr::filter(
      (mean_umi_gt1 > 1 & 
       mean_umi_gt1 <= params$max_mean_umi_gt1 & 
       prop_2_3 > params$min_prop_2_3)
    ) %>%
    dplyr::pull(Variants)
  
  # Filter GTsummary: keep good variants regardless of UMI, but for non-good variants require Freq > 1
  # First, keep all rows with Freq > 2
  gt_high <- gts %>%
    dplyr::filter(Freq >= 2)
  
  # Next, for Freq == 1, keep only those with good_variants and depth <= rule_params$max_depth
  gt_one <- gts %>%
    dplyr::filter(
      Freq == 1,
      Variants %in% good_variants,
      depth <= params$max_depth
    )
  
    # Combine
  ob@GTsummary.filtered <- dplyr::bind_rows(gt_high, gt_one) %>%
    dplyr::group_by(Variants) %>%
    dplyr::filter(dplyr::n() >= min_cells_per_variant) %>%
    dplyr::ungroup()
  
  ob <- update_redeemR_from_GTsummary(ob, QualifiedTotalCts)
  
  # Print matrix dimensions
  print_redeemR_matrix_dims(ob, filter_name)
  
  ob
}

#' Subset redeemR object to a whitelist of cells and rebuild matrices
#'
#' @description
#' Filters \code{@GTsummary.filtered} to the provided set of cells, then
#' synchronizes \code{@V.fitered}, updates metadata, and rebuilds matrices and
#' the depth matrix via \code{update_redeemR_from_GTsummary}. Finally prints
#' matrix dimensions (optionally labeled) and returns the updated object.
#'
#' @param redeemR_obj A \code{redeemR} object.
#' @param cell_whitelist Character vector of cell barcodes to keep.
#' @param QualifiedTotalCts Optional “QualifiedTotalCts” table forwarded to
#'   \code{Add_DepthMatrix_filter2} during rebuild.
#' @param filter_name Optional label passed to \code{print_redeemR_matrix_dims}.
#'
#' @return The updated \code{redeemR} object restricted to the whitelist.
#'
#' @examples
#' # ob <- subset_redeem(ob, cell_whitelist = c("AAAC...", "AAAG..."))
#'
#' @export
#' @importFrom dplyr filter
subset_redeem <- function(redeemR_obj, cell_whitelist, QualifiedTotalCts = NULL, filter_name = NULL){
    redeemR_obj@GTsummary.filtered<-redeemR_obj@GTsummary.filtered %>% filter(Cell %in% cell_whitelist)
    redeemR_obj <- update_redeemR_from_GTsummary(redeemR_obj, QualifiedTotalCts)
    print_redeemR_matrix_dims(redeemR_obj, filter_name)
    return(redeemR_obj)
}