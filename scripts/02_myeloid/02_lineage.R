options(stringsAsFactors = FALSE)
set.seed(1234)
required_pkgs <- c( "Seurat", "SeuratObject", "dplyr", "ggplot2", "patchwork", "clustree", "AUCell", "GSVA", "pheatmap", "Matrix" )
missing_pkgs <- required_pkgs[ !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE) ]
if (length(missing_pkgs) > 0) {
  stop("Missing R packages: ", paste(missing_pkgs, collapse = ", "))
}
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(clustree)
  library(AUCell)
  library(GSVA)
  library(pheatmap)
  library(Matrix)
})
project_dir <- Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
step2A_dir <- file.path( project_dir, "04_results", "step2A_GSE162610_myeloid_refinement" )
step2A_rds <- file.path( step2A_dir, "GSE162610_step2A_myeloid_refined_multires.rds" )
step1F_dir <- file.path( project_dir, "04_results", "step1F_GSE162610_broad_annotation" )
step1F_rds <- file.path( step1F_dir, "GSE162610_broad_annotation_audit.rds" )
out_dir <- file.path( project_dir, "04_results", "step2B_GSE162610_lineage_adjudication" )
fig_dir <- file.path(out_dir, "figures")
tab_dir <- file.path(out_dir, "tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(step2A_rds)) {
  stop("Required input not found: ", step2A_rds)
}
if (!file.exists(step1F_rds)) {
  stop("Required input not found: ", step1F_rds)
}
obj <- readRDS(step2A_rds)
if (!inherits(obj, "Seurat")) {
  stop("Input object is not a Seurat object.")
}
DefaultAssay(obj) <- "RNA"
meta <- obj[[]]
meta_names <- colnames(meta)
if (!"myeloid_cluster_r04" %in% meta_names) {
  stop("Required input not found: ")
}
obj$myeloid_cluster_r04 <- as.character(obj$myeloid_cluster_r04)
Idents(obj) <- "myeloid_cluster_r04"
print(table(obj$myeloid_cluster_r04))
expected_clusters <- as.character(0:11)
observed_clusters <- sort(unique(as.character(obj$myeloid_cluster_r04)))
if (!all(expected_clusters %in% observed_clusters)) {
  warning( "Validation check.", paste(observed_clusters, collapse = ", ") )
}
intended_resolutions <- c(0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0)
r04_candidates <- grep( "res\\.0\\.4$", meta_names, value = TRUE )
if (length(r04_candidates) == 0) {
  stop("Required input not found: ")
}
target_r04 <- as.character(obj$myeloid_cluster_r04)
r04_match <- vapply( r04_candidates, function(cc) {
    x <- as.character(meta[[cc]])
    mean(x == target_r04, na.rm = TRUE)
  },
  numeric(1)
)
best_r04_col <- r04_candidates[which.max(r04_match)]
best_match <- max(r04_match)
if (!is.finite(best_match) || best_match < 0.95) {
  stop( "Validation failed.", "Best candidate: ", best_r04_col, " | agreement = ", round(best_match, 4) )
}
cluster_prefix <- sub("0\\.4$", "", best_r04_col)
find_resolution_col <- function(resolution_value) {
  candidates <- meta_names[startsWith(meta_names, cluster_prefix)]
  suffix <- sub(cluster_prefix, "", candidates, fixed = TRUE)
  suffix_num <- suppressWarnings(as.numeric(suffix))
  hit <- candidates[ is.finite(suffix_num) & abs(suffix_num - resolution_value) < 1e-8 ]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}
res_source_cols <- vapply( intended_resolutions, find_resolution_col, character(1) )
res_map <- data.frame( resolution = intended_resolutions, source_column = res_source_cols, stringsAsFactors = FALSE )
if (any(is.na(res_source_cols))) {
  stop( "Required input not found: ", paste( intended_resolutions[is.na(res_source_cols)], collapse = ", " ) )
}
tree_meta <- data.frame( row.names = rownames(meta) )
for (i in seq_along(intended_resolutions)) {
  r <- intended_resolutions[i]
  out_name <- paste0( "step2B_res.", format(r, trim = TRUE, scientific = FALSE) )
  tree_meta[[out_name]] <- as.factor( as.character(meta[[res_source_cols[i]]]) )
}
cluster_number_tbl <- data.frame( resolution = intended_resolutions, n_clusters = vapply( tree_meta, function(x) length(unique(as.character(x))), integer(1) ), source_column = res_source_cols, stringsAsFactors = FALSE )
write.csv( cluster_number_tbl, file.path(tab_dir, "01_corrected_cluster_number_by_resolution.csv"), row.names = FALSE )
p_tree <- clustree( tree_meta, prefix = "step2B_res." )
ggsave( file.path( fig_dir, "01_corrected_myeloid_clustree_Step2A_resolutions_only.pdf" ), p_tree, width = 11, height = 8 )
target_cells <- colnames(obj)
author_col <- "NOT_FOUND"
author_source_file <- "NOT_FOUND"
author_source_type <- "NOT_FOUND"
author_match_method <- "NOT_FOUND"
author_vec_final <- setNames( rep(NA_character_, length(target_cells)), target_cells )
detect_author_column <- function(df) {
  nm <- colnames(df)
  if (is.null(nm)) return(NA_character_)
  exact_priority <- c( "original_V30", "Original_V30", "V30", "v30" )
  hit <- exact_priority[exact_priority %in% nm]
  if (length(hit) > 0) return(hit[1])
  hit <- grep( "original.*V30|V30.*original|(^|[._])V30($|[._])", nm, value = TRUE, ignore.case = TRUE )
  if (length(hit) > 0) return(hit[1])
  annotation_terms <- c( "Microglia", "Macrophage", "Monocyte", "Astrocyte", "Oligodend", "Endothelial", "Neutrophil", "Fibroblast", "Ependymal", "T cell", "B cell" )
  candidate_cols <- nm[ vapply( df, function(x) is.character(x) || is.factor(x), logical(1) ) ]
  if (length(candidate_cols) == 0) return(NA_character_)
  semantic_score <- vapply( candidate_cols, function(cc) {
      vals <- unique(as.character(df[[cc]]))
      vals <- vals[!is.na(vals)]
      if (length(vals) == 0) return(0)
      sum( vapply( annotation_terms, function(term) {
            any( grepl( term, vals, ignore.case = TRUE, fixed = FALSE ) )
          },
          logical(1)
        )
      )
    },
    numeric(1)
  )
  if (max(semantic_score, na.rm = TRUE) >= 2) {
    return(candidate_cols[which.max(semantic_score)])
  }
  NA_character_
}
find_best_cell_id_vector <- function( df, target_cells, external_barcodes = NULL ) {
  best_ids <- NULL
  best_method <- "NONE"
  best_overlap <- 0L
  rn <- rownames(df)
  if (!is.null(rn) && length(rn) == nrow(df)) {
    rn_chr <- as.character(rn)
    ov <- sum(rn_chr %in% target_cells)
    if (ov > best_overlap) {
      best_ids <- rn_chr
      best_method <- "rownames"
      best_overlap <- ov
    }
  }
  for (cc in colnames(df)) {
    x <- df[[cc]]
    if ( !(is.character(x) || is.factor(x)) && !grepl( "barcode|cell|id", cc, ignore.case = TRUE ) ) {
      next
    }
    vals <- as.character(x)
    if (length(vals) != nrow(df)) next
    ov <- sum(vals %in% target_cells, na.rm = TRUE)
    if (ov > best_overlap) {
      best_ids <- vals
      best_method <- paste0("column:", cc)
      best_overlap <- ov
    }
  }
  if ( !is.null(external_barcodes) && length(external_barcodes) == nrow(df) ) {
    ext <- as.character(external_barcodes)
    ov <- sum(ext %in% target_cells)
    if (ov > best_overlap) {
      best_ids <- ext
      best_method <- "external_barcode_file_by_row_order"
      best_overlap <- ov
    }
  }
  list( ids = best_ids, method = best_method, overlap = best_overlap )
}
recover_author_from_df <- function( df, source_file, source_type, target_cells, external_barcodes = NULL ) {
  if (!is.data.frame(df)) {
    df <- as.data.frame(df)
  }
  ann_col <- detect_author_column(df)
  if (is.na(ann_col)) {
    return( list( success = FALSE, reason = "annotation_column_not_found" ) )
  }
  id_info <- find_best_cell_id_vector( df = df, target_cells = target_cells, external_barcodes = external_barcodes )
  if (is.null(id_info$ids) || id_info$overlap == 0) {
    return( list( success = FALSE, reason = "cell_barcode_match_not_found", annotation_col = ann_col ) )
  }
  ids <- as.character(id_info$ids)
  ann <- as.character(df[[ann_col]])
  good <- !is.na(ids) & ids != ""
  ids <- ids[good]
  ann <- ann[good]
  keep_unique <- !duplicated(ids)
  ids <- ids[keep_unique]
  ann <- ann[keep_unique]
  ann_map <- setNames(ann, ids)
  recovered <- unname( ann_map[target_cells] )
  non_na <- sum( !is.na(recovered) & recovered != "" )
  list( success = non_na > 0, recovered = recovered, annotation_col = ann_col, match_method = id_info$method, direct_overlap = id_info$overlap, recovered_nonNA = non_na, source_file = source_file, source_type = source_type )
}
safe_read_csv <- function(path) {
  tryCatch( read.csv( path, stringsAsFactors = FALSE, check.names = FALSE ), error = function(e) {
      NULL
    }
  )
}
safe_read_tsv_gz <- function(path) {
  tryCatch( read.delim( gzfile(path), header = TRUE, stringsAsFactors = FALSE, check.names = FALSE ), error = function(e) {
      NULL
    }
  )
}
raw_barcode_file <- file.path( project_dir, "01_raw_data", "GSE162610", "GSE162610_barcodes.tsv.gz" )
geo_barcodes <- NULL
if (file.exists(raw_barcode_file)) {
  geo_barcodes <- tryCatch( {
      x <- read.delim( gzfile(raw_barcode_file), header = FALSE, stringsAsFactors = FALSE, check.names = FALSE )
      as.character(x[[1]])
    },
    error = function(e) {
      NULL
    }
  )
}
recovery_attempts <- list()
obj1F <- tryCatch( readRDS(step1F_rds), error = function(e) NULL )
if (!is.null(obj1F) && inherits(obj1F, "Seurat")) {
  meta1F <- obj1F[[]]
  result_A <- recover_author_from_df( df = meta1F, source_file = step1F_rds, source_type = "Step1F_RDS_metadata", target_cells = target_cells, external_barcodes = NULL )
  recovery_attempts[["A_Step1F"]] <- result_A
  rm(obj1F, meta1F)
  gc()
} else {
  recovery_attempts[["A_Step1F"]] <- list( success = FALSE, reason = "Step1F_RDS_unreadable" )
}
step0_candidate_files <- c( file.path( project_dir, "00_admin", "step0_GSE162610_audit", "11_GSE162610_cell_metadata_with_auto_fields.csv" ), file.path( project_dir, "00_admin", "step0_GSE162610_confirmed", "01_GSE162610_cell_metadata_CONFIRMED.csv" ), file.path( project_dir, "00_admin", "step0_GSE162610_audit", "08_cell_metadata_merged.csv" ) )
step0_candidate_types <- c( "Step0_auto_fields_metadata", "Step0_confirmed_metadata", "Step0_merged_metadata" )
for (ii in seq_along(step0_candidate_files)) {
  ff <- step0_candidate_files[ii]
  if (!file.exists(ff)) {
    recovery_attempts[[paste0("Step0_", ii)]] <- list( success = FALSE, reason = "file_not_found", source_file = ff )
    next
  }
  df_tmp <- safe_read_csv(ff)
  if (is.null(df_tmp)) {
    recovery_attempts[[paste0("Step0_", ii)]] <- list( success = FALSE, reason = "read_failed", source_file = ff )
    next
  }
  result_tmp <- recover_author_from_df( df = df_tmp, source_file = ff, source_type = step0_candidate_types[ii], target_cells = target_cells, external_barcodes = geo_barcodes )
  recovery_attempts[[paste0("Step0_", ii)]] <- result_tmp
  rm(df_tmp)
  gc()
}
raw_meta_file <- file.path( project_dir, "01_raw_data", "GSE162610", "GSE162610_barcode_metadata.tsv.gz" )
if (file.exists(raw_meta_file)) {
  raw_meta <- safe_read_tsv_gz( raw_meta_file )
  if (!is.null(raw_meta)) {
    result_E <- recover_author_from_df( df = raw_meta, source_file = raw_meta_file, source_type = "GEO_raw_barcode_metadata", target_cells = target_cells, external_barcodes = geo_barcodes )
    recovery_attempts[["E_GEO_raw"]] <- result_E
    rm(raw_meta)
    gc()
  } else {
    recovery_attempts[["E_GEO_raw"]] <- list( success = FALSE, reason = "raw_metadata_read_failed" )
  }
} else {
  recovery_attempts[["E_GEO_raw"]] <- list( success = FALSE, reason = "raw_metadata_file_not_found", source_file = raw_meta_file )
}
success_idx <- which( vapply( recovery_attempts, function(x) isTRUE(x$success), logical(1) ) )
if (length(success_idx) > 0) {
  success_results <- recovery_attempts[ success_idx ]
  coverage_values <- vapply( success_results, function(x) {
      if (is.null(x$recovered_nonNA)) 0 else x$recovered_nonNA
    },
    numeric(1)
  )
  best_name <- names(success_results)[ which.max(coverage_values) ]
  best_result <- success_results[[best_name]]
  author_vec_final <- best_result$recovered
  names(author_vec_final) <- target_cells
  author_col <- best_result$annotation_col
  author_source_file <- best_result$source_file
  author_source_type <- best_result$source_type
  author_match_method <- best_result$match_method
} else {
  warning( paste0( "original_V30 could not be recovered from Step1F, Step0, or raw GEO metadata.\n", "Step 2B will CONTINUE with marker + AUCell + ssGSEA + broad-origin evidence.\n", "Author annotation cross-check will be recorded as unavailable." ) )
}
obj$original_V30_step1F <- unname( author_vec_final[colnames(obj)] )
coverage_n <- sum( !is.na(obj$original_V30_step1F) & obj$original_V30_step1F != "" )
coverage_pct <- ( coverage_n / ncol(obj) * 100 )
author_recovery_summary <- data.frame( source_type = author_source_type, source_file = author_source_file, source_column = author_col, cell_match_method = author_match_method, step2B_cells = ncol(obj), recovered_nonNA = coverage_n, recovered_percent = coverage_pct, stringsAsFactors = FALSE )
write.csv( author_recovery_summary, file.path( tab_dir, "02_originalV30_recovery_summary.csv" ), row.names = FALSE )
attempt_summary <- do.call( rbind, lapply( names(recovery_attempts), function(nm) {
      x <- recovery_attempts[[nm]]
      data.frame( attempt = nm, success = ifelse( is.null(x$success), FALSE, isTRUE(x$success) ), source_file = ifelse( is.null(x$source_file), NA_character_, x$source_file ), annotation_col = ifelse( is.null(x$annotation_col), NA_character_, x$annotation_col ), match_method = ifelse( is.null(x$match_method), NA_character_, x$match_method ), direct_overlap = ifelse( is.null(x$direct_overlap), NA_real_, x$direct_overlap ), recovered_nonNA = ifelse( is.null(x$recovered_nonNA), NA_real_, x$recovered_nonNA ), reason = ifelse( is.null(x$reason), NA_character_, x$reason ), stringsAsFactors = FALSE )
    }
  )
)
write.csv( attempt_summary, file.path( tab_dir, "02B_originalV30_recovery_attempts.csv" ), row.names = FALSE )
author_cross <- as.data.frame( table( myeloid_cluster_r04 = obj$myeloid_cluster_r04, original_V30 = obj$original_V30_step1F, useNA = "ifany" ) )
author_cross <- author_cross %>% group_by(myeloid_cluster_r04) %>% mutate( cluster_total = sum(Freq), proportion_within_cluster = Freq / cluster_total ) %>% ungroup()
write.csv( author_cross, file.path( tab_dir, "03_cluster_vs_originalV30_step1F.csv" ), row.names = FALSE )
if (coverage_n > 0) {
  author_tab <- table( obj$myeloid_cluster_r04, obj$original_V30_step1F, useNA = "no" )
  if ( nrow(author_tab) > 0 && ncol(author_tab) > 0 ) {
    author_prop <- prop.table( author_tab, margin = 1 )
    pdf( file.path( fig_dir, "02_cluster_vs_originalV30_heatmap.pdf" ), width = max( 8, ncol(author_prop) * 0.55 + 4 ), height = 7 )
    pheatmap( author_prop, cluster_rows = FALSE, cluster_cols = FALSE, main = "Step 2B cluster vs original_V30 | row proportion", fontsize = 9 )
    dev.off()
  }
} else {
}
gc()
gene_sets <- list( Resident_microglia = c( "P2ry12", "Tmem119", "Sall1", "Fcrls", "Hexb", "Cx3cr1", "Siglech", "Gpr34", "P2ry13" ), Infiltrating_mono_macro = c( "Ccr2", "Ly6c2", "Ms4a7", "Fn1", "Fcgr1", "S100a8", "S100a9", "Chil3", "Clec4n", "F13a1", "Arg1" ), BAM_like = c( "Lyve1", "Cd163", "Mrc1", "Mgl2", "Clec10a", "Pf4" ), Cycling_S = c( "Ccne2", "Hells", "Uhrf1", "Dtl", "Chek1", "Rad51", "Mcm10", "Mcm5", "Cdc6", "Mcm2" ), Cycling_G2M = c( "Mki67", "Top2a", "Ube2c", "Birc5", "Cenpf", "Ccna2", "Aurkb" ), IFN_response = c( "Ifit1", "Ifit2", "Ifit3", "Isg15", "Mx1", "Rsad2", "Usp18", "Oasl1" ), MHCII_APC = c( "H2-Aa", "H2-Ab1", "H2-Eb1", "Cd74", "Itgax" ) )
rna_genes <- rownames(obj[["RNA"]])
gene_sets_use <- lapply( gene_sets, function(gs) intersect(gs, rna_genes) )
gene_set_availability <- do.call( rbind, lapply( names(gene_sets_use), function(nm) {
      data.frame( gene_set = nm, n_requested = length(gene_sets[[nm]]), n_available = length(gene_sets_use[[nm]]), available_genes = paste(gene_sets_use[[nm]], collapse = ";"), stringsAsFactors = FALSE )
    }
  )
)
write.csv( gene_set_availability, file.path(tab_dir, "04_lineage_gene_set_availability.csv"), row.names = FALSE )
bad_sets <- names(gene_sets_use)[ vapply(gene_sets_use, length, integer(1)) < 3 ]
if (length(bad_sets) > 0) {
  stop( "Threshold check failed.", paste(bad_sets, collapse = ", ") )
}
dot_features <- unique(c( gene_sets_use$Resident_microglia, gene_sets_use$Infiltrating_mono_macro, gene_sets_use$BAM_like, gene_sets_use$Cycling_S, gene_sets_use$Cycling_G2M ))
Idents(obj) <- "myeloid_cluster_r04"
p_dot <- DotPlot( obj, features = dot_features, group.by = "myeloid_cluster_r04", assay = "RNA" ) + RotatedAxis() + ggtitle( "Step 2B lineage adjudication | resident vs infiltrating vs BAM vs cycling" )
ggsave( file.path( fig_dir, "03_expanded_lineage_marker_DotPlot.pdf" ), p_dot, width = 19, height = 7 )
get_rna_matrix <- function(seurat_obj) {
  x <- tryCatch( GetAssayData( seurat_obj, assay = "RNA", layer = "counts" ), error = function(e) NULL )
  if (is.null(x)) {
    x <- tryCatch( GetAssayData( seurat_obj, assay = "RNA", slot = "counts" ), error = function(e) NULL )
  }
  if (is.null(x)) {
    x <- tryCatch( GetAssayData( seurat_obj, assay = "RNA", layer = "data" ), error = function(e) NULL )
  }
  if (is.null(x)) {
    stop("Validation failed.")
  }
  x
}
expr_auc <- get_rna_matrix(obj)
if (!inherits(expr_auc, "dgCMatrix")) {
  expr_auc <- as(expr_auc, "dgCMatrix")
}
keep_gene_auc <- Matrix::rowSums(expr_auc) > 0
expr_auc <- expr_auc[keep_gene_auc, , drop = FALSE]
cells_rankings <- AUCell_buildRankings( expr_auc, plotStats = FALSE, nCores = 1, verbose = TRUE )
auc_max_rank <- max( 50, ceiling(0.05 * nrow(expr_auc)) )
cells_AUC <- tryCatch( AUCell_calcAUC( gene_sets_use, cells_rankings, aucMaxRank = auc_max_rank, nCores = 1 ), error = function(e1) {
    AUCell_calcAUC( gene_sets_use, cells_rankings, aucMaxRank = auc_max_rank )
  }
)
auc_mat <- as.matrix(getAUC(cells_AUC))
auc_mat <- auc_mat[, colnames(obj), drop = FALSE]
for (nm in rownames(auc_mat)) {
  col_nm <- paste0("AUCell_", nm)
  obj[[col_nm]] <- as.numeric(auc_mat[nm, ])
}
auc_cols <- paste0("AUCell_", rownames(auc_mat))
auc_summary_list <- lapply( sort(unique(as.character(obj$myeloid_cluster_r04))), function(cl) {
    idx <- which(obj$myeloid_cluster_r04 == cl)
    out <- data.frame( myeloid_cluster_r04 = cl, n_cells = length(idx), stringsAsFactors = FALSE )
    for (nm in rownames(auc_mat)) {
      vals <- auc_mat[nm, idx]
      out[[paste0(nm, "_mean")]] <- mean(vals, na.rm = TRUE)
      out[[paste0(nm, "_median")]] <- median(vals, na.rm = TRUE)
    }
    out
  }
)
auc_cluster_summary <- bind_rows(auc_summary_list)
write.csv( auc_cluster_summary, file.path(tab_dir, "05_AUCell_cluster_summary.csv"), row.names = FALSE )
cluster_levels <- sort(unique(as.character(obj$myeloid_cluster_r04)))
auc_cluster_median <- sapply( cluster_levels, function(cl) {
    idx <- which(obj$myeloid_cluster_r04 == cl)
    apply( auc_mat[, idx, drop = FALSE], 1, median, na.rm = TRUE )
  }
)
colnames(auc_cluster_median) <- cluster_levels
row_z <- function(x) {
  sx <- sd(x, na.rm = TRUE)
  if (!is.finite(sx) || sx == 0) return(rep(0, length(x)))
  (x - mean(x, na.rm = TRUE)) / sx
}
auc_cluster_z <- t( apply( auc_cluster_median, 1, row_z ) )
pdf( file.path( fig_dir, "04_AUCell_cluster_signature_heatmap.pdf" ), width = 10, height = 6 )
pheatmap( auc_cluster_z, cluster_rows = FALSE, cluster_cols = FALSE, main = "AUCell signature activity | cluster median Z-score" )
dev.off()
lineage_auc_features <- c( "AUCell_Resident_microglia", "AUCell_Infiltrating_mono_macro", "AUCell_BAM_like" )
p_auc_umap <- FeaturePlot( obj, features = lineage_auc_features, reduction = "umap.myeloid", ncol = 3, order = TRUE )
ggsave( file.path( fig_dir, "05_AUCell_lineage_signature_UMAP.pdf" ), p_auc_umap, width = 15, height = 5 )
focus_clusters <- c("0", "3", "4", "8", "1", "2", "6")
plot_meta <- obj[[]]
plot_meta <- plot_meta[ as.character(plot_meta$myeloid_cluster_r04) %in% focus_clusters, , drop = FALSE ]
build_long_auc <- function(df, cols, cluster_col) {
  out <- vector("list", length(cols))
  for (i in seq_along(cols)) {
    cc <- cols[i]
    out[[i]] <- data.frame( cluster = as.character(df[[cluster_col]]), signature = sub("^AUCell_", "", cc), AUC = as.numeric(df[[cc]]), stringsAsFactors = FALSE )
  }
  bind_rows(out)
}
focus_long <- build_long_auc( plot_meta, lineage_auc_features, "myeloid_cluster_r04" )
focus_long$cluster <- factor( focus_long$cluster, levels = focus_clusters )
p_focus <- ggplot( focus_long, aes(x = cluster, y = AUC) ) + geom_violin( scale = "width", trim = TRUE ) + geom_boxplot( width = 0.15, outlier.shape = NA ) + facet_wrap( ~ signature, scales = "free_y", ncol = 3 ) + theme_bw() + labs( title = "Focused AUCell adjudication of cluster 4 and cluster 8", x = "myeloid cluster", y = "AUCell AUC" )
ggsave( file.path( fig_dir, "06_cluster4_8_AUCell_violin.pdf" ), p_focus, width = 13, height = 5 )
avg_expr <- AverageExpression( obj, assays = "RNA", group.by = "myeloid_cluster_r04", layer = "data", verbose = FALSE )[["RNA"]]
avg_expr <- as.matrix(avg_expr)
run_ssgsea_compat <- function(expr_matrix, gene_sets_list) {
  ns_exports <- getNamespaceExports("GSVA")
  if (any(ns_exports == "ssgseaParam")) {
    param <- tryCatch( GSVA::ssgseaParam( exprData = expr_matrix, geneSets = gene_sets_list, minSize = 2, maxSize = Inf ), error = function(e1) {
        GSVA::ssgseaParam( expr_matrix, gene_sets_list, minSize = 2, maxSize = Inf )
      }
    )
    out <- GSVA::gsva( param, verbose = FALSE )
  } else {
    out <- GSVA::gsva( expr_matrix, gene_sets_list, method = "ssgsea", kcdf = "Gaussian", abs.ranking = TRUE, verbose = FALSE )
  }
  as.matrix(out)
}
ssgsea_mat <- run_ssgsea_compat( avg_expr, gene_sets_use )
write.csv( data.frame( gene_set = rownames(ssgsea_mat), ssgsea_mat, check.names = FALSE ), file.path(tab_dir, "06_pseudobulk_ssGSEA_scores.csv"), row.names = FALSE )
ssgsea_z <- t( apply( ssgsea_mat, 1, row_z ) )
pdf( file.path( fig_dir, "07_pseudobulk_ssGSEA_signature_heatmap.pdf" ), width = 10, height = 6 )
pheatmap( ssgsea_z, cluster_rows = FALSE, cluster_cols = FALSE, main = "Cluster-average ssGSEA | row Z-score" )
dev.off()
get_cluster_auc_z <- function(signature, cluster_id) {
  if (!signature %in% rownames(auc_cluster_z)) return(NA_real_)
  if (!cluster_id %in% colnames(auc_cluster_z)) return(NA_real_)
  auc_cluster_z[signature, cluster_id]
}
get_cluster_ssgsea_z <- function(signature, cluster_id) {
  if (!signature %in% rownames(ssgsea_z)) return(NA_real_)
  candidates <- c( cluster_id, paste0("g", cluster_id) )
  hit <- candidates[candidates %in% colnames(ssgsea_z)]
  if (length(hit) == 0) return(NA_real_)
  ssgsea_z[signature, hit[1]]
}
get_broad_origin_prop <- function(cluster_id, broad_id) {
  if (!"broad_origin" %in% colnames(obj[[]])) return(NA_real_)
  x <- obj[[]]
  x <- x[x$myeloid_cluster_r04 == cluster_id, , drop = FALSE]
  if (nrow(x) == 0) return(NA_real_)
  mean(as.character(x$broad_origin) == broad_id, na.rm = TRUE)
}
get_author_top <- function(cluster_id) {
  x <- obj[[]]
  x <- x[ x$myeloid_cluster_r04 == cluster_id & !is.na(x$original_V30_step1F), , drop = FALSE ]
  if (nrow(x) == 0) {
    return( list( label = NA_character_, prop = NA_real_ ) )
  }
  tt <- sort( table(x$original_V30_step1F), decreasing = TRUE )
  list( label = names(tt)[1], prop = as.numeric(tt[1]) / sum(tt) )
}
evidence_clusters <- c("4", "8")
evidence_tbl <- do.call( rbind, lapply( evidence_clusters, function(cl) {
      author_info <- get_author_top(cl)
      data.frame( cluster = cl, AUCell_resident_z = get_cluster_auc_z("Resident_microglia", cl), AUCell_infiltrating_z = get_cluster_auc_z("Infiltrating_mono_macro", cl), AUCell_BAM_z = get_cluster_auc_z("BAM_like", cl), ssGSEA_resident_z = get_cluster_ssgsea_z("Resident_microglia", cl), ssGSEA_infiltrating_z = get_cluster_ssgsea_z("Infiltrating_mono_macro", cl), ssGSEA_BAM_z = get_cluster_ssgsea_z("BAM_like", cl), broad_origin_0_prop = get_broad_origin_prop(cl, "0"), broad_origin_4_prop = get_broad_origin_prop(cl, "4"), broad_origin_2_prop = get_broad_origin_prop(cl, "2"), top_original_V30 = author_info$label, top_original_V30_prop = author_info$prop, final_lineage_call = NA_character_, reviewer_note = NA_character_, stringsAsFactors = FALSE )
    }
  )
)
write.csv( evidence_tbl, file.path( tab_dir, "07_cluster4_8_adjudication_evidence.csv" ), row.names = FALSE )
provisional_lineage <- data.frame( cluster = as.character(0:11), step2A_provisional_call = c( "clean_resident_microglia", "infiltrating_macrophage", "infiltrating_monocyte", "clean_activated_resident", "AMBIGUOUS_review_in_Step2B", "clean_cycling_S_microglia", "macrophage_biased_infiltrating", "clean_cycling_G2M_microglia", "AMBIGUOUS_BAM_like_review", "clean_IFN_resident_microglia", "IFN_infiltrating_mono_macro", "MHCII_APC_infiltrating_macro" ), final_step2_call = NA_character_, stringsAsFactors = FALSE )
write.csv( provisional_lineage, file.path( tab_dir, "08_provisional_lineage_framework_for_review.csv" ), row.names = FALSE )
step2B_rds <- file.path( out_dir, "GSE162610_step2B_lineage_adjudication.rds" )
saveRDS( obj, step2B_rds )
report_file <- file.path( out_dir, "09_STEP2B_lineage_adjudication_report.txt" )
sink(report_file)
print(table(obj$myeloid_cluster_r04))
print(cluster_number_tbl)
print(author_recovery_summary)
print(gene_set_availability)
print(evidence_tbl)
print(provisional_lineage)
sink()
