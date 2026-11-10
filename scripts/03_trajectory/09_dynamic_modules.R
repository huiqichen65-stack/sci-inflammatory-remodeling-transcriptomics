options(stringsAsFactors = FALSE)
set.seed(20260801)
PROJECT_ROOT <- Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
STEP5CR_DIR <- file.path( PROJECT_ROOT, "04_results", "step5CR_GSE162610_Y29_prebranch_progressive_divergence" )
STEP6A_DIR <- file.path( PROJECT_ROOT, "04_results", "step6A_GSE162610_7dpi_branch_specific_DEG" )
STEP5BR2_DIR <- file.path( PROJECT_ROOT, "04_results", "step5BR2_GSE162610_single_root_branch_audit" )
INPUT_STEP5CR <- file.path( STEP5CR_DIR, "GSE162610_step5CR_results.rds" )
INPUT_STEP6A <- file.path( STEP6A_DIR, "GSE162610_step6A_7dpi_branch_DEG_results.rds" )
INPUT_SEURAT <- file.path( STEP5BR2_DIR, "GSE162610_step5BR2_single_root_seurat.rds" )
OUTDIR <- file.path( PROJECT_ROOT, "04_results", "step6B_GSE162610_dynamic_branch_modules" )
FIGDIR <- file.path(OUTDIR, "figures")
TABDIR <- file.path(OUTDIR, "tables")
dir.create(FIGDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABDIR, recursive = TRUE, showWarnings = FALSE)
OUT_FIG1 <- file.path( FIGDIR, "Fig6B_01_dynamic_module_heatmap.pdf" )
OUT_FIG2 <- file.path( FIGDIR, "Fig6B_02_module_centroid_trajectories.pdf" )
OUT_FIG3 <- file.path( FIGDIR, "Fig6B_03_7dpi_pseudobulk_module_scores.pdf" )
OUT_FIG4 <- file.path( FIGDIR, "Fig6B_04_module_evidence_and_replicate_effects.pdf" )
OUT_XLSX <- file.path( OUTDIR, "Table6B_dynamic_branch_modules.xlsx" )
OUT_RDS <- file.path( OUTDIR, "GSE162610_step6B_dynamic_modules_results.rds" )
OUT_REP <- file.path( OUTDIR, "Step6B_report.txt" )
OUT_SESS <- file.path( OUTDIR, "Step6B_sessionInfo.txt" )
N_PRE_BINS <- 3L
N_POST_BINS <- 8L
MIN_DETECTION_FRACTION <- 0.05
K_CANDIDATES <- 4:8
MIN_MODULE_SIZE <- 15L
KMEANS_NSTART <- 50L
KMEANS_ITER_MAX <- 200L
TOP_HEAT_GENES_PER_MODULE <- 25L
TOP_REPORT_GENES_PER_MODULE <- 15L
MIN_ABS_MODULE_EFFECT_PER_SAMPLE <- 0.15
required_pkgs <- c( "Seurat", "SeuratObject", "Matrix", "dplyr", "tidyr", "tibble", "ggplot2", "patchwork", "cluster", "openxlsx" )
missing_pkgs <- required_pkgs[ !vapply( required_pkgs, requireNamespace, logical(1), quietly = TRUE ) ]
if (length(missing_pkgs) > 0) {
  stop( paste0( "缺少以下R包：", paste(missing_pkgs, collapse = ", ") ), call. = FALSE )
}
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(cluster)
})
capture_table <- function(x) {
  paste( capture.output( print( tibble::as_tibble(x), n = Inf, width = Inf ) ), collapse = "\n" )
}
get_norm <- function(obj) {
  assay <- DefaultAssay(obj)
  x <- tryCatch( Seurat::GetAssayData( obj, assay = assay, layer = "data" ), error = function(e) NULL )
  if (!is.null(x)) return(x)
  Seurat::GetAssayData( obj, assay = assay, slot = "data" )
}
safe_sheet <- function(wb, sheet_name, df) {
  openxlsx::addWorksheet( wb, sheetName = sheet_name )
  ddf <- tryCatch( as.data.frame(df), error = function(e) NULL )
  if ( is.null(ddf) || nrow(ddf) == 0 || ncol(ddf) == 0 ) {
    openxlsx::writeData( wb, sheet = sheet_name, x = data.frame( note = "No result under the predefined Step6B criterion.", stringsAsFactors = FALSE ) )
  } else {
    openxlsx::writeDataTable( wb, sheet = sheet_name, x = ddf, tableStyle = "TableStyleMedium2" )
  }
  openxlsx::freezePane( wb, sheet = sheet_name, firstRow = TRUE )
}
row_zscore <- function(mat) {
  mu <- rowMeans(mat, na.rm = TRUE)
  ss <- apply(mat, 1, stats::sd, na.rm = TRUE)
  keep <- is.finite(ss) & ss > 0
  z <- matrix( NA_real_, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat) )
  z[keep, ] <- sweep( sweep( mat[keep, , drop = FALSE], 1, mu[keep], "-" ), 1, ss[keep], "/" )
  z
}
balanced_rank_bin <- function( primary_order, secondary_order, cell_id, n_bins) {
  n <- length(primary_order)
  if (n == 0) {
    return(integer(0))
  }
  if (n < n_bins) {
    stop( paste0( "A sample x segment group contains only ", n, " cells but requires ", n_bins, " bins." ), call. = FALSE )
  }
  ord <- order( primary_order, secondary_order, cell_id, na.last = TRUE )
  out <- rep(NA_integer_, n)
  out[ord] <- floor( (seq_along(ord) - 1) * n_bins / n ) + 1L
  as.integer(out)
}
mean_sparse_by_group <- function(expr_mat, group_vector) {
  group_vector <- factor(group_vector)
  design <- Matrix::sparse.model.matrix( ~ 0 + group_vector )
  colnames(design) <- levels(group_vector)
  sums <- expr_mat %*% design
  counts <- Matrix::colSums(design)
  means <- sweep( as.matrix(sums), 2, counts, "/" )
  means
}
module_label_from_centroid <- function( centroid, pre_cols, arm1_cols, arm2_cols) {
  arm1_early <- head(arm1_cols, 2)
  arm2_early <- head(arm2_cols, 2)
  arm1_mid <- arm1_cols[ seq.int( max(2, floor(length(arm1_cols) / 3)), max(3, ceiling(2 * length(arm1_cols) / 3)) ) ]
  arm2_mid <- arm2_cols[ seq.int( max(2, floor(length(arm2_cols) / 3)), max(3, ceiling(2 * length(arm2_cols) / 3)) ) ]
  arm1_late <- tail(arm1_cols, 3)
  arm2_late <- tail(arm2_cols, 3)
  pre_score <- mean(centroid[pre_cols], na.rm = TRUE)
  early_score <- mean( centroid[c(arm1_early, arm2_early)], na.rm = TRUE )
  mid_score <- mean( centroid[c(arm1_mid, arm2_mid)], na.rm = TRUE )
  late1 <- mean( centroid[arm1_late], na.rm = TRUE )
  late2 <- mean( centroid[arm2_late], na.rm = TRUE )
  scores <- c( Shared_prebranch = pre_score, Early_shared = early_score, Intermediate_shared = mid_score, Arm1_terminal = late1, Arm2_terminal = late2 )
  if ( max(late1, late2) == max(scores, na.rm = TRUE) && abs(late1 - late2) < 0.20 ) {
    label <- "Shared_late"
  } else {
    label <- names(which.max(scores))
  }
  list( label = label, pre_score = pre_score, early_score = early_score, mid_score = mid_score, arm1_late_score = late1, arm2_late_score = late2, arm1_minus_arm2_late = late1 - late2 )
}
safe_cor <- function(x, y) {
  if ( length(x) != length(y) || stats::sd(x, na.rm = TRUE) == 0 || stats::sd(y, na.rm = TRUE) == 0 ) {
    return(NA_real_)
  }
  suppressWarnings( stats::cor( x, y, method = "pearson", use = "pairwise.complete.obs" ) )
}
for (ff in c(INPUT_STEP5CR, INPUT_STEP6A, INPUT_SEURAT)) {
  if (!file.exists(ff)) {
    stop( paste0("找不到输入文件：", ff), call. = FALSE )
  }
}
step5cr <- readRDS(INPUT_STEP5CR)
step6a <- readRDS(INPUT_STEP6A)
obj <- readRDS(INPUT_SEURAT)
DefaultAssay(obj) <- "RNA"
required_step5 <- c( "prebranch_cells", "arm1_cells", "arm2_cells", "analysis_metadata", "progressive_shape_genes", "strong_progressive_genes" )
missing_step5 <- required_step5[ !required_step5 %in% names(step5cr) ]
if (length(missing_step5) > 0) {
  stop( paste0( "Step5C-R结果缺少：", paste(missing_step5, collapse = ", ") ), call. = FALSE )
}
required_step6a <- c( "high_confidence_branch_markers", "master_evidence", "pseudobulk_logCPM", "pseudobulk_metadata" )
missing_step6a <- required_step6a[ !required_step6a %in% names(step6a) ]
if (length(missing_step6a) > 0) {
  stop( paste0( "Step6A结果缺少：", paste(missing_step6a, collapse = ", ") ), call. = FALSE )
}
step5_progressive <- unique( as.character( step5cr$progressive_shape_genes ) )
step5_strong <- unique( as.character( step5cr$strong_progressive_genes ) )
step6a_highconf <- unique( as.character( step6a$high_confidence_branch_markers$gene_id ) )
candidate_genes <- union( step5_progressive, step6a_highconf )
analysis_meta <- as.data.frame( step5cr$analysis_metadata )
required_meta <- c( "cell", "state", "sample", "signed_progress", "post_progress" )
missing_meta <- required_meta[ !required_meta %in% colnames(analysis_meta) ]
if (length(missing_meta) > 0) {
  stop( paste0( "Step5C-R analysis_metadata缺少：", paste(missing_meta, collapse = ", ") ), call. = FALSE )
}
pt_candidates <- c( "retained_step5BR_pseudotime", "step5BR_pseudotime", "pseudotime" )
PT_TIE_COL <- pt_candidates[ pt_candidates %in% colnames(analysis_meta) ][1]
if ( length(PT_TIE_COL) == 0 || is.na(PT_TIE_COL) ) {
  warning( "No retained pseudotime column found; cell ID will be used to break graph-distance ties." )
  analysis_meta$step6B_tie_pseudotime <- 0
} else {
  analysis_meta$step6B_tie_pseudotime <- as.numeric( analysis_meta[[PT_TIE_COL]] )
}
analysis_meta <- analysis_meta |> dplyr::filter( state %in% c( "Shared_prebranch", "Arm_1", "Arm_2" ), cell %in% colnames(obj) )
analysis_cells <- analysis_meta$cell
if (length(unique(analysis_cells)) != nrow(analysis_meta)) {
  stop( "Step5C-R analysis_metadata存在重复cell。", call. = FALSE )
}
norm_mat <- get_norm(obj)
candidate_genes <- intersect( candidate_genes, rownames(norm_mat) )
expr_candidate <- norm_mat[ candidate_genes, analysis_cells, drop = FALSE ]
detection_fraction <- Matrix::rowMeans( expr_candidate > 0 )
keep_genes <- names(detection_fraction)[ is.finite(detection_fraction) & detection_fraction >= MIN_DETECTION_FRACTION ]
if (length(keep_genes) < 50) {
  stop( "表达过滤后候选动态基因少于50个。", call. = FALSE )
}
expr_use <- norm_mat[ keep_genes, analysis_cells, drop = FALSE ]
analysis_meta <- analysis_meta |> dplyr::mutate( segment = dplyr::case_when( state == "Shared_prebranch" ~ "Pre", state == "Arm_1" ~ "Arm1", state == "Arm_2" ~ "Arm2", TRUE ~ NA_character_ ), local_progress = dplyr::case_when( segment == "Pre" ~ signed_progress, segment %in% c("Arm1", "Arm2") ~ post_progress, TRUE ~ NA_real_ ) )
group_size_audit <- analysis_meta |> dplyr::count( sample, segment, name = "n_cells" ) |> dplyr::mutate( required_bins = ifelse( segment == "Pre", N_PRE_BINS, N_POST_BINS ), sufficient_cells = n_cells >= required_bins )
if (any(!group_size_audit$sufficient_cells)) {
  print(group_size_audit)
  stop( "At least one biological sample x trajectory segment has fewer cells than bins.", call. = FALSE )
}
analysis_meta <- analysis_meta |> dplyr::group_by( sample, segment ) |> dplyr::group_modify( ~{
      segment_current <- as.character(.y$segment[[1]])
      nbin <- if ( segment_current == "Pre" ) {
        N_PRE_BINS
      } else {
        N_POST_BINS
      }
      .x$bin <- balanced_rank_bin( primary_order = .x$local_progress, secondary_order = .x$step6B_tie_pseudotime, cell_id = .x$cell, n_bins = nbin )
      .x
    }
  ) |> dplyr::ungroup() |> dplyr::mutate(
    sample_segment_bin = paste( sample, segment, sprintf("%02d", bin), sep = "__" ), trajectory_bin = paste( segment, sprintf("%02d", bin), sep = "_" )
  )
sample_bin_summary <- analysis_meta |> dplyr::count( sample, segment, bin, trajectory_bin, name = "n_cells" )
expected_group_count <- sum( group_size_audit$required_bins )
if (nrow(sample_bin_summary) != expected_group_count) {
  print(sample_bin_summary)
  stop( "Sample-specific trajectory binning did not produce all expected bins.", call. = FALSE )
}
sample_group_means <- mean_sparse_by_group( expr_use, analysis_meta$sample_segment_bin )
sample_group_meta <- tibble::tibble( sample_segment_bin = colnames(sample_group_means) ) |> tidyr::separate( sample_segment_bin, into = c( "sample", "segment", "bin" ), sep = "__", remove = FALSE ) |> dplyr::mutate( bin = as.integer(bin), trajectory_bin = paste( segment, sprintf("%02d", bin), sep = "_" ) )
trajectory_order <- c( paste0( "Pre_", sprintf("%02d", seq_len(N_PRE_BINS)) ), paste0( "Arm1_", sprintf("%02d", seq_len(N_POST_BINS)) ), paste0( "Arm2_", sprintf("%02d", seq_len(N_POST_BINS)) ) )
profile_mat <- matrix( NA_real_, nrow = length(keep_genes), ncol = length(trajectory_order), dimnames = list( keep_genes, trajectory_order ) )
for (tb in trajectory_order) {
  cols_tb <- sample_group_meta$sample_segment_bin[ sample_group_meta$trajectory_bin == tb ]
  cols_tb <- intersect( cols_tb, colnames(sample_group_means) )
  if (length(cols_tb) == 0) next
  profile_mat[, tb] <- rowMeans( sample_group_means[ , cols_tb, drop = FALSE ], na.rm = TRUE )
}
empty_profile_bins <- colnames(profile_mat)[ colSums(is.finite(profile_mat)) == 0 ]
if (length(empty_profile_bins) > 0) {
  stop( paste0( "The following trajectory bins are empty after sample balancing: ", paste(empty_profile_bins, collapse = ", ") ), call. = FALSE )
}
valid_profile_genes <- rownames(profile_mat)[ apply( profile_mat, 1, function(v) {
      all(is.finite(v)) && stats::sd(v) > 0
    }
  )
]
profile_mat <- profile_mat[ valid_profile_genes, , drop = FALSE ]
profile_z <- row_zscore(profile_mat)
profile_z <- profile_z[ apply(profile_z, 1, function(v) all(is.finite(v))), , drop = FALSE ]
if (nrow(profile_z) < 50) {
  stop( "可用于动态模块聚类的完整表达轨迹基因少于50个。", call. = FALSE )
}
dist_mat <- stats::dist( profile_z, method = "euclidean" )
k_search <- list()
fit_by_k <- list()
for (kk in K_CANDIDATES) {
  km <- stats::kmeans( profile_z, centers = kk, nstart = KMEANS_NSTART, iter.max = KMEANS_ITER_MAX )
  module_sizes <- table(km$cluster)
  min_size <- min(module_sizes)
  sil <- cluster::silhouette( km$cluster, dist_mat )
  avg_sil <- mean( sil[, "sil_width"], na.rm = TRUE )
  valid_size <- min_size >= MIN_MODULE_SIZE
  k_search[[as.character(kk)]] <- tibble::tibble( K = kk, average_silhouette = avg_sil, minimum_module_size = min_size, all_modules_ge_minimum = valid_size )
  fit_by_k[[as.character(kk)]] <- km
}
k_search_tbl <- dplyr::bind_rows( k_search )
valid_k <- k_search_tbl |> dplyr::filter( all_modules_ge_minimum )
if (nrow(valid_k) > 0) {
  selected_k <- valid_k |> dplyr::arrange( dplyr::desc(average_silhouette), K ) |> dplyr::slice(1) |> dplyr::pull(K)
} else {
  selected_k <- k_search_tbl |> dplyr::arrange( dplyr::desc(average_silhouette), dplyr::desc(minimum_module_size), K ) |> dplyr::slice(1) |> dplyr::pull(K)
  warning( "No K satisfied the predefined minimum module size; selected the best available K." )
}
selected_k <- as.integer(selected_k)
selected_fit <- fit_by_k[[as.character(selected_k)]]
raw_module <- selected_fit$cluster
module_centroids <- lapply( sort(unique(raw_module)), function(mm) {
    colMeans( profile_z[ raw_module == mm, , drop = FALSE ], na.rm = TRUE )
  }
)
names(module_centroids) <- sort( unique(raw_module) )
pre_cols <- grep( "^Pre_", colnames(profile_z), value = TRUE )
arm1_cols <- grep( "^Arm1_", colnames(profile_z), value = TRUE )
arm2_cols <- grep( "^Arm2_", colnames(profile_z), value = TRUE )
module_annotation_list <- lapply( names(module_centroids), function(mm) {
    ann <- module_label_from_centroid( centroid = module_centroids[[mm]], pre_cols = pre_cols, arm1_cols = arm1_cols, arm2_cols = arm2_cols )
    tibble::tibble( raw_module = as.integer(mm), temporal_label = ann$label, pre_score = ann$pre_score, early_score = ann$early_score, intermediate_score = ann$mid_score, arm1_late_score = ann$arm1_late_score, arm2_late_score = ann$arm2_late_score, arm1_minus_arm2_late = ann$arm1_minus_arm2_late )
  }
)
module_annotation <- dplyr::bind_rows( module_annotation_list )
module_annotation <- module_annotation |> dplyr::mutate( module_name = paste0( "M", raw_module, "_", temporal_label ) )
module_name_map <- stats::setNames( module_annotation$module_name, module_annotation$raw_module )
module_assignment <- tibble::tibble( gene_id = rownames(profile_z), raw_module = as.integer(raw_module), module_name = unname( module_name_map[ as.character(raw_module) ] ) )
module_assignment$module_membership_cor <- vapply( seq_len(nrow(module_assignment)), function(ii) {
    gg <- module_assignment$gene_id[ii]
    mm <- as.character( module_assignment$raw_module[ii] )
    safe_cor( profile_z[gg, ], module_centroids[[mm]] )
  },
  numeric(1)
)
step5_master_raw <- as.data.frame( step5cr$master_table )
required_step5_master_cols <- c( "gene_id", "progressive_shape_q_value", "anchored_branch_q_value", "terminal_direction" )
missing_step5_master_cols <- required_step5_master_cols[ !required_step5_master_cols %in% colnames(step5_master_raw) ]
if (length(missing_step5_master_cols) > 0) {
  stop( paste0( "Step5C-R master_table缺少：", paste( missing_step5_master_cols, collapse = ", " ) ), call. = FALSE )
}
strong_col_candidates <- c( "strong_progressive_pattern", "strong_progressive" )
STRONG_PROGRESSIVE_COL <- strong_col_candidates[ strong_col_candidates %in% colnames(step5_master_raw) ][1]
if ( length(STRONG_PROGRESSIVE_COL) == 0 || is.na(STRONG_PROGRESSIVE_COL) ) {
  warning( "Step5C-R master_table没有strong-progressive标记列；将根据strong_progressive_genes重新构建。" )
  step5_master_raw$strong_progressive_pattern <- step5_master_raw$gene_id %in%
    step5_strong
} else {
  step5_master_raw$strong_progressive_pattern <- as.logical( step5_master_raw[[STRONG_PROGRESSIVE_COL]] )
}
step5_master <- step5_master_raw |> dplyr::select( gene_id, progressive_shape_q_value, anchored_branch_q_value, strong_progressive_pattern, terminal_direction )
step6a_master <- step6a$master_evidence |> dplyr::select( gene_id, cell_avg_log2FC, cell_BH_FDR, high_confidence_branch_marker, strict_pseudobulk_supported, final_direction )
gene_modules <- module_assignment |> dplyr::left_join( step5_master, by = "gene_id" ) |> dplyr::left_join( step6a_master, by = "gene_id" ) |> dplyr::mutate( source_Step5_progressive = gene_id %in% step5_progressive, source_Step5_strong = gene_id %in% step5_strong, source_Step6A_highconf = gene_id %in% step6a_highconf, evidence_source = dplyr::case_when( source_Step5_progressive & source_Step6A_highconf ~ "Both", source_Step5_progressive ~ "Step5C-R_only", source_Step6A_highconf ~ "Step6A_only", TRUE ~ "Other" ) )
module_sizes <- gene_modules |> dplyr::count( raw_module, module_name, name = "n_genes" )
module_annotation <- module_annotation |> dplyr::left_join( module_sizes, by = c( "raw_module", "module_name" ) )
pb_logcpm <- as.matrix( step6a$pseudobulk_logCPM )
pb_meta <- as.data.frame( step6a$pseudobulk_metadata )
if ( !"pb_group" %in% colnames(pb_meta) || !"sample" %in% colnames(pb_meta) || !"branch" %in% colnames(pb_meta) ) {
  stop( "Step6A pseudobulk_metadata缺少pb_group/sample/branch。", call. = FALSE )
}
if (!all(pb_meta$pb_group %in% colnames(pb_logcpm))) {
  stop( "Step6A pseudobulk metadata与logCPM列不匹配。", call. = FALSE )
}
pb_meta <- pb_meta[ match( colnames(pb_logcpm), pb_meta$pb_group ), , drop = FALSE ]
module_pb_list <- lapply( unique(gene_modules$module_name), function(mod) {
    genes_mod <- gene_modules$gene_id[ gene_modules$module_name == mod ]
    genes_mod <- intersect( genes_mod, rownames(pb_logcpm) )
    if (length(genes_mod) < 2) {
      return(NULL)
    }
    mm <- pb_logcpm[ genes_mod, , drop = FALSE ]
    mm_z <- row_zscore(mm)
    mm_z <- mm_z[ apply(mm_z, 1, function(v) all(is.finite(v))), , drop = FALSE ]
    if (nrow(mm_z) < 2) {
      return(NULL)
    }
    tibble::tibble( module_name = mod, pb_group = colnames(mm_z), module_score = colMeans( mm_z, na.rm = TRUE ) )
  }
)
module_pb_scores <- dplyr::bind_rows( module_pb_list ) |> dplyr::left_join( tibble::as_tibble(pb_meta), by = "pb_group" )
module_effects <- module_pb_scores |> dplyr::select( module_name, sample, branch, module_score ) |> tidyr::pivot_wider( names_from = branch, values_from = module_score ) |> dplyr::mutate( Arm1_minus_Arm2 = Arm_1 - Arm_2 )
module_robustness <- module_effects |> dplyr::group_by( module_name ) |> dplyr::summarise( n_samples = dplyr::n(), mean_Arm1_minus_Arm2 = mean( Arm1_minus_Arm2, na.rm = TRUE ), min_abs_effect = min( abs(Arm1_minus_Arm2), na.rm = TRUE ), same_direction_all_samples = all(Arm1_minus_Arm2 > 0) || all(Arm1_minus_Arm2 < 0), robust_terminal_branch_effect = same_direction_all_samples && min_abs_effect >= MIN_ABS_MODULE_EFFECT_PER_SAMPLE, .groups = "drop" )
module_annotation <- module_annotation |> dplyr::left_join( module_robustness, by = "module_name" )
module_source_counts <- gene_modules |> dplyr::count( module_name, evidence_source, name = "n_genes" )
top_genes_by_module <- gene_modules |> dplyr::group_by( module_name ) |> dplyr::arrange( dplyr::desc(module_membership_cor), progressive_shape_q_value, cell_BH_FDR, .by_group = TRUE ) |> dplyr::slice_head( n = TOP_REPORT_GENES_PER_MODULE ) |> dplyr::ungroup()
heat_genes <- gene_modules |> dplyr::group_by( module_name ) |> dplyr::arrange( dplyr::desc(module_membership_cor), .by_group = TRUE ) |> dplyr::slice_head( n = TOP_HEAT_GENES_PER_MODULE ) |> dplyr::ungroup()
heat_gene_order <- heat_genes |> dplyr::arrange( module_name, module_membership_cor ) |> dplyr::pull(gene_id)
heat_long <- as.data.frame( profile_z[ heat_gene_order, , drop = FALSE ] ) |> tibble::rownames_to_column( "gene_id" ) |> tidyr::pivot_longer( cols = -gene_id, names_to = "trajectory_bin", values_to = "row_z" ) |> dplyr::left_join( heat_genes |> dplyr::select( gene_id, module_name ), by = "gene_id" ) |> tidyr::separate( trajectory_bin, into = c( "segment", "bin" ), sep = "_", remove = FALSE ) |> dplyr::mutate( bin = as.integer(bin), segment = factor( segment, levels = c( "Pre", "Arm1", "Arm2" ) ), gene_id = factor( gene_id, levels = heat_gene_order ) )
p1 <- ggplot2::ggplot( heat_long, ggplot2::aes( x = bin, y = gene_id, fill = row_z ) ) + ggplot2::geom_tile() + ggplot2::facet_grid( module_name ~ segment, scales = "free", space = "free" ) + ggplot2::scale_fill_gradient2( midpoint = 0, name = "Row z-score" ) + ggplot2::labs( title = "Sample-balanced dynamic branch modules", subtitle = paste0( "Union of Step5C-R progressive genes and Step6A high-confidence markers | K = ", selected_k ), x = "Trajectory bin", y = "Gene" ) + ggplot2::theme_classic( base_size = 8 ) + ggplot2::theme( axis.text.y = ggplot2::element_text( size = 5.2 ), strip.text.y = ggplot2::element_text( angle = 0 ) )
grDevices::pdf( OUT_FIG1, width = 13, height = max( 11, 2.0 * selected_k + 5 ), useDingbats = FALSE )
print(p1)
grDevices::dev.off()
centroid_long <- dplyr::bind_rows( lapply( names(module_centroids), function(mm) {
      tibble::tibble( raw_module = as.integer(mm), trajectory_bin = names(module_centroids[[mm]]), centroid_z = as.numeric(module_centroids[[mm]]) )
    }
  )
) |> dplyr::left_join(
    module_annotation |> dplyr::select( raw_module, module_name ), by = "raw_module"
  ) |> tidyr::separate(
    trajectory_bin, into = c( "segment", "bin" ), sep = "_", remove = FALSE
  ) |> dplyr::mutate(
    bin = as.integer(bin)
  )
pre_centroid <- centroid_long |> dplyr::filter( segment == "Pre" ) |> dplyr::mutate( x = -1 + (bin - 1) / max(N_PRE_BINS - 1, 1), lineage = "Shared_prebranch" )
arm1_centroid <- centroid_long |> dplyr::filter( segment == "Arm1" ) |> dplyr::mutate( x = (bin - 1) / max(N_POST_BINS - 1, 1), lineage = "Arm_1" )
arm2_centroid <- centroid_long |> dplyr::filter( segment == "Arm2" ) |> dplyr::mutate( x = (bin - 1) / max(N_POST_BINS - 1, 1), lineage = "Arm_2" )
pre_endpoint <- pre_centroid |> dplyr::group_by( module_name ) |> dplyr::slice_max( order_by = x, n = 1, with_ties = FALSE ) |> dplyr::ungroup()
centroid_plot_tbl <- dplyr::bind_rows( pre_centroid, pre_endpoint |> dplyr::mutate( lineage = "Arm_1" ), pre_endpoint |> dplyr::mutate( lineage = "Arm_2" ), arm1_centroid, arm2_centroid )
p2 <- ggplot2::ggplot( centroid_plot_tbl, ggplot2::aes( x = x, y = centroid_z, group = lineage, linetype = lineage ) ) + ggplot2::geom_line( linewidth = 0.9 ) + ggplot2::geom_point( size = 1.4 ) + ggplot2::geom_vline( xintercept = 0, linetype = 2, linewidth = 0.4 ) + ggplot2::facet_wrap( ~module_name, scales = "free_y", ncol = 2 ) + ggplot2::labs( title = "Dynamic module centroid trajectories", subtitle = "Negative x = shared prebranch; positive x = progress away from Y_29", x = "Branchpoint-centered progression", y = "Module centroid z-score", linetype = "Lineage" ) + ggplot2::theme_classic( base_size = 10 )
grDevices::pdf( OUT_FIG2, width = 11, height = max( 7, 2.5 * ceiling(selected_k / 2) ), useDingbats = FALSE )
print(p2)
grDevices::dev.off()
p3 <- ggplot2::ggplot( module_pb_scores, ggplot2::aes( x = branch, y = module_score, group = sample, shape = sample ) ) + ggplot2::geom_line( linewidth = 0.65, alpha = 0.8 ) + ggplot2::geom_point( size = 2.7 ) + ggplot2::facet_wrap( ~module_name, scales = "free_y", ncol = 2 ) + ggplot2::labs( title = "7 dpi pseudobulk dynamic-module scores", subtitle = "Each line connects Arm_1 and Arm_2 within the same biological sample", x = "Terminal branch", y = "Gene-centered module score", shape = "Biological sample" ) + ggplot2::theme_classic( base_size = 10 )
grDevices::pdf( OUT_FIG3, width = 11, height = max( 7, 2.5 * ceiling(selected_k / 2) ), useDingbats = FALSE )
print(p3)
grDevices::dev.off()
p4a <- ggplot2::ggplot( module_source_counts, ggplot2::aes( x = module_name, y = n_genes, fill = evidence_source ) ) + ggplot2::geom_col() + ggplot2::labs( title = "Module evidence composition", subtitle = "Step5C-R dynamic evidence vs Step6A terminal evidence", x = NULL, y = "Genes", fill = "Evidence source" ) + ggplot2::theme_classic( base_size = 10 ) + ggplot2::theme( axis.text.x = ggplot2::element_text( angle = 30, hjust = 1 ) )
p4b <- ggplot2::ggplot( module_effects, ggplot2::aes( x = module_name, y = Arm1_minus_Arm2, group = sample, shape = sample ) ) + ggplot2::geom_hline( yintercept = 0, linetype = 2 ) + ggplot2::geom_line( linewidth = 0.65, alpha = 0.8 ) + ggplot2::geom_point( size = 2.7 ) + ggplot2::labs( title = "7 dpi module branch effect by biological sample", subtitle = "Positive = Arm_1 higher; negative = Arm_2 higher", x = NULL, y = "Arm_1 - Arm_2 module score", shape = "Biological sample" ) + ggplot2::theme_classic( base_size = 10 ) + ggplot2::theme( axis.text.x = ggplot2::element_text( angle = 30, hjust = 1 ) )
grDevices::pdf( OUT_FIG4, width = 13, height = 6.5, useDingbats = FALSE )
print(p4a | p4b)
grDevices::dev.off()
utils::write.csv( gene_modules, file.path( TABDIR, "01_dynamic_module_gene_assignments.csv" ), row.names = FALSE )
utils::write.csv( module_annotation, file.path( TABDIR, "02_dynamic_module_summary.csv" ), row.names = FALSE )
utils::write.csv( top_genes_by_module, file.path( TABDIR, "03_top_genes_by_module.csv" ), row.names = FALSE )
utils::write.csv( module_effects, file.path( TABDIR, "04_module_7dpi_replicate_effects.csv" ), row.names = FALSE )
profile_table <- as.data.frame( profile_mat ) |> tibble::rownames_to_column( "gene_id" )
profile_z_table <- as.data.frame( profile_z ) |> tibble::rownames_to_column( "gene_id" )
parameter_tbl <- tibble::tibble( parameter = c( "Step5C-R progressive genes", "Step5C-R strong progressive genes", "Step6A high-confidence markers", "Union candidate genes", "Genes after detection filter", "Genes entering clustering", "Prebranch bins", "Postbranch bins per arm", "Selected K", "Selected average silhouette", "Minimum module size threshold", "Minimum 7 dpi per-sample module effect" ), value = c( length(step5_progressive), length(step5_strong), length(step6a_highconf), length(candidate_genes), length(keep_genes), nrow(profile_z), N_PRE_BINS, N_POST_BINS, selected_k, k_search_tbl$average_silhouette[ k_search_tbl$K == selected_k ], MIN_MODULE_SIZE, MIN_ABS_MODULE_EFFECT_PER_SAMPLE ) )
wb <- openxlsx::createWorkbook()
safe_sheet(wb, "Parameters", parameter_tbl)
safe_sheet(wb, "K_selection", k_search_tbl)
safe_sheet(wb, "Module_summary", module_annotation)
safe_sheet(wb, "Gene_assignments", gene_modules)
safe_sheet(wb, "Top_genes", top_genes_by_module)
safe_sheet(wb, "Sample_segment_audit", group_size_audit)
safe_sheet(wb, "Sample_bin_counts", sample_bin_summary)
safe_sheet(wb, "Profile_means", profile_table)
safe_sheet(wb, "Profile_zscores", profile_z_table)
safe_sheet(wb, "PB_module_scores", module_pb_scores)
safe_sheet(wb, "Module_effects", module_effects)
safe_sheet(wb, "Source_counts", module_source_counts)
openxlsx::saveWorkbook( wb, OUT_XLSX, overwrite = TRUE )
result_object <- list( candidate_genes = candidate_genes, genes_entering_clustering = rownames(profile_z), selected_k = selected_k, k_search = k_search_tbl, sample_segment_audit = group_size_audit, sample_bin_summary = sample_bin_summary, sample_balanced_profile = profile_mat, sample_balanced_profile_z = profile_z, module_annotation = module_annotation, gene_modules = gene_modules, top_genes_by_module = top_genes_by_module, module_pseudobulk_scores = module_pb_scores, module_7dpi_effects = module_effects, module_robustness = module_robustness, module_source_counts = module_source_counts )
saveRDS( result_object, OUT_RDS, compress = TRUE )
module_report_tbl <- module_annotation |> dplyr::select( module_name, n_genes, temporal_label, pre_score, early_score, intermediate_score, arm1_late_score, arm2_late_score, arm1_minus_arm2_late, mean_Arm1_minus_Arm2, min_abs_effect, same_direction_all_samples, robust_terminal_branch_effect ) |> dplyr::arrange( module_name )
top_gene_report <- top_genes_by_module |> dplyr::select( module_name, gene_id, module_membership_cor, evidence_source, progressive_shape_q_value, cell_avg_log2FC, cell_BH_FDR, final_direction, strict_pseudobulk_supported ) |> dplyr::arrange( module_name, dplyr::desc(module_membership_cor) )
writeLines( capture.output(sessionInfo()), OUT_SESS, useBytes = TRUE )
print( module_annotation |> dplyr::select( module_name, n_genes, arm1_minus_arm2_late, mean_Arm1_minus_Arm2, same_direction_all_samples, robust_terminal_branch_effect ) )
