options(stringsAsFactors = FALSE)
set.seed(20260731)
PROJECT_ROOT <- Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
INPUT_RDS <- file.path( PROJECT_ROOT, "04_results", "step5A_GSE162610_balanced_root_trajectory_input", "GSE162610_step5A_balanced_root_trajectory_input.rds" )
OUTDIR <- file.path( PROJECT_ROOT, "04_results", "step5AR_GSE162610_consensus_homeostatic_root_refinement" )
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
OUT_FIG1 <- file.path(OUTDIR, "Fig5AR_01_cluster0_consensus_space_audit.pdf")
OUT_FIG2 <- file.path(OUTDIR, "Fig5AR_02_consensus_core_selection.pdf")
OUT_FIG3 <- file.path(OUTDIR, "Fig5AR_03_final_consensus_root.pdf")
OUT_XLSX <- file.path(OUTDIR, "Table5AR_consensus_root_audit.xlsx")
OUT_RDS  <- file.path(OUTDIR, "GSE162610_step5AR_consensus_root_trajectory_input.rds")
OUT_REP  <- file.path(OUTDIR, "Step5AR_report.txt")
OUT_SESS <- file.path(OUTDIR, "Step5AR_sessionInfo.txt")
TIME_COL       <- "timepoint"
SAMPLE_COL     <- "sample_id"
CLUSTER_COL    <- "step3A_cluster_r04"
UMAP_REDUCTION <- "umap.microglia.step3A"
HARMONY_REDUCTION <- "harmony.step3A"
HARMONY_DIMS <- 1:30
EXPECTED_UNINJ_SAMPLES <- c( "uninj_sample1", "uninj_sample2", "uninj_sample3" )
ROOT_CLUSTER <- "0"
ROOT_SCORE_PCT_MIN       <- 0.40
CENTRALITY_PCT_MIN       <- 0.40
CROSS_SAMPLE_PCT_MIN     <- 0.25
TARGET_ROOT_PER_SAMPLE <- 60
MIN_ROOT_PER_SAMPLE    <- 20
W_ROOT    <- 0.50
W_CENTRAL <- 0.30
W_CROSS   <- 0.20
HOMEOSTATIC_GENES <- c( "P2ry12", "P2ry13", "Tmem119", "Sall1", "Fcrls", "Hexb", "Cx3cr1", "Gpr34", "Siglech", "Olfml3", "Slc2a5", "Adgrg1" )
ACTIVATED_DAM_GENES <- c( "Apoe", "Lpl", "Spp1", "Gpnmb", "Lgals3", "Fabp5", "Cst7", "Itgax", "Axl", "Clec7a", "Ccl2", "Ccl3", "Ccl4", "Il1b", "Tnf", "Nfkbia" )
required_pkgs <- c( "Seurat", "SeuratObject", "Matrix", "dplyr", "tidyr", "tibble", "ggplot2", "patchwork", "openxlsx", "RANN" )
missing_pkgs <- required_pkgs[ !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE) ]
if (length(missing_pkgs) > 0) {
  stop(paste0("Missing R packages: ", paste(missing_pkgs, collapse = ", ")), call. = FALSE)
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
})
SAMPLE_COLORS <- c( "uninj_sample1" = "#E76F8A", "uninj_sample2" = "#76B947", "uninj_sample3" = "#27A9D8" )
TIME_COLORS <- c( "Uninjured" = "#66C2A5", "1 dpi"     = "#FC8D62", "3 dpi"     = "#8DA0CB", "7 dpi"     = "#E78AC3" )
ROOT_SCORE_COLORS <- c( "#313695", "#74ADD1", "#FFFFBF", "#FDAE61", "#A50026" )
CENTRALITY_COLORS <- c( "#FFF5F0", "#FB6A4A", "#99000D" )
CROSS_COLORS <- c( "#F7FCF5", "#74C476", "#00441B" )
ROOT_COL <- "#D73027"
BG_COL <- "#D9D9D9"
capture_table <- function(x) {
  paste( capture.output( print(tibble::as_tibble(x), n = Inf, width = Inf) ), collapse = "\n" )
}
get_data_matrix <- function(obj) {
  assay <- DefaultAssay(obj)
  x <- tryCatch( Seurat::GetAssayData(obj, assay = assay, layer = "data"), error = function(e) NULL )
  if (!is.null(x)) return(x)
  Seurat::GetAssayData(obj, assay = assay, slot = "data")
}
calc_gene_z_program <- function(mat, genes) {
  genes <- intersect(genes, rownames(mat))
  if (length(genes) == 0) {
    stop("Requested program has zero available genes.", call. = FALSE)
  }
  xx <- as.matrix(mat[genes, , drop = FALSE])
  mu <- rowMeans(xx, na.rm = TRUE)
  ss <- apply(xx, 1, stats::sd, na.rm = TRUE)
  ss[!is.finite(ss) | ss == 0] <- 1
  zz <- sweep(xx, 1, mu, "-")
  zz <- sweep(zz, 1, ss, "/")
  list( score = colMeans(zz, na.rm = TRUE), genes = genes )
}
euclidean_to_vector <- function(mat, vec) {
  sqrt(rowSums((mat - matrix( vec, nrow = nrow(mat), ncol = ncol(mat), byrow = TRUE )) ^ 2))
}
if (!file.exists(INPUT_RDS)) {
  stop("找不到Step5A正式trajectory input RDS：\n", INPUT_RDS)
}
obj <- readRDS(INPUT_RDS)
if (!inherits(obj, "Seurat")) {
  stop("Input is not a Seurat object.", call. = FALSE)
}
DefaultAssay(obj) <- "RNA"
required_meta <- c( TIME_COL, SAMPLE_COL, CLUSTER_COL, "transition_status", "shared_transition_3dpi_highconf" )
missing_meta <- required_meta[ !required_meta %in% colnames(obj[[]]) ]
if (length(missing_meta) > 0) {
  stop( "缺少metadata：", paste(missing_meta, collapse = ", "), call. = FALSE )
}
required_reductions <- c( UMAP_REDUCTION, HARMONY_REDUCTION )
missing_reductions <- setdiff( required_reductions, Seurat::Reductions(obj) )
if (length(missing_reductions) > 0) {
  stop( paste0( "缺少正式reduction：", paste(missing_reductions, collapse = ", "), "\n重新运行Harmony或UMAP。" ), call. = FALSE )
}
time_vec <- as.character(obj[[TIME_COL, drop = TRUE]])
sample_vec <- as.character(obj[[SAMPLE_COL, drop = TRUE]])
cluster_vec <- as.character(obj[[CLUSTER_COL, drop = TRUE]])
candidate_cells <- colnames(obj)[ time_vec == "Uninjured" & cluster_vec == ROOT_CLUSTER ]
candidate_tbl <- tibble::tibble( cell = candidate_cells, sample = sample_vec[match(candidate_cells, colnames(obj))], cluster = cluster_vec[match(candidate_cells, colnames(obj))] )
candidate_counts <- candidate_tbl |> dplyr::count(sample, name = "n_candidates") |> tidyr::complete( sample = EXPECTED_UNINJ_SAMPLES, fill = list(n_candidates = 0) )
print(candidate_counts)
if (any(candidate_counts$n_candidates < MIN_ROOT_PER_SAMPLE)) {
  stop( "至少一个Uninjured sample在cluster 0中的候选细胞 < ", MIN_ROOT_PER_SAMPLE, "，不能继续。", call. = FALSE )
}
expr <- get_data_matrix(obj)
if (!identical(colnames(expr), colnames(obj))) {
  expr <- expr[, colnames(obj), drop = FALSE]
}
expr_candidate <- expr[, candidate_cells, drop = FALSE]
homeo_obj <- calc_gene_z_program( expr_candidate, HOMEOSTATIC_GENES )
dam_obj <- calc_gene_z_program( expr_candidate, ACTIVATED_DAM_GENES )
candidate_tbl$homeostatic_program_score <- as.numeric(homeo_obj$score)
candidate_tbl$activated_DAM_penalty_score <- as.numeric(dam_obj$score)
candidate_tbl$root_score <- candidate_tbl$homeostatic_program_score - candidate_tbl$activated_DAM_penalty_score
harmony <- Seurat::Embeddings( obj, reduction = HARMONY_REDUCTION )
if (ncol(harmony) < max(HARMONY_DIMS)) {
  stop( paste0( HARMONY_REDUCTION, " only has ", ncol(harmony), " dimensions; expected >=30." ), call. = FALSE )
}
harm_candidate <- harmony[ candidate_cells, HARMONY_DIMS, drop = FALSE ]
harm_mu <- colMeans(harm_candidate)
harm_sd <- apply(harm_candidate, 2, stats::sd)
harm_sd[!is.finite(harm_sd) | harm_sd == 0] <- 1
harm_z <- sweep(harm_candidate, 2, harm_mu, "-")
harm_z <- sweep(harm_z, 2, harm_sd, "/")
sample_centroids <- lapply( EXPECTED_UNINJ_SAMPLES, function(smp) {
    cells <- candidate_tbl$cell[candidate_tbl$sample == smp]
    colMeans(harm_z[cells, , drop = FALSE])
  }
)
names(sample_centroids) <- EXPECTED_UNINJ_SAMPLES
consensus_centroid <- Reduce("+", sample_centroids) / length(sample_centroids)
candidate_tbl$distance_to_consensus_centroid <- euclidean_to_vector( harm_z[candidate_tbl$cell, , drop = FALSE], consensus_centroid )
cross_metric_list <- lapply( EXPECTED_UNINJ_SAMPLES, function(target_sample) {
    target_cells <- candidate_tbl$cell[ candidate_tbl$sample == target_sample ]
    other_samples <- setdiff( EXPECTED_UNINJ_SAMPLES, target_sample )
    target_mat <- harm_z[target_cells, , drop = FALSE]
    nearest_by_other <- lapply( other_samples, function(other_sample) {
        other_cells <- candidate_tbl$cell[ candidate_tbl$sample == other_sample ]
        other_mat <- harm_z[ other_cells, , drop = FALSE ]
        nn <- RANN::nn2( data = other_mat, query = target_mat, k = 1 )
        as.numeric(nn$nn.dists[, 1])
      }
    )
    nearest_mat <- do.call( cbind, nearest_by_other )
    tibble::tibble( cell = target_cells, nearest_other_sample_mean = rowMeans(nearest_mat), nearest_other_sample_max = apply(nearest_mat, 1, max) )
  }
) |> dplyr::bind_rows()
candidate_tbl <- candidate_tbl |> dplyr::left_join( cross_metric_list, by = "cell" )
candidate_tbl <- candidate_tbl |> dplyr::group_by(sample) |> dplyr::mutate( root_score_percentile = dplyr::percent_rank(root_score), consensus_centrality_percentile = 1 - dplyr::percent_rank( distance_to_consensus_centroid ), cross_sample_proximity_percentile = 1 - dplyr::percent_rank( nearest_other_sample_max ) ) |> dplyr::ungroup() |> dplyr::mutate( pass_root_score = root_score_percentile >= ROOT_SCORE_PCT_MIN, pass_centrality = consensus_centrality_percentile >= CENTRALITY_PCT_MIN, pass_cross_sample = cross_sample_proximity_percentile >= CROSS_SAMPLE_PCT_MIN, highconf_consensus = pass_root_score & pass_centrality & pass_cross_sample, consensus_priority = W_ROOT * root_score_percentile + W_CENTRAL * consensus_centrality_percentile + W_CROSS * cross_sample_proximity_percentile )
highconf_counts <- candidate_tbl |> dplyr::filter(highconf_consensus) |> dplyr::count(sample, name = "n_highconf") |> tidyr::complete( sample = EXPECTED_UNINJ_SAMPLES, fill = list(n_highconf = 0) )
print(highconf_counts)
min_highconf <- min(highconf_counts$n_highconf)
if (min_highconf < MIN_ROOT_PER_SAMPLE) {
  stop( paste0( "最少样本仅有 ", min_highconf, " 个高置信Harmony-space root candidates (<", MIN_ROOT_PER_SAMPLE, ")。\n", "停止运行并检查 root-candidate 阈值。" ), call. = FALSE )
}
final_n_per_sample <- min( TARGET_ROOT_PER_SAMPLE, min_highconf )
final_root_tbl <- candidate_tbl |> dplyr::filter(highconf_consensus) |> dplyr::group_by(sample) |> dplyr::arrange( dplyr::desc(consensus_priority), dplyr::desc(root_score), .by_group = TRUE ) |> dplyr::mutate( within_sample_rank = dplyr::row_number() ) |> dplyr::filter( within_sample_rank <= final_n_per_sample ) |> dplyr::ungroup()
final_root_cells <- final_root_tbl$cell
final_counts <- final_root_tbl |> dplyr::count(sample, name = "n_root")
if ( nrow(final_counts) != length(EXPECTED_UNINJ_SAMPLES) || length(unique(final_counts$n_root)) != 1 ) {
  stop( "Internal error: final root is not exactly sample-balanced.", call. = FALSE )
}
obj$step5AR_v2_root_candidate <- FALSE
obj$step5AR_v2_highconf_consensus <- FALSE
obj$step5AR_v2_consensus_root <- FALSE
obj$step5AR_v2_root_score <- NA_real_
obj$step5AR_v2_root_score_percentile <- NA_real_
obj$step5AR_v2_distance_to_consensus <- NA_real_
obj$step5AR_v2_centrality_percentile <- NA_real_
obj$step5AR_v2_cross_sample_distance <- NA_real_
obj$step5AR_v2_cross_sample_percentile <- NA_real_
obj$step5AR_v2_consensus_priority <- NA_real_
idx <- match(candidate_tbl$cell, colnames(obj))
obj$step5AR_v2_root_candidate[idx] <- TRUE
obj$step5AR_v2_highconf_consensus[ match( candidate_tbl$cell[candidate_tbl$highconf_consensus], colnames(obj) ) ] <- TRUE
obj$step5AR_v2_consensus_root[ match(final_root_cells, colnames(obj)) ] <- TRUE
obj$step5AR_v2_root_score[idx] <- candidate_tbl$root_score
obj$step5AR_v2_root_score_percentile[idx] <- candidate_tbl$root_score_percentile
obj$step5AR_v2_distance_to_consensus[idx] <- candidate_tbl$distance_to_consensus_centroid
obj$step5AR_v2_centrality_percentile[idx] <- candidate_tbl$consensus_centrality_percentile
obj$step5AR_v2_cross_sample_distance[idx] <- candidate_tbl$nearest_other_sample_max
obj$step5AR_v2_cross_sample_percentile[idx] <- candidate_tbl$cross_sample_proximity_percentile
obj$step5AR_v2_consensus_priority[idx] <- candidate_tbl$consensus_priority
centroid_table <- function(cells, label) {
  tmp <- tibble::tibble( cell = cells, sample = sample_vec[match(cells, colnames(obj))] )
  cent <- lapply( EXPECTED_UNINJ_SAMPLES, function(smp) {
      cc <- tmp$cell[tmp$sample == smp]
      colMeans(harm_z[cc, , drop = FALSE])
    }
  )
  names(cent) <- EXPECTED_UNINJ_SAMPLES
  pairs <- combn( EXPECTED_UNINJ_SAMPLES, 2, simplify = FALSE )
  lapply( pairs, function(pp) {
      tibble::tibble( stage = label, sample_A = pp[1], sample_B = pp[2], centroid_distance = sqrt( sum( (cent[[pp[1]]] - cent[[pp[2]]]) ^ 2 ) ) )
    }
  ) |> dplyr::bind_rows()
}
candidate_centroid_dist <- centroid_table( candidate_cells, "All_cluster0_candidates" )
final_centroid_dist <- centroid_table( final_root_cells, "Final_consensus_root" )
centroid_distance_compare <- dplyr::bind_rows( candidate_centroid_dist, final_centroid_dist )
mean_candidate_centroid_distance <- mean(candidate_centroid_dist$centroid_distance)
mean_final_centroid_distance <- mean(final_centroid_dist$centroid_distance)
centroid_distance_reduction <- 1 - mean_final_centroid_distance / mean_candidate_centroid_distance
final_summary <- final_root_tbl |> dplyr::group_by(sample) |> dplyr::summarise( n_root = dplyr::n(), root_score_mean = mean(root_score, na.rm = TRUE), root_score_median = median(root_score, na.rm = TRUE), centrality_pct_mean = mean( consensus_centrality_percentile, na.rm = TRUE ), cross_sample_pct_mean = mean( cross_sample_proximity_percentile, na.rm = TRUE ), consensus_priority_mean = mean( consensus_priority, na.rm = TRUE ), .groups = "drop" )
umap <- Seurat::Embeddings( obj, reduction = UMAP_REDUCTION )
plot_tbl <- tibble::tibble( cell = rownames(umap), UMAP_1 = umap[, 1], UMAP_2 = umap[, 2], timepoint = time_vec[ match(rownames(umap), colnames(obj)) ], sample = sample_vec[ match(rownames(umap), colnames(obj)) ], cluster = cluster_vec[ match(rownames(umap), colnames(obj)) ] )
candidate_plot <- plot_tbl |> dplyr::filter(cell %in% candidate_cells) |> dplyr::left_join( candidate_tbl, by = c("cell", "sample", "cluster") )
trajectory_plot <- plot_tbl |> dplyr::mutate( final_root = cell %in% final_root_cells )
theme5 <- ggplot2::theme_classic(base_size = 11) + ggplot2::theme( axis.text = ggplot2::element_text(color = "#222222"), axis.title = ggplot2::element_text(color = "#222222"), plot.title = ggplot2::element_text(face = "bold", color = "#111111"), plot.subtitle = ggplot2::element_text(color = "#333333") )
p1a <- ggplot2::ggplot( candidate_plot, ggplot2::aes( UMAP_1, UMAP_2, color = factor(sample) ) ) + ggplot2::geom_point(size = 0.65, alpha = 0.75) + ggplot2::scale_color_manual( values = SAMPLE_COLORS, name = "Uninjured sample" ) + ggplot2::labs( title = "Shared root state: Uninjured Step3A cluster 0", subtitle = "All three biological samples are represented" ) + theme5
p1b <- ggplot2::ggplot( candidate_plot, ggplot2::aes( UMAP_1, UMAP_2, color = root_score ) ) + ggplot2::geom_point(size = 0.65, alpha = 0.80) + ggplot2::scale_color_gradientn( colors = ROOT_SCORE_COLORS, name = "Root score" ) + ggplot2::labs( title = "Biological homeostatic-root score" ) + theme5
p1c <- ggplot2::ggplot( candidate_plot, ggplot2::aes( UMAP_1, UMAP_2, color = consensus_centrality_percentile ) ) + ggplot2::geom_point(size = 0.65, alpha = 0.80) + ggplot2::scale_color_gradientn( colors = CENTRALITY_COLORS, limits = c(0, 1), name = "Consensus\ncentrality" ) + ggplot2::labs( title = "Closeness to equal-weight Harmony consensus centroid", subtitle = "Calculated in harmony.step3A dimensions 1:30" ) + theme5
p1d <- ggplot2::ggplot( candidate_plot, ggplot2::aes( UMAP_1, UMAP_2, color = cross_sample_proximity_percentile ) ) + ggplot2::geom_point(size = 0.65, alpha = 0.80) + ggplot2::scale_color_gradientn( colors = CROSS_COLORS, limits = c(0, 1), name = "Cross-sample\nproximity" ) + ggplot2::labs( title = "Proximity to cells from the other biological samples", subtitle = "Calculated in standardized Harmony space" ) + theme5
grDevices::pdf( OUT_FIG1, width = 14, height = 10, useDingbats = FALSE )
print( (p1a | p1b) / (p1c | p1d) + patchwork::plot_annotation( title = "Step 5A-R v2 | Harmony-space consensus-root audit" ) )
grDevices::dev.off()
candidate_plot$core_status <- ifelse( candidate_plot$highconf_consensus, "High-confidence consensus", "Not selected" )
p2a <- ggplot2::ggplot( candidate_plot, ggplot2::aes(UMAP_1, UMAP_2) ) + ggplot2::geom_point( data = candidate_plot |> dplyr::filter(core_status == "Not selected"), color = BG_COL, size = 0.55, alpha = 0.35 ) + ggplot2::geom_point( data = candidate_plot |> dplyr::filter(core_status == "High-confidence consensus"), color = "#2C7BB6", size = 0.80, alpha = 0.90 ) + ggplot2::labs( title = "High-confidence consensus core", subtitle = "Biological root quality + Harmony centrality + cross-sample proximity" ) + theme5
priority_long <- candidate_tbl |> dplyr::select( cell, sample, root_score_percentile, consensus_centrality_percentile, cross_sample_proximity_percentile, consensus_priority ) |> tidyr::pivot_longer( cols = c( root_score_percentile, consensus_centrality_percentile, cross_sample_proximity_percentile, consensus_priority ), names_to = "metric", values_to = "value" ) |> dplyr::mutate( metric = dplyr::recode( metric, root_score_percentile = "Root-score percentile", consensus_centrality_percentile = "Consensus centrality", cross_sample_proximity_percentile = "Cross-sample proximity", consensus_priority = "Composite consensus priority" ) )
p2b <- ggplot2::ggplot( priority_long, ggplot2::aes( x = factor(sample), y = value, fill = factor(sample) ) ) + ggplot2::geom_violin( trim = TRUE, scale = "width", alpha = 0.75, linewidth = 0.25 ) + ggplot2::geom_boxplot( width = 0.12, outlier.shape = NA, fill = "#FFFFFF", linewidth = 0.25 ) + ggplot2::facet_wrap( ~metric, ncol = 2, scales = "free_y" ) + ggplot2::scale_fill_manual( values = SAMPLE_COLORS, guide = "none" ) + ggplot2::labs( title = "Consensus metrics by biological sample", x = "Uninjured sample", y = "Percentile / priority" ) + theme5 + ggplot2::theme( axis.text.x = ggplot2::element_text(angle = 25, hjust = 1) )
highconf_plot <- highconf_counts |> dplyr::mutate( sample = factor( sample, levels = EXPECTED_UNINJ_SAMPLES ) )
p2c <- ggplot2::ggplot( highconf_plot, ggplot2::aes( x = sample, y = n_highconf, fill = sample ) ) + ggplot2::geom_col(width = 0.70) + ggplot2::geom_hline( yintercept = MIN_ROOT_PER_SAMPLE, linetype = 2, linewidth = 0.4, color = "#555555" ) + ggplot2::scale_fill_manual( values = SAMPLE_COLORS, guide = "none" ) + ggplot2::labs( title = "High-confidence cells before final balancing", subtitle = paste0( "Dashed line = minimum ", MIN_ROOT_PER_SAMPLE, "/sample" ), x = "Uninjured sample", y = "Cells" ) + theme5 + ggplot2::theme( axis.text.x = ggplot2::element_text(angle = 25, hjust = 1) )
p2d <- ggplot2::ggplot( final_root_tbl, ggplot2::aes( x = factor(sample), y = consensus_priority, fill = factor(sample) ) ) + ggplot2::geom_boxplot( width = 0.55, outlier.shape = NA, alpha = 0.80 ) + ggplot2::scale_fill_manual( values = SAMPLE_COLORS, guide = "none" ) + ggplot2::labs( title = "Final selected-root consensus quality", subtitle = paste0( final_n_per_sample, " cells/sample" ), x = "Uninjured sample", y = "Consensus priority" ) + theme5 + ggplot2::theme( axis.text.x = ggplot2::element_text(angle = 25, hjust = 1) )
grDevices::pdf( OUT_FIG2, width = 14, height = 10, useDingbats = FALSE )
print( (p2a | p2b) / (p2c | p2d) + patchwork::plot_annotation( title = "Step 5A-R v2 | Consensus-core filtering and sample balancing" ) )
grDevices::dev.off()
final_plot <- candidate_plot |> dplyr::mutate( final_root = cell %in% final_root_cells )
p3a <- ggplot2::ggplot( final_plot, ggplot2::aes(UMAP_1, UMAP_2) ) + ggplot2::geom_point( color = BG_COL, size = 0.55, alpha = 0.30 ) + ggplot2::geom_point( data = final_plot |> dplyr::filter(final_root), color = ROOT_COL, size = 0.95, alpha = 0.95 ) + ggplot2::labs( title = "Final consensus root within shared cluster 0", subtitle = paste0( final_n_per_sample, " cells/sample; red = future Monocle root cells" ) ) + theme5
p3b <- ggplot2::ggplot( final_plot |> dplyr::filter(final_root), ggplot2::aes( UMAP_1, UMAP_2, color = factor(sample) ) ) + ggplot2::geom_point( size = 1.00, alpha = 0.85 ) + ggplot2::scale_color_manual( values = SAMPLE_COLORS, name = "Root sample" ) + ggplot2::labs( title = "Final root cells by biological sample" ) + theme5
p3c <- ggplot2::ggplot( trajectory_plot, ggplot2::aes( UMAP_1, UMAP_2, color = factor(timepoint) ) ) + ggplot2::geom_point( size = 0.35, alpha = 0.42 ) + ggplot2::scale_color_manual( values = TIME_COLORS, name = "True time" ) + ggplot2::geom_point( data = trajectory_plot |> dplyr::filter(final_root), inherit.aes = FALSE, ggplot2::aes(UMAP_1, UMAP_2), color = ROOT_COL, size = 0.95, alpha = 0.95 ) + ggplot2::labs( title = "Final root on unchanged 17,863-cell trajectory input", subtitle = "Graph and pseudotime are NOT learned until Step 5B" ) + theme5
centroid_plot <- centroid_distance_compare |> dplyr::mutate( comparison = paste( sample_A, sample_B, sep = " vs " ) )
p3d <- ggplot2::ggplot( centroid_plot, ggplot2::aes( x = comparison, y = centroid_distance, fill = stage ) ) + ggplot2::geom_col( position = "dodge", width = 0.70 ) + ggplot2::scale_fill_manual( values = c( "All_cluster0_candidates" = "#9ECAE1", "Final_consensus_root" = "#D73027" ), name = NULL ) + ggplot2::labs( title = "Cross-sample centroid convergence in Harmony space", subtitle = paste0( "Mean distance reduction = ", round( 100 * centroid_distance_reduction, 1 ), "%" ), x = NULL, y = "Standardized Harmony centroid distance" ) + theme5 + ggplot2::theme( axis.text.x = ggplot2::element_text( angle = 25, hjust = 1 ) )
grDevices::pdf( OUT_FIG3, width = 14, height = 10, useDingbats = FALSE )
print( (p3a | p3b) / (p3c | p3d) + patchwork::plot_annotation( title = "Step 5A-R v2 | Final cross-sample consensus homeostatic root" ) )
grDevices::dev.off()
parameters <- tibble::tibble( parameter = c( "Root state", "Harmony reduction", "Harmony dimensions", "Root-score percentile minimum", "Consensus-centrality percentile minimum", "Cross-sample proximity percentile minimum", "Target root cells/sample", "Minimum root cells/sample", "Actual root cells/sample", "Final root total", "Mean candidate centroid distance", "Mean final-root centroid distance", "Centroid distance reduction" ), value = c( "Uninjured + Step3A cluster 0", HARMONY_REDUCTION, "1:30", ROOT_SCORE_PCT_MIN, CENTRALITY_PCT_MIN, CROSS_SAMPLE_PCT_MIN, TARGET_ROOT_PER_SAMPLE, MIN_ROOT_PER_SAMPLE, final_n_per_sample, length(final_root_cells), mean_candidate_centroid_distance, mean_final_centroid_distance, centroid_distance_reduction ) )
gene_audit <- tibble::tibble( program = c( rep("Homeostatic", length(HOMEOSTATIC_GENES)), rep("Activated_DAM_penalty", length(ACTIVATED_DAM_GENES)) ), gene = c( HOMEOSTATIC_GENES, ACTIVATED_DAM_GENES ) ) |> dplyr::mutate( available = gene %in% rownames(expr) )
wb <- openxlsx::createWorkbook()
add_sheet <- function(name, df) {
  openxlsx::addWorksheet(wb, name)
  openxlsx::writeDataTable( wb, sheet = name, x = df, tableStyle = "TableStyleMedium2" )
  openxlsx::freezePane( wb, sheet = name, firstRow = TRUE )
}
add_sheet("Parameters", parameters)
add_sheet("Candidate_counts", candidate_counts)
add_sheet("Gene_audit", gene_audit)
add_sheet("Candidate_metrics", candidate_tbl)
add_sheet("Highconf_counts", highconf_counts)
add_sheet("Final_root_cells", final_root_tbl)
add_sheet("Final_root_summary", final_summary)
add_sheet("Centroid_distances", centroid_distance_compare)
openxlsx::saveWorkbook( wb, OUT_XLSX, overwrite = TRUE )
saveRDS( obj, OUT_RDS, compress = TRUE )
writeLines( capture.output(sessionInfo()), OUT_SESS, useBytes = TRUE )
