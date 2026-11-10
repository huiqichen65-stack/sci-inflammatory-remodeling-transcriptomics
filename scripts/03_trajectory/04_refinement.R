options(stringsAsFactors = FALSE)
set.seed(20260731)
PROJECT_ROOT <- Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
INPUT_RDS <- file.path( PROJECT_ROOT, "04_results", "step5AR_GSE162610_consensus_homeostatic_root_refinement", "GSE162610_step5AR_consensus_root_trajectory_input.rds" )
OUTDIR <- file.path( PROJECT_ROOT, "04_results", "step5BR_GSE162610_injury_response_trajectory_refinement" )
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
OUT_FIG1 <- file.path(OUTDIR, "Fig5BR_01_injury_input_and_root.pdf")
OUT_FIG2 <- file.path(OUTDIR, "Fig5BR_02_Monocle3_injury_pseudotime.pdf")
OUT_FIG3 <- file.path(OUTDIR, "Fig5BR_03_branchpoint_and_arms.pdf")
OUT_FIG4 <- file.path(OUTDIR, "Fig5BR_04_posthoc_repair_scar_projection.pdf")
OUT_FIG5 <- file.path(OUTDIR, "Fig5BR_05_root_sensitivity.pdf")
OUT_XLSX <- file.path(OUTDIR, "Table5BR_injury_trajectory_results.xlsx")
OUT_CDS  <- file.path(OUTDIR, "GSE162610_step5BR_injury_monocle3_cds.rds")
OUT_RDS  <- file.path(OUTDIR, "GSE162610_step5BR_injury_trajectory_seurat.rds")
OUT_REP  <- file.path(OUTDIR, "Step5BR_report.txt")
OUT_SESS <- file.path(OUTDIR, "Step5BR_sessionInfo.txt")
TIME_COL       <- "timepoint"
SAMPLE_COL     <- "sample_id"
CLUSTER_COL    <- "step3A_cluster_r04"
UMAP_REDUCTION <- "umap.microglia.step3A"
HARMONY_REDUCTION <- "harmony.step3A"
HARMONY_DIMS <- 1:30
REPAIR_SCORE_COL <- "repair_candidate_score"
SCAR_SCORE_COL   <- "scar_candidate_score"
EARLY_ROOT_CLUSTERS <- c("2", "5", "10")
DAM_CLUSTERS <- c("0", "1", "4")
CENTRALITY_PCT_MIN <- 0.40
CROSS_SAMPLE_PCT_MIN <- 0.25
TARGET_ROOT_PER_SAMPLE <- 75
MIN_ROOT_PER_SAMPLE <- 30
W_CENTRAL <- 0.65
W_CROSS   <- 0.35
NUM_DIM <- 50
MONOCLE_K <- 30
LEARN_GRAPH_EUCLIDEAN_DISTANCE_RATIO <- 0.8
LEARN_GRAPH_CLOSE_LOOP <- FALSE
MIN_TOTAL_CELLS_PER_DOWNSTREAM_ARM <- 100
MIN_7DPI_CELLS_PER_DOWNSTREAM_ARM <- 100
BRANCHPOINT_GRAPH_RADIUS <- 2
required_pkgs <- c( "Seurat", "SeuratObject", "Matrix", "monocle3", "SingleCellExperiment", "SummarizedExperiment", "S4Vectors", "igraph", "RANN", "dplyr", "tidyr", "tibble", "ggplot2", "patchwork", "openxlsx" )
missing_pkgs <- required_pkgs[ !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE) ]
if (length(missing_pkgs) > 0) {
  stop(paste0("Missing R packages: ", paste(missing_pkgs, collapse = ", ")), call. = FALSE)
}
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(monocle3)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(S4Vectors)
  library(igraph)
  library(RANN)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
})
TIME_COLORS <- c( "1 dpi" = "#FC8D62", "3 dpi" = "#8DA0CB", "7 dpi" = "#E78AC3" )
CLUSTER_COLORS <- c( "0"  = "#4E79A7", "1"  = "#E15759", "2"  = "#59A14F", "4"  = "#B07AA1", "5"  = "#76B7B2", "10" = "#2F4B7C" )
ROOT_COL <- "#D73027"
BG_COL <- "#D9D9D9"
capture_table <- function(x) {
  paste( capture.output( print( tibble::as_tibble(x), n = Inf, width = Inf ) ), collapse = "\n" )
}
get_counts_matrix <- function(obj) {
  assay <- DefaultAssay(obj)
  x <- tryCatch( Seurat::GetAssayData( obj, assay = assay, layer = "counts" ), error = function(e) NULL )
  if (!is.null(x)) return(x)
  Seurat::GetAssayData( obj, assay = assay, slot = "counts" )
}
get_pseudotime_safe <- function(cds) {
  pt <- tryCatch( monocle3::pseudotime( cds, reduction_method = "UMAP" ), error = function(e) NULL )
  if (!is.null(pt)) return(pt)
  aux <- cds@principal_graph_aux[["UMAP"]]
  if ( !is.null(aux) && "pseudotime" %in% names(aux) ) {
    return(aux$pseudotime)
  }
  stop( "Could not extract Monocle3 pseudotime.", call. = FALSE )
}
get_clusters_safe <- function(cds) {
  x <- tryCatch( monocle3::clusters( cds, reduction_method = "UMAP" ), error = function(e) NULL )
  if (!is.null(x)) return(x)
  cds@clusters[["UMAP"]]$clusters
}
get_partitions_safe <- function(cds) {
  x <- tryCatch( monocle3::partitions( cds, reduction_method = "UMAP" ), error = function(e) NULL )
  if (!is.null(x)) return(x)
  cds@clusters[["UMAP"]]$partitions
}
get_closest_vertex <- function(cds, pg) {
  aux <- cds@principal_graph_aux[["UMAP"]]
  raw <- aux$pr_graph_cell_proj_closest_vertex
  if (is.null(raw)) {
    stop( "Monocle3 closest-principal-vertex mapping not found.", call. = FALSE )
  }
  if (is.matrix(raw) || is.data.frame(raw)) {
    vv <- raw[, 1]
    names(vv) <- rownames(raw)
  } else {
    vv <- raw
  }
  if (is.null(names(vv))) {
    names(vv) <- colnames(cds)
  }
  vv_chr <- as.character(vv)
  pg_names <- igraph::V(pg)$name
  if (all(vv_chr %in% pg_names)) {
    return(stats::setNames(vv_chr, names(vv)))
  }
  y_names <- paste0("Y_", vv_chr)
  if (all(y_names %in% pg_names)) {
    return(stats::setNames(y_names, names(vv)))
  }
  idx <- suppressWarnings(as.integer(vv_chr))
  if ( all(!is.na(idx)) && all(idx >= 1) && all(idx <= igraph::vcount(pg)) ) {
    return( stats::setNames( pg_names[idx], names(vv) ) )
  }
  stop( "Could not translate closest-vertex mapping to principal graph vertex names.", call. = FALSE )
}
euclidean_to_vector <- function(mat, vec) {
  sqrt( rowSums( ( mat - matrix( vec, nrow = nrow(mat), ncol = ncol(mat), byrow = TRUE ) ) ^ 2 ) )
}
cells_near_graph_vertex <- function( pg, closest_vertex, vertex_name, radius = 2) {
  dd <- igraph::distances( pg, v = vertex_name, to = igraph::V(pg), mode = "all" )
  near_vertices <- colnames(dd)[ as.numeric(dd[1, ]) <= radius ]
  names(closest_vertex)[ closest_vertex %in% near_vertices ]
}
add_sheet_safe <- function(wb, name, df) {
  openxlsx::addWorksheet( wb, sheetName = name )
  ddf <- tryCatch( as.data.frame(df), error = function(e) NULL )
  if ( is.null(ddf) || nrow(ddf) == 0 || ncol(ddf) == 0 ) {
    openxlsx::writeData( wb, sheet = name, x = data.frame( note = "No results generated under the predefined Step 5B-R criteria.", stringsAsFactors = FALSE ) )
  } else {
    openxlsx::writeDataTable( wb, sheet = name, x = ddf, tableStyle = "TableStyleMedium2" )
  }
  openxlsx::freezePane( wb, sheet = name, firstRow = TRUE )
}
if (!file.exists(INPUT_RDS)) {
  stop( paste0( "找不到Step5A-R正式输入：\n", INPUT_RDS ), call. = FALSE )
}
obj_full <- readRDS(INPUT_RDS)
if (!inherits(obj_full, "Seurat")) {
  stop( "Step5A-R input is not a Seurat object.", call. = FALSE )
}
DefaultAssay(obj_full) <- "RNA"
required_meta <- c( TIME_COL, SAMPLE_COL, CLUSTER_COL, "transition_status", "shared_transition_3dpi_highconf", REPAIR_SCORE_COL, SCAR_SCORE_COL )
missing_meta <- required_meta[ !required_meta %in% colnames(obj_full[[]]) ]
if (length(missing_meta) > 0) {
  stop( paste0( "缺少Step5B-R所需metadata：", paste(missing_meta, collapse = ", ") ), call. = FALSE )
}
required_reductions <- c( UMAP_REDUCTION, HARMONY_REDUCTION )
missing_reductions <- setdiff( required_reductions, Seurat::Reductions(obj_full) )
if (length(missing_reductions) > 0) {
  stop( paste0( "缺少正式reduction：", paste(missing_reductions, collapse = ", ") ), call. = FALSE )
}
time_full <- as.character( obj_full[[TIME_COL, drop = TRUE]] )
sample_full <- as.character( obj_full[[SAMPLE_COL, drop = TRUE]] )
cluster_full <- as.character( obj_full[[CLUSTER_COL, drop = TRUE]] )
is_1dpi_early <- ( time_full == "1 dpi" & cluster_full %in% EARLY_ROOT_CLUSTERS )
is_3dpi_transition_landscape <- ( time_full == "3 dpi" & ( cluster_full %in% DAM_CLUSTERS | as.logical( obj_full$shared_transition_3dpi_highconf ) ) )
is_7dpi_DAM_landscape <- ( time_full == "7 dpi" & cluster_full %in% DAM_CLUSTERS )
injury_keep <- ( is_1dpi_early | is_3dpi_transition_landscape | is_7dpi_DAM_landscape )
injury_cells <- colnames(obj_full)[ injury_keep ]
obj <- subset( obj_full, cells = injury_cells )
DefaultAssay(obj) <- "RNA"
time_vec <- as.character( obj[[TIME_COL, drop = TRUE]] )
sample_vec <- as.character( obj[[SAMPLE_COL, drop = TRUE]] )
cluster_vec <- as.character( obj[[CLUSTER_COL, drop = TRUE]] )
injury_input_summary <- obj[[]] |> dplyr::count( .data[[TIME_COL]], .data[[CLUSTER_COL]], name = "n_cells" )
colnames(injury_input_summary)[ colnames(injury_input_summary) == TIME_COL ] <- "timepoint"
colnames(injury_input_summary)[ colnames(injury_input_summary) == CLUSTER_COL ] <- "cluster"
time_counts <- tibble::tibble( timepoint = time_vec ) |> dplyr::count( timepoint, name = "n_cells" )
print(time_counts)
if ( !all( c("1 dpi", "3 dpi", "7 dpi") %in% time_counts$timepoint ) ) {
  stop( "Primary injury trajectory does not contain all 1/3/7 dpi timepoints.", call. = FALSE )
}
root_candidate_cells <- colnames(obj)[ time_vec == "1 dpi" & cluster_vec %in% EARLY_ROOT_CLUSTERS ]
root_candidate_tbl <- tibble::tibble( cell = root_candidate_cells, sample = sample_vec[ match( root_candidate_cells, colnames(obj) ) ], cluster = cluster_vec[ match( root_candidate_cells, colnames(obj) ) ] )
root_samples <- sort( unique( root_candidate_tbl$sample ) )
if (length(root_samples) < 2) {
  stop( "1 dpi early-activated root candidates contain fewer than two biological samples.", call. = FALSE )
}
root_candidate_counts <- root_candidate_tbl |> dplyr::count( sample, cluster, name = "n_cells" )
harmony <- Seurat::Embeddings( obj_full, reduction = HARMONY_REDUCTION )
if (ncol(harmony) < max(HARMONY_DIMS)) {
  stop( paste0( HARMONY_REDUCTION, " has fewer than 30 dimensions." ), call. = FALSE )
}
harm_root <- harmony[ root_candidate_cells, HARMONY_DIMS, drop = FALSE ]
harm_mu <- colMeans(harm_root)
harm_sd <- apply( harm_root, 2, stats::sd )
harm_sd[ !is.finite(harm_sd) | harm_sd == 0 ] <- 1
harm_z <- sweep( harm_root, 2, harm_mu, "-" )
harm_z <- sweep( harm_z, 2, harm_sd, "/" )
sample_centroids <- lapply( root_samples, function(smp) {
    cc <- root_candidate_tbl$cell[ root_candidate_tbl$sample == smp ]
    colMeans( harm_z[ cc, , drop = FALSE ] )
  }
)
names(sample_centroids) <- root_samples
consensus_centroid <- Reduce("+", sample_centroids) / length(sample_centroids)
root_candidate_tbl$distance_to_consensus <- euclidean_to_vector( harm_z[ root_candidate_tbl$cell, , drop = FALSE ], consensus_centroid )
cross_list <- lapply( root_samples, function(target_sample) {
    target_cells <- root_candidate_tbl$cell[ root_candidate_tbl$sample == target_sample ]
    other_samples <- setdiff( root_samples, target_sample )
    target_mat <- harm_z[ target_cells, , drop = FALSE ]
    other_dist <- lapply( other_samples, function(other_sample) {
        other_cells <- root_candidate_tbl$cell[ root_candidate_tbl$sample == other_sample ]
        other_mat <- harm_z[ other_cells, , drop = FALSE ]
        nn <- RANN::nn2( data = other_mat, query = target_mat, k = 1 )
        as.numeric( nn$nn.dists[, 1] )
      }
    )
    dist_mat <- do.call( cbind, other_dist )
    tibble::tibble( cell = target_cells, nearest_other_sample_mean = rowMeans(dist_mat), nearest_other_sample_max = apply( dist_mat, 1, max ) )
  }
) |> dplyr::bind_rows()
root_candidate_tbl <- root_candidate_tbl |> dplyr::left_join( cross_list, by = "cell" ) |> dplyr::group_by(sample) |> dplyr::mutate( consensus_centrality_percentile = 1 - dplyr::percent_rank( distance_to_consensus ), cross_sample_proximity_percentile = 1 - dplyr::percent_rank( nearest_other_sample_max ) ) |> dplyr::ungroup() |> dplyr::mutate( highconf_root = consensus_centrality_percentile >= CENTRALITY_PCT_MIN & cross_sample_proximity_percentile >= CROSS_SAMPLE_PCT_MIN, root_priority = W_CENTRAL * consensus_centrality_percentile + W_CROSS * cross_sample_proximity_percentile )
highconf_root_counts <- root_candidate_tbl |> dplyr::filter( highconf_root ) |> dplyr::count( sample, name = "n_highconf" ) |> tidyr::complete( sample = root_samples, fill = list( n_highconf = 0 ) )
print(highconf_root_counts)
min_highconf_root <- min( highconf_root_counts$n_highconf )
if ( min_highconf_root < MIN_ROOT_PER_SAMPLE ) {
  stop( paste0( "最少1 dpi样本只有 ", min_highconf_root, " 个高置信root candidates (<", MIN_ROOT_PER_SAMPLE, ")。\n", "停止运行并检查 root-candidate 阈值。" ), call. = FALSE )
}
root_n_per_sample <- min( TARGET_ROOT_PER_SAMPLE, min_highconf_root )
primary_root_tbl <- root_candidate_tbl |> dplyr::filter( highconf_root ) |> dplyr::group_by(sample) |> dplyr::arrange( dplyr::desc(root_priority), .by_group = TRUE ) |> dplyr::mutate( within_sample_rank = dplyr::row_number() ) |> dplyr::filter( within_sample_rank <= root_n_per_sample ) |> dplyr::ungroup()
primary_root_cells <- primary_root_tbl$cell
sensitivity_root_cells <- root_candidate_tbl$cell[ root_candidate_tbl$highconf_root ]
formal_umap_full <- Seurat::Embeddings( obj_full, reduction = UMAP_REDUCTION )
formal_umap <- formal_umap_full[ colnames(obj), 1:2, drop = FALSE ]
root_plot_df <- tibble::tibble( cell = root_candidate_cells, UMAP_1 = formal_umap_full[ root_candidate_cells, 1 ], UMAP_2 = formal_umap_full[ root_candidate_cells, 2 ] ) |> dplyr::left_join( root_candidate_tbl, by = "cell" ) |> dplyr::mutate( primary_root = cell %in% primary_root_cells )
injury_plot_df <- tibble::tibble( cell = colnames(obj), UMAP_1 = formal_umap[, 1], UMAP_2 = formal_umap[, 2], timepoint = time_vec, sample = sample_vec, cluster = cluster_vec, primary_root = colnames(obj) %in% primary_root_cells )
injury_samples <- sort( unique( injury_plot_df$sample ) )
sample_palette <- grDevices::hcl.colors( length(injury_samples), palette = "Dark 3" )
names(sample_palette) <- injury_samples
root_sample_palette <- grDevices::hcl.colors( length(root_samples), palette = "Dark 3" )
names(root_sample_palette) <- root_samples
theme5 <- ggplot2::theme_classic( base_size = 11 ) + ggplot2::theme( axis.text = ggplot2::element_text( color = "#222222" ), axis.title = ggplot2::element_text( color = "#222222" ), plot.title = ggplot2::element_text( face = "bold", color = "#111111" ) )
p1a <- ggplot2::ggplot( injury_plot_df, ggplot2::aes( UMAP_1, UMAP_2, color = factor(timepoint) ) ) + ggplot2::geom_point( size = 0.45, alpha = 0.60 ) + ggplot2::scale_color_manual( values = TIME_COLORS, name = "True time" ) + ggplot2::labs( title = "Primary injury-response trajectory input", subtitle = "1 dpi early activated -> 3 dpi transition/DAM -> 7 dpi DAM landscape" ) + theme5
p1b <- ggplot2::ggplot( injury_plot_df, ggplot2::aes( UMAP_1, UMAP_2, color = factor(cluster) ) ) + ggplot2::geom_point( size = 0.45, alpha = 0.60 ) + ggplot2::scale_color_manual( values = CLUSTER_COLORS, name = "Step3A cluster" ) + ggplot2::labs( title = "Step3A states in injury trajectory" ) + theme5
p1c <- ggplot2::ggplot( root_plot_df, ggplot2::aes( UMAP_1, UMAP_2, color = factor(sample) ) ) + ggplot2::geom_point( size = 0.65, alpha = 0.70 ) + ggplot2::scale_color_manual( values = root_sample_palette, name = "1 dpi sample" ) + ggplot2::geom_point( data = root_plot_df |> dplyr::filter( primary_root ), inherit.aes = FALSE, ggplot2::aes( UMAP_1, UMAP_2 ), color = ROOT_COL, size = 0.95, alpha = 0.95 ) + ggplot2::labs( title = "Sample-balanced 1 dpi consensus root", subtitle = paste0( root_n_per_sample, " cells/sample; red = primary root" ) ) + theme5
root_metric_long <- root_candidate_tbl |> dplyr::select( sample, consensus_centrality_percentile, cross_sample_proximity_percentile, root_priority ) |> tidyr::pivot_longer( cols = c( consensus_centrality_percentile, cross_sample_proximity_percentile, root_priority ), names_to = "metric", values_to = "value" ) |> dplyr::mutate( metric = dplyr::recode( metric, consensus_centrality_percentile = "Harmony consensus centrality", cross_sample_proximity_percentile = "Cross-sample proximity", root_priority = "Composite root priority" ) )
p1d <- ggplot2::ggplot( root_metric_long, ggplot2::aes( x = factor(sample), y = value, fill = factor(sample) ) ) + ggplot2::geom_violin( trim = TRUE, scale = "width", alpha = 0.75, linewidth = 0.25 ) + ggplot2::geom_boxplot( width = 0.12, outlier.shape = NA, fill = "#FFFFFF", linewidth = 0.25 ) + ggplot2::facet_wrap( ~metric, nrow = 1 ) + ggplot2::scale_fill_manual( values = root_sample_palette, guide = "none" ) + ggplot2::labs( title = "1 dpi root-consensus metrics", x = "Biological sample", y = "Percentile / priority" ) + theme5 + ggplot2::theme( axis.text.x = ggplot2::element_text( angle = 25, hjust = 1 ) )
grDevices::pdf( OUT_FIG1, width = 14, height = 10, useDingbats = FALSE )
print( (p1a | p1b) / (p1c | p1d) + patchwork::plot_annotation( title = "Step 5B-R | Injury-response input and 1 dpi root" ) )
grDevices::dev.off()
counts <- get_counts_matrix(obj)
if (!identical( colnames(counts), colnames(obj) )) {
  counts <- counts[ , colnames(obj), drop = FALSE ]
}
cell_metadata <- obj[[]]
cell_metadata <- cell_metadata[ colnames(counts), , drop = FALSE ]
gene_metadata <- data.frame( gene_short_name = rownames(counts), row.names = rownames(counts), stringsAsFactors = FALSE )
cds <- monocle3::new_cell_data_set( counts, cell_metadata = cell_metadata, gene_metadata = gene_metadata )
cds <- monocle3::preprocess_cds( cds, num_dim = NUM_DIM, method = "PCA" )
colnames(formal_umap) <- c( "UMAP_1", "UMAP_2" )
assign_ok <- tryCatch( {
    SingleCellExperiment::reducedDim( cds, "UMAP" ) <- formal_umap
    TRUE
  },
  error = function(e) FALSE
)
if (!assign_ok) {
  cds <- monocle3::reduce_dimension( cds, reduction_method = "UMAP", preprocess_method = "PCA" )
  cds@int_colData$reducedDims$UMAP <- formal_umap
}
imported_umap <- SingleCellExperiment::reducedDim( cds, "UMAP" )
umap_max_abs_diff <- max( abs( imported_umap - formal_umap ), na.rm = TRUE )
if (umap_max_abs_diff > 1e-8) {
  stop( "Imported Monocle UMAP does not exactly match formal Step3A UMAP.", call. = FALSE )
}
cds <- monocle3::cluster_cells( cds, reduction_method = "UMAP", k = MONOCLE_K, cluster_method = "leiden", random_seed = 20260731, verbose = FALSE )
mono_cluster <- get_clusters_safe(cds)
mono_partition <- get_partitions_safe(cds)
colData(cds)$monocle3_cluster <- as.character(mono_cluster)
colData(cds)$monocle3_partition <- as.character(mono_partition)
cds <- monocle3::learn_graph( cds, use_partition = FALSE, close_loop = LEARN_GRAPH_CLOSE_LOOP, learn_graph_control = list( euclidean_distance_ratio = LEARN_GRAPH_EUCLIDEAN_DISTANCE_RATIO ) )
pg <- monocle3::principal_graph( cds )[["UMAP"]]
if (is.null(pg)) {
  stop( "Monocle3 principal graph was not generated.", call. = FALSE )
}
graph_components <- igraph::components(pg)$no
cds <- monocle3::order_cells( cds, reduction_method = "UMAP", root_cells = primary_root_cells )
pt_primary <- get_pseudotime_safe( cds )
pt_primary <- pt_primary[ colnames(cds) ]
finite_fraction <- mean( is.finite( pt_primary ) )
cds_sens <- monocle3::order_cells( cds, reduction_method = "UMAP", root_cells = sensitivity_root_cells )
pt_sens <- get_pseudotime_safe( cds_sens )
pt_sens <- pt_sens[ colnames(cds_sens) ]
both_finite <- ( is.finite(pt_primary) & is.finite(pt_sens) )
root_sensitivity_spearman <- suppressWarnings( stats::cor( pt_primary[both_finite], pt_sens[both_finite], method = "spearman" ) )
colData(cds)$pseudotime_primary <- as.numeric( pt_primary )
colData(cds)$pseudotime_sensitivity <- as.numeric( pt_sens )
cell_df <- tibble::tibble( cell = colnames(cds), timepoint = as.character( colData(cds)[[TIME_COL]] ), sample = as.character( colData(cds)[[SAMPLE_COL]] ), cluster = as.character( colData(cds)[[CLUSTER_COL]] ), shared_transition_3dpi_highconf = as.logical( colData(cds)$shared_transition_3dpi_highconf ), transition_status = as.character( colData(cds)$transition_status ), repair_candidate_score = as.numeric( colData(cds)[[REPAIR_SCORE_COL]] ), scar_candidate_score = as.numeric( colData(cds)[[SCAR_SCORE_COL]] ), primary_root = colnames(cds) %in% primary_root_cells, sensitivity_root = colnames(cds) %in% sensitivity_root_cells, pseudotime = as.numeric( pt_primary ), pseudotime_sensitivity = as.numeric( pt_sens ), monocle3_cluster = as.character( mono_cluster ), monocle3_partition = as.character( mono_partition ) )
cell_df$repair_minus_scar <- cell_df$repair_candidate_score - cell_df$scar_candidate_score
injury_time_levels <- c( "1 dpi", "3 dpi", "7 dpi" )
injury_day <- c( "1 dpi" = 1, "3 dpi" = 3, "7 dpi" = 7 )
cell_df$true_day <- unname( injury_day[ cell_df$timepoint ] )
time_summary <- cell_df |> dplyr::filter( is.finite( pseudotime ) ) |> dplyr::mutate( timepoint = factor( timepoint, levels = injury_time_levels, ordered = TRUE ) ) |> dplyr::group_by( timepoint ) |> dplyr::summarise( n_cells = dplyr::n(), pseudotime_mean = mean(pseudotime), pseudotime_median = median(pseudotime), pseudotime_q25 = stats::quantile( pseudotime, 0.25 ), pseudotime_q75 = stats::quantile( pseudotime, 0.75 ), .groups = "drop" )
time_medians <- time_summary$pseudotime_median
names(time_medians) <- as.character( time_summary$timepoint )
time_monotonic <- all( diff( time_medians[ injury_time_levels ] ) > 0 )
cell_level_time_spearman <- suppressWarnings( stats::cor( cell_df$true_day[ is.finite( cell_df$pseudotime ) ], cell_df$pseudotime[ is.finite( cell_df$pseudotime ) ], method = "spearman" ) )
sample_pt_summary <- cell_df |> dplyr::filter( is.finite( pseudotime ) ) |> dplyr::group_by( timepoint, sample ) |> dplyr::summarise( n_cells = dplyr::n(), median_pseudotime = median( pseudotime ), mean_pseudotime = mean( pseudotime ), .groups = "drop" ) |> dplyr::mutate( true_day = unname( injury_day[ timepoint ] ) )
sample_level_time_spearman <- suppressWarnings( stats::cor( sample_pt_summary$true_day, sample_pt_summary$median_pseudotime, method = "spearman" ) )
reference_df <- cell_df |> dplyr::mutate( reference_state = dplyr::case_when( primary_root ~ "1dpi_consensus_root", timepoint == "1 dpi" & cluster %in% EARLY_ROOT_CLUSTERS ~ "1dpi_early_activated", timepoint == "3 dpi" & shared_transition_3dpi_highconf ~ "3dpi_shared_transition", timepoint == "7 dpi" & cluster == "1" & transition_status == "Mature_DAM_like" ~ "7dpi_established_DAM", TRUE ~ "Other" ) )
reference_summary <- reference_df |> dplyr::filter( reference_state != "Other", is.finite( pseudotime ) ) |> dplyr::group_by( reference_state ) |> dplyr::summarise( n_cells = dplyr::n(), pseudotime_mean = mean( pseudotime ), pseudotime_median = median( pseudotime ), pseudotime_q25 = stats::quantile( pseudotime, 0.25 ), pseudotime_q75 = stats::quantile( pseudotime, 0.75 ), .groups = "drop" )
closest_vertex <- get_closest_vertex( cds, pg )
closest_vertex <- closest_vertex[ cell_df$cell ]
cell_df$principal_vertex <- as.character( closest_vertex[ cell_df$cell ] )
root_vertex_counts <- tibble::tibble( cell = primary_root_cells, principal_vertex = closest_vertex[ primary_root_cells ] ) |> dplyr::count( principal_vertex, name = "n_root_cells" ) |> dplyr::arrange( dplyr::desc( n_root_cells ) )
if (nrow(root_vertex_counts) == 0) {
  stop( "Primary root cells could not be mapped to principal graph vertices.", call. = FALSE )
}
ROOT_ANCHOR_VERTEX <- root_vertex_counts$principal_vertex[1]
degree_vec <- igraph::degree( pg )
branch_vertices <- names( degree_vec[ degree_vec >= 3 ] )
leaf_vertices <- names( degree_vec[ degree_vec == 1 ] )
root_dist <- igraph::distances( pg, v = ROOT_ANCHOR_VERTEX, to = igraph::V(pg) )
graph_vertex_table <- tibble::tibble( vertex = names( degree_vec ), degree = as.integer( degree_vec ), is_root_anchor = names(degree_vec) == ROOT_ANCHOR_VERTEX, is_branch = names(degree_vec) %in% branch_vertices, is_leaf = names(degree_vec) %in% leaf_vertices, graph_distance_from_root = as.numeric( root_dist[ 1, names(degree_vec) ] ) )
branchpoint_rows <- list()
branch_component_cache <- list()
if (length(branch_vertices) > 0) {
  for (bp in branch_vertices) {
    g2 <- igraph::delete_vertices( pg, bp )
    comp <- igraph::components( g2 )
    membership <- comp$membership
    names(membership) <- igraph::V(g2)$name
    root_side_component <- membership[ ROOT_ANCHOR_VERTEX ]
    if ( length(root_side_component) == 0 || is.na(root_side_component) ) {
      branchpoint_rows[[bp]] <- tibble::tibble( branch_vertex = bp, degree = igraph::degree( pg, v = bp ), downstream_components = NA_integer_, downstream_major_arms = NA_integer_, min_top2_7dpi_cells = NA_integer_, top2_7dpi_sum = NA_integer_, near_cells = NA_integer_, near_3dpi_cells = NA_integer_, near_transition_cells = NA_integer_, near_transition_fraction = NA_real_ )
      next
    }
    cv <- closest_vertex[ cell_df$cell ]
    comp_for_cell <- membership[ cv ]
    tmp <- cell_df |> dplyr::mutate( graph_component = as.integer( comp_for_cell[ cell ] ) )
    downstream <- tmp |> dplyr::filter( !is.na( graph_component ), graph_component != root_side_component ) |> dplyr::group_by( graph_component ) |> dplyr::summarise( n_cells = dplyr::n(), n_1dpi = sum( timepoint == "1 dpi" ), n_3dpi = sum( timepoint == "3 dpi" ), n_7dpi = sum( timepoint == "7 dpi" ), n_transition = sum( shared_transition_3dpi_highconf, na.rm = TRUE ), median_pseudotime = median( pseudotime, na.rm = TRUE ), .groups = "drop" ) |> dplyr::arrange( dplyr::desc( n_7dpi ) )
    major <- downstream |> dplyr::filter( n_cells >= MIN_TOTAL_CELLS_PER_DOWNSTREAM_ARM, n_7dpi >= MIN_7DPI_CELLS_PER_DOWNSTREAM_ARM )
    top2 <- head( major$n_7dpi, 2 )
    near_cells <- cells_near_graph_vertex( pg, closest_vertex, bp, radius = BRANCHPOINT_GRAPH_RADIUS )
    near_df <- cell_df |> dplyr::filter( cell %in% near_cells )
    branchpoint_rows[[bp]] <- tibble::tibble( branch_vertex = bp, degree = igraph::degree( pg, v = bp ), downstream_components = nrow( downstream ), downstream_major_arms = nrow( major ), min_top2_7dpi_cells = ifelse( length(top2) >= 2, min(top2), 0 ), top2_7dpi_sum = sum(top2), near_cells = nrow( near_df ), near_3dpi_cells = sum( near_df$timepoint == "3 dpi" ), near_transition_cells = sum( near_df$shared_transition_3dpi_highconf, na.rm = TRUE ), near_transition_fraction = ifelse( nrow(near_df) > 0, mean( near_df$shared_transition_3dpi_highconf, na.rm = TRUE ), NA_real_ ) )
    branch_component_cache[[bp]] <- list( membership = membership, root_side_component = root_side_component, downstream = downstream, major = major )
  }
}
branchpoint_summary <- if ( length(branchpoint_rows) > 0 ) {
  dplyr::bind_rows( branchpoint_rows ) |> dplyr::left_join( graph_vertex_table |> dplyr::select( vertex, graph_distance_from_root ), by = c( "branch_vertex" = "vertex" ) ) |> dplyr::mutate( candidate_bifurcation = !is.na( downstream_major_arms ) & downstream_major_arms >= 2 ) |> dplyr::arrange( dplyr::desc( candidate_bifurcation ), dplyr::desc( min_top2_7dpi_cells ), dplyr::desc( near_transition_cells ), graph_distance_from_root )
} else {
  tibble::tibble()
}
primary_branchpoint <- NA_character_
if ( nrow(branchpoint_summary) > 0 && any( branchpoint_summary$candidate_bifurcation, na.rm = TRUE ) ) {
  primary_branchpoint <- branchpoint_summary$branch_vertex[ which( branchpoint_summary$candidate_bifurcation )[1] ]
} else {
}
cell_df$branch_arm <- NA_character_
cell_df$major_late_arm <- FALSE
arm_summary <- tibble::tibble()
branch_program_summary <- tibble::tibble()
branch_program_sample_summary <- tibble::tibble()
if (!is.na(primary_branchpoint)) {
  cache <- branch_component_cache[[ primary_branchpoint ]]
  membership <- cache$membership
  root_side_component <- cache$root_side_component
  cv <- closest_vertex[ cell_df$cell ]
  comp_for_cell <- membership[ cv ]
  cell_df$branch_component <- as.integer( comp_for_cell[ cell_df$cell ] )
  downstream_raw <- cell_df |> dplyr::filter( !is.na( branch_component ), branch_component != root_side_component ) |> dplyr::group_by( branch_component ) |> dplyr::summarise( n_cells = dplyr::n(), n_1dpi = sum( timepoint == "1 dpi" ), n_3dpi = sum( timepoint == "3 dpi" ), n_7dpi = sum( timepoint == "7 dpi" ), n_transition = sum( shared_transition_3dpi_highconf, na.rm = TRUE ), pseudotime_median = median( pseudotime, na.rm = TRUE ), .groups = "drop" ) |> dplyr::arrange( dplyr::desc( n_7dpi ) ) |> dplyr::mutate( arm_rank = dplyr::row_number(), branch_arm = paste0( "Arm_", arm_rank ), major_late_arm = n_cells >= MIN_TOTAL_CELLS_PER_DOWNSTREAM_ARM & n_7dpi >= MIN_7DPI_CELLS_PER_DOWNSTREAM_ARM )
  arm_map <- stats::setNames( downstream_raw$branch_arm, downstream_raw$branch_component )
  major_map <- stats::setNames( downstream_raw$major_late_arm, downstream_raw$branch_component )
  root_side_idx <- which( !is.na( cell_df$branch_component ) & cell_df$branch_component == root_side_component )
  cell_df$branch_arm[ root_side_idx ] <- "Pre_branch_root_side"
  downstream_idx <- which( !is.na( cell_df$branch_component ) & cell_df$branch_component != root_side_component )
  cell_df$branch_arm[ downstream_idx ] <- unname( arm_map[ as.character( cell_df$branch_component[ downstream_idx ] ) ] )
  cell_df$major_late_arm[ downstream_idx ] <- unname( major_map[ as.character( cell_df$branch_component[ downstream_idx ] ) ] )
  cell_df$branch_arm[ cell_df$principal_vertex == primary_branchpoint ] <- "Branchpoint_vertex"
  arm_summary <- downstream_raw
  branch_program_summary <- cell_df |> dplyr::filter( timepoint == "7 dpi", !is.na( branch_arm ), grepl( "^Arm_", branch_arm ) ) |> dplyr::group_by( branch_arm ) |> dplyr::summarise( n_cells = dplyr::n(), repair_mean = mean( repair_candidate_score, na.rm = TRUE ), repair_median = median( repair_candidate_score, na.rm = TRUE ), scar_mean = mean( scar_candidate_score, na.rm = TRUE ), scar_median = median( scar_candidate_score, na.rm = TRUE ), repair_minus_scar_mean = mean( repair_minus_scar, na.rm = TRUE ), repair_minus_scar_median = median( repair_minus_scar, na.rm = TRUE ), pseudotime_median = median( pseudotime, na.rm = TRUE ), .groups = "drop" )
  branch_program_sample_summary <- cell_df |> dplyr::filter( timepoint == "7 dpi", !is.na( branch_arm ), grepl( "^Arm_", branch_arm ) ) |> dplyr::group_by( sample, branch_arm ) |> dplyr::summarise( n_cells = dplyr::n(), repair_median = median( repair_candidate_score, na.rm = TRUE ), scar_median = median( scar_candidate_score, na.rm = TRUE ), repair_minus_scar_median = median( repair_minus_scar, na.rm = TRUE ), .groups = "drop" )
}
colData(cds)$branch_arm <- cell_df$branch_arm
colData(cds)$major_late_arm <- cell_df$major_late_arm
obj_full$step5BR_injury_input <- FALSE
obj_full$step5BR_primary_root <- FALSE
obj_full$step5BR_sensitivity_root <- FALSE
obj_full$step5BR_pseudotime <- NA_real_
obj_full$step5BR_pseudotime_sensitivity <- NA_real_
obj_full$step5BR_monocle3_partition <- NA_character_
obj_full$step5BR_monocle3_cluster <- NA_character_
obj_full$step5BR_principal_vertex <- NA_character_
obj_full$step5BR_branch_arm <- NA_character_
obj_full$step5BR_major_late_arm <- FALSE
obj_full$step5BR_injury_input[ match( colnames(obj), colnames(obj_full) ) ] <- TRUE
obj_full$step5BR_primary_root[ match( primary_root_cells, colnames(obj_full) ) ] <- TRUE
obj_full$step5BR_sensitivity_root[ match( sensitivity_root_cells, colnames(obj_full) ) ] <- TRUE
idx_full <- match( cell_df$cell, colnames(obj_full) )
obj_full$step5BR_pseudotime[ idx_full ] <- cell_df$pseudotime
obj_full$step5BR_pseudotime_sensitivity[ idx_full ] <- cell_df$pseudotime_sensitivity
obj_full$step5BR_monocle3_partition[ idx_full ] <- cell_df$monocle3_partition
obj_full$step5BR_monocle3_cluster[ idx_full ] <- cell_df$monocle3_cluster
obj_full$step5BR_principal_vertex[ idx_full ] <- cell_df$principal_vertex
obj_full$step5BR_branch_arm[ idx_full ] <- cell_df$branch_arm
obj_full$step5BR_major_late_arm[ idx_full ] <- cell_df$major_late_arm
obj_full$step5BR_primary_branchpoint <- ifelse( is.na( primary_branchpoint ), NA_character_, primary_branchpoint )
p2a <- monocle3::plot_cells( cds, reduction_method = "UMAP", color_cells_by = TIME_COL, show_trajectory_graph = TRUE, label_cell_groups = FALSE, label_leaves = TRUE, label_branch_points = TRUE, graph_label_size = 3, cell_size = 0.45 ) + ggplot2::scale_color_manual( values = TIME_COLORS ) + ggplot2::ggtitle( "Injury-response principal graph by true time" )
p2b <- monocle3::plot_cells( cds, reduction_method = "UMAP", color_cells_by = "pseudotime", show_trajectory_graph = TRUE, label_cell_groups = FALSE, label_leaves = TRUE, label_branch_points = TRUE, graph_label_size = 3, cell_size = 0.45 ) + ggplot2::ggtitle( "Primary injury-response pseudotime" )
pt_violin_df <- cell_df |> dplyr::filter( is.finite( pseudotime ) ) |> dplyr::mutate( timepoint = factor( timepoint, levels = injury_time_levels, ordered = TRUE ) )
p2c <- ggplot2::ggplot( pt_violin_df, ggplot2::aes( x = timepoint, y = pseudotime, fill = timepoint ) ) + ggplot2::geom_violin( trim = TRUE, scale = "width", alpha = 0.75, linewidth = 0.25 ) + ggplot2::geom_boxplot( width = 0.12, outlier.shape = NA, fill = "#FFFFFF", linewidth = 0.25 ) + ggplot2::scale_fill_manual( values = TIME_COLORS, guide = "none" ) + ggplot2::labs( title = "Pseudotime versus true injury time", subtitle = paste0( "Cell-level Spearman = ", round( cell_level_time_spearman, 3 ), " | Monotonic medians = ", time_monotonic ), x = "True time", y = "Pseudotime" ) + theme5
reference_order <- c( "1dpi_consensus_root", "1dpi_early_activated", "3dpi_shared_transition", "7dpi_established_DAM" )
reference_plot_df <- reference_df |> dplyr::mutate( reference_state = factor( reference_state, levels = reference_order, ordered = TRUE ) ) |> dplyr::filter( !is.na( reference_state ), reference_state != "Other", is.finite( pseudotime ) )
p2d <- ggplot2::ggplot( reference_plot_df, ggplot2::aes( x = reference_state, y = pseudotime ) ) + ggplot2::geom_violin( fill = "#9ECAE1", trim = TRUE, scale = "width", linewidth = 0.25 ) + ggplot2::geom_boxplot( width = 0.12, outlier.shape = NA, fill = "#FFFFFF", linewidth = 0.25 ) + ggplot2::labs( title = "Formal biological references along injury pseudotime", x = NULL, y = "Pseudotime" ) + ggplot2::theme_classic( base_size = 10 ) + ggplot2::theme( axis.text.x = ggplot2::element_text( angle = 25, hjust = 1 ) )
grDevices::pdf( OUT_FIG2, width = 14, height = 10, useDingbats = FALSE )
print( (p2a | p2b) / (p2c | p2d) + patchwork::plot_annotation( title = "Step 5B-R | Injury-response Monocle3 pseudotime" ) )
grDevices::dev.off()
if (!is.na(primary_branchpoint)) {
  colData(cds)$branch_arm_display <- ifelse( is.na( cell_df$branch_arm ), "Unassigned", cell_df$branch_arm )
  p3a <- monocle3::plot_cells( cds, reduction_method = "UMAP", color_cells_by = "branch_arm_display", show_trajectory_graph = TRUE, label_cell_groups = FALSE, label_leaves = TRUE, label_branch_points = TRUE, graph_label_size = 3, cell_size = 0.45 ) + ggplot2::ggtitle( paste0( "Objective downstream arms from ", primary_branchpoint ) )
  arm_time_df <- cell_df |> dplyr::filter( !is.na( branch_arm ), grepl( "^Arm_", branch_arm ) ) |> dplyr::count( branch_arm, timepoint, name = "n_cells" ) |> dplyr::group_by( branch_arm ) |> dplyr::mutate( fraction = n_cells / sum(n_cells) ) |> dplyr::ungroup()
  p3b <- ggplot2::ggplot( arm_time_df, ggplot2::aes( x = branch_arm, y = fraction, fill = timepoint ) ) + ggplot2::geom_col( width = 0.70 ) + ggplot2::scale_fill_manual( values = TIME_COLORS, name = "True time" ) + ggplot2::labs( title = "True-time composition of graph arms", x = "Graph arm", y = "Fraction" ) + theme5
  near_bp_cells <- cells_near_graph_vertex( pg, closest_vertex, primary_branchpoint, radius = BRANCHPOINT_GRAPH_RADIUS )
  near_bp_df <- cell_df |> dplyr::filter( cell %in% near_bp_cells )
  p3c <- ggplot2::ggplot( near_bp_df, ggplot2::aes( x = factor( timepoint, levels = injury_time_levels ), fill = shared_transition_3dpi_highconf ) ) + ggplot2::geom_bar() + ggplot2::scale_fill_manual( values = c( "FALSE" = "#D9D9D9", "TRUE" = "#D73027" ), name = "High-conf\n3 dpi transition" ) + ggplot2::labs( title = paste0( "Cells near ", primary_branchpoint ), x = "True time", y = "Cells" ) + theme5
} else {
  p3a <- monocle3::plot_cells( cds, reduction_method = "UMAP", color_cells_by = "pseudotime", show_trajectory_graph = TRUE, label_cell_groups = FALSE, label_leaves = TRUE, label_branch_points = TRUE, graph_label_size = 3, cell_size = 0.45 ) + ggplot2::ggtitle( "No objective two-major-7dpi-arm bifurcation detected" )
  placeholder <- tibble::tibble( x = 1, y = 1 )
  p3b <- ggplot2::ggplot( placeholder, ggplot2::aes(x, y) ) + ggplot2::annotate( "text", x = 1, y = 1, label = "No branch-arm assignment was forced.", size = 5 ) + ggplot2::theme_void()
  p3c <- p3b
}
if (nrow(branchpoint_summary) > 0) {
  bp_plot_df <- branchpoint_summary |> dplyr::mutate( branch_vertex = factor( branch_vertex, levels = rev( branch_vertex ) ) )
  p3d <- ggplot2::ggplot( bp_plot_df, ggplot2::aes( x = branch_vertex, y = min_top2_7dpi_cells, fill = candidate_bifurcation ) ) + ggplot2::geom_col( width = 0.70 ) + ggplot2::coord_flip() + ggplot2::scale_fill_manual( values = c( "FALSE" = "#BDBDBD", "TRUE" = "#2C7BB6" ), name = "Two major\n7 dpi arms" ) + ggplot2::labs( title = "Objective branch-point ranking", subtitle = "Minimum 7 dpi support across two largest downstream arms", x = "Principal-graph branch vertex", y = "Minimum top-two 7 dpi cells" ) + ggplot2::theme_classic( base_size = 10 )
} else {
  p3d <- ggplot2::ggplot( tibble::tibble( x = 1, y = 1 ), ggplot2::aes(x, y) ) + ggplot2::annotate( "text", x = 1, y = 1, label = "No degree >=3 branch vertices detected.", size = 5 ) + ggplot2::theme_void()
}
grDevices::pdf( OUT_FIG3, width = 14, height = 10, useDingbats = FALSE )
print( (p3a | p3b) / (p3c | p3d) + patchwork::plot_annotation( title = "Step 5B-R | Injury-response branch-point discovery" ) )
grDevices::dev.off()
p4a <- monocle3::plot_cells( cds, reduction_method = "UMAP", color_cells_by = REPAIR_SCORE_COL, show_trajectory_graph = TRUE, label_cell_groups = FALSE, label_leaves = FALSE, label_branch_points = TRUE, graph_label_size = 3, cell_size = 0.45 ) + ggplot2::ggtitle( "Post-hoc repair-related candidate program" )
p4b <- monocle3::plot_cells( cds, reduction_method = "UMAP", color_cells_by = SCAR_SCORE_COL, show_trajectory_graph = TRUE, label_cell_groups = FALSE, label_leaves = FALSE, label_branch_points = TRUE, graph_label_size = 3, cell_size = 0.45 ) + ggplot2::ggtitle( "Post-hoc scar-related candidate program" )
colData(cds)$repair_minus_scar <- cell_df$repair_minus_scar
p4c <- monocle3::plot_cells( cds, reduction_method = "UMAP", color_cells_by = "repair_minus_scar", show_trajectory_graph = TRUE, label_cell_groups = FALSE, label_leaves = FALSE, label_branch_points = TRUE, graph_label_size = 3, cell_size = 0.45 ) + ggplot2::ggtitle( "Post-hoc repair minus scar coordinate" )
if ( !is.na( primary_branchpoint ) && nrow( branch_program_sample_summary ) > 0 ) {
  branch_long <- branch_program_sample_summary |> dplyr::select( sample, branch_arm, repair_median, scar_median ) |> tidyr::pivot_longer( cols = c( repair_median, scar_median ), names_to = "program", values_to = "sample_median_score" ) |> dplyr::mutate( program = dplyr::recode( program, repair_median = "Repair-related", scar_median = "Scar-related" ) )
  p4d <- ggplot2::ggplot( branch_long, ggplot2::aes( x = branch_arm, y = sample_median_score, shape = sample ) ) + ggplot2::geom_point( size = 3 ) + ggplot2::facet_wrap( ~program, nrow = 1, scales = "free_y" ) + ggplot2::labs( title = "7 dpi arm-level programs by biological sample", subtitle = "Post-hoc only; topology was fixed before these scores were read", x = "Graph arm", y = "Sample median score" ) + ggplot2::theme_classic( base_size = 10 )
} else {
  p4d <- ggplot2::ggplot( tibble::tibble( x = 1, y = 1 ), ggplot2::aes(x, y) ) + ggplot2::annotate( "text", x = 1, y = 1, label = "No arm-level repair/scar comparison:\nno objective bifurcation was forced.", size = 4.5 ) + ggplot2::theme_void()
}
grDevices::pdf( OUT_FIG4, width = 14, height = 10, useDingbats = FALSE )
print( (p4a | p4b) / (p4c | p4d) + patchwork::plot_annotation( title = "Step 5B-R | Step4 programs projected after injury-graph learning" ) )
grDevices::dev.off()
sens_df <- cell_df |> dplyr::filter( is.finite( pseudotime ), is.finite( pseudotime_sensitivity ) )
p5a <- ggplot2::ggplot( sens_df, ggplot2::aes( x = pseudotime, y = pseudotime_sensitivity ) ) + ggplot2::geom_point( size = 0.45, alpha = 0.30 ) + ggplot2::labs( title = "Primary versus broad-root injury pseudotime", subtitle = paste0( "Spearman = ", round( root_sensitivity_spearman, 4 ) ), x = "Primary pseudotime", y = "Sensitivity pseudotime" ) + theme5
sens_time <- sens_df |> dplyr::group_by( timepoint ) |> dplyr::summarise( primary_median = median( pseudotime ), sensitivity_median = median( pseudotime_sensitivity ), .groups = "drop" ) |> tidyr::pivot_longer( cols = c( primary_median, sensitivity_median ), names_to = "root_definition", values_to = "median_pseudotime" ) |> dplyr::mutate( timepoint = factor( timepoint, levels = injury_time_levels, ordered = TRUE ), root_definition = dplyr::recode( root_definition, primary_median = "Primary balanced 1 dpi root", sensitivity_median = "Broad high-confidence 1 dpi root" ) )
p5b <- ggplot2::ggplot( sens_time, ggplot2::aes( x = timepoint, y = median_pseudotime, group = root_definition, linetype = root_definition ) ) + ggplot2::geom_line( linewidth = 0.8 ) + ggplot2::geom_point( size = 2.5 ) + ggplot2::labs( title = "True-time ordering under two 1 dpi root definitions", x = "True time", y = "Median pseudotime", linetype = "Root definition" ) + theme5
root_set_plot <- injury_plot_df |> dplyr::mutate( sensitivity_root = cell %in% sensitivity_root_cells )
p5c <- ggplot2::ggplot( root_set_plot, ggplot2::aes( UMAP_1, UMAP_2 ) ) + ggplot2::geom_point( color = BG_COL, size = 0.35, alpha = 0.25 ) + ggplot2::geom_point( data = root_set_plot |> dplyr::filter( sensitivity_root ), color = "#2C7BB6", size = 0.65, alpha = 0.55 ) + ggplot2::geom_point( data = root_set_plot |> dplyr::filter( primary_root ), color = ROOT_COL, size = 0.95, alpha = 0.95 ) + ggplot2::labs( title = "Primary and sensitivity 1 dpi root sets", subtitle = "Red = primary balanced root; blue = broader high-confidence root" ) + theme5
normalized_diff <- sens_df |> dplyr::mutate( primary_scaled = ( pseudotime - min( pseudotime ) ) / ( max( pseudotime ) - min( pseudotime ) ), sensitivity_scaled = ( pseudotime_sensitivity - min( pseudotime_sensitivity ) ) / ( max( pseudotime_sensitivity ) - min( pseudotime_sensitivity ) ), scaled_difference = sensitivity_scaled - primary_scaled )
p5d <- ggplot2::ggplot( normalized_diff, ggplot2::aes( x = scaled_difference ) ) + ggplot2::geom_histogram( bins = 60 ) + ggplot2::labs( title = "Normalized root-sensitivity pseudotime difference", x = "Sensitivity - primary scaled pseudotime", y = "Cells" ) + theme5
grDevices::pdf( OUT_FIG5, width = 14, height = 10, useDingbats = FALSE )
print( (p5a | p5b) / (p5c | p5d) + patchwork::plot_annotation( title = "Step 5B-R | Injury-root sensitivity" ) )
grDevices::dev.off()
parameter_tbl <- tibble::tibble( parameter = c( "Primary injury trajectory cells", "1 dpi early-root clusters", "3 dpi inclusion rule", "7 dpi DAM clusters", "1 dpi root samples", "Primary root cells/sample", "Primary root total", "Sensitivity root total", "Formal UMAP", "Formal Harmony", "Principal graph vertices", "Principal graph edges", "Principal graph connected components", "Finite pseudotime fraction", "Cell-level true-time Spearman", "Sample-level true-time Spearman", "True-time medians monotonic", "Root sensitivity Spearman", "Root principal vertices", "Root anchor vertex", "Branch vertices", "Leaf vertices", "Objective primary branchpoint" ), value = c( ncol(obj), paste( EARLY_ROOT_CLUSTERS, collapse = ";" ), "cluster 0/1/4 OR shared_transition_3dpi_highconf", paste( DAM_CLUSTERS, collapse = ";" ), paste( root_samples, collapse = ";" ), root_n_per_sample, length( primary_root_cells ), length( sensitivity_root_cells ), UMAP_REDUCTION, paste0( HARMONY_REDUCTION, " dims 1:30" ), igraph::vcount(pg), igraph::ecount(pg), graph_components, finite_fraction, cell_level_time_spearman, sample_level_time_spearman, time_monotonic, root_sensitivity_spearman, nrow( root_vertex_counts ), ROOT_ANCHOR_VERTEX, length( branch_vertices ), length( leaf_vertices ), ifelse( is.na( primary_branchpoint ), "None", primary_branchpoint ) ) )
wb <- openxlsx::createWorkbook()
add_sheet_safe( wb, "Parameters", parameter_tbl )
add_sheet_safe( wb, "Injury_input", injury_input_summary )
add_sheet_safe( wb, "Time_counts", time_counts )
add_sheet_safe( wb, "Root_candidate_counts", root_candidate_counts )
add_sheet_safe( wb, "Highconf_root_counts", highconf_root_counts )
add_sheet_safe( wb, "Root_candidate_metrics", root_candidate_tbl )
add_sheet_safe( wb, "Primary_root_cells", primary_root_tbl )
add_sheet_safe( wb, "Root_vertex_counts", root_vertex_counts )
add_sheet_safe( wb, "Time_pseudotime", time_summary )
add_sheet_safe( wb, "Sample_pseudotime", sample_pt_summary )
add_sheet_safe( wb, "Reference_states", reference_summary )
add_sheet_safe( wb, "Graph_vertices", graph_vertex_table )
add_sheet_safe( wb, "Branchpoint_ranking", branchpoint_summary )
add_sheet_safe( wb, "Arm_summary", arm_summary )
add_sheet_safe( wb, "Branch_program_7dpi", branch_program_summary )
add_sheet_safe( wb, "Branch_program_sample", branch_program_sample_summary )
add_sheet_safe( wb, "Cell_results", cell_df )
openxlsx::saveWorkbook( wb, OUT_XLSX, overwrite = TRUE )
saveRDS( cds, OUT_CDS, compress = TRUE )
saveRDS( obj_full, OUT_RDS, compress = TRUE )
writeLines( capture.output( sessionInfo() ), OUT_SESS, useBytes = TRUE )
