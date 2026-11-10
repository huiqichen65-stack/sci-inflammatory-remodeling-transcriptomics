options(stringsAsFactors = FALSE)
set.seed(20260731)
PROJECT_ROOT <- Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
INPUT_RDS <- file.path( PROJECT_ROOT, "04_results", "step5AR_GSE162610_consensus_homeostatic_root_refinement", "GSE162610_step5AR_consensus_root_trajectory_input.rds" )
OUTDIR <- file.path( PROJECT_ROOT, "04_results", "step5B_GSE162610_Monocle3_fate_bifurcation" )
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
OUT_FIG1 <- file.path(OUTDIR, "Fig5B_01_Monocle3_graph_topology.pdf")
OUT_FIG2 <- file.path(OUTDIR, "Fig5B_02_primary_pseudotime_validation.pdf")
OUT_FIG3 <- file.path(OUTDIR, "Fig5B_03_branchpoint_and_arm_assignment.pdf")
OUT_FIG4 <- file.path(OUTDIR, "Fig5B_04_Step4_program_projection.pdf")
OUT_FIG5 <- file.path(OUTDIR, "Fig5B_05_root_sensitivity.pdf")
OUT_XLSX <- file.path(OUTDIR, "Table5B_Monocle3_trajectory_results.xlsx")
OUT_CDS  <- file.path(OUTDIR, "GSE162610_step5B_monocle3_primary_cds.rds")
OUT_RDS  <- file.path(OUTDIR, "GSE162610_step5B_trajectory_scored_seurat.rds")
OUT_REP  <- file.path(OUTDIR, "Step5B_report.txt")
OUT_SESS <- file.path(OUTDIR, "Step5B_sessionInfo.txt")
TIME_COL       <- "timepoint"
SAMPLE_COL     <- "sample_id"
CLUSTER_COL    <- "step3A_cluster_r04"
UMAP_REDUCTION <- "umap.microglia.step3A"
PRIMARY_ROOT_COL <- "step5AR_v2_consensus_root"
SENS_ROOT_COL    <- "step5AR_v2_highconf_consensus"
REPAIR_SCORE_COL <- "repair_candidate_score"
SCAR_SCORE_COL   <- "scar_candidate_score"
NUM_DIM <- 50
MONOCLE_K <- 30
LEARN_GRAPH_EUCLIDEAN_DISTANCE_RATIO <- 0.8
LEARN_GRAPH_CLOSE_LOOP <- FALSE
MIN_TOTAL_CELLS_PER_DOWNSTREAM_ARM <- 100
MIN_7DPI_CELLS_PER_DOWNSTREAM_ARM <- 100
BRANCHPOINT_NEIGHBOR_GRAPH_RADIUS <- 2
required_pkgs <- c( "Seurat", "SeuratObject", "Matrix", "monocle3", "SingleCellExperiment", "SummarizedExperiment", "S4Vectors", "igraph", "dplyr", "tidyr", "tibble", "ggplot2", "patchwork", "openxlsx" )
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
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
})
TIME_COLORS <- c( "Uninjured" = "#66C2A5", "1 dpi"     = "#FC8D62", "3 dpi"     = "#8DA0CB", "7 dpi"     = "#E78AC3" )
CLUSTER_COLORS <- c( "0"  = "#4E79A7", "1"  = "#E15759", "2"  = "#59A14F", "3"  = "#F28E2B", "4"  = "#B07AA1", "5"  = "#76B7B2", "10" = "#2F4B7C" )
ROOT_COL <- "#D73027"
BG_COL   <- "#D9D9D9"
REPAIR_COLORS <- c("#F7FBFF", "#6BAED6", "#08306B")
SCAR_COLORS   <- c("#FFF5F0", "#FB6A4A", "#99000D")
DIFF_COLORS   <- c("#B2182B", "#F7F7F7", "#2166AC")
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
  stop( "Could not translate Monocle3 closest-vertex mapping to graph vertex names.", call. = FALSE )
}
cells_near_graph_vertex <- function( pg, closest_vertex, vertex_name, radius = 2) {
  dd <- igraph::distances( pg, v = vertex_name, to = igraph::V(pg), mode = "all" )
  near_vertices <- colnames(dd)[ as.numeric(dd[1, ]) <= radius ]
  names(closest_vertex)[ closest_vertex %in% near_vertices ]
}
summarize_branch_vertex <- function( bp, pg, closest_vertex, root_vertices, cell_df, radius = 2) {
  g2 <- igraph::delete_vertices( pg, bp )
  comp <- igraph::components(g2)
  membership <- comp$membership
  names(membership) <- igraph::V(g2)$name
  roots_remaining <- setdiff( root_vertices, bp )
  root_comp_ids <- unique( membership[ roots_remaining[ roots_remaining %in% names(membership) ] ] )
  root_comp_ids <- unique( root_comp_ids[ !is.na(root_comp_ids) ] )
  if (length(root_comp_ids) != 1) {
    return( tibble::tibble( branch_vertex = bp, degree = igraph::degree(pg, v = bp), root_side_components = length(root_comp_ids), downstream_components = NA_integer_, downstream_major_arms = NA_integer_, min_top2_7dpi_cells = NA_integer_, top2_7dpi_sum = NA_integer_, near_cells = NA_integer_, near_3dpi_cells = NA_integer_, near_transition_cells = NA_integer_, near_transition_fraction = NA_real_ ) )
  }
  root_comp <- root_comp_ids[1]
  cv <- closest_vertex[cell_df$cell]
  comp_for_cell <- membership[cv]
  tmp <- cell_df |> dplyr::mutate( graph_component = as.integer( comp_for_cell[cell] ) )
  downstream <- tmp |> dplyr::filter( !is.na(graph_component), graph_component != root_comp ) |> dplyr::group_by(graph_component) |> dplyr::summarise( n_cells = dplyr::n(), n_7dpi = sum(timepoint == "7 dpi"), n_3dpi = sum(timepoint == "3 dpi"), n_1dpi = sum(timepoint == "1 dpi"), n_uninj = sum(timepoint == "Uninjured"), n_transition = sum( shared_transition_3dpi_highconf, na.rm = TRUE ), median_pseudotime = median( pseudotime, na.rm = TRUE ), .groups = "drop" )
  major <- downstream |> dplyr::filter( n_cells >= MIN_TOTAL_CELLS_PER_DOWNSTREAM_ARM, n_7dpi >= MIN_7DPI_CELLS_PER_DOWNSTREAM_ARM ) |> dplyr::arrange( dplyr::desc(n_7dpi) )
  top2 <- head( major$n_7dpi, 2 )
  near_cells_vec <- cells_near_graph_vertex( pg, closest_vertex, bp, radius = radius )
  near_df <- cell_df |> dplyr::filter( cell %in% near_cells_vec )
  tibble::tibble( branch_vertex = bp, degree = igraph::degree(pg, v = bp), root_side_components = length(root_comp_ids), downstream_components = nrow(downstream), downstream_major_arms = nrow(major), min_top2_7dpi_cells = ifelse( length(top2) >= 2, min(top2), 0 ), top2_7dpi_sum = ifelse( length(top2) >= 2, sum(top2), sum(top2) ), near_cells = nrow(near_df), near_3dpi_cells = sum( near_df$timepoint == "3 dpi" ), near_transition_cells = sum( near_df$shared_transition_3dpi_highconf, na.rm = TRUE ), near_transition_fraction = ifelse( nrow(near_df) > 0, mean( near_df$shared_transition_3dpi_highconf, na.rm = TRUE ), NA_real_ ) )
}
if (!file.exists(INPUT_RDS)) {
  stop( paste0( "Required input not found: ", INPUT_RDS ), call. = FALSE )
}
obj <- readRDS(INPUT_RDS)
if (!inherits(obj, "Seurat")) {
  stop( "Step5A-R input is not a Seurat object.", call. = FALSE )
}
DefaultAssay(obj) <- "RNA"
required_meta <- c( TIME_COL, SAMPLE_COL, CLUSTER_COL, PRIMARY_ROOT_COL, SENS_ROOT_COL, "transition_status", "shared_transition_3dpi_highconf", REPAIR_SCORE_COL, SCAR_SCORE_COL )
missing_meta <- required_meta[ !required_meta %in% colnames(obj[[]]) ]
if (length(missing_meta) > 0) {
  stop( paste0( "Missing metadata: ", paste(missing_meta, collapse = ", ") ), call. = FALSE )
}
if (!UMAP_REDUCTION %in% Seurat::Reductions(obj)) {
  stop( paste0( "Required input not found: ", UMAP_REDUCTION ), call. = FALSE )
}
primary_root_cells <- colnames(obj)[ as.logical( obj[[PRIMARY_ROOT_COL, drop = TRUE]] ) ]
sens_root_cells <- colnames(obj)[ as.logical( obj[[SENS_ROOT_COL, drop = TRUE]] ) ]
if (length(primary_root_cells) != 72) {
  warning( paste0( "Primary root contains ", length(primary_root_cells), " cells rather than expected 72. ", "The script will continue using the stored frozen root metadata." ) )
}
if (length(primary_root_cells) < 30) {
  stop( "Primary root contains <30 cells.", call. = FALSE )
}
counts <- get_counts_matrix(obj)
if (!identical( colnames(counts), colnames(obj) )) {
  counts <- counts[ , colnames(obj), drop = FALSE ]
}
cell_metadata <- obj[[]]
cell_metadata <- cell_metadata[ colnames(counts), , drop = FALSE ]
gene_metadata <- data.frame( gene_short_name = rownames(counts), row.names = rownames(counts), stringsAsFactors = FALSE )
cds <- monocle3::new_cell_data_set( counts, cell_metadata = cell_metadata, gene_metadata = gene_metadata )
cds <- monocle3::preprocess_cds( cds, num_dim = NUM_DIM, method = "PCA" )
formal_umap <- Seurat::Embeddings( obj, reduction = UMAP_REDUCTION )
formal_umap <- formal_umap[ colnames(cds), 1:2, drop = FALSE ]
colnames(formal_umap) <- c( "UMAP_1", "UMAP_2" )
assign_ok <- tryCatch( {
    SingleCellExperiment::reducedDim( cds, "UMAP" ) <- formal_umap
    TRUE
  },
  error = function(e) {
    FALSE
  }
)
if (!assign_ok) {
  cds <- monocle3::reduce_dimension( cds, reduction_method = "UMAP", preprocess_method = "PCA" )
  cds@int_colData$reducedDims$UMAP <- formal_umap
}
imported_umap <- SingleCellExperiment::reducedDim( cds, "UMAP" )
if ( nrow(imported_umap) != nrow(formal_umap) ) {
  stop( "Imported UMAP dimensions match formal Seurat UMAP.", call. = FALSE )
}
umap_max_abs_diff <- max( abs( imported_umap - formal_umap ), na.rm = TRUE )
if (umap_max_abs_diff > 1e-8) {
  stop( "Monocle3 UMAP does not exactly match formal Seurat UMAP.", call. = FALSE )
}
cds <- monocle3::cluster_cells( cds, reduction_method = "UMAP", k = MONOCLE_K, cluster_method = "leiden", random_seed = 20260731, verbose = FALSE )
mono_cluster <- get_clusters_safe(cds)
mono_partition <- get_partitions_safe(cds)
colData(cds)$monocle3_cluster <- as.character(mono_cluster)
colData(cds)$monocle3_partition <- as.character(mono_partition)
partition_summary <- tibble::tibble( cell = colnames(cds), partition = as.character(mono_partition), timepoint = as.character( colData(cds)[[TIME_COL]] ), sample = as.character( colData(cds)[[SAMPLE_COL]] ) ) |> dplyr::count( partition, timepoint, name = "n_cells" )
cds <- monocle3::learn_graph( cds, use_partition = FALSE, close_loop = LEARN_GRAPH_CLOSE_LOOP, learn_graph_control = list( euclidean_distance_ratio = LEARN_GRAPH_EUCLIDEAN_DISTANCE_RATIO ) )
pg <- monocle3::principal_graph( cds )[["UMAP"]]
if (is.null(pg)) {
  stop( "Monocle3 principal graph was not generated.", call. = FALSE )
}
graph_components <- igraph::components(pg)$no
if (graph_components != 1) {
  warning( paste0( "Principal graph has ", graph_components, " connected components. ", "Interpretation requires extra caution." ) )
}
cds <- monocle3::order_cells( cds, reduction_method = "UMAP", root_cells = primary_root_cells )
pt_primary <- get_pseudotime_safe(cds)
pt_primary <- pt_primary[ colnames(cds) ]
finite_primary <- is.finite(pt_primary)
finite_fraction <- mean( finite_primary )
if (finite_fraction < 0.90) {
  warning( "Less than 90% of trajectory cells have finite primary pseudotime." )
}
colData(cds)$pseudotime_primary <- as.numeric(pt_primary)
cds_sens <- monocle3::order_cells( cds, reduction_method = "UMAP", root_cells = sens_root_cells )
pt_sens <- get_pseudotime_safe( cds_sens )
pt_sens <- pt_sens[ colnames(cds_sens) ]
both_finite <- is.finite(pt_primary) & is.finite(pt_sens)
pseudotime_root_sensitivity_spearman <- suppressWarnings( stats::cor( pt_primary[both_finite], pt_sens[both_finite], method = "spearman" ) )
colData(cds)$pseudotime_sensitivity <- as.numeric(pt_sens)
cell_df <- tibble::tibble( cell = colnames(cds), timepoint = as.character( colData(cds)[[TIME_COL]] ), sample = as.character( colData(cds)[[SAMPLE_COL]] ), step3A_cluster = as.character( colData(cds)[[CLUSTER_COL]] ), transition_status = as.character( colData(cds)$transition_status ), shared_transition_3dpi_highconf = as.logical( colData(cds)$shared_transition_3dpi_highconf ), primary_root = cell %in% primary_root_cells, sensitivity_root = cell %in% sens_root_cells, repair_candidate_score = as.numeric( colData(cds)[[REPAIR_SCORE_COL]] ), scar_candidate_score = as.numeric( colData(cds)[[SCAR_SCORE_COL]] ), pseudotime = as.numeric(pt_primary), pseudotime_sensitivity = as.numeric(pt_sens), monocle3_cluster = as.character(mono_cluster), monocle3_partition = as.character(mono_partition) )
cell_df$repair_minus_scar <- cell_df$repair_candidate_score - cell_df$scar_candidate_score
time_levels <- c( "Uninjured", "1 dpi", "3 dpi", "7 dpi" )
time_days <- c( "Uninjured" = 0, "1 dpi" = 1, "3 dpi" = 3, "7 dpi" = 7 )
time_summary <- cell_df |> dplyr::filter( is.finite(pseudotime) ) |> dplyr::mutate( timepoint = factor( timepoint, levels = time_levels, ordered = TRUE ) ) |> dplyr::group_by(timepoint) |> dplyr::summarise( n_cells = dplyr::n(), pseudotime_mean = mean(pseudotime), pseudotime_median = median(pseudotime), pseudotime_q25 = stats::quantile( pseudotime, 0.25 ), pseudotime_q75 = stats::quantile( pseudotime, 0.75 ), .groups = "drop" )
time_medians <- time_summary$pseudotime_median
names(time_medians) <- as.character(time_summary$timepoint)
time_monotonic <- all( diff( time_medians[ time_levels ] ) > 0 )
sample_pt_summary <- cell_df |> dplyr::filter( is.finite(pseudotime) ) |> dplyr::group_by( timepoint, sample ) |> dplyr::summarise( n_cells = dplyr::n(), median_pseudotime = median(pseudotime), mean_pseudotime = mean(pseudotime), .groups = "drop" ) |> dplyr::mutate( true_day = unname( time_days[timepoint] ) )
sample_level_time_spearman <- suppressWarnings( stats::cor( sample_pt_summary$true_day, sample_pt_summary$median_pseudotime, method = "spearman" ) )
cell_df$true_day <- unname( time_days[ cell_df$timepoint ] )
cell_level_time_spearman <- suppressWarnings( stats::cor( cell_df$true_day[ is.finite( cell_df$pseudotime ) ], cell_df$pseudotime[ is.finite( cell_df$pseudotime ) ], method = "spearman" ) )
reference_state_summary <- cell_df |> dplyr::mutate( reference_state = dplyr::case_when( primary_root ~ "Balanced_homeostatic_root", timepoint == "1 dpi" & step3A_cluster %in% c("2", "5", "10") ~ "Early_activated_reference", shared_transition_3dpi_highconf ~ "Shared_transition_3dpi", timepoint == "7 dpi" & step3A_cluster == "1" & transition_status == "Mature_DAM_like" ~ "Established_DAM_7dpi", TRUE ~ "Other" ) ) |> dplyr::filter( reference_state != "Other", is.finite(pseudotime) ) |> dplyr::group_by( reference_state ) |> dplyr::summarise( n_cells = dplyr::n(), pseudotime_mean = mean(pseudotime), pseudotime_median = median(pseudotime), pseudotime_q25 = stats::quantile( pseudotime, 0.25 ), pseudotime_q75 = stats::quantile( pseudotime, 0.75 ), .groups = "drop" )
closest_vertex <- get_closest_vertex( cds, pg )
closest_vertex <- closest_vertex[ cell_df$cell ]
cell_df$principal_vertex <- as.character( closest_vertex[ cell_df$cell ] )
root_vertices <- unique( closest_vertex[ primary_root_cells ] )
root_vertices <- unique( root_vertices[ !is.na(root_vertices) ] )
degree_vec <- igraph::degree( pg )
branch_vertices <- names( degree_vec[ degree_vec >= 3 ] )
leaf_vertices <- names( degree_vec[ degree_vec == 1 ] )
graph_vertex_table <- tibble::tibble( vertex = names(degree_vec), degree = as.integer(degree_vec), is_root_vertex = names(degree_vec) %in% root_vertices, is_branch_vertex = names(degree_vec) %in% branch_vertices, is_leaf_vertex = names(degree_vec) %in% leaf_vertices )
root_dist_mat <- igraph::distances( pg, v = igraph::V(pg), to = root_vertices )
root_dist_min <- apply( root_dist_mat, 1, min )
graph_vertex_table$graph_distance_from_root <- as.numeric( root_dist_min[ graph_vertex_table$vertex ] )
branchpoint_summary <- tibble::tibble()
if (length(branch_vertices) > 0) {
  branchpoint_summary <- lapply( branch_vertices, function(bp) {
      summarize_branch_vertex( bp = bp, pg = pg, closest_vertex = closest_vertex, root_vertices = root_vertices, cell_df = cell_df, radius = BRANCHPOINT_NEIGHBOR_GRAPH_RADIUS )
    }
  ) |> dplyr::bind_rows() |> dplyr::left_join(
      graph_vertex_table |> dplyr::select( vertex, graph_distance_from_root ), by = c( "branch_vertex" = "vertex" )
    ) |> dplyr::mutate(
      candidate_bifurcation = !is.na( downstream_major_arms ) & downstream_major_arms >= 2
    ) |> dplyr::arrange(
      dplyr::desc( candidate_bifurcation ), dplyr::desc( min_top2_7dpi_cells ), dplyr::desc( near_transition_cells ), dplyr::desc( top2_7dpi_sum )
    )
}
primary_branchpoint <- NA_character_
if ( nrow(branchpoint_summary) > 0 && any( branchpoint_summary$candidate_bifurcation, na.rm = TRUE ) ) {
  primary_branchpoint <- branchpoint_summary$branch_vertex[ which( branchpoint_summary$candidate_bifurcation )[1] ]
} else {
}
arm_summary <- tibble::tibble()
leaf_summary <- tibble::tibble()
branch_program_summary <- tibble::tibble()
branch_program_sample_summary <- tibble::tibble()
cell_df$branch_arm <- NA_character_
cell_df$major_late_arm <- FALSE
if (!is.na(primary_branchpoint)) {
  g2 <- igraph::delete_vertices( pg, primary_branchpoint )
  comp <- igraph::components(g2)
  membership <- comp$membership
  names(membership) <- igraph::V(g2)$name
  roots_remaining <- setdiff( root_vertices, primary_branchpoint )
  root_comp <- unique( membership[ roots_remaining[ roots_remaining %in% names(membership) ] ] )
  root_comp <- unique( root_comp[ !is.na(root_comp) ] )
  if (length(root_comp) != 1) {
    warning( "Primary branch point does not yield one unambiguous root-side component." )
  } else {
    root_comp <- root_comp[1]
    cv <- closest_vertex[ cell_df$cell ]
    comp_for_cell <- membership[ cv ]
    cell_df$branch_component <- as.integer( comp_for_cell[ cell_df$cell ] )
    at_bp <- cell_df$principal_vertex == primary_branchpoint
    arm_raw <- cell_df |> dplyr::filter( !is.na(branch_component), branch_component != root_comp ) |> dplyr::group_by( branch_component ) |> dplyr::summarise( n_cells = dplyr::n(), n_uninj = sum( timepoint == "Uninjured" ), n_1dpi = sum( timepoint == "1 dpi" ), n_3dpi = sum( timepoint == "3 dpi" ), n_7dpi = sum( timepoint == "7 dpi" ), n_transition = sum( shared_transition_3dpi_highconf, na.rm = TRUE ), pseudotime_median = median( pseudotime, na.rm = TRUE ), .groups = "drop" ) |> dplyr::arrange( dplyr::desc(n_7dpi), dplyr::desc(n_cells) ) |> dplyr::mutate( arm_rank = dplyr::row_number(), branch_arm = paste0( "Arm_", arm_rank ), major_late_arm = n_cells >= MIN_TOTAL_CELLS_PER_DOWNSTREAM_ARM & n_7dpi >= MIN_7DPI_CELLS_PER_DOWNSTREAM_ARM )
    arm_map <- stats::setNames( arm_raw$branch_arm, arm_raw$branch_component )
    major_map <- stats::setNames( arm_raw$major_late_arm, arm_raw$branch_component )
    cell_df$branch_arm[ !is.na( cell_df$branch_component ) & cell_df$branch_component == root_comp ] <- "Pre_branch_root_side"
    downstream_idx <- which( !is.na( cell_df$branch_component ) & cell_df$branch_component != root_comp )
    cell_df$branch_arm[ downstream_idx ] <- unname( arm_map[ as.character( cell_df$branch_component[ downstream_idx ] ) ] )
    cell_df$major_late_arm[ downstream_idx ] <- unname( major_map[ as.character( cell_df$branch_component[ downstream_idx ] ) ] )
    cell_df$branch_arm[ at_bp ] <- "Branchpoint_vertex"
    arm_summary <- arm_raw
    leaf_summary <- lapply( leaf_vertices, function(lf) {
        comp_id <- membership[lf]
        lf_cells <- names( closest_vertex )[ closest_vertex == lf ]
        tibble::tibble( leaf_vertex = lf, graph_component = ifelse( length(comp_id) == 0 || is.na(comp_id), NA_integer_, as.integer(comp_id) ), branch_arm = ifelse( is.na(comp_id), NA_character_, ifelse( comp_id == root_comp, "Pre_branch_root_side", unname( arm_map[ as.character( comp_id ) ] ) ) ), n_cells_mapped_exactly = length(lf_cells), n_7dpi_mapped_exactly = sum( cell_df$cell %in% lf_cells & cell_df$timepoint == "7 dpi" ) )
      }
    ) |> dplyr::bind_rows()
    branch_program_summary <- cell_df |> dplyr::filter( timepoint == "7 dpi", !is.na(branch_arm), grepl( "^Arm_", branch_arm ) ) |> dplyr::group_by( branch_arm ) |> dplyr::summarise( n_cells = dplyr::n(), repair_mean = mean( repair_candidate_score, na.rm = TRUE ), repair_median = median( repair_candidate_score, na.rm = TRUE ), scar_mean = mean( scar_candidate_score, na.rm = TRUE ), scar_median = median( scar_candidate_score, na.rm = TRUE ), repair_minus_scar_mean = mean( repair_minus_scar, na.rm = TRUE ), repair_minus_scar_median = median( repair_minus_scar, na.rm = TRUE ), pseudotime_median = median( pseudotime, na.rm = TRUE ), .groups = "drop" )
    branch_program_sample_summary <- cell_df |> dplyr::filter( timepoint == "7 dpi", !is.na(branch_arm), grepl( "^Arm_", branch_arm ) ) |> dplyr::group_by( sample, branch_arm ) |> dplyr::summarise( n_cells = dplyr::n(), repair_median = median( repair_candidate_score, na.rm = TRUE ), scar_median = median( scar_candidate_score, na.rm = TRUE ), repair_minus_scar_median = median( repair_minus_scar, na.rm = TRUE ), .groups = "drop" )
  }
}
colData(cds)$branch_arm <- cell_df$branch_arm
colData(cds)$major_late_arm <- cell_df$major_late_arm
obj$step5B_pseudotime_primary <- NA_real_
obj$step5B_pseudotime_sensitivity <- NA_real_
obj$step5B_monocle3_partition <- NA_character_
obj$step5B_monocle3_cluster <- NA_character_
obj$step5B_principal_vertex <- NA_character_
obj$step5B_branch_arm <- NA_character_
obj$step5B_major_late_arm <- FALSE
idx_obj <- match( cell_df$cell, colnames(obj) )
obj$step5B_pseudotime_primary[idx_obj] <- cell_df$pseudotime
obj$step5B_pseudotime_sensitivity[idx_obj] <- cell_df$pseudotime_sensitivity
obj$step5B_monocle3_partition[idx_obj] <- cell_df$monocle3_partition
obj$step5B_monocle3_cluster[idx_obj] <- cell_df$monocle3_cluster
obj$step5B_principal_vertex[idx_obj] <- cell_df$principal_vertex
obj$step5B_branch_arm[idx_obj] <- cell_df$branch_arm
obj$step5B_major_late_arm[idx_obj] <- cell_df$major_late_arm
obj$step5B_primary_branchpoint <- ifelse( is.na(primary_branchpoint), NA_character_, primary_branchpoint )
p1a <- monocle3::plot_cells( cds, reduction_method = "UMAP", color_cells_by = TIME_COL, show_trajectory_graph = FALSE, label_cell_groups = FALSE, label_leaves = FALSE, label_branch_points = FALSE, cell_size = 0.5 ) + ggplot2::scale_color_manual( values = TIME_COLORS ) + ggplot2::ggtitle( "True time on formal Step3A UMAP" )
p1b <- monocle3::plot_cells( cds, reduction_method = "UMAP", color_cells_by = CLUSTER_COL, show_trajectory_graph = FALSE, label_cell_groups = FALSE, label_leaves = FALSE, label_branch_points = FALSE, cell_size = 0.5 ) + ggplot2::scale_color_manual( values = CLUSTER_COLORS ) + ggplot2::ggtitle( "Formal Step3A states" )
p1c <- monocle3::plot_cells( cds, reduction_method = "UMAP", color_cells_by = "monocle3_partition", show_trajectory_graph = FALSE, label_cell_groups = FALSE, label_leaves = FALSE, label_branch_points = FALSE, cell_size = 0.5 ) + ggplot2::ggtitle( "Monocle3 partitions" )
p1d <- monocle3::plot_cells( cds, reduction_method = "UMAP", color_cells_by = TIME_COL, show_trajectory_graph = TRUE, label_cell_groups = FALSE, label_leaves = TRUE, label_branch_points = TRUE, graph_label_size = 3, cell_size = 0.4 ) + ggplot2::scale_color_manual( values = TIME_COLORS ) + ggplot2::ggtitle( "Principal graph learned independently" )
grDevices::pdf( OUT_FIG1, width = 14, height = 10, useDingbats = FALSE )
print( (p1a | p1b) / (p1c | p1d) + patchwork::plot_annotation( title = "Step 5B | Monocle3 graph topology" ) )
grDevices::dev.off()
p2a <- monocle3::plot_cells( cds, reduction_method = "UMAP", color_cells_by = "pseudotime", show_trajectory_graph = TRUE, label_cell_groups = FALSE, label_leaves = TRUE, label_branch_points = TRUE, graph_label_size = 3, cell_size = 0.45 ) + ggplot2::ggtitle( "Primary pseudotime | frozen 72-cell consensus root" )
pt_plot_df <- cell_df |> dplyr::filter( is.finite(pseudotime) ) |> dplyr::mutate( timepoint = factor( timepoint, levels = time_levels, ordered = TRUE ) )
p2b <- ggplot2::ggplot( pt_plot_df, ggplot2::aes( x = timepoint, y = pseudotime, fill = timepoint ) ) + ggplot2::geom_violin( trim = TRUE, scale = "width", alpha = 0.75, linewidth = 0.25 ) + ggplot2::geom_boxplot( width = 0.12, outlier.shape = NA, fill = "#FFFFFF", linewidth = 0.25 ) + ggplot2::scale_fill_manual( values = TIME_COLORS, guide = "none" ) + ggplot2::labs( title = "Pseudotime versus true time", subtitle = paste0( "Cell-level Spearman = ", round( cell_level_time_spearman, 3 ) ), x = "True time", y = "Primary pseudotime" ) + ggplot2::theme_classic(base_size = 11)
p2c <- ggplot2::ggplot( sample_pt_summary, ggplot2::aes( x = true_day, y = median_pseudotime, label = sample ) ) + ggplot2::geom_point( size = 2.7 ) + ggplot2::geom_smooth( method = "lm", se = FALSE, linewidth = 0.7 ) + ggplot2::geom_text( size = 3, vjust = -0.6, check_overlap = TRUE ) + ggplot2::labs( title = "Biological-sample median pseudotime", subtitle = paste0( "Sample-level Spearman = ", round( sample_level_time_spearman, 3 ) ), x = "True day after injury", y = "Median pseudotime" ) + ggplot2::theme_classic(base_size = 11)
ref_order <- c( "Balanced_homeostatic_root", "Early_activated_reference", "Shared_transition_3dpi", "Established_DAM_7dpi" )
ref_cell_plot <- cell_df |> dplyr::mutate( reference_state = dplyr::case_when( primary_root ~ "Balanced_homeostatic_root", timepoint == "1 dpi" & step3A_cluster %in% c("2", "5", "10") ~ "Early_activated_reference", shared_transition_3dpi_highconf ~ "Shared_transition_3dpi", timepoint == "7 dpi" & step3A_cluster == "1" & transition_status == "Mature_DAM_like" ~ "Established_DAM_7dpi", TRUE ~ NA_character_ ), reference_state = factor( reference_state, levels = ref_order, ordered = TRUE ) ) |> dplyr::filter( !is.na(reference_state), is.finite(pseudotime) )
p2d <- ggplot2::ggplot( ref_cell_plot, ggplot2::aes( x = reference_state, y = pseudotime ) ) + ggplot2::geom_violin( fill = "#9ECAE1", trim = TRUE, scale = "width", linewidth = 0.25 ) + ggplot2::geom_boxplot( width = 0.12, outlier.shape = NA, fill = "#FFFFFF", linewidth = 0.25 ) + ggplot2::labs( title = "Pseudotime across pre-defined Step3 biological references", x = NULL, y = "Primary pseudotime" ) + ggplot2::theme_classic(base_size = 10) + ggplot2::theme( axis.text.x = ggplot2::element_text( angle = 25, hjust = 1 ) )
grDevices::pdf( OUT_FIG2, width = 14, height = 10, useDingbats = FALSE )
print( (p2a | p2b) / (p2c | p2d) + patchwork::plot_annotation( title = "Step 5B | Primary pseudotime validation" ) )
grDevices::dev.off()
if (!is.na(primary_branchpoint)) {
  colData(cds)$branch_arm_display <- ifelse( is.na(cell_df$branch_arm), "Unassigned", cell_df$branch_arm )
  p3a <- monocle3::plot_cells( cds, reduction_method = "UMAP", color_cells_by = "branch_arm_display", show_trajectory_graph = TRUE, label_cell_groups = FALSE, label_leaves = TRUE, label_branch_points = TRUE, graph_label_size = 3, cell_size = 0.45 ) + ggplot2::ggtitle( paste0( "Graph arms after objective branch point ", primary_branchpoint ) )
  arm_time_plot <- cell_df |> dplyr::filter( !is.na(branch_arm), grepl( "^Arm_", branch_arm ) ) |> dplyr::count( branch_arm, timepoint, name = "n_cells" ) |> dplyr::group_by( branch_arm ) |> dplyr::mutate( fraction = n_cells / sum(n_cells) ) |> dplyr::ungroup()
  p3b <- ggplot2::ggplot( arm_time_plot, ggplot2::aes( x = branch_arm, y = fraction, fill = timepoint ) ) + ggplot2::geom_col( width = 0.70 ) + ggplot2::scale_fill_manual( values = TIME_COLORS, name = "True time" ) + ggplot2::labs( title = "True-time composition of downstream graph arms", x = "Graph arm", y = "Fraction of cells" ) + ggplot2::theme_classic(base_size = 11)
  near_primary <- cells_near_graph_vertex( pg, closest_vertex, primary_branchpoint, radius = BRANCHPOINT_NEIGHBOR_GRAPH_RADIUS )
  near_plot_df <- cell_df |> dplyr::filter( cell %in% near_primary )
  p3c <- ggplot2::ggplot( near_plot_df, ggplot2::aes( x = factor(timepoint), fill = shared_transition_3dpi_highconf ) ) + ggplot2::geom_bar( position = "stack" ) + ggplot2::scale_fill_manual( values = c( "FALSE" = "#D9D9D9", "TRUE" = "#D73027" ), name = "High-conf\n3 dpi transition" ) + ggplot2::labs( title = paste0( "Cells near ", primary_branchpoint, " (graph radius ", BRANCHPOINT_NEIGHBOR_GRAPH_RADIUS, ")" ), x = "True time", y = "Cells" ) + ggplot2::theme_classic(base_size = 11)
  bp_rank_plot <- branchpoint_summary |> dplyr::mutate( branch_vertex = factor( branch_vertex, levels = rev( branch_vertex ) ) )
  p3d <- ggplot2::ggplot( bp_rank_plot, ggplot2::aes( x = branch_vertex, y = min_top2_7dpi_cells, fill = candidate_bifurcation ) ) + ggplot2::geom_col( width = 0.70 ) + ggplot2::coord_flip() + ggplot2::scale_fill_manual( values = c( "FALSE" = "#BDBDBD", "TRUE" = "#2C7BB6" ), name = "Two major\nlate arms" ) + ggplot2::labs( title = "Objective branch-point ranking", subtitle = "Minimum 7 dpi support across the two largest downstream arms", x = "Principal-graph branch vertex", y = "Minimum top-two 7 dpi cells" ) + ggplot2::theme_classic(base_size = 10)
} else {
  p3a <- monocle3::plot_cells( cds, reduction_method = "UMAP", color_cells_by = "pseudotime", show_trajectory_graph = TRUE, label_cell_groups = FALSE, label_leaves = TRUE, label_branch_points = TRUE, cell_size = 0.45 ) + ggplot2::ggtitle( "No objective two-major-late-arm branch point detected" )
  placeholder <- tibble::tibble( x = 1, y = 1 )
  p3b <- ggplot2::ggplot( placeholder, ggplot2::aes(x, y) ) + ggplot2::annotate( "text", x = 1, y = 1, label = "No branch arm assignment was forced.", size = 5 ) + ggplot2::theme_void()
  p3c <- p3b
  p3d <- p3b
}
grDevices::pdf( OUT_FIG3, width = 14, height = 10, useDingbats = FALSE )
print( (p3a | p3b) / (p3c | p3d) + patchwork::plot_annotation( title = "Step 5B | Objective branch-point and arm discovery" ) )
grDevices::dev.off()
p4a <- monocle3::plot_cells( cds, reduction_method = "UMAP", color_cells_by = REPAIR_SCORE_COL, show_trajectory_graph = TRUE, label_cell_groups = FALSE, label_leaves = FALSE, label_branch_points = TRUE, graph_label_size = 3, cell_size = 0.45 ) + ggplot2::scale_color_gradientn( colors = REPAIR_COLORS ) + ggplot2::ggtitle( "Post-hoc repair-related candidate program" )
p4b <- monocle3::plot_cells( cds, reduction_method = "UMAP", color_cells_by = SCAR_SCORE_COL, show_trajectory_graph = TRUE, label_cell_groups = FALSE, label_leaves = FALSE, label_branch_points = TRUE, graph_label_size = 3, cell_size = 0.45 ) + ggplot2::scale_color_gradientn( colors = SCAR_COLORS ) + ggplot2::ggtitle( "Post-hoc scar-related candidate program" )
colData(cds)$repair_minus_scar <- cell_df$repair_minus_scar
p4c <- monocle3::plot_cells( cds, reduction_method = "UMAP", color_cells_by = "repair_minus_scar", show_trajectory_graph = TRUE, label_cell_groups = FALSE, label_leaves = FALSE, label_branch_points = TRUE, graph_label_size = 3, cell_size = 0.45 ) + ggplot2::scale_color_gradient2( low = DIFF_COLORS[1], mid = DIFF_COLORS[2], high = DIFF_COLORS[3], midpoint = 0 ) + ggplot2::ggtitle( "Post-hoc repair minus scar coordinate" )
if ( !is.na(primary_branchpoint) && nrow( branch_program_sample_summary ) > 0 ) {
  branch_long <- branch_program_sample_summary |> dplyr::select( sample, branch_arm, repair_median, scar_median ) |> tidyr::pivot_longer( cols = c( repair_median, scar_median ), names_to = "program", values_to = "sample_median_score" ) |> dplyr::mutate( program = dplyr::recode( program, repair_median = "Repair-related", scar_median = "Scar-related" ) )
  p4d <- ggplot2::ggplot( branch_long, ggplot2::aes( x = branch_arm, y = sample_median_score, shape = sample ) ) + ggplot2::geom_point( size = 3 ) + ggplot2::facet_wrap( ~program, nrow = 1, scales = "free_y" ) + ggplot2::labs( title = "7 dpi branch-arm program scores by biological sample", subtitle = "Post-hoc only; branch topology was learned without these scores", x = "Graph arm", y = "Sample-level median score" ) + ggplot2::theme_classic(base_size = 10) + ggplot2::theme( axis.text.x = ggplot2::element_text( angle = 25, hjust = 1 ) )
} else {
  p4d <- ggplot2::ggplot( tibble::tibble( x = 1, y = 1 ), ggplot2::aes(x, y) ) + ggplot2::annotate( "text", x = 1, y = 1, label = "No branch-arm program comparison:\nno objective bifurcation was forced.", size = 4.5 ) + ggplot2::theme_void()
}
grDevices::pdf( OUT_FIG4, width = 14, height = 10, useDingbats = FALSE )
print( (p4a | p4b) / (p4c | p4d) + patchwork::plot_annotation( title = "Step 5B | Step4 programs projected only after independent graph learning" ) )
grDevices::dev.off()
sens_df <- cell_df |> dplyr::filter( is.finite(pseudotime), is.finite( pseudotime_sensitivity ) )
p5a <- ggplot2::ggplot( sens_df, ggplot2::aes( x = pseudotime, y = pseudotime_sensitivity ) ) + ggplot2::geom_point( size = 0.45, alpha = 0.30 ) + ggplot2::labs( title = "Primary vs broad-root pseudotime", subtitle = paste0( "Spearman = ", round( pseudotime_root_sensitivity_spearman, 4 ) ), x = "Primary 72-cell-root pseudotime", y = "Broad high-confidence-root pseudotime" ) + ggplot2::theme_classic(base_size = 11)
sens_time <- sens_df |> dplyr::group_by( timepoint ) |> dplyr::summarise( primary_median = median(pseudotime), sensitivity_median = median( pseudotime_sensitivity ), .groups = "drop" ) |> tidyr::pivot_longer( cols = c( primary_median, sensitivity_median ), names_to = "root_definition", values_to = "median_pseudotime" ) |> dplyr::mutate( timepoint = factor( timepoint, levels = time_levels, ordered = TRUE ), root_definition = dplyr::recode( root_definition, primary_median = "Primary 72-cell root", sensitivity_median = "Broad high-confidence root" ) )
p5b <- ggplot2::ggplot( sens_time, ggplot2::aes( x = timepoint, y = median_pseudotime, group = root_definition, linetype = root_definition ) ) + ggplot2::geom_line( linewidth = 0.8 ) + ggplot2::geom_point( size = 2.5 ) + ggplot2::labs( title = "True-time ordering under two root definitions", x = "True time", y = "Median pseudotime", linetype = "Root definition" ) + ggplot2::theme_classic(base_size = 11)
root_umap_df <- tibble::tibble( cell = rownames(formal_umap), UMAP_1 = formal_umap[, 1], UMAP_2 = formal_umap[, 2], primary_root = rownames(formal_umap) %in% primary_root_cells, sensitivity_root = rownames(formal_umap) %in% sens_root_cells )
p5c <- ggplot2::ggplot( root_umap_df, ggplot2::aes( UMAP_1, UMAP_2 ) ) + ggplot2::geom_point( color = BG_COL, size = 0.35, alpha = 0.25 ) + ggplot2::geom_point( data = root_umap_df |> dplyr::filter( sensitivity_root ), color = "#2C7BB6", size = 0.70, alpha = 0.55 ) + ggplot2::geom_point( data = root_umap_df |> dplyr::filter( primary_root ), color = ROOT_COL, size = 0.95, alpha = 0.95 ) + ggplot2::labs( title = "Root sensitivity sets", subtitle = "Red = primary 72; blue = broader high-confidence root" ) + ggplot2::theme_classic(base_size = 11)
pt_diff <- sens_df |> dplyr::mutate( primary_scaled = (pseudotime - min(pseudotime)) / ( max(pseudotime) - min(pseudotime) ), sensitivity_scaled = ( pseudotime_sensitivity - min( pseudotime_sensitivity ) ) / ( max( pseudotime_sensitivity ) - min( pseudotime_sensitivity ) ), scaled_difference = sensitivity_scaled - primary_scaled )
p5d <- ggplot2::ggplot( pt_diff, ggplot2::aes( x = scaled_difference ) ) + ggplot2::geom_histogram( bins = 60 ) + ggplot2::labs( title = "Normalized pseudotime difference", x = "Sensitivity - primary scaled pseudotime", y = "Cells" ) + ggplot2::theme_classic(base_size = 11)
grDevices::pdf( OUT_FIG5, width = 14, height = 10, useDingbats = FALSE )
print( (p5a | p5b) / (p5c | p5d) + patchwork::plot_annotation( title = "Step 5B | Root-orientation sensitivity analysis" ) )
grDevices::dev.off()
parameter_tbl <- tibble::tibble( parameter = c( "Trajectory input cells", "Primary root cells", "Sensitivity root cells", "Formal UMAP", "PCA dimensions", "Monocle k", "learn_graph use_partition", "learn_graph close_loop", "euclidean_distance_ratio", "Principal graph vertices", "Principal graph edges", "Principal graph connected components", "Finite pseudotime fraction", "Cell-level true-time Spearman", "Sample-level true-time Spearman", "True-time median monotonic", "Primary-vs-sensitivity pseudotime Spearman", "Branch vertices", "Leaf vertices", "Primary objective branch point" ), value = c( ncol(obj), length(primary_root_cells), length(sens_root_cells), UMAP_REDUCTION, NUM_DIM, MONOCLE_K, FALSE, LEARN_GRAPH_CLOSE_LOOP, LEARN_GRAPH_EUCLIDEAN_DISTANCE_RATIO, igraph::vcount(pg), igraph::ecount(pg), graph_components, finite_fraction, cell_level_time_spearman, sample_level_time_spearman, time_monotonic, pseudotime_root_sensitivity_spearman, length(branch_vertices), length(leaf_vertices), ifelse( is.na( primary_branchpoint ), "None", primary_branchpoint ) ) )
root_cells_tbl <- tibble::tibble( cell = primary_root_cells, sample = as.character( colData(cds)[ primary_root_cells, SAMPLE_COL ] ), step3A_cluster = as.character( colData(cds)[ primary_root_cells, CLUSTER_COL ] ), pseudotime = pt_primary[ primary_root_cells ] )
wb <- openxlsx::createWorkbook()
add_sheet <- function(name, df) {
  openxlsx::addWorksheet( wb, sheetName = name )
  if ( is.null(df) || ncol(as.data.frame(df)) == 0 || nrow(as.data.frame(df)) == 0 ) {
    openxlsx::writeData( wb, sheet = name, x = data.frame( note = "No results generated for this table under the predefined Step 5B criteria.", stringsAsFactors = FALSE ) )
  } else {
    openxlsx::writeDataTable( wb, sheet = name, x = as.data.frame(df), tableStyle = "TableStyleMedium2" )
  }
  openxlsx::freezePane( wb, sheet = name, firstRow = TRUE )
}
add_sheet("Parameters", parameter_tbl)
add_sheet("Primary_root_cells", root_cells_tbl)
add_sheet("Time_pseudotime", time_summary)
add_sheet("Sample_pseudotime", sample_pt_summary)
add_sheet("Reference_states", reference_state_summary)
add_sheet("Partitions", partition_summary)
add_sheet("Graph_vertices", graph_vertex_table)
add_sheet("Branchpoint_ranking", branchpoint_summary)
add_sheet("Arm_summary", arm_summary)
add_sheet("Leaf_summary", leaf_summary)
add_sheet("Branch_program_7dpi", branch_program_summary)
add_sheet( "Branch_program_sample", branch_program_sample_summary )
add_sheet("Cell_results", cell_df)
openxlsx::saveWorkbook( wb, OUT_XLSX, overwrite = TRUE )
saveRDS( cds, OUT_CDS, compress = TRUE )
saveRDS( obj, OUT_RDS, compress = TRUE )
writeLines( capture.output( sessionInfo() ), OUT_SESS, useBytes = TRUE )
