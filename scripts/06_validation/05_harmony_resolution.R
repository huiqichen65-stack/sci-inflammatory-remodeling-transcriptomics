options( stringsAsFactors = FALSE, warn = 1, Seurat.object.assay.version = "v5", future.globals.maxSize = 16 * 1024^3 )
SEED <- 20260803L
PCS_USE <- 1:40
RESOLUTIONS <- c(0.1, 0.2, 0.3, 0.4)
DISPLAY_RESOLUTION <- 0.2
K_GRAPH <- 30L
K_MIXING <- 30L
N_TREES <- 50L
SNN_PRUNE <- 1 / 15
UMAP_N_NEIGHBORS <- 30L
UMAP_MIN_DIST <- 0.3
HCON_THETA <- 1
HCON_LAMBDA <- 2
HDEF_THETA <- 2
HDEF_LAMBDA <- 1
HARMONY_MAX_ITER <- 20L
HARMONY_NCORES <- 4L
LIBRARY_DOMINANCE_THRESHOLD <- 0.80
TIMEPOINT_DOMINANCE_THRESHOLD <- 0.95
MIN_CLUSTER_SIZE_FOR_DOMINANCE <- 500L
MIN_ACCEPTABLE_MEAN_KNN_OVERLAP <- 0.30
expected_timepoints <- c( "Uninjured", "1 dpi", "4 dpi", "7 dpi", "14 dpi", "1 month", "2 months" )
branch_specs <- data.frame( branch = c( "pca", "harmony_conservative", "harmony_default" ), reduction = c( "pca", "harmony_conservative", "harmony_default" ), graph_prefix = c( "cmp_pca", "cmp_hcon", "cmp_hdef" ), cluster_prefix = c( "pca", "hcon", "hdef" ), umap = c( "umap_pca_unintegrated", "umap_harmony_conservative", "umap_harmony_default" ), stringsAsFactors = FALSE )
project_root <- Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
input_object_path <- file.path( project_root, "04_results", "step10B_GSE234774_normalization_PCA_UMAP_clustering", "RDS", "GSE234774_step10B_LogNormalize_PCA_UMAP_clusters_BPCells.rds" )
output_dir <- file.path( project_root, "04_results", "step10C_GSE234774_Harmony_resolution_sensitivity" )
tables_dir <- file.path(output_dir, "tables")
figures_dir <- file.path(output_dir, "figures")
rds_dir <- file.path(output_dir, "RDS")
for (d in c( output_dir, tables_dir, figures_dir, rds_dir )) {
  dir.create( d, recursive = TRUE, showWarnings = FALSE )
}
run_log <- file.path( output_dir, "00_run_log.txt" )
if (file.exists(run_log)) {
  file.remove(run_log)
}
output_object_path <- file.path( rds_dir, "GSE234774_step10C_PCA_Harmony_resolution_comparison_BPCells.rds" )
log_message <- function(...) {
  message_text <- paste0(...)
  timestamp <- format( Sys.time(), "%Y-%m-%d %H:%M:%S" )
  line <- paste0( "[", timestamp, "] ", message_text )
  cat(line, "\n")
  cat( line, "\n", file = run_log, append = TRUE )
}
write_csv_safe <- function( x, filename ) {
  path <- file.path( tables_dir, filename )
  utils::write.csv( x, path, row.names = FALSE, na = "" )
  log_message( "Wrote: ", path )
  invisible(path)
}
stop_with_report <- function(reason) {
  report_path <- file.path( output_dir, "10_STEP10C_FINAL_REPORT.txt" )
  report <- c( "GSE234774 STEP 10C HARMONY / RESOLUTION SENSITIVITY", "============================================================", paste0( "Run time: ", format( Sys.time(), "%Y-%m-%d %H:%M:%S" ) ), paste0( "Output directory: ", output_dir ), "", "FINAL DECISION: STOP", "", paste0( "Reason: ", reason ), "", "Do not freeze a reduction or clustering resolution." )
  writeLines( report, report_path, useBytes = TRUE )
  stop( reason, call. = FALSE )
}
safe_package_version <- function(pkg) {
  if (!requireNamespace( pkg, quietly = TRUE )) {
    return(NA_character_)
  }
  as.character( utils::packageVersion(pkg) )
}
resolution_tag <- function(resolution) {
  gsub( "\\.", "_", format( resolution, scientific = FALSE, trim = TRUE ) )
}
cluster_column_name <- function( prefix, resolution ) {
  paste0( prefix, "_res_", resolution_tag(resolution) )
}
adjusted_rand_index <- function( x, y ) {
  x <- as.character(x)
  y <- as.character(y)
  ok <- !is.na(x) & !is.na(y)
  x <- x[ok]
  y <- y[ok]
  contingency <- table(x, y)
  n <- sum(contingency)
  if (n < 2L) {
    return(NA_real_)
  }
  choose2 <- function(z) {
    z * (z - 1) / 2
  }
  sum_cells <- sum( choose2(contingency) )
  sum_rows <- sum( choose2( rowSums(contingency) ) )
  sum_cols <- sum( choose2( colSums(contingency) ) )
  total_pairs <- choose2(n)
  expected_index <- ( sum_rows * sum_cols ) / total_pairs
  maximum_index <- 0.5 * ( sum_rows + sum_cols )
  denominator <- maximum_index - expected_index
  if (abs(denominator) < .Machine$double.eps) {
    return(1)
  }
  ( sum_cells - expected_index ) / denominator
}
clean_neighbor_indices <- function( index_matrix, k ) {
  n <- nrow(index_matrix)
  cleaned <- matrix( NA_integer_, nrow = n, ncol = k )
  for (i in seq_len(n)) {
    values <- index_matrix[i, ]
    values <- values[ !is.na(values) & values != i ]
    take <- min( k, length(values) )
    if (take > 0L) {
      cleaned[ i, seq_len(take) ] <- values[ seq_len(take) ]
    }
  }
  cleaned
}
calculate_within_timepoint_mixing <- function( embeddings, metadata, reduction_name, k = 30L ) {
  results <- list()
  neighbor_cache <- list()
  for (timepoint_name in expected_timepoints) {
    cells <- rownames(metadata)[ as.character( metadata$timepoint ) == timepoint_name ]
    sub_meta <- metadata[ cells, , drop = FALSE ]
    sub_embeddings <- embeddings[ cells, , drop = FALSE ]
    libraries <- as.character( sub_meta$biological_library )
    library_levels <- unique(libraries)
    if ( length(library_levels) < 2L || nrow(sub_embeddings) <= k + 1L ) {
      stop( paste0( "Insufficient libraries/cells for mixing metric at ", timepoint_name, "." ) )
    }
    log_message( "Finding within-timepoint neighbors: ", reduction_name, " / ", timepoint_name, "." )
    neighbor_object <- Seurat::FindNeighbors( object = sub_embeddings, k.param = min( k + 1L, nrow(sub_embeddings) - 1L ), return.neighbor = TRUE, compute.SNN = FALSE, nn.method = "annoy", n.trees = N_TREES, annoy.metric = "euclidean", verbose = FALSE )
    raw_indices <- SeuratObject::Indices( neighbor_object )
    clean_indices <- clean_neighbor_indices( raw_indices, k = k )
    neighbor_libraries <- matrix( libraries[ clean_indices ], nrow = nrow(clean_indices), ncol = ncol(clean_indices) )
    valid_neighbor_counts <- rowSums( !is.na( neighbor_libraries ) )
    same_library_matrix <- sweep( neighbor_libraries, MARGIN = 1L, STATS = libraries, FUN = "==" )
    same_library_fraction <- rowSums( same_library_matrix, na.rm = TRUE ) / valid_neighbor_counts
    other_library_fraction <- 1 - same_library_fraction
    has_other_library_neighbor <- rowSums( !same_library_matrix & !is.na( same_library_matrix ), na.rm = TRUE ) > 0L
    probability_matrix <- vapply( library_levels, function(library_name) {
        rowSums( neighbor_libraries == library_name, na.rm = TRUE ) / valid_neighbor_counts
      },
      numeric( nrow(neighbor_libraries) )
    )
    if (is.null( dim(probability_matrix) )) {
      probability_matrix <- matrix( probability_matrix, ncol = 1L )
    }
    entropy_terms <- ifelse( probability_matrix > 0, probability_matrix * log( probability_matrix ), 0 )
    normalized_entropy <- -rowSums( entropy_terms ) / log( length( library_levels ) )
    cell_metrics <- data.frame( barcode = cells, reduction = reduction_name, timepoint = timepoint_name, biological_library = libraries, other_library_fraction = other_library_fraction, normalized_library_entropy = normalized_entropy, has_other_library_neighbor = has_other_library_neighbor, valid_neighbor_count = valid_neighbor_counts, stringsAsFactors = FALSE )
    summary_metrics <- data.table::as.data.table( cell_metrics )[ , .( n_cells = .N, n_libraries = uniqueN( biological_library ), mean_other_library_fraction = mean( other_library_fraction, na.rm = TRUE ), median_other_library_fraction = median( other_library_fraction, na.rm = TRUE ), mean_normalized_library_entropy = mean( normalized_library_entropy, na.rm = TRUE ), median_normalized_library_entropy = median( normalized_library_entropy, na.rm = TRUE ), fraction_with_other_library_neighbor = mean( has_other_library_neighbor, na.rm = TRUE ) ), by = .( reduction, timepoint ) ]
    results[[timepoint_name]] <- list( cell_metrics = cell_metrics, summary_metrics = as.data.frame( summary_metrics ) )
    neighbor_cache[[timepoint_name]] <- list( cells = cells, indices = clean_indices )
    rm( neighbor_object, raw_indices, clean_indices, neighbor_libraries, probability_matrix )
    invisible(gc())
  }
  list( summaries = do.call( rbind, lapply( results, `[[`, "summary_metrics" ) ), cell_metrics = do.call( rbind, lapply( results, `[[`, "cell_metrics" ) ), neighbor_cache = neighbor_cache )
}
calculate_neighbor_overlap <- function( reference_cache, comparison_cache, comparison_name ) {
  rows <- list()
  for (timepoint_name in expected_timepoints) {
    ref <- reference_cache[[timepoint_name]]
    cmp <- comparison_cache[[timepoint_name]]
    if (!identical( ref$cells, cmp$cells )) {
      stop( paste0( "Cell order differs in neighbor caches at ", timepoint_name, "." ) )
    }
    ref_idx <- ref$indices
    cmp_idx <- cmp$indices
    overlap_fraction <- vapply( seq_len( nrow(ref_idx) ), function(i) {
        ref_values <- ref_idx[i, ]
        cmp_values <- cmp_idx[i, ]
        ref_values <- ref_values[ !is.na(ref_values) ]
        cmp_values <- cmp_values[ !is.na(cmp_values) ]
        denominator <- min( length(ref_values), length(cmp_values) )
        if (denominator == 0L) {
          return(NA_real_)
        }
        length( intersect( ref_values, cmp_values ) ) / denominator
      },
      numeric(1)
    )
    rows[[timepoint_name]] <- data.frame( comparison = comparison_name, timepoint = timepoint_name, n_cells = length( overlap_fraction ), mean_KNN_overlap = mean( overlap_fraction, na.rm = TRUE ), median_KNN_overlap = median( overlap_fraction, na.rm = TRUE ), q10_KNN_overlap = as.numeric( stats::quantile( overlap_fraction, 0.10, na.rm = TRUE ) ), q90_KNN_overlap = as.numeric( stats::quantile( overlap_fraction, 0.90, na.rm = TRUE ) ), stringsAsFactors = FALSE )
  }
  do.call( rbind, rows )
}
run_harmony_matrix_compat <- function( pca_embeddings, metadata, theta, lambda, branch_label ) {
  harmony_version <- utils::packageVersion( "harmony" )
  log_message( "Running ", branch_label, " with theta=", theta, ", lambda=", lambda, ", harmony version ", harmony_version, "." )
  set.seed(SEED)
  if ( harmony_version >= "2.0.0" ) {
    output <- harmony::RunHarmony( data_mat = pca_embeddings, meta_data = metadata, vars_use = "biological_library", theta = theta, lambda = lambda, sigma = 0.1, nclust = NULL, max_iter = HARMONY_MAX_ITER, early_stop = TRUE, ncores = HARMONY_NCORES, plot_convergence = FALSE, return_object = FALSE, verbose = TRUE )
  } else if (
    exists( "HarmonyMatrix", where = asNamespace( "harmony" ), inherits = FALSE )
  ) {
    output <- harmony::HarmonyMatrix( data_mat = pca_embeddings, meta_data = metadata, vars_use = "biological_library", do_pca = FALSE, theta = theta, lambda = lambda, sigma = 0.1, nclust = NULL, tau = 0, block.size = 0.05, max.iter.harmony = HARMONY_MAX_ITER, plot_convergence = FALSE, return_object = FALSE, verbose = TRUE )
  } else {
    stop( paste0( "Unsupported harmony version: ", harmony_version, ". Update harmony from CRAN." ) )
  }
  output <- as.matrix( output )
  if ( nrow(output) != nrow(pca_embeddings) ) {
    if ( ncol(output) == nrow(pca_embeddings) ) {
      output <- t(output)
    } else {
      stop( paste0( branch_label, " returned unexpected dimensions." ) )
    }
  }
  rownames(output) <- rownames( pca_embeddings )
  output
}
summarize_clustering <- function( metadata, cluster_column, branch, resolution ) {
  clusters <- as.character( metadata[ , cluster_column ] )
  cluster_dt <- data.table::data.table( cluster = clusters, biological_library = as.character( metadata$biological_library ), timepoint = as.character( metadata$timepoint ) )
  sizes <- cluster_dt[ , .( n_cells = .N, n_libraries = uniqueN( biological_library ), n_timepoints = uniqueN( timepoint ) ), by = cluster ]
  library_composition <- cluster_dt[ , .( n_cells = .N ), by = .( cluster, biological_library ) ]
  library_composition <- merge( library_composition, sizes[ , c( "cluster", "n_cells" ), with = FALSE ], by = "cluster", suffixes = c( "_library", "_cluster" ), all.x = TRUE, sort = FALSE )
  library_composition$fraction <- library_composition$n_cells_library / library_composition$n_cells_cluster
  max_library <- library_composition[ order( cluster, -fraction ), .SD[1L], by = cluster ]
  timepoint_composition <- cluster_dt[ , .( n_cells = .N ), by = .( cluster, timepoint ) ]
  timepoint_composition <- merge( timepoint_composition, sizes[ , c( "cluster", "n_cells" ), with = FALSE ], by = "cluster", suffixes = c( "_timepoint", "_cluster" ), all.x = TRUE, sort = FALSE )
  timepoint_composition$fraction <- timepoint_composition$n_cells_timepoint / timepoint_composition$n_cells_cluster
  max_timepoint <- timepoint_composition[ order( cluster, -fraction ), .SD[1L], by = cluster ]
  sizes <- merge( sizes, max_library[ , .( cluster, dominant_library = biological_library, maximum_library_fraction = fraction ) ], by = "cluster", all.x = TRUE, sort = FALSE )
  sizes <- merge( sizes, max_timepoint[ , .( cluster, dominant_timepoint = timepoint, maximum_timepoint_fraction = fraction ) ], by = "cluster", all.x = TRUE, sort = FALSE )
  sizes$branch <- branch
  sizes$resolution <- resolution
  sizes$cluster_column <- cluster_column
  sizes$library_dominance_review <- sizes$n_cells >= MIN_CLUSTER_SIZE_FOR_DOMINANCE & sizes$maximum_library_fraction >= LIBRARY_DOMINANCE_THRESHOLD
  sizes$timepoint_dominance_review <- sizes$n_cells >= MIN_CLUSTER_SIZE_FOR_DOMINANCE & sizes$maximum_timepoint_fraction >= TIMEPOINT_DOMINANCE_THRESHOLD
  list( cluster_detail = as.data.frame( sizes ), summary = data.frame( branch = branch, resolution = resolution, cluster_column = cluster_column, n_clusters = nrow(sizes), minimum_cluster_size = min( sizes$n_cells ), median_cluster_size = median( sizes$n_cells ), maximum_cluster_size = max( sizes$n_cells ), n_clusters_lt_100 = sum( sizes$n_cells < 100L ), n_clusters_lt_500 = sum( sizes$n_cells < 500L ), n_library_dominant_clusters = sum( sizes$library_dominance_review ), n_timepoint_dominant_clusters = sum( sizes$timepoint_dominance_review ), stringsAsFactors = FALSE ) )
}
log_message( "Starting package preflight." )
required_packages <- c( "Matrix", "data.table", "ggplot2", "patchwork", "scales", "Seurat", "SeuratObject", "harmony" )
missing_packages <- required_packages[ !vapply( required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1) ) ]
if ( "harmony" %in% missing_packages ) {
  log_message( "Installing harmony from CRAN." )
  try( utils::install.packages( "harmony", repos = "https://cloud.r-project.org", dependencies = TRUE ), silent = TRUE )
}
missing_after_install <- required_packages[ !vapply( required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1) ) ]
if ( length( missing_after_install ) > 0L ) {
  stop_with_report( paste0( "Missing required package(s): ", paste( missing_after_install, collapse = ", " ), "." ) )
}
suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(Seurat)
  library(SeuratObject)
  library(harmony)
})
package_versions <- data.frame( package = c( "R", "Matrix", "data.table", "ggplot2", "patchwork", "scales", "Seurat", "SeuratObject", "harmony" ), version = c( paste0( R.version$major, ".", R.version$minor ), safe_package_version("Matrix"), safe_package_version("data.table"), safe_package_version("ggplot2"), safe_package_version("patchwork"), safe_package_version("scales"), safe_package_version("Seurat"), safe_package_version("SeuratObject"), safe_package_version("harmony") ), stringsAsFactors = FALSE )
write_csv_safe( package_versions, "00_package_versions.csv" )
if ( !file.exists( input_object_path ) ) {
  stop_with_report( paste0( "Input object not found: ", input_object_path ) )
}
log_message( "Loading Step 10B object." )
obj <- readRDS( input_object_path )
if ( !inherits( obj, "Seurat" ) ) {
  stop_with_report( "Input RDS is not a Seurat object." )
}
if ( nrow(obj) != 26918L || ncol(obj) != 144359L ) {
  stop_with_report( paste0( "Expected 26,918 features x 144,359 singlets; detected ", nrow(obj), " x ", ncol(obj), "." ) )
}
required_metadata <- c( "biological_library", "timepoint", "replicate_id", "nCount_RNA", "nFeature_RNA", "percent.mt", "scDblFinder.class" )
missing_metadata <- setdiff( required_metadata, names( obj@meta.data ) )
if ( length( missing_metadata ) > 0L ) {
  stop_with_report( paste0( "Missing metadata: ", paste( missing_metadata, collapse = ", " ), "." ) )
}
if ( !all( as.character( obj$scDblFinder.class ) == "singlet" ) ) {
  stop_with_report( "Input contains non-singlet nuclei." )
}
if ( length( unique( as.character( obj$biological_library ) ) ) != 20L ) {
  stop_with_report( "Input does not contain exactly 20 biological libraries." )
}
if ( !all( expected_timepoints %in% unique( as.character( obj$timepoint ) ) ) ) {
  stop_with_report( "Expected timepoints are incomplete." )
}
if ( !"pca" %in% names( obj@reductions ) ) {
  stop_with_report( "PCA reduction is absent." )
}
pca_embeddings <- Seurat::Embeddings( object = obj, reduction = "pca" )[ , PCS_USE, drop = FALSE ]
if ( nrow(pca_embeddings) != ncol(obj) || ncol(pca_embeddings) != length(PCS_USE) ) {
  stop_with_report( "PCA embedding dimensions are incorrect." )
}
metadata <- obj@meta.data
metadata$timepoint <- factor( as.character( metadata$timepoint ), levels = expected_timepoints, ordered = TRUE )
metadata$biological_library <- factor( as.character( metadata$biological_library ) )
rownames(metadata) <- colnames(obj)
if ( !identical( rownames(pca_embeddings), rownames(metadata) ) ) {
  stop_with_report( "PCA embeddings and metadata are not aligned." )
}
pca_copy_for_integrity <- pca_embeddings
cell_order_before <- colnames(obj)
input_audit <- data.frame( item = c( "Input object", "Features", "Singlet nuclei", "Biological libraries", "Timepoints", "PCA dimensions used", "Existing unintegrated UMAP", "RNA assay retained" ), value = c( normalizePath( input_object_path, winslash = "/", mustWork = TRUE ), nrow(obj), ncol(obj), length( unique( metadata$biological_library ) ), length( unique( metadata$timepoint ) ), paste0( min(PCS_USE), "-", max(PCS_USE) ), "umap_pca_unintegrated" %in% names( obj@reductions ), "RNA" %in% names( obj@assays ) ), stringsAsFactors = FALSE )
write_csv_safe( input_audit, "01_input_object_audit.csv" )
harmony_conservative_embeddings <- run_harmony_matrix_compat( pca_embeddings = pca_embeddings, metadata = metadata, theta = HCON_THETA, lambda = HCON_LAMBDA, branch_label = "conservative Harmony" )
harmony_default_embeddings <- run_harmony_matrix_compat( pca_embeddings = pca_embeddings, metadata = metadata, theta = HDEF_THETA, lambda = HDEF_LAMBDA, branch_label = "default-strength Harmony" )
if ( !identical( rownames( harmony_conservative_embeddings ), rownames( pca_embeddings ) ) || !identical( rownames( harmony_default_embeddings ), rownames( pca_embeddings ) ) ) {
  stop_with_report( "Harmony embedding cell names are not aligned." )
}
colnames( harmony_conservative_embeddings ) <- paste0( "HCON_", seq_len( ncol( harmony_conservative_embeddings ) ) )
colnames( harmony_default_embeddings ) <- paste0( "HDEF_", seq_len( ncol( harmony_default_embeddings ) ) )
obj[["harmony_conservative"]] <- SeuratObject::CreateDimReducObject( embeddings = harmony_conservative_embeddings, key = "HCON_", assay = "RNA" )
obj[["harmony_default"]] <- SeuratObject::CreateDimReducObject( embeddings = harmony_default_embeddings, key = "HDEF_", assay = "RNA" )
cluster_summary_rows <- list()
cluster_detail_rows <- list()
for ( branch_index in seq_len( nrow(branch_specs) ) ) {
  branch <- branch_specs$branch[ branch_index ]
  reduction_name <- branch_specs$reduction[ branch_index ]
  graph_prefix <- branch_specs$graph_prefix[ branch_index ]
  cluster_prefix <- branch_specs$cluster_prefix[ branch_index ]
  umap_name <- branch_specs$umap[ branch_index ]
  graph_names <- c( paste0( graph_prefix, "_nn" ), paste0( graph_prefix, "_snn" ) )
  log_message( "Building graph for ", branch, "." )
  obj <- Seurat::FindNeighbors( object = obj, reduction = reduction_name, dims = PCS_USE, k.param = K_GRAPH, compute.SNN = TRUE, prune.SNN = SNN_PRUNE, nn.method = "annoy", n.trees = N_TREES, annoy.metric = "euclidean", graph.name = graph_names, verbose = TRUE )
  for ( resolution_index in seq_along( RESOLUTIONS ) ) {
    resolution <- RESOLUTIONS[ resolution_index ]
    cluster_column <- cluster_column_name( cluster_prefix, resolution )
    log_message( "Clustering ", branch, " at resolution ", resolution, "." )
    obj <- Seurat::FindClusters( object = obj, graph.name = graph_names[[2L]], cluster.name = cluster_column, resolution = resolution, algorithm = 1, n.start = 10L, n.iter = 10L, random.seed = SEED + branch_index * 100L + resolution_index, group.singletons = TRUE, verbose = TRUE )
    clustering_summary <- summarize_clustering( metadata = obj@meta.data, cluster_column = cluster_column, branch = branch, resolution = resolution )
    cluster_summary_rows[[ length( cluster_summary_rows ) + 1L ]] <- clustering_summary$summary
    cluster_detail_rows[[ length( cluster_detail_rows ) + 1L ]] <- clustering_summary$cluster_detail
  }
  if ( branch != "pca" || !"umap_pca_unintegrated" %in% names( obj@reductions ) ) {
    log_message( "Running UMAP for ", branch, "." )
    obj <- Seurat::RunUMAP( object = obj, reduction = reduction_name, dims = PCS_USE, reduction.name = umap_name, reduction.key = paste0( toupper( substr( cluster_prefix, 1, 4 ) ), "UMAP_" ), umap.method = "uwot", n.neighbors = UMAP_N_NEIGHBORS, n.components = 2L, metric = "cosine", min.dist = UMAP_MIN_DIST, uwot.sgd = TRUE, seed.use = SEED + branch_index * 1000L, verbose = TRUE )
  }
}
cluster_summary <- do.call( rbind, cluster_summary_rows )
cluster_detail <- do.call( rbind, cluster_detail_rows )
write_csv_safe( cluster_detail, "02_cluster_detail_all_branches_resolutions.csv" )
stability_rows <- list()
for ( branch_index in seq_len( nrow(branch_specs) ) ) {
  branch <- branch_specs$branch[ branch_index ]
  prefix <- branch_specs$cluster_prefix[ branch_index ]
  for ( i in seq_len( length(RESOLUTIONS) - 1L ) ) {
    res_a <- RESOLUTIONS[[i]]
    res_b <- RESOLUTIONS[[i + 1L]]
    col_a <- cluster_column_name( prefix, res_a )
    col_b <- cluster_column_name( prefix, res_b )
    stability_rows[[ length( stability_rows ) + 1L ]] <- data.frame( comparison_type = "adjacent_resolution_within_branch", branch_a = branch, resolution_a = res_a, branch_b = branch, resolution_b = res_b, adjusted_rand_index = adjusted_rand_index( obj@meta.data[ , col_a ], obj@meta.data[ , col_b ] ), stringsAsFactors = FALSE )
  }
}
for (resolution in RESOLUTIONS) {
  pca_col <- cluster_column_name( "pca", resolution )
  hcon_col <- cluster_column_name( "hcon", resolution )
  hdef_col <- cluster_column_name( "hdef", resolution )
  stability_rows[[ length( stability_rows ) + 1L ]] <- data.frame( comparison_type = "same_resolution_between_branches", branch_a = "pca", resolution_a = resolution, branch_b = "harmony_conservative", resolution_b = resolution, adjusted_rand_index = adjusted_rand_index( obj@meta.data[ , pca_col ], obj@meta.data[ , hcon_col ] ), stringsAsFactors = FALSE )
  stability_rows[[ length( stability_rows ) + 1L ]] <- data.frame( comparison_type = "same_resolution_between_branches", branch_a = "pca", resolution_a = resolution, branch_b = "harmony_default", resolution_b = resolution, adjusted_rand_index = adjusted_rand_index( obj@meta.data[ , pca_col ], obj@meta.data[ , hdef_col ] ), stringsAsFactors = FALSE )
  stability_rows[[ length( stability_rows ) + 1L ]] <- data.frame( comparison_type = "same_resolution_between_branches", branch_a = "harmony_conservative", resolution_a = resolution, branch_b = "harmony_default", resolution_b = resolution, adjusted_rand_index = adjusted_rand_index( obj@meta.data[ , hcon_col ], obj@meta.data[ , hdef_col ] ), stringsAsFactors = FALSE )
}
stability_table <- do.call( rbind, stability_rows )
resolution_and_stability <- merge( cluster_summary, aggregate( adjusted_rand_index ~ branch_a + resolution_a, data = stability_table[ stability_table$comparison_type == "adjacent_resolution_within_branch", , drop = FALSE ], FUN = mean ), by.x = c( "branch", "resolution" ), by.y = c( "branch_a", "resolution_a" ), all.x = TRUE, sort = FALSE )
names( resolution_and_stability )[ names( resolution_and_stability ) == "adjusted_rand_index" ] <- "mean_adjacent_resolution_ARI"
write_csv_safe( resolution_and_stability, "03_resolution_cluster_summary.csv" )
write_csv_safe( stability_table, "03B_clustering_ARI_comparisons.csv" )
log_message( "Calculating PCA within-timepoint library mixing." )
mixing_pca <- calculate_within_timepoint_mixing( embeddings = pca_embeddings, metadata = metadata, reduction_name = "pca", k = K_MIXING )
log_message( "Calculating conservative-Harmony within-timepoint library mixing." )
mixing_hcon <- calculate_within_timepoint_mixing( embeddings = harmony_conservative_embeddings, metadata = metadata, reduction_name = "harmony_conservative", k = K_MIXING )
log_message( "Calculating default-Harmony within-timepoint library mixing." )
mixing_hdef <- calculate_within_timepoint_mixing( embeddings = harmony_default_embeddings, metadata = metadata, reduction_name = "harmony_default", k = K_MIXING )
mixing_summary <- do.call( rbind, list( mixing_pca$summaries, mixing_hcon$summaries, mixing_hdef$summaries ) )
mixing_summary$timepoint <- factor( mixing_summary$timepoint, levels = expected_timepoints, ordered = TRUE )
mixing_pca_reference <- mixing_summary[ mixing_summary$reduction == "pca", c( "timepoint", "mean_other_library_fraction", "mean_normalized_library_entropy" ), drop = FALSE ]
names( mixing_pca_reference )[ names( mixing_pca_reference ) == "mean_other_library_fraction" ] <- "PCA_mean_other_library_fraction"
names( mixing_pca_reference )[ names( mixing_pca_reference ) == "mean_normalized_library_entropy" ] <- "PCA_mean_normalized_library_entropy"
mixing_summary <- merge( mixing_summary, mixing_pca_reference, by = "timepoint", all.x = TRUE, sort = FALSE )
mixing_summary$gain_vs_PCA_other_library_fraction <- mixing_summary$mean_other_library_fraction - mixing_summary$PCA_mean_other_library_fraction
mixing_summary$gain_vs_PCA_normalized_entropy <- mixing_summary$mean_normalized_library_entropy - mixing_summary$PCA_mean_normalized_library_entropy
mixing_summary <- mixing_summary[ order( mixing_summary$timepoint, match( mixing_summary$reduction, branch_specs$branch ) ), , drop = FALSE ]
neighbor_overlap_hcon <- calculate_neighbor_overlap( reference_cache = mixing_pca$neighbor_cache, comparison_cache = mixing_hcon$neighbor_cache, comparison_name = "harmony_conservative_vs_pca" )
neighbor_overlap_hdef <- calculate_neighbor_overlap( reference_cache = mixing_pca$neighbor_cache, comparison_cache = mixing_hdef$neighbor_cache, comparison_name = "harmony_default_vs_pca" )
neighbor_overlap <- rbind( neighbor_overlap_hcon, neighbor_overlap_hdef )
mixing_summary <- merge( mixing_summary, transform( neighbor_overlap, reduction = ifelse( comparison == "harmony_conservative_vs_pca", "harmony_conservative", "harmony_default" ) )[ , c( "reduction", "timepoint", "mean_KNN_overlap", "median_KNN_overlap" ), drop = FALSE ], by = c( "reduction", "timepoint" ), all.x = TRUE, sort = FALSE )
mixing_summary$mean_KNN_overlap[ mixing_summary$reduction == "pca" ] <- 1
mixing_summary$median_KNN_overlap[ mixing_summary$reduction == "pca" ] <- 1
write_csv_safe( mixing_summary, "04_within_timepoint_library_mixing.csv" )
write_csv_safe( neighbor_overlap, "04B_local_structure_preservation.csv" )
log_message( "Creating comparison figures." )
display_cluster_columns <- c( cluster_column_name( "pca", DISPLAY_RESOLUTION ), cluster_column_name( "hcon", DISPLAY_RESOLUTION ), cluster_column_name( "hdef", DISPLAY_RESOLUTION ) )
global_plot_list <- list()
for ( branch_index in seq_len( nrow(branch_specs) ) ) {
  branch <- branch_specs$branch[ branch_index ]
  umap_name <- branch_specs$umap[ branch_index ]
  cluster_column <- display_cluster_columns[ branch_index ]
  p_cluster <- Seurat::DimPlot( object = obj, reduction = umap_name, group.by = cluster_column, label = TRUE, repel = TRUE, raster = TRUE ) + ggplot2::labs( title = paste0( branch, ": clusters at resolution ", DISPLAY_RESOLUTION ) ) + Seurat::NoLegend()
  p_timepoint <- Seurat::DimPlot( object = obj, reduction = umap_name, group.by = "timepoint", raster = TRUE ) + ggplot2::labs( title = paste0( branch, ": timepoint" ) )
  p_library <- Seurat::DimPlot( object = obj, reduction = umap_name, group.by = "biological_library", raster = TRUE ) + ggplot2::labs( title = paste0( branch, ": biological library" ) )
  global_plot_list[[branch]] <- p_cluster | p_timepoint | p_library
}
global_comparison_plot <- global_plot_list[["pca"]] / global_plot_list[["harmony_conservative"]] / global_plot_list[["harmony_default"]]
ggplot2::ggsave( filename = file.path( figures_dir, "01_global_PCA_Harmony_UMAP_comparison.pdf" ), plot = global_comparison_plot, width = 22, height = 18, limitsize = FALSE )
seven_day_cells <- rownames(metadata)[ as.character( metadata$timepoint ) == "7 dpi" ]
seven_day_plot_list <- list()
for ( branch_index in seq_len( nrow(branch_specs) ) ) {
  branch <- branch_specs$branch[ branch_index ]
  umap_name <- branch_specs$umap[ branch_index ]
  coordinates <- Seurat::Embeddings( object = obj, reduction = umap_name )[ seven_day_cells, , drop = FALSE ]
  seven_day_df <- data.frame( dim1 = coordinates[ , 1 ], dim2 = coordinates[ , 2 ], biological_library = as.character( metadata[ seven_day_cells, "biological_library" ] ), stringsAsFactors = FALSE )
  seven_day_plot_list[[branch]] <- ggplot2::ggplot( seven_day_df, ggplot2::aes( x = dim1, y = dim2, colour = biological_library ) ) + ggplot2::geom_point( size = 0.25, alpha = 0.60 ) + ggplot2::coord_equal() + ggplot2::labs( x = "UMAP 1", y = "UMAP 2", colour = NULL, title = paste0( branch, ": 7 dpi libraries" ) ) + ggplot2::theme_bw()
}
seven_day_comparison_plot <- seven_day_plot_list[["pca"]] | seven_day_plot_list[["harmony_conservative"]] | seven_day_plot_list[["harmony_default"]]
ggplot2::ggsave( filename = file.path( figures_dir, "02_7dpi_library_mixing_comparison.pdf" ), plot = seven_day_comparison_plot, width = 18, height = 6, limitsize = FALSE )
resolution_pdf <- file.path( figures_dir, "03_resolution_0.1_to_0.4_UMAP_comparison.pdf" )
grDevices::pdf( resolution_pdf, width = 16, height = 12, onefile = TRUE )
for ( branch_index in seq_len( nrow(branch_specs) ) ) {
  branch <- branch_specs$branch[ branch_index ]
  prefix <- branch_specs$cluster_prefix[ branch_index ]
  umap_name <- branch_specs$umap[ branch_index ]
  branch_plots <- lapply( RESOLUTIONS, function(resolution) {
      cluster_column <- cluster_column_name( prefix, resolution )
      Seurat::DimPlot( object = obj, reduction = umap_name, group.by = cluster_column, label = TRUE, repel = TRUE, raster = TRUE ) + ggplot2::labs( title = paste0( branch, " | resolution ", resolution ) ) + Seurat::NoLegend()
    }
  )
  print( patchwork::wrap_plots( branch_plots, ncol = 2 ) )
}
grDevices::dev.off()
log_message( "Wrote: ", resolution_pdf )
obj@misc$GSE234774_step10C <- list( input_object = normalizePath( input_object_path, winslash = "/", mustWork = TRUE ), reductions_compared = branch_specs$reduction, resolutions_compared = RESOLUTIONS, display_resolution = DISPLAY_RESOLUTION, harmony_conservative = list( theta = HCON_THETA, lambda = HCON_LAMBDA ), harmony_default = list( theta = HDEF_THETA, lambda = HDEF_LAMBDA ), harmony_grouping_variable = "biological_library", mixing_evaluation = "within each timepoint only", k_graph = K_GRAPH, k_mixing = K_MIXING, pca_dimensions = PCS_USE, RNA_expression_changed = FALSE, final_reduction_frozen = FALSE, final_resolution_frozen = FALSE, biological_unit = "biological_library", completed_at = as.character( Sys.time() ) )
saveRDS( obj, output_object_path, compress = FALSE )
log_message( "Saved comparison object: ", output_object_path )
reload_obj <- tryCatch( readRDS( output_object_path ), error = function(e) e )
required_reductions <- c( "pca", "harmony_conservative", "harmony_default", "umap_pca_unintegrated", "umap_harmony_conservative", "umap_harmony_default" )
required_cluster_columns <- unlist( lapply( seq_len( nrow(branch_specs) ), function(branch_index) {
      vapply( RESOLUTIONS, function(resolution) {
          cluster_column_name( branch_specs$cluster_prefix[ branch_index ], resolution )
        },
        character(1)
      )
    }
  ), use.names = FALSE
)
reload_pass <- !inherits( reload_obj, "error" ) && inherits( reload_obj, "Seurat" ) && nrow(reload_obj) == nrow(obj) && ncol(reload_obj) == ncol(obj) && all( required_reductions %in% names( reload_obj@reductions ) ) && all( required_cluster_columns %in% names( reload_obj@meta.data ) )
rm( reload_obj )
pca_integrity_pass <- identical( pca_copy_for_integrity, Seurat::Embeddings( object = obj, reduction = "pca" )[ , PCS_USE, drop = FALSE ] )
cell_order_pass <- identical( cell_order_before, colnames(obj) )
seven_day_summary <- mixing_summary[ as.character( mixing_summary$timepoint ) == "7 dpi", , drop = FALSE ]
seven_day_pca_mixing <- seven_day_summary$mean_other_library_fraction[ seven_day_summary$reduction == "pca" ]
seven_day_hcon_mixing <- seven_day_summary$mean_other_library_fraction[ seven_day_summary$reduction == "harmony_conservative" ]
seven_day_hdef_mixing <- seven_day_summary$mean_other_library_fraction[ seven_day_summary$reduction == "harmony_default" ]
mean_hcon_overlap <- mean( neighbor_overlap$mean_KNN_overlap[ neighbor_overlap$comparison == "harmony_conservative_vs_pca" ], na.rm = TRUE )
mean_hdef_overlap <- mean( neighbor_overlap$mean_KNN_overlap[ neighbor_overlap$comparison == "harmony_default_vs_pca" ], na.rm = TRUE )
audit_checks <- data.frame( check = c( "Input contains 26,918 features and 144,359 singlets", "PCA cell order and coordinates were preserved", "Conservative Harmony contains all nuclei and 40 dimensions", "Default-strength Harmony contains all nuclei and 40 dimensions", "Three UMAPs contain all nuclei", "All 12 branch-resolution cluster columns were created", "Within-timepoint library-mixing metrics were calculated for all branches and timepoints", "Local-structure preservation metrics were calculated for both Harmony branches", "Saved object reloads with all reductions and cluster columns", "RNA counts and normalized expression were not replaced", "Final reduction has been frozen", "Final clustering resolution has been frozen", "Harmony and resolution figures have been manually reviewed" ), required_for_pass = c( TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE ), status = c( ifelse( nrow(obj) == 26918L && ncol(obj) == 144359L, "PASS", "STOP" ), ifelse( pca_integrity_pass && cell_order_pass, "PASS", "STOP" ), ifelse( nrow( harmony_conservative_embeddings ) == ncol(obj) && ncol( harmony_conservative_embeddings ) == length(PCS_USE), "PASS", "STOP" ), ifelse( nrow( harmony_default_embeddings ) == ncol(obj) && ncol( harmony_default_embeddings ) == length(PCS_USE), "PASS", "STOP" ), ifelse( all( vapply( branch_specs$umap, function(umap_name) {
            nrow( Seurat::Embeddings( obj, umap_name ) ) == ncol(obj)
          },
          logical(1)
        )
      ), "PASS", "STOP"
    ), ifelse(
      all( required_cluster_columns %in% names( obj@meta.data ) ), "PASS", "STOP"
    ), ifelse(
      nrow(mixing_summary) == nrow(branch_specs) * length(expected_timepoints) && !anyNA( mixing_summary$mean_other_library_fraction ), "PASS", "STOP"
    ), ifelse(
      nrow(neighbor_overlap) == 2L * length(expected_timepoints) && !anyNA( neighbor_overlap$mean_KNN_overlap ), "PASS", "STOP"
    ), ifelse(
      reload_pass, "PASS", "STOP"
    ), "CONFIRMED_NOT_CHANGED", "NOT_FROZEN", "NOT_FROZEN", "REVIEW"
  ), detail = c(
    paste0( nrow(obj), " features x ", ncol(obj), " nuclei" ), paste0( "PCA integrity=", pca_integrity_pass, "; cell order=", cell_order_pass ), paste0( dim( harmony_conservative_embeddings ), collapse = " x " ), paste0( dim( harmony_default_embeddings ), collapse = " x " ), paste( branch_specs$umap, collapse = ", " ), paste0( length( required_cluster_columns ), " columns" ), paste0( nrow(mixing_summary), " branch-timepoint rows" ), paste0( nrow(neighbor_overlap), " comparison-timepoint rows" ), ifelse( reload_pass, "RDS reload PASS", "RDS reload failed" ), "Harmony was applied to PCA embeddings only.", "Pending review.", "Pending review.", paste0( "7 dpi other-library fraction: PCA=", round( seven_day_pca_mixing, 3 ), "; conservative=", round( seven_day_hcon_mixing, 3 ), "; default=", round( seven_day_hdef_mixing, 3 ), ". Mean KNN overlap: conservative=", round( mean_hcon_overlap, 3 ), "; default=", round( mean_hdef_overlap, 3 ), "." )
  ), stringsAsFactors = FALSE
)
write_csv_safe( audit_checks, "09_STEP10C_audit_checks.csv" )
required_fail <- audit_checks$required_for_pass & audit_checks$status != "PASS"
final_decision <- if ( any( required_fail ) ) {
  "STOP"
} else {
  "PASS_WITH_MANUAL_METHOD_SELECTION"
}
cluster_lines <- apply( resolution_and_stability[ order( match( resolution_and_stability$branch, branch_specs$branch ), resolution_and_stability$resolution ), , drop = FALSE ], 1L, function(row) {
    paste0( row[["branch"]], " resolution ", row[["resolution"]], ": ", row[["n_clusters"]], " clusters; min/median/max=", row[["minimum_cluster_size"]], "/", round( as.numeric( row[["median_cluster_size"]] ), 1 ), "/", row[["maximum_cluster_size"]], "; <500=", row[["n_clusters_lt_500"]], "; library-dominant=", row[["n_library_dominant_clusters"]], "; timepoint-dominant=", row[["n_timepoint_dominant_clusters"]] )
  }
)
mixing_lines <- apply( mixing_summary[ order( mixing_summary$timepoint, match( mixing_summary$reduction, branch_specs$branch ) ), , drop = FALSE ], 1L, function(row) {
    paste0( row[["timepoint"]], " | ", row[["reduction"]], ": other-library fraction=", round( as.numeric( row[["mean_other_library_fraction"]] ), 3 ), "; normalized entropy=", round( as.numeric( row[["mean_normalized_library_entropy"]] ), 3 ), "; KNN overlap with PCA=", round( as.numeric( row[["mean_KNN_overlap"]] ), 3 ) )
  }
)
report_lines <- c( "GSE234774 STEP 10C HARMONY / RESOLUTION SENSITIVITY", "============================================================", paste0( "Run time: ", format( Sys.time(), "%Y-%m-%d %H:%M:%S" ) ), paste0( "Input object: ", input_object_path ), paste0( "Output directory: ", output_dir ), "", paste0( "FINAL DECISION: ", final_decision ), "", "DESIGN", "Unintegrated PCA was preserved.", paste0( "Conservative Harmony: theta=", HCON_THETA, ", lambda=", HCON_LAMBDA, "." ), paste0( "Default-strength Harmony: theta=", HDEF_THETA, ", lambda=", HDEF_LAMBDA, "." ), "Harmony grouping variable: biological_library.", "Library mixing was evaluated separately within each timepoint.", paste0( "Resolutions compared: ", paste( RESOLUTIONS, collapse = ", " ), "." ), "", "7 DPI MIXING", paste0( "PCA mean other-library neighbor fraction: ", round( seven_day_pca_mixing, 3 ) ), paste0( "Conservative Harmony: ", round( seven_day_hcon_mixing, 3 ) ), paste0( "Default-strength Harmony: ", round( seven_day_hdef_mixing, 3 ) ), "", "LOCAL STRUCTURE PRESERVATION", paste0( "Mean conservative-Harmony versus PCA KNN overlap: ", round( mean_hcon_overlap, 3 ) ), paste0( "Mean default-Harmony versus PCA KNN overlap: ", round( mean_hdef_overlap, 3 ) ), paste0( "Review threshold: mean overlap < ", MIN_ACCEPTABLE_MEAN_KNN_OVERLAP, " suggests strong structural alteration." ), "", "CLUSTERING SUMMARY", cluster_lines, "", "WITHIN-TIMEPOINT MIXING SUMMARY", mixing_lines, "", "INTERPRETATION BOUNDARIES", "1. Better library mixing alone is not sufficient to choose Harmony.", "2. The selected branch must retain recognizable biological structure and avoid collapsing timepoint-associated states.", "3. Resolution is not selected solely by the number of clusters.", "4. Whole-atlas annotation should favor a stable coarse solution; finer myeloid states will be resolved after myeloid extraction.", "5. RNA counts and LogNormalize expression remain unchanged.", "6. No reduction or resolution is frozen by this script.", "", "PRIMARY COMPARISON OBJECT", output_object_path, "", "REQUIRED REVIEW", "Review all three exported PDFs together with the mixing, resolution, audit, and final-report tables.", "Then select one reduction branch and one coarse whole-atlas resolution before marker-first annotation." )
if ( final_decision == "STOP" ) {
  failed <- audit_checks[ audit_checks$required_for_pass & audit_checks$status != "PASS", , drop = FALSE ]
  report_lines <- c( report_lines, "", "STOP REASONS" )
  for ( i in seq_len( nrow(failed) ) ) {
    report_lines <- c( report_lines, paste0( i, ". ", failed$check[[i]], " -- ", failed$detail[[i]] ) )
  }
}
report_path <- file.path( output_dir, "10_STEP10C_FINAL_REPORT.txt" )
writeLines( report_lines, report_path, useBytes = TRUE )
log_message( "Final decision: ", final_decision )
log_message( "Final report: ", report_path )
cat("\n")
cat("============================================================\n")
cat("GSE234774 STEP 10C COMPLETED\n")
cat( "FINAL DECISION:", final_decision, "\n" )
cat( "Output:", output_dir, "\n" )
cat("Comparison object:\n")
cat( " ", output_object_path, "\n", sep = "" )
cat("Key outputs:\n")
cat(" 1) figures/01_global_PCA_Harmony_UMAP_comparison.pdf\n")
cat(" 2) figures/02_7dpi_library_mixing_comparison.pdf\n")
cat(" 3) figures/03_resolution_0.1_to_0.4_UMAP_comparison.pdf\n")
cat(" 4) tables/03_resolution_cluster_summary.csv\n")
cat(" 5) tables/04_within_timepoint_library_mixing.csv\n")
cat(" 6) tables/09_STEP10C_audit_checks.csv\n")
cat(" 7) 10_STEP10C_FINAL_REPORT.txt\n")
cat("============================================================\n")
if ( final_decision == "STOP" ) {
  stop( "GSE234774 Step 10C decision is STOP. Review the final report.", call. = FALSE )
}
