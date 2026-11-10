options( stringsAsFactors = FALSE, warn = 1, Seurat.object.assay.version = "v5", future.globals.maxSize = 16 * 1024^3 )
SEED <- 20260803L
FROZEN_REDUCTION <- "harmony_myeloid_conservative"
FROZEN_UMAP <- "umap_myeloid_harmony"
FROZEN_CLUSTER_COLUMN <- "myhcon_res_0_6"
EXPECTED_FEATURES <- 26918L
EXPECTED_MYLOID_NUCLEI <- 25205L
EXPECTED_CLUSTERS <- as.character(0:8)
MAX_CELLS_PER_CLUSTER_FOR_MARKERS <- 3000L
MARKER_LOGFC_THRESHOLD <- 0.25
MARKER_MIN_PCT <- 0.10
MARKER_MIN_DIFF_PCT <- 0.05
MARKER_RETURN_P_ADJ <- 0.05
TOP_MARKERS_PER_CLUSTER <- 30L
REUSE_COMPLETED_MARKERS <- TRUE
LIBRARY_DOMINANCE_REVIEW_THRESHOLD <- 0.80
TIMEPOINT_DOMINANCE_REVIEW_THRESHOLD <- 0.95
MIN_CLUSTER_SIZE_FOR_DOMINANCE_REVIEW <- 200L
expected_timepoints <- c( "Uninjured", "1 dpi", "4 dpi", "7 dpi", "14 dpi", "1 month", "2 months" )
working_label_map <- c( "0" = "Homeostatic_resident_biased_microglia", "1" = "Activated_transitional_microglia", "2" = "Phagocytic_DAM_like_microglia", "3" = "Intermediate_resolving_like_microglia", "4" = "MHCII_BAM_like_ambiguous_myeloid", "5" = "Macrophage_like_activated_ambiguous_myeloid", "6" = "Cycling_myeloid", "7" = "IFN_responsive_microglia", "8" = "Inflammatory_proliferative_ambiguous_myeloid" )
resident_eligibility_map <- c( "0" = "LIKELY_RESIDENT_LINEAGE", "1" = "LIKELY_RESIDENT_LINEAGE", "2" = "LIKELY_RESIDENT_LINEAGE", "3" = "LIKELY_RESIDENT_LINEAGE", "4" = "AMBIGUOUS_RESIDENT_VS_INFILTRATING", "5" = "AMBIGUOUS_RESIDENT_VS_INFILTRATING", "6" = "LIKELY_RESIDENT_LINEAGE", "7" = "LIKELY_RESIDENT_LINEAGE", "8" = "AMBIGUOUS_RESIDENT_VS_INFILTRATING" )
expected_primary_panel_map <- c( "0" = "Resident_homeostatic", "1" = "Activated_DAM_phagolysosomal", "2" = "Activated_DAM_phagolysosomal", "3" = "Resident_homeostatic", "4" = "MHCII_antigen_presentation", "5" = "Macrophage_monocyte_like", "6" = "Cell_cycle", "7" = "Interferon_response", "8" = "Inflammatory_chemokine" )
if ( !identical( names(working_label_map), EXPECTED_CLUSTERS ) || !identical( names(resident_eligibility_map), EXPECTED_CLUSTERS ) || !identical( names(expected_primary_panel_map), EXPECTED_CLUSTERS ) ) {
  stop( "The frozen cluster maps match clusters 0-8.", call. = FALSE )
}
marker_panels_requested <- list( Myeloid_general = c( "Ptprc", "Aif1", "Tyrobp", "Lyz2", "Csf1r" ), Resident_homeostatic = c( "P2ry12", "Tmem119", "Sall1", "Fcrls", "Hexb", "Cx3cr1", "Gpr34", "Olfml3" ), Activated_DAM_phagolysosomal = c( "Apoe", "Cst7", "Trem2", "Itgax", "Lgals3", "Cd68", "Ctsb", "Ctsd", "Lpl", "Spp1", "Gpnmb" ), Macrophage_monocyte_like = c( "Ccr2", "Ly6c2", "Plac8", "Ms4a7", "Fcgr3", "S100a8", "S100a9", "Retnlg" ), BAM = c( "Mrc1", "Lyve1", "Pf4", "Folr2", "Cd163", "Cbr2" ), MHCII_antigen_presentation = c( "Cd74", "H2-Ab1", "H2-Aa", "Ciita", "H2-Eb1", "Fcgr2b" ), Interferon_response = c( "Isg15", "Ifit1", "Ifit3", "Stat1", "Rsad2", "Irf7" ), Inflammatory_chemokine = c( "Il1b", "Tnf", "Nfkbia", "Ccl2", "Ccl3", "Ccl4" ), Cell_cycle = c( "Mki67", "Top2a", "Stmn1", "Tubb5", "Cenpf", "Ube2c", "Bub1", "Kif11" ) )
project_root <- Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
input_object_path <- file.path( project_root, "04_results", "step10E_GSE234774_myeloid_reclustering_Harmony_resolution", "RDS", "GSE234774_step10E_myeloid_PCA_Harmony_resolution_comparison_BPCells.rds" )
output_dir <- file.path( project_root, "04_results", "step10F_GSE234774_myeloid_marker_first_annotation_evidence" )
tables_dir <- file.path( output_dir, "tables" )
figures_dir <- file.path( output_dir, "figures" )
rds_dir <- file.path( output_dir, "RDS" )
for (directory in c( output_dir, tables_dir, figures_dir, rds_dir )) {
  dir.create( directory, recursive = TRUE, showWarnings = FALSE )
}
run_log <- file.path( output_dir, "00_run_log.txt" )
if (file.exists( run_log )) {
  file.remove( run_log )
}
output_object_path <- file.path( rds_dir, "GSE234774_step10F_myeloid_annotation_evidence_BPCells.rds" )
write_csv_safe <- function( x, filename ) {
  path <- file.path( tables_dir, filename )
  utils::write.csv( x, path, row.names = FALSE, na = "" )
  invisible( path )
}
stop_with_report <- function(reason) {
  report_path <- file.path( output_dir, "10_STEP10F_FINAL_REPORT.txt" )
  stop( reason, call. = FALSE )
}
safe_package_version <- function( package_name ) {
  if (!requireNamespace( package_name, quietly = TRUE )) {
    return( NA_character_ )
  }
  as.character( utils::packageVersion( package_name ) )
}
resolve_features <- function( requested_features, available_features ) {
  available_upper <- toupper( available_features )
  first_index <- !duplicated( available_upper )
  lookup <- setNames( available_features[ first_index ], available_upper[ first_index ] )
  resolved <- unname( lookup[ toupper( requested_features ) ] )
  names( resolved ) <- requested_features
  resolved
}
deduplicate_plot_panels <- function( panel_list ) {
  used <- character(0)
  output <- list()
  for (panel_name in names( panel_list )) {
    genes <- panel_list[[panel_name]]
    genes <- genes[ !genes %in% used ]
    if (length( genes ) > 0L) {
      output[[panel_name]] <- genes
      used <- c( used, genes )
    }
  }
  output
}
safe_zscore <- function(x) {
  value_sd <- stats::sd( x, na.rm = TRUE )
  if ( !is.finite( value_sd ) || value_sd == 0 ) {
    return( rep( 0, length( x ) ) )
  }
  output <- ( x - mean( x, na.rm = TRUE ) ) / value_sd
  output[ !is.finite( output ) ] <- 0
  output
}
dominant_composition <- function( metadata, cluster_column, value_column, value_name, fraction_name ) {
  composition <- data.table::data.table( cluster = as.character( metadata[ , cluster_column ] ), value = as.character( metadata[ , value_column ] ) )[ , .( n_cells = .N ), by = .( cluster, value ) ]
  totals <- composition[ , .( cluster_total = sum( n_cells ) ), by = cluster ]
  composition <- merge( composition, totals, by = "cluster", all.x = TRUE, sort = FALSE )
  composition$fraction <- composition$n_cells / composition$cluster_total
  top <- composition[ order( cluster, -fraction, value ), .SD[ 1L ], by = cluster ]
  output <- data.frame( cluster = top$cluster, dominant_value = top$value, dominant_fraction = top$fraction, stringsAsFactors = FALSE )
  names( output )[ names( output ) == "dominant_value" ] <- value_name
  names( output )[ names( output ) == "dominant_fraction" ] <- fraction_name
  output
}
plot_score_umap <- function( plot_data, score_column, title_text ) {
  current_data <- plot_data
  current_data$score <- current_data[[score_column]]
  ggplot2::ggplot( current_data, ggplot2::aes( x = UMAP_1, y = UMAP_2 ) ) + ggplot2::geom_point( ggplot2::aes( colour = score ), size = 0.25, alpha = 0.70 ) + ggplot2::scale_colour_viridis_c( option = "magma" ) + ggplot2::coord_equal() + ggplot2::labs( title = title_text, x = "Myeloid Harmony UMAP 1", y = "Myeloid Harmony UMAP 2", colour = "Z score" ) + ggplot2::theme_bw()
}
required_packages <- c( "Matrix", "data.table", "ggplot2", "patchwork", "scales", "Seurat", "SeuratObject" )
missing_packages <- required_packages[ !vapply( required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1) ) ]
if ( length( missing_packages ) > 0L ) {
  stop_with_report( paste0( "Missing required package(s): ", paste( missing_packages, collapse = ", " ), "." ) )
}
suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(Seurat)
  library(SeuratObject)
})
package_versions <- data.frame( package = c( "R", "Matrix", "data.table", "ggplot2", "patchwork", "scales", "Seurat", "SeuratObject" ), version = c( paste0( R.version$major, ".", R.version$minor ), safe_package_version( "Matrix" ), safe_package_version( "data.table" ), safe_package_version( "ggplot2" ), safe_package_version( "patchwork" ), safe_package_version( "scales" ), safe_package_version( "Seurat" ), safe_package_version( "SeuratObject" ) ), stringsAsFactors = FALSE )
write_csv_safe( package_versions, "00_package_versions.csv" )
if (!file.exists( input_object_path )) {
  stop_with_report( paste0( "Input Step 10E object not found: ", input_object_path ) )
}
obj <- readRDS( input_object_path )
if (!inherits( obj, "Seurat" )) {
  stop_with_report( "Input RDS is not a Seurat object." )
}
if ( nrow( obj ) != EXPECTED_FEATURES || ncol( obj ) != EXPECTED_MYLOID_NUCLEI ) {
  stop_with_report( paste0( "Expected ", EXPECTED_FEATURES, " features x ", EXPECTED_MYLOID_NUCLEI, " nuclei; detected ", nrow( obj ), " x ", ncol( obj ), "." ) )
}
required_metadata <- c( "biological_library", "timepoint", "replicate_id", "parent_atlas_cluster", "Phase", "S.Score", "G2M.Score", "nCount_RNA", "nFeature_RNA", "percent.mt", FROZEN_CLUSTER_COLUMN )
missing_metadata <- setdiff( required_metadata, names( obj@meta.data ) )
if ( length( missing_metadata ) > 0L ) {
  stop_with_report( paste0( "Missing required metadata column(s): ", paste( missing_metadata, collapse = ", " ), "." ) )
}
required_reductions <- c( FROZEN_REDUCTION, FROZEN_UMAP )
missing_reductions <- setdiff( required_reductions, names( obj@reductions ) )
if ( length( missing_reductions ) > 0L ) {
  stop_with_report( paste0( "Missing frozen reduction(s): ", paste( missing_reductions, collapse = ", " ), "." ) )
}
cluster_values <- as.character( obj@meta.data[ , FROZEN_CLUSTER_COLUMN ] )
observed_clusters <- sort( unique( cluster_values ), method = "radix" )
if ( !setequal( observed_clusters, EXPECTED_CLUSTERS ) ) {
  stop_with_report( paste0( "Expected frozen clusters 0-8; detected ", paste( observed_clusters, collapse = ", " ), "." ) )
}
obj@meta.data[ , FROZEN_CLUSTER_COLUMN ] <- factor( cluster_values, levels = EXPECTED_CLUSTERS )
Seurat::DefaultAssay( obj ) <- "RNA"
Seurat::Idents( obj ) <- FROZEN_CLUSTER_COLUMN
rna_layers <- SeuratObject::Layers( object = obj, assay = "RNA", search = NA )
counts_layers <- rna_layers[ grepl( "^counts($|\\.)", rna_layers ) ]
data_layers <- rna_layers[ grepl( "^data($|\\.)", rna_layers ) ]
if ( length( counts_layers ) != 1L || length( data_layers ) != 1L ) {
  stop_with_report( paste0( "Expected one RNA counts layer and one data layer; detected ", length( counts_layers ), " and ", length( data_layers ), "." ) )
}
counts_layer_name <- counts_layers[[1L]]
data_layer_name <- data_layers[[1L]]
input_audit <- data.frame( item = c( "Input object", "Features", "Myeloid nuclei", "Biological libraries", "Timepoints", "Frozen reduction", "Frozen UMAP", "Frozen cluster column", "Frozen cluster count", "RNA counts layer", "RNA data layer" ), value = c( normalizePath( input_object_path, winslash = "/", mustWork = TRUE ), nrow( obj ), ncol( obj ), length( unique( as.character( obj$biological_library ) ) ), length( unique( as.character( obj$timepoint ) ) ), FROZEN_REDUCTION, FROZEN_UMAP, FROZEN_CLUSTER_COLUMN, length( observed_clusters ), counts_layer_name, data_layer_name ), stringsAsFactors = FALSE )
write_csv_safe( input_audit, "01_input_object_audit.csv" )
all_requested_features <- unique( unlist( marker_panels_requested, use.names = FALSE ) )
feature_resolution <- resolve_features( requested_features = all_requested_features, available_features = rownames( obj ) )
feature_resolution_table <- data.frame( requested_gene = names( feature_resolution ), resolved_feature = unname( feature_resolution ), present = !is.na( feature_resolution ), stringsAsFactors = FALSE )
write_csv_safe( feature_resolution_table, "02_marker_feature_resolution.csv" )
marker_panels_present <- lapply( marker_panels_requested, function(panel_genes) {
    resolved <- feature_resolution[ panel_genes ]
    unique( unname( resolved[ !is.na( resolved ) ] ) )
  }
)
panel_availability <- data.frame( panel = names( marker_panels_requested ), requested_n = vapply( marker_panels_requested, length, integer(1) ), present_n = vapply( marker_panels_present, length, integer(1) ), present_features = vapply( marker_panels_present, function(features) {
      paste( features, collapse = ";" )
    },
    character(1)
  ), stringsAsFactors = FALSE
)
panel_availability$present_fraction <- panel_availability$present_n / panel_availability$requested_n
write_csv_safe( panel_availability, "02B_marker_panel_availability.csv" )
if ( any( panel_availability$present_n < 2L ) ) {
  failed_panels <- panel_availability$panel[ panel_availability$present_n < 2L ]
  stop_with_report( paste0( "Fewer than two genes were resolved for panel(s): ", paste( failed_panels, collapse = ", " ), "." ) )
}
all_present_panel_features <- unique( unlist( marker_panels_present, use.names = FALSE ) )
panel_data <- SeuratObject::LayerData( object = obj, assay = "RNA", layer = data_layer_name, features = all_present_panel_features, fast = FALSE )
panel_counts <- SeuratObject::LayerData( object = obj, assay = "RNA", layer = counts_layer_name, features = all_present_panel_features, fast = FALSE )
panel_data <- methods::as( panel_data, "dgCMatrix" )
panel_counts <- methods::as( panel_counts, "dgCMatrix" )
if ( !identical( colnames( panel_data ), colnames( obj ) ) || !identical( colnames( panel_counts ), colnames( obj ) ) ) {
  stop_with_report( "Panel expression matrices match the Seurat cell order." )
}
score_columns_raw <- character(0)
score_columns_z <- character(0)
for (panel_name in names( marker_panels_present )) {
  panel_features <- marker_panels_present[[panel_name]]
  raw_score <- Matrix::colMeans( panel_data[ panel_features, , drop = FALSE ] )
  z_score <- safe_zscore( raw_score )
  raw_column <- paste0( "panel_mean_", panel_name )
  z_column <- paste0( "panel_z_", panel_name )
  obj@meta.data[ , raw_column ] <- as.numeric( raw_score )
  obj@meta.data[ , z_column ] <- as.numeric( z_score )
  score_columns_raw <- c( score_columns_raw, raw_column )
  score_columns_z <- c( score_columns_z, z_column )
}
cell_clusters <- as.character( obj@meta.data[ , FROZEN_CLUSTER_COLUMN ] )
obj$myeloid_working_label <- unname( working_label_map[ cell_clusters ] )
obj$resident_lineage_eligibility <- unname( resident_eligibility_map[ cell_clusters ] )
obj$myeloid_annotation_status <- "WORKING_LABELS_NOT_FROZEN"
if ( anyNA( obj$myeloid_working_label ) || anyNA( obj$resident_lineage_eligibility ) ) {
  stop_with_report( "Working labels were not mapped to every myeloid nucleus." )
}
all_markers_path <- file.path( tables_dir, "03_all_positive_myeloid_cluster_markers.csv.gz" )
top_markers_path <- file.path( tables_dir, "04_top30_positive_markers_by_myeloid_cluster.csv" )
reuse_marker_table <- isTRUE( REUSE_COMPLETED_MARKERS ) && file.exists( top_markers_path )
if (reuse_marker_table) {
  top_markers <- data.table::fread( top_markers_path, data.table = FALSE, check.names = FALSE )
  if ( !all( c( "cluster", "gene" ) %in% names( top_markers ) ) ) {
    stop_with_report( "The completed top-marker table lacks cluster or gene columns." )
  }
  fold_change_column <- intersect( c( "avg_log2FC", "avg_logFC", "avg_diff" ), names( top_markers ) )
  if ( length( fold_change_column ) == 0L ) {
    stop_with_report( "Could not identify the fold-change column in the completed top-marker table." )
  }
  fold_change_column <- fold_change_column[[1L]]
  top_markers$cluster <- as.character( top_markers$cluster )
} else {
  set.seed( SEED )
  all_markers <- Seurat::FindAllMarkers( object = obj, assay = "RNA", slot = "data", only.pos = TRUE, test.use = "wilcox", logfc.threshold = MARKER_LOGFC_THRESHOLD, min.pct = MARKER_MIN_PCT, min.diff.pct = MARKER_MIN_DIFF_PCT, return.thresh = MARKER_RETURN_P_ADJ, max.cells.per.ident = MAX_CELLS_PER_CLUSTER_FOR_MARKERS, random.seed = SEED, densify = FALSE, verbose = TRUE )
  if ( is.null( all_markers ) || nrow( all_markers ) == 0L ) {
    stop_with_report( "FindAllMarkers returned no positive markers." )
  }
  if ( !all( c( "cluster", "gene" ) %in% names( all_markers ) ) ) {
    stop_with_report( "FindAllMarkers output lacks cluster or gene columns." )
  }
  fold_change_column <- intersect( c( "avg_log2FC", "avg_logFC", "avg_diff" ), names( all_markers ) )
  if ( length( fold_change_column ) == 0L ) {
    stop_with_report( "Could not identify the marker fold-change column." )
  }
  fold_change_column <- fold_change_column[[1L]]
  all_markers$cluster <- as.character( all_markers$cluster )
  data.table::fwrite( all_markers, file = all_markers_path, sep = ",", quote = TRUE, na = "NA", compress = "gzip" )
  marker_table <- data.table::as.data.table( all_markers )
  marker_table[ , marker_is_technical := grepl( "^(mt-|Rpl|Rps|Hb[ab])", gene, ignore.case = TRUE ) ]
  rankable_markers <- marker_table[ marker_is_technical == FALSE ]
  data.table::setorderv( rankable_markers, cols = c( "cluster", fold_change_column ), order = c( 1L, -1L ), na.last = TRUE )
  top_markers <- rankable_markers[ , head( .SD, TOP_MARKERS_PER_CLUSTER ), by = cluster ]
  top_markers <- as.data.frame( top_markers )
  write_csv_safe( top_markers, "04_top30_positive_markers_by_myeloid_cluster.csv" )
}
marker_cluster_coverage <- sort( unique( as.character( top_markers$cluster ) ), method = "radix" )
if ( !setequal( marker_cluster_coverage, EXPECTED_CLUSTERS ) ) {
  stop_with_report( paste0( "Top-marker table does not cover all clusters 0-8. Detected: ", paste( marker_cluster_coverage, collapse = ", " ), "." ) )
}
top_marker_strings <- data.table::as.data.table( top_markers )[ order( cluster, -get( fold_change_column ) ), .( top10_markers = paste( head( gene, 10L ), collapse = ";" ) ), by = cluster ]
cluster_metadata <- obj@meta.data
cluster_metadata$cluster <- as.character( cluster_metadata[ , FROZEN_CLUSTER_COLUMN ] )
cluster_summary <- data.table::as.data.table( cluster_metadata )[ , .( n_cells = .N, n_libraries = uniqueN( biological_library ), n_timepoints = uniqueN( timepoint ), median_nFeature_RNA = median( nFeature_RNA ), median_nCount_RNA = median( nCount_RNA ), median_percent_mt = median( percent.mt ), fraction_G1 = mean( Phase == "G1" ), fraction_S = mean( Phase == "S" ), fraction_G2M = mean( Phase == "G2M" ), fraction_parent_cluster_21 = mean( parent_atlas_cluster == "21" ) ), by = cluster ]
dominant_library <- dominant_composition( metadata = cluster_metadata, cluster_column = "cluster", value_column = "biological_library", value_name = "dominant_library", fraction_name = "maximum_library_fraction" )
dominant_timepoint <- dominant_composition( metadata = cluster_metadata, cluster_column = "cluster", value_column = "timepoint", value_name = "dominant_timepoint", fraction_name = "maximum_timepoint_fraction" )
cluster_summary <- merge( cluster_summary, dominant_library, by = "cluster", all.x = TRUE, sort = FALSE )
cluster_summary <- merge( cluster_summary, dominant_timepoint, by = "cluster", all.x = TRUE, sort = FALSE )
for (panel_name in names( marker_panels_present )) {
  raw_column <- paste0( "panel_mean_", panel_name )
  z_column <- paste0( "panel_z_", panel_name )
  score_summary <- data.table::as.data.table( cluster_metadata )[ , .( panel_raw_mean = mean( get( raw_column ) ), panel_raw_median = median( get( raw_column ) ), panel_z_mean = mean( get( z_column ) ), panel_z_median = median( get( z_column ) ) ), by = cluster ]
  names( score_summary )[ names( score_summary ) == "panel_raw_mean" ] <- paste0( "mean_", raw_column )
  names( score_summary )[ names( score_summary ) == "panel_raw_median" ] <- paste0( "median_", raw_column )
  names( score_summary )[ names( score_summary ) == "panel_z_mean" ] <- paste0( "mean_", z_column )
  names( score_summary )[ names( score_summary ) == "panel_z_median" ] <- paste0( "median_", z_column )
  cluster_summary <- merge( cluster_summary, score_summary, by = "cluster", all.x = TRUE, sort = FALSE )
}
cluster_factor <- factor( cluster_metadata$cluster, levels = EXPECTED_CLUSTERS )
cluster_index <- as.integer( cluster_factor )
cluster_sizes <- tabulate( cluster_index, nbins = length( EXPECTED_CLUSTERS ) )
names( cluster_sizes ) <- EXPECTED_CLUSTERS
cluster_indicator <- Matrix::sparseMatrix( i = seq_len( ncol( obj ) ), j = cluster_index, x = 1, dims = c( ncol( obj ), length( EXPECTED_CLUSTERS ) ), dimnames = list( colnames( obj ), EXPECTED_CLUSTERS ) )
detection_counts <- as.matrix( ( panel_counts > 0 ) %*% cluster_indicator )
detection_fraction <- sweep( detection_counts, MARGIN = 2L, STATS = cluster_sizes, FUN = "/" )
for (panel_name in names( marker_panels_present )) {
  panel_features <- marker_panels_present[[panel_name]]
  panel_detection_mean <- colMeans( detection_fraction[ panel_features, , drop = FALSE ], na.rm = TRUE )
  cluster_summary[ , paste0( "mean_detection_", panel_name ) ] <- panel_detection_mean[ cluster_summary$cluster ]
}
cluster_summary$working_label <- unname( working_label_map[ cluster_summary$cluster ] )
cluster_summary$resident_lineage_eligibility <- unname( resident_eligibility_map[ cluster_summary$cluster ] )
cluster_summary$expected_primary_panel <- unname( expected_primary_panel_map[ cluster_summary$cluster ] )
cluster_summary <- as.data.frame( cluster_summary )
panel_z_mean_columns <- paste0( "mean_panel_z_", names( marker_panels_present ) )
cluster_panel_matrix <- as.matrix( cluster_summary[ , panel_z_mean_columns, drop = FALSE ] )
rownames( cluster_panel_matrix ) <- cluster_summary$cluster
colnames( cluster_panel_matrix ) <- names( marker_panels_present )
top_panel_index <- apply( cluster_panel_matrix, 1L, which.max )
cluster_summary$observed_top_panel <- colnames( cluster_panel_matrix )[ top_panel_index ]
sorted_panel_scores <- t( apply( cluster_panel_matrix, 1L, function(values) {
      sort( values, decreasing = TRUE )
    }
  )
)
cluster_summary$top_panel_score <- sorted_panel_scores[ , 1L ]
cluster_summary$second_panel_score <- sorted_panel_scores[ , 2L ]
cluster_summary$top_panel_margin <- cluster_summary$top_panel_score - cluster_summary$second_panel_score
cluster_summary$expected_panel_is_top <- cluster_summary$expected_primary_panel == cluster_summary$observed_top_panel
cluster_summary$library_dominance_review <- cluster_summary$n_cells >= MIN_CLUSTER_SIZE_FOR_DOMINANCE_REVIEW & cluster_summary$maximum_library_fraction >= LIBRARY_DOMINANCE_REVIEW_THRESHOLD
cluster_summary$timepoint_dominance_review <- cluster_summary$n_cells >= MIN_CLUSTER_SIZE_FOR_DOMINANCE_REVIEW & cluster_summary$maximum_timepoint_fraction >= TIMEPOINT_DOMINANCE_REVIEW_THRESHOLD
cluster_summary <- merge( cluster_summary, as.data.frame( top_marker_strings ), by = "cluster", all.x = TRUE, sort = FALSE )
cluster_summary <- as.data.frame( cluster_summary )
cluster_summary <- cluster_summary[ match( EXPECTED_CLUSTERS, cluster_summary$cluster ), , drop = FALSE ]
write_csv_safe( cluster_summary, "05_myeloid_cluster_annotation_evidence.csv" )
obj$myeloid_working_label <- factor( obj$myeloid_working_label, levels = unname( working_label_map ) )
obj$resident_lineage_eligibility <- factor( obj$resident_lineage_eligibility, levels = c( "LIKELY_RESIDENT_LINEAGE", "AMBIGUOUS_RESIDENT_VS_INFILTRATING" ) )
p_cluster <- Seurat::DimPlot( object = obj, reduction = FROZEN_UMAP, group.by = FROZEN_CLUSTER_COLUMN, label = TRUE, repel = TRUE, raster = TRUE ) + ggplot2::labs( title = "Frozen myeloid Harmony clusters", subtitle = "myhcon_res_0_6; resolution 0.6" ) + Seurat::NoLegend()
p_working_label <- Seurat::DimPlot( object = obj, reduction = FROZEN_UMAP, group.by = "myeloid_working_label", raster = TRUE ) + ggplot2::labs( title = "Step 10E working labels", subtitle = "Review-only; not frozen" )
p_eligibility <- Seurat::DimPlot( object = obj, reduction = FROZEN_UMAP, group.by = "resident_lineage_eligibility", raster = TRUE ) + ggplot2::labs( title = "Resident-lineage eligibility", subtitle = "Ambiguous clusters remain separate" )
p_timepoint <- Seurat::DimPlot( object = obj, reduction = FROZEN_UMAP, group.by = "timepoint", raster = TRUE ) + ggplot2::labs( title = "Myeloid distribution by timepoint" )
ggplot2::ggsave( filename = file.path( figures_dir, "01_myeloid_working_annotation_UMAP.pdf" ), plot = ( p_cluster | p_working_label ) / ( p_eligibility | p_timepoint ), width = 19, height = 14, limitsize = FALSE )
umap_embeddings <- Seurat::Embeddings( object = obj, reduction = FROZEN_UMAP )
score_plot_data <- data.frame( UMAP_1 = umap_embeddings[ , 1L ], UMAP_2 = umap_embeddings[ , 2L ], obj@meta.data[ , score_columns_z, drop = FALSE ], check.names = FALSE, stringsAsFactors = FALSE )
score_plot_names <- c( "Resident_homeostatic", "Activated_DAM_phagolysosomal", "Macrophage_monocyte_like", "BAM", "MHCII_antigen_presentation", "Interferon_response", "Inflammatory_chemokine", "Cell_cycle" )
score_plots <- lapply( score_plot_names, function(panel_name) {
    plot_score_umap( plot_data = score_plot_data, score_column = paste0( "panel_z_", panel_name ), title_text = panel_name )
  }
)
ggplot2::ggsave( filename = file.path( figures_dir, "02_myeloid_panel_score_UMAPs.pdf" ), plot = patchwork::wrap_plots( score_plots, ncol = 4 ), width = 20, height = 10, limitsize = FALSE )
marker_plot_panels <- deduplicate_plot_panels( marker_panels_present )
marker_review_pdf <- file.path( figures_dir, "03_myeloid_marker_DotPlot_and_panel_heatmap.pdf" )
grDevices::pdf( marker_review_pdf, width = 21, height = 10, onefile = TRUE )
print( Seurat::DotPlot( object = obj, features = marker_plot_panels, assay = "RNA", group.by = FROZEN_CLUSTER_COLUMN, cluster.idents = FALSE, dot.scale = 6 ) + Seurat::RotatedAxis() + ggplot2::labs( title = "Marker-first myeloid subtype evidence", x = NULL, y = "Frozen myeloid Harmony cluster" ) + ggplot2::theme( axis.text.x = ggplot2::element_text( size = 7 ) ) )
heatmap_data <- data.frame( cluster = rep( cluster_summary$cluster, each = length( marker_panels_present ) ), panel = rep( names( marker_panels_present ), times = nrow( cluster_summary ) ), mean_z = as.vector( t( cluster_panel_matrix[ cluster_summary$cluster, names( marker_panels_present ), drop = FALSE ] ) ), stringsAsFactors = FALSE )
heatmap_data$cluster <- factor( heatmap_data$cluster, levels = rev( EXPECTED_CLUSTERS ) )
heatmap_data$panel <- factor( heatmap_data$panel, levels = names( marker_panels_present ) )
print( ggplot2::ggplot( heatmap_data, ggplot2::aes( x = panel, y = cluster, fill = mean_z ) ) + ggplot2::geom_tile() + ggplot2::scale_fill_gradient2( midpoint = 0 ) + ggplot2::labs( title = "Cluster mean standardized marker-panel scores", x = NULL, y = "Frozen myeloid Harmony cluster", fill = "Mean Z" ) + ggplot2::theme_bw() + ggplot2::theme( axis.text.x = ggplot2::element_text( angle = 45, hjust = 1 ) ) )
grDevices::dev.off()
obj@misc$GSE234774_step10F <- list( input_object = normalizePath( input_object_path, winslash = "/", mustWork = TRUE ), frozen_reduction = FROZEN_REDUCTION, frozen_umap = FROZEN_UMAP, frozen_cluster_column = FROZEN_CLUSTER_COLUMN, frozen_resolution = 0.6, marker_assay = "RNA", marker_layer = data_layer_name, marker_test = "wilcox", max_cells_per_cluster = MAX_CELLS_PER_CLUSTER_FOR_MARKERS, scoring_method = "sparse mean LogNormalize expression plus across-cell z score", whole_transcriptome_AUCell_run = FALSE, whole_transcriptome_ssGSEA_run = FALSE, working_labels_frozen = FALSE, clean_resident_object_frozen = FALSE, Arm_1_Arm_2_scoring_run = FALSE, completed_at = as.character( Sys.time() ) )
saveRDS( obj, output_object_path, compress = FALSE )
reload_obj <- tryCatch( readRDS( output_object_path ), error = function(error_object) {
    error_object
  }
)
reload_pass <- !inherits( reload_obj, "error" ) && inherits( reload_obj, "Seurat" ) && nrow( reload_obj ) == EXPECTED_FEATURES && ncol( reload_obj ) == EXPECTED_MYLOID_NUCLEI && all( c( "myeloid_working_label", "resident_lineage_eligibility", "myeloid_annotation_status" ) %in% names( reload_obj@meta.data ) ) && all( score_columns_z %in% names( reload_obj@meta.data ) ) && all( required_reductions %in% names( reload_obj@reductions ) )
rm( reload_obj )
invisible( gc() )
required_figure_paths <- file.path( figures_dir, c( "01_myeloid_working_annotation_UMAP.pdf", "02_myeloid_panel_score_UMAPs.pdf", "03_myeloid_marker_DotPlot_and_panel_heatmap.pdf" ) )
n_expected_panel_matches <- sum( cluster_summary$expected_panel_is_top, na.rm = TRUE )
n_ambiguous_clusters <- sum( cluster_summary$resident_lineage_eligibility == "AMBIGUOUS_RESIDENT_VS_INFILTRATING" )
audit_checks <- data.frame( check = c( "Input contains 26,918 features and 25,205 myeloid nuclei", "Frozen Harmony reduction, UMAP, and myhcon_res_0_6 clusters 0-8 are present", "Every curated marker panel contains at least two resolved genes", "Sparse per-cell panel scores were mapped to every nucleus", "Positive markers are available for all nine clusters", "Cluster annotation evidence table contains all nine clusters", "Working labels and resident-lineage eligibility were mapped to every nucleus", "All three required review PDFs were created", "Saved evidence object reloads with scores, labels, and reductions", "Whole-transcriptome dense AUCell/ssGSEA was run", "Working myeloid subtype labels have been frozen", "A clean resident-microglia object has been frozen", "GSE162610 Arm_1/Arm_2 scoring was run" ), required_for_pass = c( TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE ), status = c( "PASS", "PASS", ifelse( all( panel_availability$present_n >= 2L ), "PASS", "STOP" ), ifelse( all( score_columns_z %in% names( obj@meta.data ) ) && !anyNA( obj@meta.data[ , score_columns_z, drop = FALSE ] ), "PASS", "STOP" ), ifelse( setequal( marker_cluster_coverage, EXPECTED_CLUSTERS ), "PASS", "STOP" ), ifelse( nrow( cluster_summary ) == length( EXPECTED_CLUSTERS ), "PASS", "STOP" ), ifelse( !anyNA( obj$myeloid_working_label ) && !anyNA( obj$resident_lineage_eligibility ), "PASS", "STOP" ), ifelse( all( file.exists( required_figure_paths ) ), "PASS", "STOP" ), ifelse( reload_pass, "PASS", "STOP" ), "CONFIRMED_NOT_RUN", "NOT_FROZEN", "NOT_FROZEN", "CONFIRMED_NOT_RUN" ), detail = c( paste0( nrow( obj ), " features x ", ncol( obj ), " nuclei" ), paste( FROZEN_REDUCTION, FROZEN_UMAP, FROZEN_CLUSTER_COLUMN, sep = "; " ), paste0( "minimum present genes per panel=", min( panel_availability$present_n ) ), paste0( length( score_columns_z ), " z-score columns; all nuclei scored" ), paste0( "marker clusters=", length( marker_cluster_coverage ), "; reused completed markers=", reuse_marker_table ), paste0( "evidence rows=", nrow( cluster_summary ), "; expected panel ranked first in ", n_expected_panel_matches, "/9 clusters" ), paste0( "likely resident clusters=", paste( cluster_summary$cluster[ cluster_summary$resident_lineage_eligibility == "LIKELY_RESIDENT_LINEAGE" ], collapse = "," ), "; ambiguous clusters=", paste( cluster_summary$cluster[ cluster_summary$resident_lineage_eligibility == "AMBIGUOUS_RESIDENT_VS_INFILTRATING" ], collapse = "," ) ), paste( basename( required_figure_paths ), collapse = "; " ), ifelse( reload_pass, "RDS reload PASS", "RDS reload failed" ), "Sparse panel means were used; full-transcriptome dense scoring was intentionally deferred.", "Manual review of seven outputs remains required.", paste0( "Ambiguous clusters requiring explicit decision=", n_ambiguous_clusters ), "Arm_1/Arm_2 remain frozen from GSE162610." ), stringsAsFactors = FALSE )
write_csv_safe( audit_checks, "09_STEP10F_audit_checks.csv" )
required_failures <- audit_checks$required_for_pass & audit_checks$status != "PASS"
final_decision <- if ( any( required_failures ) ) {
  "STOP"
} else {
  "PASS_WITH_MYELOID_ANNOTATION_REVIEW"
}
cluster_lines <- apply( cluster_summary, 1L, function(row) {
    paste0( "cluster ", row[["cluster"]], ": ", row[["working_label"]], "; n=", row[["n_cells"]], "; eligibility=", row[["resident_lineage_eligibility"]], "; expected panel=", row[["expected_primary_panel"]], "; observed top panel=", row[["observed_top_panel"]], "; expected_is_top=", row[["expected_panel_is_top"]], "; dominant timepoint=", row[["dominant_timepoint"]], " (", round( as.numeric( row[["maximum_timepoint_fraction"]] ), 3 ), "); dominant library fraction=", round( as.numeric( row[["maximum_library_fraction"]] ), 3 ), "; top markers=", row[["top10_markers"]] )
  }
)
if ( final_decision == "STOP" ) {
  failed <- audit_checks[ audit_checks$required_for_pass & audit_checks$status != "PASS", , drop = FALSE ]
  for ( index in seq_len( nrow( failed ) ) ) {
  }
}
report_path <- file.path( output_dir, "10_STEP10F_FINAL_REPORT.txt" )
if ( final_decision == "STOP" ) {
  stop( "GSE234774 Step 10F decision is STOP. Review the final report.", call. = FALSE )
}
