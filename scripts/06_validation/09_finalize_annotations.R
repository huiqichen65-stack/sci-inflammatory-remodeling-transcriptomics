options( stringsAsFactors = FALSE, warn = 1, Seurat.object.assay.version = "v5", future.globals.maxSize = 16 * 1024^3 )
FROZEN_REDUCTION <- "harmony_myeloid_conservative"
FROZEN_UMAP <- "umap_myeloid_harmony"
FROZEN_CLUSTER_COLUMN <- "myhcon_res_0_6"
EXPECTED_FEATURES <- 26918L
EXPECTED_MYLOID_NUCLEI <- 25205L
EXPECTED_CLUSTERS <- as.character(0:8)
PRIMARY_CLUSTERS <- c("0", "1", "2", "3", "7")
NONIFN_CLUSTERS <- c("0", "1", "2", "3")
CYCLING_CLUSTERS <- c("6")
EXCLUDED_CLUSTERS <- c("4", "5", "8")
COMPLETE_RESIDENT_LINEAGE_CLUSTERS <- c( PRIMARY_CLUSTERS, CYCLING_CLUSTERS )
EXPECTED_CLUSTER_COUNTS <- c( "0" = 6871L, "1" = 5692L, "2" = 3651L, "3" = 3321L, "4" = 2034L, "5" = 1963L, "6" = 650L, "7" = 644L, "8" = 379L )
EXPECTED_OBJECT_COUNTS <- c( FULL_ANNOTATED_MYELOID = 25205L, PRIMARY_NONCYCLING_RESIDENT = 20179L, NONIFN_SENSITIVITY_RESIDENT = 19535L, CYCLING_RESIDENT = 650L, EXCLUDED_NONMICROGLIAL_OR_CONTAMINANT = 4376L )
expected_timepoints <- c( "Uninjured", "1 dpi", "4 dpi", "7 dpi", "14 dpi", "1 month", "2 months" )
final_label_map <- c( "0" = "Homeostatic_microglia", "1" = "Activated_transitional_microglia", "2" = "Phagocytic_DAM_like_microglia", "3" = "Resident_associated_inflammatory_intermediate_microglia", "4" = "Mrc1_MHCII_BAM_like_macrophage", "5" = "Arg1_F7_macrophage_like_non_microglial_myeloid", "6" = "Cycling_microglia", "7" = "IFN_responsive_microglia", "8" = "Neuronal_contaminant_residual_doublet_like" )
resident_status_map <- c( "0" = "RESIDENT_LINEAGE_INCLUDED", "1" = "RESIDENT_LINEAGE_INCLUDED", "2" = "RESIDENT_LINEAGE_INCLUDED", "3" = "RESIDENT_LINEAGE_INCLUDED", "4" = "NON_MICROGLIAL_MYELOID_EXCLUDED", "5" = "NON_MICROGLIAL_MYELOID_EXCLUDED", "6" = "RESIDENT_LINEAGE_INCLUDED", "7" = "RESIDENT_LINEAGE_INCLUDED", "8" = "CONTAMINANT_EXCLUDED" )
arm_analysis_set_map <- c( "0" = "PRIMARY_NONCYCLING_RESIDENT", "1" = "PRIMARY_NONCYCLING_RESIDENT", "2" = "PRIMARY_NONCYCLING_RESIDENT", "3" = "PRIMARY_NONCYCLING_RESIDENT", "4" = "EXCLUDED_NONMICROGLIAL_OR_CONTAMINANT", "5" = "EXCLUDED_NONMICROGLIAL_OR_CONTAMINANT", "6" = "CYCLING_RESIDENT_SEPARATE", "7" = "PRIMARY_NONCYCLING_RESIDENT", "8" = "EXCLUDED_NONMICROGLIAL_OR_CONTAMINANT" )
nonIFN_sensitivity_map <- c( "0" = "INCLUDED_NONIFN_SENSITIVITY", "1" = "INCLUDED_NONIFN_SENSITIVITY", "2" = "INCLUDED_NONIFN_SENSITIVITY", "3" = "INCLUDED_NONIFN_SENSITIVITY", "4" = "NOT_ELIGIBLE", "5" = "NOT_ELIGIBLE", "6" = "NOT_ELIGIBLE", "7" = "EXCLUDED_IFN_SENSITIVITY", "8" = "NOT_ELIGIBLE" )
for (mapping in list( final_label_map, resident_status_map, arm_analysis_set_map, nonIFN_sensitivity_map )) {
  if (!identical( names(mapping), EXPECTED_CLUSTERS )) {
    stop( "At least one frozen cluster map does not match clusters 0-8.", call. = FALSE )
  }
}
project_root <- Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
input_object_path <- file.path( project_root, "04_results", "step10F_GSE234774_myeloid_marker_first_annotation_evidence", "RDS", "GSE234774_step10F_myeloid_annotation_evidence_BPCells.rds" )
output_dir <- file.path( project_root, "04_results", "step10G_GSE234774_frozen_myeloid_annotations_resident_objects" )
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
object_paths <- c( FULL_ANNOTATED_MYELOID = file.path( rds_dir, "GSE234774_step10G_full_myeloid_final_annotations_BPCells.rds" ), PRIMARY_NONCYCLING_RESIDENT = file.path( rds_dir, "GSE234774_step10G_primary_noncycling_resident_microglia_BPCells.rds" ), NONIFN_SENSITIVITY_RESIDENT = file.path( rds_dir, "GSE234774_step10G_nonIFN_sensitivity_resident_microglia_BPCells.rds" ), CYCLING_RESIDENT = file.path( rds_dir, "GSE234774_step10G_cycling_resident_microglia_BPCells.rds" ), EXCLUDED_NONMICROGLIAL_OR_CONTAMINANT = file.path( rds_dir, "GSE234774_step10G_excluded_nonmicroglial_contaminant_BPCells.rds" ) )
log_message <- function(...) {
  message_text <- paste0(...)
  timestamp <- format( Sys.time(), "%Y-%m-%d %H:%M:%S" )
  line <- paste0( "[", timestamp, "] ", message_text )
  cat( line, "\n" )
  cat( line, "\n", file = run_log, append = TRUE )
}
write_csv_safe <- function( x, filename ) {
  path <- file.path( tables_dir, filename )
  utils::write.csv( x, path, row.names = FALSE, na = "" )
  log_message( "Wrote: ", path )
  invisible( path )
}
stop_with_report <- function(reason) {
  report_path <- file.path( output_dir, "10_STEP10G_FINAL_REPORT.txt" )
  report_lines <- c( "GSE234774 STEP 10G FROZEN MYELOID ANNOTATIONS / RESIDENT OBJECTS", "============================================================", paste0( "Run time: ", format( Sys.time(), "%Y-%m-%d %H:%M:%S" ) ), paste0( "Output directory: ", output_dir ), "", "FINAL DECISION: STOP", "", paste0( "Reason: ", reason ), "", "Do not start Arm_1/Arm_2 projection." )
  writeLines( report_lines, report_path, useBytes = TRUE )
  stop( reason, call. = FALSE )
}
safe_package_version <- function( package_name ) {
  if (!requireNamespace( package_name, quietly = TRUE )) {
    return( NA_character_ )
  }
  as.character( utils::packageVersion( package_name ) )
}
safe_valid_object <- function(object) {
  tryCatch( {
      methods::validObject( object )
      TRUE
    },
    error = function(error_object) {
      FALSE
    }
  )
}
make_subset_object <- function( object, cells, object_name, expected_nuclei ) {
  if ( length( cells ) != expected_nuclei ) {
    stop( paste0( object_name, " cell-vector count mismatch before subsetting: expected ", expected_nuclei, ", detected ", length( cells ), "." ), call. = FALSE )
  }
  subset_object <- subset( x = object, cells = cells )
  if ( nrow( subset_object ) != EXPECTED_FEATURES || ncol( subset_object ) != expected_nuclei ) {
    stop( paste0( object_name, " dimensions are incorrect after subsetting: ", nrow( subset_object ), " x ", ncol( subset_object ), "." ), call. = FALSE )
  }
  if (!identical( colnames( subset_object ), cells )) {
    stop( paste0( object_name, " barcode order is incorrect after subsetting." ), call. = FALSE )
  }
  required_subset_metadata <- c( "myeloid_final_label", "resident_lineage_status", "arm_analysis_set", "nonIFN_sensitivity_status", "step10G_annotation_status", FROZEN_CLUSTER_COLUMN )
  if (!all( required_subset_metadata %in% names( subset_object@meta.data ) )) {
    stop( paste0( object_name, " is missing frozen annotation metadata." ), call. = FALSE )
  }
  if (!all( c( FROZEN_REDUCTION, FROZEN_UMAP ) %in% names( subset_object@reductions ) )) {
    stop( paste0( object_name, " is missing the frozen Harmony reduction or UMAP." ), call. = FALSE )
  }
  if (!safe_valid_object( subset_object )) {
    stop( paste0( object_name, " failed Seurat object validity checking." ), call. = FALSE )
  }
  subset_object
}
sample_layer_values <- function( object, layer_name, features, cells ) {
  SeuratObject::LayerData( object = object, assay = "RNA", layer = layer_name, features = features, cells = cells, fast = FALSE )
}
compare_matrix_values <- function( matrix_a, matrix_b, tolerance ) {
  if (!identical( dim( matrix_a ), dim( matrix_b ) )) {
    return( FALSE )
  }
  if (!identical( rownames( matrix_a ), rownames( matrix_b ) )) {
    return( FALSE )
  }
  if (!identical( colnames( matrix_a ), colnames( matrix_b ) )) {
    return( FALSE )
  }
  isTRUE( all.equal( as.matrix( matrix_a ), as.matrix( matrix_b ), tolerance = tolerance, check.attributes = TRUE ) )
}
deduplicate_marker_panels <- function( panel_list ) {
  used <- character(0)
  output <- list()
  for (panel_name in names( panel_list )) {
    genes <- panel_list[[panel_name]]
    genes <- genes[ !genes %in% used ]
    if ( length( genes ) > 0L ) {
      output[[panel_name]] <- genes
      used <- c( used, genes )
    }
  }
  output
}
resolve_features <- function( requested_features, available_features ) {
  available_upper <- toupper( available_features )
  first_index <- !duplicated( available_upper )
  lookup <- setNames( available_features[ first_index ], available_upper[ first_index ] )
  resolved <- unname( lookup[ toupper( requested_features ) ] )
  names( resolved ) <- requested_features
  resolved
}
log_message( "Starting package preflight." )
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
  stop_with_report( paste0( "Input Step 10F object not found: ", input_object_path ) )
}
log_message( "Loading Step 10F myeloid annotation-evidence object." )
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
if (!all( c( FROZEN_REDUCTION, FROZEN_UMAP ) %in% names( obj@reductions ) )) {
  stop_with_report( "Frozen myeloid Harmony reduction or UMAP is absent." )
}
cluster_values <- as.character( obj@meta.data[ , FROZEN_CLUSTER_COLUMN ] )
observed_clusters <- sort( unique( cluster_values ), method = "radix" )
if (!setequal( observed_clusters, EXPECTED_CLUSTERS )) {
  stop_with_report( paste0( "Expected clusters 0-8; detected ", paste( observed_clusters, collapse = ", " ), "." ) )
}
observed_cluster_counts <- table( factor( cluster_values, levels = EXPECTED_CLUSTERS ) )
observed_cluster_counts <- as.integer( observed_cluster_counts )
names( observed_cluster_counts ) <- EXPECTED_CLUSTERS
if (!identical( observed_cluster_counts, EXPECTED_CLUSTER_COUNTS )) {
  stop_with_report( paste0( "Frozen cluster counts do not match the reviewed Step 10F counts. ", "Observed=", paste( paste0( names( observed_cluster_counts ), ":", observed_cluster_counts ), collapse = ";" ), "." ) )
}
Seurat::DefaultAssay( obj ) <- "RNA"
Seurat::Idents( obj ) <- FROZEN_CLUSTER_COLUMN
rna_layers <- SeuratObject::Layers( object = obj, assay = "RNA", search = NA )
counts_layers <- rna_layers[ grepl( "^counts($|\\.)", rna_layers ) ]
data_layers <- rna_layers[ grepl( "^data($|\\.)", rna_layers ) ]
if ( length( counts_layers ) != 1L || length( data_layers ) != 1L ) {
  stop_with_report( paste0( "Expected one RNA counts layer and one RNA data layer; detected ", length( counts_layers ), " and ", length( data_layers ), "." ) )
}
counts_layer_name <- counts_layers[[1L]]
data_layer_name <- data_layers[[1L]]
set.seed( 20260803L )
integrity_features <- rownames( obj )[ sort( sample( seq_len( nrow( obj ) ), size = 30L, replace = FALSE ) ) ]
integrity_cells <- colnames( obj )[ sort( sample( seq_len( ncol( obj ) ), size = 120L, replace = FALSE ) ) ]
input_counts_sample <- sample_layer_values( object = obj, layer_name = counts_layer_name, features = integrity_features, cells = integrity_cells )
input_data_sample <- sample_layer_values( object = obj, layer_name = data_layer_name, features = integrity_features, cells = integrity_cells )
input_audit <- data.frame( item = c( "Input object", "Features", "Myeloid nuclei", "Frozen reduction", "Frozen UMAP", "Frozen cluster column", "Frozen cluster count", "RNA counts layer", "RNA data layer" ), value = c( normalizePath( input_object_path, winslash = "/", mustWork = TRUE ), nrow( obj ), ncol( obj ), FROZEN_REDUCTION, FROZEN_UMAP, FROZEN_CLUSTER_COLUMN, length( observed_clusters ), counts_layer_name, data_layer_name ), stringsAsFactors = FALSE )
write_csv_safe( input_audit, "01_input_object_audit.csv" )
log_message( "Freezing reviewed cluster-level annotations." )
obj$myeloid_final_label <- unname( final_label_map[ cluster_values ] )
obj$resident_lineage_status <- unname( resident_status_map[ cluster_values ] )
obj$arm_analysis_set <- unname( arm_analysis_set_map[ cluster_values ] )
obj$nonIFN_sensitivity_status <- unname( nonIFN_sensitivity_map[ cluster_values ] )
obj$step10G_annotation_status <- "FROZEN_CLUSTER_LEVEL_ANNOTATION"
obj$primary_noncycling_resident_member <- cluster_values %in%
  PRIMARY_CLUSTERS
obj$nonIFN_sensitivity_member <- cluster_values %in%
  NONIFN_CLUSTERS
obj$cycling_resident_member <- cluster_values %in%
  CYCLING_CLUSTERS
obj$excluded_nonmicroglial_or_contaminant_member <- cluster_values %in%
  EXCLUDED_CLUSTERS
obj$complete_resident_lineage_member <- cluster_values %in%
  COMPLETE_RESIDENT_LINEAGE_CLUSTERS
if ( anyNA( obj$myeloid_final_label ) || anyNA( obj$resident_lineage_status ) || anyNA( obj$arm_analysis_set ) || anyNA( obj$nonIFN_sensitivity_status ) ) {
  stop_with_report( "Frozen annotations were not mapped to every nucleus." )
}
obj$myeloid_final_label <- factor( obj$myeloid_final_label, levels = unname( final_label_map ) )
obj$resident_lineage_status <- factor( obj$resident_lineage_status, levels = c( "RESIDENT_LINEAGE_INCLUDED", "NON_MICROGLIAL_MYELOID_EXCLUDED", "CONTAMINANT_EXCLUDED" ) )
obj$arm_analysis_set <- factor( obj$arm_analysis_set, levels = c( "PRIMARY_NONCYCLING_RESIDENT", "CYCLING_RESIDENT_SEPARATE", "EXCLUDED_NONMICROGLIAL_OR_CONTAMINANT" ) )
obj$nonIFN_sensitivity_status <- factor( obj$nonIFN_sensitivity_status, levels = c( "INCLUDED_NONIFN_SENSITIVITY", "EXCLUDED_IFN_SENSITIVITY", "NOT_ELIGIBLE" ) )
exclusive_membership_count <- ( as.integer( obj$primary_noncycling_resident_member ) + as.integer( obj$cycling_resident_member ) + as.integer( obj$excluded_nonmicroglial_or_contaminant_member ) )
exclusive_partition_pass <- all( exclusive_membership_count == 1L )
if (!exclusive_partition_pass) {
  stop_with_report( "Primary, cycling, and excluded membership do not form an exact one-way partition." )
}
nonIFN_subset_pass <- all( !obj$nonIFN_sensitivity_member | obj$primary_noncycling_resident_member )
if (!nonIFN_subset_pass) {
  stop_with_report( "The non-IFN sensitivity set is not a strict subset of the primary resident set." )
}
observed_object_counts <- c( FULL_ANNOTATED_MYELOID = ncol( obj ), PRIMARY_NONCYCLING_RESIDENT = sum( obj$primary_noncycling_resident_member ), NONIFN_SENSITIVITY_RESIDENT = sum( obj$nonIFN_sensitivity_member ), CYCLING_RESIDENT = sum( obj$cycling_resident_member ), EXCLUDED_NONMICROGLIAL_OR_CONTAMINANT = sum( obj$excluded_nonmicroglial_or_contaminant_member ) )
if (!identical( as.integer( observed_object_counts ), as.integer( EXPECTED_OBJECT_COUNTS[ names( observed_object_counts ) ] ) )) {
  stop_with_report( paste0( "Object membership counts do not match expectations. Observed=", paste( paste0( names( observed_object_counts ), ":", observed_object_counts ), collapse = ";" ), "." ) )
}
primary_cells <- colnames( obj )[ obj$primary_noncycling_resident_member ]
nonIFN_cells <- colnames( obj )[ obj$nonIFN_sensitivity_member ]
cycling_cells <- colnames( obj )[ obj$cycling_resident_member ]
excluded_cells <- colnames( obj )[ obj$excluded_nonmicroglial_or_contaminant_member ]
source_library_manifest <- unique( data.table::data.table( timepoint = as.character( obj$timepoint ), biological_library = as.character( obj$biological_library ) ) )
source_library_manifest <- source_library_manifest[ !is.na( timepoint ) & nzchar( timepoint ) & !is.na( biological_library ) & nzchar( biological_library ) ]
all_source_libraries <- sort( unique( source_library_manifest$biological_library ) )
all_source_timepoints <- unique( source_library_manifest$timepoint )
if ( length( all_source_libraries ) == 0L || length( all_source_timepoints ) == 0L ) {
  stop_with_report( "The Step 10F myeloid source cohort has no valid library/timepoint metadata." )
}
unexpected_source_timepoints <- setdiff( all_source_timepoints, expected_timepoints )
if ( length( unexpected_source_timepoints ) > 0L ) {
  stop_with_report( paste0( "Unexpected timepoint label(s) in the Step 10F source cohort: ", paste( unexpected_source_timepoints, collapse = ", " ), "." ) )
}
library_timepoint_mapping <- source_library_manifest[ , .( n_timepoint_labels = uniqueN( timepoint ), timepoint_labels = paste( sort( unique( timepoint ) ), collapse = ";" ) ), by = biological_library ]
if ( any( library_timepoint_mapping$n_timepoint_labels != 1L ) ) {
  problematic_libraries <- library_timepoint_mapping$biological_library[ library_timepoint_mapping$n_timepoint_labels != 1L ]
  stop_with_report( paste0( "At least one biological_library maps to multiple timepoints: ", paste( problematic_libraries, collapse = ", " ), "." ) )
}
missing_source_timepoints <- setdiff( expected_timepoints, all_source_timepoints )
log_message( "Step 10F myeloid source coverage: ", length( all_source_libraries ), " libraries and ", length( all_source_timepoints ), "/", length( expected_timepoints ), " expected timepoints." )
if ( length( missing_source_timepoints ) > 0L ) {
  log_message( "Source-coverage review: expected timepoints absent from the Step 10F myeloid compartment: ", paste( missing_source_timepoints, collapse = ", " ), "." )
}
primary_coverage_counts <- data.table::data.table( timepoint = as.character( obj$timepoint[ obj$primary_noncycling_resident_member ] ), biological_library = as.character( obj$biological_library[ obj$primary_noncycling_resident_member ] ) )[ , .( primary_n_cells = .N ), by = .( timepoint, biological_library ) ]
primary_library_coverage <- merge( source_library_manifest, primary_coverage_counts, by = c( "timepoint", "biological_library" ), all.x = TRUE, sort = FALSE )
primary_library_coverage$primary_n_cells[ is.na( primary_library_coverage$primary_n_cells ) ] <- 0L
primary_library_coverage$retained_in_primary <- primary_library_coverage$primary_n_cells > 0L
primary_library_coverage$timepoint <- factor( primary_library_coverage$timepoint, levels = expected_timepoints, ordered = TRUE )
data.table::setorder( primary_library_coverage, timepoint, biological_library )
primary_timepoint_coverage <- primary_library_coverage[ , .( source_libraries = .N, retained_libraries = sum( retained_in_primary ), missing_libraries = sum( !retained_in_primary ), primary_n_cells = sum( primary_n_cells ), missing_library_names = paste( biological_library[ !retained_in_primary ], collapse = ";" ) ), by = timepoint ]
primary_present_libraries <- sort( primary_library_coverage$biological_library[ primary_library_coverage$retained_in_primary ] )
primary_present_libraries <- unique( primary_present_libraries )
missing_primary_libraries <- setdiff( all_source_libraries, primary_present_libraries )
primary_present_timepoints <- as.character( primary_timepoint_coverage$timepoint[ primary_timepoint_coverage$primary_n_cells > 0L ] )
missing_primary_from_source_timepoints <- setdiff( all_source_timepoints, primary_present_timepoints )
missing_primary_expected_timepoints <- setdiff( expected_timepoints, primary_present_timepoints )
primary_all_libraries_retained <- length( missing_primary_libraries ) == 0L
primary_all_source_timepoints_retained <- length( missing_primary_from_source_timepoints ) == 0L
source_has_all_expected_timepoints <- length( missing_source_timepoints ) == 0L
low_replication_timepoints <- as.character( primary_timepoint_coverage$timepoint[ primary_timepoint_coverage$retained_libraries < 2L ] )
primary_coverage_count_pass <- sum( primary_library_coverage$primary_n_cells ) == length( primary_cells )
if (!primary_coverage_count_pass) {
  stop_with_report( "Primary library/timepoint coverage counts do not sum to the frozen primary-object count." )
}
source_coverage_summary_rows <- rbind( data.frame( record_type = "SOURCE_COHORT_SUMMARY", timepoint = "ALL", biological_library = "ALL", source_libraries = length( all_source_libraries ), retained_libraries = length( primary_present_libraries ), missing_libraries = length( missing_primary_libraries ), primary_n_cells = length( primary_cells ), retained_in_primary = NA, missing_library_names = paste( missing_primary_libraries, collapse = ";" ), stringsAsFactors = FALSE ), data.frame( record_type = "SOURCE_EXPECTED_TIMEPOINT_COVERAGE", timepoint = paste( expected_timepoints, collapse = ";" ), biological_library = "ALL", source_libraries = length( all_source_libraries ), retained_libraries = NA_integer_, missing_libraries = length( missing_source_timepoints ), primary_n_cells = length( primary_cells ), retained_in_primary = source_has_all_expected_timepoints, missing_library_names = paste( missing_source_timepoints, collapse = ";" ), stringsAsFactors = FALSE ) )
primary_coverage_precheck <- rbind( source_coverage_summary_rows, data.frame( record_type = "SOURCE_LIBRARY", timepoint = as.character( primary_library_coverage$timepoint ), biological_library = primary_library_coverage$biological_library, source_libraries = NA_integer_, retained_libraries = NA_integer_, missing_libraries = NA_integer_, primary_n_cells = primary_library_coverage$primary_n_cells, retained_in_primary = primary_library_coverage$retained_in_primary, missing_library_names = "", stringsAsFactors = FALSE ), data.frame( record_type = "TIMEPOINT_SUMMARY", timepoint = as.character( primary_timepoint_coverage$timepoint ), biological_library = "ALL", source_libraries = primary_timepoint_coverage$source_libraries, retained_libraries = primary_timepoint_coverage$retained_libraries, missing_libraries = primary_timepoint_coverage$missing_libraries, primary_n_cells = primary_timepoint_coverage$primary_n_cells, retained_in_primary = primary_timepoint_coverage$primary_n_cells > 0L, missing_library_names = primary_timepoint_coverage$missing_library_names, stringsAsFactors = FALSE ) )
write_csv_safe( primary_coverage_precheck, "02_primary_resident_library_timepoint_coverage_precheck.csv" )
log_message( "Primary resident coverage: ", length( primary_present_libraries ), "/", length( all_source_libraries ), " source libraries retained; ", length( primary_present_timepoints ), "/", length( expected_timepoints ), " project timepoints retained in the primary object." )
if ( length( missing_primary_libraries ) > 0L ) {
  log_message( "Coverage review: source libraries with zero primary resident nuclei: ", paste( missing_primary_libraries, collapse = ", " ), "." )
}
if ( length( missing_primary_expected_timepoints ) > 0L ) {
  log_message( "Coverage review: expected project timepoints absent from the primary resident object: ", paste( missing_primary_expected_timepoints, collapse = ", " ), "." )
}
if ( length( missing_primary_from_source_timepoints ) > 0L ) {
  log_message( "Coverage review: Step 10F source timepoints lost after primary resident filtering: ", paste( missing_primary_from_source_timepoints, collapse = ", " ), "." )
}
if ( length( low_replication_timepoints ) > 0L ) {
  log_message( "Coverage review: timepoints represented by fewer than two retained libraries: ", paste( low_replication_timepoints, collapse = ", " ), "." )
}
log_message( "Creating Step 10G annotation and object-membership figures." )
p_cluster <- Seurat::DimPlot( object = obj, reduction = FROZEN_UMAP, group.by = FROZEN_CLUSTER_COLUMN, label = TRUE, repel = TRUE, raster = TRUE ) + ggplot2::labs( title = "Frozen myeloid Harmony clusters", subtitle = "myhcon_res_0_6; resolution 0.6" ) + Seurat::NoLegend()
p_label <- Seurat::DimPlot( object = obj, reduction = FROZEN_UMAP, group.by = "myeloid_final_label", raster = TRUE ) + ggplot2::labs( title = "Frozen final myeloid annotation" )
p_analysis_set <- Seurat::DimPlot( object = obj, reduction = FROZEN_UMAP, group.by = "arm_analysis_set", raster = TRUE ) + ggplot2::labs( title = "Frozen Arm-analysis membership" )
p_timepoint <- Seurat::DimPlot( object = obj, reduction = FROZEN_UMAP, group.by = "timepoint", raster = TRUE ) + ggplot2::labs( title = "Timepoint distribution after annotation freeze" )
ggplot2::ggsave( filename = file.path( figures_dir, "01_frozen_annotation_and_analysis_set_UMAP.pdf" ), plot = ( p_cluster | p_label ) / ( p_analysis_set | p_timepoint ), width = 20, height = 15, limitsize = FALSE )
composition_data <- data.table::as.data.table( obj@meta.data )[ , .( n_cells = .N ), by = .( arm_analysis_set, timepoint ) ]
composition_data[ , within_timepoint_fraction := n_cells / sum( n_cells ), by = timepoint ]
p_counts <- ggplot2::ggplot( composition_data, ggplot2::aes( x = timepoint, y = n_cells, fill = arm_analysis_set ) ) + ggplot2::geom_col() + ggplot2::scale_y_continuous( labels = scales::comma ) + ggplot2::labs( title = "Analysis-set nuclei by timepoint", x = NULL, y = "Nuclei", fill = "Analysis set" ) + ggplot2::theme_bw() + ggplot2::theme( axis.text.x = ggplot2::element_text( angle = 45, hjust = 1 ) )
p_fraction <- ggplot2::ggplot( composition_data, ggplot2::aes( x = timepoint, y = within_timepoint_fraction, fill = arm_analysis_set ) ) + ggplot2::geom_col() + ggplot2::scale_y_continuous( labels = scales::percent ) + ggplot2::labs( title = "Analysis-set fraction within each timepoint", x = NULL, y = "Fraction", fill = "Analysis set" ) + ggplot2::theme_bw() + ggplot2::theme( axis.text.x = ggplot2::element_text( angle = 45, hjust = 1 ) )
ggplot2::ggsave( filename = file.path( figures_dir, "02_analysis_set_timepoint_composition.pdf" ), plot = p_counts / p_fraction, width = 13, height = 11, limitsize = FALSE )
review_marker_panels_requested <- list( Resident = c( "P2ry12", "Tmem119", "Sall1", "Fcrls", "Hexb", "Cx3cr1", "Siglech" ), Activated_DAM = c( "Apoe", "Cst7", "Trem2", "Itgax", "Lgals3", "Cd68", "Mmp12" ), BAM_Macrophage = c( "Mrc1", "Cd74", "H2-Ab1", "H2-Aa", "Arg1", "F7", "F10" ), Cycling = c( "Mki67", "Top2a", "Kif11", "Bub1", "Kif23" ), IFN = c( "Isg15", "Ifit1", "Ifit3", "Rsad2", "Irf7" ), Neuronal_contamination = c( "Syt1", "Nrcam", "Rims2", "Nxph1", "Kcnd2" ) )
all_review_markers <- unique( unlist( review_marker_panels_requested, use.names = FALSE ) )
resolved_review_markers <- resolve_features( requested_features = all_review_markers, available_features = rownames( obj ) )
review_marker_panels_present <- lapply( review_marker_panels_requested, function(panel_genes) {
    resolved <- resolved_review_markers[ panel_genes ]
    unique( unname( resolved[ !is.na( resolved ) ] ) )
  }
)
review_marker_panels_present <- deduplicate_marker_panels( review_marker_panels_present )
if ( any( vapply( review_marker_panels_present, length, integer(1) ) < 2L ) ) {
  stop_with_report( "At least one final-review marker panel contains fewer than two resolved genes." )
}
p_marker <- Seurat::DotPlot( object = obj, features = review_marker_panels_present, assay = "RNA", group.by = "myeloid_final_label", cluster.idents = FALSE, dot.scale = 6 ) + Seurat::RotatedAxis() + ggplot2::labs( title = "Marker evidence after final annotation freeze", x = NULL, y = "Frozen final myeloid label" ) + ggplot2::theme( axis.text.x = ggplot2::element_text( size = 7 ), axis.text.y = ggplot2::element_text( size = 8 ) )
ggplot2::ggsave( filename = file.path( figures_dir, "03_final_annotation_marker_DotPlot.pdf" ), plot = p_marker, width = 19, height = 10, limitsize = FALSE )
cluster_annotation_table <- data.frame( cluster = EXPECTED_CLUSTERS, expected_n_cells = as.integer( EXPECTED_CLUSTER_COUNTS[ EXPECTED_CLUSTERS ] ), observed_n_cells = as.integer( observed_cluster_counts[ EXPECTED_CLUSTERS ] ), myeloid_final_label = unname( final_label_map[ EXPECTED_CLUSTERS ] ), resident_lineage_status = unname( resident_status_map[ EXPECTED_CLUSTERS ] ), arm_analysis_set = unname( arm_analysis_set_map[ EXPECTED_CLUSTERS ] ), nonIFN_sensitivity_status = unname( nonIFN_sensitivity_map[ EXPECTED_CLUSTERS ] ), primary_noncycling_resident_member = EXPECTED_CLUSTERS %in% PRIMARY_CLUSTERS, nonIFN_sensitivity_member = EXPECTED_CLUSTERS %in% NONIFN_CLUSTERS, cycling_resident_member = EXPECTED_CLUSTERS %in% CYCLING_CLUSTERS, excluded_nonmicroglial_or_contaminant_member = EXPECTED_CLUSTERS %in% EXCLUDED_CLUSTERS, complete_resident_lineage_member = EXPECTED_CLUSTERS %in% COMPLETE_RESIDENT_LINEAGE_CLUSTERS, count_match = as.integer( EXPECTED_CLUSTER_COUNTS[ EXPECTED_CLUSTERS ] ) == as.integer( observed_cluster_counts[ EXPECTED_CLUSTERS ] ), stringsAsFactors = FALSE )
write_csv_safe( cluster_annotation_table, "04_frozen_cluster_annotation_and_object_membership.csv" )
object_cell_sets <- list( FULL_ANNOTATED_MYELOID = colnames( obj ), PRIMARY_NONCYCLING_RESIDENT = primary_cells, NONIFN_SENSITIVITY_RESIDENT = nonIFN_cells, CYCLING_RESIDENT = cycling_cells, EXCLUDED_NONMICROGLIAL_OR_CONTAMINANT = excluded_cells )
object_composition_rows <- list()
for (object_name in names( object_cell_sets )) {
  cells <- object_cell_sets[[object_name]]
  metadata_subset <- obj@meta.data[ cells, , drop = FALSE ]
  overall_row <- data.frame( record_type = "OBJECT_TOTAL", object_name = object_name, timepoint = "ALL", biological_library = "ALL", myeloid_final_label = "ALL", n_cells = nrow( metadata_subset ), object_fraction = 1, stringsAsFactors = FALSE )
  timepoint_rows <- data.table::as.data.table( metadata_subset )[ , .( n_cells = .N ), by = .( timepoint ) ]
  timepoint_rows$record_type <- "TIMEPOINT"
  timepoint_rows$object_name <- object_name
  timepoint_rows$biological_library <- "ALL"
  timepoint_rows$myeloid_final_label <- "ALL"
  timepoint_rows$object_fraction <- timepoint_rows$n_cells / nrow( metadata_subset )
  library_rows <- data.table::as.data.table( metadata_subset )[ , .( n_cells = .N ), by = .( timepoint, biological_library ) ]
  library_rows$record_type <- "TIMEPOINT_LIBRARY"
  library_rows$object_name <- object_name
  library_rows$myeloid_final_label <- "ALL"
  library_rows$object_fraction <- library_rows$n_cells / nrow( metadata_subset )
  label_rows <- data.table::as.data.table( metadata_subset )[ , .( n_cells = .N ), by = .( myeloid_final_label ) ]
  label_rows$record_type <- "FINAL_LABEL"
  label_rows$object_name <- object_name
  label_rows$timepoint <- "ALL"
  label_rows$biological_library <- "ALL"
  label_rows$object_fraction <- label_rows$n_cells / nrow( metadata_subset )
  object_composition_rows[[object_name]] <- rbind( overall_row, as.data.frame( timepoint_rows[ , c( "record_type", "object_name", "timepoint", "biological_library", "myeloid_final_label", "n_cells", "object_fraction" ), with = FALSE ] ), as.data.frame( library_rows[ , c( "record_type", "object_name", "timepoint", "biological_library", "myeloid_final_label", "n_cells", "object_fraction" ), with = FALSE ] ), as.data.frame( label_rows[ , c( "record_type", "object_name", "timepoint", "biological_library", "myeloid_final_label", "n_cells", "object_fraction" ), with = FALSE ] ) )
}
object_composition_table <- do.call( rbind, object_composition_rows )
rownames( object_composition_table ) <- NULL
primary_source_coverage_rows <- data.frame( record_type = "PRIMARY_SOURCE_LIBRARY_COVERAGE", object_name = "PRIMARY_NONCYCLING_RESIDENT", timepoint = as.character( primary_library_coverage$timepoint ), biological_library = primary_library_coverage$biological_library, myeloid_final_label = "ALL", n_cells = primary_library_coverage$primary_n_cells, object_fraction = primary_library_coverage$primary_n_cells / length( primary_cells ), stringsAsFactors = FALSE )
object_composition_table <- rbind( object_composition_table, primary_source_coverage_rows )
write_csv_safe( object_composition_table, "05_frozen_object_composition_summary.csv" )
obj@misc$GSE234774_step10G <- list( input_object = normalizePath( input_object_path, winslash = "/", mustWork = TRUE ), frozen_reduction = FROZEN_REDUCTION, frozen_umap = FROZEN_UMAP, frozen_cluster_column = FROZEN_CLUSTER_COLUMN, frozen_resolution = 0.6, final_label_map = final_label_map, resident_status_map = resident_status_map, arm_analysis_set_map = arm_analysis_set_map, primary_clusters = PRIMARY_CLUSTERS, nonIFN_clusters = NONIFN_CLUSTERS, cycling_clusters = CYCLING_CLUSTERS, excluded_clusters = EXCLUDED_CLUSTERS, primary_noncycling_resident_n = length( primary_cells ), nonIFN_sensitivity_n = length( nonIFN_cells ), cycling_resident_n = length( cycling_cells ), excluded_n = length( excluded_cells ), RNA_layers_modified = FALSE, reprocessing_run = FALSE, Arm_1_Arm_2_scoring_run = FALSE, completed_at = as.character( Sys.time() ) )
log_message( "Saving the full annotated myeloid object." )
saveRDS( obj, object_paths[[ "FULL_ANNOTATED_MYELOID" ]], compress = FALSE )
subset_definitions <- list( PRIMARY_NONCYCLING_RESIDENT = list( cells = primary_cells, expected_n = EXPECTED_OBJECT_COUNTS[[ "PRIMARY_NONCYCLING_RESIDENT" ]], object_role = "PRIMARY_ARM_ANALYSIS_INCLUDES_IFN" ), NONIFN_SENSITIVITY_RESIDENT = list( cells = nonIFN_cells, expected_n = EXPECTED_OBJECT_COUNTS[[ "NONIFN_SENSITIVITY_RESIDENT" ]], object_role = "SENSITIVITY_EXCLUDES_IFN_CLUSTER_7" ), CYCLING_RESIDENT = list( cells = cycling_cells, expected_n = EXPECTED_OBJECT_COUNTS[[ "CYCLING_RESIDENT" ]], object_role = "CYCLING_RESIDENT_SEPARATE_LAYER" ), EXCLUDED_NONMICROGLIAL_OR_CONTAMINANT = list( cells = excluded_cells, expected_n = EXPECTED_OBJECT_COUNTS[[ "EXCLUDED_NONMICROGLIAL_OR_CONTAMINANT" ]], object_role = "EXCLUDED_FROM_RESIDENT_ARM_ANALYSIS" ) )
save_summary_rows <- list()
full_object_valid <- safe_valid_object( obj )
save_summary_rows[[ "FULL_ANNOTATED_MYELOID" ]] <- data.frame( object_name = "FULL_ANNOTATED_MYELOID", expected_n_cells = EXPECTED_OBJECT_COUNTS[[ "FULL_ANNOTATED_MYELOID" ]], observed_n_cells = ncol( obj ), expected_n_features = EXPECTED_FEATURES, observed_n_features = nrow( obj ), valid_object = full_object_valid, object_role = "FROZEN_FULL_MYELOID_ANNOTATION", output_path = object_paths[[ "FULL_ANNOTATED_MYELOID" ]], stringsAsFactors = FALSE )
for (object_name in names( subset_definitions )) {
  definition <- subset_definitions[[object_name]]
  log_message( "Creating and saving object: ", object_name, "." )
  subset_object <- make_subset_object( object = obj, cells = definition$cells, object_name = object_name, expected_nuclei = as.integer( definition$expected_n ) )
  subset_object@misc$GSE234774_step10G_object_role <- list( object_name = object_name, object_role = definition$object_role, frozen_cluster_column = FROZEN_CLUSTER_COLUMN, source_full_object = object_paths[[ "FULL_ANNOTATED_MYELOID" ]], completed_at = as.character( Sys.time() ) )
  saveRDS( subset_object, object_paths[[ object_name ]], compress = FALSE )
  save_summary_rows[[object_name]] <- data.frame( object_name = object_name, expected_n_cells = as.integer( definition$expected_n ), observed_n_cells = ncol( subset_object ), expected_n_features = EXPECTED_FEATURES, observed_n_features = nrow( subset_object ), valid_object = safe_valid_object( subset_object ), object_role = definition$object_role, output_path = object_paths[[ object_name ]], stringsAsFactors = FALSE )
  rm( subset_object )
  invisible( gc() )
}
save_summary <- do.call( rbind, save_summary_rows )
rownames( save_summary ) <- NULL
log_message( "Reloading and auditing all five Step 10G objects." )
reload_audit_rows <- list()
for (object_name in names( object_cell_sets )) {
  expected_cells <- object_cell_sets[[object_name]]
  reloaded_object <- tryCatch( readRDS( object_paths[[ object_name ]] ), error = function(error_object) {
      error_object
    }
  )
  reload_success <- !inherits( reloaded_object, "error" )
  if (reload_success) {
    dimension_pass <- nrow( reloaded_object ) == EXPECTED_FEATURES && ncol( reloaded_object ) == length( expected_cells )
    barcode_pass <- identical( colnames( reloaded_object ), expected_cells )
    metadata_pass <- all( c( "myeloid_final_label", "resident_lineage_status", "arm_analysis_set", "nonIFN_sensitivity_status", "step10G_annotation_status" ) %in% names( reloaded_object@meta.data ) )
    reduction_pass <- all( c( FROZEN_REDUCTION, FROZEN_UMAP ) %in% names( reloaded_object@reductions ) )
    validity_pass <- safe_valid_object( reloaded_object )
  } else {
    dimension_pass <- FALSE
    barcode_pass <- FALSE
    metadata_pass <- FALSE
    reduction_pass <- FALSE
    validity_pass <- FALSE
  }
  reload_audit_rows[[object_name]] <- data.frame( object_name = object_name, reload_success = reload_success, dimension_pass = dimension_pass, barcode_pass = barcode_pass, metadata_pass = metadata_pass, reduction_pass = reduction_pass, validity_pass = validity_pass, stringsAsFactors = FALSE )
  rm( reloaded_object )
  invisible( gc() )
}
reload_audit <- do.call( rbind, reload_audit_rows )
rownames( reload_audit ) <- NULL
all_objects_reload_pass <- all( reload_audit$reload_success & reload_audit$dimension_pass & reload_audit$barcode_pass & reload_audit$metadata_pass & reload_audit$reduction_pass & reload_audit$validity_pass )
annotated_counts_sample <- sample_layer_values( object = obj, layer_name = counts_layer_name, features = integrity_features, cells = integrity_cells )
annotated_data_sample <- sample_layer_values( object = obj, layer_name = data_layer_name, features = integrity_features, cells = integrity_cells )
counts_unchanged <- compare_matrix_values( matrix_a = input_counts_sample, matrix_b = annotated_counts_sample, tolerance = 0 )
data_unchanged <- compare_matrix_values( matrix_a = input_data_sample, matrix_b = annotated_data_sample, tolerance = 1e-12 )
if ( !counts_unchanged || !data_unchanged ) {
  stop_with_report( paste0( "RNA integrity verification failed. counts_unchanged=", counts_unchanged, "; data_unchanged=", data_unchanged, "." ) )
}
primary_metadata <- obj@meta.data[ primary_cells, , drop = FALSE ]
audit_checks <- data.frame( check = c( "Input contains 26,918 features and 25,205 myeloid nuclei", "Frozen myhcon_res_0_6 cluster counts match the reviewed Step 10F counts", "Every nucleus received a frozen final label and analysis-set assignment", "Primary, cycling, and excluded groups form an exact one-way partition", "The non-IFN sensitivity group is a strict subset of the primary resident group", "All five object counts match the frozen expected counts", "Primary resident library/timepoint coverage was quantified and reported", "All five RDS objects reload with exact barcode order, metadata, reductions, and valid structure", "RNA counts and LogNormalize data remain unchanged", "Cluster 7 is included in the primary resident object", "Cluster 6 is excluded from the primary object and retained in the cycling object", "Clusters 4, 5, and 8 are excluded from both resident analysis objects", "All three required figures and two core tables were created", "Arm_1/Arm_2 scoring was run" ), required_for_pass = c( TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE ), status = c( "PASS", ifelse( identical( observed_cluster_counts, EXPECTED_CLUSTER_COUNTS ), "PASS", "STOP" ), ifelse( !anyNA( obj$myeloid_final_label ) && !anyNA( obj$arm_analysis_set ), "PASS", "STOP" ), ifelse( exclusive_partition_pass, "PASS", "STOP" ), ifelse( nonIFN_subset_pass && all( nonIFN_cells %in% primary_cells ) && !any( primary_cells[ cluster_values[ match( primary_cells, colnames( obj ) ) ] == "7" ] %in% nonIFN_cells ), "PASS", "STOP" ), ifelse( identical( as.integer( observed_object_counts ), as.integer( EXPECTED_OBJECT_COUNTS[ names( observed_object_counts ) ] ) ), "PASS", "STOP" ), ifelse( primary_coverage_count_pass && nrow( primary_library_coverage ) == nrow( source_library_manifest ), "PASS", "STOP" ), ifelse( all_objects_reload_pass, "PASS", "STOP" ), ifelse( counts_unchanged && data_unchanged, "PASS", "STOP" ), ifelse( all( colnames( obj )[ cluster_values == "7" ] %in% primary_cells ), "PASS", "STOP" ), ifelse( !any( colnames( obj )[ cluster_values == "6" ] %in% primary_cells ) && setequal( colnames( obj )[ cluster_values == "6" ], cycling_cells ), "PASS", "STOP" ), ifelse( !any( colnames( obj )[ cluster_values %in% EXCLUDED_CLUSTERS ] %in% primary_cells ) && !any( colnames( obj )[ cluster_values %in% EXCLUDED_CLUSTERS ] %in% nonIFN_cells ) && setequal( colnames( obj )[ cluster_values %in% EXCLUDED_CLUSTERS ], excluded_cells ), "PASS", "STOP" ), ifelse( all( file.exists( file.path( figures_dir, c( "01_frozen_annotation_and_analysis_set_UMAP.pdf", "02_analysis_set_timepoint_composition.pdf", "03_final_annotation_marker_DotPlot.pdf" ) ) ) ) && all( file.exists( file.path( tables_dir, c( "04_frozen_cluster_annotation_and_object_membership.csv", "05_frozen_object_composition_summary.csv" ) ) ) ), "PASS", "STOP" ), "CONFIRMED_NOT_RUN" ), detail = c( paste0( nrow( obj ), " features x ", ncol( obj ), " nuclei" ), paste( paste0( EXPECTED_CLUSTERS, ":", observed_cluster_counts[ EXPECTED_CLUSTERS ] ), collapse = "; " ), paste0( length( unique( as.character( obj$myeloid_final_label ) ) ), " final labels; ", length( unique( as.character( obj$arm_analysis_set ) ) ), " mutually exclusive analysis sets" ), paste0( "exclusive membership count range=", min( exclusive_membership_count ), "-", max( exclusive_membership_count ) ), paste0( "non-IFN n=", length( nonIFN_cells ), "; primary n=", length( primary_cells ), "; cluster 7 retained only in primary" ), paste( paste0( names( observed_object_counts ), "=", observed_object_counts ), collapse = "; " ), paste0( "retained libraries=", length( primary_present_libraries ), "/", length( all_source_libraries ), "; Step10F source timepoints=", length( all_source_timepoints ), "/", length( expected_timepoints ), "; primary retained source timepoints=", length( primary_present_timepoints ), "/", length( all_source_timepoints ), "; primary retained project timepoints=", length( primary_present_timepoints ), "/", length( expected_timepoints ), "; missing libraries=", ifelse( length( missing_primary_libraries ) == 0L, "none", paste( missing_primary_libraries, collapse = "," ) ), "; source missing expected timepoints=", ifelse( length( missing_source_timepoints ) == 0L, "none", paste( missing_source_timepoints, collapse = "," ) ), "; primary missing expected timepoints=", ifelse( length( missing_primary_expected_timepoints ) == 0L, "none", paste( missing_primary_expected_timepoints, collapse = "," ) ), "; source timepoints lost after primary filtering=", ifelse( length( missing_primary_from_source_timepoints ) == 0L, "none", paste( missing_primary_from_source_timepoints, collapse = "," ) ), "; timepoints with <2 retained libraries=", ifelse( length( low_replication_timepoints ) == 0L, "none", paste( low_replication_timepoints, collapse = "," ) ) ), paste0( "objects passing reload audit=", sum( reload_audit$reload_success & reload_audit$dimension_pass & reload_audit$barcode_pass & reload_audit$metadata_pass & reload_audit$reduction_pass & reload_audit$validity_pass ), "/", nrow( reload_audit ) ), paste0( "30 features x 120 nuclei sampled; counts unchanged=", counts_unchanged, "; LogNormalize data unchanged=", data_unchanged ), paste0( "cluster 7 nuclei=", EXPECTED_CLUSTER_COUNTS[[ "7" ]], "; primary includes IFN-responsive microglia" ), paste0( "cluster 6 nuclei=", length( cycling_cells ), "; separated from primary Arm analysis" ), paste0( "excluded clusters=", paste( EXCLUDED_CLUSTERS, collapse = "," ), "; excluded nuclei=", length( excluded_cells ) ), "Three PDFs and two core CSV tables present.", "Arm_1/Arm_2 signatures remain frozen from GSE162610 and were not scored." ), stringsAsFactors = FALSE )
write_csv_safe( audit_checks, "09_STEP10G_audit_checks.csv" )
required_failures <- audit_checks$required_for_pass & audit_checks$status != "PASS"
coverage_review_needed <- !primary_all_libraries_retained || !primary_all_source_timepoints_retained || !source_has_all_expected_timepoints || length( low_replication_timepoints ) > 0L
final_decision <- if ( any( required_failures ) ) {
  "STOP"
} else if (
  coverage_review_needed
) {
  "PASS_FROZEN_RESIDENT_OBJECTS_WITH_COVERAGE_REVIEW"
} else {
  "PASS_FROZEN_RESIDENT_OBJECTS"
}
object_lines <- paste0( save_summary$object_name, ": ", save_summary$observed_n_cells, " nuclei; ", save_summary$observed_n_features, " features; valid=", save_summary$valid_object, "; role=", save_summary$object_role, "; path=", save_summary$output_path )
cluster_lines <- paste0( "cluster ", cluster_annotation_table$cluster, ": ", cluster_annotation_table$myeloid_final_label, "; n=", cluster_annotation_table$observed_n_cells, "; resident_status=", cluster_annotation_table$resident_lineage_status, "; arm_set=", cluster_annotation_table$arm_analysis_set, "; nonIFN_status=", cluster_annotation_table$nonIFN_sensitivity_status )
report_lines <- c( "GSE234774 STEP 10G FROZEN MYELOID ANNOTATIONS / RESIDENT OBJECTS", "============================================================", paste0( "Run time: ", format( Sys.time(), "%Y-%m-%d %H:%M:%S" ) ), paste0( "Input object: ", input_object_path ), paste0( "Output directory: ", output_dir ), "", paste0( "FINAL DECISION: ", final_decision ), "", "FROZEN CLUSTER ANNOTATIONS", cluster_lines, "", "OBJECTS", object_lines, "", "PRIMARY ARM ANALYSIS DESIGN", "Primary non-cycling resident clusters: 0, 1, 2, 3, 7.", paste0( "Primary non-cycling resident nuclei: ", length( primary_cells ), "." ), "Cluster 7 IFN-responsive microglia is retained in the primary analysis.", paste0( "Step 10F myeloid source libraries: ", length( all_source_libraries ), "." ), paste0( "Step 10F myeloid source timepoints: ", length( all_source_timepoints ), "/", length( expected_timepoints ), "; absent expected timepoints=", ifelse( length( missing_source_timepoints ) == 0L, "none", paste( missing_source_timepoints, collapse = ", " ) ), "." ), paste0( "Primary retained source libraries: ", length( primary_present_libraries ), "/", length( all_source_libraries ), "." ), paste0( "Source libraries with zero retained primary resident nuclei: ", ifelse( length( missing_primary_libraries ) == 0L, "none", paste( missing_primary_libraries, collapse = ", " ) ), "." ), paste0( "Primary retained Step 10F source timepoints: ", length( primary_present_timepoints ), "/", length( all_source_timepoints ), "." ), paste0( "Primary retained expected project timepoints: ", length( primary_present_timepoints ), "/", length( expected_timepoints ), "." ), paste0( "Absent primary timepoints: ", ifelse( length( missing_primary_expected_timepoints ) == 0L, "none", paste( missing_primary_expected_timepoints, collapse = ", " ) ), "." ), paste0( "Step 10F source timepoints lost after primary filtering: ", ifelse( length( missing_primary_from_source_timepoints ) == 0L, "none", paste( missing_primary_from_source_timepoints, collapse = ", " ) ), "." ), paste0( "Timepoints with fewer than two retained biological libraries: ", ifelse( length( low_replication_timepoints ) == 0L, "none", paste( low_replication_timepoints, collapse = ", " ) ), "." ), "A zero-cell source library is a retained-coverage limitation, not an object-construction error.", "Non-IFN sensitivity clusters: 0, 1, 2, 3.", paste0( "Non-IFN sensitivity nuclei: ", length( nonIFN_cells ), "." ), "Cycling cluster 6 is retained as resident lineage but analyzed separately.", paste0( "Cycling resident nuclei: ", length( cycling_cells ), "." ), "Clusters 4, 5, and 8 are excluded from resident Arm analysis.", paste0( "Excluded nuclei: ", length( excluded_cells ), "." ), "", "INTEGRITY", paste0( "Exclusive primary/cycling/excluded partition: ", exclusive_partition_pass ), paste0( "Non-IFN is a strict primary subset: ", nonIFN_subset_pass ), paste0( "All five objects reload and validate: ", all_objects_reload_pass ), paste0( "RNA counts unchanged: ", counts_unchanged ), paste0( "RNA LogNormalize data unchanged: ", data_unchanged ), "", "INTERPRETATION BOUNDARIES", "1. Cluster-level myeloid annotation and object membership are frozen in this step.", "2. Cycling microglia remains a state, not a terminal lineage.", "3. IFN-responsive microglia remains in the primary resident analysis and is removed only in a dedicated sensitivity object.", "4. BAM-like, macrophage-like non-microglial, and neuronal-contaminant clusters are excluded from resident Arm analysis.", "5. Harmony is retained only for visualization and cluster context.", "6. RNA counts and LogNormalize data were not modified.", "7. Arm_1/Arm_2 scoring was not run.", "", "STATUS", "Project the frozen GSE162610 Arm_1 and Arm_2 signatures onto the primary non-cycling resident object.", "Repeat the projection in the non-IFN sensitivity object.", "Use only retained biological libraries for sample-level summaries and explicitly report any zero-cell libraries or low-replication timepoints.", "Report cycling resident microglia separately without forcing cluster 6 into either Arm." )
if ( final_decision == "STOP" ) {
  failed <- audit_checks[ audit_checks$required_for_pass & audit_checks$status != "PASS", , drop = FALSE ]
  report_lines <- c( report_lines, "", "STOP REASONS" )
  for ( index in seq_len( nrow( failed ) ) ) {
    report_lines <- c( report_lines, paste0( index, ". ", failed$check[ index ], " -- ", failed$detail[ index ] ) )
  }
}
report_path <- file.path( output_dir, "10_STEP10G_FINAL_REPORT.txt" )
writeLines( report_lines, report_path, useBytes = TRUE )
log_message( "Final decision: ", final_decision )
log_message( "Final report: ", report_path )
cat("\n")
cat("============================================================\n")
cat("GSE234774 STEP 10G COMPLETED\n")
cat( "FINAL DECISION:", final_decision, "\n" )
cat( "Output:", output_dir, "\n" )
cat("Primary Arm-analysis object:\n")
cat( " ", object_paths[[ "PRIMARY_NONCYCLING_RESIDENT" ]], "\n", sep = "" )
cat("Non-IFN sensitivity object:\n")
cat( " ", object_paths[[ "NONIFN_SENSITIVITY_RESIDENT" ]], "\n", sep = "" )
cat("Key outputs:\n")
cat(" 1) figures/01_frozen_annotation_and_analysis_set_UMAP.pdf\n")
cat(" 2) figures/02_analysis_set_timepoint_composition.pdf\n")
cat(" 3) figures/03_final_annotation_marker_DotPlot.pdf\n")
cat(" 4) tables/04_frozen_cluster_annotation_and_object_membership.csv\n")
cat(" 5) tables/05_frozen_object_composition_summary.csv\n")
cat(" 6) tables/09_STEP10G_audit_checks.csv\n")
cat(" 7) 10_STEP10G_FINAL_REPORT.txt\n")
cat("============================================================\n")
if ( final_decision == "STOP" ) {
  stop( "GSE234774 Step 10G decision is STOP. Review the final report.", call. = FALSE )
}
