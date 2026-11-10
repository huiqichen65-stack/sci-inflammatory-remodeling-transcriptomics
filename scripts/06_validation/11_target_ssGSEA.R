options( stringsAsFactors = FALSE, warn = 1, Seurat.object.assay.version = "v5", future.globals.maxSize = 24 * 1024^3 )
set.seed( 20260803L )
PROJECT_ROOT <- Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
MAIN_TIMEPOINTS <- c( "Uninjured", "7 dpi", "14 dpi", "1 month", "2 months" )
TIMEPOINT_LEVELS <- MAIN_TIMEPOINTS
FROZEN_CLUSTER_COLUMN <- "myhcon_res_0_6"
FROZEN_UMAP <- "umap_myeloid_harmony"
PRIMARY_CLUSTERS <- c( "0", "1", "2", "3", "7" )
EXPECTED_PRIMARY_INPUT_FEATURES <- 26918L
EXPECTED_PRIMARY_INPUT_NUCLEI <- 20179L
MIN_CELLS_PER_LIBRARY <- 50L
MIN_LIBRARIES_PER_MAIN_TIMEPOINT <- 2L
MIN_SIGNATURE_GENES_PRESENT <- 15L
MIN_SIGNATURE_MAPPING_FRACTION <- 0.70
AUC_MAX_RANK_FRACTION <- 0.05
SSGSEA_ALPHA <- 0.25
PSEUDOBULK_PRIOR_COUNT <- 2
SHARED_STEP6B_MODULE <- "M3_Shared_prebranch"
ARM2_STEP6B_MODULES <- c( "M1_Arm2_terminal", "M2_Arm2_terminal" )
TECHNICAL_GENE_PATTERN <- "^(mt-|Rpl|Rps|Hb[ab])"
SCORED_PROGRAM_NAMES <- c( "Arm_1_specific", "Arm_2_specific", "Shared_prebranch" )
PRIMARY_OBJECT_PATH <- file.path( PROJECT_ROOT, "04_results", "step10G_GSE234774_frozen_myeloid_annotations_resident_objects", "RDS", "GSE234774_step10G_primary_noncycling_resident_microglia_BPCells.rds" )
STEP6A_DIR <- file.path( PROJECT_ROOT, "04_results", "step6A_GSE162610_7dpi_branch_specific_DEG" )
STEP6A_ARM1_FILE <- file.path( STEP6A_DIR, "tables", "03_Arm1_biased_high_confidence_genes.csv" )
STEP6A_ARM2_FILE <- file.path( STEP6A_DIR, "tables", "04_Arm2_biased_high_confidence_genes.csv" )
STEP6B_ASSIGNMENT_FILE <- file.path( PROJECT_ROOT, "04_results", "step6B_GSE162610_dynamic_branch_modules", "tables", "01_dynamic_module_gene_assignments.csv" )
OUTPUT_DIR <- file.path( PROJECT_ROOT, "04_results", "step11B_GSE234774_refined_main_Arm_specific_shared_AUCell_ssGSEA" )
FIGURES_DIR <- file.path( OUTPUT_DIR, "figures" )
TABLES_DIR <- file.path( OUTPUT_DIR, "tables" )
RDS_DIR <- file.path( OUTPUT_DIR, "RDS" )
for (directory in c( OUTPUT_DIR, FIGURES_DIR, TABLES_DIR, RDS_DIR )) {
  dir.create( directory, recursive = TRUE, showWarnings = FALSE )
}
RUN_LOG <- file.path( OUTPUT_DIR, "00_run_log.txt" )
if (file.exists( RUN_LOG )) {
  file.remove( RUN_LOG )
}
SCORED_OBJECT_PATH <- file.path( RDS_DIR, "GSE234774_step11B_refined_main_Arm_scored_BPCells.rds" )
RESULTS_RDS_PATH <- file.path( RDS_DIR, "GSE234774_step11B_refined_Arm_AUCell_ssGSEA_results.rds" )
log_message <- function(...) {
  message_text <- paste0(...)
  timestamp <- format( Sys.time(), "%Y-%m-%d %H:%M:%S" )
  line <- paste0( "[", timestamp, "] ", message_text )
  cat( line, "\n" )
  cat( line, "\n", file = RUN_LOG, append = TRUE )
}
write_csv_safe <- function( x, filename ) {
  path <- file.path( TABLES_DIR, filename )
  utils::write.csv( x, path, row.names = FALSE, na = "" )
  log_message( "Wrote: ", path )
  invisible( path )
}
stop_with_report <- function(reason) {
  report_path <- file.path( OUTPUT_DIR, "10_STEP11B_FINAL_REPORT.txt" )
  report_lines <- c( "GSE234774 STEP 11B REFINED MAIN ARM VALIDATION", "============================================================", paste0( "Run time: ", format( Sys.time(), "%Y-%m-%d %H:%M:%S" ) ), paste0( "Output directory: ", OUTPUT_DIR ), "", "FINAL DECISION: STOP", "", paste0( "Reason: ", reason ), "", "Do not interpret refined Arm validation from this run." )
  writeLines( report_lines, report_path, useBytes = TRUE )
  stop( reason, call. = FALSE )
}
safe_package_version <- function( package_name ) {
  if (!requireNamespace( package_name, quietly = TRUE )) {
    return( NA_character_ )
  }
  as.character( utils::packageVersion( package_name ) )
}
require_columns <- function( data, columns, label ) {
  missing_columns <- setdiff( columns, colnames( data ) )
  if ( length( missing_columns ) > 0L ) {
    stop_with_report( paste0( label, " is missing required column(s): ", paste( missing_columns, collapse = ", " ), "." ) )
  }
}
detect_gene_column <- function( data, label ) {
  candidates <- c( "gene_id", "gene", "Gene", "symbol", "SYMBOL" )
  detected <- intersect( candidates, colnames( data ) )
  if ( length( detected ) == 0L ) {
    stop_with_report( paste0( "Could not identify a gene-symbol column in ", label, "." ) )
  }
  detected[[1L]]
}
clean_gene_vector <- function(x) {
  output <- trimws( as.character( x ) )
  output <- output[ !is.na( output ) & nzchar( output ) ]
  unique( output )
}
resolve_features <- function( requested_features, available_features ) {
  available_upper <- toupper( available_features )
  first_index <- !duplicated( available_upper )
  lookup <- setNames( available_features[ first_index ], available_upper[ first_index ] )
  resolved <- unname( lookup[ toupper( requested_features ) ] )
  names( resolved ) <- requested_features
  resolved
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
to_dgCMatrix <- function( matrix_object, label ) {
  if (inherits( matrix_object, "dgCMatrix" )) {
    return( matrix_object )
  }
  converted <- tryCatch( methods::as( matrix_object, "dgCMatrix" ), error = function(error_object) {
      error_object
    }
  )
  if (inherits( converted, "error" )) {
    intermediate <- tryCatch( methods::as( matrix_object, "CsparseMatrix" ), error = function(error_object) {
        error_object
      }
    )
    if (!inherits( intermediate, "error" )) {
      converted <- tryCatch( methods::as( intermediate, "dgCMatrix" ), error = function(error_object) {
          error_object
        }
      )
    }
  }
  if (inherits( converted, "error" )) {
    stop_with_report( paste0( "Could not convert ", label, " to dgCMatrix without dense whole-cell conversion: ", conditionMessage( converted ) ) )
  }
  converted
}
run_ssgsea_version_safe <- function( expression_matrix, gene_sets ) {
  gsva_exports <- getNamespaceExports( "GSVA" )
  if ( "ssgseaParam" %in% gsva_exports ) {
    log_message( "Running ssGSEA with the current GSVA ssgseaParam() interface." )
    parameter_arguments <- list( exprData = expression_matrix, geneSets = gene_sets, minSize = MIN_SIGNATURE_GENES_PRESENT, maxSize = Inf, alpha = SSGSEA_ALPHA, normalize = TRUE, checkNA = "yes", use = "all.obs", ondisk = "no", verbose = TRUE )
    parameter_formals <- names( formals( GSVA::ssgseaParam ) )
    if ( !"..." %in% parameter_formals ) {
      parameter_arguments <- parameter_arguments[ names( parameter_arguments ) %in% parameter_formals ]
    }
    parameter_object <- do.call( GSVA::ssgseaParam, parameter_arguments )
    result <- GSVA::gsva( parameter_object, verbose = TRUE, BPPARAM = BiocParallel::SerialParam() )
  } else {
    log_message( "Running ssGSEA with the legacy GSVA::gsva(method='ssgsea') interface." )
    legacy_arguments <- list( expr = expression_matrix, gset.idx.list = gene_sets, method = "ssgsea", kcdf = "Gaussian", abs.ranking = TRUE, ssgsea.norm = TRUE, min.sz = MIN_SIGNATURE_GENES_PRESENT, max.sz = Inf, parallel.sz = 1L, verbose = TRUE )
    legacy_formals <- names( formals( GSVA::gsva ) )
    if ( !"..." %in% legacy_formals ) {
      legacy_arguments <- legacy_arguments[ names( legacy_arguments ) %in% legacy_formals ]
    }
    result <- do.call( GSVA::gsva, legacy_arguments )
  }
  as.matrix( result )
}
plot_umap_score <- function( plot_data, score_column, title_text, legend_title ) {
  current_data <- plot_data
  current_data$score <- current_data[[score_column]]
  ggplot2::ggplot( current_data, ggplot2::aes( x = UMAP_1, y = UMAP_2 ) ) + ggplot2::geom_point( ggplot2::aes( colour = score ), size = 0.28, alpha = 0.72 ) + ggplot2::scale_colour_viridis_c() + ggplot2::coord_equal() + ggplot2::labs( title = title_text, x = "Myeloid Harmony UMAP 1", y = "Myeloid Harmony UMAP 2", colour = legend_title ) + ggplot2::theme_bw()
}
safe_cor_test <- function( x, y ) {
  complete <- stats::complete.cases( x, y )
  x_complete <- x[ complete ]
  y_complete <- y[ complete ]
  if ( length( x_complete ) < 4L || stats::sd( x_complete ) == 0 || stats::sd( y_complete ) == 0 ) {
    return( list( rho = NA_real_, p_value = NA_real_, n = length( x_complete ) ) )
  }
  test <- suppressWarnings( stats::cor.test( x_complete, y_complete, method = "spearman", exact = FALSE ) )
  list( rho = unname( test$estimate ), p_value = test$p.value, n = length( x_complete ) )
}
log_message( "Starting package preflight." )
required_packages <- c( "Matrix", "data.table", "ggplot2", "patchwork", "scales", "Seurat", "SeuratObject", "AUCell", "GSVA", "edgeR", "BiocParallel" )
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
  library(AUCell)
  library(GSVA)
  library(edgeR)
  library(BiocParallel)
})
package_versions <- data.frame( package = c( "R", required_packages ), version = c( paste0( R.version$major, ".", R.version$minor ), vapply( required_packages, safe_package_version, FUN.VALUE = character(1) ) ), stringsAsFactors = FALSE )
write_csv_safe( package_versions, "00_package_versions.csv" )
input_files <- c( primary_resident_object = PRIMARY_OBJECT_PATH, Step6A_Arm1_highconf = STEP6A_ARM1_FILE, Step6A_Arm2_highconf = STEP6A_ARM2_FILE, Step6B_module_assignments = STEP6B_ASSIGNMENT_FILE )
input_file_audit <- data.frame( input_name = names( input_files ), path = unname( input_files ), exists = file.exists( input_files ), size_bytes = vapply( input_files, function(path) {
      if ( file.exists( path ) ) {
        file.info( path )$size
      } else {
        NA_real_
      }
    },
    numeric(1)
  ), stringsAsFactors = FALSE
)
write_csv_safe( input_file_audit, "01_input_file_audit.csv" )
if ( !all( input_file_audit$exists ) ) {
  stop_with_report( paste0( "Missing required input file(s): ", paste( input_file_audit$input_name[ !input_file_audit$exists ], collapse = ", " ), "." ) )
}
log_message( "Loading the frozen Step 10G primary non-cycling resident object." )
primary_object <- readRDS( PRIMARY_OBJECT_PATH )
if (!inherits( primary_object, "Seurat" )) {
  stop_with_report( "The Step 10G primary input is not a Seurat object." )
}
if ( nrow( primary_object ) != EXPECTED_PRIMARY_INPUT_FEATURES || ncol( primary_object ) != EXPECTED_PRIMARY_INPUT_NUCLEI ) {
  stop_with_report( paste0( "Unexpected Step 10G primary-object dimensions. Expected ", EXPECTED_PRIMARY_INPUT_FEATURES, " x ", EXPECTED_PRIMARY_INPUT_NUCLEI, "; detected ", nrow( primary_object ), " x ", ncol( primary_object ), "." ) )
}
required_metadata <- c( "biological_library", "timepoint", "myeloid_final_label", "arm_analysis_set", FROZEN_CLUSTER_COLUMN )
missing_metadata <- setdiff( required_metadata, colnames( primary_object@meta.data ) )
if ( length( missing_metadata ) > 0L ) {
  stop_with_report( paste0( "Primary object is missing metadata column(s): ", paste( missing_metadata, collapse = ", " ), "." ) )
}
if (!FROZEN_UMAP %in% names( primary_object@reductions ) ) {
  stop_with_report( paste0( "Frozen UMAP reduction is missing: ", FROZEN_UMAP, "." ) )
}
input_cluster_values <- sort( unique( as.character( primary_object@meta.data[ , FROZEN_CLUSTER_COLUMN ] ) ), method = "radix" )
if (!setequal( input_cluster_values, PRIMARY_CLUSTERS )) {
  stop_with_report( paste0( "Primary object clusters do not match 0,1,2,3,7. Detected: ", paste( input_cluster_values, collapse = ", " ), "." ) )
}
if ( any( as.character( primary_object$arm_analysis_set ) != "PRIMARY_NONCYCLING_RESIDENT" ) ) {
  stop_with_report( "The Step 10G primary object contains cells outside PRIMARY_NONCYCLING_RESIDENT." )
}
primary_timepoint_values <- as.character( primary_object$timepoint )
main_timepoint_cells <- colnames( primary_object )[ primary_timepoint_values %in% MAIN_TIMEPOINTS ]
if ( length( main_timepoint_cells ) == 0L ) {
  stop_with_report( "No cells were retained for the five main validation timepoints." )
}
main_preQC_object <- subset( x = primary_object, cells = main_timepoint_cells )
main_preQC_object$timepoint <- factor( as.character( main_preQC_object$timepoint ), levels = TIMEPOINT_LEVELS, ordered = TRUE )
preQC_library_counts <- data.table::data.table( biological_library = as.character( main_preQC_object$biological_library ), timepoint = as.character( main_preQC_object$timepoint ) )[ , .( n_cells_before_library_QC = .N ), by = .( biological_library, timepoint ) ]
preQC_library_counts$timepoint <- factor( preQC_library_counts$timepoint, levels = TIMEPOINT_LEVELS, ordered = TRUE )
data.table::setorder( preQC_library_counts, timepoint, biological_library )
preQC_library_counts$minimum_cell_threshold <- MIN_CELLS_PER_LIBRARY
preQC_library_counts$library_QC_status <- ifelse( preQC_library_counts$n_cells_before_library_QC >= MIN_CELLS_PER_LIBRARY, "INCLUDED", "EXCLUDED_LOW_CELL_COUNT" )
included_libraries <- preQC_library_counts$biological_library[ preQC_library_counts$library_QC_status == "INCLUDED" ]
excluded_libraries <- preQC_library_counts$biological_library[ preQC_library_counts$library_QC_status == "EXCLUDED_LOW_CELL_COUNT" ]
included_libraries <- sort( unique( included_libraries ) )
excluded_libraries <- sort( unique( excluded_libraries ) )
if ( length( included_libraries ) == 0L ) {
  stop_with_report( "All main-timepoint libraries failed the minimum-cell threshold." )
}
main_cells <- colnames( main_preQC_object )[ as.character( main_preQC_object$biological_library ) %in% included_libraries ]
main_object <- subset( x = main_preQC_object, cells = main_cells )
main_object$timepoint <- factor( as.character( main_object$timepoint ), levels = TIMEPOINT_LEVELS, ordered = TRUE )
if ( any( as.character( main_object$timepoint ) %in% c( "1 dpi", "4 dpi" ) ) ) {
  stop_with_report( "The refined main analysis still contains 1 dpi or 4 dpi cells." )
}
if ( "Cycling_microglia" %in% as.character( main_object$myeloid_final_label ) ) {
  stop_with_report( "Cycling microglia is present in the refined main object." )
}
if ( !"IFN_responsive_microglia" %in% as.character( main_object$myeloid_final_label ) ) {
  stop_with_report( "IFN-responsive microglia cluster 7 is absent from the refined main object." )
}
postQC_library_counts <- data.table::data.table( biological_library = as.character( main_object$biological_library ), timepoint = as.character( main_object$timepoint ) )[ , .( n_cells_after_library_QC = .N ), by = .( biological_library, timepoint ) ]
library_QC_audit <- merge( preQC_library_counts, postQC_library_counts, by = c( "biological_library", "timepoint" ), all.x = TRUE, sort = FALSE )
library_QC_audit$n_cells_after_library_QC[ is.na( library_QC_audit$n_cells_after_library_QC ) ] <- 0L
library_QC_audit$timepoint <- factor( library_QC_audit$timepoint, levels = TIMEPOINT_LEVELS, ordered = TRUE )
data.table::setorder( library_QC_audit, timepoint, biological_library )
postQC_timepoint_replication <- data.table::data.table( biological_library = as.character( main_object$biological_library ), timepoint = as.character( main_object$timepoint ) )[ , .( n_cells = .N ), by = .( biological_library, timepoint ) ][ , .( n_libraries = uniqueN( biological_library ), n_cells = sum( n_cells ), minimum_library_cells = min( n_cells ), maximum_library_cells = max( n_cells ) ), by = timepoint ]
postQC_timepoint_replication$timepoint <- factor( postQC_timepoint_replication$timepoint, levels = TIMEPOINT_LEVELS, ordered = TRUE )
data.table::setorder( postQC_timepoint_replication, timepoint )
missing_main_timepoints <- setdiff( MAIN_TIMEPOINTS, as.character( postQC_timepoint_replication$timepoint ) )
if ( length( missing_main_timepoints ) > 0L ) {
  stop_with_report( paste0( "A main validation timepoint was lost after library QC: ", paste( missing_main_timepoints, collapse = ", " ), "." ) )
}
low_replication_timepoints <- as.character( postQC_timepoint_replication$timepoint[ postQC_timepoint_replication$n_libraries < MIN_LIBRARIES_PER_MAIN_TIMEPOINT ] )
if ( length( low_replication_timepoints ) > 0L ) {
  stop_with_report( paste0( "A refined main timepoint has fewer than ", MIN_LIBRARIES_PER_MAIN_TIMEPOINT, " retained libraries: ", paste( low_replication_timepoints, collapse = ", " ), "." ) )
}
library_QC_output <- rbind( data.frame( record_type = "LIBRARY", biological_library = library_QC_audit$biological_library, timepoint = as.character( library_QC_audit$timepoint ), n_cells_before_library_QC = library_QC_audit$n_cells_before_library_QC, minimum_cell_threshold = library_QC_audit$minimum_cell_threshold, library_QC_status = library_QC_audit$library_QC_status, n_cells_after_library_QC = library_QC_audit$n_cells_after_library_QC, n_libraries_after_QC = NA_integer_, minimum_retained_library_cells = NA_integer_, maximum_retained_library_cells = NA_integer_, stringsAsFactors = FALSE ), data.frame( record_type = "TIMEPOINT_SUMMARY", biological_library = "ALL", timepoint = as.character( postQC_timepoint_replication$timepoint ), n_cells_before_library_QC = NA_integer_, minimum_cell_threshold = MIN_CELLS_PER_LIBRARY, library_QC_status = "TIMEPOINT_RETAINED", n_cells_after_library_QC = postQC_timepoint_replication$n_cells, n_libraries_after_QC = postQC_timepoint_replication$n_libraries, minimum_retained_library_cells = postQC_timepoint_replication$minimum_library_cells, maximum_retained_library_cells = postQC_timepoint_replication$maximum_library_cells, stringsAsFactors = FALSE ) )
write_csv_safe( library_QC_output, "02_main_library_minimum_cell_QC.csv" )
log_message( "Library QC retained ", length( included_libraries ), " libraries and excluded ", length( excluded_libraries ), " low-cell library/libraries." )
if ( length( excluded_libraries ) > 0L ) {
  log_message( "Excluded low-cell libraries: ", paste( excluded_libraries, collapse = ", " ), "." )
}
log_message( "Loading Step6A and Step6B source-defined gene sets." )
step6a_arm1 <- data.table::fread( STEP6A_ARM1_FILE, data.table = FALSE, check.names = FALSE, showProgress = FALSE )
step6a_arm2 <- data.table::fread( STEP6A_ARM2_FILE, data.table = FALSE, check.names = FALSE, showProgress = FALSE )
step6b_assignments <- data.table::fread( STEP6B_ASSIGNMENT_FILE, data.table = FALSE, check.names = FALSE, showProgress = FALSE )
arm1_gene_column <- detect_gene_column( step6a_arm1, "Step6A Arm_1 high-confidence marker table" )
arm2_gene_column <- detect_gene_column( step6a_arm2, "Step6A Arm_2 high-confidence marker table" )
require_columns( step6b_assignments, c( "gene_id", "module_name" ), "Step6B dynamic module assignment table" )
step6a_arm1_genes <- clean_gene_vector( step6a_arm1[ , arm1_gene_column ] )
step6a_arm2_genes <- clean_gene_vector( step6a_arm2[ , arm2_gene_column ] )
step6b_assignments$gene_id <- trimws( as.character( step6b_assignments$gene_id ) )
step6b_assignments$module_name <- trimws( as.character( step6b_assignments$module_name ) )
step6b_assignments <- step6b_assignments[ !is.na( step6b_assignments$gene_id ) & nzchar( step6b_assignments$gene_id ) & !is.na( step6b_assignments$module_name ) & nzchar( step6b_assignments$module_name ), , drop = FALSE ]
required_modules <- c( SHARED_STEP6B_MODULE, ARM2_STEP6B_MODULES )
missing_modules <- setdiff( required_modules, unique( step6b_assignments$module_name ) )
if ( length( missing_modules ) > 0L ) {
  stop_with_report( paste0( "Step6B assignment table is missing module(s): ", paste( missing_modules, collapse = ", " ), "." ) )
}
shared_source_genes <- clean_gene_vector( step6b_assignments$gene_id[ step6b_assignments$module_name == SHARED_STEP6B_MODULE ] )
arm2_module_genes <- clean_gene_vector( step6b_assignments$gene_id[ step6b_assignments$module_name %in% ARM2_STEP6B_MODULES ] )
arm1_source_genes <- unique( step6a_arm1_genes )
arm2_source_genes <- unique( c( step6a_arm2_genes, arm2_module_genes ) )
source_overlap_before_shared_separation <- intersect( arm1_source_genes, arm2_source_genes )
residual_overlap_after_shared_separation <- setdiff( source_overlap_before_shared_separation, shared_source_genes )
arm1_specific_requested <- setdiff( arm1_source_genes, union( shared_source_genes, arm2_source_genes ) )
arm2_specific_requested <- setdiff( arm2_source_genes, union( shared_source_genes, arm1_source_genes ) )
shared_prebranch_requested <- unique( shared_source_genes )
if ( length( arm1_specific_requested ) < MIN_SIGNATURE_GENES_PRESENT || length( arm2_specific_requested ) < MIN_SIGNATURE_GENES_PRESENT || length( shared_prebranch_requested ) < MIN_SIGNATURE_GENES_PRESENT ) {
  stop_with_report( paste0( "A refined source-defined program is too small before external mapping. ", "Arm_1_specific=", length( arm1_specific_requested ), "; Arm_2_specific=", length( arm2_specific_requested ), "; Shared_prebranch=", length( shared_prebranch_requested ), "." ) )
}
if ( length( intersect( arm1_specific_requested, arm2_specific_requested ) ) > 0L || length( intersect( arm1_specific_requested, shared_prebranch_requested ) ) > 0L || length( intersect( arm2_specific_requested, shared_prebranch_requested ) ) > 0L ) {
  stop_with_report( "The refined Arm-specific and shared source sets are not mutually exclusive." )
}
available_features <- rownames( main_object )
arm1_resolution <- resolve_features( requested_features = arm1_specific_requested, available_features = available_features )
arm2_resolution <- resolve_features( requested_features = arm2_specific_requested, available_features = available_features )
shared_resolution <- resolve_features( requested_features = shared_prebranch_requested, available_features = available_features )
arm1_present <- unique( unname( arm1_resolution[ !is.na( arm1_resolution ) ] ) )
arm2_present <- unique( unname( arm2_resolution[ !is.na( arm2_resolution ) ] ) )
shared_present <- unique( unname( shared_resolution[ !is.na( shared_resolution ) ] ) )
arm1_mapping_fraction <- length( arm1_present ) / length( arm1_specific_requested )
arm2_mapping_fraction <- length( arm2_present ) / length( arm2_specific_requested )
shared_mapping_fraction <- length( shared_present ) / length( shared_prebranch_requested )
if ( length( arm1_present ) < MIN_SIGNATURE_GENES_PRESENT || length( arm2_present ) < MIN_SIGNATURE_GENES_PRESENT || length( shared_present ) < MIN_SIGNATURE_GENES_PRESENT || arm1_mapping_fraction < MIN_SIGNATURE_MAPPING_FRACTION || arm2_mapping_fraction < MIN_SIGNATURE_MAPPING_FRACTION || shared_mapping_fraction < MIN_SIGNATURE_MAPPING_FRACTION ) {
  stop_with_report( paste0( "Refined-program mapping failed. Arm_1_specific=", length( arm1_present ), "/", length( arm1_specific_requested ), " (", round( arm1_mapping_fraction, 3 ), "); Arm_2_specific=", length( arm2_present ), "/", length( arm2_specific_requested ), " (", round( arm2_mapping_fraction, 3 ), "); Shared_prebranch=", length( shared_present ), "/", length( shared_prebranch_requested ), " (", round( shared_mapping_fraction, 3 ), ")." ) )
}
if ( length( intersect( arm1_present, arm2_present ) ) > 0L || length( intersect( arm1_present, shared_present ) ) > 0L || length( intersect( arm2_present, shared_present ) ) > 0L ) {
  stop_with_report( "Mapped refined programs are not mutually exclusive." )
}
scored_gene_sets <- list( Arm_1_specific = arm1_present, Arm_2_specific = arm2_present, Shared_prebranch = shared_present )
all_source_genes <- sort( unique( c( arm1_source_genes, arm2_source_genes, shared_source_genes ) ) )
all_source_resolution <- resolve_features( requested_features = all_source_genes, available_features = available_features )
signature_gene_audit <- data.frame( record_type = "GENE", gene = all_source_genes, in_Step6A_Arm1_source = all_source_genes %in% step6a_arm1_genes, in_Step6A_Arm2_source = all_source_genes %in% step6a_arm2_genes, in_Step6B_Arm2_M1_M2_source = all_source_genes %in% arm2_module_genes, in_M3_Shared_prebranch = all_source_genes %in% shared_source_genes, in_source_Arm1_Arm2_overlap = all_source_genes %in% source_overlap_before_shared_separation, in_residual_overlap_not_scored = all_source_genes %in% residual_overlap_after_shared_separation, in_Arm1_specific = all_source_genes %in% arm1_specific_requested, in_Arm2_specific = all_source_genes %in% arm2_specific_requested, in_Shared_prebranch = all_source_genes %in% shared_prebranch_requested, resolved_external_feature = unname( all_source_resolution ), present_in_GSE234774 = !is.na( all_source_resolution ), technical_gene_pattern = grepl( TECHNICAL_GENE_PATTERN, all_source_genes, ignore.case = TRUE ), input_gene_count = NA_integer_, present_gene_count = NA_integer_, mapping_fraction = NA_real_, stringsAsFactors = FALSE )
signature_summary <- data.frame( record_type = "SUMMARY", gene = c( "Step6A_Arm1_source", "Step6A_Arm2_source", "Step6B_Arm2_M1_M2_source", "M3_Shared_prebranch", "Arm1_Arm2_source_overlap_before_shared_separation", "Residual_overlap_not_scored", "Arm_1_specific", "Arm_2_specific", "Shared_prebranch_scored" ), in_Step6A_Arm1_source = NA, in_Step6A_Arm2_source = NA, in_Step6B_Arm2_M1_M2_source = NA, in_M3_Shared_prebranch = NA, in_source_Arm1_Arm2_overlap = NA, in_residual_overlap_not_scored = NA, in_Arm1_specific = NA, in_Arm2_specific = NA, in_Shared_prebranch = NA, resolved_external_feature = "", present_in_GSE234774 = NA, technical_gene_pattern = NA, input_gene_count = c( length( step6a_arm1_genes ), length( step6a_arm2_genes ), length( arm2_module_genes ), length( shared_source_genes ), length( source_overlap_before_shared_separation ), length( residual_overlap_after_shared_separation ), length( arm1_specific_requested ), length( arm2_specific_requested ), length( shared_prebranch_requested ) ), present_gene_count = c( sum( !is.na( resolve_features( step6a_arm1_genes, available_features ) ) ), sum( !is.na( resolve_features( step6a_arm2_genes, available_features ) ) ), sum( !is.na( resolve_features( arm2_module_genes, available_features ) ) ), length( shared_present ), sum( !is.na( resolve_features( source_overlap_before_shared_separation, available_features ) ) ), sum( !is.na( resolve_features( residual_overlap_after_shared_separation, available_features ) ) ), length( arm1_present ), length( arm2_present ), length( shared_present ) ), mapping_fraction = c( mean( !is.na( resolve_features( step6a_arm1_genes, available_features ) ) ), mean( !is.na( resolve_features( step6a_arm2_genes, available_features ) ) ), mean( !is.na( resolve_features( arm2_module_genes, available_features ) ) ), shared_mapping_fraction, ifelse( length( source_overlap_before_shared_separation ) == 0L, NA_real_, mean( !is.na( resolve_features( source_overlap_before_shared_separation, available_features ) ) ) ), ifelse( length( residual_overlap_after_shared_separation ) == 0L, NA_real_, mean( !is.na( resolve_features( residual_overlap_after_shared_separation, available_features ) ) ) ), arm1_mapping_fraction, arm2_mapping_fraction, shared_mapping_fraction ), stringsAsFactors = FALSE )
signature_gene_audit <- rbind( signature_summary, signature_gene_audit )
write_csv_safe( signature_gene_audit, "04_refined_Arm_specific_shared_signature_audit.csv" )
log_message( "Refined mapping PASS: Arm_1_specific ", length( arm1_present ), "/", length( arm1_specific_requested ), "; Arm_2_specific ", length( arm2_present ), "/", length( arm2_specific_requested ), "; Shared_prebranch ", length( shared_present ), "/", length( shared_prebranch_requested ), "; residual overlap not scored=", length( residual_overlap_after_shared_separation ), "." )
Seurat::DefaultAssay( main_object ) <- "RNA"
rna_layers <- SeuratObject::Layers( object = main_object, assay = "RNA", search = NA )
counts_layers <- rna_layers[ grepl( "^counts($|\\.)", rna_layers ) ]
data_layers <- rna_layers[ grepl( "^data($|\\.)", rna_layers ) ]
if ( length( counts_layers ) != 1L || length( data_layers ) != 1L ) {
  stop_with_report( paste0( "Expected one RNA counts layer and one RNA data layer; detected ", length( counts_layers ), " and ", length( data_layers ), "." ) )
}
counts_layer_name <- counts_layers[[1L]]
data_layer_name <- data_layers[[1L]]
log_message( "Reading the refined main-analysis sparse RNA counts." )
counts_matrix_raw <- SeuratObject::LayerData( object = main_object, assay = "RNA", layer = counts_layer_name, fast = FALSE )
counts_matrix <- to_dgCMatrix( counts_matrix_raw, "refined main-analysis RNA counts" )
rm( counts_matrix_raw )
invisible( gc() )
if ( nrow( counts_matrix ) != nrow( main_object ) || ncol( counts_matrix ) != ncol( main_object ) || !identical( rownames( counts_matrix ), rownames( main_object ) ) || !identical( colnames( counts_matrix ), colnames( main_object ) ) ) {
  stop_with_report( "Sparse counts matrix does not match the refined main Seurat object." )
}
log_message( "Building AUCell rankings on the refined sparse counts matrix." )
auc_max_rank <- max( 1L, ceiling( AUC_MAX_RANK_FRACTION * nrow( counts_matrix ) ) )
set.seed( 20260803L )
cell_rankings <- AUCell::AUCell_buildRankings( exprMat = counts_matrix, featureType = "genes", plotStats = FALSE, splitByBlocks = TRUE, BPPARAM = BiocParallel::SerialParam(), keepZeroesAsNA = FALSE, verbose = TRUE )
log_message( "Calculating normalized AUCell AUC for Arm_1_specific, Arm_2_specific, and Shared_prebranch." )
cell_auc_object <- AUCell::AUCell_calcAUC( geneSets = scored_gene_sets, rankings = cell_rankings, nCores = 1L, normAUC = TRUE, aucMaxRank = auc_max_rank, verbose = TRUE )
cell_auc_matrix <- as.matrix( AUCell::getAUC( cell_auc_object ) )
if ( !all( SCORED_PROGRAM_NAMES %in% rownames( cell_auc_matrix ) ) || !identical( colnames( cell_auc_matrix ), colnames( main_object ) ) || any( !is.finite( cell_auc_matrix[ SCORED_PROGRAM_NAMES, , drop = FALSE ] ) ) ) {
  stop_with_report( "Refined AUCell output failed program, cell-order, or finite-value checks." )
}
main_object$Arm_1_specific_AUCell <- as.numeric( cell_auc_matrix[ "Arm_1_specific", colnames( main_object ) ] )
main_object$Arm_2_specific_AUCell <- as.numeric( cell_auc_matrix[ "Arm_2_specific", colnames( main_object ) ] )
main_object$Shared_prebranch_AUCell <- as.numeric( cell_auc_matrix[ "Shared_prebranch", colnames( main_object ) ] )
main_object$Arm_specific_balance_AUCell <- main_object$Arm_1_specific_AUCell - main_object$Arm_2_specific_AUCell
rm( cell_rankings, cell_auc_object )
invisible( gc() )
cell_score_table <- data.table::data.table( cell = colnames( main_object ), biological_library = as.character( main_object$biological_library ), timepoint = as.character( main_object$timepoint ), myeloid_final_label = as.character( main_object$myeloid_final_label ), Arm_1_specific_AUCell = main_object$Arm_1_specific_AUCell, Arm_2_specific_AUCell = main_object$Arm_2_specific_AUCell, Shared_prebranch_AUCell = main_object$Shared_prebranch_AUCell, Arm_specific_balance_AUCell = main_object$Arm_specific_balance_AUCell )
library_auc_summary <- cell_score_table[ , .( n_cells = .N, Arm_1_specific_AUCell_median = median( Arm_1_specific_AUCell ), Arm_1_specific_AUCell_mean = mean( Arm_1_specific_AUCell ), Arm_2_specific_AUCell_median = median( Arm_2_specific_AUCell ), Arm_2_specific_AUCell_mean = mean( Arm_2_specific_AUCell ), Shared_prebranch_AUCell_median = median( Shared_prebranch_AUCell ), Shared_prebranch_AUCell_mean = mean( Shared_prebranch_AUCell ), Arm_specific_balance_AUCell_median = median( Arm_specific_balance_AUCell ), Arm_specific_balance_AUCell_mean = mean( Arm_specific_balance_AUCell ) ), by = .( biological_library, timepoint ) ]
log_message( "Aggregating raw counts by retained biological_library." )
library_order_table <- unique( library_auc_summary[ , c( "biological_library", "timepoint" ), with = FALSE ] )
library_order_table$timepoint <- factor( library_order_table$timepoint, levels = TIMEPOINT_LEVELS, ordered = TRUE )
data.table::setorder( library_order_table, timepoint, biological_library )
library_order <- library_order_table$biological_library
cell_library_factor <- factor( as.character( main_object$biological_library ), levels = library_order )
if ( anyNA( cell_library_factor ) ) {
  stop_with_report( "At least one refined main-analysis cell could not be mapped to a retained pseudobulk library." )
}
library_indicator <- Matrix::sparseMatrix( i = seq_len( ncol( main_object ) ), j = as.integer( cell_library_factor ), x = 1, dims = c( ncol( main_object ), length( library_order ) ), dimnames = list( colnames( main_object ), library_order ) )
pseudobulk_counts_sparse <- counts_matrix %*%
    library_indicator
if ( ncol( pseudobulk_counts_sparse ) != length( library_order ) || !identical( colnames( pseudobulk_counts_sparse ), library_order ) ) {
  stop_with_report( "Refined pseudobulk aggregation failed the library-order check." )
}
pseudobulk_library_totals <- Matrix::colSums( pseudobulk_counts_sparse )
if ( any( pseudobulk_library_totals <= 0 ) ) {
  stop_with_report( "At least one retained pseudobulk library has zero total counts." )
}
pseudobulk_keep_genes <- Matrix::rowSums( pseudobulk_counts_sparse ) > 0
pseudobulk_counts <- as.matrix( pseudobulk_counts_sparse[ pseudobulk_keep_genes, , drop = FALSE ] )
rm( pseudobulk_counts_sparse, library_indicator )
invisible( gc() )
sample_manifest <- data.frame( biological_library = library_order, timepoint = as.character( library_order_table$timepoint[ match( library_order, library_order_table$biological_library ) ] ), stringsAsFactors = FALSE )
rownames( sample_manifest ) <- sample_manifest$biological_library
dge <- edgeR::DGEList( counts = pseudobulk_counts, samples = sample_manifest )
dge <- edgeR::calcNormFactors( dge, method = "TMM" )
pseudobulk_logcpm <- edgeR::cpm( dge, log = TRUE, prior.count = PSEUDOBULK_PRIOR_COUNT )
if ( any( !is.finite( pseudobulk_logcpm ) ) || !identical( colnames( pseudobulk_logcpm ), library_order ) ) {
  stop_with_report( "Refined TMM-normalized pseudobulk logCPM failed validation." )
}
ssgsea_gene_sets <- lapply( scored_gene_sets, function(genes) {
    intersect( genes, rownames( pseudobulk_logcpm ) )
  }
)
if ( any( vapply( ssgsea_gene_sets, length, integer(1) ) < MIN_SIGNATURE_GENES_PRESENT ) ) {
  stop_with_report( "A refined program has too few genes in the pseudobulk logCPM matrix." )
}
ssgsea_scores <- run_ssgsea_version_safe( expression_matrix = pseudobulk_logcpm, gene_sets = ssgsea_gene_sets )
if ( !all( SCORED_PROGRAM_NAMES %in% rownames( ssgsea_scores ) ) || !all( library_order %in% colnames( ssgsea_scores ) ) ) {
  stop_with_report( "Refined ssGSEA output does not contain all three programs and retained libraries." )
}
ssgsea_scores <- ssgsea_scores[ SCORED_PROGRAM_NAMES, library_order, drop = FALSE ]
if ( any( !is.finite( ssgsea_scores ) ) ) {
  stop_with_report( "Refined ssGSEA returned non-finite values." )
}
library_scores <- merge( as.data.frame( library_auc_summary ), data.frame( biological_library = library_order, Arm_1_specific_ssGSEA = as.numeric( ssgsea_scores[ "Arm_1_specific", library_order ] ), Arm_2_specific_ssGSEA = as.numeric( ssgsea_scores[ "Arm_2_specific", library_order ] ), Shared_prebranch_ssGSEA = as.numeric( ssgsea_scores[ "Shared_prebranch", library_order ] ), stringsAsFactors = FALSE ), by = "biological_library", all.x = TRUE, sort = FALSE )
library_scores <- library_scores[ match( library_order, library_scores$biological_library ), , drop = FALSE ]
library_scores$Arm_specific_balance_ssGSEA <- library_scores$Arm_1_specific_ssGSEA - library_scores$Arm_2_specific_ssGSEA
library_scores$timepoint <- factor( library_scores$timepoint, levels = TIMEPOINT_LEVELS, ordered = TRUE )
score_columns_to_standardize <- c( "Arm_1_specific_AUCell_median", "Arm_2_specific_AUCell_median", "Shared_prebranch_AUCell_median", "Arm_specific_balance_AUCell_median", "Arm_1_specific_ssGSEA", "Arm_2_specific_ssGSEA", "Shared_prebranch_ssGSEA", "Arm_specific_balance_ssGSEA" )
for (score_column in score_columns_to_standardize) {
  library_scores[ , paste0( score_column, "_z" ) ] <- safe_zscore( library_scores[[score_column]] )
}
concordance_inputs <- list( Arm_1_specific = c( "Arm_1_specific_AUCell_median", "Arm_1_specific_ssGSEA" ), Arm_2_specific = c( "Arm_2_specific_AUCell_median", "Arm_2_specific_ssGSEA" ), Shared_prebranch = c( "Shared_prebranch_AUCell_median", "Shared_prebranch_ssGSEA" ), Arm_specific_balance = c( "Arm_specific_balance_AUCell_median", "Arm_specific_balance_ssGSEA" ) )
concordance_rows <- lapply( names( concordance_inputs ), function(program_name) {
    columns <- concordance_inputs[[program_name]]
    test_result <- safe_cor_test( library_scores[[columns[[1L]]]], library_scores[[columns[[2L]]]] )
    data.frame( program = program_name, AUCell_library_statistic = "median_per_cell_normalized_AUC", ssGSEA_library_statistic = "TMM_logCPM_pseudobulk_ssGSEA", spearman_rho = test_result$rho, p_value = test_result$p_value, n_libraries = test_result$n, stringsAsFactors = FALSE )
  }
)
concordance_table <- do.call( rbind, concordance_rows )
rownames( concordance_table ) <- NULL
library_long_parts <- list( data.frame( biological_library = library_scores$biological_library, timepoint = as.character( library_scores$timepoint ), n_cells = library_scores$n_cells, method = "AUCell_library_median", program = "Arm_1_specific", score = library_scores$Arm_1_specific_AUCell_median, score_z = library_scores$Arm_1_specific_AUCell_median_z, stringsAsFactors = FALSE ), data.frame( biological_library = library_scores$biological_library, timepoint = as.character( library_scores$timepoint ), n_cells = library_scores$n_cells, method = "AUCell_library_median", program = "Arm_2_specific", score = library_scores$Arm_2_specific_AUCell_median, score_z = library_scores$Arm_2_specific_AUCell_median_z, stringsAsFactors = FALSE ), data.frame( biological_library = library_scores$biological_library, timepoint = as.character( library_scores$timepoint ), n_cells = library_scores$n_cells, method = "AUCell_library_median", program = "Shared_prebranch", score = library_scores$Shared_prebranch_AUCell_median, score_z = library_scores$Shared_prebranch_AUCell_median_z, stringsAsFactors = FALSE ), data.frame( biological_library = library_scores$biological_library, timepoint = as.character( library_scores$timepoint ), n_cells = library_scores$n_cells, method = "AUCell_library_median", program = "Arm_specific_balance", score = library_scores$Arm_specific_balance_AUCell_median, score_z = library_scores$Arm_specific_balance_AUCell_median_z, stringsAsFactors = FALSE ), data.frame( biological_library = library_scores$biological_library, timepoint = as.character( library_scores$timepoint ), n_cells = library_scores$n_cells, method = "ssGSEA_pseudobulk", program = "Arm_1_specific", score = library_scores$Arm_1_specific_ssGSEA, score_z = library_scores$Arm_1_specific_ssGSEA_z, stringsAsFactors = FALSE ), data.frame( biological_library = library_scores$biological_library, timepoint = as.character( library_scores$timepoint ), n_cells = library_scores$n_cells, method = "ssGSEA_pseudobulk", program = "Arm_2_specific", score = library_scores$Arm_2_specific_ssGSEA, score_z = library_scores$Arm_2_specific_ssGSEA_z, stringsAsFactors = FALSE ), data.frame( biological_library = library_scores$biological_library, timepoint = as.character( library_scores$timepoint ), n_cells = library_scores$n_cells, method = "ssGSEA_pseudobulk", program = "Shared_prebranch", score = library_scores$Shared_prebranch_ssGSEA, score_z = library_scores$Shared_prebranch_ssGSEA_z, stringsAsFactors = FALSE ), data.frame( biological_library = library_scores$biological_library, timepoint = as.character( library_scores$timepoint ), n_cells = library_scores$n_cells, method = "ssGSEA_pseudobulk", program = "Arm_specific_balance", score = library_scores$Arm_specific_balance_ssGSEA, score_z = library_scores$Arm_specific_balance_ssGSEA_z, stringsAsFactors = FALSE ) )
library_score_long <- do.call( rbind, library_long_parts )
library_score_long$timepoint <- factor( library_score_long$timepoint, levels = TIMEPOINT_LEVELS, ordered = TRUE )
timepoint_score_summary <- data.table::as.data.table( library_score_long )[ , .( n_libraries = uniqueN( biological_library ), total_cells = sum( n_cells ), mean_score = mean( score ), median_score = median( score ), minimum_score = min( score ), maximum_score = max( score ), sd_score = stats::sd( score ), mean_score_z = mean( score_z ), median_score_z = median( score_z ) ), by = .( timepoint, method, program ) ]
library_rows <- data.frame( record_type = "LIBRARY_SCORE", biological_library = library_score_long$biological_library, timepoint = as.character( library_score_long$timepoint ), n_cells = library_score_long$n_cells, n_libraries = 1L, method = library_score_long$method, program = library_score_long$program, score = library_score_long$score, score_z = library_score_long$score_z, mean_score = NA_real_, median_score = NA_real_, minimum_score = NA_real_, maximum_score = NA_real_, sd_score = NA_real_, mean_score_z = NA_real_, median_score_z = NA_real_, stringsAsFactors = FALSE )
timepoint_rows <- data.frame( record_type = "TIMEPOINT_SUMMARY", biological_library = "ALL", timepoint = as.character( timepoint_score_summary$timepoint ), n_cells = timepoint_score_summary$total_cells, n_libraries = timepoint_score_summary$n_libraries, method = timepoint_score_summary$method, program = timepoint_score_summary$program, score = timepoint_score_summary$median_score, score_z = timepoint_score_summary$median_score_z, mean_score = timepoint_score_summary$mean_score, median_score = timepoint_score_summary$median_score, minimum_score = timepoint_score_summary$minimum_score, maximum_score = timepoint_score_summary$maximum_score, sd_score = timepoint_score_summary$sd_score, mean_score_z = timepoint_score_summary$mean_score_z, median_score_z = timepoint_score_summary$median_score_z, stringsAsFactors = FALSE )
combined_score_table <- rbind( library_rows, timepoint_rows )
combined_score_table$timepoint <- factor( combined_score_table$timepoint, levels = TIMEPOINT_LEVELS, ordered = TRUE )
combined_score_table <- combined_score_table[ order( combined_score_table$record_type, combined_score_table$method, combined_score_table$program, combined_score_table$timepoint, combined_score_table$biological_library ), , drop = FALSE ]
combined_score_table$timepoint <- as.character( combined_score_table$timepoint )
write_csv_safe( combined_score_table, "05_refined_main_library_and_timepoint_scores.csv" )
write_csv_safe( library_scores, "06_refined_library_scores_wide.csv" )
write_csv_safe( concordance_table, "07_refined_AUCell_ssGSEA_concordance.csv" )
log_message( "Creating refined AUCell UMAP and library-level figures." )
umap_embeddings <- Seurat::Embeddings( object = main_object, reduction = FROZEN_UMAP )
umap_plot_data <- data.frame( UMAP_1 = umap_embeddings[ , 1L ], UMAP_2 = umap_embeddings[ , 2L ], Arm_1_specific_AUCell = main_object$Arm_1_specific_AUCell, Arm_2_specific_AUCell = main_object$Arm_2_specific_AUCell, Shared_prebranch_AUCell = main_object$Shared_prebranch_AUCell, Arm_specific_balance_AUCell = main_object$Arm_specific_balance_AUCell, timepoint = main_object$timepoint, stringsAsFactors = FALSE )
p_arm1_umap <- plot_umap_score( plot_data = umap_plot_data, score_column = "Arm_1_specific_AUCell", title_text = "Arm_1-specific AUCell", legend_title = "Normalized AUC" )
p_arm2_umap <- plot_umap_score( plot_data = umap_plot_data, score_column = "Arm_2_specific_AUCell", title_text = "Arm_2-specific AUCell", legend_title = "Normalized AUC" )
p_shared_umap <- plot_umap_score( plot_data = umap_plot_data, score_column = "Shared_prebranch_AUCell", title_text = "Shared prebranch AUCell", legend_title = "Normalized AUC" )
p_balance_umap <- plot_umap_score( plot_data = umap_plot_data, score_column = "Arm_specific_balance_AUCell", title_text = "Arm-specific AUCell balance", legend_title = "Arm_1 - Arm_2" )
ggplot2::ggsave( filename = file.path( FIGURES_DIR, "01_refined_AUCell_specific_shared_UMAP.pdf" ), plot = ( p_arm1_umap | p_arm2_umap ) / ( p_shared_umap | p_balance_umap ), width = 19, height = 14, limitsize = FALSE )
plot_library_long <- library_score_long
plot_library_long$timepoint <- factor( as.character( plot_library_long$timepoint ), levels = TIMEPOINT_LEVELS, ordered = TRUE )
p_library_timecourse <- ggplot2::ggplot( plot_library_long, ggplot2::aes( x = timepoint, y = score, group = biological_library ) ) + ggplot2::geom_point( ggplot2::aes( shape = biological_library ), size = 2.4, position = ggplot2::position_jitter( width = 0.08, height = 0 ) ) + ggplot2::stat_summary( ggplot2::aes( group = 1 ), fun = median, geom = "line", linewidth = 0.8 ) + ggplot2::stat_summary( ggplot2::aes( group = 1 ), fun = median, geom = "point", size = 3.0 ) + ggplot2::facet_grid( rows = ggplot2::vars( method ), cols = ggplot2::vars( program ), scales = "free_y" ) + ggplot2::labs( title = "Refined Arm-specific and shared scores by biological library", subtitle = "Points are retained biological libraries; line is the timepoint median", x = NULL, y = "Method-specific score", shape = "Biological library" ) + ggplot2::theme_bw() + ggplot2::theme( axis.text.x = ggplot2::element_text( angle = 45, hjust = 1 ), legend.position = "none" )
ggplot2::ggsave( filename = file.path( FIGURES_DIR, "02_refined_library_AUCell_ssGSEA_timecourse.pdf" ), plot = p_library_timecourse, width = 22, height = 10, limitsize = FALSE )
concordance_plot_parts <- lapply( names( concordance_inputs ), function(program_name) {
    columns <- concordance_inputs[[program_name]]
    data.frame( biological_library = library_scores$biological_library, timepoint = library_scores$timepoint, program = program_name, AUCell_z = library_scores[[ paste0( columns[[1L]], "_z" ) ]], ssGSEA_z = library_scores[[ paste0( columns[[2L]], "_z" ) ]], stringsAsFactors = FALSE )
  }
)
concordance_plot_data <- do.call( rbind, concordance_plot_parts )
concordance_labels <- data.frame( program = concordance_table$program, label = paste0( "Spearman rho=", ifelse( is.na( concordance_table$spearman_rho ), "NA", format( round( concordance_table$spearman_rho, 3 ), nsmall = 3 ) ), "\nn=", concordance_table$n_libraries ), x = -Inf, y = Inf, stringsAsFactors = FALSE )
p_concordance <- ggplot2::ggplot( concordance_plot_data, ggplot2::aes( x = AUCell_z, y = ssGSEA_z, shape = timepoint ) ) + ggplot2::geom_hline( yintercept = 0, linetype = 3 ) + ggplot2::geom_vline( xintercept = 0, linetype = 3 ) + ggplot2::geom_point( size = 2.8 ) + ggplot2::geom_smooth( method = "lm", se = FALSE, linewidth = 0.7 ) + ggplot2::geom_text( data = concordance_labels, ggplot2::aes( x = x, y = y, label = label ), inherit.aes = FALSE, hjust = -0.05, vjust = 1.1, size = 3.2 ) + ggplot2::facet_wrap( ggplot2::vars( program ), scales = "free" ) + ggplot2::labs( title = "Refined library-level AUCell and pseudobulk ssGSEA concordance", subtitle = "Scores are standardized separately within each method", x = "AUCell library-median Z score", y = "Pseudobulk ssGSEA Z score", shape = "Timepoint" ) + ggplot2::theme_bw()
ggplot2::ggsave( filename = file.path( FIGURES_DIR, "03_refined_AUCell_ssGSEA_concordance.pdf" ), plot = p_concordance, width = 18, height = 7, limitsize = FALSE )
main_object@misc$GSE234774_step11B <- list( main_timepoints = MAIN_TIMEPOINTS, minimum_cells_per_library = MIN_CELLS_PER_LIBRARY, included_libraries = included_libraries, excluded_low_cell_libraries = excluded_libraries, statistical_unit = "biological_library", source_based_refinement = "M3 separated as shared; Arm-specific sets made mutually exclusive", Arm1_source = "Step6A Arm1 high-confidence", Arm2_source = "Step6A Arm2 high-confidence union Step6B M1/M2 Arm2 terminal", Shared_source = "Step6B M3 Shared prebranch", Arm1_specific_requested = arm1_specific_requested, Arm2_specific_requested = arm2_specific_requested, Shared_prebranch_requested = shared_prebranch_requested, residual_overlap_not_scored = residual_overlap_after_shared_separation, Arm1_specific_present = arm1_present, Arm2_specific_present = arm2_present, Shared_prebranch_present = shared_present, AUCell_aucMaxRank = auc_max_rank, AUCell_normAUC = TRUE, pseudobulk_normalization = "edgeR TMM logCPM", pseudobulk_prior_count = PSEUDOBULK_PRIOR_COUNT, ssGSEA_alpha = SSGSEA_ALPHA, cell_level_p_values_run = FALSE, sensitivity_analysis_run = FALSE, cycling_analysis_run = FALSE, completed_at = as.character( Sys.time() ) )
log_message( "Saving the refined scored Seurat object." )
saveRDS( main_object, SCORED_OBJECT_PATH, compress = FALSE )
results_object <- list( main_timepoints = MAIN_TIMEPOINTS, minimum_cells_per_library = MIN_CELLS_PER_LIBRARY, library_QC = as.data.frame( library_QC_audit ), included_libraries = included_libraries, excluded_libraries = excluded_libraries, source_sets = list( Arm1_source = arm1_source_genes, Arm2_source = arm2_source_genes, Shared_source = shared_source_genes, residual_overlap_not_scored = residual_overlap_after_shared_separation ), scored_gene_sets_requested = list( Arm_1_specific = arm1_specific_requested, Arm_2_specific = arm2_specific_requested, Shared_prebranch = shared_prebranch_requested ), scored_gene_sets_present = scored_gene_sets, AUCell_aucMaxRank = auc_max_rank, cell_AUCell = cell_auc_matrix, library_scores = library_scores, pseudobulk_counts = pseudobulk_counts, pseudobulk_logCPM = pseudobulk_logcpm, ssGSEA_scores = ssgsea_scores, concordance = concordance_table, timepoint_summary = as.data.frame( timepoint_score_summary ), completed_at = as.character( Sys.time() ) )
saveRDS( results_object, RESULTS_RDS_PATH, compress = TRUE )
reloaded_scored_object <- tryCatch( readRDS( SCORED_OBJECT_PATH ), error = function(error_object) {
    error_object
  }
)
scored_object_reload_pass <- !inherits( reloaded_scored_object, "error" ) && inherits( reloaded_scored_object, "Seurat" ) && ncol( reloaded_scored_object ) == ncol( main_object ) && identical( colnames( reloaded_scored_object ), colnames( main_object ) ) && all( c( "Arm_1_specific_AUCell", "Arm_2_specific_AUCell", "Shared_prebranch_AUCell", "Arm_specific_balance_AUCell" ) %in% colnames( reloaded_scored_object@meta.data ) )
rm( reloaded_scored_object )
invisible( gc() )
reloaded_results <- tryCatch( readRDS( RESULTS_RDS_PATH ), error = function(error_object) {
    error_object
  }
)
results_reload_pass <- !inherits( reloaded_results, "error" ) && is.list( reloaded_results ) && all( c( "library_scores", "ssGSEA_scores", "cell_AUCell", "scored_gene_sets_present", "library_QC" ) %in% names( reloaded_results ) )
rm( reloaded_results )
invisible( gc() )
required_figure_paths <- file.path( FIGURES_DIR, c( "01_refined_AUCell_specific_shared_UMAP.pdf", "02_refined_library_AUCell_ssGSEA_timecourse.pdf", "03_refined_AUCell_ssGSEA_concordance.pdf" ) )
required_table_paths <- file.path( TABLES_DIR, c( "04_refined_Arm_specific_shared_signature_audit.csv", "05_refined_main_library_and_timepoint_scores.csv" ) )
program_mapping_pass <- length( arm1_present ) >= MIN_SIGNATURE_GENES_PRESENT && length( arm2_present ) >= MIN_SIGNATURE_GENES_PRESENT && length( shared_present ) >= MIN_SIGNATURE_GENES_PRESENT && arm1_mapping_fraction >= MIN_SIGNATURE_MAPPING_FRACTION && arm2_mapping_fraction >= MIN_SIGNATURE_MAPPING_FRACTION && shared_mapping_fraction >= MIN_SIGNATURE_MAPPING_FRACTION
mutual_exclusivity_pass <- length( intersect( arm1_present, arm2_present ) ) == 0L && length( intersect( arm1_present, shared_present ) ) == 0L && length( intersect( arm2_present, shared_present ) ) == 0L
audit_checks <- data.frame( check = c( "Frozen Step 10G primary object dimensions and cluster membership are correct", "Only the five prespecified main validation timepoints are included", "Libraries below the 50-cell threshold were excluded before scoring", "Every retained main timepoint has at least two biological libraries", "IFN-responsive cluster 7 is retained and cycling cluster 6 is absent", "M3_Shared_prebranch is scored independently rather than merged into Arm_1", "Arm_1-specific, Arm_2-specific, and shared mapped programs are mutually exclusive", "All three refined programs meet mapping thresholds", "AUCell returned finite scores for every retained main-analysis cell", "Pseudobulk produced one nonzero profile per retained biological library", "ssGSEA returned finite scores for all programs and retained libraries", "Library-level AUCell-ssGSEA concordance statistics were produced", "All three required figures and two core tables were created", "Scored Seurat object and complete results RDS reload successfully", "Cell-level P values were calculated", "Non-IFN sensitivity, 1 dpi, 4 dpi, or cycling analyses were run" ), required_for_pass = c( TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE ), status = c( "PASS", ifelse( setequal( unique( as.character( main_object$timepoint ) ), MAIN_TIMEPOINTS ), "PASS", "STOP" ), ifelse( all( library_QC_audit$n_cells_after_library_QC[ library_QC_audit$library_QC_status == "EXCLUDED_LOW_CELL_COUNT" ] == 0L ) && all( library_QC_audit$n_cells_before_library_QC[ library_QC_audit$library_QC_status == "INCLUDED" ] >= MIN_CELLS_PER_LIBRARY ), "PASS", "STOP" ), ifelse( all( postQC_timepoint_replication$n_libraries >= MIN_LIBRARIES_PER_MAIN_TIMEPOINT ), "PASS", "STOP" ), ifelse( "IFN_responsive_microglia" %in% as.character( main_object$myeloid_final_label ) && !"Cycling_microglia" %in% as.character( main_object$myeloid_final_label ), "PASS", "STOP" ), ifelse( identical( sort( shared_prebranch_requested ), sort( shared_source_genes ) ) && length( intersect( arm1_specific_requested, shared_prebranch_requested ) ) == 0L, "PASS", "STOP" ), ifelse( mutual_exclusivity_pass, "PASS", "STOP" ), ifelse( program_mapping_pass, "PASS", "STOP" ), ifelse( ncol( cell_auc_matrix ) == ncol( main_object ) && all( is.finite( cell_auc_matrix[ SCORED_PROGRAM_NAMES, , drop = FALSE ] ) ), "PASS", "STOP" ), ifelse( ncol( pseudobulk_counts ) == length( library_order ) && all( colSums( pseudobulk_counts ) > 0 ), "PASS", "STOP" ), ifelse( ncol( ssgsea_scores ) == length( library_order ) && all( is.finite( ssgsea_scores ) ), "PASS", "STOP" ), ifelse( nrow( concordance_table ) == 4L && all( concordance_table$n_libraries == length( library_order ) ), "PASS", "STOP" ), ifelse( all( file.exists( required_figure_paths ) ) && all( file.exists( required_table_paths ) ), "PASS", "STOP" ), ifelse( scored_object_reload_pass && results_reload_pass, "PASS", "STOP" ), "CONFIRMED_NOT_RUN", "CONFIRMED_NOT_RUN" ), detail = c( paste0( nrow( primary_object ), " features x ", ncol( primary_object ), " nuclei; clusters=", paste( input_cluster_values, collapse = "," ) ), paste( MAIN_TIMEPOINTS, collapse = ", " ), paste0( "threshold=", MIN_CELLS_PER_LIBRARY, "; included libraries=", length( included_libraries ), "; excluded libraries=", ifelse( length( excluded_libraries ) == 0L, "none", paste( excluded_libraries, collapse = "," ) ) ), paste( paste0( as.character( postQC_timepoint_replication$timepoint ), "=", postQC_timepoint_replication$n_libraries, " libraries" ), collapse = "; " ), paste0( "labels=", paste( sort( unique( as.character( main_object$myeloid_final_label ) ) ), collapse = ";" ) ), paste0( "Shared source=", SHARED_STEP6B_MODULE, "; shared genes=", length( shared_prebranch_requested ) ), paste0( "mapped intersections: Arm1/Arm2=", length( intersect( arm1_present, arm2_present ) ), "; Arm1/shared=", length( intersect( arm1_present, shared_present ) ), "; Arm2/shared=", length( intersect( arm2_present, shared_present ) ) ), paste0( "Arm1 ", length( arm1_present ), "/", length( arm1_specific_requested ), " (", round( arm1_mapping_fraction, 3 ), "); Arm2 ", length( arm2_present ), "/", length( arm2_specific_requested ), " (", round( arm2_mapping_fraction, 3 ), "); Shared ", length( shared_present ), "/", length( shared_prebranch_requested ), " (", round( shared_mapping_fraction, 3 ), ")" ), paste0( ncol( cell_auc_matrix ), " cells; aucMaxRank=", auc_max_rank ), paste0( ncol( pseudobulk_counts ), " libraries; total-count range=", min( colSums( pseudobulk_counts ) ), "-", max( colSums( pseudobulk_counts ) ) ), paste0( ncol( ssgsea_scores ), " libraries; GSVA version=", safe_package_version( "GSVA" ) ), paste( paste0( concordance_table$program, " rho=", round( concordance_table$spearman_rho, 3 ) ), collapse = "; " ), paste( c( basename( required_figure_paths ), basename( required_table_paths ) ), collapse = "; " ), paste0( "scored object reload=", scored_object_reload_pass, "; results RDS reload=", results_reload_pass ), "No cell-level inferential test was run.", "Main timepoints only; no sensitivity or cycling analysis." ), stringsAsFactors = FALSE )
write_csv_safe( audit_checks, "09_STEP11B_audit_checks.csv" )
required_failures <- audit_checks$required_for_pass & audit_checks$status != "PASS"
final_decision <- if ( any( required_failures ) ) {
  "STOP"
} else {
  "PASS_WITH_REFINED_MAIN_ARM_VALIDATION_REVIEW"
}
timepoint_lines <- apply( as.data.frame( postQC_timepoint_replication ), 1L, function(row) {
    paste0( row[["timepoint"]], ": ", row[["n_libraries"]], " retained libraries; ", row[["n_cells"]], " nuclei; library-cell range ", row[["minimum_library_cells"]], "-", row[["maximum_library_cells"]] )
  }
)
concordance_lines <- apply( concordance_table, 1L, function(row) {
    paste0( row[["program"]], ": Spearman rho=", ifelse( is.na( as.numeric( row[["spearman_rho"]] ) ), "NA", round( as.numeric( row[["spearman_rho"]] ), 3 ) ), "; p=", ifelse( is.na( as.numeric( row[["p_value"]] ) ), "NA", signif( as.numeric( row[["p_value"]] ), 4 ) ), "; n=", row[["n_libraries"]] )
  }
)
report_lines <- c( "GSE234774 STEP 11B REFINED MAIN ARM VALIDATION", "============================================================", paste0( "Run time: ", format( Sys.time(), "%Y-%m-%d %H:%M:%S" ) ), paste0( "Primary input: ", PRIMARY_OBJECT_PATH ), paste0( "Output directory: ", OUTPUT_DIR ), "", paste0( "FINAL DECISION: ", final_decision ), "", "FORMAL SCOPE", "Main validation timepoints only: Uninjured, 7 dpi, 14 dpi, 1 month, 2 months.", "Primary non-cycling resident microglia only; IFN-responsive cluster 7 retained.", paste0( "Minimum retained cells per biological library: ", MIN_CELLS_PER_LIBRARY, "." ), paste0( "Excluded low-cell libraries: ", ifelse( length( excluded_libraries ) == 0L, "none", paste( excluded_libraries, collapse = ", " ) ), "." ), "No 1 dpi, 4 dpi, non-IFN sensitivity, or cycling analysis was run.", "The statistical unit is biological_library.", "", "REFINED SOURCE-DEFINED PROGRAMS", "Shared_prebranch = Step6B M3_Shared_prebranch.", "Arm_1_specific = Step6A Arm_1-biased high-confidence genes minus shared and Arm_2-source genes.", "Arm_2_specific = Step6A Arm_2-biased high-confidence genes union Step6B M1/M2 Arm_2-terminal genes, minus shared and Arm_1-source genes.", "Residual Arm_1/Arm_2 source overlap after shared separation was audited and not scored.", paste0( "Residual overlap not scored: ", length( residual_overlap_after_shared_separation ), " genes." ), paste0( "Arm_1_specific mapped genes: ", length( arm1_present ), "/", length( arm1_specific_requested ), " (", round( 100 * arm1_mapping_fraction, 2 ), "%)." ), paste0( "Arm_2_specific mapped genes: ", length( arm2_present ), "/", length( arm2_specific_requested ), " (", round( 100 * arm2_mapping_fraction, 2 ), "%)." ), paste0( "Shared_prebranch mapped genes: ", length( shared_present ), "/", length( shared_prebranch_requested ), " (", round( 100 * shared_mapping_fraction, 2 ), "%)." ), "", "MAIN ANALYSIS COVERAGE AFTER LIBRARY QC", timepoint_lines, paste0( "Total retained nuclei: ", ncol( main_object ) ), paste0( "Total retained biological libraries: ", length( library_order ) ), "", "AUCELL", "Input: sparse single-cell RNA counts.", paste0( "AUCell aucMaxRank: ", auc_max_rank, " genes (", 100 * AUC_MAX_RANK_FRACTION, "% of ranked features)." ), "Normalized AUC was calculated per cell for Arm_1_specific, Arm_2_specific, and Shared_prebranch.", "Formal library-level AUCell summary: median per-cell AUC.", "", "PSEUDOBULK SSGSEA", "Raw counts were summed within each retained biological_library.", paste0( "Pseudobulk profiles: ", ncol( pseudobulk_counts ), "." ), paste0( "Normalization: edgeR TMM logCPM with prior.count=", PSEUDOBULK_PRIOR_COUNT, "." ), paste0( "ssGSEA alpha=", SSGSEA_ALPHA, "; final ssGSEA normalization enabled." ), "", "LIBRARY-LEVEL METHOD CONCORDANCE", concordance_lines, "", "INTERPRETATION BOUNDARIES", "1. Technical PASS does not prove biological replication; time directions and between-library consistency require review.", "2. Cell-level AUCell distributions are descriptive; biological interpretation uses library summaries.", "3. AUCell and ssGSEA raw magnitudes are method-specific and are not directly compared.", "4. Arm-specific balance is Arm_1_specific minus Arm_2_specific within each method.", "5. M3_Shared_prebranch is interpreted independently and is not evidence for either specific Arm.", "6. The refinement is based only on frozen source labels and the prespecified library-cell threshold.", "7. Arm_1 remains resolution/homeostasis-biased metabolic/lipid handling.", "8. Arm_2 remains inflammatory immune/vascular remodeling and scar-modulating; it is not a terminal-fate claim.", "", "PRIMARY OUTPUT OBJECTS", paste0( "Scored Seurat object: ", SCORED_OBJECT_PATH ), paste0( "Complete results RDS: ", RESULTS_RDS_PATH ), "", "NEXT DECISION AFTER SEVEN-FILE REVIEW", "Determine whether the refined Arm_1-specific and Arm_2-specific programs show concordant AUCell and pseudobulk ssGSEA directions.", "Evaluate Shared_prebranch separately and do not use it to support either Arm.", "Judge the external validation as replicated, partially replicated, sample-sensitive, or not replicated." )
if ( final_decision == "STOP" ) {
  failed <- audit_checks[ audit_checks$required_for_pass & audit_checks$status != "PASS", , drop = FALSE ]
  report_lines <- c( report_lines, "", "STOP REASONS" )
  for ( index in seq_len( nrow( failed ) ) ) {
    report_lines <- c( report_lines, paste0( index, ". ", failed$check[ index ], " -- ", failed$detail[ index ] ) )
  }
}
report_path <- file.path( OUTPUT_DIR, "10_STEP11B_FINAL_REPORT.txt" )
writeLines( report_lines, report_path, useBytes = TRUE )
log_message( "Final decision: ", final_decision )
log_message( "Final report: ", report_path )
cat("\n")
cat("============================================================\n")
cat("GSE234774 STEP 11B COMPLETED\n")
cat( "FINAL DECISION:", final_decision, "\n" )
cat( "Output:", OUTPUT_DIR, "\n" )
cat("Key outputs:\n")
cat(" 1) figures/01_refined_AUCell_specific_shared_UMAP.pdf\n")
cat(" 2) figures/02_refined_library_AUCell_ssGSEA_timecourse.pdf\n")
cat(" 3) figures/03_refined_AUCell_ssGSEA_concordance.pdf\n")
cat(" 4) tables/04_refined_Arm_specific_shared_signature_audit.csv\n")
cat(" 5) tables/05_refined_main_library_and_timepoint_scores.csv\n")
cat(" 6) tables/09_STEP11B_audit_checks.csv\n")
cat(" 7) 10_STEP11B_FINAL_REPORT.txt\n")
cat("============================================================\n")
if ( final_decision == "STOP" ) {
  stop( "GSE234774 Step 11B decision is STOP. Review the final report.", call. = FALSE )
}
