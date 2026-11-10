options( stringsAsFactors = FALSE, warn = 1, width = 180, future.globals.maxSize = 24 * 1024^3 )
PROJECT_ROOT <- Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
CAPTURE_AREAS_TO_RUN <- "ALL"
MAX_CORES <- 6L
RESUME_COMPLETED <- TRUE
FORCE_RERUN_CAPTURE_AREAS <- character()
SAVE_FULL_RCTD_OBJECT <- FALSE
WEIGHT_SUM_TOLERANCE <- 1e-6
NEGATIVE_WEIGHT_TOLERANCE <- 1e-10
STEP12C_DIR <- file.path( PROJECT_ROOT, "04_results", "step12C_GSE234774_snRNA_reference_RCTD_input_audit" )
STEP12C_REPORT <- file.path( STEP12C_DIR, "10_STEP12C_FINAL_REPORT.txt" )
REFERENCE_RDS <- file.path( STEP12C_DIR, "RDS", "GSE234774_step12C_RCTD_reference.rds" )
REFERENCE_METADATA_RDS <- file.path( STEP12C_DIR, "RDS", "GSE234774_step12C_RCTD_reference_cell_metadata.rds" )
STEP12C_QUERY_INDEX <- file.path( STEP12C_DIR, "tables", "05_RCTD_query_capture_area_index.csv" )
STEP12B_DIR <- file.path( PROJECT_ROOT, "04_results", "step12B_fast_BPCells_spatial_section_QC" )
SPATIAL_BPCELLS_DIR <- file.path( STEP12B_DIR, "BPCells", "GSE234774_spatial_2d_full_counts" )
STEP12B_SPOT_QC <- file.path( STEP12B_DIR, "tables", "04_spot_QC_flags.csv.gz" )
STEP12A2_ASSET_DIR <- file.path( PROJECT_ROOT, "04_results", "step12A2_GSE234774_spatial_metadata_asset_linkage_repair", "audit_extracted_2d_assets" )
OUTPUT_DIR <- file.path( PROJECT_ROOT, "04_results", "step12D_GSE234774_full_mode_RCTD_per_capture_area" )
CAPTURE_ROOT <- file.path( OUTPUT_DIR, "capture_areas" )
TABLE_DIR <- file.path( OUTPUT_DIR, "tables" )
RDS_DIR <- file.path( OUTPUT_DIR, "RDS" )
RUN_LOG_FILE <- file.path( OUTPUT_DIR, "00_run_log.txt" )
SESSION_INFO_FILE <- file.path( OUTPUT_DIR, "00_sessionInfo.txt" )
OUT_INPUT_AUDIT <- file.path( TABLE_DIR, "01_step12D_input_and_parameter_audit.csv" )
OUT_RUN_SUMMARY <- file.path( TABLE_DIR, "02_capture_area_RCTD_run_summary.csv" )
OUT_NORMALIZED_WEIGHTS <- file.path( TABLE_DIR, "03_all_capture_areas_normalized_weights.csv.gz" )
OUT_RAW_WEIGHTS <- file.path( TABLE_DIR, "04_all_capture_areas_raw_weights.csv.gz" )
OUT_SECTION_SUMMARY <- file.path( TABLE_DIR, "05_cell_type_weight_summary_by_section.csv" )
OUT_CHECKS <- file.path( TABLE_DIR, "09_STEP12D_audit_checks.csv" )
OUT_FINAL_REPORT <- file.path( OUTPUT_DIR, "10_STEP12D_FINAL_REPORT.txt" )
OUT_NORMALIZED_WEIGHTS_RDS <- file.path( RDS_DIR, "GSE234774_step12D_all_spots_normalized_RCTD_weights.rds" )
OUT_RAW_WEIGHTS_RDS <- file.path( RDS_DIR, "GSE234774_step12D_all_spots_raw_RCTD_weights.rds" )
OUT_HANDOFF_RDS <- file.path( RDS_DIR, "GSE234774_step12D_RCTD_weights_handoff.rds" )
for (directory in c( OUTPUT_DIR, CAPTURE_ROOT, TABLE_DIR, RDS_DIR )) {
  dir.create( directory, recursive = TRUE, showWarnings = FALSE )
}
required_packages <- c( "data.table", "Matrix", "BPCells", "spacexr" )
missing_packages <- required_packages[ !vapply( required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1) ) ]
if (length(missing_packages) > 0L) {
  stop( "Missing required package(s): ", paste( missing_packages, collapse = ", " ), call. = FALSE )
}
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(BPCells)
  library(spacexr)
})
run_log <- character()
log_message <- function(...) {
  line <- paste0( "[", format( Sys.time(), "%Y-%m-%d %H:%M:%S" ), "] ", paste0(...) )
  cat(line, "\n")
  run_log <<- c( run_log, line )
}
write_run_log <- function() {
  writeLines( run_log, RUN_LOG_FILE, useBytes = TRUE )
}
collapse_unique <- function( x, max_values = 40L) {
  x <- trimws( as.character(x) )
  x <- unique( x[ !is.na(x) & nzchar(x) ] )
  if (length(x) == 0L) {
    return(NA_character_)
  }
  if (length(x) > max_values) {
    x <- c( x[seq_len(max_values)], paste0( "... +", length(x) - max_values, " more" ) )
  }
  paste( x, collapse = " | " )
}
required_columns <- function( table, columns, label) {
  missing <- setdiff( columns, names(table) )
  if (length(missing) > 0L) {
    stop( label, " is missing required column(s): ", paste( missing, collapse = ", " ), call. = FALSE )
  }
  invisible(TRUE)
}
as_logical_strict <- function( x, label) {
  if (is.logical(x)) {
    return(x)
  }
  value <- toupper( trimws( as.character(x) ) )
  output <- rep( NA, length(value) )
  output[ value %in% c( "TRUE", "T", "1" ) ] <- TRUE
  output[ value %in% c( "FALSE", "F", "0" ) ] <- FALSE
  bad <- is.na(output) & !is.na(value) & nzchar(value)
  if (any(bad)) {
    stop( label, " contains invalid logical values: ", paste( unique(value[bad]), collapse = " | " ), call. = FALSE )
  }
  output
}
regex_escape <- function(x) {
  gsub( "([][{}()+*^$|\\\\?.])", "\\\\\\1", x )
}
file_size_safe <- function(path) {
  if (!file.exists(path)) {
    return(NA_real_)
  }
  as.numeric( file.info(path)$size )
}
elapsed_minutes <- function(start_time) {
  round( as.numeric( difftime( Sys.time(), start_time, units = "mins" ) ), 3 )
}
read_position_file <- function(path) {
  position_table <- data.table::fread( path, header = FALSE, data.table = TRUE, showProgress = FALSE )
  if (ncol(position_table) < 6L) {
    stop( "Position file has fewer than six columns: ", basename(path), call. = FALSE )
  }
  first_value <- tolower( trimws( as.character( position_table[[1L]][1L] ) ) )
  if (first_value %in% c( "barcode", "barcodes" )) {
    position_table <- data.table::fread( path, header = TRUE, data.table = TRUE, showProgress = FALSE )
    normalized_names <- tolower( gsub( "[^a-z0-9]+", "_", names(position_table) ) )
    setnames( position_table, normalized_names )
    barcode_candidates <- c( "barcode", "barcodes" )
    barcode_column <- barcode_candidates[ barcode_candidates %in% names(position_table) ][1L]
    if ( length(barcode_column) == 0L || is.na(barcode_column) ) {
      stop( "Position file header is present but no barcode column was found: ", basename(path), call. = FALSE )
    }
    requested_columns <- c( barcode_column, "in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres" )
    missing_columns <- setdiff( requested_columns, names(position_table) )
    if (length(missing_columns) > 0L) {
      stop( "Position file is missing column(s): ", paste( missing_columns, collapse = ", " ), call. = FALSE )
    }
    position_table <- position_table[ , ..requested_columns ]
    setnames( position_table, c( "raw_spot_barcode", "in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres" ) )
  } else {
    position_table <- position_table[ , 1:6 ]
    setnames( position_table, c( "raw_spot_barcode", "in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres" ) )
  }
  position_table[ , raw_spot_barcode := trimws( as.character( raw_spot_barcode ) ) ]
  position_table[ , in_tissue := suppressWarnings( as.integer( in_tissue ) ) ]
  position_table[ , array_row := suppressWarnings( as.numeric( array_row ) ) ]
  position_table[ , array_col := suppressWarnings( as.numeric( array_col ) ) ]
  position_table[ , pxl_row_in_fullres := suppressWarnings( as.numeric( pxl_row_in_fullres ) ) ]
  position_table[ , pxl_col_in_fullres := suppressWarnings( as.numeric( pxl_col_in_fullres ) ) ]
  position_table
}
to_dgCMatrix <- function( matrix_object, label) {
  if (inherits( matrix_object, "dgCMatrix" )) {
    return(matrix_object)
  }
  converted <- tryCatch( methods::as( matrix_object, "dgCMatrix" ), error = function(e) e )
  if (inherits( converted, "error" )) {
    intermediate <- tryCatch( methods::as( matrix_object, "CsparseMatrix" ), error = function(e) e )
    if (!inherits( intermediate, "error" )) {
      converted <- tryCatch( methods::as( intermediate, "dgCMatrix" ), error = function(e) e )
    }
  }
  if (inherits( converted, "error" )) {
    stop( "Could not convert ", label, " to dgCMatrix without dense conversion: ", conditionMessage( converted ), call. = FALSE )
  }
  converted
}
get_reference_counts <- function( reference_object) {
  if ( base::isS4( reference_object ) && "counts" %in% methods::slotNames( reference_object ) ) {
    return( methods::slot( reference_object, "counts" ) )
  }
  if ( is.list( reference_object ) && !is.null( reference_object$counts ) ) {
    return( reference_object$counts )
  }
  stop( "The saved Reference object does not expose a counts component.", call. = FALSE )
}
get_rctd_results <- function( rctd_object) {
  if ( base::isS4( rctd_object ) && "results" %in% methods::slotNames( rctd_object ) ) {
    return( methods::slot( rctd_object, "results" ) )
  }
  if ( is.list( rctd_object ) && !is.null( rctd_object$results ) ) {
    return( rctd_object$results )
  }
  stop( "The completed RCTD object has no accessible results component.", call. = FALSE )
}
orient_weight_matrix <- function( weights, spot_names, label) {
  if (inherits( weights, "Matrix" )) {
    weights <- as.matrix( weights )
  } else {
    weights <- as.matrix( weights )
  }
  if ( !is.null( rownames(weights) ) && setequal( rownames(weights), spot_names ) && nrow(weights) == length(spot_names) ) {
    weights <- weights[ spot_names, , drop = FALSE ]
    return(weights)
  }
  if ( !is.null( colnames(weights) ) && setequal( colnames(weights), spot_names ) && ncol(weights) == length(spot_names) ) {
    weights <- t( weights[ , spot_names, drop = FALSE ] )
    return(weights)
  }
  stop( label, " cannot be oriented to the expected ", length(spot_names), " spatial spots.", call. = FALSE )
}
manual_normalize_weights <- function( weights) {
  row_sums <- rowSums( weights )
  if (any( !is.finite(row_sums) | row_sums <= 0 )) {
    stop( "Cannot normalize weights because at least one row sum is nonpositive or non-finite.", call. = FALSE )
  }
  weights / row_sums
}
capture_success_files <- function( capture_area) {
  area_dir <- file.path( CAPTURE_ROOT, capture_area )
  list( area_dir = area_dir, success_flag = file.path( area_dir, "SUCCESS.flag" ), raw_weights_rds = file.path( area_dir, paste0( capture_area, "_RCTD_full_raw_weights.rds" ) ), normalized_weights_rds = file.path( area_dir, paste0( capture_area, "_RCTD_full_normalized_weights.rds" ) ), raw_weights_csv = file.path( area_dir, paste0( capture_area, "_RCTD_full_raw_weights.csv.gz" ) ), normalized_weights_csv = file.path( area_dir, paste0( capture_area, "_RCTD_full_normalized_weights.csv.gz" ) ), spot_metadata_csv = file.path( area_dir, paste0( capture_area, "_spot_metadata.csv.gz" ) ), audit_csv = file.path( area_dir, paste0( capture_area, "_RCTD_run_audit.csv" ) ), results_components_txt = file.path( area_dir, paste0( capture_area, "_RCTD_results_components.txt" ) ), warnings_txt = file.path( area_dir, paste0( capture_area, "_RCTD_warnings.txt" ) ), report_txt = file.path( area_dir, paste0( capture_area, "_RCTD_FINAL_REPORT.txt" ) ), full_rctd_rds = file.path( area_dir, paste0( capture_area, "_RCTD_full_object.rds" ) ) )
}
capture_is_complete <- function( capture_area) {
  files <- capture_success_files( capture_area )
  required <- c( files$success_flag, files$raw_weights_rds, files$normalized_weights_rds, files$spot_metadata_csv, files$audit_csv, files$report_txt )
  all( file.exists( required ) ) && all( file.info( required )$size > 0 )
}
write_capture_failure <- function( capture_area, reason, warning_lines, elapsed) {
  files <- capture_success_files( capture_area )
  dir.create( files$area_dir, recursive = TRUE, showWarnings = FALSE )
  failure_report <- c( "============================================================", paste0( "GSE234774 Step 12D RCTD | ", capture_area ), "============================================================", paste0( "Generated: ", format( Sys.time(), "%Y-%m-%d %H:%M:%S" ) ), "", "FINAL_STATUS: FAIL", paste0( "ERROR: ", reason ), paste0( "ELAPSED_MINUTES: ", elapsed ), "", "WARNINGS", if ( length(warning_lines) == 0L ) {
      "- None captured before failure."
    } else {
      paste0( "- ", warning_lines )
    },
    "", "BOUNDARIES", "- Other successful capture areas remain valid.", "- Rerun the same script; successful capture areas will be skipped.", "- No target-program score was calculated.", "============================================================"
  )
  writeLines( failure_report, files$report_txt, useBytes = TRUE )
  writeLines( if ( length(warning_lines) == 0L ) {
      "No warnings captured."
    } else {
      warning_lines
    },
    files$warnings_txt, useBytes = TRUE
  )
  if (file.exists( files$success_flag )) {
    unlink( files$success_flag, force = TRUE )
  }
}
checks <- data.table( check_id = character(), category = character(), description = character(), status = character(), observed = character(), required = character(), critical = logical(), action_if_not_pass = character() )
add_check <- function( check_id, category, description, passed, observed, required, critical = TRUE, action_if_not_pass = "") {
  checks <<- rbind( checks, data.table( check_id = as.character( check_id ), category = as.character( category ), description = as.character( description ), status = if ( isTRUE( passed ) ) {
        "PASS"
      } else {
        "FAIL"
      },
      observed = as.character( observed ), required = as.character( required ), critical = isTRUE( critical ), action_if_not_pass = as.character( action_if_not_pass )
    ), fill = TRUE
  )
}
add_review <- function( check_id, category, description, observed, required, action_if_not_pass) {
  checks <<- rbind( checks, data.table( check_id = as.character( check_id ), category = as.character( category ), description = as.character( description ), status = "REVIEW", observed = as.character( observed ), required = as.character( required ), critical = FALSE, action_if_not_pass = as.character( action_if_not_pass ) ), fill = TRUE )
}
write_fatal_global_report <- function( reason) {
  report <- c( "============================================================", "GSE234774 Step 12D full-mode RCTD", "============================================================", paste0( "Generated: ", format( Sys.time(), "%Y-%m-%d %H:%M:%S" ) ), paste0( "Output directory: ", OUTPUT_DIR ), "", "FINAL_STATUS: STOP", "", "FATAL PREFLIGHT ERROR", paste0( "- ", reason ), "", "BOUNDARIES", "- No target-program score was calculated.", "- Existing successful capture-area checkpoints were not deleted.", "", "NEXT ACTION", "- review this report and 00_run_log.txt.", "============================================================" )
  writeLines( report, OUT_FINAL_REPORT, useBytes = TRUE )
  log_message( "FATAL PREFLIGHT ERROR: ", reason )
  write_run_log()
}
run_one_capture_area <- function( capture_area, query_index, spot_qc, spatial_counts, reference_object, reference_genes, reference_gene_indices, SpatialRNA_constructor, create_RCTD_function, run_RCTD_function, normalize_weights_function) {
  area_start <- Sys.time()
  warning_lines <- character()
  files <- capture_success_files( capture_area )
  dir.create( files$area_dir, recursive = TRUE, showWarnings = FALSE )
  if ( capture_area %in% FORCE_RERUN_CAPTURE_AREAS ) {
    unlink( files$area_dir, recursive = TRUE, force = TRUE )
    dir.create( files$area_dir, recursive = TRUE, showWarnings = FALSE )
  }
  current_query_row <- query_index[ capture_area_id == capture_area ]
  if (nrow( current_query_row ) != 1L) {
    stop( "Expected one Step 12C query-index row for ", capture_area, "; observed ", nrow( current_query_row ), ".", call. = FALSE )
  }
  current_spots <- spot_qc[ capture_area_id == capture_area ]
  expected_spot_names <- as.character( current_spots$barcode )
  if ( length( expected_spot_names ) != current_query_row$n_spots ) {
    stop( "Spot count differs from the Step 12C query index for ", capture_area, ".", call. = FALSE )
  }
  spatial_column_indices <- match( expected_spot_names, colnames( spatial_counts ) )
  if (anyNA( spatial_column_indices )) {
    stop( "Some spatial barcodes are absent from the BPCells matrix for ", capture_area, ".", call. = FALSE )
  }
  log_message( capture_area, ": materializing ", length( reference_genes ), " genes x ", length( expected_spot_names ), " spots." )
  query_view <- spatial_counts[ reference_gene_indices, spatial_column_indices, drop = FALSE ]
  query_counts <- to_dgCMatrix( query_view, paste0( capture_area, " spatial query" ) )
  rm( query_view )
  invisible( gc() )
  if ( !identical( rownames( query_counts ), reference_genes ) ) {
    rownames( query_counts ) <- reference_genes
  }
  if ( !identical( colnames( query_counts ), expected_spot_names ) ) {
    colnames( query_counts ) <- expected_spot_names
  }
  query_nUMI <- as.numeric( Matrix::colSums( query_counts ) )
  names( query_nUMI ) <- colnames( query_counts )
  if (any( !is.finite( query_nUMI ) | query_nUMI <= 0 )) {
    stop( "At least one ", capture_area, " query spot has zero/non-finite nUMI in the shared gene space.", call. = FALSE )
  }
  position_file <- file.path( STEP12A2_ASSET_DIR, current_query_row$position_member )
  if (!file.exists( position_file )) {
    stop( "Position file is missing for ", capture_area, ": ", position_file, call. = FALSE )
  }
  position_table <- read_position_file( position_file )
  processed_prefix <- paste0( capture_area, "-" )
  if (!all( startsWith( expected_spot_names, processed_prefix ) )) {
    stop( "Processed barcode prefix mismatch for ", capture_area, ".", call. = FALSE )
  }
  raw_spot_barcodes <- sub( paste0( "^", regex_escape( processed_prefix ) ), "", expected_spot_names, perl = TRUE )
  position_match <- match( raw_spot_barcodes, position_table$raw_spot_barcode )
  if (anyNA( position_match )) {
    stop( "Position coordinates are missing for ", sum( is.na( position_match ) ), " spots in ", capture_area, ".", call. = FALSE )
  }
  matched_positions <- position_table[ position_match ]
  coordinate_table <- data.frame( x = matched_positions$pxl_col_in_fullres, y = matched_positions$pxl_row_in_fullres, row.names = expected_spot_names, check.names = FALSE )
  coordinate_ok <- nrow( coordinate_table ) == ncol( query_counts ) && identical( rownames( coordinate_table ), colnames( query_counts ) ) && all( is.finite( coordinate_table$x ) ) && all( is.finite( coordinate_table$y ) )
  if (!coordinate_ok) {
    stop( "Coordinate-count alignment failed for ", capture_area, ".", call. = FALSE )
  }
  current_spots[ , query_nUMI_shared_genes := query_nUMI[ barcode ] ]
  current_spots[ , coordinate_x_fullres := coordinate_table[ barcode, "x" ] ]
  current_spots[ , coordinate_y_fullres := coordinate_table[ barcode, "y" ] ]
  current_spots[ , position_in_tissue := matched_positions$in_tissue ]
  data.table::fwrite( current_spots, files$spot_metadata_csv, bom = TRUE, compress = "gzip" )
  query_object <- SpatialRNA_constructor( coordinate_table, query_counts, query_nUMI )
  log_message( capture_area, ": create.RCTD(max_cores=", MAX_CORES, ")." )
  rctd_object <- withCallingHandlers( create_RCTD_function( query_object, reference_object, max_cores = MAX_CORES ), warning = function(w) {
      warning_lines <<- c( warning_lines, paste0( "create.RCTD: ", conditionMessage( w ) ) )
      invokeRestart( "muffleWarning" )
    }
  )
  log_message( capture_area, ": run.RCTD(doublet_mode='full')." )
  rctd_object <- withCallingHandlers( run_RCTD_function( rctd_object, doublet_mode = "full" ), warning = function(w) {
      warning_lines <<- c( warning_lines, paste0( "run.RCTD: ", conditionMessage( w ) ) )
      invokeRestart( "muffleWarning" )
    }
  )
  rctd_results <- get_rctd_results( rctd_object )
  result_components <- names( rctd_results )
  if ( is.null( result_components ) ) {
    result_components <- character()
  }
  writeLines( c( paste0( "spacexr version: ", as.character( utils::packageVersion( "spacexr" ) ) ), paste0( "results components: ", paste( result_components, collapse = " | " ) ) ), files$results_components_txt, useBytes = TRUE )
  if ( is.null( rctd_results$weights ) ) {
    stop( "run.RCTD completed without a results$weights component for ", capture_area, ".", call. = FALSE )
  }
  raw_weights <- orient_weight_matrix( rctd_results$weights, expected_spot_names, paste0( capture_area, " raw RCTD weights" ) )
  if ( is.null( colnames( raw_weights ) ) || anyDuplicated( colnames( raw_weights ) ) ) {
    stop( "Raw RCTD weights have missing or duplicated cell-type columns for ", capture_area, ".", call. = FALSE )
  }
  normalization_method <- "spacexr::normalize_weights"
  normalized_weights <- tryCatch( normalize_weights_function( raw_weights ), error = function(e) {
      normalization_method <<- "manual_row_sum_fallback"
      manual_normalize_weights( raw_weights )
    }
  )
  normalized_weights <- orient_weight_matrix( normalized_weights, expected_spot_names, paste0( capture_area, " normalized RCTD weights" ) )
  normalized_weights <- normalized_weights[ , colnames( raw_weights ), drop = FALSE ]
  raw_finite <- all( is.finite( raw_weights ) )
  normalized_finite <- all( is.finite( normalized_weights ) )
  raw_nonnegative <- min( raw_weights, na.rm = TRUE ) >= -NEGATIVE_WEIGHT_TOLERANCE
  normalized_nonnegative <- min( normalized_weights, na.rm = TRUE ) >= -NEGATIVE_WEIGHT_TOLERANCE
  raw_row_sums <- rowSums( raw_weights )
  normalized_row_sums <- rowSums( normalized_weights )
  raw_positive_sums <- all( is.finite( raw_row_sums ) & raw_row_sums > 0 )
  normalized_sum_deviation <- max( abs( normalized_row_sums - 1 ), na.rm = TRUE )
  normalized_sums_ok <- is.finite( normalized_sum_deviation ) && normalized_sum_deviation <= WEIGHT_SUM_TOLERANCE
  weight_dimensions_ok <- nrow( raw_weights ) == length( expected_spot_names ) && nrow( normalized_weights ) == length( expected_spot_names ) && ncol( raw_weights ) == ncol( normalized_weights ) && identical( rownames( raw_weights ), expected_spot_names ) && identical( rownames( normalized_weights ), expected_spot_names )
  weight_audit_pass <- weight_dimensions_ok && raw_finite && normalized_finite && raw_nonnegative && normalized_nonnegative && raw_positive_sums && normalized_sums_ok
  if (!weight_audit_pass) {
    stop( "RCTD weight audit failed for ", capture_area, ": dimensions=", weight_dimensions_ok, ", raw_finite=", raw_finite, ", normalized_finite=", normalized_finite, ", raw_nonnegative=", raw_nonnegative, ", normalized_nonnegative=", normalized_nonnegative, ", raw_positive_sums=", raw_positive_sums, ", max_normalized_sum_deviation=", normalized_sum_deviation, ".", call. = FALSE )
  }
  saveRDS( raw_weights, files$raw_weights_rds, compress = FALSE )
  saveRDS( normalized_weights, files$normalized_weights_rds, compress = FALSE )
  raw_weight_table <- data.table( barcode = rownames( raw_weights ), raw_weights, check.names = FALSE )
  normalized_weight_table <- data.table( barcode = rownames( normalized_weights ), normalized_weights, check.names = FALSE )
  data.table::fwrite( raw_weight_table, files$raw_weights_csv, bom = TRUE, compress = "gzip" )
  data.table::fwrite( normalized_weight_table, files$normalized_weights_csv, bom = TRUE, compress = "gzip" )
  if ( !is.null( rctd_results$results_df ) ) {
    results_df <- as.data.table( rctd_results$results_df, keep.rownames = "barcode" )
    data.table::fwrite( results_df, file.path( files$area_dir, paste0( capture_area, "_RCTD_results_df.csv.gz" ) ), bom = TRUE, compress = "gzip" )
  }
  if (SAVE_FULL_RCTD_OBJECT) {
    saveRDS( rctd_object, files$full_rctd_rds, compress = FALSE )
  }
  warning_lines <- unique( warning_lines )
  writeLines( if ( length( warning_lines ) == 0L ) {
      "No warnings captured."
    } else {
      warning_lines
    },
    files$warnings_txt, useBytes = TRUE
  )
  area_elapsed <- elapsed_minutes( area_start )
  area_audit <- data.table( capture_area_id = capture_area, status = "PASS", elapsed_minutes = area_elapsed, n_spots = nrow( raw_weights ), n_reference_genes = nrow( query_counts ), n_cell_types = ncol( raw_weights ), cell_types = collapse_unique( colnames( raw_weights ) ), query_nonzero_entries = length( query_counts@x ), min_query_nUMI = min( query_nUMI ), median_query_nUMI = median( query_nUMI ), max_query_nUMI = max( query_nUMI ), result_components = collapse_unique( result_components ), normalization_method = normalization_method, raw_weight_min = min( raw_weights ), raw_weight_max = max( raw_weights ), raw_row_sum_min = min( raw_row_sums ), raw_row_sum_median = median( raw_row_sums ), raw_row_sum_max = max( raw_row_sums ), normalized_weight_min = min( normalized_weights ), normalized_weight_max = max( normalized_weights ), normalized_row_sum_min = min( normalized_row_sums ), normalized_row_sum_median = median( normalized_row_sums ), normalized_row_sum_max = max( normalized_row_sums ), max_normalized_row_sum_deviation = normalized_sum_deviation, warning_count = length( warning_lines ), full_RCTD_object_saved = SAVE_FULL_RCTD_OBJECT, raw_weights_rds = files$raw_weights_rds, normalized_weights_rds = files$normalized_weights_rds )
  data.table::fwrite( area_audit, files$audit_csv, bom = TRUE )
  area_report <- c( "============================================================", paste0( "GSE234774 Step 12D full-mode RCTD | ", capture_area ), "============================================================", paste0( "Generated: ", format( Sys.time(), "%Y-%m-%d %H:%M:%S" ) ), "", "FINAL_STATUS: PASS", paste0( "ELAPSED_MINUTES: ", area_elapsed ), "", "QUERY", paste0( "- Shared reference genes: ", nrow( query_counts ) ), paste0( "- Spatial spots: ", ncol( query_counts ) ), paste0( "- Query nonzero entries: ", length( query_counts@x ) ), paste0( "- Coordinate rows: ", nrow( coordinate_table ) ), "- All Step 12B spots were retained.", "", "RCTD", "- Mode: full", paste0( "- max_cores: ", MAX_CORES ), paste0( "- Cell types: ", ncol( raw_weights ), " [", paste( colnames( raw_weights ), collapse = ", " ), "]" ), paste0( "- Result components: ", paste( result_components, collapse = ", " ) ), "", "WEIGHT AUDIT", paste0( "- Raw weights finite/nonnegative: ", raw_finite && raw_nonnegative ), paste0( "- Raw row sums positive: ", raw_positive_sums ), paste0( "- Normalized weights finite/nonnegative: ", normalized_finite && normalized_nonnegative ), paste0( "- Maximum normalized row-sum deviation from 1: ", format( normalized_sum_deviation, scientific = TRUE ) ), paste0( "- Normalization method: ", normalization_method ), "", "WARNINGS", if ( length( warning_lines ) == 0L ) {
      "- None captured."
    } else {
      paste0( "- ", warning_lines )
    },
    "", "BOUNDARIES", "- No target-program score was calculated.", "- No spatial enrichment test was performed.", "- Full RCTD object saving was optional; weights are the frozen Step 12D outputs.", "============================================================"
  )
  writeLines( area_report, files$report_txt, useBytes = TRUE )
  writeLines( c( paste0( "capture_area_id=", capture_area ), paste0( "completed_at=", format( Sys.time(), "%Y-%m-%d %H:%M:%S" ) ), "status=PASS" ), files$success_flag, useBytes = TRUE )
  rm( query_counts, query_object, rctd_object, rctd_results, raw_weights, normalized_weights, raw_weight_table, normalized_weight_table )
  invisible( gc() )
  area_audit
}
main <- function() {
  workflow_start <- Sys.time()
  log_message( "============================================================" )
  log_message( "GSE234774 Step 12D full-mode RCTD started." )
  log_message( "============================================================" )
  required_inputs <- c( Step12C_report = STEP12C_REPORT, Reference_RDS = REFERENCE_RDS, Reference_metadata_RDS = REFERENCE_METADATA_RDS, Step12C_query_index = STEP12C_QUERY_INDEX, Step12B_spot_QC = STEP12B_SPOT_QC )
  file_status <- file.exists( required_inputs )
  cache_status <- dir.exists( SPATIAL_BPCELLS_DIR )
  asset_status <- dir.exists( STEP12A2_ASSET_DIR )
  add_check( "C01", "input", "All Step 12D frozen inputs, the spatial BPCells cache, and coordinate assets exist.", all( file_status ) && cache_status && asset_status, paste0( paste( paste0( names( required_inputs ), "=", file_status ), collapse = "; " ), "; spatial_BPCells=", cache_status, "; coordinate_assets=", asset_status ), "All TRUE", TRUE, "Restore frozen upstream outputs; do not silently substitute files." )
  if ( !all( file_status ) || !cache_status || !asset_status ) {
    stop( "At least one required Step 12D input is missing.", call. = FALSE )
  }
  step12C_text <- readLines( STEP12C_REPORT, warn = FALSE, encoding = "UTF-8" )
  step12C_status <- step12C_text[ grepl( "^FINAL_STATUS:", step12C_text ) ]
  step12C_pass <- length( step12C_status ) == 1L && grepl( "^FINAL_STATUS: (PASS|PASS_WITH_REVIEW)$", step12C_status )
  add_check( "C02", "upstream", "Step 12C releases the saved Reference and query manifests.", step12C_pass, if ( length( step12C_status ) > 0L ) {
      step12C_status[[1L]]
    } else {
      "No FINAL_STATUS line"
    },
    "PASS or PASS_WITH_REVIEW", TRUE, "Do not run RCTD until Step 12C passes."
  )
  if (!step12C_pass) {
    stop( "Step 12C does not release Step 12D.", call. = FALSE )
  }
  query_index <- data.table::fread( STEP12C_QUERY_INDEX, data.table = TRUE, showProgress = FALSE )
  spot_qc <- data.table::fread( STEP12B_SPOT_QC, data.table = TRUE, showProgress = FALSE )
  required_columns( query_index, c( "capture_area_id", "n_spots", "n_sections", "position_member", "common_reference_genes", "spot_retention_rule", "planned_RCTD_mode" ), "Step 12C query index" )
  required_columns( spot_qc, c( "barcode", "capture_area_id", "section_id", "timepoint", "replicate_index", "biological_replicate_id", "injury_condition", "nCount_Spatial_BPCells", "qc_review_flag" ), "Step 12B spot-QC table" )
  query_index[ , capture_area_id := toupper( trimws( as.character( capture_area_id ) ) ) ]
  spot_qc[ , capture_area_id := toupper( trimws( as.character( capture_area_id ) ) ) ]
  spot_qc[ , qc_review_flag := as_logical_strict( qc_review_flag, "qc_review_flag" ) ]
  query_design_ok <- nrow( query_index ) == 12L && sum( query_index$n_spots ) == 33941L && all( query_index$n_sections == 3L ) && all( query_index$spot_retention_rule == "PRIMARY_RETAIN_ALL_33941_SPOTS" ) && all( query_index$planned_RCTD_mode == "full" ) && nrow( spot_qc ) == 33941L && !anyDuplicated( spot_qc$barcode )
  add_check( "C03", "query_design", "The frozen design contains 12 separate full-mode queries and all 33,941 spots.", query_design_ok, paste0( "queries=", nrow( query_index ), "; query spots=", sum( query_index$n_spots ), "; spot-QC rows=", nrow( spot_qc ), "; retained review spots=", sum( spot_qc$qc_review_flag ) ), "12 queries, 33,941 spots, full mode, all spots retained", TRUE, "Do not merge areas or remove review spots before RCTD." )
  if (!query_design_ok) {
    stop( "The frozen Step 12D query design was not reproduced.", call. = FALSE )
  }
  reference_object <- readRDS( REFERENCE_RDS )
  reference_metadata <- readRDS( REFERENCE_METADATA_RDS )
  reference_counts <- get_reference_counts( reference_object )
  reference_genes <- rownames( reference_counts )
  reference_ok <- nrow( reference_counts ) == 24321L && ncol( reference_counts ) == 19593L && length( reference_genes ) == 24321L && !anyDuplicated( reference_genes ) && nrow( reference_metadata ) == 19593L
  add_check( "C04", "reference", "The frozen Step 12C Reference reloads with expected dimensions.", reference_ok, paste0( "counts=", nrow( reference_counts ), " x ", ncol( reference_counts ), "; metadata rows=", nrow( reference_metadata ), "; class=", paste( class( reference_object ), collapse = " | " ) ), "24,321 genes x 19,593 nuclei", TRUE, "Do not rebuild or alter the frozen Step 12C Reference." )
  if (!reference_ok) {
    stop( "The saved Step 12C Reference does not match the frozen design.", call. = FALSE )
  }
  spatial_counts <- BPCells::open_matrix_dir( SPATIAL_BPCELLS_DIR )
  spatial_gene_indices <- match( reference_genes, rownames( spatial_counts ) )
  reference_spatial_gene_alignment <- !anyNA( spatial_gene_indices ) && length( unique( spatial_gene_indices ) ) == length( reference_genes )
  add_check( "C05", "gene_alignment", "Every frozen Reference gene is present exactly once in the spatial BPCells matrix.", reference_spatial_gene_alignment, paste0( "matched=", sum( !is.na( spatial_gene_indices ) ), "/", length( reference_genes ), "; unique indices=", length( unique( spatial_gene_indices[ !is.na( spatial_gene_indices ) ] ) ) ), "24,321/24,321 exact unique matches", TRUE, "Do not run RCTD with a changed feature space." )
  if (!reference_spatial_gene_alignment) {
    stop( "Reference-spatial gene alignment failed.", call. = FALSE )
  }
  spacexr_namespace <- asNamespace( "spacexr" )
  required_spacexr_functions <- c( "SpatialRNA", "create.RCTD", "run.RCTD" )
  function_status <- vapply( required_spacexr_functions, exists, logical(1), envir = spacexr_namespace, inherits = FALSE )
  normalize_available <- exists( "normalize_weights", envir = spacexr_namespace, inherits = FALSE )
  add_check( "C06", "software", "Installed spacexr provides the required RCTD constructors and run function.", all( function_status ), paste( paste0( names( function_status ), "=", function_status ), collapse = "; " ), "SpatialRNA, create.RCTD, and run.RCTD all TRUE", TRUE, "Use the installed spacexr 2.2.1 API; do not mix incompatible wrappers." )
  if (!all( function_status )) {
    stop( "The installed spacexr namespace lacks a required function.", call. = FALSE )
  }
  if (!normalize_available) {
    add_review( "R01", "software", "spacexr::normalize_weights is unavailable.", "Manual row-sum normalization will be used.", "Normalized rows must still sum to 1 within tolerance.", "Record the fallback in every capture-area audit." )
  }
  SpatialRNA_constructor <- get( "SpatialRNA", envir = spacexr_namespace )
  create_RCTD_function <- get( "create.RCTD", envir = spacexr_namespace )
  run_RCTD_function <- get( "run.RCTD", envir = spacexr_namespace )
  normalize_weights_function <- if ( normalize_available ) {
    get( "normalize_weights", envir = spacexr_namespace )
  } else {
    manual_normalize_weights
  }
  detected_cores <- parallel::detectCores( logical = TRUE )
  if ( !is.finite( detected_cores ) || detected_cores < 1L ) {
    detected_cores <- 1L
  }
  MAX_CORES <<- max( 1L, min( as.integer( MAX_CORES ), as.integer( detected_cores ) ) )
  all_capture_areas <- as.character( query_index$capture_area_id )
  capture_areas_to_run <- if ( length( CAPTURE_AREAS_TO_RUN ) == 1L && identical( toupper( CAPTURE_AREAS_TO_RUN ), "ALL" ) ) {
    all_capture_areas
  } else {
    toupper( as.character( CAPTURE_AREAS_TO_RUN ) )
  }
  unknown_requested <- setdiff( capture_areas_to_run, all_capture_areas )
  if (length( unknown_requested ) > 0L) {
    stop( "Unknown requested capture area(s): ", paste( unknown_requested, collapse = " | " ), call. = FALSE )
  }
  parameter_audit <- rbindlist( list( data.table( item_type = "package", item = required_packages, value = vapply( required_packages, function(package_name) {
              as.character( utils::packageVersion( package_name ) )
            },
            FUN.VALUE = character(1)
          ), note = "Version used in Step 12D"
      ), data.table(
        item_type = "parameter", item = c( "CAPTURE_AREAS_TO_RUN", "MAX_CORES", "detected_logical_cores", "RESUME_COMPLETED", "FORCE_RERUN_CAPTURE_AREAS", "SAVE_FULL_RCTD_OBJECT", "WEIGHT_SUM_TOLERANCE", "NEGATIVE_WEIGHT_TOLERANCE", "RCTD_mode", "spatial_spot_rule", "reference_dimensions" ), value = c( paste( capture_areas_to_run, collapse = " | " ), MAX_CORES, detected_cores, RESUME_COMPLETED, collapse_unique( FORCE_RERUN_CAPTURE_AREAS ), SAVE_FULL_RCTD_OBJECT, WEIGHT_SUM_TOLERANCE, NEGATIVE_WEIGHT_TOLERANCE, "full", "retain_all_33941_spots", paste0( nrow( reference_counts ), " genes x ", ncol( reference_counts ), " nuclei" ) ), note = c( "Sequential capture-area run order", "Internal spacexr cores per area", "Windows logical cores", "Skip SUCCESS.flag checkpoints", "Intentional replacements only", "Weights are saved regardless", "Normalized row-sum audit", "Allow tiny floating-point negatives only", "Visium broad-mixture mode", "Includes five QC-review spots", "Frozen Step 12C Reference" )
      )
    ), fill = TRUE
  )
  data.table::fwrite( parameter_audit, OUT_INPUT_AUDIT, bom = TRUE )
  run_summary_list <- vector( mode = "list", length = length( capture_areas_to_run ) )
  for ( capture_index in seq_along( capture_areas_to_run ) ) {
    capture_area <- capture_areas_to_run[ capture_index ]
    log_message( "[", capture_index, "/", length( capture_areas_to_run ), "] ", capture_area )
    if ( RESUME_COMPLETED && !capture_area %in% FORCE_RERUN_CAPTURE_AREAS && capture_is_complete( capture_area ) ) {
      log_message( capture_area, ": existing successful checkpoint found; skipping." )
      files <- capture_success_files( capture_area )
      prior_audit <- data.table::fread( files$audit_csv, data.table = TRUE, showProgress = FALSE )
      prior_audit[ , run_action := "SKIPPED_EXISTING_SUCCESS" ]
      run_summary_list[[ capture_index ]] <- prior_audit
      next
    }
    area_warnings <- character()
    area_start <- Sys.time()
    area_result <- tryCatch( run_one_capture_area( capture_area = capture_area, query_index = query_index, spot_qc = spot_qc, spatial_counts = spatial_counts, reference_object = reference_object, reference_genes = reference_genes, reference_gene_indices = spatial_gene_indices, SpatialRNA_constructor = SpatialRNA_constructor, create_RCTD_function = create_RCTD_function, run_RCTD_function = run_RCTD_function, normalize_weights_function = normalize_weights_function ), error = function(e) {
        write_capture_failure( capture_area = capture_area, reason = conditionMessage( e ), warning_lines = area_warnings, elapsed = elapsed_minutes( area_start ) )
        data.table( capture_area_id = capture_area, status = "FAIL", elapsed_minutes = elapsed_minutes( area_start ), error_message = conditionMessage( e ), run_action = "ATTEMPTED_FAILED" )
      }
    )
    if ( "run_action" %in% names( area_result ) == FALSE ) {
      area_result[ , run_action := "ATTEMPTED_PASS" ]
    }
    run_summary_list[[ capture_index ]] <- area_result
    data.table::fwrite( rbindlist( run_summary_list[ !vapply( run_summary_list, is.null, logical(1) ) ], fill = TRUE ), OUT_RUN_SUMMARY, bom = TRUE )
    invisible( gc() )
  }
  run_summary <- rbindlist( run_summary_list[ !vapply( run_summary_list, is.null, logical(1) ) ], fill = TRUE )
  all_completed <- all_capture_areas[ vapply( all_capture_areas, capture_is_complete, logical(1) ) ]
  unrepresented_completed <- setdiff( all_completed, run_summary$capture_area_id )
  if (length( unrepresented_completed ) > 0L) {
    extra_rows <- rbindlist( lapply( unrepresented_completed, function(capture_area) {
          files <- capture_success_files( capture_area )
          row <- data.table::fread( files$audit_csv, data.table = TRUE, showProgress = FALSE )
          row[ , run_action := "INCLUDED_EXISTING_SUCCESS" ]
          row
        }
      ), fill = TRUE
    )
    run_summary <- rbind( run_summary, extra_rows, fill = TRUE )
  }
  run_summary[ , step12D_capture_order__ := match( capture_area_id, all_capture_areas ) ]
  data.table::setorder( run_summary, step12D_capture_order__ )
  run_summary[ , step12D_capture_order__ := NULL ]
  data.table::fwrite( run_summary, OUT_RUN_SUMMARY, bom = TRUE )
  completed_capture_areas <- all_capture_areas[ vapply( all_capture_areas, capture_is_complete, logical(1) ) ]
  failed_capture_areas <- setdiff( all_capture_areas, completed_capture_areas )
  all_areas_complete <- length( failed_capture_areas ) == 0L
  add_check( "C07", "RCTD_completion", "All 12 capture areas have successful checkpoint outputs.", all_areas_complete, paste0( "complete=", length( completed_capture_areas ), "/12; incomplete=", paste( failed_capture_areas, collapse = " | " ) ), "12/12 complete", TRUE, "Rerun the script; completed areas will be skipped automatically." )
  combined_raw <- NULL
  combined_normalized <- NULL
  combined_metadata <- NULL
  section_summary <- NULL
  common_cell_types_ok <- FALSE
  all_weights_valid <- FALSE
  if (length( completed_capture_areas ) > 0L) {
    raw_list <- list()
    normalized_list <- list()
    metadata_list <- list()
    cell_type_sets <- list()
    for ( capture_area in completed_capture_areas ) {
      files <- capture_success_files( capture_area )
      raw_matrix <- readRDS( files$raw_weights_rds )
      normalized_matrix <- readRDS( files$normalized_weights_rds )
      current_metadata <- data.table::fread( files$spot_metadata_csv, data.table = TRUE, showProgress = FALSE )
      raw_matrix <- orient_weight_matrix( raw_matrix, as.character( current_metadata$barcode ), paste0( capture_area, " saved raw weights" ) )
      normalized_matrix <- orient_weight_matrix( normalized_matrix, as.character( current_metadata$barcode ), paste0( capture_area, " saved normalized weights" ) )
      raw_list[[ capture_area ]] <- raw_matrix
      normalized_list[[ capture_area ]] <- normalized_matrix
      metadata_list[[ capture_area ]] <- current_metadata
      cell_type_sets[[ capture_area ]] <- colnames( normalized_matrix )
    }
    first_cell_types <- cell_type_sets[[ 1L ]]
    common_cell_types_ok <- all( vapply( cell_type_sets, function(x) {
          setequal( x, first_cell_types )
        },
        logical(1)
      )
    )
    if (common_cell_types_ok) {
      raw_list <- lapply( raw_list, function(x) {
          x[ , first_cell_types, drop = FALSE ]
        }
      )
      normalized_list <- lapply( normalized_list, function(x) {
          x[ , first_cell_types, drop = FALSE ]
        }
      )
      combined_raw <- do.call( rbind, raw_list )
      combined_normalized <- do.call( rbind, normalized_list )
      combined_metadata <- rbindlist( metadata_list, fill = TRUE )
      combined_metadata <- combined_metadata[ match( rownames( combined_normalized ), barcode ) ]
      all_weights_valid <- nrow( combined_normalized ) == sum( query_index[ capture_area_id %in% completed_capture_areas, n_spots ] ) && !anyDuplicated( rownames( combined_normalized ) ) && all( is.finite( combined_normalized ) ) && min( combined_normalized, na.rm = TRUE ) >= -NEGATIVE_WEIGHT_TOLERANCE && max( abs( rowSums( combined_normalized ) - 1 ), na.rm = TRUE ) <= WEIGHT_SUM_TOLERANCE && identical( combined_metadata$barcode, rownames( combined_normalized ) )
    }
  }
  add_check( "C08", "weight_schema", "All completed capture areas have the same RCTD cell-type columns.", if ( length( completed_capture_areas ) > 0L ) {
      common_cell_types_ok
    } else {
      FALSE
    },
    if ( length( completed_capture_areas ) > 0L ) {
      paste0( "completed areas=", length( completed_capture_areas ), "; shared cell types=", if ( common_cell_types_ok ) {
          paste( colnames( combined_normalized ), collapse = " | " )
        } else {
          "MISMATCH"
        }
      )
    } else {
      "No completed areas"
    },
    "Identical broad-cell-type set", TRUE, "Do not combine or compare weights with inconsistent cell-type columns."
  )
  add_check( "C09", "weight_audit", "Combined normalized weights are finite, nonnegative, uniquely indexed, and sum to 1.", if ( length( completed_capture_areas ) > 0L ) {
      all_weights_valid
    } else {
      FALSE
    },
    if ( length( completed_capture_areas ) > 0L && common_cell_types_ok ) {
      paste0( "spots=", nrow( combined_normalized ), "; duplicated barcodes=", anyDuplicated( rownames( combined_normalized ) ), "; min weight=", min( combined_normalized ), "; max row-sum deviation=", max( abs( rowSums( combined_normalized ) - 1 ) ) )
    } else {
      "Combined matrix unavailable"
    },
    "All completed spots pass", TRUE, "Resolve failed or inconsistent capture-area outputs before downstream scoring."
  )
  if ( length( completed_capture_areas ) > 0L && common_cell_types_ok && all_weights_valid ) {
    saveRDS( combined_raw, OUT_RAW_WEIGHTS_RDS, compress = FALSE )
    saveRDS( combined_normalized, OUT_NORMALIZED_WEIGHTS_RDS, compress = FALSE )
    raw_table <- data.table( barcode = rownames( combined_raw ), combined_raw, check.names = FALSE )
    normalized_table <- data.table( barcode = rownames( combined_normalized ), combined_normalized, check.names = FALSE )
    raw_table <- cbind( combined_metadata[ , .( barcode, capture_area_id, section_id, timepoint, replicate_index, biological_replicate_id, injury_condition, qc_review_flag ) ], raw_table[ , !"barcode" ] )
    normalized_table <- cbind( combined_metadata[ , .( barcode, capture_area_id, section_id, timepoint, replicate_index, biological_replicate_id, injury_condition, qc_review_flag ) ], normalized_table[ , !"barcode" ] )
    data.table::fwrite( raw_table, OUT_RAW_WEIGHTS, bom = TRUE, compress = "gzip" )
    data.table::fwrite( normalized_table, OUT_NORMALIZED_WEIGHTS, bom = TRUE, compress = "gzip" )
    weight_long <- data.table::melt( normalized_table, id.vars = c( "barcode", "capture_area_id", "section_id", "timepoint", "replicate_index", "biological_replicate_id", "injury_condition", "qc_review_flag" ), variable.name = "rctd_cell_type", value.name = "normalized_weight" )
    section_summary <- weight_long[ , .( n_spots = .N, mean_weight = mean( normalized_weight ), median_weight = median( normalized_weight ), sd_weight = stats::sd( normalized_weight ), q10_weight = as.numeric( quantile( normalized_weight, 0.10, names = FALSE ) ), q90_weight = as.numeric( quantile( normalized_weight, 0.90, names = FALSE ) ), min_weight = min( normalized_weight ), max_weight = max( normalized_weight ) ), by = .( capture_area_id, section_id, timepoint, replicate_index, biological_replicate_id, injury_condition, rctd_cell_type ) ]
    data.table::fwrite( section_summary, OUT_SECTION_SUMMARY, bom = TRUE )
    handoff <- list( step = "Step12D", completed_at = as.character( Sys.time() ), final_all_areas_complete = all_areas_complete, completed_capture_areas = completed_capture_areas, incomplete_capture_areas = failed_capture_areas, normalized_weights_rds = OUT_NORMALIZED_WEIGHTS_RDS, raw_weights_rds = OUT_RAW_WEIGHTS_RDS, normalized_weights_csv = OUT_NORMALIZED_WEIGHTS, section_summary_csv = OUT_SECTION_SUMMARY, cell_types = colnames( combined_normalized ), n_spots = nrow( combined_normalized ), reference_rds = REFERENCE_RDS, spatial_BPCells_path = SPATIAL_BPCELLS_DIR, frozen_spot_rule = "retain_all_33941_spots", target_program_scored = FALSE, next_step = paste( "Step 12E computes the target-program score and adjusts it", "for Microglia RCTD weight. Do not treat the target program", "as an RCTD cell type." ) )
    saveRDS( handoff, OUT_HANDOFF_RDS, compress = FALSE )
  }
  critical_failures <- checks[ critical == TRUE & status == "FAIL" ]
  review_items <- checks[ status == "REVIEW" ]
  final_status <- if ( all_areas_complete && nrow( critical_failures ) == 0L ) {
    if ( nrow( review_items ) > 0L ) {
      "PASS_WITH_REVIEW"
    } else {
      "PASS"
    }
  } else if (
    length( completed_capture_areas ) > 0L
  ) {
    "PARTIAL_CHECKPOINT"
  } else {
    "STOP"
  }
  data.table::fwrite( checks, OUT_CHECKS, bom = TRUE )
  total_elapsed <- elapsed_minutes( workflow_start )
  report <- c( "============================================================", "GSE234774 Step 12D full-mode RCTD per capture area", "============================================================", paste0( "Generated: ", format( Sys.time(), "%Y-%m-%d %H:%M:%S" ) ), paste0( "Output directory: ", OUTPUT_DIR ), "", paste0( "FINAL_STATUS: ", final_status ), "", "DESIGN", "- Twelve capture areas are processed separately.", "- RCTD mode: full.", paste0( "- spacexr max_cores per capture area: ", MAX_CORES, "." ), "- All 33,941 Step 12B spots, including five review spots, remain eligible.", "- Target inflammatory-remodeling program is not an RCTD cell type.", "", "COMPLETION", paste0( "- Successful capture areas: ", length( completed_capture_areas ), "/12 [", paste( completed_capture_areas, collapse = ", " ), "]." ), paste0( "- Incomplete capture areas: ", length( failed_capture_areas ), " [", paste( failed_capture_areas, collapse = ", " ), "]." ), paste0( "- Successful spots represented: ", if ( !is.null( combined_normalized ) ) {
        nrow( combined_normalized )
      } else {
        0L
      },
      "/33,941."
    ), "", "WEIGHT OUTPUT", if (
      !is.null( combined_normalized ) && common_cell_types_ok
    ) {
      paste0( "- Broad cell types: ", ncol( combined_normalized ), " [", paste( colnames( combined_normalized ), collapse = ", " ), "]." )
    } else {
      "- Combined broad-cell-type schema is not yet available."
    },
    if ( !is.null( combined_normalized ) && all_weights_valid ) {
      paste0( "- Maximum normalized row-sum deviation from 1: ", format( max( abs( rowSums( combined_normalized ) - 1 ) ), scientific = TRUE ), "." )
    } else {
      "- Combined normalized weight audit is incomplete."
    },
    paste0( "- Full RCTD objects saved: ", SAVE_FULL_RCTD_OBJECT, "." ), "", "CRITICAL CHECKS", paste0( "- PASS: ", nrow( checks[ critical == TRUE & status == "PASS" ] ), "; FAIL: ", nrow( critical_failures ), "." )
  )
  if (nrow( critical_failures ) > 0L) {
    report <- c( report, "", "FAILED CHECKS", paste0( "- ", critical_failures$check_id, ": ", critical_failures$description, " Observed: ", critical_failures$observed, " Action: ", critical_failures$action_if_not_pass ) )
  }
  if (nrow( review_items ) > 0L) {
    report <- c( report, "", "REVIEW ITEMS", paste0( "- ", review_items$check_id, ": ", review_items$description, " Observed: ", review_items$observed, " Action: ", review_items$action_if_not_pass ) )
  }
  report <- c( report, "", "INTERPRETATION", if ( final_status %in% c( "PASS", "PASS_WITH_REVIEW" ) ) {
      paste( "- Step 12D is technically complete.", "All capture areas have full-mode RCTD weights with valid", "normalized row sums." )
    } else if (
      final_status == "PARTIAL_CHECKPOINT"
    ) {
      paste( "- Step 12D is partially complete.", "Existing successful capture areas are preserved.", "Rerun this script to continue the remaining areas." )
    } else {
      "- Step 12D is not released."
    },
    "- These weights estimate broad cell-type composition only.", "- No target-program score or lesion-zone test was performed.", "", "NEXT ACTION", if ( final_status %in% c( "PASS", "PASS_WITH_REVIEW" ) ) {
      paste( "- Step 12E: compute the frozen target-program spatial score,", "then adjust/stratify it by the Microglia RCTD weight.", "Do not redefine the program as a cell type." )
    } else {
      "- Rerun Step 12D; successful capture-area checkpoints will be skipped."
    },
    "", "KEY OUTPUTS:", "- tables/01_step12D_input_and_parameter_audit.csv", "- tables/02_capture_area_RCTD_run_summary.csv", "- tables/03_all_capture_areas_normalized_weights.csv.gz (when all areas complete)", "- tables/04_all_capture_areas_raw_weights.csv.gz (when all areas complete)", "- tables/05_cell_type_weight_summary_by_section.csv (when all areas complete)", "- tables/09_STEP12D_audit_checks.csv", "- 10_STEP12D_FINAL_REPORT.txt", "", "LARGE INTERMEDIATE FILES", paste0( "- ", OUT_NORMALIZED_WEIGHTS_RDS ), paste0( "- ", OUT_RAW_WEIGHTS_RDS ), paste0( "- ", OUT_HANDOFF_RDS ), "- capture_areas/* full RCTD objects, if optional saving was enabled.", "", paste0( "TOTAL ELAPSED TIME THIS INVOCATION: ", total_elapsed, " minutes" ), "============================================================"
  )
  writeLines( report, OUT_FINAL_REPORT, useBytes = TRUE )
  writeLines( capture.output( sessionInfo() ), SESSION_INFO_FILE, useBytes = TRUE )
  log_message( "FINAL_STATUS: ", final_status )
  log_message( "Successful capture areas: ", length( completed_capture_areas ), "/12" )
  log_message( "Incomplete capture areas: ", paste( failed_capture_areas, collapse = " | " ) )
  log_message( "Open first: ", OUT_FINAL_REPORT )
  write_run_log()
  cat("\n")
  cat( "============================================================\n" )
  cat( "Step 12D completed: ", final_status, "\n", sep = "" )
  cat( "Successful capture areas: ", length( completed_capture_areas ), "/12\n", sep = "" )
  cat( "Open first: ", OUT_FINAL_REPORT, "\n", sep = "" )
  cat( "============================================================\n" )
  invisible( list( final_status = final_status, completed_capture_areas = completed_capture_areas, incomplete_capture_areas = failed_capture_areas, run_summary = run_summary, checks = checks, normalized_weights_rds = if ( file.exists( OUT_NORMALIZED_WEIGHTS_RDS ) ) {
          OUT_NORMALIZED_WEIGHTS_RDS
        } else {
          NA_character_
        }
    )
  )
}
tryCatch( main(), error = function(e) {
    write_fatal_global_report( conditionMessage( e ) )
    message( "Step 12D stopped during preflight/global processing: ", conditionMessage( e ) )
  }
)
