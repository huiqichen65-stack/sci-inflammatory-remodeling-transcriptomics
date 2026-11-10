options( stringsAsFactors = FALSE, warn = 1, width = 180, future.globals.maxSize = 24 * 1024^3 )
PROJECT_ROOT <- Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
CAPTURE_AREAS_TO_RUN <- "ALL"
RESUME_COMPLETED <- TRUE
FORCE_RERUN_CAPTURE_AREAS <- character()
AUC_MAX_RANK_FRACTION <- 0.05
MIN_TARGET_MAPPING_FRACTION <- 0.90
MIN_SHARED_MAPPING_FRACTION <- 0.90
MIN_TARGET_GENES_PRESENT <- 100L
MIN_SHARED_GENES_PRESENT <- 30L
LOGCP10K_SCALE_FACTOR <- 10000
AUCELL_NCORES <- 1L
POINT_CEX <- 0.42
PLOT_QUANTILE_LOW <- 0.01
PLOT_QUANTILE_HIGH <- 0.99
EXPECTED_SPATIAL_GENES <- 25966L
EXPECTED_SPATIAL_SPOTS <- 33941L
EXPECTED_CAPTURE_AREAS <- 12L
EXPECTED_SECTIONS <- 36L
EXPECTED_TARGET_STEP11B_GENES <- 496L
EXPECTED_SHARED_STEP11B_GENES <- 111L
STEP11B_DIR <- file.path( PROJECT_ROOT, "04_results", "step11B_GSE234774_refined_main_Arm_specific_shared_AUCell_ssGSEA" )
STEP11B_REPORT <- file.path( STEP11B_DIR, "10_STEP11B_FINAL_REPORT.txt" )
STEP11B_RESULTS_RDS <- file.path( STEP11B_DIR, "RDS", "GSE234774_step11B_refined_Arm_AUCell_ssGSEA_results.rds" )
STEP12D_DIR <- file.path( PROJECT_ROOT, "04_results", "step12D_GSE234774_full_mode_RCTD_per_capture_area" )
STEP12D_REPORT <- file.path( STEP12D_DIR, "10_STEP12D_FINAL_REPORT.txt" )
STEP12D_NORMALIZED_WEIGHTS_RDS <- file.path( STEP12D_DIR, "RDS", "GSE234774_step12D_all_spots_normalized_RCTD_weights.rds" )
STEP12D_CAPTURE_ROOT <- file.path( STEP12D_DIR, "capture_areas" )
STEP12B_DIR <- file.path( PROJECT_ROOT, "04_results", "step12B_fast_BPCells_spatial_section_QC" )
STEP12B_REPORT <- file.path( STEP12B_DIR, "10_STEP12B_FAST_FINAL_REPORT.txt" )
STEP12B_SPOT_QC <- file.path( STEP12B_DIR, "tables", "04_spot_QC_flags.csv.gz" )
STEP12B_CAPTURE_AUDIT <- file.path( STEP12B_DIR, "tables", "02_capture_area_image_coordinate_audit.csv" )
SPATIAL_BPCELLS_DIR <- file.path( STEP12B_DIR, "BPCells", "GSE234774_spatial_2d_full_counts" )
LOWRES_PREVIEW_ROOT <- file.path( STEP12B_DIR, "lowres_asset_preview" )
OUTPUT_DIR <- file.path( PROJECT_ROOT, "04_results", "step12E_GSE234774_target_program_spatial_scoring_microglia_adjustment" )
CAPTURE_SCORE_ROOT <- file.path( OUTPUT_DIR, "capture_area_scores" )
TABLE_DIR <- file.path( OUTPUT_DIR, "tables" )
FIGURE_DIR <- file.path( OUTPUT_DIR, "figures" )
RDS_DIR <- file.path( OUTPUT_DIR, "RDS" )
RUN_LOG_FILE <- file.path( OUTPUT_DIR, "00_run_log.txt" )
SESSION_INFO_FILE <- file.path( OUTPUT_DIR, "00_sessionInfo.txt" )
OUT_INPUT_GENE_AUDIT <- file.path( TABLE_DIR, "01_input_gene_set_mapping_audit.csv" )
OUT_SPOT_SCORES <- file.path( TABLE_DIR, "02_spot_target_shared_microglia_scores.csv.gz" )
OUT_SECTION_SUMMARY <- file.path( TABLE_DIR, "03_section_target_shared_score_summary.csv" )
OUT_ADJUSTMENT_AUDIT <- file.path( TABLE_DIR, "04_microglia_adjustment_audit.csv" )
OUT_SPATIAL_MAPS <- file.path( FIGURE_DIR, "01_target_microglia_adjusted_shared_spatial_maps.pdf" )
OUT_CHECKS <- file.path( TABLE_DIR, "09_STEP12E_audit_checks.csv" )
OUT_FINAL_REPORT <- file.path( OUTPUT_DIR, "10_STEP12E_FINAL_REPORT.txt" )
OUT_SPOT_SCORES_RDS <- file.path( RDS_DIR, "GSE234774_step12E_spot_target_shared_microglia_scores.rds" )
OUT_HANDOFF_RDS <- file.path( RDS_DIR, "GSE234774_step12E_target_program_spatial_score_handoff.rds" )
for (directory in c( OUTPUT_DIR, CAPTURE_SCORE_ROOT, TABLE_DIR, FIGURE_DIR, RDS_DIR )) {
  dir.create( directory, recursive = TRUE, showWarnings = FALSE )
}
required_packages <- c( "data.table", "Matrix", "BPCells", "AUCell", "png" )
missing_packages <- required_packages[ !vapply( required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1) ) ]
if (length(missing_packages) > 0L) {
  stop( "Missing required package(s): ", paste( missing_packages, collapse = ", " ), ". Install the missing package(s) before rerunning Step 12E.", call. = FALSE )
}
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(BPCells)
  library(AUCell)
  library(png)
})
run_log <- character()
write_run_log <- function() {
  writeLines( run_log, RUN_LOG_FILE, useBytes = TRUE )
}
elapsed_minutes <- function(start_time) {
  round( as.numeric( difftime( Sys.time(), start_time, units = "mins" ) ), 3 )
}
required_columns <- function( table, columns, label) {
  missing <- setdiff( columns, names(table) )
  if (length(missing) > 0L) {
    stop( label, " is missing required column(s): ", paste( missing, collapse = ", " ), call. = FALSE )
  }
  invisible(TRUE)
}
collapse_unique <- function( x, max_values = 50L) {
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
    stop( label, " contains invalid logical value(s): ", paste( unique(value[bad]), collapse = " | " ), call. = FALSE )
  }
  output
}
clean_gene_vector <- function(x) {
  output <- trimws( as.character(x) )
  output <- output[ !is.na(output) & nzchar(output) ]
  unique(output)
}
resolve_features <- function( requested_features, available_features) {
  available_upper <- toupper( available_features )
  first_index <- !duplicated( available_upper )
  lookup <- setNames( available_features[first_index], available_upper[first_index] )
  resolved <- unname( lookup[ toupper( requested_features ) ] )
  names(resolved) <- requested_features
  resolved
}
safe_zscore <- function(x) {
  output <- rep( NA_real_, length(x) )
  valid <- is.finite(x)
  if (sum(valid) < 2L) {
    output[valid] <- 0
    return(output)
  }
  value_sd <- stats::sd( x[valid] )
  if ( !is.finite(value_sd) || value_sd == 0 ) {
    output[valid] <- 0
    return(output)
  }
  output[valid] <- ( x[valid] - mean( x[valid] ) ) / value_sd
  output
}
safe_spearman <- function( x, y) {
  valid <- is.finite(x) & is.finite(y)
  if ( sum(valid) < 4L || stats::sd(x[valid]) == 0 || stats::sd(y[valid]) == 0 ) {
    return(NA_real_)
  }
  suppressWarnings( stats::cor( x[valid], y[valid], method = "spearman" ) )
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
    stop( "Could not convert ", label, " to dgCMatrix without dense conversion: ", conditionMessage(converted), call. = FALSE )
  }
  converted
}
file_size_safe <- function(path) {
  if (!file.exists(path)) {
    return(NA_real_)
  }
  as.numeric( file.info(path)$size )
}
build_rankings_version_safe <- function( expression_matrix) {
  function_object <- AUCell::AUCell_buildRankings
  formal_names <- names( formals( function_object ) )
  arguments <- list( exprMat = expression_matrix, plotStats = FALSE, verbose = TRUE )
  if ("splitByBlocks" %in% formal_names) {
    arguments$splitByBlocks <- TRUE
  }
  if ("nCores" %in% formal_names) {
    arguments$nCores <- AUCELL_NCORES
  }
  if (!"..." %in% formal_names) {
    arguments <- arguments[ names(arguments) %in% formal_names ]
  }
  do.call( function_object, arguments )
}
calculate_auc_version_safe <- function( gene_sets, rankings, auc_max_rank) {
  function_object <- AUCell::AUCell_calcAUC
  formal_names <- names( formals( function_object ) )
  arguments <- list( geneSets = gene_sets, rankings = rankings, aucMaxRank = auc_max_rank, normAUC = TRUE, verbose = TRUE )
  if ("nCores" %in% formal_names) {
    arguments$nCores <- AUCELL_NCORES
  }
  if (!"..." %in% formal_names) {
    arguments <- arguments[ names(arguments) %in% formal_names ]
  }
  do.call( function_object, arguments )
}
compute_program_expression_metrics <- function( counts_matrix, genes, total_counts) {
  gene_indices <- match( genes, rownames( counts_matrix ) )
  if (anyNA(gene_indices)) {
    stop( "A mapped program gene is absent from the capture-area matrix.", call. = FALSE )
  }
  subset_counts <- counts_matrix[ gene_indices, , drop = FALSE ]
  normalization_factors <- LOGCP10K_SCALE_FACTOR / total_counts
  normalized_sparse <- subset_counts %*%
    Matrix::Diagonal( x = normalization_factors )
  normalized_sparse@x <- log1p( normalized_sparse@x )
  mean_logCP10K <- as.numeric( Matrix::colMeans( normalized_sparse ) )
  detected_fraction <- as.numeric( Matrix::colMeans( subset_counts > 0 ) )
  list( mean_logCP10K = mean_logCP10K, detected_fraction = detected_fraction )
}
capture_files <- function( capture_area) {
  area_dir <- file.path( CAPTURE_SCORE_ROOT, capture_area )
  list( area_dir = area_dir, success_flag = file.path( area_dir, "SUCCESS.flag" ), scores_rds = file.path( area_dir, paste0( capture_area, "_target_shared_scores.rds" ) ), scores_csv = file.path( area_dir, paste0( capture_area, "_target_shared_scores.csv.gz" ) ), audit_csv = file.path( area_dir, paste0( capture_area, "_score_audit.csv" ) ), report_txt = file.path( area_dir, paste0( capture_area, "_STEP12E_REPORT.txt" ) ) )
}
capture_is_complete <- function( capture_area) {
  files <- capture_files( capture_area )
  required <- c( files$success_flag, files$scores_rds, files$scores_csv, files$audit_csv, files$report_txt )
  all( file.exists(required) ) && all( file.info(required)$size > 0 )
}
write_capture_failure <- function( capture_area, reason, elapsed) {
  files <- capture_files( capture_area )
  dir.create( files$area_dir, recursive = TRUE, showWarnings = FALSE )
  if (file.exists( files$success_flag )) {
    unlink( files$success_flag, force = TRUE )
  }
}
fit_section_fixed_microglia_adjustment <- function( score, microglia_weight, section_id, score_label) {
  model_data <- data.frame( score = as.numeric(score), microglia_weight = as.numeric( microglia_weight ), section_id = factor( as.character( section_id ) ), stringsAsFactors = FALSE )
  valid <- stats::complete.cases( model_data ) & is.finite( model_data$score ) & is.finite( model_data$microglia_weight )
  if (sum(valid) < 100L) {
    stop( "Too few complete spots for ", score_label, " microglia adjustment.", call. = FALSE )
  }
  model <- stats::lm( score ~ microglia_weight + section_id, data = model_data[valid, , drop = FALSE] )
  residual_vector <- rep( NA_real_, nrow( model_data ) )
  residual_vector[valid] <- stats::residuals( model )
  coefficient_names <- names( stats::coef( model ) )
  microglia_coefficient <- if ( "microglia_weight" %in% coefficient_names ) {
    unname( stats::coef( model )[ "microglia_weight" ] )
  } else {
    NA_real_
  }
  model_summary <- summary( model )
  list( residual = residual_vector, residual_z = safe_zscore( residual_vector ), microglia_coefficient = microglia_coefficient, r_squared = unname( model_summary$r.squared ), adjusted_r_squared = unname( model_summary$adj.r.squared ), n_spots = sum(valid), n_sections = nlevels( droplevels( model_data$section_id[valid] ) ) )
}
add_within_section_columns <- function( score_table) {
  data <- data.table::copy( score_table )
  data[ , target_program_adjusted_z_within_section := safe_zscore( target_program_microglia_adjusted_residual ), by = section_id ]
  data[ , shared_prebranch_adjusted_z_within_section := safe_zscore( shared_prebranch_microglia_adjusted_residual ), by = section_id ]
  data[ , microglia_weight_quartile_within_section := {
      ranks <- data.table::frank( Microglia_RCTD_weight, ties.method = "average", na.last = "keep" )
      valid_n <- sum( is.finite( Microglia_RCTD_weight ) )
      output <- rep( NA_integer_, .N )
      if (valid_n > 0L) {
        scaled_rank <- ranks / valid_n
        output[ is.finite( scaled_rank ) ] <- pmin( 4L, pmax( 1L, ceiling( scaled_rank[ is.finite( scaled_rank ) ] * 4 ) ) )
      }
      output
    },
    by = section_id
  ]
  data[ , target_adjusted_quartile_within_section := {
      ranks <- data.table::frank( target_program_adjusted_z_within_section, ties.method = "average", na.last = "keep" )
      valid_n <- sum( is.finite( target_program_adjusted_z_within_section ) )
      output <- rep( NA_integer_, .N )
      if (valid_n > 0L) {
        scaled_rank <- ranks / valid_n
        output[ is.finite( scaled_rank ) ] <- pmin( 4L, pmax( 1L, ceiling( scaled_rank[ is.finite( scaled_rank ) ] * 4 ) ) )
      }
      output
    },
    by = section_id
  ]
  data
}
resolve_hcl_palette <- function( requested, fallback = "Viridis") {
  palettes <- grDevices::hcl.pals()
  exact <- palettes[ tolower(palettes) == tolower(requested) ]
  if (length(exact) > 0L) {
    return(exact[[1L]])
  }
  fallback_match <- palettes[ tolower(palettes) == tolower(fallback) ]
  if (length(fallback_match) > 0L) {
    return(fallback_match[[1L]])
  }
  palettes[[1L]]
}
continuous_colors <- function( values, palette = "Viridis") {
  values <- as.numeric(values)
  finite_values <- values[ is.finite(values) ]
  if (length(finite_values) == 0L) {
    return( rep( "#808080AA", length(values) ) )
  }
  lower <- as.numeric( stats::quantile( finite_values, probs = PLOT_QUANTILE_LOW, names = FALSE, na.rm = TRUE, type = 8 ) )
  upper <- as.numeric( stats::quantile( finite_values, probs = PLOT_QUANTILE_HIGH, names = FALSE, na.rm = TRUE, type = 8 ) )
  if ( !is.finite(lower) || !is.finite(upper) || upper <= lower ) {
    lower <- min( finite_values )
    upper <- max( finite_values )
  }
  if ( !is.finite(lower) || !is.finite(upper) || upper <= lower ) {
    return( rep( "#3B528BAA", length(values) ) )
  }
  clipped <- pmin( upper, pmax( lower, values ) )
  scaled <- ( clipped - lower ) / ( upper - lower )
  palette_name <- resolve_hcl_palette( palette )
  palette_values <- grDevices::hcl.colors( 100L, palette = palette_name )
  indices <- floor( scaled * 99 ) + 1L
  output <- rep( "#808080AA", length(values) )
  valid <- is.finite( indices )
  output[valid] <- grDevices::adjustcolor( palette_values[ pmax( 1L, pmin( 100L, indices[valid] ) ) ], alpha.f = 0.80 )
  output
}
draw_image_panel <- function( image_array, image_width, image_height, x, y, values, title, palette = "Viridis") {
  graphics::plot.new()
  graphics::plot.window( xlim = c( 0, image_width ), ylim = c( image_height, 0 ), asp = 1 )
  graphics::rasterImage( image_array, xleft = 0, ybottom = image_height, xright = image_width, ytop = 0, interpolate = FALSE )
  graphics::points( x, y, pch = 16, cex = POINT_CEX, col = continuous_colors( values, palette = palette ) )
  graphics::box( col = "#666666" )
  graphics::title( main = title, cex.main = 0.92 )
  finite_values <- values[ is.finite(values) ]
  if (length(finite_values) > 0L) {
    graphics::mtext( paste0( "median=", round( stats::median( finite_values ), 4 ), "; q01–q99=", round( stats::quantile( finite_values, PLOT_QUANTILE_LOW, names = FALSE ), 4 ), "–", round( stats::quantile( finite_values, PLOT_QUANTILE_HIGH, names = FALSE ), 4 ) ), side = 1, line = -1.1, cex = 0.61 )
  }
}
checks <- data.table( check_id = character(), category = character(), description = character(), status = character(), observed = character(), required = character(), critical = logical(), action_if_not_pass = character() )
add_check <- function( check_id, category, description, passed, observed, required, critical = TRUE, action_if_not_pass = "") {
  checks <<- rbind( checks, data.table( check_id = as.character( check_id ), category = as.character( category ), description = as.character( description ), status = if ( isTRUE(passed) ) {
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
write_fatal_report <- function(reason) {
  write_run_log()
}
score_one_capture_area <- function( capture_area, spatial_counts, spot_qc, target_genes_spatial, shared_genes_spatial, auc_max_rank) {
  area_start <- Sys.time()
  files <- capture_files( capture_area )
  dir.create( files$area_dir, recursive = TRUE, showWarnings = FALSE )
  if ( capture_area %in% FORCE_RERUN_CAPTURE_AREAS ) {
    unlink( files$area_dir, recursive = TRUE, force = TRUE )
    dir.create( files$area_dir, recursive = TRUE, showWarnings = FALSE )
  }
  current_metadata <- spot_qc[ capture_area_id == capture_area ]
  expected_spots <- as.character( current_metadata$barcode )
  column_indices <- match( expected_spots, colnames( spatial_counts ) )
  if (anyNA( column_indices )) {
    stop( "Spatial BPCells matrix is missing one or more ", capture_area, " spot barcodes.", call. = FALSE )
  }
  area_view <- spatial_counts[ , column_indices, drop = FALSE ]
  area_counts <- to_dgCMatrix( area_view, paste0( capture_area, " full spatial counts" ) )
  rm(area_view)
  invisible(gc())
  if ( !identical( colnames( area_counts ), expected_spots ) ) {
    area_counts <- area_counts[ , match( expected_spots, colnames( area_counts ) ), drop = FALSE ]
  }
  if ( !identical( colnames( area_counts ), expected_spots ) ) {
    stop( "Capture-area count columns could not be aligned for ", capture_area, ".", call. = FALSE )
  }
  total_counts <- as.numeric( Matrix::colSums( area_counts ) )
  if (any( !is.finite(total_counts) | total_counts <= 0 )) {
    stop( capture_area, " contains zero/non-finite total-count spots.", call. = FALSE )
  }
  gene_sets <- list( target_program = target_genes_spatial, Shared_prebranch = shared_genes_spatial )
  rankings <- build_rankings_version_safe( area_counts )
  auc_object <- calculate_auc_version_safe( gene_sets = gene_sets, rankings = rankings, auc_max_rank = auc_max_rank )
  auc_matrix <- as.matrix( AUCell::getAUC( auc_object ) )
  expected_programs <- names( gene_sets )
  if ( !all( expected_programs %in% rownames( auc_matrix ) ) || !all( expected_spots %in% colnames( auc_matrix ) ) ) {
    stop( "AUCell output is missing a program or spot for ", capture_area, ".", call. = FALSE )
  }
  auc_matrix <- auc_matrix[ expected_programs, expected_spots, drop = FALSE ]
  target_expression_metrics <- compute_program_expression_metrics( counts_matrix = area_counts, genes = target_genes_spatial, total_counts = total_counts )
  shared_expression_metrics <- compute_program_expression_metrics( counts_matrix = area_counts, genes = shared_genes_spatial, total_counts = total_counts )
  score_table <- data.table( barcode = expected_spots, capture_area_id = capture_area, target_program_AUCell_raw = as.numeric( auc_matrix[ "target_program", expected_spots ] ), shared_prebranch_AUCell_raw = as.numeric( auc_matrix[ "Shared_prebranch", expected_spots ] ), target_program_mean_logCP10K = target_expression_metrics$mean_logCP10K, shared_prebranch_mean_logCP10K = shared_expression_metrics$mean_logCP10K, target_program_detected_gene_fraction = target_expression_metrics$detected_fraction, shared_prebranch_detected_gene_fraction = shared_expression_metrics$detected_fraction, total_spatial_counts = total_counts, aucMaxRank = auc_max_rank )
  score_valid <- nrow( score_table ) == length( expected_spots ) && !anyDuplicated( score_table$barcode ) && all( is.finite( score_table$target_program_AUCell_raw ) ) && all( is.finite( score_table$shared_prebranch_AUCell_raw ) ) && all( score_table$target_program_AUCell_raw >= 0 ) && all( score_table$shared_prebranch_AUCell_raw >= 0 ) && all( score_table$target_program_detected_gene_fraction >= 0 & score_table$target_program_detected_gene_fraction <= 1 ) && all( score_table$shared_prebranch_detected_gene_fraction >= 0 & score_table$shared_prebranch_detected_gene_fraction <= 1 )
  if (!score_valid) {
    stop( "Capture-area score validation failed for ", capture_area, ".", call. = FALSE )
  }
  saveRDS( score_table, files$scores_rds, compress = FALSE )
  data.table::fwrite( score_table, files$scores_csv, bom = TRUE, compress = "gzip" )
  area_elapsed <- elapsed_minutes( area_start )
  area_audit <- data.table( capture_area_id = capture_area, status = "PASS", elapsed_minutes = area_elapsed, n_genes_ranked = nrow( area_counts ), n_spots = ncol( area_counts ), n_nonzero = length( area_counts@x ), aucMaxRank = auc_max_rank, target_genes = length( target_genes_spatial ), shared_genes = length( shared_genes_spatial ), target_AUCell_min = min( score_table$target_program_AUCell_raw ), target_AUCell_median = median( score_table$target_program_AUCell_raw ), target_AUCell_max = max( score_table$target_program_AUCell_raw ), shared_AUCell_min = min( score_table$shared_prebranch_AUCell_raw ), shared_AUCell_median = median( score_table$shared_prebranch_AUCell_raw ), shared_AUCell_max = max( score_table$shared_prebranch_AUCell_raw ), target_detected_fraction_median = median( score_table$target_program_detected_gene_fraction ), shared_detected_fraction_median = median( score_table$shared_prebranch_detected_gene_fraction ), score_RDS = files$scores_rds )
  data.table::fwrite( area_audit, files$audit_csv, bom = TRUE )
  area_report <- c( "============================================================", paste0( "GSE234774 Step 12E | ", capture_area ), "============================================================", paste0( "Generated: ", format( Sys.time(), "%Y-%m-%d %H:%M:%S" ) ), "", "FINAL_STATUS: PASS", paste0( "ELAPSED_MINUTES: ", area_elapsed ), "", "AUCELL", paste0( "- Genes ranked: ", nrow( area_counts ) ), paste0( "- Spots scored: ", ncol( area_counts ) ), paste0( "- aucMaxRank: ", auc_max_rank, " (", AUC_MAX_RANK_FRACTION * 100, "% of spatial genes)" ), paste0( "- Target genes: ", length( target_genes_spatial ) ), paste0( "- Shared_prebranch genes: ", length( shared_genes_spatial ) ), "", "BOUNDARIES", "- Scores are descriptive spot-level measurements.", "- Microglia adjustment is performed only after all capture areas are combined.", "- No lesion-zone or recipient-neighborhood test was performed.", "============================================================" )
  writeLines( area_report, files$report_txt, useBytes = TRUE )
  writeLines( c( paste0( "capture_area_id=", capture_area ), paste0( "completed_at=", format( Sys.time(), "%Y-%m-%d %H:%M:%S" ) ), "status=PASS" ), files$success_flag, useBytes = TRUE )
  rm( rankings, auc_object, auc_matrix, area_counts )
  invisible(gc())
  area_audit
}
main <- function() {
  workflow_start <- Sys.time()
  required_inputs <- c( Step11B_report = STEP11B_REPORT, Step11B_results_RDS = STEP11B_RESULTS_RDS, Step12D_report = STEP12D_REPORT, Step12D_normalized_weights_RDS = STEP12D_NORMALIZED_WEIGHTS_RDS, Step12B_report = STEP12B_REPORT, Step12B_spot_QC = STEP12B_SPOT_QC, Step12B_capture_audit = STEP12B_CAPTURE_AUDIT )
  input_exists <- file.exists( required_inputs )
  spatial_cache_exists <- dir.exists( SPATIAL_BPCELLS_DIR )
  add_check( "C01", "input", "All frozen Step 11B/12B/12D inputs and the spatial BPCells matrix exist.", all(input_exists) && spatial_cache_exists, paste0( paste( paste0( names( required_inputs ), "=", input_exists ), collapse = "; " ), "; spatial_BPCells=", spatial_cache_exists ), "All TRUE", TRUE, "Restore frozen upstream outputs; reconstruct gene sets from memory." )
  if ( !all(input_exists) || !spatial_cache_exists ) {
    stop( "At least one required Step 12E input is missing.", call. = FALSE )
  }
  step11B_text <- readLines( STEP11B_REPORT, warn = FALSE, encoding = "UTF-8" )
  step11B_status_lines <- step11B_text[ grepl( "^FINAL DECISION:", step11B_text ) | grepl( "^FINAL_STATUS:", step11B_text ) ]
  step11B_pass <- length( step11B_status_lines ) >= 1L && any( grepl( "PASS", step11B_status_lines ) )
  step12B_text <- readLines( STEP12B_REPORT, warn = FALSE, encoding = "UTF-8" )
  step12B_status_lines <- step12B_text[ grepl( "^FINAL_STATUS:", step12B_text ) ]
  step12B_pass <- length( step12B_status_lines ) == 1L && grepl( "^FINAL_STATUS: (PASS|PASS_WITH_REVIEW)$", step12B_status_lines )
  step12D_text <- readLines( STEP12D_REPORT, warn = FALSE, encoding = "UTF-8" )
  step12D_status_lines <- step12D_text[ grepl( "^FINAL_STATUS:", step12D_text ) ]
  step12D_pass <- length( step12D_status_lines ) == 1L && grepl( "^FINAL_STATUS: (PASS|PASS_WITH_REVIEW)$", step12D_status_lines )
  add_check( "C02", "upstream", "Step 11B, Step 12B, and Step 12D release target scoring.", step11B_pass && step12B_pass && step12D_pass, paste0( "Step11B=", collapse_unique( step11B_status_lines ), "; Step12B=", collapse_unique( step12B_status_lines ), "; Step12D=", collapse_unique( step12D_status_lines ) ), "All upstream stages PASS/PASS_WITH_REVIEW", TRUE, "score spatial target programs before the frozen stages pass." )
  if ( !step11B_pass || !step12B_pass || !step12D_pass ) {
    stop( "An upstream frozen stage does not release Step 12E.", call. = FALSE )
  }
  step11B_results <- readRDS( STEP11B_RESULTS_RDS )
  required_step11B_fields <- c( "scored_gene_sets_present", "scored_gene_sets_requested", "source_sets" )
  missing_step11B_fields <- setdiff( required_step11B_fields, names( step11B_results ) )
  if (length( missing_step11B_fields ) > 0L) {
    stop( "Step 11B results RDS is missing field(s): ", paste( missing_step11B_fields, collapse = ", " ), call. = FALSE )
  }
  frozen_gene_sets <- step11B_results$scored_gene_sets_present
  required_gene_sets <- c( "Arm_2_specific", "Shared_prebranch" )
  if ( !is.list( frozen_gene_sets ) || !all( required_gene_sets %in% names( frozen_gene_sets ) ) ) {
    stop( "Step 11B scored_gene_sets_present does not contain Arm_2_specific and Shared_prebranch.", call. = FALSE )
  }
  target_genes_step11B <- clean_gene_vector( frozen_gene_sets$Arm_2_specific )
  shared_genes_step11B <- clean_gene_vector( frozen_gene_sets$Shared_prebranch )
  frozen_sizes_ok <- length( target_genes_step11B ) == EXPECTED_TARGET_STEP11B_GENES && length( shared_genes_step11B ) == EXPECTED_SHARED_STEP11B_GENES && length( intersect( target_genes_step11B, shared_genes_step11B ) ) == 0L
  add_check( "C03", "gene_set", "Frozen Step 11B target and Shared_prebranch gene sets have expected sizes and are disjoint.", frozen_sizes_ok, paste0( "target=", length( target_genes_step11B ), "; shared=", length( shared_genes_step11B ), "; overlap=", length( intersect( target_genes_step11B, shared_genes_step11B ) ) ), paste0( "target=", EXPECTED_TARGET_STEP11B_GENES, "; shared=", EXPECTED_SHARED_STEP11B_GENES, "; overlap=0" ), TRUE, "redefine the frozen programs in the spatial dataset." )
  if (!frozen_sizes_ok) {
    stop( "Frozen Step 11B gene-set integrity check failed.", call. = FALSE )
  }
  spatial_counts <- BPCells::open_matrix_dir( SPATIAL_BPCELLS_DIR )
  spatial_dimensions_ok <- nrow( spatial_counts ) == EXPECTED_SPATIAL_GENES && ncol( spatial_counts ) == EXPECTED_SPATIAL_SPOTS && !anyDuplicated( rownames( spatial_counts ) ) && !anyDuplicated( colnames( spatial_counts ) )
  add_check( "C04", "spatial_matrix", "The frozen spatial BPCells matrix has expected dimensions and unique names.", spatial_dimensions_ok, paste0( "dimensions=", nrow( spatial_counts ), " x ", ncol( spatial_counts ), "; duplicated genes=", anyDuplicated( rownames( spatial_counts ) ), "; duplicated spots=", anyDuplicated( colnames( spatial_counts ) ) ), "25,966 genes x 33,941 spots with unique names", TRUE, "score a changed spatial matrix." )
  if (!spatial_dimensions_ok) {
    stop( "Spatial BPCells matrix integrity check failed.", call. = FALSE )
  }
  target_resolution <- resolve_features( requested_features = target_genes_step11B, available_features = rownames( spatial_counts ) )
  shared_resolution <- resolve_features( requested_features = shared_genes_step11B, available_features = rownames( spatial_counts ) )
  target_genes_spatial <- unique( unname( target_resolution[ !is.na( target_resolution ) ] ) )
  shared_genes_spatial <- unique( unname( shared_resolution[ !is.na( shared_resolution ) ] ) )
  target_mapping_fraction <- length( target_genes_spatial ) / length( target_genes_step11B )
  shared_mapping_fraction <- length( shared_genes_spatial ) / length( shared_genes_step11B )
  spatial_program_overlap <- intersect( target_genes_spatial, shared_genes_spatial )
  mapping_ok <- length( target_genes_spatial ) >= MIN_TARGET_GENES_PRESENT && length( shared_genes_spatial ) >= MIN_SHARED_GENES_PRESENT && target_mapping_fraction >= MIN_TARGET_MAPPING_FRACTION && shared_mapping_fraction >= MIN_SHARED_MAPPING_FRACTION && length( spatial_program_overlap ) == 0L
  add_check( "C05", "gene_set", "Target and Shared_prebranch map sufficiently to spatial genes and remain separate.", mapping_ok, paste0( "target=", length( target_genes_spatial ), "/", length( target_genes_step11B ), " (", round( target_mapping_fraction, 4 ), "); shared=", length( shared_genes_spatial ), "/", length( shared_genes_step11B ), " (", round( shared_mapping_fraction, 4 ), "); spatial overlap=", length( spatial_program_overlap ) ), "Both mapping fractions >=0.90 and overlap=0", TRUE, "optimize or replace missing genes using spatial expression." )
  if (!mapping_ok) {
    stop( "Frozen-program spatial mapping failed.", call. = FALSE )
  }
  target_gene_audit <- data.table( program = "target_program", Step11B_name = "Arm_2_specific", requested_gene = names( target_resolution ), spatial_gene = unname( target_resolution ), mapping_status = ifelse( is.na( target_resolution ), "NOT_PRESENT", "PRESENT" ) )
  shared_gene_audit <- data.table( program = "Shared_prebranch", Step11B_name = "Shared_prebranch", requested_gene = names( shared_resolution ), spatial_gene = unname( shared_resolution ), mapping_status = ifelse( is.na( shared_resolution ), "NOT_PRESENT", "PRESENT" ) )
  package_input_audit <- data.table( program = "AUDIT_METADATA", Step11B_name = c( "R_version", required_packages, "Step11B_results_RDS", "Step12D_normalized_weights_RDS", "spatial_BPCells_path", "target_mapping_fraction", "shared_mapping_fraction", "AUCell_aucMaxRank_fraction" ), requested_gene = c( R.version.string, vapply( required_packages, function(package_name) {
          as.character( utils::packageVersion( package_name ) )
        },
        FUN.VALUE = character(1)
      ), STEP11B_RESULTS_RDS, STEP12D_NORMALIZED_WEIGHTS_RDS, SPATIAL_BPCELLS_DIR, as.character(
        target_mapping_fraction
      ), as.character(
        shared_mapping_fraction
      ), as.character(
        AUC_MAX_RANK_FRACTION
      )
    ), spatial_gene = NA_character_, mapping_status = "AUDIT"
  )
  input_gene_audit <- rbindlist( list( package_input_audit, target_gene_audit, shared_gene_audit ), fill = TRUE )
  data.table::fwrite( input_gene_audit, OUT_INPUT_GENE_AUDIT, bom = TRUE )
  spot_qc <- data.table::fread( STEP12B_SPOT_QC, data.table = TRUE, showProgress = FALSE )
  required_columns( spot_qc, c( "barcode", "capture_area_id", "section_id", "timepoint", "replicate_index", "biological_replicate_id", "injury_condition", "nCount_Spatial_BPCells", "nFeature_Spatial_BPCells", "percent.mt", "qc_review_flag" ), "Step 12B spot-QC table" )
  spot_qc[ , capture_area_id := toupper( trimws( as.character( capture_area_id ) ) ) ]
  spot_qc[ , qc_review_flag := as_logical_strict( qc_review_flag, "qc_review_flag" ) ]
  normalized_weights <- readRDS( STEP12D_NORMALIZED_WEIGHTS_RDS )
  normalized_weights <- as.matrix( normalized_weights )
  normalized_weight_ok <- nrow( normalized_weights ) == EXPECTED_SPATIAL_SPOTS && "Microglia" %in%
      colnames( normalized_weights ) && !anyDuplicated( rownames( normalized_weights ) ) && all( is.finite( normalized_weights ) ) && max( abs( rowSums( normalized_weights ) - 1 ) ) <= 1e-6
  metadata_weight_alignment <- nrow( spot_qc ) == EXPECTED_SPATIAL_SPOTS && !anyDuplicated( spot_qc$barcode ) && setequal( spot_qc$barcode, rownames( normalized_weights ) )
  add_check( "C06", "RCTD", "Normalized RCTD weights are valid and align exactly to all spatial barcodes.", normalized_weight_ok && metadata_weight_alignment, paste0( "weights=", nrow( normalized_weights ), " x ", ncol( normalized_weights ), "; Microglia present=", "Microglia" %in% colnames( normalized_weights ), "; metadata rows=", nrow( spot_qc ), "; barcode set equal=", setequal( spot_qc$barcode, rownames( normalized_weights ) ) ), "33,941 aligned spots and Microglia weight", TRUE, "adjust target scores with misaligned RCTD weights." )
  if ( !normalized_weight_ok || !metadata_weight_alignment ) {
    stop( "Step 12D normalized-weight alignment failed.", call. = FALSE )
  }
  spot_qc <- spot_qc[ match( rownames( normalized_weights ), barcode ) ]
  if ( anyNA( spot_qc$barcode ) || !identical( spot_qc$barcode, rownames( normalized_weights ) ) ) {
    stop( "Spot metadata could not be reordered to normalized RCTD weights.", call. = FALSE )
  }
  spot_qc[ , Microglia_RCTD_weight := normalized_weights[ barcode, "Microglia" ] ]
  all_capture_areas <- sort( unique( spot_qc$capture_area_id ) )
  capture_structure_ok <- length( all_capture_areas ) == EXPECTED_CAPTURE_AREAS && data.table::uniqueN( spot_qc$section_id ) == EXPECTED_SECTIONS
  add_check( "C07", "metadata", "Spot metadata contains 12 capture areas and 36 sections.", capture_structure_ok, paste0( "capture areas=", length( all_capture_areas ), "; sections=", uniqueN( spot_qc$section_id ) ), "12 capture areas and 36 sections", TRUE, "perform section-fixed adjustment with incomplete metadata." )
  if (!capture_structure_ok) {
    stop( "Spatial capture-area/section structure is incomplete.", call. = FALSE )
  }
  auc_max_rank <- max( 1L, floor( nrow( spatial_counts ) * AUC_MAX_RANK_FRACTION ) )
  capture_areas_to_run <- if ( length( CAPTURE_AREAS_TO_RUN ) == 1L && identical( toupper( CAPTURE_AREAS_TO_RUN ), "ALL" ) ) {
    all_capture_areas
  } else {
    toupper( as.character( CAPTURE_AREAS_TO_RUN ) )
  }
  unknown_capture_areas <- setdiff( capture_areas_to_run, all_capture_areas )
  if (length( unknown_capture_areas ) > 0L) {
    stop( "Unknown capture area(s): ", paste( unknown_capture_areas, collapse = " | " ), call. = FALSE )
  }
  area_audit_list <- vector( "list", length( capture_areas_to_run ) )
  for ( area_index in seq_along( capture_areas_to_run ) ) {
    capture_area <- capture_areas_to_run[ area_index ]
    if ( RESUME_COMPLETED && !capture_area %in% FORCE_RERUN_CAPTURE_AREAS && capture_is_complete( capture_area ) ) {
      files <- capture_files( capture_area )
      prior_audit <- data.table::fread( files$audit_csv, data.table = TRUE, showProgress = FALSE )
      prior_audit[ , run_action := "SKIPPED_EXISTING_SUCCESS" ]
      area_audit_list[[ area_index ]] <- prior_audit
      next
    }
    area_start <- Sys.time()
    area_result <- tryCatch( score_one_capture_area( capture_area = capture_area, spatial_counts = spatial_counts, spot_qc = spot_qc, target_genes_spatial = target_genes_spatial, shared_genes_spatial = shared_genes_spatial, auc_max_rank = auc_max_rank ), error = function(e) {
        write_capture_failure( capture_area = capture_area, reason = conditionMessage( e ), elapsed = elapsed_minutes( area_start ) )
        data.table( capture_area_id = capture_area, status = "FAIL", elapsed_minutes = elapsed_minutes( area_start ), error_message = conditionMessage( e ), run_action = "ATTEMPTED_FAILED" )
      }
    )
    if (!"run_action" %in% names( area_result )) {
      area_result[ , run_action := "ATTEMPTED_PASS" ]
    }
    area_audit_list[[ area_index ]] <- area_result
    invisible(gc())
  }
  completed_capture_areas <- all_capture_areas[ vapply( all_capture_areas, capture_is_complete, logical(1) ) ]
  incomplete_capture_areas <- setdiff( all_capture_areas, completed_capture_areas )
  all_scores_complete <- length( incomplete_capture_areas ) == 0L
  add_check( "C08", "scoring_completion", "All 12 capture areas have successful AUCell score checkpoints.", all_scores_complete, paste0( "complete=", length( completed_capture_areas ), "/12; incomplete=", paste( incomplete_capture_areas, collapse = " | " ) ), "12/12 complete", TRUE, "Rerun the script; completed capture areas will be skipped." )
  if (!all_scores_complete) {
    data.table::fwrite( checks, OUT_CHECKS, bom = TRUE )
    stop( "Step 12E has incomplete capture-area score checkpoints: ", paste( incomplete_capture_areas, collapse = " | " ), call. = FALSE )
  }
  score_list <- lapply( completed_capture_areas, function(capture_area) {
      readRDS( capture_files( capture_area )$scores_rds )
    }
  )
  combined_scores <- rbindlist( score_list, fill = TRUE )
  combined_scores <- combined_scores[ match( spot_qc$barcode, barcode ) ]
  combined_score_ok <- nrow( combined_scores ) == EXPECTED_SPATIAL_SPOTS && !anyDuplicated( combined_scores$barcode ) && identical( combined_scores$barcode, spot_qc$barcode ) && all( is.finite( combined_scores$target_program_AUCell_raw ) ) && all( is.finite( combined_scores$shared_prebranch_AUCell_raw ) )
  add_check( "C09", "combined_score", "Combined target and Shared_prebranch scores cover all spots exactly once.", combined_score_ok, paste0( "rows=", nrow( combined_scores ), "; duplicates=", anyDuplicated( combined_scores$barcode ), "; exact metadata order=", identical( combined_scores$barcode, spot_qc$barcode ) ), "33,941 unique aligned spots", TRUE, "fit Microglia adjustment with incomplete score coverage." )
  if (!combined_score_ok) {
    stop( "Combined capture-area score validation failed.", call. = FALSE )
  }
  spot_scores <- cbind( spot_qc[ , .( barcode, capture_area_id, section_id, slide_id, timepoint, replicate_index, biological_replicate_id, injury_condition, section_order, x, y, z, nCount_Spatial_BPCells, nFeature_Spatial_BPCells, percent.mt, qc_review_flag, Microglia_RCTD_weight ) ], combined_scores[ , !c( "barcode", "capture_area_id" ) ] )
  target_adjustment <- fit_section_fixed_microglia_adjustment( score = spot_scores$target_program_AUCell_raw, microglia_weight = spot_scores$Microglia_RCTD_weight, section_id = spot_scores$section_id, score_label = "target program" )
  shared_adjustment <- fit_section_fixed_microglia_adjustment( score = spot_scores$shared_prebranch_AUCell_raw, microglia_weight = spot_scores$Microglia_RCTD_weight, section_id = spot_scores$section_id, score_label = "Shared_prebranch" )
  spot_scores[ , target_program_microglia_adjusted_residual := target_adjustment$residual ]
  spot_scores[ , target_program_microglia_adjusted_z := target_adjustment$residual_z ]
  spot_scores[ , shared_prebranch_microglia_adjusted_residual := shared_adjustment$residual ]
  spot_scores[ , shared_prebranch_microglia_adjusted_z := shared_adjustment$residual_z ]
  spot_scores <- add_within_section_columns( spot_scores )
  adjusted_scores_valid <- all( is.finite( spot_scores$target_program_microglia_adjusted_z ) ) && all( is.finite( spot_scores$shared_prebranch_microglia_adjusted_z ) ) && abs( mean( spot_scores$target_program_microglia_adjusted_z ) ) <
      1e-8 && abs( stats::sd( spot_scores$target_program_microglia_adjusted_z ) - 1 ) <
      1e-6
  add_check( "C10", "microglia_adjustment", "Target and Shared_prebranch Microglia-adjusted scores are finite for every spot.", adjusted_scores_valid, paste0( "target finite=", sum( is.finite( spot_scores$target_program_microglia_adjusted_z ) ), "/", nrow( spot_scores ), "; shared finite=", sum( is.finite( spot_scores$shared_prebranch_microglia_adjusted_z ) ), "/", nrow( spot_scores ), "; target adjusted mean=", mean( spot_scores$target_program_microglia_adjusted_z ), "; target adjusted SD=", stats::sd( spot_scores$target_program_microglia_adjusted_z ) ), "33,941 finite adjusted scores; target adjusted mean~0 and SD~1", TRUE, "proceed to lesion-zone analysis with failed residualization." )
  if (!adjusted_scores_valid) {
    stop( "Microglia-adjusted score validation failed.", call. = FALSE )
  }
  global_adjustment_audit <- data.table( record_type = "GLOBAL_MODEL", section_id = "ALL", capture_area_id = "ALL", timepoint = "ALL", n_spots = nrow( spot_scores ), target_raw_vs_Microglia_spearman = safe_spearman( spot_scores$target_program_AUCell_raw, spot_scores$Microglia_RCTD_weight ), target_adjusted_vs_Microglia_spearman = safe_spearman( spot_scores$target_program_microglia_adjusted_z, spot_scores$Microglia_RCTD_weight ), shared_raw_vs_Microglia_spearman = safe_spearman( spot_scores$shared_prebranch_AUCell_raw, spot_scores$Microglia_RCTD_weight ), shared_adjusted_vs_Microglia_spearman = safe_spearman( spot_scores$shared_prebranch_microglia_adjusted_z, spot_scores$Microglia_RCTD_weight ), target_microglia_model_coefficient = target_adjustment$microglia_coefficient, target_model_R2 = target_adjustment$r_squared, target_model_adjusted_R2 = target_adjustment$adjusted_r_squared, shared_microglia_model_coefficient = shared_adjustment$microglia_coefficient, shared_model_R2 = shared_adjustment$r_squared, shared_model_adjusted_R2 = shared_adjustment$adjusted_r_squared, inferential_p_values_run = FALSE )
  section_adjustment_audit <- spot_scores[ , .( record_type = "SECTION_DESCRIPTIVE", capture_area_id = as.character( capture_area_id[1L] ), timepoint = as.character( timepoint[1L] ), n_spots = .N, target_raw_vs_Microglia_spearman = safe_spearman( target_program_AUCell_raw, Microglia_RCTD_weight ), target_adjusted_vs_Microglia_spearman = safe_spearman( target_program_microglia_adjusted_z, Microglia_RCTD_weight ), shared_raw_vs_Microglia_spearman = safe_spearman( shared_prebranch_AUCell_raw, Microglia_RCTD_weight ), shared_adjusted_vs_Microglia_spearman = safe_spearman( shared_prebranch_microglia_adjusted_z, Microglia_RCTD_weight ), target_microglia_model_coefficient = NA_real_, target_model_R2 = NA_real_, target_model_adjusted_R2 = NA_real_, shared_microglia_model_coefficient = NA_real_, shared_model_R2 = NA_real_, shared_model_adjusted_R2 = NA_real_, inferential_p_values_run = FALSE ), by = section_id ]
  adjustment_audit <- rbind( global_adjustment_audit, section_adjustment_audit, fill = TRUE )
  section_summary <- spot_scores[ , .( n_spots = .N, n_QC_review_spots = sum( qc_review_flag ), Microglia_weight_mean = mean( Microglia_RCTD_weight ), Microglia_weight_median = median( Microglia_RCTD_weight ), target_AUCell_mean = mean( target_program_AUCell_raw ), target_AUCell_median = median( target_program_AUCell_raw ), target_AUCell_q10 = as.numeric( quantile( target_program_AUCell_raw, 0.10, names = FALSE ) ), target_AUCell_q90 = as.numeric( quantile( target_program_AUCell_raw, 0.90, names = FALSE ) ), target_adjusted_z_mean = mean( target_program_microglia_adjusted_z ), target_adjusted_z_median = median( target_program_microglia_adjusted_z ), target_adjusted_z_q10 = as.numeric( quantile( target_program_microglia_adjusted_z, 0.10, names = FALSE ) ), target_adjusted_z_q90 = as.numeric( quantile( target_program_microglia_adjusted_z, 0.90, names = FALSE ) ), shared_AUCell_mean = mean( shared_prebranch_AUCell_raw ), shared_AUCell_median = median( shared_prebranch_AUCell_raw ), shared_adjusted_z_mean = mean( shared_prebranch_microglia_adjusted_z ), shared_adjusted_z_median = median( shared_prebranch_microglia_adjusted_z ), target_detected_gene_fraction_median = median( target_program_detected_gene_fraction ), shared_detected_gene_fraction_median = median( shared_prebranch_detected_gene_fraction ), target_raw_vs_Microglia_spearman = safe_spearman( target_program_AUCell_raw, Microglia_RCTD_weight ), target_adjusted_vs_Microglia_spearman = safe_spearman( target_program_microglia_adjusted_z, Microglia_RCTD_weight ) ), by = .( capture_area_id, section_id, timepoint, replicate_index, biological_replicate_id, injury_condition ) ]
  section_summary_ok <- nrow( section_summary ) == EXPECTED_SECTIONS && sum( section_summary$n_spots ) == EXPECTED_SPATIAL_SPOTS
  add_check( "C11", "section_summary", "Section-level descriptive summaries preserve 36 sections and all spots.", section_summary_ok, paste0( "sections=", nrow( section_summary ), "; summed spots=", sum( section_summary$n_spots ) ), "36 sections and 33,941 summed spots", TRUE, "use an incomplete section summary in Step 12F." )
  if (!section_summary_ok) {
    stop( "Section summary integrity check failed.", call. = FALSE )
  }
  data.table::fwrite( spot_scores, OUT_SPOT_SCORES, bom = TRUE, compress = "gzip" )
  data.table::fwrite( section_summary, OUT_SECTION_SUMMARY, bom = TRUE )
  data.table::fwrite( adjustment_audit, OUT_ADJUSTMENT_AUDIT, bom = TRUE )
  saveRDS( spot_scores, OUT_SPOT_SCORES_RDS, compress = FALSE )
  capture_audit <- data.table::fread( STEP12B_CAPTURE_AUDIT, data.table = TRUE, showProgress = FALSE )
  required_columns( capture_audit, c( "capture_area_id", "tissue_lowres_scalef", "lowres_image_width", "lowres_image_height", "asset_linkage_status" ), "Step 12B capture-area image audit" )
  capture_audit[ , capture_area_id := toupper( trimws( as.character( capture_area_id ) ) ) ]
  grDevices::pdf( OUT_SPATIAL_MAPS, width = 11, height = 8.5, onefile = TRUE, useDingbats = FALSE )
  pdf_open <- TRUE
  on.exit( {
      if (isTRUE(pdf_open)) {
        try( grDevices::dev.off(), silent = TRUE )
      }
    },
    add = TRUE
  )
  plot_area_pass <- logical( length( all_capture_areas ) )
  for ( area_index in seq_along( all_capture_areas ) ) {
    capture_area <- all_capture_areas[ area_index ]
    current_scores <- spot_scores[ capture_area_id == capture_area ]
    coordinate_file <- file.path( STEP12D_CAPTURE_ROOT, capture_area, paste0( capture_area, "_spot_metadata.csv.gz" ) )
    image_file <- file.path( LOWRES_PREVIEW_ROOT, capture_area, "tissue_lowres_image.png" )
    audit_row <- capture_audit[ capture_area_id == capture_area ]
    if ( nrow( audit_row ) != 1L || !file.exists( coordinate_file ) || !file.exists( image_file ) ) {
      stop( "Missing plotting asset for ", capture_area, ".", call. = FALSE )
    }
    coordinates <- data.table::fread( coordinate_file, data.table = TRUE, showProgress = FALSE )
    required_columns( coordinates, c( "barcode", "coordinate_x_fullres", "coordinate_y_fullres" ), paste0( capture_area, " coordinate table" ) )
    coordinates <- coordinates[ match( current_scores$barcode, barcode ) ]
    if ( anyNA( coordinates$barcode ) || !identical( coordinates$barcode, current_scores$barcode ) ) {
      stop( "Coordinate-score barcode alignment failed for ", capture_area, ".", call. = FALSE )
    }
    image_array <- png::readPNG( image_file, native = FALSE, info = TRUE )
    image_dimensions <- dim( image_array )
    image_height <- image_dimensions[1L]
    image_width <- image_dimensions[2L]
    lowres_scale <- as.numeric( audit_row$tissue_lowres_scalef )
    image_x <- as.numeric( coordinates$coordinate_x_fullres ) * lowres_scale
    image_y <- as.numeric( coordinates$coordinate_y_fullres ) * lowres_scale
    plot_alignment_ok <- all( is.finite( image_x ) ) && all( is.finite( image_y ) ) && all( image_x >= 0 & image_x <= image_width + 2 ) && all( image_y >= 0 & image_y <= image_height + 2 )
    plot_area_pass[ area_index ] <- plot_alignment_ok
    if (!plot_alignment_ok) {
      stop( "Image-coordinate bounds failed for ", capture_area, ".", call. = FALSE )
    }
    graphics::par( mfrow = c( 2, 2 ), mar = c( 0.5, 0.5, 2.4, 0.5 ), oma = c( 0.4, 0.4, 2.0, 0.4 ), xaxs = "i", yaxs = "i" )
    draw_image_panel( image_array = image_array, image_width = image_width, image_height = image_height, x = image_x, y = image_y, values = current_scores$target_program_AUCell_raw, title = "Target program AUCell (raw)", palette = "Inferno" )
    draw_image_panel( image_array = image_array, image_width = image_width, image_height = image_height, x = image_x, y = image_y, values = current_scores$Microglia_RCTD_weight, title = "Microglia RCTD weight", palette = "Viridis" )
    draw_image_panel( image_array = image_array, image_width = image_width, image_height = image_height, x = image_x, y = image_y, values = current_scores$target_program_microglia_adjusted_z, title = "Target score adjusted for Microglia + section", palette = "Blue-Red 3" )
    draw_image_panel( image_array = image_array, image_width = image_width, image_height = image_height, x = image_x, y = image_y, values = current_scores$shared_prebranch_AUCell_raw, title = "Shared_prebranch AUCell (separate)", palette = "Mako" )
    graphics::mtext( paste0( "GSE234774 Step 12E | ", capture_area, " | no lesion-zone inference" ), side = 3, outer = TRUE, line = 0.35, cex = 1.02, font = 2 )
    rm( image_array, coordinates )
    invisible(gc())
  }
  grDevices::dev.off()
  pdf_open <- FALSE
  plot_output_ok <- all( plot_area_pass ) && file.exists( OUT_SPATIAL_MAPS ) && file_size_safe( OUT_SPATIAL_MAPS ) > 0
  add_check( "C12", "spatial_plot", "All 12 capture areas produced aligned four-panel spatial score maps.", plot_output_ok, paste0( "aligned areas=", sum( plot_area_pass ), "/12; PDF bytes=", file_size_safe( OUT_SPATIAL_MAPS ) ), "12/12 aligned and nonempty PDF", TRUE, "interpret maps with coordinate or image misalignment." )
  if (!plot_output_ok) {
    stop( "Spatial score map generation failed.", call. = FALSE )
  }
  global_before <- global_adjustment_audit$
    target_raw_vs_Microglia_spearman
  global_after <- global_adjustment_audit$
    target_adjusted_vs_Microglia_spearman
  if ( is.finite( global_before ) && is.finite( global_after ) && abs( global_after ) >= abs( global_before ) ) {
    add_review( "R01", "microglia_adjustment", "Global target–Microglia Spearman correlation was not reduced in magnitude after adjustment.", paste0( "before=", round( global_before, 4 ), "; after=", round( global_after, 4 ) ), "Review scatterplots/section patterns before program-specific localization claims.", "interpret high raw target scores as independent of Microglia abundance." )
  }
  handoff <- list( step = "Step12E", completed_at = as.character( Sys.time() ), frozen_target_definition = "Step11B Arm_2_specific renamed persistent inflammatory-remodeling target program", shared_definition = "Step11B Shared_prebranch retained separately", target_genes_Step11B = target_genes_step11B, target_genes_spatial = target_genes_spatial, shared_genes_Step11B = shared_genes_step11B, shared_genes_spatial = shared_genes_spatial, primary_raw_score = "target_program_AUCell_raw", primary_adjusted_score = "target_program_microglia_adjusted_z", adjustment_model = "target_AUCell ~ Microglia_RCTD_weight + section fixed effects", spot_score_RDS = OUT_SPOT_SCORES_RDS, spot_score_CSV = OUT_SPOT_SCORES, section_summary_CSV = OUT_SECTION_SUMMARY, spatial_maps_PDF = OUT_SPATIAL_MAPS, n_spots = nrow( spot_scores ), n_sections = uniqueN( spot_scores$section_id ), lesion_zone_analysis_run = FALSE, target_name_upgraded_to_scar_associated = FALSE, next_step = paste( "Step 12F defines lesion core/border/perilesional/remote zones", "using section-aware spatial evidence and evaluates raw and", "Microglia-adjusted target scores with section/animal as the", "biological unit." ) )
  saveRDS( handoff, OUT_HANDOFF_RDS, compress = FALSE )
  critical_failures <- checks[ critical == TRUE & status == "FAIL" ]
  review_items <- checks[ status == "REVIEW" ]
  final_status <- if ( nrow( critical_failures ) > 0L ) {
    "STOP"
  } else if (
    nrow( review_items ) > 0L
  ) {
    "PASS_WITH_REVIEW"
  } else {
    "PASS"
  }
  data.table::fwrite( checks, OUT_CHECKS, bom = TRUE )
  total_elapsed <- elapsed_minutes( workflow_start )
  if (nrow( critical_failures ) > 0L) {
  }
  if (nrow( review_items ) > 0L) {
  }
  writeLines( capture.output( sessionInfo() ), SESSION_INFO_FILE, useBytes = TRUE )
  write_run_log()
  invisible( list( final_status = final_status, spot_scores_RDS = OUT_SPOT_SCORES_RDS, handoff_RDS = OUT_HANDOFF_RDS, section_summary = section_summary, adjustment_audit = adjustment_audit, checks = checks ) )
}
tryCatch( main(), error = function(e) {
    write_fatal_report( conditionMessage( e ) )
  }
)
