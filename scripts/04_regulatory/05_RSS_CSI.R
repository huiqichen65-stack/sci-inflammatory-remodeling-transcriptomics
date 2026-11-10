rm(list = ls())
gc()
options( stringsAsFactors = FALSE, width = 240, warn = 1, timeout = 86400 )
set.seed(20260802)
PROJECT_ROOT <- Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
EXPECTED_CELLS <- 4203L
EXPECTED_REGULONS <- 292L
EXPECTED_PRIMARY_REGULONS <- 215L
EXPECTED_INFORMATIVE_REGULONS <- 284L
STATE_ORDER <- c( "Shared_prebranch", "Arm_1", "Arm_2" )
EXPECTED_STATE_COUNTS <- c( Shared_prebranch = 116L, Arm_1 = 2289L, Arm_2 = 1798L )
PRIMARY_MIN_TARGETS <- 10L
RSS_TOP_N_PER_STATE <- 20L
RSS_SENSITIVITY_TOP_N <- 20L
MIN_SAMPLE_BRANCH_CELLS <- 20L
REQUIRED_7DPI_SAMPLES <- 2L
CSI_K_CANDIDATES <- 3:12
CSI_MIN_MODULE_SIZE <- 3L
CSI_FIXED_K_SENSITIVITY <- 10L
KEY_TFS <- unique( c( "Spi1", "Irf8", "Atf3", "Fos", "Junb", "Egr1", "Cebpb", "Mafb", "Nfkb1", "Nfkb2", "Rela", "Relb", "Stat1", "Stat3", "Irf1", "Irf7", "Tfeb", "Srebf1", "Nfe2l2", "Bhlhe41", "Mitf", "Nr1h3", "Pparg" ) )
STEP7A_INPUT_DIR <- file.path( PROJECT_ROOT, "04_results", "step7A_GSE162610_pySCENIC_input_environment_audit", "input" )
METADATA_FILE <- file.path( STEP7A_INPUT_DIR, "GSE162610_step7A_frozen_4203_cell_metadata.csv" )
STEP7C_DIR <- file.path( PROJECT_ROOT, "04_results", "step7C_GSE162610_cisTarget_motif_pruning_regulon_inference" )
REGULON_SUMMARY_FILE <- file.path( STEP7C_DIR, "tables", "04_regulon_summary.tsv" )
STEP7D_DIR <- file.path( PROJECT_ROOT, "04_results", "step7D_GSE162610_AUCell_regulon_activity_scoring" )
AUC_RDS_FILE <- file.path( STEP7D_DIR, "RDS", "GSE162610_step7D_AUC_regulons_by_cells.rds" )
AUC_DISTRIBUTION_FILE <- file.path( STEP7D_DIR, "tables", "04_regulon_AUC_distribution_summary.tsv" )
STEP7D_FINAL_DECISION_FILE <- file.path( STEP7D_DIR, "tables", "08_final_decision.csv" )
OUTDIR <- file.path( PROJECT_ROOT, "04_results", "step7E_GSE162610_RSS_CSI_regulon_specificity_modules" )
TABDIR <- file.path(OUTDIR, "tables")
FIGDIR <- file.path(OUTDIR, "figures")
RDSDIR <- file.path(OUTDIR, "RDS")
REPORTDIR <- file.path(OUTDIR, "reports")
for (directory in c( OUTDIR, TABDIR, FIGDIR, RDSDIR, REPORTDIR )) {
  dir.create( directory, recursive = TRUE, showWarnings = FALSE )
}
INPUT_AUDIT_FILE <- file.path( TABDIR, "00_input_file_audit.csv" )
METADATA_AUDIT_FILE <- file.path( TABDIR, "01_cell_metadata_audit.csv" )
ANALYSIS_SET_FILE <- file.path( TABDIR, "02_regulon_analysis_sets.tsv" )
RSS_PRIMARY_LONG_FILE <- file.path( TABDIR, "03_RSS_primary_long.tsv" )
RSS_PRIMARY_WIDE_FILE <- file.path( TABDIR, "03b_RSS_primary_wide.tsv" )
RSS_SENSITIVITY_LONG_FILE <- file.path( TABDIR, "03c_RSS_sensitivity_informative_long.tsv" )
RSS_SENSITIVITY_AUDIT_FILE <- file.path( TABDIR, "03d_RSS_primary_vs_sensitivity_audit.tsv" )
RSS_TOP_FILE <- file.path( TABDIR, "04_RSS_top_regulons_by_state.tsv" )
SAMPLE_DELTA_FILE <- file.path( TABDIR, "05_sample_Arm1_vs_Arm2_regulon_deltas.tsv" )
SAMPLE_ROBUSTNESS_FILE <- file.path( TABDIR, "05b_sample_branch_direction_robustness.tsv" )
PCC_MATRIX_FILE <- file.path( TABDIR, "06_regulon_Pearson_correlation_matrix.tsv" )
CSI_MATRIX_FILE <- file.path( TABDIR, "06b_CSI_matrix_primary.tsv" )
CSI_PAIR_FILE <- file.path( TABDIR, "06c_CSI_pairs_primary.tsv" )
MODULE_ASSIGNMENT_FILE <- file.path( TABDIR, "07_CSI_module_assignments.tsv" )
MODULE_STATE_FILE <- file.path( TABDIR, "08_CSI_module_state_summary.tsv" )
KEY_TF_FILE <- file.path( TABDIR, "09_key_TF_RSS_CSI_audit.tsv" )
K_SELECTION_FILE <- file.path( TABDIR, "10_CSI_k_selection.tsv" )
MODULE_SUMMARY_FILE <- file.path( TABDIR, "11_CSI_module_summary.tsv" )
FIGURE_AUDIT_FILE <- file.path( TABDIR, "11b_figure_generation_audit.csv" )
FINAL_DECISION_FILE <- file.path( TABDIR, "12_final_decision.csv" )
CSI_MATRIX_RDS <- file.path( RDSDIR, "GSE162610_step7E_CSI_matrix_primary.rds" )
MODULE_SCORE_RDS <- file.path( RDSDIR, "GSE162610_step7E_CSI_module_scores.rds" )
STEP7E_RESULTS_RDS <- file.path( RDSDIR, "GSE162610_step7E_RSS_CSI_results.rds" )
REPORT_FILE <- file.path( REPORTDIR, "Step7E_report.txt" )
ERROR_FILE <- file.path( REPORTDIR, "Step7E_ERROR.txt" )
SESSION_FILE <- file.path( REPORTDIR, "Step7E_sessionInfo.txt" )
RSS_HEATMAP_FILE <- file.path( FIGDIR, "01_RSS_top_regulons_heatmap.pdf" )
RSS_RANK_FILE <- file.path( FIGDIR, "02_RSS_rankings.pdf" )
CSI_HEATMAP_FILE <- file.path( FIGDIR, "03_CSI_primary_heatmap.pdf" )
MODULE_HEATMAP_FILE <- file.path( FIGDIR, "04_CSI_module_state_activity_heatmap.pdf" )
required_packages <- c( "data.table", "cluster" )
missing_packages <- required_packages[ !vapply( required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1) ) ]
if (length(missing_packages) > 0L) {
  stop( paste0( "Required R packages are missing: ", paste( missing_packages, collapse = ", " ) ), call. = FALSE )
}
timestamp <- function() {
  format( Sys.time(), "%Y-%m-%d %H:%M:%S" )
}
msg <- function(...) {
  cat( sprintf( "[%s] %s\n", timestamp(), paste0(..., collapse = "") ) )
  flush.console()
}
safe_path <- function(path) {
  normalizePath( path, winslash = "/", mustWork = FALSE )
}
format_duration <- function(seconds) {
  seconds <- as.numeric(seconds)
  if (!is.finite(seconds)) {
    return(NA_character_)
  }
  if (seconds < 60) {
    return(paste0(round(seconds, 1), " sec"))
  }
  if (seconds < 3600) {
    return(paste0(round(seconds / 60, 1), " min"))
  }
  paste0(round(seconds / 3600, 2), " h")
}
format_bytes <- function(bytes) {
  bytes <- as.numeric(bytes)
  if (!is.finite(bytes)) {
    return(NA_character_)
  }
  units <- c("B", "KB", "MB", "GB", "TB")
  index <- 1L
  while ( bytes >= 1024 && index < length(units) ) {
    bytes <- bytes / 1024
    index <- index + 1L
  }
  paste0( format( round(bytes, 3), nsmall = 3, trim = TRUE ), " ", units[[index]] )
}
write_csv <- function(data, path) {
  utils::write.csv( data, file = path, row.names = FALSE, fileEncoding = "UTF-8", na = "" )
}
write_tsv <- function(data, path) {
  data.table::fwrite( data, file = path, sep = "\t", quote = FALSE, na = "" )
}
abort_step <- function(text) {
  writeLines( c( "======================================================================", "GSE162610 | STEP 7E ERROR", "======================================================================", paste0("Time: ", timestamp()), "", text ), con = ERROR_FILE, useBytes = TRUE )
  stop( text, call. = FALSE )
}
normalize_name <- function(x) {
  gsub( "[^a-z0-9]", "", tolower( as.character(x) ) )
}
find_column <- function( data, candidates, required = TRUE ) {
  normalized_columns <- normalize_name( colnames(data) )
  normalized_candidates <- normalize_name( candidates )
  positions <- match( normalized_candidates, normalized_columns, nomatch = 0L )
  positions <- positions[ positions > 0L ]
  if (length(positions) > 0L) {
    return( colnames(data)[ positions[[1]] ] )
  }
  if (required) {
    stop( paste0( "Could not resolve required metadata column. Candidates: ", paste( candidates, collapse = ", " ) ), call. = FALSE )
  }
  NA_character_
}
standardize_state <- function(x) {
  normalized <- normalize_name( trimws( as.character(x) ) )
  output <- rep( NA_character_, length(normalized) )
  output[ normalized %in% c( "sharedprebranch", "prebranch", "shared" ) ] <- "Shared_prebranch"
  output[ normalized %in% c( "arm1", "branch1" ) ] <- "Arm_1"
  output[ normalized %in% c( "arm2", "branch2" ) ] <- "Arm_2"
  output
}
extract_timepoint_day <- function(x) {
  x <- tolower( trimws( as.character(x) ) )
  output <- rep( NA_real_, length(x) )
  zero_index <- grepl( "uninj|control|baseline|naive|sham", x )
  output[ zero_index ] <- 0
  numeric_text <- sub( "^.*?([0-9]+\\.?[0-9]*).*$", "\\1", x )
  has_numeric <- grepl( "[0-9]", x )
  output[ has_numeric ] <- suppressWarnings( as.numeric( numeric_text[ has_numeric ] ) )
  output
}
matrix_to_table <- function( matrix_object, row_name_column ) {
  output <- data.frame( matrix_object, check.names = FALSE, stringsAsFactors = FALSE )
  output <- cbind( setNames( data.frame( rownames(matrix_object), stringsAsFactors = FALSE ), row_name_column ), output )
  rownames(output) <- NULL
  output
}
row_second_largest <- function(matrix_object) {
  apply( matrix_object, 1L, function(values) {
      values <- sort( values, decreasing = TRUE, na.last = NA )
      if (length(values) < 2L) {
        return(NA_real_)
      }
      values[[2]]
    }
  )
}
safe_row_sd <- function(matrix_object) {
  apply( matrix_object, 1L, stats::sd, na.rm = TRUE )
}
calculate_js_distance <- function( p, q ) {
  p <- as.numeric(p)
  q <- as.numeric(q)
  if ( any(!is.finite(p)) || any(!is.finite(q)) || any(p < 0) || any(q < 0) ) {
    stop( "Jensen-Shannon inputs must be finite and non-negative.", call. = FALSE )
  }
  p_sum <- sum(p)
  q_sum <- sum(q)
  if ( p_sum <= 0 || q_sum <= 0 ) {
    return(NA_real_)
  }
  p <- p / p_sum
  q <- q / q_sum
  midpoint <- 0.5 * ( p + q )
  kl_divergence <- function( distribution, reference ) {
    keep <- distribution > 0
    sum( distribution[keep] * log( distribution[keep] / reference[keep] ) )
  }
  js_divergence <- 0.5 * ( kl_divergence( p, midpoint ) + kl_divergence( q, midpoint ) )
  sqrt( max( js_divergence, 0 ) )
}
calculate_rss <- function( activity_regulons_by_cells, labels, state_order ) {
  if ( ncol(activity_regulons_by_cells) != length(labels) ) {
    stop( "RSS activity matrix and labels have different cell counts.", call. = FALSE )
  }
  rss <- matrix( NA_real_, nrow = nrow(activity_regulons_by_cells), ncol = length(state_order), dimnames = list( rownames(activity_regulons_by_cells), state_order ) )
  for (state in state_order) {
    ideal_state <- as.numeric( labels == state )
    for ( regulon_index in seq_len( nrow(activity_regulons_by_cells) ) ) {
      rss[ regulon_index, state ] <- 1 - calculate_js_distance( activity_regulons_by_cells[ regulon_index, , drop = TRUE ], ideal_state )
    }
  }
  rss
}
calculate_csi <- function( correlation_matrix ) {
  number_regulons <- nrow( correlation_matrix )
  csi <- matrix( 1, nrow = number_regulons, ncol = number_regulons, dimnames = dimnames( correlation_matrix ) )
  if (number_regulons < 3L) {
    return(csi)
  }
  all_indices <- seq_len( number_regulons )
  for ( first_index in seq_len( number_regulons - 1L ) ) {
    for ( second_index in seq.int( first_index + 1L, number_regulons ) ) {
      pair_correlation <- correlation_matrix[ first_index, second_index ]
      other_indices <- all_indices[ all_indices != first_index & all_indices != second_index ]
      neither_higher <- correlation_matrix[ first_index, other_indices ] <= pair_correlation + 1e-12 & correlation_matrix[ second_index, other_indices ] <= pair_correlation + 1e-12
      csi_value <- mean( neither_higher )
      csi[ first_index, second_index ] <- csi_value
      csi[ second_index, first_index ] <- csi_value
    }
  }
  csi
}
relabel_modules <- function( raw_clusters, dendrogram_order ) {
  ordered_regulons <- names( raw_clusters )[ dendrogram_order ]
  cluster_order <- unique( raw_clusters[ ordered_regulons ] )
  mapping <- setNames( paste0( "M", seq_along(cluster_order) ), as.character(cluster_order) )
  output <- unname( mapping[ as.character(raw_clusters) ] )
  names(output) <- names( raw_clusters )
  output
}
mean_upper_triangle <- function(matrix_object) {
  if (nrow(matrix_object) < 2L) {
    return(NA_real_)
  }
  mean( matrix_object[ upper.tri(matrix_object) ], na.rm = TRUE )
}
collapse_top <- function( values, n = 10L ) {
  values <- values[ !is.na(values) & nzchar(values) ]
  paste( utils::head( values, n ), collapse = " | " )
}
overall_start <- Sys.time()
msg("==============================================================")
msg("GSE162610 Step 7E | RSS and CSI regulon modules")
msg("==============================================================")
msg("Primary analysis: 215 large informative regulons.")
msg("Sensitivity analysis: 284 informative regulons.")
required_files <- c( METADATA_FILE, REGULON_SUMMARY_FILE, AUC_RDS_FILE, AUC_DISTRIBUTION_FILE, STEP7D_FINAL_DECISION_FILE )
missing_files <- required_files[ !file.exists(required_files) ]
if (length(missing_files) > 0L) {
  abort_step( paste0( "Required files are missing:\n", paste( safe_path(missing_files), collapse = "\n" ) ) )
}
input_audit <- data.frame( role = c( "Step7A_frozen_cell_metadata", "Step7C_regulon_summary", "Step7D_AUC_matrix_RDS", "Step7D_AUC_distribution", "Step7D_final_decision" ), path = safe_path(required_files), size_bytes = file.info(required_files)$size, size_human = vapply( file.info(required_files)$size, format_bytes, FUN.VALUE = character(1) ), exists = file.exists(required_files), stringsAsFactors = FALSE )
write_csv( input_audit, INPUT_AUDIT_FILE )
step7d_decision <- data.table::fread( STEP7D_FINAL_DECISION_FILE, data.table = FALSE, showProgress = FALSE )
step7d_pass <- any( step7d_decision$domain == "Overall Step 7D" & tolower( as.character( step7d_decision$status ) ) == "pass" )
if (!step7d_pass) {
  abort_step( "Step 7D does not contain an Overall PASS decision." )
}
msg("Loading Step 7D AUCell matrix...")
auc_all <- readRDS( AUC_RDS_FILE )
if ( !is.matrix(auc_all) || nrow(auc_all) != EXPECTED_REGULONS || ncol(auc_all) != EXPECTED_CELLS ) {
  abort_step( paste0( "Unexpected AUC dimensions: ", paste( dim(auc_all), collapse = " x " ) ) )
}
storage.mode(auc_all) <- "double"
if ( is.null(rownames(auc_all)) || is.null(colnames(auc_all)) || anyDuplicated(rownames(auc_all)) > 0L || anyDuplicated(colnames(auc_all)) > 0L || any(!is.finite(auc_all)) || any(auc_all < 0) ) {
  abort_step( "The AUCell matrix failed names, duplication, finite-value or non-negative checks." )
}
auc_distribution <- data.table::fread( AUC_DISTRIBUTION_FILE, data.table = FALSE, showProgress = FALSE )
required_distribution_columns <- c( "regulon", "TF", "n_targets", "low_information_flag" )
missing_distribution_columns <- setdiff( required_distribution_columns, colnames(auc_distribution) )
if (length(missing_distribution_columns) > 0L) {
  abort_step( paste0( "Step 7D distribution table is missing: ", paste( missing_distribution_columns, collapse = ", " ) ) )
}
auc_distribution$regulon <- as.character( auc_distribution$regulon )
auc_distribution$TF <- as.character( auc_distribution$TF )
auc_distribution$n_targets <- as.integer( auc_distribution$n_targets )
auc_distribution$low_information_flag <- as.logical( auc_distribution$low_information_flag )
if ( !setequal( rownames(auc_all), auc_distribution$regulon ) ) {
  abort_step( "Regulon names differ between the AUC matrix and Step 7D information table." )
}
auc_distribution <- auc_distribution[ match( rownames(auc_all), auc_distribution$regulon ), , drop = FALSE ]
primary_regulons <- auc_distribution$regulon[ auc_distribution$n_targets >= PRIMARY_MIN_TARGETS & !auc_distribution$low_information_flag ]
informative_regulons <- auc_distribution$regulon[ !auc_distribution$low_information_flag ]
if ( length(primary_regulons) != EXPECTED_PRIMARY_REGULONS || length(informative_regulons) != EXPECTED_INFORMATIVE_REGULONS ) {
  abort_step( paste0( "Unexpected analysis-set sizes. Primary=", length(primary_regulons), "; informative=", length(informative_regulons), "." ) )
}
analysis_sets <- auc_distribution
analysis_sets$in_primary_215 <- analysis_sets$regulon %in%
  primary_regulons
analysis_sets$in_informative_284 <- analysis_sets$regulon %in%
  informative_regulons
write_tsv( analysis_sets, ANALYSIS_SET_FILE )
auc_primary <- auc_all[ primary_regulons, , drop = FALSE ]
auc_informative <- auc_all[ informative_regulons, , drop = FALSE ]
msg("Loading frozen cell metadata...")
msg("Expected Step 7A columns: cell, step5CR_state, sample_id, timepoint, retained_step5BR_pseudotime, Y29_signed_progress.")
metadata_raw <- data.table::fread( METADATA_FILE, data.table = FALSE, check.names = FALSE, showProgress = FALSE )
cell_id_column <- find_column( metadata_raw, c( "CellID", "cell_id", "cell", "barcode", "barcodes", "rownames" ), required = FALSE )
if (is.na(cell_id_column)) {
  overlap_fraction <- vapply( metadata_raw, function(values) {
      mean( as.character(values) %in% colnames(auc_all) )
    },
    FUN.VALUE = numeric(1)
  )
  best_column <- which.max( overlap_fraction )
  if ( length(best_column) == 0L || overlap_fraction[ best_column ] < 0.90 ) {
    abort_step( "Could not identify the cell-ID metadata column." )
  }
  cell_id_column <- colnames( metadata_raw )[ best_column ]
}
state_column <- find_column( metadata_raw, c( "step5CR_state", "step5cr_state", "BranchState", "branch_state", "state", "branch", "trajectory_state" ) )
sample_column <- find_column( metadata_raw, c( "SampleID", "sample_id", "sample", "orig.ident", "orig_ident", "replicate" ) )
timepoint_column <- find_column( metadata_raw, c( "Timepoint", "time_point", "time", "dpi", "day" ) )
pseudotime_column <- find_column( metadata_raw, c( "retained_step5BR_pseudotime", "retained_step5br_pseudotime", "Pseudotime", "retained_pseudotime", "injury_pseudotime", "monocle3_pseudotime" ), required = FALSE )
y29_column <- find_column( metadata_raw, c( "Y29SignedProgress", "Y29_signed_progress", "Y_29_signed_progress", "Y29Progress", "Y_29_progress" ), required = FALSE )
metadata <- data.frame( CellID = as.character( metadata_raw[[ cell_id_column ]] ), BranchState = standardize_state( metadata_raw[[ state_column ]] ), SampleID = as.character( metadata_raw[[ sample_column ]] ), Timepoint = as.character( metadata_raw[[ timepoint_column ]] ), stringsAsFactors = FALSE )
metadata$TimepointDay <- extract_timepoint_day( metadata$Timepoint )
metadata$Pseudotime <- if ( !is.na(pseudotime_column) ) {
  suppressWarnings( as.numeric( metadata_raw[[ pseudotime_column ]] ) )
} else {
  NA_real_
}
metadata$Y29SignedProgress <- if ( !is.na(y29_column) ) {
  suppressWarnings( as.numeric( metadata_raw[[ y29_column ]] ) )
} else {
  NA_real_
}
if ( nrow(metadata) != EXPECTED_CELLS || anyDuplicated(metadata$CellID) > 0L || any(is.na(metadata$BranchState)) ) {
  abort_step( "Metadata failed cell-count, duplicate-ID or state checks." )
}
metadata_index <- match( colnames(auc_all), metadata$CellID )
if (anyNA(metadata_index)) {
  abort_step( "Some AUCell cells are absent from metadata." )
}
metadata <- metadata[ metadata_index, , drop = FALSE ]
if ( !identical( metadata$CellID, colnames(auc_all) ) ) {
  abort_step( "Metadata could not be aligned exactly to the AUCell matrix." )
}
metadata$BranchState <- factor( metadata$BranchState, levels = STATE_ORDER )
observed_state_counts <- table( metadata$BranchState )
state_count_pass <- all( as.integer( observed_state_counts[ names(EXPECTED_STATE_COUNTS) ] ) == as.integer( EXPECTED_STATE_COUNTS ) )
if (!state_count_pass) {
  abort_step( paste0( "Unexpected state counts: ", paste( names(observed_state_counts), observed_state_counts, sep = "=", collapse = "; " ) ) )
}
metadata_audit <- data.frame( metric = c( "metadata_rows", "unique_cell_IDs", "exact_AUC_cell_order", "Shared_prebranch_cells", "Arm_1_cells", "Arm_2_cells", "unique_samples", "unique_timepoint_labels", "seven_dpi_cells", "finite_pseudotime_cells", "finite_Y29SignedProgress_cells", "resolved_cell_ID_column", "resolved_state_column", "resolved_sample_column", "resolved_timepoint_column", "resolved_pseudotime_column", "resolved_Y29_column" ), value = c( nrow(metadata), length(unique(metadata$CellID)), identical( metadata$CellID, colnames(auc_all) ), observed_state_counts[ "Shared_prebranch" ], observed_state_counts[ "Arm_1" ], observed_state_counts[ "Arm_2" ], length(unique(metadata$SampleID)), length(unique(metadata$Timepoint)), sum( metadata$TimepointDay == 7, na.rm = TRUE ), sum(is.finite(metadata$Pseudotime)), sum(is.finite(metadata$Y29SignedProgress)), cell_id_column, state_column, sample_column, timepoint_column, pseudotime_column, y29_column ), stringsAsFactors = FALSE )
write_csv( metadata_audit, METADATA_AUDIT_FILE )
rss_start <- Sys.time()
msg("Calculating RSS for 215 primary regulons...")
synthetic_activity <- matrix( c( 1, 1, 0, 0 ), nrow = 1L, dimnames = list( "synthetic", paste0( "cell", 1:4 ) ) )
synthetic_labels <- c( "A", "A", "B", "B" )
synthetic_rss <- calculate_rss( synthetic_activity, synthetic_labels, c( "A", "B" ) )
synthetic_rss_pass <- abs( synthetic_rss[ "synthetic", "A" ] - 1 ) < 1e-12
rss_primary <- calculate_rss( auc_primary, as.character( metadata$BranchState ), STATE_ORDER )
rss_informative <- calculate_rss( auc_informative, as.character( metadata$BranchState ), STATE_ORDER )
rss_primary_pass <- all(is.finite(rss_primary)) && all( rss_primary >= 0 & rss_primary <= 1 )
rss_sensitivity_pass <- all(is.finite(rss_informative)) && all( rss_informative >= 0 & rss_informative <= 1 )
if ( !synthetic_rss_pass || !rss_primary_pass || !rss_sensitivity_pass ) {
  abort_step( "RSS failed the synthetic control or finite/range checks." )
}
primary_info <- auc_distribution[ match( rownames(rss_primary), auc_distribution$regulon ), , drop = FALSE ]
informative_info <- auc_distribution[ match( rownames(rss_informative), auc_distribution$regulon ), , drop = FALSE ]
preferred_index <- max.col( rss_primary, ties.method = "first" )
preferred_state <- colnames( rss_primary )[ preferred_index ]
preferred_rss <- rss_primary[ cbind( seq_len(nrow(rss_primary)), preferred_index ) ]
second_rss <- row_second_largest( rss_primary )
rss_margin <- preferred_rss - second_rss
state_mean_auc <- sapply( STATE_ORDER, function(state) {
    rowMeans( auc_primary[ , metadata$BranchState == state, drop = FALSE ] )
  }
)
state_median_auc <- sapply( STATE_ORDER, function(state) {
    apply( auc_primary[ , metadata$BranchState == state, drop = FALSE ], 1L, stats::median )
  }
)
colnames(state_mean_auc) <- STATE_ORDER
colnames(state_median_auc) <- STATE_ORDER
arm1_auc <- auc_primary[ , metadata$BranchState == "Arm_1", drop = FALSE ]
arm2_auc <- auc_primary[ , metadata$BranchState == "Arm_2", drop = FALSE ]
arm1_mean <- rowMeans(arm1_auc)
arm2_mean <- rowMeans(arm2_auc)
arm_delta <- arm1_mean - arm2_mean
arm1_sd <- safe_row_sd(arm1_auc)
arm2_sd <- safe_row_sd(arm2_auc)
pooled_sd <- sqrt( ( ( ncol(arm1_auc) - 1 ) * arm1_sd^2 + ( ncol(arm2_auc) - 1 ) * arm2_sd^2 ) / ( ncol(arm1_auc) + ncol(arm2_auc) - 2 ) )
arm_standardized_difference <- ifelse( pooled_sd > 0, arm_delta / pooled_sd, NA_real_ )
rss_primary_wide <- data.frame( regulon = rownames(rss_primary), TF = primary_info$TF, n_targets = primary_info$n_targets, Shared_prebranch_RSS = rss_primary[ , "Shared_prebranch" ], Arm_1_RSS = rss_primary[ , "Arm_1" ], Arm_2_RSS = rss_primary[ , "Arm_2" ], preferred_state = preferred_state, preferred_RSS = preferred_rss, second_highest_RSS = second_rss, RSS_margin = rss_margin, Shared_prebranch_mean_AUC = state_mean_auc[ , "Shared_prebranch" ], Arm_1_mean_AUC = arm1_mean, Arm_2_mean_AUC = arm2_mean, Arm1_minus_Arm2_mean_AUC = arm_delta, Arm1_vs_Arm2_standardized_difference = arm_standardized_difference, stringsAsFactors = FALSE )
write_tsv( rss_primary_wide, RSS_PRIMARY_WIDE_FILE )
rss_primary_long <- do.call( rbind, lapply( STATE_ORDER, function(state) {
      other_states <- setdiff( STATE_ORDER, state )
      data.frame( regulon = rownames(rss_primary), TF = primary_info$TF, n_targets = primary_info$n_targets, state = state, RSS = rss_primary[ , state ], RSS_rank_in_state = rank( -rss_primary[ , state ], ties.method = "min" ), mean_AUC_in_state = state_mean_auc[ , state ], median_AUC_in_state = state_median_auc[ , state ], mean_AUC_other_states = rowMeans( state_mean_auc[ , other_states, drop = FALSE ] ), preferred_state = preferred_state, preferred_RSS = preferred_rss, second_highest_RSS = second_rss, RSS_margin = rss_margin, Arm1_minus_Arm2_mean_AUC = arm_delta, Arm1_vs_Arm2_standardized_difference = arm_standardized_difference, stringsAsFactors = FALSE )
    }
  )
)
rss_primary_long <- rss_primary_long[ order( match( rss_primary_long$state, STATE_ORDER ), rss_primary_long$RSS_rank_in_state, -rss_primary_long$RSS ), , drop = FALSE ]
write_tsv( rss_primary_long, RSS_PRIMARY_LONG_FILE )
rss_sensitivity_long <- do.call( rbind, lapply( STATE_ORDER, function(state) {
      data.frame( regulon = rownames(rss_informative), TF = informative_info$TF, n_targets = informative_info$n_targets, state = state, RSS = rss_informative[ , state ], RSS_rank_in_state = rank( -rss_informative[ , state ], ties.method = "min" ), stringsAsFactors = FALSE )
    }
  )
)
write_tsv( rss_sensitivity_long, RSS_SENSITIVITY_LONG_FILE )
rss_sensitivity_audit <- do.call( rbind, lapply( STATE_ORDER, function(state) {
      primary_rank <- rank( -rss_primary[ , state ], ties.method = "average" )
      sensitivity_rank <- rank( -rss_informative[ , state ], ties.method = "average" )
      sensitivity_rank_on_primary <- sensitivity_rank[ match( rownames(rss_primary), rownames(rss_informative) ) ]
      primary_top <- rownames(rss_primary)[ order( rss_primary[ , state ], decreasing = TRUE ) ][ seq_len( min( RSS_SENSITIVITY_TOP_N, nrow(rss_primary) ) ) ]
      sensitivity_top <- rownames(rss_informative)[ order( rss_informative[ , state ], decreasing = TRUE ) ][ seq_len( min( RSS_SENSITIVITY_TOP_N, nrow(rss_informative) ) ) ]
      data.frame( state = state, primary_regulons = nrow(rss_primary), sensitivity_regulons = nrow(rss_informative), Spearman_rank_correlation_on_primary = suppressWarnings( stats::cor( primary_rank, sensitivity_rank_on_primary, method = "spearman" ) ), top20_overlap_n = length( intersect( primary_top, sensitivity_top ) ), top20_overlap_fraction = length( intersect( primary_top, sensitivity_top ) ) / RSS_SENSITIVITY_TOP_N, stringsAsFactors = FALSE )
    }
  )
)
write_tsv( rss_sensitivity_audit, RSS_SENSITIVITY_AUDIT_FILE )
rss_elapsed <- as.numeric( difftime( Sys.time(), rss_start, units = "secs" ) )
msg( "RSS PASS in ", format_duration(rss_elapsed), "." )
sample_start <- Sys.time()
msg("Calculating sample-level Arm_1 versus Arm_2 directions...")
metadata_dt <- data.table::as.data.table( metadata )
sample_branch_counts <- metadata_dt[ BranchState %in% c( "Arm_1", "Arm_2" ), .N, by = .( SampleID, Timepoint, TimepointDay, BranchState ) ]
sample_branch_count_wide <- data.table::dcast( sample_branch_counts, SampleID + Timepoint + TimepointDay ~ BranchState, value.var = "N", fill = 0 )
if (!"Arm_1" %in% colnames(sample_branch_count_wide)) {
  sample_branch_count_wide$Arm_1 <- 0L
}
if (!"Arm_2" %in% colnames(sample_branch_count_wide)) {
  sample_branch_count_wide$Arm_2 <- 0L
}
evaluable_samples <- sample_branch_count_wide[ Arm_1 >= MIN_SAMPLE_BRANCH_CELLS & Arm_2 >= MIN_SAMPLE_BRANCH_CELLS ]
sample_delta_parts <- list()
if (nrow(evaluable_samples) > 0L) {
  for ( sample_index in seq_len(nrow(evaluable_samples)) ) {
    sample_id <- as.character( evaluable_samples$SampleID[ sample_index ] )
    sample_timepoint <- as.character( evaluable_samples$Timepoint[ sample_index ] )
    sample_timepoint_day <- as.numeric( evaluable_samples$TimepointDay[ sample_index ] )
    arm1_cells <- metadata$SampleID == sample_id & metadata$BranchState == "Arm_1"
    arm2_cells <- metadata$SampleID == sample_id & metadata$BranchState == "Arm_2"
    current <- data.frame( SampleID = sample_id, Timepoint = sample_timepoint, TimepointDay = sample_timepoint_day, n_Arm_1 = sum(arm1_cells), n_Arm_2 = sum(arm2_cells), regulon = rownames(auc_primary), TF = primary_info$TF, mean_AUC_Arm_1 = rowMeans( auc_primary[ , arm1_cells, drop = FALSE ] ), mean_AUC_Arm_2 = rowMeans( auc_primary[ , arm2_cells, drop = FALSE ] ), stringsAsFactors = FALSE )
    current$Arm1_minus_Arm2 <- current$mean_AUC_Arm_1 - current$mean_AUC_Arm_2
    sample_delta_parts[[ length(sample_delta_parts) + 1L ]] <- current
  }
}
if (length(sample_delta_parts) > 0L) {
  sample_delta <- data.table::rbindlist( sample_delta_parts, use.names = TRUE, fill = TRUE )
} else {
  sample_delta <- data.table::data.table( SampleID = character(0), Timepoint = character(0), TimepointDay = numeric(0), n_Arm_1 = integer(0), n_Arm_2 = integer(0), regulon = character(0), TF = character(0), mean_AUC_Arm_1 = numeric(0), mean_AUC_Arm_2 = numeric(0), Arm1_minus_Arm2 = numeric(0) )
}
overall_delta_lookup <- setNames( arm_delta, rownames(auc_primary) )
if (nrow(sample_delta) > 0L) {
  sample_delta[ , overall_Arm1_minus_Arm2 := overall_delta_lookup[ regulon ] ]
  sample_delta[ , same_direction_as_overall := sign(Arm1_minus_Arm2) == sign(overall_Arm1_minus_Arm2) & sign(overall_Arm1_minus_Arm2) != 0 ]
  sample_delta[ , is_7dpi := is.finite(TimepointDay) & abs(TimepointDay - 7) < 1e-8 ]
} else {
  sample_delta[ , overall_Arm1_minus_Arm2 := numeric(0) ]
  sample_delta[ , same_direction_as_overall := logical(0) ]
  sample_delta[ , is_7dpi := logical(0) ]
}
write_tsv( sample_delta, SAMPLE_DELTA_FILE )
if (nrow(sample_delta) > 0L) {
  sample_robustness <- sample_delta[ , .( TF = TF[[1]], overall_Arm1_minus_Arm2 = overall_Arm1_minus_Arm2[[1]], overall_preferred_arm = ifelse( overall_Arm1_minus_Arm2[[1]] > 0, "Arm_1", ifelse( overall_Arm1_minus_Arm2[[1]] < 0, "Arm_2", "Tie" ) ), n_evaluable_samples = .N, n_samples_same_direction = sum( same_direction_as_overall, na.rm = TRUE ), sample_direction_consistency_fraction = mean( same_direction_as_overall, na.rm = TRUE ), median_sample_delta = stats::median( Arm1_minus_Arm2, na.rm = TRUE ), minimum_sample_delta = min( Arm1_minus_Arm2, na.rm = TRUE ), maximum_sample_delta = max( Arm1_minus_Arm2, na.rm = TRUE ), n_7dpi_samples = sum( is_7dpi, na.rm = TRUE ), n_7dpi_samples_same_direction = sum( is_7dpi & same_direction_as_overall, na.rm = TRUE ), seven_dpi_direction_consistent = sum( is_7dpi, na.rm = TRUE ) >= REQUIRED_7DPI_SAMPLES & sum( is_7dpi & same_direction_as_overall, na.rm = TRUE ) == sum( is_7dpi, na.rm = TRUE ) ), by = regulon ]
} else {
  sample_robustness <- data.table::data.table( regulon = rownames(auc_primary), TF = primary_info$TF, overall_Arm1_minus_Arm2 = arm_delta, overall_preferred_arm = ifelse( arm_delta > 0, "Arm_1", ifelse( arm_delta < 0, "Arm_2", "Tie" ) ), n_evaluable_samples = 0L, n_samples_same_direction = 0L, sample_direction_consistency_fraction = NA_real_, median_sample_delta = NA_real_, minimum_sample_delta = NA_real_, maximum_sample_delta = NA_real_, n_7dpi_samples = 0L, n_7dpi_samples_same_direction = 0L, seven_dpi_direction_consistent = FALSE )
}
write_tsv( sample_robustness, SAMPLE_ROBUSTNESS_FILE )
sample_robustness_df <- as.data.frame( sample_robustness )
rownames(sample_robustness_df) <- sample_robustness_df$regulon
rss_top <- do.call( rbind, lapply( STATE_ORDER, function(state) {
      current <- rss_primary_long[ rss_primary_long$state == state, , drop = FALSE ]
      current <- current[ order( current$RSS_rank_in_state, -current$RSS ), , drop = FALSE ]
      current <- utils::head( current, RSS_TOP_N_PER_STATE )
      current$is_preferred_state <- current$preferred_state == state
      if ( state %in% c( "Arm_1", "Arm_2" ) ) {
        current$seven_dpi_direction_consistent <- sample_robustness_df[ current$regulon, "seven_dpi_direction_consistent" ]
        current$sample_direction_consistency_fraction <- sample_robustness_df[ current$regulon, "sample_direction_consistency_fraction" ]
      } else {
        current$seven_dpi_direction_consistent <- NA
        current$sample_direction_consistency_fraction <- NA_real_
      }
      current
    }
  )
)
write_tsv( rss_top, RSS_TOP_FILE )
observed_7dpi_samples <- length( unique( sample_delta$SampleID[ sample_delta$is_7dpi ] ) )
sample_audit_pass <- observed_7dpi_samples >= REQUIRED_7DPI_SAMPLES
sample_elapsed <- as.numeric( difftime( Sys.time(), sample_start, units = "secs" ) )
msg( "Sample audit complete in ", format_duration(sample_elapsed), "; evaluable 7 dpi samples=", observed_7dpi_samples, "." )
csi_start <- Sys.time()
msg("Calculating Pearson correlation matrix...")
regulon_correlation <- stats::cor( t(auc_primary), method = "pearson", use = "pairwise.complete.obs" )
if ( any(!is.finite(regulon_correlation)) || max( abs( regulon_correlation - t(regulon_correlation) ) ) > 1e-10 ) {
  abort_step( "Regulon correlation matrix is non-finite or asymmetric." )
}
diag(regulon_correlation) <- 1
msg("Calculating CSI matrix...")
csi_matrix <- calculate_csi( regulon_correlation )
csi_pass <- all(is.finite(csi_matrix)) && all( csi_matrix >= 0 & csi_matrix <= 1 ) && max( abs( csi_matrix - t(csi_matrix) ) ) <
    1e-12 && max( abs( diag(csi_matrix) - 1 ) ) <
    1e-12
if (!csi_pass) {
  abort_step( "CSI matrix failed finite, range, symmetry or diagonal checks." )
}
saveRDS( csi_matrix, CSI_MATRIX_RDS, compress = TRUE )
write_tsv( matrix_to_table( regulon_correlation, "regulon" ), PCC_MATRIX_FILE )
write_tsv( matrix_to_table( csi_matrix, "regulon" ), CSI_MATRIX_FILE )
upper_index <- which( upper.tri(csi_matrix), arr.ind = TRUE )
csi_pairs <- data.frame( regulon_1 = rownames(csi_matrix)[ upper_index[ , 1 ] ], regulon_2 = colnames(csi_matrix)[ upper_index[ , 2 ] ], Pearson_correlation = regulon_correlation[ upper_index ], CSI = csi_matrix[ upper_index ], stringsAsFactors = FALSE )
csi_pairs <- csi_pairs[ order( -csi_pairs$CSI, -csi_pairs$Pearson_correlation ), , drop = FALSE ]
write_tsv( csi_pairs, CSI_PAIR_FILE )
msg("Selecting CSI module count...")
csi_distance <- stats::dist( csi_matrix, method = "euclidean" )
csi_hclust <- stats::hclust( csi_distance, method = "average" )
valid_k <- CSI_K_CANDIDATES[ CSI_K_CANDIDATES >= 2L & CSI_K_CANDIDATES < nrow(csi_matrix) ]
k_rows <- vector( "list", length(valid_k) )
for (index in seq_along(valid_k)) {
  current_k <- valid_k[[index]]
  clusters <- stats::cutree( csi_hclust, k = current_k )
  silhouette_object <- cluster::silhouette( clusters, csi_distance )
  cluster_sizes <- table( clusters )
  k_rows[[index]] <- data.frame( k = current_k, average_silhouette_width = mean( silhouette_object[ , "sil_width" ] ), median_silhouette_width = stats::median( silhouette_object[ , "sil_width" ] ), minimum_silhouette_width = min( silhouette_object[ , "sil_width" ] ), minimum_module_size = min( cluster_sizes ), maximum_module_size = max( cluster_sizes ), singleton_modules = sum( cluster_sizes == 1L ), minimum_size_constraint_pass = min(cluster_sizes) >= CSI_MIN_MODULE_SIZE, stringsAsFactors = FALSE )
}
k_selection <- do.call( rbind, k_rows )
eligible_k <- k_selection[ k_selection$minimum_size_constraint_pass, , drop = FALSE ]
if (nrow(eligible_k) > 0L) {
  eligible_k <- eligible_k[ order( -eligible_k$average_silhouette_width, eligible_k$k ), , drop = FALSE ]
  chosen_k <- eligible_k$k[[1]]
  chosen_k_reason <- paste0( "Highest average silhouette among k with minimum module size >= ", CSI_MIN_MODULE_SIZE, "." )
} else {
  fallback_k <- k_selection[ order( -k_selection$average_silhouette_width, k_selection$k ), , drop = FALSE ]
  chosen_k <- fallback_k$k[[1]]
  chosen_k_reason <- paste0( "Fallback highest average silhouette; no candidate met minimum module size ", CSI_MIN_MODULE_SIZE, "." )
}
k_selection$chosen_primary_k <- k_selection$k == chosen_k
write_tsv( k_selection, K_SELECTION_FILE )
primary_raw_clusters <- stats::cutree( csi_hclust, k = chosen_k )
primary_modules <- relabel_modules( primary_raw_clusters, csi_hclust$order )
fixed_k <- min( CSI_FIXED_K_SENSITIVITY, nrow(csi_matrix) - 1L )
fixed_raw_clusters <- stats::cutree( csi_hclust, k = fixed_k )
fixed_modules <- relabel_modules( fixed_raw_clusters, csi_hclust$order )
module_assignment <- rss_primary_wide[ match( names(primary_modules), rss_primary_wide$regulon ), , drop = FALSE ]
module_assignment$CSI_module <- primary_modules[ module_assignment$regulon ]
module_assignment$CSI_module_fixed_k10 <- fixed_modules[ module_assignment$regulon ]
module_assignment$seven_dpi_direction_consistent <- sample_robustness_df[ module_assignment$regulon, "seven_dpi_direction_consistent" ]
module_assignment$sample_direction_consistency_fraction <- sample_robustness_df[ module_assignment$regulon, "sample_direction_consistency_fraction" ]
module_assignment <- module_assignment[ order( as.integer( sub( "^M", "", module_assignment$CSI_module ) ), -module_assignment$preferred_RSS ), , drop = FALSE ]
write_tsv( module_assignment, MODULE_ASSIGNMENT_FILE )
csi_elapsed <- as.numeric( difftime( Sys.time(), csi_start, units = "secs" ) )
msg( "CSI PASS in ", format_duration(csi_elapsed), "; selected k=", chosen_k, "." )
module_start <- Sys.time()
msg("Calculating CSI module activity...")
auc_row_mean <- rowMeans( auc_primary )
auc_row_sd <- safe_row_sd( auc_primary )
auc_z <- sweep( auc_primary, 1L, auc_row_mean, "-" )
auc_z <- sweep( auc_z, 1L, auc_row_sd, "/" )
auc_z[ !is.finite(auc_z) ] <- 0
module_ids <- unique( module_assignment$CSI_module )
module_ids <- module_ids[ order( as.integer( sub( "^M", "", module_ids ) ) ) ]
module_raw_activity <- matrix( NA_real_, nrow = length(module_ids), ncol = EXPECTED_CELLS, dimnames = list( module_ids, colnames(auc_primary) ) )
module_z_activity <- module_raw_activity
for (module_id in module_ids) {
  member_regulons <- module_assignment$regulon[ module_assignment$CSI_module == module_id ]
  module_raw_activity[ module_id, ] <- colMeans( auc_primary[ member_regulons, , drop = FALSE ] )
  module_z_activity[ module_id, ] <- colMeans( auc_z[ member_regulons, , drop = FALSE ] )
}
module_rss <- calculate_rss( module_raw_activity, as.character( metadata$BranchState ), STATE_ORDER )
module_state_parts <- list()
for (module_id in module_ids) {
  for (state in STATE_ORDER) {
    state_cells <- metadata$BranchState == state
    module_state_parts[[ length(module_state_parts) + 1L ]] <- data.frame( CSI_module = module_id, module_size = sum( module_assignment$CSI_module == module_id ), state = state, module_RSS = module_rss[ module_id, state ], module_RSS_rank_in_state = rank( -module_rss[ , state ], ties.method = "min" )[ match( module_id, rownames(module_rss) ) ], mean_module_raw_AUC = mean( module_raw_activity[ module_id, state_cells ] ), median_module_raw_AUC = stats::median( module_raw_activity[ module_id, state_cells ] ), mean_module_z_activity = mean( module_z_activity[ module_id, state_cells ] ), median_module_z_activity = stats::median( module_z_activity[ module_id, state_cells ] ), stringsAsFactors = FALSE )
  }
}
module_state_summary <- do.call( rbind, module_state_parts )
write_tsv( module_state_summary, MODULE_STATE_FILE )
module_summary_parts <- list()
for (module_id in module_ids) {
  members <- module_assignment[ module_assignment$CSI_module == module_id, , drop = FALSE ]
  member_regulons <- members$regulon
  within_csi <- mean_upper_triangle( csi_matrix[ member_regulons, member_regulons, drop = FALSE ] )
  outside_regulons <- setdiff( rownames(csi_matrix), member_regulons )
  between_csi <- if ( length(outside_regulons) > 0L ) {
    mean( csi_matrix[ member_regulons, outside_regulons, drop = FALSE ], na.rm = TRUE )
  } else {
    NA_real_
  }
  preferred_module_state <- colnames( module_rss )[ which.max( module_rss[ module_id, ] ) ]
  preferred_column <- paste0( preferred_module_state, "_RSS" )
  members <- members[ order( -members[[ preferred_column ]], -members$preferred_RSS ), , drop = FALSE ]
  key_tfs <- intersect( members$TF, KEY_TFS )
  module_summary_parts[[ length(module_summary_parts) + 1L ]] <- data.frame( CSI_module = module_id, module_size = nrow(members), preferred_state_by_module_RSS = preferred_module_state, preferred_module_RSS = max( module_rss[ module_id, ] ), Shared_prebranch_module_RSS = module_rss[ module_id, "Shared_prebranch" ], Arm_1_module_RSS = module_rss[ module_id, "Arm_1" ], Arm_2_module_RSS = module_rss[ module_id, "Arm_2" ], Arm1_minus_Arm2_mean_module_z = mean( module_z_activity[ module_id, metadata$BranchState == "Arm_1" ] ) - mean( module_z_activity[ module_id, metadata$BranchState == "Arm_2" ] ), mean_within_module_CSI = within_csi, mean_between_module_CSI = between_csi, within_minus_between_CSI = within_csi - between_csi, mean_regulon_target_count = mean( members$n_targets ), median_regulon_target_count = stats::median( members$n_targets ), number_key_TFs = length(key_tfs), key_TFs = paste( key_tfs, collapse = " | " ), top_regulons_for_preferred_state = collapse_top( members$regulon, 10L ), stringsAsFactors = FALSE )
}
module_summary <- do.call( rbind, module_summary_parts )
write_tsv( module_summary, MODULE_SUMMARY_FILE )
saveRDS( list( raw_mean_AUC_modules_by_cells = module_raw_activity, z_mean_modules_by_cells = module_z_activity, module_RSS = module_rss, module_assignment = module_assignment ), MODULE_SCORE_RDS, compress = TRUE )
module_elapsed <- as.numeric( difftime( Sys.time(), module_start, units = "secs" ) )
msg( "Module activity complete in ", format_duration(module_elapsed), "." )
key_tf_parts <- list()
for (current_tf in KEY_TFS) {
  matching <- module_assignment[ module_assignment$TF == current_tf, , drop = FALSE ]
  if (nrow(matching) == 0L) {
    key_tf_parts[[ length(key_tf_parts) + 1L ]] <- data.frame( TF = current_tf, regulon = "", present_in_primary_215 = FALSE, n_targets = 0L, Shared_prebranch_RSS = NA_real_, Arm_1_RSS = NA_real_, Arm_2_RSS = NA_real_, preferred_state = NA_character_, preferred_RSS = NA_real_, RSS_margin = NA_real_, Arm1_minus_Arm2_mean_AUC = NA_real_, Arm1_vs_Arm2_standardized_difference = NA_real_, seven_dpi_direction_consistent = NA, sample_direction_consistency_fraction = NA_real_, CSI_module = NA_character_, module_preferred_state = NA_character_, module_preferred_RSS = NA_real_, stringsAsFactors = FALSE )
  } else {
    matching <- matching[ 1L, , drop = FALSE ]
    current_module <- matching$CSI_module[ 1L ]
    current_module_summary <- module_summary[ module_summary$CSI_module == current_module, , drop = FALSE ]
    key_tf_parts[[ length(key_tf_parts) + 1L ]] <- data.frame( TF = current_tf, regulon = matching$regulon, present_in_primary_215 = TRUE, n_targets = matching$n_targets, Shared_prebranch_RSS = matching$Shared_prebranch_RSS, Arm_1_RSS = matching$Arm_1_RSS, Arm_2_RSS = matching$Arm_2_RSS, preferred_state = matching$preferred_state, preferred_RSS = matching$preferred_RSS, RSS_margin = matching$RSS_margin, Arm1_minus_Arm2_mean_AUC = matching$Arm1_minus_Arm2_mean_AUC, Arm1_vs_Arm2_standardized_difference = matching$Arm1_vs_Arm2_standardized_difference, seven_dpi_direction_consistent = matching$seven_dpi_direction_consistent, sample_direction_consistency_fraction = matching$sample_direction_consistency_fraction, CSI_module = current_module, module_preferred_state = current_module_summary$preferred_state_by_module_RSS, module_preferred_RSS = current_module_summary$preferred_module_RSS, stringsAsFactors = FALSE )
  }
}
key_tf_audit <- do.call( rbind, key_tf_parts )
write_tsv( key_tf_audit, KEY_TF_FILE )
figure_audit <- data.frame( figure = c( "RSS_top_heatmap", "RSS_rankings", "CSI_heatmap", "module_state_heatmap" ), created = FALSE, path = c( RSS_HEATMAP_FILE, RSS_RANK_FILE, CSI_HEATMAP_FILE, MODULE_HEATMAP_FILE ), stringsAsFactors = FALSE )
try( {
    grDevices::pdf( RSS_RANK_FILE, width = 12, height = 4.5 )
    old_par <- graphics::par( mfrow = c(1, 3), mar = c(4, 4, 3, 1) )
    for (state in STATE_ORDER) {
      ordered <- order( rss_primary[ , state ], decreasing = TRUE )
      values <- rss_primary[ ordered, state ]
      graphics::plot( seq_along(values), values, pch = 16, cex = 0.45, xlab = "Regulon rank", ylab = "RSS", main = state )
      top_indices <- seq_len( min( RSS_TOP_N_PER_STATE, length(values) ) )
      graphics::text( top_indices, values[ top_indices ], labels = rownames(rss_primary)[ ordered[ top_indices ] ], pos = 4, cex = 0.42 )
    }
    graphics::par(old_par)
    grDevices::dev.off()
    figure_audit$created[ figure_audit$figure == "RSS_rankings" ] <- TRUE
  },
  silent = TRUE
)
if ( requireNamespace( "pheatmap", quietly = TRUE ) ) {
  try( {
      top_regulons <- unique( rss_top$regulon )
      pheatmap::pheatmap( rss_primary[ top_regulons, STATE_ORDER, drop = FALSE ], cluster_rows = TRUE, cluster_cols = FALSE, scale = "none", fontsize_row = 5, filename = RSS_HEATMAP_FILE, width = 6, height = max( 7, length(top_regulons) * 0.12 ), main = "Top state-associated regulons by RSS" )
      figure_audit$created[ figure_audit$figure == "RSS_top_heatmap" ] <- TRUE
    },
    silent = TRUE
  )
  try( {
      annotation <- data.frame( CSI_module = module_assignment$CSI_module, row.names = module_assignment$regulon, stringsAsFactors = FALSE )
      pheatmap::pheatmap( csi_matrix, cluster_rows = csi_hclust, cluster_cols = csi_hclust, show_rownames = FALSE, show_colnames = FALSE, annotation_row = annotation, annotation_col = annotation, filename = CSI_HEATMAP_FILE, width = 9, height = 9, main = paste0( "CSI regulon modules (k=", chosen_k, ")" ) )
      figure_audit$created[ figure_audit$figure == "CSI_heatmap" ] <- TRUE
    },
    silent = TRUE
  )
  try( {
      module_state_matrix <- matrix( NA_real_, nrow = length(module_ids), ncol = length(STATE_ORDER), dimnames = list( module_ids, STATE_ORDER ) )
      for (module_id in module_ids) {
        for (state in STATE_ORDER) {
          module_state_matrix[ module_id, state ] <- module_state_summary$mean_module_z_activity[ module_state_summary$CSI_module == module_id & module_state_summary$state == state ]
        }
      }
      pheatmap::pheatmap( module_state_matrix, cluster_rows = TRUE, cluster_cols = FALSE, scale = "none", filename = MODULE_HEATMAP_FILE, width = 6, height = max( 4, nrow(module_state_matrix) * 0.45 ), main = "CSI module mean standardized activity" )
      figure_audit$created[ figure_audit$figure == "module_state_heatmap" ] <- TRUE
    },
    silent = TRUE
  )
}
write_csv( figure_audit, FIGURE_AUDIT_FILE )
module_assignment_pass <- nrow(module_assignment) == EXPECTED_PRIMARY_REGULONS && all(!is.na(module_assignment$CSI_module)) && length( unique(module_assignment$CSI_module) ) == chosen_k
module_activity_pass <- nrow(module_raw_activity) == chosen_k && ncol(module_raw_activity) == EXPECTED_CELLS && all(is.finite(module_raw_activity)) && all(is.finite(module_z_activity))
final_decision <- data.frame( domain = c( "Required input files", "Step 7D Overall PASS", "AUC matrix dimensions and values", "Frozen metadata cell order", "Frozen branch-state counts", "Primary regulon set equals 215", "Informative sensitivity set equals 284", "Synthetic RSS control", "Primary RSS finite and within range", "Sensitivity RSS finite and within range", "At least two evaluable 7 dpi samples", "CSI finite symmetric and within range", "CSI diagonal equals one", "CSI module count selected", "All 215 regulons assigned to a CSI module", "Module activity matrices complete", "Key TF audit produced", "Complete Step 7E RDS produced", "Overall Step 7E" ), pass = c( all(input_audit$exists), step7d_pass, nrow(auc_all) == EXPECTED_REGULONS && ncol(auc_all) == EXPECTED_CELLS && all(is.finite(auc_all)) && all(auc_all >= 0), identical( metadata$CellID, colnames(auc_all) ), state_count_pass, length(primary_regulons) == EXPECTED_PRIMARY_REGULONS, length(informative_regulons) == EXPECTED_INFORMATIVE_REGULONS, synthetic_rss_pass, rss_primary_pass, rss_sensitivity_pass, sample_audit_pass, csi_pass, max( abs( diag(csi_matrix) - 1 ) ) < 1e-12, is.finite(chosen_k) && chosen_k >= 2L, module_assignment_pass, module_activity_pass, file.exists(KEY_TF_FILE) && nrow(key_tf_audit) == length(KEY_TFS), FALSE, FALSE ), stringsAsFactors = FALSE )
final_decision$status <- ifelse( final_decision$pass, "PASS", "FAIL" )
step7e_results <- list( step = "Step7E", generated_at = timestamp(), parameters = list( state_order = STATE_ORDER, primary_min_targets = PRIMARY_MIN_TARGETS, primary_regulons = EXPECTED_PRIMARY_REGULONS, informative_regulons = EXPECTED_INFORMATIVE_REGULONS, min_sample_branch_cells = MIN_SAMPLE_BRANCH_CELLS, csi_k_candidates = CSI_K_CANDIDATES, chosen_k = chosen_k, chosen_k_reason = chosen_k_reason, fixed_k_sensitivity = fixed_k ), metadata = metadata, primary_regulons = primary_regulons, informative_regulons = informative_regulons, RSS_primary = rss_primary, RSS_informative_sensitivity = rss_informative, RSS_primary_long = rss_primary_long, RSS_sensitivity_audit = rss_sensitivity_audit, sample_delta = sample_delta, sample_robustness = sample_robustness, Pearson_correlation = regulon_correlation, CSI = csi_matrix, CSI_hclust = csi_hclust, CSI_k_selection = k_selection, CSI_module_assignment = module_assignment, CSI_module_summary = module_summary, CSI_module_state_summary = module_state_summary, CSI_module_raw_activity = module_raw_activity, CSI_module_z_activity = module_z_activity, CSI_module_RSS = module_rss, key_TF_audit = key_tf_audit )
saveRDS( step7e_results, STEP7E_RESULTS_RDS, compress = TRUE )
rds_pass <- file.exists( STEP7E_RESULTS_RDS )
final_decision$pass[ final_decision$domain == "Complete Step 7E RDS produced" ] <- rds_pass
non_overall_pass <- final_decision$pass[ final_decision$domain != "Overall Step 7E" ]
final_decision$pass[ final_decision$domain == "Overall Step 7E" ] <- all(non_overall_pass)
final_decision$status <- ifelse( final_decision$pass, "PASS", "FAIL" )
write_csv( final_decision, FINAL_DECISION_FILE )
FINAL_STATUS <- if ( isTRUE( final_decision$pass[ final_decision$domain == "Overall Step 7E" ][[1]] ) ) {
  "PASS"
} else {
  "FAIL"
}
step7e_results$status <- FINAL_STATUS
step7e_results$final_decision <- final_decision
saveRDS( step7e_results, STEP7E_RESULTS_RDS, compress = TRUE )
overall_elapsed <- as.numeric( difftime( Sys.time(), overall_start, units = "secs" ) )
top_text <- unlist( lapply( STATE_ORDER, function(state) {
      current <- rss_top[ rss_top$state == state, c( "RSS_rank_in_state", "regulon", "TF", "RSS", "preferred_state", "RSS_margin", "Arm1_minus_Arm2_mean_AUC", "seven_dpi_direction_consistent" ), drop = FALSE ]
      c( paste0( "[", state, "]" ), paste( capture.output( print( current, row.names = FALSE ) ), collapse = "\n" ), "" )
    }
  )
)
report_lines <- c( "======================================================================", "GSE162610 | STEP 7E REPORT", "RSS state specificity and CSI regulon modules", "======================================================================", "", paste0( "Generated: ", timestamp() ), paste0( "Project root: ", PROJECT_ROOT ), paste0( "Overall Step 7E status: ", FINAL_STATUS ), "", "FORMAL SCOPE:", paste0( "Step 7E calculated RSS for Shared_prebranch, Arm_1 and Arm_2, ", "then calculated regulon-regulon CSI and CSI-based modules." ), paste0( "Arm_1 and Arm_2 remain neutral trajectory labels. No fate-driver ", "claim was made." ), "", "FROZEN INPUT:", paste0( "Cells: ", EXPECTED_CELLS, " (Shared_prebranch=116; Arm_1=2,289; Arm_2=1,798)" ), paste0( "All Step 7D regulons: ", EXPECTED_REGULONS ), paste0( "Primary regulons: ", length(primary_regulons), " (targets >= ", PRIMARY_MIN_TARGETS, "; informative)" ), paste0( "Sensitivity regulons: ", length(informative_regulons), " (all informative)" ), "", "RSS METHOD:", paste0( "RSS = 1 - Jensen-Shannon distance between each continuous AUCell ", "profile and each binary state-membership profile." ), paste0( "Synthetic perfect-specificity control RSS: ", synthetic_rss[ "synthetic", "A" ], " (PASS=", synthetic_rss_pass, ")" ), paste0( "Primary RSS range: ", min(rss_primary), " to ", max(rss_primary) ), "", "SAMPLE-LEVEL DIRECTION AUDIT:", paste0( "Minimum cells per arm per sample: ", MIN_SAMPLE_BRANCH_CELLS ), paste0( "Evaluable sample-timepoint units: ", nrow(evaluable_samples) ), paste0( "Evaluable 7 dpi biological samples: ", observed_7dpi_samples ), paste0( "This is a direction-consistency audit; cells were not treated as ", "independent biological replicates." ), "", "CSI METHOD:", paste0( "Pearson correlations were calculated across all 4,203 cells for ", "the 215 primary regulons." ), paste0( "CSI(A,B) is the fraction of other regulons k for which neither ", "cor(A,k) nor cor(B,k) exceeds cor(A,B)." ), paste0( "Clustering: Euclidean distance between CSI profiles; average linkage." ), paste0( "Chosen module count: ", chosen_k, ". ", chosen_k_reason ), paste0( "Fixed k=10 labels were exported as a sensitivity column." ), "", "TIMING:", paste0( "RSS: ", format_duration(rss_elapsed) ), paste0( "Sample audit: ", format_duration(sample_elapsed) ), paste0( "CSI and clustering: ", format_duration(csi_elapsed) ), paste0( "Module activity: ", format_duration(module_elapsed) ), paste0( "Total Step 7E: ", format_duration(overall_elapsed) ), "", "TOP RSS REGULONS BY STATE:", top_text, "RSS PRIMARY-VERSUS-SENSITIVITY AUDIT:", paste( capture.output( print( rss_sensitivity_audit, row.names = FALSE ) ), collapse = "\n" ), "", "CSI MODULE SUMMARY:", paste( capture.output( print( module_summary, row.names = FALSE ) ), collapse = "\n" ), "", "KEY TF RSS/CSI AUDIT:", paste( capture.output( print( key_tf_audit, row.names = FALSE ) ), collapse = "\n" ), "", "FINAL DECISION:", paste( capture.output( print( final_decision, row.names = FALSE ) ), collapse = "\n" ), "", "INTERPRETATION BOUNDARY:", paste0( "RSS-high regulons are state-associated, and CSI modules are ", "co-varying regulatory programs. Neither result proves pre-branch ", "activation or causal fate control." ), "", "NEXT ACTION:", if (FINAL_STATUS == "PASS") {
    paste0( "Proceed to Step 7F: smooth regulon AUC along retained pseudotime ", "and Y_29 signed progress. A fate-driving candidate must separate ", "before or early after the branch point and before terminal programs." )
  } else {
    paste0( "Do not proceed to fate-timing analysis until the failed Step 7E ", "audit domain is resolved." )
  },
  "", "PRIMARY OUTPUTS:", paste0( "RSS primary table: ", safe_path(RSS_PRIMARY_LONG_FILE) ), paste0( "RSS top-state table: ", safe_path(RSS_TOP_FILE) ), paste0( "CSI module assignments: ", safe_path(MODULE_ASSIGNMENT_FILE) ), paste0( "CSI module-state summary: ", safe_path(MODULE_STATE_FILE) ), paste0( "Complete Step 7E RDS: ", safe_path(STEP7E_RESULTS_RDS) ), "", "KEY OUTPUTS:", "1) reports/Step7E_report.txt", "2) tables/04_RSS_top_regulons_by_state.tsv", "3) tables/07_CSI_module_assignments.tsv", "4) tables/08_CSI_module_state_summary.tsv", "5) tables/12_final_decision.csv"
)
writeLines( report_lines, con = REPORT_FILE, useBytes = TRUE )
writeLines( capture.output( sessionInfo() ), con = SESSION_FILE, useBytes = TRUE )
gc()
msg("==============================================================")
msg("STEP 7E COMPLETE")
msg("==============================================================")
msg("Step 7E status: ", FINAL_STATUS)
msg( "RSS: ", nrow(rss_primary), " primary regulons x ", ncol(rss_primary), " states." )
msg( "CSI: ", nrow(csi_matrix), " regulons; ", chosen_k, " selected modules." )
msg( "Evaluable 7 dpi samples: ", observed_7dpi_samples, "." )
msg( "Total time: ", format_duration(overall_elapsed) )
msg( "Report: ", safe_path(REPORT_FILE) )
msg("")
msg("Key outputs:")
msg("1) reports/Step7E_report.txt")
msg("2) tables/04_RSS_top_regulons_by_state.tsv")
msg("3) tables/07_CSI_module_assignments.tsv")
msg("4) tables/08_CSI_module_state_summary.tsv")
msg("5) tables/12_final_decision.csv")
