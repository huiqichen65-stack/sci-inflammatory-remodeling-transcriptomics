rm(list = ls())
gc()
options( stringsAsFactors = FALSE, width = 240, warn = 1, timeout = 86400 )
cat("R parser check PASS: the complete Step 7F file was parsed by source().\n")
set.seed(20260802)
PROJECT_ROOT <- Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
EXPECTED_CELLS <- 4203L
EXPECTED_REGULONS <- 292L
EXPECTED_PRIMARY_REGULONS <- 215L
EXPECTED_STATE_COUNTS <- c( Shared_prebranch = 116L, Arm_1 = 2289L, Arm_2 = 1798L )
STATE_ORDER <- c( "Shared_prebranch", "Arm_1", "Arm_2" )
STATE_COLORS <- c( Shared_prebranch = "#6A3D9A", Arm_1 = "#2F6B9A", Arm_2 = "#B64B45" )
MODEL_DF_SIGNED <- 3L
MODEL_DF_POST <- 2L
POST_GRID_N <- 101L
PRE_GRID_N <- 41L
MIN_CONSECUTIVE_GRID_POINTS <- 5L
PROGRESSIVE_Q_MAX <- 0.05
BRANCH_Q_MAX <- 0.05
MIN_TERMINAL_EFFECT_Z <- 0.25
STRONG_TERMINAL_EFFECT_Z <- 0.50
ONSET_RELATIVE_TERMINAL_FRACTION <- 0.25
ONSET_RESIDUAL_SD_FRACTION <- 0.15
EARLY_PROGRESS_MAX <- 0.25
MIDDLE_PROGRESS_MAX <- 0.50
REFERENCE_ONSET_TOLERANCE <- 0.02
FORMAL_TERMINAL_MODULES <- c( "M1_Arm2_terminal", "M2_Arm2_terminal" )
ARM2_PRIORITY_TFS <- c( "Nfkb1", "Rel", "Nfe2l2", "Rela", "Relb", "Nfkb2", "Irf1", "Egr1", "Atf4", "Mafk", "Maff", "Fosb" )
TRANSITION_STRESS_TFS <- c( "Atf3", "Fos", "Junb" )
ARM1_IDENTITY_TFS <- c( "Mafb", "Spi1", "Irf8", "Irf7", "Stat1", "Cebpb", "Tfeb", "Srebf1", "Bhlhe40", "Borcs8", "Cux1", "Gabpa" )
KEY_TFS <- unique( c( ARM2_PRIORITY_TFS, TRANSITION_STRESS_TFS, ARM1_IDENTITY_TFS ) )
STEP7A_INPUT_DIR <- file.path( PROJECT_ROOT, "04_results", "step7A_GSE162610_pySCENIC_input_environment_audit", "input" )
METADATA_FILE <- file.path( STEP7A_INPUT_DIR, "GSE162610_step7A_frozen_4203_cell_metadata.csv" )
STEP7D_DIR <- file.path( PROJECT_ROOT, "04_results", "step7D_GSE162610_AUCell_regulon_activity_scoring" )
AUC_RDS_FILE <- file.path( STEP7D_DIR, "RDS", "GSE162610_step7D_AUC_regulons_by_cells.rds" )
STEP7D_DECISION_FILE <- file.path( STEP7D_DIR, "tables", "08_final_decision.csv" )
STEP7E_DIR <- file.path( PROJECT_ROOT, "04_results", "step7E_GSE162610_RSS_CSI_regulon_specificity_modules" )
STEP7E_RDS_FILE <- file.path( STEP7E_DIR, "RDS", "GSE162610_step7E_RSS_CSI_results.rds" )
STEP7E_DECISION_FILE <- file.path( STEP7E_DIR, "tables", "12_final_decision.csv" )
OUTDIR <- file.path( PROJECT_ROOT, "04_results", "step7F_GSE162610_regulon_fate_timing" )
TABDIR <- file.path(OUTDIR, "tables")
FIGDIR <- file.path(OUTDIR, "figures")
RDSDIR <- file.path(OUTDIR, "RDS")
REPORTDIR <- file.path(OUTDIR, "reports")
for (directory in c( OUTDIR, TABDIR, FIGDIR, RDSDIR, REPORTDIR )) {
  dir.create( directory, recursive = TRUE, showWarnings = FALSE )
}
INPUT_AUDIT_FILE <- file.path( TABDIR, "00_input_file_audit.csv" )
METADATA_AUDIT_FILE <- file.path( TABDIR, "01_trajectory_metadata_audit.csv" )
MODEL_AUDIT_FILE <- file.path( TABDIR, "02_model_fit_audit.csv" )
REGULON_MASTER_FILE <- file.path( TABDIR, "03_regulon_fate_timing_master.tsv" )
KEY_CANDIDATE_FILE <- file.path( TABDIR, "04_key_candidate_fate_timing.tsv" )
TERMINAL_MODULE_FILE <- file.path( TABDIR, "05_Step6B_terminal_module_timing.tsv" )
CSI_MODULE_FILE <- file.path( TABDIR, "05b_CSI_regulatory_module_timing.tsv" )
PSEUDOTIME_COR_FILE <- file.path( TABDIR, "06_retained_pseudotime_correlations.tsv" )
BINNED_TRAJECTORY_FILE <- file.path( TABDIR, "07_key_analyte_sample_balanced_bins.tsv.gz" )
FINAL_DECISION_FILE <- file.path( TABDIR, "08_final_decision.csv" )
RESULT_RDS_FILE <- file.path( RDSDIR, "GSE162610_step7F_regulon_fate_timing_results.rds" )
REPORT_FILE <- file.path( REPORTDIR, "Step7F_report.txt" )
SESSION_FILE <- file.path( REPORTDIR, "Step7F_sessionInfo.txt" )
ERROR_FILE <- file.path( REPORTDIR, "Step7F_ERROR.txt" )
FIG_KEY_CURVES <- file.path( FIGDIR, "Fig7F_01_key_regulon_Y29_timing_curves.pdf" )
FIG_DIFFERENCE_HEATMAP <- file.path( FIGDIR, "Fig7F_02_key_regulon_Arm1_minus_Arm2_heatmap.pdf" )
FIG_MODULE_CURVES <- file.path( FIGDIR, "Fig7F_03_terminal_and_CSI_module_timing_curves.pdf" )
FIG_ONSET_SCATTER <- file.path( FIGDIR, "Fig7F_04_regulon_onset_vs_terminal_reference.pdf" )
if (!requireNamespace("data.table", quietly = TRUE)) {
  stop( "Required R package 'data.table' is not installed.", call. = FALSE )
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
write_csv <- function(data, path) {
  utils::write.csv( data, file = path, row.names = FALSE, fileEncoding = "UTF-8", na = "" )
}
write_tsv <- function(data, path) {
  data.table::fwrite( data, file = path, sep = "\t", quote = FALSE, na = "" )
}
abort_step <- function(text) {
  writeLines( c( "======================================================================", "GSE162610 | STEP 7F ERROR", "======================================================================", paste0("Time: ", timestamp()), "", text ), con = ERROR_FILE, useBytes = TRUE )
  stop( text, call. = FALSE )
}
decision_pass <- function( decision_table, domain_name ) {
  if ( is.null(decision_table) || !"domain" %in% colnames(decision_table) ) {
    return(FALSE)
  }
  index <- which( as.character(decision_table$domain) == domain_name )
  if (length(index) == 0L) {
    return(FALSE)
  }
  if ("pass" %in% colnames(decision_table)) {
    value <- decision_table$pass[ index[[1]] ]
    return( isTRUE(value) || identical( toupper(as.character(value)), "TRUE" ) || identical( as.character(value), "1" ) )
  }
  if ("status" %in% colnames(decision_table)) {
    return( identical( toupper( as.character( decision_table$status[ index[[1]] ] ) ), "PASS" ) )
  }
  FALSE
}
standardize_state <- function(x) {
  x <- as.character(x)
  x[ x %in% c( "Shared_prebranch", "shared_prebranch", "Pre", "prebranch" ) ] <- "Shared_prebranch"
  x[ x %in% c( "Arm_1", "Arm1", "arm1" ) ] <- "Arm_1"
  x[ x %in% c( "Arm_2", "Arm2", "arm2" ) ] <- "Arm_2"
  x
}
extract_nested_p <- function( smaller_model, larger_model ) {
  comparison <- tryCatch( stats::anova( smaller_model, larger_model ), error = function(e) {
      NULL
    }
  )
  if ( is.null(comparison) || nrow(comparison) < 2L || !"Pr(>F)" %in% colnames(comparison) ) {
    return(NA_real_)
  }
  as.numeric( comparison[ 2L, "Pr(>F)" ] )
}
first_consecutive_true <- function( logical_vector, run_length ) {
  if ( length(logical_vector) < run_length || run_length < 1L ) {
    return(NA_integer_)
  }
  run_object <- rle( logical_vector )
  ends <- cumsum( run_object$lengths )
  starts <- ends - run_object$lengths + 1L
  valid_runs <- which( run_object$values & run_object$lengths >= run_length )
  if (length(valid_runs) == 0L) {
    return(NA_integer_)
  }
  starts[ valid_runs[[1]] ]
}
onset_category <- function(onset_progress) {
  if (!is.finite(onset_progress)) {
    return("not_detected")
  }
  if (onset_progress <= 0.05) {
    return("Y29_immediate")
  }
  if (onset_progress <= EARLY_PROGRESS_MAX) {
    return("early_downstream")
  }
  if (onset_progress <= MIDDLE_PROGRESS_MAX) {
    return("middle_downstream")
  }
  "late_downstream"
}
mean_predictions_over_samples <- function( model, newdata_template, sample_levels ) {
  prediction_matrix <- vapply( sample_levels, function(sample_level) {
      current_newdata <- newdata_template
      current_newdata$sample_factor <- factor( rep( sample_level, nrow(current_newdata) ), levels = sample_levels )
      as.numeric( stats::predict( model, newdata = current_newdata ) )
    },
    FUN.VALUE = numeric( nrow(newdata_template) )
  )
  if (is.null(dim(prediction_matrix))) {
    return( as.numeric(prediction_matrix) )
  }
  rowMeans( prediction_matrix )
}
fit_y29_timing_model <- function( y, metadata, pre_grid, post_grid ) {
  model_data <- data.frame( y = as.numeric(y), sample_factor = factor( metadata$sample_id ), signed_progress = as.numeric( metadata$Y29_signed_progress ), post_progress = as.numeric( metadata$Y29_postbranch_progress ), branch_contrast = as.numeric( metadata$branch_contrast ), stringsAsFactors = FALSE )
  complete_index <- stats::complete.cases( model_data )
  model_data <- model_data[ complete_index, , drop = FALSE ]
  if ( nrow(model_data) < 100L || stats::sd(model_data$y) <= 0 ) {
    return( list( success = FALSE, error = "insufficient_complete_rows_or_zero_variance" ) )
  }
  full_formula <- stats::as.formula( paste0( "y ~ sample_factor + ", "splines::ns(signed_progress, df = ", MODEL_DF_SIGNED, ") + branch_contrast + ", "branch_contrast:splines::ns(post_progress, df = ", MODEL_DF_POST, ")" ) )
  additive_formula <- stats::as.formula( paste0( "y ~ sample_factor + ", "splines::ns(signed_progress, df = ", MODEL_DF_SIGNED, ") + branch_contrast" ) )
  null_formula <- stats::as.formula( paste0( "y ~ sample_factor + ", "splines::ns(signed_progress, df = ", MODEL_DF_SIGNED, ")" ) )
  full_model <- tryCatch( stats::lm( full_formula, data = model_data ), error = function(e) {
      e
    }
  )
  additive_model <- tryCatch( stats::lm( additive_formula, data = model_data ), error = function(e) {
      e
    }
  )
  null_model <- tryCatch( stats::lm( null_formula, data = model_data ), error = function(e) {
      e
    }
  )
  if ( inherits(full_model, "error") || inherits(additive_model, "error") || inherits(null_model, "error") ) {
    error_text <- paste( c( if (inherits(full_model, "error")) {
          conditionMessage(full_model)
        } else {
          NULL
        },
        if (inherits(additive_model, "error")) {
          conditionMessage(additive_model)
        } else {
          NULL
        },
        if (inherits(null_model, "error")) {
          conditionMessage(null_model)
        } else {
          NULL
        }
      ), collapse = " | "
    )
    return( list( success = FALSE, error = error_text ) )
  }
  sample_levels <- levels( model_data$sample_factor )
  pre_template <- data.frame( sample_factor = factor( rep( sample_levels[[1]], length(pre_grid) ), levels = sample_levels ), signed_progress = pre_grid, post_progress = rep( 0, length(pre_grid) ), branch_contrast = rep( 0, length(pre_grid) ) )
  arm1_template <- data.frame( sample_factor = factor( rep( sample_levels[[1]], length(post_grid) ), levels = sample_levels ), signed_progress = post_grid, post_progress = post_grid, branch_contrast = rep( 0.5, length(post_grid) ) )
  arm2_template <- data.frame( sample_factor = factor( rep( sample_levels[[1]], length(post_grid) ), levels = sample_levels ), signed_progress = post_grid, post_progress = post_grid, branch_contrast = rep( -0.5, length(post_grid) ) )
  pre_prediction <- tryCatch( mean_predictions_over_samples( full_model, pre_template, sample_levels ), error = function(e) {
      rep( NA_real_, length(pre_grid) )
    }
  )
  arm1_prediction <- tryCatch( mean_predictions_over_samples( full_model, arm1_template, sample_levels ), error = function(e) {
      rep( NA_real_, length(post_grid) )
    }
  )
  arm2_prediction <- tryCatch( mean_predictions_over_samples( full_model, arm2_template, sample_levels ), error = function(e) {
      rep( NA_real_, length(post_grid) )
    }
  )
  if ( any(!is.finite(pre_prediction)) || any(!is.finite(arm1_prediction)) || any(!is.finite(arm2_prediction)) ) {
    return( list( success = FALSE, error = "non_finite_model_predictions" ) )
  }
  difference_curve <- arm1_prediction - arm2_prediction
  late_index <- post_grid >= stats::quantile( post_grid, probs = 0.80, names = FALSE )
  late_difference <- mean( difference_curve[ late_index ] )
  residual_sd <- stats::sigma( full_model )
  if ( !is.finite(residual_sd) || residual_sd <= 0 ) {
    residual_sd <- stats::sd( model_data$y )
  }
  terminal_effect_z <- late_difference / residual_sd
  onset_threshold <- max( ONSET_RELATIVE_TERMINAL_FRACTION * abs(late_difference), ONSET_RESIDUAL_SD_FRACTION * residual_sd, .Machine$double.eps )
  direction_sign <- sign( late_difference )
  onset_condition <- if ( direction_sign == 0 ) {
    rep( FALSE, length(post_grid) )
  } else {
    sign(difference_curve) == direction_sign & abs(difference_curve) >= onset_threshold
  }
  onset_index <- first_consecutive_true( onset_condition, MIN_CONSECUTIVE_GRID_POINTS )
  onset_progress_raw <- if ( is.na(onset_index) ) {
    NA_real_
  } else {
    post_grid[ onset_index ]
  }
  prebranch_delta <- tail( pre_prediction, 1L ) - pre_prediction[[1]]
  list( success = TRUE, error = "", progressive_shape_p = extract_nested_p( additive_model, full_model ), branch_deviation_p = extract_nested_p( null_model, full_model ), adjusted_r_squared = summary( full_model )$adj.r.squared, residual_sd = residual_sd, prebranch_common_delta = prebranch_delta, prebranch_common_delta_z = prebranch_delta / residual_sd, Y29_difference = difference_curve[[1]], Y29_difference_z = difference_curve[[1]] / residual_sd, late_difference = late_difference, terminal_effect_z = terminal_effect_z, late_direction = if ( late_difference > 0 ) {
      "Arm_1"
    } else if (
      late_difference < 0
    ) {
      "Arm_2"
    } else {
      "Tie"
    },
    onset_threshold = onset_threshold, onset_progress_raw = onset_progress_raw, pre_grid = pre_grid, post_grid = post_grid, pre_prediction = pre_prediction, arm1_prediction = arm1_prediction, arm2_prediction = arm2_prediction, difference_curve = difference_curve
  )
}
rank_bin <- function( values, number_bins ) {
  if (length(values) == 0L) {
    return(integer(0))
  }
  ranked <- rank( values, ties.method = "average", na.last = "keep" )
  bins <- ceiling( ranked / sum(is.finite(ranked)) * number_bins )
  bins[ bins < 1L ] <- 1L
  bins[ bins > number_bins ] <- number_bins
  as.integer(bins)
}
build_sample_balanced_bins <- function( analyte_matrix, metadata, analyte_type ) {
  if ( ncol(analyte_matrix) != nrow(metadata) ) {
    stop( "Analyte matrix and metadata have different cell counts.", call. = FALSE )
  }
  metadata_work <- metadata
  metadata_work$segment <- as.character( metadata_work$step5CR_state )
  metadata_work$Y29_local_axis <- ifelse( metadata_work$segment == "Shared_prebranch", metadata_work$Y29_signed_progress, metadata_work$Y29_postbranch_progress )
  metadata_work$Y29_bin <- NA_integer_
  metadata_work$PT_bin <- NA_integer_
  group_keys <- unique( paste( metadata_work$sample_id, metadata_work$segment, sep = "||" ) )
  for (group_key in group_keys) {
    split_key <- strsplit( group_key, "\\|\\|" )[[1]]
    sample_value <- split_key[[1]]
    state_value <- split_key[[2]]
    index <- which( metadata_work$sample_id == sample_value & metadata_work$segment == state_value )
    y29_bins <- if ( state_value == "Shared_prebranch" ) {
      8L
    } else {
      20L
    }
    metadata_work$Y29_bin[ index ] <- rank_bin( metadata_work$Y29_local_axis[ index ], y29_bins )
    metadata_work$PT_bin[ index ] <- rank_bin( metadata_work$retained_step5BR_pseudotime[ index ], 20L )
  }
  output_parts <- list()
  for ( analyte_index in seq_len( nrow(analyte_matrix) ) ) {
    analyte_name <- rownames( analyte_matrix )[ analyte_index ]
    analyte_values <- as.numeric( analyte_matrix[ analyte_index, , drop = TRUE ] )
    working <- data.table::data.table( sample_id = metadata_work$sample_id, state = metadata_work$segment, Y29_local_axis = metadata_work$Y29_local_axis, Y29_bin = metadata_work$Y29_bin, retained_pseudotime = metadata_work$retained_step5BR_pseudotime, PT_bin = metadata_work$PT_bin, value = analyte_values )
    y29_sample <- working[ is.finite(value) & is.finite(Y29_local_axis) & !is.na(Y29_bin), .( axis_value = mean( Y29_local_axis ), sample_mean = mean( value ), n_cells = .N ), by = .( sample_id, state, bin = Y29_bin ) ]
    y29_balanced <- y29_sample[ , .( axis_value = mean( axis_value ), balanced_mean = mean( sample_mean ), balanced_sd = stats::sd( sample_mean ), n_samples = .N, total_cells = sum( n_cells ) ), by = .( state, bin ) ]
    y29_balanced[ , `:=`( analyte = analyte_name, analyte_type = analyte_type, axis = "Y29_progress" ) ]
    pt_sample <- working[ is.finite(value) & is.finite(retained_pseudotime) & !is.na(PT_bin), .( axis_value = mean( retained_pseudotime ), sample_mean = mean( value ), n_cells = .N ), by = .( sample_id, state, bin = PT_bin ) ]
    pt_balanced <- pt_sample[ , .( axis_value = mean( axis_value ), balanced_mean = mean( sample_mean ), balanced_sd = stats::sd( sample_mean ), n_samples = .N, total_cells = sum( n_cells ) ), by = .( state, bin ) ]
    pt_balanced[ , `:=`( analyte = analyte_name, analyte_type = analyte_type, axis = "retained_pseudotime" ) ]
    output_parts[[ length(output_parts) + 1L ]] <- data.table::rbindlist( list( y29_balanced, pt_balanced ), use.names = TRUE, fill = TRUE )
  }
  data.table::rbindlist( output_parts, use.names = TRUE, fill = TRUE )
}
classify_regulon_timing <- function( row, formal_terminal_reference ) {
  if (!isTRUE(row$model_success)) {
    return("model_failed")
  }
  if (!isTRUE(row$seven_dpi_direction_consistent)) {
    return("sample_sensitive_regulon")
  }
  if ( is.na(row$progressive_shape_q) || row$progressive_shape_q > PROGRESSIVE_Q_MAX ) {
    return("state_associated_without_progressive_timing")
  }
  if ( is.na(row$branch_deviation_q) || row$branch_deviation_q > BRANCH_Q_MAX ) {
    return("progressive_shape_without_stable_branch_deviation")
  }
  if ( is.na(row$onset_progress_qualified) ) {
    return("significant_but_onset_not_detected")
  }
  if ( !is.na(row$RSS_preferred_state) && row$RSS_preferred_state %in% c( "Arm_1", "Arm_2" ) && row$RSS_preferred_state != row$late_direction ) {
    return("RSS_dynamic_direction_discordant")
  }
  if ( row$late_direction == "Arm_2" ) {
    if ( is.finite( formal_terminal_reference ) && row$onset_progress_qualified <= formal_terminal_reference + REFERENCE_ONSET_TOLERANCE ) {
      if ( abs( row$terminal_effect_z ) >= STRONG_TERMINAL_EFFECT_Z ) {
        return( "early_Arm2_fate_timing_candidate_requires_Step7G" )
      }
      return( "early_weak_Arm2_timing_candidate_requires_Step7G" )
    }
    return( "late_Arm2_branch_associated_regulon" )
  }
  if ( row$late_direction == "Arm_1" ) {
    if ( row$onset_progress_qualified <= EARLY_PROGRESS_MAX ) {
      return( "early_Arm1_identity_retention_candidate_not_fate_driver" )
    }
    return( "late_Arm1_branch_associated_regulon" )
  }
  "ambiguous_timing"
}
overall_start <- Sys.time()
msg("==============================================================")
msg("GSE162610 Step 7F | regulon fate-timing analysis")
msg("==============================================================")
msg("Formal timing model is centered on Y_29.")
msg("Retained Step5B-R pseudotime is used as an orthogonal descriptive audit.")
required_files <- c( METADATA_FILE, AUC_RDS_FILE, STEP7D_DECISION_FILE, STEP7E_RDS_FILE, STEP7E_DECISION_FILE )
missing_files <- required_files[ !file.exists( required_files ) ]
if (length(missing_files) > 0L) {
  abort_step( paste0( "Required files are missing:\n", paste( safe_path( missing_files ), collapse = "\n" ) ) )
}
input_audit <- data.frame( role = c( "Step7A_frozen_metadata", "Step7D_AUC_matrix", "Step7D_final_decision", "Step7E_complete_results", "Step7E_final_decision" ), path = safe_path( required_files ), exists = file.exists( required_files ), size_bytes = file.info( required_files )$size, stringsAsFactors = FALSE )
write_csv( input_audit, INPUT_AUDIT_FILE )
step7d_decision <- data.table::fread( STEP7D_DECISION_FILE, data.table = FALSE, showProgress = FALSE )
step7e_decision <- data.table::fread( STEP7E_DECISION_FILE, data.table = FALSE, showProgress = FALSE )
step7d_pass <- decision_pass( step7d_decision, "Overall Step 7D" )
step7e_pass <- decision_pass( step7e_decision, "Overall Step 7E" )
if (!step7d_pass) {
  abort_step( "Step 7D does not have an Overall PASS decision." )
}
if (!step7e_pass) {
  abort_step( "Step 7E does not have an Overall PASS decision." )
}
msg("Loading the continuous regulon AUC matrix...")
auc_all <- readRDS( AUC_RDS_FILE )
if ( !is.matrix(auc_all) || nrow(auc_all) != EXPECTED_REGULONS || ncol(auc_all) != EXPECTED_CELLS ) {
  abort_step( paste0( "Unexpected AUCell dimensions: ", paste( dim(auc_all), collapse = " x " ) ) )
}
storage.mode(auc_all) <- "double"
if ( any(!is.finite(auc_all)) || any(auc_all < 0) || anyDuplicated(rownames(auc_all)) > 0L || anyDuplicated(colnames(auc_all)) > 0L ) {
  abort_step( "AUCell matrix failed finite, non-negative or duplicate-name checks." )
}
step7e <- readRDS( STEP7E_RDS_FILE )
required_step7e_fields <- c( "primary_regulons", "CSI_module_assignment", "sample_robustness", "CSI_module_raw_activity", "final_decision" )
missing_step7e_fields <- setdiff( required_step7e_fields, names(step7e) )
if (length(missing_step7e_fields) > 0L) {
  abort_step( paste0( "Step 7E RDS is missing: ", paste( missing_step7e_fields, collapse = ", " ) ) )
}
if ( !decision_pass( step7e$final_decision, "Overall Step 7E" ) ) {
  abort_step( "The Step 7E RDS itself does not contain an Overall PASS decision." )
}
primary_regulons <- as.character( step7e$primary_regulons )
if ( length(primary_regulons) != EXPECTED_PRIMARY_REGULONS || !all(primary_regulons %in% rownames(auc_all)) ) {
  abort_step( "Step 7E primary-regulon set is incomplete or does not match Step 7D." )
}
auc_primary <- auc_all[ primary_regulons, , drop = FALSE ]
metadata <- data.table::fread( METADATA_FILE, data.table = FALSE, check.names = FALSE, showProgress = FALSE )
required_metadata_columns <- c( "cell", "sample_id", "timepoint", "step5CR_state", "graph_distance_from_Y29", "Y29_signed_progress", "Y29_postbranch_progress", "retained_step5BR_pseudotime" )
missing_metadata_columns <- setdiff( required_metadata_columns, colnames(metadata) )
if (length(missing_metadata_columns) > 0L) {
  abort_step( paste0( "Step 7A metadata is missing: ", paste( missing_metadata_columns, collapse = ", " ) ) )
}
module_score_columns <- grep( "^step6B_mean_norm_score__", colnames(metadata), value = TRUE )
if (length(module_score_columns) < 5L) {
  abort_step( paste0( "Fewer than five Step6B module-score columns were found. Found: ", paste( module_score_columns, collapse = ", " ) ) )
}
metadata$cell <- as.character( metadata$cell )
metadata$step5CR_state <- standardize_state( metadata$step5CR_state )
metadata_index <- match( colnames(auc_primary), metadata$cell )
if (anyNA(metadata_index)) {
  abort_step( "Some AUCell cells are absent from Step 7A metadata." )
}
metadata <- metadata[ metadata_index, , drop = FALSE ]
if ( !identical( metadata$cell, colnames(auc_primary) ) ) {
  abort_step( "Metadata order does not exactly match AUCell columns." )
}
metadata$sample_id <- as.character( metadata$sample_id )
metadata$timepoint <- as.character( metadata$timepoint )
numeric_metadata_columns <- c( "graph_distance_from_Y29", "Y29_signed_progress", "Y29_postbranch_progress", "retained_step5BR_pseudotime", module_score_columns )
for (column_name in numeric_metadata_columns) {
  metadata[[column_name]] <- suppressWarnings( as.numeric( metadata[[column_name]] ) )
}
metadata$branch_contrast <- 0
metadata$branch_contrast[ metadata$step5CR_state == "Arm_1" ] <- 0.5
metadata$branch_contrast[ metadata$step5CR_state == "Arm_2" ] <- -0.5
observed_state_counts <- table( factor( metadata$step5CR_state, levels = STATE_ORDER ) )
state_count_pass <- all( as.integer( observed_state_counts[ names( EXPECTED_STATE_COUNTS ) ] ) == as.integer( EXPECTED_STATE_COUNTS ) )
if (!state_count_pass) {
  abort_step( paste0( "Frozen state counts differ: ", paste( names(observed_state_counts), observed_state_counts, sep = "=", collapse = "; " ) ) )
}
critical_numeric_pass <- all( vapply( c( "Y29_signed_progress", "Y29_postbranch_progress", "retained_step5BR_pseudotime" ), function(column_name) {
      all( is.finite( metadata[[column_name]] ) )
    },
    FUN.VALUE = logical(1)
  )
)
module_numeric_pass <- all( vapply( module_score_columns, function(column_name) {
      all( is.finite( metadata[[column_name]] ) )
    },
    FUN.VALUE = logical(1)
  )
)
if ( !critical_numeric_pass || !module_numeric_pass ) {
  abort_step( "Trajectory or Step6B module-score metadata contains non-finite values." )
}
if ( any( metadata$Y29_signed_progress[ metadata$step5CR_state == "Shared_prebranch" ] > 1e-10 ) || any( metadata$Y29_postbranch_progress[ metadata$step5CR_state == "Shared_prebranch" ] > 1e-10 ) || any( metadata$Y29_signed_progress[ metadata$step5CR_state %in% c( "Arm_1", "Arm_2" ) ] < -1e-10 ) ) {
  abort_step( "Y_29 progress signs are inconsistent with the frozen Step5C-R definition." )
}
metadata_audit <- data.frame( metric = c( "metadata_rows", "exact_AUC_cell_order", "Shared_prebranch_cells", "Arm_1_cells", "Arm_2_cells", "unique_samples", "unique_timepoints", "minimum_Y29_signed_progress", "maximum_Y29_signed_progress", "maximum_Y29_postbranch_progress", "minimum_retained_pseudotime", "maximum_retained_pseudotime", "Step6B_module_score_columns", "all_critical_trajectory_values_finite", "all_module_scores_finite" ), value = c( nrow(metadata), identical( metadata$cell, colnames(auc_primary) ), observed_state_counts[ "Shared_prebranch" ], observed_state_counts[ "Arm_1" ], observed_state_counts[ "Arm_2" ], length( unique( metadata$sample_id ) ), length( unique( metadata$timepoint ) ), min( metadata$Y29_signed_progress ), max( metadata$Y29_signed_progress ), max( metadata$Y29_postbranch_progress ), min( metadata$retained_step5BR_pseudotime ), max( metadata$retained_step5BR_pseudotime ), length( module_score_columns ), critical_numeric_pass, module_numeric_pass ), stringsAsFactors = FALSE )
write_csv( metadata_audit, METADATA_AUDIT_FILE )
pre_grid <- seq( min( metadata$Y29_signed_progress[ metadata$step5CR_state == "Shared_prebranch" ] ), 0, length.out = PRE_GRID_N )
post_grid <- seq( 0, max( metadata$Y29_postbranch_progress ), length.out = POST_GRID_N )
model_start <- Sys.time()
msg( "Fitting Y_29 timing models for ", length(primary_regulons), " primary regulons..." )
regulon_fit_list <- vector( "list", length(primary_regulons) )
names(regulon_fit_list) <- primary_regulons
regulon_rows <- vector( "list", length(primary_regulons) )
for ( regulon_index in seq_along(primary_regulons) ) {
  regulon_name <- primary_regulons[ regulon_index ]
  fit_result <- fit_y29_timing_model( y = auc_primary[ regulon_name, , drop = TRUE ], metadata = metadata, pre_grid = pre_grid, post_grid = post_grid )
  regulon_fit_list[[regulon_name]] <- fit_result
  if (isTRUE(fit_result$success)) {
    regulon_rows[[regulon_index]] <- data.frame( regulon = regulon_name, model_success = TRUE, model_error = "", progressive_shape_p = fit_result$progressive_shape_p, branch_deviation_p = fit_result$branch_deviation_p, adjusted_r_squared = fit_result$adjusted_r_squared, residual_sd = fit_result$residual_sd, prebranch_common_delta = fit_result$prebranch_common_delta, prebranch_common_delta_z = fit_result$prebranch_common_delta_z, Y29_difference = fit_result$Y29_difference, Y29_difference_z = fit_result$Y29_difference_z, late_difference = fit_result$late_difference, terminal_effect_z = fit_result$terminal_effect_z, late_direction = fit_result$late_direction, onset_threshold = fit_result$onset_threshold, onset_progress_raw = fit_result$onset_progress_raw, stringsAsFactors = FALSE )
  } else {
    regulon_rows[[regulon_index]] <- data.frame( regulon = regulon_name, model_success = FALSE, model_error = fit_result$error, progressive_shape_p = NA_real_, branch_deviation_p = NA_real_, adjusted_r_squared = NA_real_, residual_sd = NA_real_, prebranch_common_delta = NA_real_, prebranch_common_delta_z = NA_real_, Y29_difference = NA_real_, Y29_difference_z = NA_real_, late_difference = NA_real_, terminal_effect_z = NA_real_, late_direction = NA_character_, onset_threshold = NA_real_, onset_progress_raw = NA_real_, stringsAsFactors = FALSE )
  }
  if ( regulon_index %% 25L == 0L ) {
    msg( "  completed ", regulon_index, "/", length(primary_regulons) )
  }
}
regulon_master <- do.call( rbind, regulon_rows )
regulon_master$progressive_shape_q <- stats::p.adjust( regulon_master$progressive_shape_p, method = "BH" )
regulon_master$branch_deviation_q <- stats::p.adjust( regulon_master$branch_deviation_p, method = "BH" )
regulon_master$onset_progress_qualified <- ifelse( regulon_master$model_success & is.finite( regulon_master$onset_progress_raw ) & is.finite( regulon_master$progressive_shape_q ) & regulon_master$progressive_shape_q <= PROGRESSIVE_Q_MAX & is.finite( regulon_master$branch_deviation_q ) & regulon_master$branch_deviation_q <= BRANCH_Q_MAX & abs( regulon_master$terminal_effect_z ) >= MIN_TERMINAL_EFFECT_Z, regulon_master$onset_progress_raw, NA_real_ )
regulon_master$onset_category <- vapply( regulon_master$onset_progress_qualified, onset_category, FUN.VALUE = character(1) )
model_elapsed <- as.numeric( difftime( Sys.time(), model_start, units = "secs" ) )
model_success_fraction <- mean( regulon_master$model_success )
model_audit <- data.frame( metric = c( "primary_regulons", "successful_models", "failed_models", "model_success_fraction", "progressive_q_le_0.05", "branch_deviation_q_le_0.05", "qualified_onset_detected", "Arm_1_late_direction", "Arm_2_late_direction", "model_elapsed_seconds" ), value = c( nrow(regulon_master), sum(regulon_master$model_success), sum(!regulon_master$model_success), model_success_fraction, sum( regulon_master$progressive_shape_q <= PROGRESSIVE_Q_MAX, na.rm = TRUE ), sum( regulon_master$branch_deviation_q <= BRANCH_Q_MAX, na.rm = TRUE ), sum( is.finite( regulon_master$onset_progress_qualified ) ), sum( regulon_master$late_direction == "Arm_1", na.rm = TRUE ), sum( regulon_master$late_direction == "Arm_2", na.rm = TRUE ), model_elapsed ), stringsAsFactors = FALSE )
write_csv( model_audit, MODEL_AUDIT_FILE )
msg( "Regulon models complete in ", format_duration(model_elapsed), "; success fraction=", round( model_success_fraction, 4 ), "." )
msg("Fitting Step6B gene-program timing models...")
step6b_module_names <- sub( "^step6B_mean_norm_score__", "", module_score_columns )
step6b_score_matrix <- t( as.matrix( metadata[ , module_score_columns, drop = FALSE ] ) )
rownames(step6b_score_matrix) <- step6b_module_names
colnames(step6b_score_matrix) <- metadata$cell
step6b_fit_list <- vector( "list", length(step6b_module_names) )
names(step6b_fit_list) <- step6b_module_names
step6b_rows <- vector( "list", length(step6b_module_names) )
for ( module_index in seq_along(step6b_module_names) ) {
  module_name <- step6b_module_names[ module_index ]
  fit_result <- fit_y29_timing_model( y = step6b_score_matrix[ module_name, , drop = TRUE ], metadata = metadata, pre_grid = pre_grid, post_grid = post_grid )
  step6b_fit_list[[module_name]] <- fit_result
  step6b_rows[[module_index]] <- data.frame( module_name = module_name, model_success = isTRUE( fit_result$success ), model_error = if ( isTRUE(fit_result$success) ) {
      ""
    } else {
      fit_result$error
    },
    progressive_shape_p = if ( isTRUE(fit_result$success) ) {
      fit_result$progressive_shape_p
    } else {
      NA_real_
    },
    branch_deviation_p = if ( isTRUE(fit_result$success) ) {
      fit_result$branch_deviation_p
    } else {
      NA_real_
    },
    residual_sd = if ( isTRUE(fit_result$success) ) {
      fit_result$residual_sd
    } else {
      NA_real_
    },
    prebranch_common_delta_z = if ( isTRUE(fit_result$success) ) {
      fit_result$prebranch_common_delta_z
    } else {
      NA_real_
    },
    Y29_difference_z = if ( isTRUE(fit_result$success) ) {
      fit_result$Y29_difference_z
    } else {
      NA_real_
    },
    late_difference = if ( isTRUE(fit_result$success) ) {
      fit_result$late_difference
    } else {
      NA_real_
    },
    terminal_effect_z = if ( isTRUE(fit_result$success) ) {
      fit_result$terminal_effect_z
    } else {
      NA_real_
    },
    late_direction = if ( isTRUE(fit_result$success) ) {
      fit_result$late_direction
    } else {
      NA_character_
    },
    onset_progress_raw = if ( isTRUE(fit_result$success) ) {
      fit_result$onset_progress_raw
    } else {
      NA_real_
    },
    formal_terminal_reference = module_name %in%
      FORMAL_TERMINAL_MODULES, stringsAsFactors = FALSE
  )
}
step6b_timing <- do.call( rbind, step6b_rows )
step6b_timing$progressive_shape_q <- stats::p.adjust( step6b_timing$progressive_shape_p, method = "BH" )
step6b_timing$branch_deviation_q <- stats::p.adjust( step6b_timing$branch_deviation_p, method = "BH" )
step6b_timing$onset_progress_qualified <- ifelse( step6b_timing$model_success & is.finite( step6b_timing$onset_progress_raw ) & step6b_timing$progressive_shape_q <= PROGRESSIVE_Q_MAX & step6b_timing$branch_deviation_q <= BRANCH_Q_MAX & abs( step6b_timing$terminal_effect_z ) >= MIN_TERMINAL_EFFECT_Z, step6b_timing$onset_progress_raw, NA_real_ )
step6b_timing$onset_category <- vapply( step6b_timing$onset_progress_qualified, onset_category, FUN.VALUE = character(1) )
formal_reference_values <- step6b_timing$onset_progress_qualified[ step6b_timing$formal_terminal_reference & step6b_timing$late_direction == "Arm_2" ]
formal_reference_values <- formal_reference_values[ is.finite( formal_reference_values ) ]
terminal_reference_source <- "M1_Arm2_terminal_and_M2_Arm2_terminal"
if (length(formal_reference_values) > 0L) {
  formal_terminal_reference <- min( formal_reference_values )
} else {
  fallback_index <- grepl( "Arm2_terminal", step6b_timing$module_name, fixed = TRUE ) & step6b_timing$late_direction == "Arm_2" & is.finite( step6b_timing$onset_progress_qualified )
  fallback_values <- step6b_timing$onset_progress_qualified[ fallback_index ]
  if (length(fallback_values) > 0L) {
    formal_terminal_reference <- min( fallback_values )
    terminal_reference_source <- "fallback_all_Arm2_terminal_modules"
  } else {
    formal_terminal_reference <- NA_real_
    terminal_reference_source <- "no_valid_terminal_reference"
  }
}
step6b_timing$formal_terminal_reference_progress <- formal_terminal_reference
step6b_timing$terminal_reference_source <- terminal_reference_source
write_tsv( step6b_timing, TERMINAL_MODULE_FILE )
msg("Fitting Step7E CSI regulatory-module timing models...")
csi_module_matrix <- as.matrix( step7e$CSI_module_raw_activity )
if ( ncol(csi_module_matrix) != EXPECTED_CELLS ) {
  abort_step( "Step 7E CSI module activity has an unexpected cell count." )
}
csi_cell_index <- match( metadata$cell, colnames(csi_module_matrix) )
if (anyNA(csi_cell_index)) {
  abort_step( "CSI module activity is missing one or more frozen cells." )
}
csi_module_matrix <- csi_module_matrix[ , csi_cell_index, drop = FALSE ]
if ( !identical( colnames(csi_module_matrix), metadata$cell ) ) {
  abort_step( "CSI module activity could not be placed in exact metadata order." )
}
csi_module_names <- rownames( csi_module_matrix )
csi_fit_list <- vector( "list", length(csi_module_names) )
names(csi_fit_list) <- csi_module_names
csi_rows <- vector( "list", length(csi_module_names) )
for ( module_index in seq_along(csi_module_names) ) {
  module_name <- csi_module_names[ module_index ]
  fit_result <- fit_y29_timing_model( y = csi_module_matrix[ module_name, , drop = TRUE ], metadata = metadata, pre_grid = pre_grid, post_grid = post_grid )
  csi_fit_list[[module_name]] <- fit_result
  csi_rows[[module_index]] <- data.frame( CSI_module = module_name, model_success = isTRUE( fit_result$success ), progressive_shape_p = if ( isTRUE(fit_result$success) ) {
      fit_result$progressive_shape_p
    } else {
      NA_real_
    },
    branch_deviation_p = if ( isTRUE(fit_result$success) ) {
      fit_result$branch_deviation_p
    } else {
      NA_real_
    },
    prebranch_common_delta_z = if ( isTRUE(fit_result$success) ) {
      fit_result$prebranch_common_delta_z
    } else {
      NA_real_
    },
    Y29_difference_z = if ( isTRUE(fit_result$success) ) {
      fit_result$Y29_difference_z
    } else {
      NA_real_
    },
    terminal_effect_z = if ( isTRUE(fit_result$success) ) {
      fit_result$terminal_effect_z
    } else {
      NA_real_
    },
    late_direction = if ( isTRUE(fit_result$success) ) {
      fit_result$late_direction
    } else {
      NA_character_
    },
    onset_progress_raw = if ( isTRUE(fit_result$success) ) {
      fit_result$onset_progress_raw
    } else {
      NA_real_
    },
    stringsAsFactors = FALSE
  )
}
csi_timing <- do.call( rbind, csi_rows )
csi_timing$progressive_shape_q <- stats::p.adjust( csi_timing$progressive_shape_p, method = "BH" )
csi_timing$branch_deviation_q <- stats::p.adjust( csi_timing$branch_deviation_p, method = "BH" )
csi_timing$onset_progress_qualified <- ifelse( csi_timing$model_success & is.finite( csi_timing$onset_progress_raw ) & csi_timing$progressive_shape_q <= PROGRESSIVE_Q_MAX & csi_timing$branch_deviation_q <= BRANCH_Q_MAX & abs( csi_timing$terminal_effect_z ) >= MIN_TERMINAL_EFFECT_Z, csi_timing$onset_progress_raw, NA_real_ )
csi_timing$onset_category <- vapply( csi_timing$onset_progress_qualified, onset_category, FUN.VALUE = character(1) )
write_tsv( csi_timing, CSI_MODULE_FILE )
module_assignment <- as.data.frame( step7e$CSI_module_assignment )
sample_robustness <- as.data.frame( step7e$sample_robustness )
required_assignment_columns <- c( "regulon", "TF", "n_targets", "preferred_state", "preferred_RSS", "RSS_margin", "CSI_module" )
missing_assignment_columns <- setdiff( required_assignment_columns, colnames(module_assignment) )
if (length(missing_assignment_columns) > 0L) {
  abort_step( paste0( "Step 7E CSI_module_assignment is missing: ", paste( missing_assignment_columns, collapse = ", " ) ) )
}
assignment_index <- match( regulon_master$regulon, module_assignment$regulon )
if (anyNA(assignment_index)) {
  abort_step( "Some primary regulons are missing from the Step 7E assignment table." )
}
regulon_master$TF <- module_assignment$TF[ assignment_index ]
regulon_master$n_targets <- module_assignment$n_targets[ assignment_index ]
regulon_master$RSS_preferred_state <- module_assignment$preferred_state[ assignment_index ]
regulon_master$preferred_RSS <- module_assignment$preferred_RSS[ assignment_index ]
regulon_master$RSS_margin <- module_assignment$RSS_margin[ assignment_index ]
regulon_master$CSI_module <- module_assignment$CSI_module[ assignment_index ]
robustness_index <- match( regulon_master$regulon, sample_robustness$regulon )
regulon_master$seven_dpi_direction_consistent <- sample_robustness$seven_dpi_direction_consistent[ robustness_index ]
regulon_master$sample_direction_consistency_fraction <- sample_robustness$sample_direction_consistency_fraction[ robustness_index ]
regulon_master$direction_matches_RSS <- regulon_master$late_direction == regulon_master$RSS_preferred_state
regulon_master$lead_before_formal_terminal_program <- if ( is.finite( formal_terminal_reference ) ) {
    formal_terminal_reference - regulon_master$onset_progress_qualified
  } else {
    NA_real_
  }
regulon_master$prebranch_common_trend <- ifelse( regulon_master$prebranch_common_delta_z >= MIN_TERMINAL_EFFECT_Z, "increasing_before_Y29", ifelse( regulon_master$prebranch_common_delta_z <= -MIN_TERMINAL_EFFECT_Z, "decreasing_before_Y29", "weak_or_flat_before_Y29" ) )
regulon_master$timing_classification <- vapply( seq_len( nrow(regulon_master) ), function(row_index) {
      classify_regulon_timing( regulon_master[ row_index, , drop = FALSE ], formal_terminal_reference )
    },
    FUN.VALUE = character(1)
  )
regulon_master$formal_terminal_reference_progress <- formal_terminal_reference
regulon_master$terminal_reference_source <- terminal_reference_source
regulon_master <- regulon_master[ order( is.na( regulon_master$onset_progress_qualified ), regulon_master$onset_progress_qualified, regulon_master$progressive_shape_q, -abs( regulon_master$terminal_effect_z ) ), , drop = FALSE ]
write_tsv( regulon_master, REGULON_MASTER_FILE )
key_candidate_table <- regulon_master[ regulon_master$TF %in% KEY_TFS, , drop = FALSE ]
key_candidate_table$priority_group <- ifelse( key_candidate_table$TF %in% ARM2_PRIORITY_TFS, "Arm2_priority", ifelse( key_candidate_table$TF %in% TRANSITION_STRESS_TFS, "transition_stress", "Arm1_identity_retention" ) )
key_candidate_table <- key_candidate_table[ order( match( key_candidate_table$priority_group, c( "Arm2_priority", "transition_stress", "Arm1_identity_retention" ) ), is.na( key_candidate_table$onset_progress_qualified ), key_candidate_table$onset_progress_qualified, key_candidate_table$progressive_shape_q ), , drop = FALSE ]
write_tsv( key_candidate_table, KEY_CANDIDATE_FILE )
msg("Calculating retained-pseudotime correlations...")
pseudotime_rows <- list()
for ( regulon_index in seq_along(primary_regulons) ) {
  regulon_name <- primary_regulons[ regulon_index ]
  for (state_name in STATE_ORDER) {
    state_index <- metadata$step5CR_state == state_name
    correlation_value <- suppressWarnings( stats::cor( auc_primary[ regulon_name, state_index, drop = TRUE ], metadata$retained_step5BR_pseudotime[ state_index ], method = "spearman", use = "complete.obs" ) )
    pseudotime_rows[[ length(pseudotime_rows) + 1L ]] <- data.frame( regulon = regulon_name, state = state_name, spearman_rho = correlation_value, n_cells = sum( state_index ), stringsAsFactors = FALSE )
  }
}
pseudotime_correlations <- do.call( rbind, pseudotime_rows )
write_tsv( pseudotime_correlations, PSEUDOTIME_COR_FILE )
msg("Building sample-balanced trajectory bins for key analytes...")
key_regulons <- key_candidate_table$regulon
key_regulon_matrix <- auc_primary[ key_regulons, , drop = FALSE ]
key_regulon_bins <- build_sample_balanced_bins( key_regulon_matrix, metadata, "regulon_AUC" )
step6b_bins <- build_sample_balanced_bins( step6b_score_matrix, metadata, "Step6B_gene_program" )
csi_bins <- build_sample_balanced_bins( csi_module_matrix, metadata, "CSI_regulatory_module" )
all_bins <- data.table::rbindlist( list( key_regulon_bins, step6b_bins, csi_bins ), use.names = TRUE, fill = TRUE )
data.table::fwrite( all_bins, file = BINNED_TRAJECTORY_FILE, sep = "\t", quote = FALSE, na = "", compress = "gzip" )
msg("Generating Step 7F figures...")
key_plot_table <- key_candidate_table[ key_candidate_table$model_success, , drop = FALSE ]
key_plot_table <- key_plot_table[ order( match( key_plot_table$priority_group, c( "Arm2_priority", "transition_stress", "Arm1_identity_retention" ) ), is.na( key_plot_table$onset_progress_qualified ), key_plot_table$onset_progress_qualified ), , drop = FALSE ]
grDevices::pdf( FIG_KEY_CURVES, width = 12, height = 9, onefile = TRUE )
old_par <- graphics::par( mfrow = c( 3, 3 ), mar = c( 3.4, 3.6, 3.2, 1.0 ), oma = c( 0, 0, 1.0, 0 ) )
for ( row_index in seq_len( nrow(key_plot_table) ) ) {
  regulon_name <- key_plot_table$regulon[ row_index ]
  fit_result <- regulon_fit_list[[regulon_name]]
  x_all <- c( fit_result$pre_grid, fit_result$post_grid, fit_result$post_grid )
  y_all <- c( fit_result$pre_prediction, fit_result$arm1_prediction, fit_result$arm2_prediction )
  y_range <- range( y_all, finite = TRUE )
  graphics::plot( fit_result$pre_grid, fit_result$pre_prediction, type = "l", lwd = 2, col = STATE_COLORS[ "Shared_prebranch" ], xlim = range( c( fit_result$pre_grid, fit_result$post_grid ) ), ylim = y_range, xlab = "Y_29 signed/postbranch progress", ylab = "Adjusted regulon AUC", main = paste0( key_plot_table$TF[ row_index ], " | ", key_plot_table$late_direction[ row_index ] ) )
  graphics::lines( fit_result$post_grid, fit_result$arm1_prediction, lwd = 2, col = STATE_COLORS[ "Arm_1" ] )
  graphics::lines( fit_result$post_grid, fit_result$arm2_prediction, lwd = 2, col = STATE_COLORS[ "Arm_2" ] )
  graphics::abline( v = 0, lty = 2, col = "grey40" )
  if ( is.finite( key_plot_table$onset_progress_qualified[ row_index ] ) ) {
    graphics::abline( v = key_plot_table$onset_progress_qualified[ row_index ], lty = 3, col = "black" )
  }
  if ( is.finite( formal_terminal_reference ) ) {
    graphics::abline( v = formal_terminal_reference, lty = 4, col = "grey20" )
  }
  graphics::mtext( paste0( "q_prog=", format( key_plot_table$progressive_shape_q[ row_index ], digits = 2, scientific = TRUE ), "; onset=", ifelse( is.finite( key_plot_table$onset_progress_qualified[ row_index ] ), round( key_plot_table$onset_progress_qualified[ row_index ], 3 ), "NA" ), "; z_late=", round( key_plot_table$terminal_effect_z[ row_index ], 2 ) ), side = 3, line = 0.25, cex = 0.62 )
  graphics::legend( "topleft", legend = c( "Shared", "Arm_1", "Arm_2" ), col = STATE_COLORS[ c( "Shared_prebranch", "Arm_1", "Arm_2" ) ], lty = 1, lwd = 2, cex = 0.55, bty = "n" )
}
graphics::par( old_par )
grDevices::dev.off()
heatmap_regulons <- key_plot_table$regulon
difference_matrix <- do.call( rbind, lapply( heatmap_regulons, function(regulon_name) {
      fit_result <- regulon_fit_list[[regulon_name]]
      fit_result$difference_curve / fit_result$residual_sd
    }
  )
)
rownames(difference_matrix) <- key_plot_table$TF
difference_matrix[ !is.finite( difference_matrix ) ] <- 0
heatmap_order <- order( key_plot_table$late_direction, key_plot_table$onset_progress_qualified, na.last = TRUE )
difference_matrix <- difference_matrix[ heatmap_order, , drop = FALSE ]
heatmap_limit <- max( abs( difference_matrix ), na.rm = TRUE )
if ( !is.finite(heatmap_limit) || heatmap_limit <= 0 ) {
  heatmap_limit <- 1
}
heatmap_breaks <- seq( -heatmap_limit, heatmap_limit, length.out = 101L )
heatmap_palette <- grDevices::colorRampPalette( c( STATE_COLORS[ "Arm_2" ], "white", STATE_COLORS[ "Arm_1" ] ) )( 100L )
grDevices::pdf( FIG_DIFFERENCE_HEATMAP, width = 10, height = max( 7, nrow(difference_matrix) * 0.25 ) )
graphics::par( mar = c( 4, 8, 3, 2 ) )
graphics::image( x = post_grid, y = seq_len( nrow(difference_matrix) ), z = t( difference_matrix[ nrow(difference_matrix):1, , drop = FALSE ] ), col = heatmap_palette, breaks = heatmap_breaks, xlab = "Post-Y_29 progress", ylab = "", yaxt = "n", main = "Standardized Arm_1 minus Arm_2 regulon-AUC difference" )
graphics::axis( side = 2, at = seq_len( nrow(difference_matrix) ), labels = rev( rownames(difference_matrix) ), las = 2, cex.axis = 0.62 )
if ( is.finite( formal_terminal_reference ) ) {
  graphics::abline( v = formal_terminal_reference, lty = 3, lwd = 1.2 )
}
grDevices::dev.off()
combined_module_names <- c( paste0( "Step6B::", step6b_module_names ), paste0( "CSI::", csi_module_names ) )
grDevices::pdf( FIG_MODULE_CURVES, width = 12, height = 8, onefile = TRUE )
old_par <- graphics::par( mfrow = c( 2, 3 ), mar = c( 3.5, 3.7, 3.0, 1.0 ) )
for (combined_name in combined_module_names) {
  if ( startsWith( combined_name, "Step6B::" ) ) {
    module_name <- sub( "^Step6B::", "", combined_name )
    fit_result <- step6b_fit_list[[module_name]]
    module_title <- module_name
  } else {
    module_name <- sub( "^CSI::", "", combined_name )
    fit_result <- csi_fit_list[[module_name]]
    module_title <- paste0( "CSI ", module_name )
  }
  if (!isTRUE(fit_result$success)) {
    graphics::plot.new()
    graphics::title( main = paste0( module_title, " | model failed" ) )
    next
  }
  y_range <- range( c( fit_result$pre_prediction, fit_result$arm1_prediction, fit_result$arm2_prediction ), finite = TRUE )
  graphics::plot( fit_result$pre_grid, fit_result$pre_prediction, type = "l", lwd = 2, col = STATE_COLORS[ "Shared_prebranch" ], xlim = range( c( fit_result$pre_grid, fit_result$post_grid ) ), ylim = y_range, xlab = "Y_29 progress", ylab = "Adjusted module score", main = module_title )
  graphics::lines( fit_result$post_grid, fit_result$arm1_prediction, lwd = 2, col = STATE_COLORS[ "Arm_1" ] )
  graphics::lines( fit_result$post_grid, fit_result$arm2_prediction, lwd = 2, col = STATE_COLORS[ "Arm_2" ] )
  graphics::abline( v = 0, lty = 2 )
  if ( is.finite( formal_terminal_reference ) ) {
    graphics::abline( v = formal_terminal_reference, lty = 4 )
  }
}
graphics::par( old_par )
grDevices::dev.off()
plot_table <- regulon_master[ is.finite( regulon_master$onset_progress_qualified ), , drop = FALSE ]
grDevices::pdf( FIG_ONSET_SCATTER, width = 9, height = 7 )
point_colors <- ifelse( plot_table$late_direction == "Arm_1", STATE_COLORS[ "Arm_1" ], STATE_COLORS[ "Arm_2" ] )
graphics::plot( plot_table$onset_progress_qualified, plot_table$terminal_effect_z, pch = 16, col = point_colors, xlab = "Qualified regulon onset after Y_29", ylab = "Late Arm_1 minus Arm_2 effect / residual SD", main = "Regulon timing relative to Step6B terminal-program onset" )
graphics::abline( h = 0, lty = 2 )
if ( is.finite( formal_terminal_reference ) ) {
  graphics::abline( v = formal_terminal_reference, lty = 3, lwd = 2 )
}
label_index <- which( plot_table$TF %in% KEY_TFS )
graphics::text( plot_table$onset_progress_qualified[ label_index ], plot_table$terminal_effect_z[ label_index ], labels = plot_table$TF[ label_index ], pos = 4, cex = 0.65 )
graphics::legend( "topright", legend = c( "Arm_1 late direction", "Arm_2 late direction", "Step6B terminal reference" ), col = c( STATE_COLORS[ "Arm_1" ], STATE_COLORS[ "Arm_2" ], "black" ), pch = c( 16, 16, NA ), lty = c( NA, NA, 3 ), bty = "n" )
grDevices::dev.off()
formal_reference_pass <- is.finite( formal_terminal_reference )
key_output_pass <- nrow(key_candidate_table) >= 20L && file.exists( KEY_CANDIDATE_FILE )
figure_pass <- all( file.exists( c( FIG_KEY_CURVES, FIG_DIFFERENCE_HEATMAP, FIG_MODULE_CURVES, FIG_ONSET_SCATTER ) ) )
final_decision <- data.frame( domain = c( "Required input files", "Step 7D Overall PASS", "Step 7E Overall PASS", "AUC dimensions and finite values", "Frozen metadata exact cell order", "Frozen branch-state counts", "Y29 trajectory metadata finite and signed correctly", "Step6B module-score metadata", "At least 95 percent regulon models successful", "Regulon q values and timing outputs produced", "Formal Step6B terminal-program reference available", "CSI regulatory-module timing complete", "Key candidate timing table produced", "Step7F figures produced", "Complete Step7F RDS produced", "Overall Step 7F" ), pass = c( all( input_audit$exists ), step7d_pass, step7e_pass, nrow(auc_all) == EXPECTED_REGULONS && ncol(auc_all) == EXPECTED_CELLS && all( is.finite( auc_all ) ), identical( metadata$cell, colnames(auc_primary) ), state_count_pass, critical_numeric_pass, module_numeric_pass && length(module_score_columns) >= 5L, model_success_fraction >= 0.95, file.exists( REGULON_MASTER_FILE ) && nrow(regulon_master) == EXPECTED_PRIMARY_REGULONS, formal_reference_pass, all( csi_timing$model_success ), key_output_pass, figure_pass, FALSE, FALSE ), stringsAsFactors = FALSE )
final_decision$status <- ifelse( final_decision$pass, "PASS", "FAIL" )
step7f_results <- list( step = "Step7F", generated_at = timestamp(), parameters = list( model_df_signed = MODEL_DF_SIGNED, model_df_post = MODEL_DF_POST, progressive_q_max = PROGRESSIVE_Q_MAX, branch_q_max = BRANCH_Q_MAX, minimum_terminal_effect_z = MIN_TERMINAL_EFFECT_Z, strong_terminal_effect_z = STRONG_TERMINAL_EFFECT_Z, minimum_consecutive_grid_points = MIN_CONSECUTIVE_GRID_POINTS, formal_terminal_modules = FORMAL_TERMINAL_MODULES, formal_terminal_reference = formal_terminal_reference, terminal_reference_source = terminal_reference_source ), metadata = metadata, primary_regulons = primary_regulons, regulon_timing_master = regulon_master, key_candidate_timing = key_candidate_table, Step6B_terminal_module_timing = step6b_timing, CSI_module_timing = csi_timing, retained_pseudotime_correlations = pseudotime_correlations, regulon_fit_curves = regulon_fit_list, Step6B_fit_curves = step6b_fit_list, CSI_fit_curves = csi_fit_list, final_decision = final_decision )
saveRDS( step7f_results, RESULT_RDS_FILE, compress = TRUE )
rds_pass <- file.exists( RESULT_RDS_FILE )
final_decision$pass[ final_decision$domain == "Complete Step7F RDS produced" ] <- rds_pass
non_overall_pass <- final_decision$pass[ final_decision$domain != "Overall Step 7F" ]
final_decision$pass[ final_decision$domain == "Overall Step 7F" ] <- all( non_overall_pass )
final_decision$status <- ifelse( final_decision$pass, "PASS", "FAIL" )
write_csv( final_decision, FINAL_DECISION_FILE )
FINAL_STATUS <- if ( isTRUE( final_decision$pass[ final_decision$domain == "Overall Step 7F" ][[1]] ) ) {
  "PASS"
} else {
  "FAIL"
}
step7f_results$status <- FINAL_STATUS
step7f_results$final_decision <- final_decision
saveRDS( step7f_results, RESULT_RDS_FILE, compress = TRUE )
overall_elapsed <- as.numeric( difftime( Sys.time(), overall_start, units = "secs" ) )
classification_counts <- as.data.frame( table( regulon_master$timing_classification, useNA = "ifany" ), stringsAsFactors = FALSE )
colnames(classification_counts) <- c( "timing_classification", "n_regulons" )
top_early_arm2 <- regulon_master[ grepl( "^early_.*Arm2", regulon_master$timing_classification ), , drop = FALSE ]
top_early_arm2 <- top_early_arm2[ order( top_early_arm2$onset_progress_qualified, -abs( top_early_arm2$terminal_effect_z ) ), , drop = FALSE ]
top_early_arm1 <- regulon_master[ regulon_master$timing_classification == "early_Arm1_identity_retention_candidate_not_fate_driver", , drop = FALSE ]
top_early_arm1 <- top_early_arm1[ order( top_early_arm1$onset_progress_qualified, -abs( top_early_arm1$terminal_effect_z ) ), , drop = FALSE ]
report_lines <- c( "======================================================================", "GSE162610 | STEP 7F REPORT", "Regulon AUC fate-timing around Y_29", "======================================================================", "", paste0( "Generated: ", timestamp() ), paste0( "Project root: ", PROJECT_ROOT ), paste0( "Overall Step 7F status: ", FINAL_STATUS ), "", "FORMAL SCOPE:", paste0( "Step 7F fitted branchpoint-centered continuous-AUC models for all ", EXPECTED_PRIMARY_REGULONS, " primary regulons." ), paste0( "The formal model matched the Step5C-R structure: sample + common ", "signed-progress spline + anchored branch contrast + branch-by-postprogress ", "spline interaction." ), paste0( "Retained Step5B-R pseudotime was used as an orthogonal descriptive ", "correlation and sample-balanced trajectory audit." ), "", "FROZEN INPUT:", paste0( "Cells: ", EXPECTED_CELLS, " (Shared_prebranch=116; Arm_1=2,289; Arm_2=1,798)" ), paste0( "Primary regulons: ", EXPECTED_PRIMARY_REGULONS ), paste0( "Step6B module-score columns: ", length( module_score_columns ) ), "", "MODEL QUALITY:", paste0( "Successful regulon models: ", sum( regulon_master$model_success ), "/", EXPECTED_PRIMARY_REGULONS, " (", round( 100 * model_success_fraction, 1 ), "%)" ), paste0( "Progressive-shape q<=0.05: ", sum( regulon_master$progressive_shape_q <= PROGRESSIVE_Q_MAX, na.rm = TRUE ) ), paste0( "Stable branch-deviation q<=0.05: ", sum( regulon_master$branch_deviation_q <= BRANCH_Q_MAX, na.rm = TRUE ) ), paste0( "Qualified onset detected: ", sum( is.finite( regulon_master$onset_progress_qualified ) ) ), "", "ONSET DEFINITION:", paste0( "The onset was the earliest run of at least ", MIN_CONSECUTIVE_GRID_POINTS, " post-Y_29 grid points with the terminal direction and an absolute ", "Arm_1-minus-Arm_2 difference exceeding max(", ONSET_RELATIVE_TERMINAL_FRACTION, " x late difference, ", ONSET_RESIDUAL_SD_FRACTION, " x residual SD)." ), paste0( "A qualified onset additionally required progressive q<=", PROGRESSIVE_Q_MAX, ", branch-deviation q<=", BRANCH_Q_MAX, " and absolute late effect >=", MIN_TERMINAL_EFFECT_Z, " residual SD." ), "", "TERMINAL PROGRAM REFERENCE:", paste0( "Reference source: ", terminal_reference_source ), paste0( "Earliest formal terminal-program onset: ", ifelse( is.finite( formal_terminal_reference ), round( formal_terminal_reference, 4 ), "NA" ) ), paste( capture.output( print( step6b_timing, row.names = FALSE ) ), collapse = "\n" ), "", "TIMING CLASSIFICATION COUNTS:", paste( capture.output( print( classification_counts, row.names = FALSE ) ), collapse = "\n" ), "", "EARLY ARM_2 FATE-TIMING CANDIDATES:", if ( nrow(top_early_arm2) > 0L ) {
    paste( capture.output( print( utils::head( top_early_arm2[ , c( "TF", "regulon", "onset_progress_qualified", "terminal_effect_z", "prebranch_common_trend", "seven_dpi_direction_consistent", "CSI_module", "timing_classification" ), drop = FALSE ], 30L ), row.names = FALSE ) ), collapse = "\n" )
  } else {
    "None met all Step 7F timing criteria."
  },
  "", "EARLY ARM_1 / IDENTITY-RETENTION CANDIDATES:", if ( nrow(top_early_arm1) > 0L ) {
    paste( capture.output( print( utils::head( top_early_arm1[ , c( "TF", "regulon", "onset_progress_qualified", "terminal_effect_z", "prebranch_common_trend", "seven_dpi_direction_consistent", "CSI_module", "timing_classification" ), drop = FALSE ], 30L ), row.names = FALSE ) ), collapse = "\n" )
  } else {
    "None met the early identity-retention timing criteria."
  },
  "", "KEY CANDIDATE AUDIT:", paste( capture.output( print( key_candidate_table[ , c( "priority_group", "TF", "regulon", "RSS_preferred_state", "late_direction", "progressive_shape_q", "branch_deviation_q", "onset_progress_qualified", "onset_category", "terminal_effect_z", "prebranch_common_trend", "lead_before_formal_terminal_program", "seven_dpi_direction_consistent", "sample_direction_consistency_fraction", "CSI_module", "timing_classification" ), drop = FALSE ], row.names = FALSE ) ), collapse = "\n" ), "", "FINAL DECISION:", paste( capture.output( print( final_decision, row.names = FALSE ) ), collapse = "\n" ), "", "INTERPRETATION BOUNDARY:", paste0( "An early Step 7F candidate has timing evidence only. It must not yet be ", "called a fate driver. Step 7G must show direction-consistent overlap of ", "its motif-supported targets with progressive genes, strong progressive ", "genes, dynamic modules and branch markers." ), paste0( "Arm_1 candidates remain identity-retention or metabolic/lipid-handling ", "candidates unless a distinct Arm_1 terminal program is demonstrated." ), "", "NEXT ACTION:", if ( FINAL_STATUS == "PASS" ) {
    paste0( "Proceed to Step 7G: intersect regulon targets with the 212 progressive ", "genes, 172 strong progressive genes, Step6B dynamic modules and Step6A ", "high-confidence branch markers." )
  } else {
    paste0( "Do not proceed to Step 7G until the failed Step 7F audit domain is resolved." )
  },
  "", "PRIMARY OUTPUTS:", paste0( "Regulon timing master: ", safe_path( REGULON_MASTER_FILE ) ), paste0( "Key candidate timing: ", safe_path( KEY_CANDIDATE_FILE ) ), paste0( "Terminal module timing: ", safe_path( TERMINAL_MODULE_FILE ) ), paste0( "Complete Step 7F RDS: ", safe_path( RESULT_RDS_FILE ) ), "", "KEY OUTPUTS:", "1) reports/Step7F_report.txt", "2) tables/04_key_candidate_fate_timing.tsv", "3) tables/05_Step6B_terminal_module_timing.tsv", "4) figures/Fig7F_01_key_regulon_Y29_timing_curves.pdf", "5) tables/08_final_decision.csv"
)
writeLines( report_lines, con = REPORT_FILE, useBytes = TRUE )
writeLines( capture.output( sessionInfo() ), con = SESSION_FILE, useBytes = TRUE )
gc()
msg("==============================================================")
msg("STEP 7F COMPLETE")
msg("==============================================================")
msg("Step 7F status: ", FINAL_STATUS)
msg( "Successful regulon models: ", sum( regulon_master$model_success ), "/", EXPECTED_PRIMARY_REGULONS )
msg( "Qualified regulon onsets: ", sum( is.finite( regulon_master$onset_progress_qualified ) ) )
msg( "Formal terminal-program reference: ", ifelse( is.finite( formal_terminal_reference ), round( formal_terminal_reference, 4 ), "NA" ) )
msg( "Total time: ", format_duration( overall_elapsed ) )
msg( "Report: ", safe_path( REPORT_FILE ) )
msg("")
msg("Key outputs:")
msg("1) reports/Step7F_report.txt")
msg("2) tables/04_key_candidate_fate_timing.tsv")
msg("3) tables/05_Step6B_terminal_module_timing.tsv")
msg("4) figures/Fig7F_01_key_regulon_Y29_timing_curves.pdf")
msg("5) tables/08_final_decision.csv")
