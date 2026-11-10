rm(list = ls())
gc()
options( stringsAsFactors = FALSE, width = 240, warn = 1, timeout = 86400, future.globals.maxSize = 32 * 1024^3 )
set.seed(20260802)
PROJECT_ROOT <- Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
FORMAL_TIMEPOINTS <- c( "1dpi", "3dpi", "7dpi" )
FORMAL_BRANCH_SENDERS <- c( "Microglia_Arm_1", "Microglia_Arm_2" )
EXCLUDED_FORMAL_TARGETS <- c( "Microglia_Arm_1", "Microglia_Arm_2", "Microglia_Shared_prebranch" )
MIN_CELLS_PER_GROUP <- 20L
MIN_TARGET_CELLS_PER_SAMPLE <- 20L
MIN_SENDER_CELLS_PER_SAMPLE <- 20L
MIN_RECEPTOR_PERCENT_EXPRESSED <- 5
CELLCHAT_P_THRESHOLD <- 0.05
MIN_BRANCH_FOLD_CHANGE <- 1.25
MAX_SAMPLE_AUDIT_INTERACTIONS <- 800L
MAX_SHORTLIST_ROWS <- 200L
MAX_FIGURE_ROWS <- 40L
FAST_DISCOVERY_MODE <- TRUE
FAST_MAX_CELLS_PER_SAMPLE_GROUP <- 150L
FAST_KEEP_ALL_GROUPS <- c( "Microglia_Shared_prebranch", "Microglia_Arm_1", "Microglia_Arm_2" )
FAST_SAMPLING_SEED <- 20260802L
FAST_NBOOT <- 20L
REUSE_VALID_CHECKPOINTS <- TRUE
N_WORKERS <- NA_integer_
CELLCHAT_PARAMETERS <- list( database = "CellChatDB.mouse", database_subset = "all_categories", probability_type = "triMean", trim = 0.10, raw.use = TRUE, population.size = FALSE, distance.use = FALSE, min.cells = MIN_CELLS_PER_GROUP, nboot = FAST_NBOOT, seed.use = 20260802L )
ATLAS_RDS_FILE <- file.path( PROJECT_ROOT, "04_results", "step1F_GSE162610_broad_annotation", "GSE162610_broad_annotation_audit.rds" )
STEP8A_DIR <- file.path( PROJECT_ROOT, "04_results", "step8A_GSE162610_CellChat_input_preflight" )
STEP8A_MANIFEST_FILE <- file.path( STEP8A_DIR, "RDS", "GSE162610_step8A_CellChat_input_manifest.rds" )
STEP8A_DECISION_FILE <- file.path( STEP8A_DIR, "tables", "08_final_decision.csv" )
STEP8A_CANDIDATE_LIGAND_FILE <- file.path( STEP8A_DIR, "tables", "06_candidate_regulon_ligand_catalog.tsv" )
OUTDIR <- file.path( PROJECT_ROOT, "04_results", "step8B_PREVIEW_GSE162610_CellChat_nboot20" )
TABDIR <- file.path(OUTDIR, "tables")
FIGDIR <- file.path(OUTDIR, "figures")
RDSDIR <- file.path(OUTDIR, "RDS")
CHECKPOINTDIR <- file.path(RDSDIR, "timepoint_checkpoints")
REPORTDIR <- file.path(OUTDIR, "reports")
for (directory in c( OUTDIR, TABDIR, FIGDIR, RDSDIR, CHECKPOINTDIR, REPORTDIR )) {
  dir.create( directory, recursive = TRUE, showWarnings = FALSE )
}
INPUT_AUDIT_FILE <- file.path( TABDIR, "00_input_file_audit.csv" )
TIMEPOINT_RUN_SUMMARY_FILE <- file.path( TABDIR, "01_timepoint_CellChat_run_summary.csv" )
FAST_SAMPLING_AUDIT_FILE <- file.path( TABDIR, "01a_fast_stratified_sampling_audit.csv" )
FAST_SAMPLING_SUMMARY_FILE <- file.path( TABDIR, "01b_fast_sampling_summary_by_timepoint.csv" )
NETWORK_SUMMARY_FILE <- file.path( TABDIR, "02_network_summary_by_timepoint.csv" )
PATHWAY_FLOW_FILE <- file.path( TABDIR, "03_pathway_information_flow_by_timepoint.tsv" )
ALL_INTERACTIONS_FILE <- file.path( TABDIR, "04_all_CellChat_interactions_by_timepoint.tsv.gz" )
PATHWAY_TARGET_DIFF_FILE <- file.path( TABDIR, "05_7dpi_branch_pathway_target_differential.tsv" )
PATHWAY_OUTGOING_DIFF_FILE <- file.path( TABDIR, "06_7dpi_branch_outgoing_pathway_differential.tsv" )
LR_DIFF_FILE <- file.path( TABDIR, "07_7dpi_branch_LR_differential.tsv" )
SAMPLE_AUDIT_LONG_FILE <- file.path( TABDIR, "08_7dpi_LR_sample_expression_audit_long.tsv.gz" )
SAMPLE_SHORTLIST_FILE <- file.path( TABDIR, "09_7dpi_sample_supported_LR_shortlist.tsv" )
CANDIDATE_INTEGRATED_FILE <- file.path( TABDIR, "10_Step7G_candidate_CellChat_integrated_shortlist.tsv" )
RECIPIENT_AUDIT_FILE <- file.path( TABDIR, "11_7dpi_recipient_sample_coverage.csv" )
FINAL_DECISION_FILE <- file.path( TABDIR, "12_final_decision.csv" )
MERGED_CELLCHAT_FILE <- file.path( RDSDIR, "GSE162610_step8B_merged_CellChat_1_3_7dpi.rds" )
RESULT_RDS_FILE <- file.path( RDSDIR, "GSE162610_step8B_PREVIEW_CellChat_results.rds" )
REPORT_FILE <- file.path( REPORTDIR, "Step8B_PREVIEW_report.txt" )
SESSION_FILE <- file.path( REPORTDIR, "Step8B_sessionInfo.txt" )
ERROR_FILE <- file.path( REPORTDIR, "Step8B_ERROR.txt" )
FIG_NETWORK_SUMMARY <- file.path( FIGDIR, "Fig8B_01_network_summary_by_timepoint.pdf" )
FIG_PATHWAY_FLOW <- file.path( FIGDIR, "Fig8B_02_pathway_information_flow_heatmap.pdf" )
FIG_BRANCH_PATHWAY <- file.path( FIGDIR, "Fig8B_03_7dpi_Arm1_vs_Arm2_outgoing_pathways.pdf" )
FIG_BRANCH_LR <- file.path( FIGDIR, "Fig8B_04_7dpi_top_branch_differential_LR.pdf" )
FIG_CANDIDATE_AXES <- file.path( FIGDIR, "Fig8B_05_Step7G_candidate_supported_axes.pdf" )
REQUIRED_PACKAGES <- c( "Seurat", "SeuratObject", "Matrix", "data.table", "ggplot2", "CellChat", "future" )
MISSING_PACKAGES <- REQUIRED_PACKAGES[ !vapply( REQUIRED_PACKAGES, requireNamespace, logical(1), quietly = TRUE ) ]
if (length(MISSING_PACKAGES) > 0L) {
  stop( paste0( "Missing required R package(s): ", paste( MISSING_PACKAGES, collapse = ", " ), "." ), call. = FALSE )
}
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(data.table)
  library(ggplot2)
  library(CellChat)
  library(future)
})
N_WORKERS <- max( 1L, min( 4L, as.integer( future::availableCores() ) - 1L ) )
if ( !is.finite(N_WORKERS) || N_WORKERS < 1L ) {
  N_WORKERS <- 1L
}
timestamp <- function() {
  format( Sys.time(), "%Y-%m-%d %H:%M:%S" )
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
    return( paste0( round(seconds, 1), " sec" ) )
  }
  if (seconds < 3600) {
    return( paste0( round( seconds / 60, 1 ), " min" ) )
  }
  paste0( round( seconds / 3600, 2 ), " h" )
}
write_csv_utf8 <- function( data, path ) {
  data.table::fwrite( data, file = path, sep = ",", quote = TRUE, na = "NA" )
}
write_tsv <- function( data, path ) {
  data.table::fwrite( data, file = path, sep = "\t", quote = FALSE, na = "NA" )
}
write_tsv_gz <- function( data, path ) {
  data.table::fwrite( data, file = path, sep = "\t", quote = FALSE, na = "NA", compress = "gzip" )
}
abort_step <- function(text) {
  writeLines( c( "======================================================================", "GSE162610 | STEP 8B PREVIEW ERROR", "======================================================================", paste0( "Time: ", timestamp() ), "", text ), con = ERROR_FILE, useBytes = TRUE )
  stop( text, call. = FALSE )
}
as_logical_safe <- function(x) {
  if (is.logical(x)) {
    return(x)
  }
  normalized <- toupper( trimws( as.character(x) ) )
  output <- rep( NA, length(normalized) )
  output[ normalized %in% c( "TRUE", "T", "1", "PASS", "YES" ) ] <- TRUE
  output[ normalized %in% c( "FALSE", "F", "0", "FAIL", "NO" ) ] <- FALSE
  output
}
decision_pass <- function( decision_table, domain_name ) {
  if ( !"domain" %in% colnames(decision_table) ) {
    return(FALSE)
  }
  index <- which( as.character( decision_table$domain ) == domain_name )
  if (length(index) == 0L) {
    return(FALSE)
  }
  if ( "pass" %in% colnames(decision_table) ) {
    value <- as_logical_safe( decision_table$pass[ index[[1L]] ] )
    return( isTRUE( value[[1L]] ) )
  }
  if ( "status" %in% colnames(decision_table) ) {
    return( identical( toupper( as.character( decision_table$status[ index[[1L]] ] ) ), "PASS" ) )
  }
  FALSE
}
add_decision_row <- function( decision_list, domain, pass, detail = "", status_if_fail = "FAIL" ) {
  decision_list[[ length(decision_list) + 1L ]] <- data.frame( domain = domain, pass = isTRUE(pass), status = if ( isTRUE(pass) ) {
      "PASS"
    } else {
      status_if_fail
    },
    detail = as.character(detail), stringsAsFactors = FALSE
  )
  decision_list
}
get_rna_data <- function(object) {
  output <- tryCatch( SeuratObject::LayerData( object = object, assay = "RNA", layer = "data" ), error = function(e) {
      NULL
    }
  )
  if (!is.null(output)) {
    return(output)
  }
  output <- tryCatch( Seurat::GetAssayData( object = object, assay = "RNA", slot = "data" ), error = function(e) {
      NULL
    }
  )
  if (is.null(output)) {
    stop( "Unable to retrieve normalized RNA data from the full atlas.", call. = FALSE )
  }
  output
}
load_cellchat_mouse_db <- function() {
  database <- tryCatch( get( "CellChatDB.mouse", envir = asNamespace( "CellChat" ), inherits = FALSE ), error = function(e) {
      NULL
    }
  )
  if (!is.null(database)) {
    return(database)
  }
  temporary_environment <- new.env( parent = emptyenv() )
  try( utils::data( "CellChatDB.mouse", package = "CellChat", envir = temporary_environment ), silent = TRUE )
  if ( exists( "CellChatDB.mouse", envir = temporary_environment, inherits = FALSE ) ) {
    return( get( "CellChatDB.mouse", envir = temporary_environment, inherits = FALSE ) )
  }
  stop( "CellChatDB.mouse could not be loaded.", call. = FALSE )
}
parameters_identical <- function( first, second ) {
  identical( first, second )
}
call_with_supported_arguments <- function( function_object, arguments ) {
  supported <- names( formals( function_object ) )
  if ( "..." %in% supported ) {
    return( do.call( function_object, arguments ) )
  }
  arguments <- arguments[ names(arguments) %in% supported ]
  do.call( function_object, arguments )
}
checkpoint_path <- function(timepoint) {
  file.path( CHECKPOINTDIR, paste0( "GSE162610_step8B_PREVIEW_CellChat_", timepoint, "_checkpoint.rds" ) )
}
checkpoint_is_valid <- function( checkpoint, timepoint, expected_cells, expected_groups ) {
  if ( !is.list(checkpoint) || !all( c( "cellchat", "parameters", "timepoint", "cells", "groups" ) %in% names(checkpoint) ) ) {
    return(FALSE)
  }
  if ( !inherits( checkpoint$cellchat, "CellChat" ) ) {
    return(FALSE)
  }
  if ( !identical( checkpoint$timepoint, timepoint ) || !identical( checkpoint$cells, expected_cells ) || !identical( checkpoint$groups, expected_groups ) || !parameters_identical( checkpoint$parameters, CELLCHAT_PARAMETERS ) ) {
    return(FALSE)
  }
  if ( !"netP" %in% slotNames( checkpoint$cellchat ) ) {
    return(FALSE)
  }
  probability_array <- tryCatch( checkpoint$cellchat@netP$prob, error = function(e) {
      NULL
    }
  )
  !is.null( probability_array ) && length( probability_array ) > 0L
}
safe_subset_communication <- function( cellchat, slot_name = "net", sources_use = NULL, targets_use = NULL, threshold = 1 ) {
  arguments <- list( object = cellchat, slot.name = slot_name, sources.use = sources_use, targets.use = targets_use, thresh = threshold )
  output <- call_with_supported_arguments( CellChat::subsetCommunication, arguments )
  as.data.frame( output, stringsAsFactors = FALSE )
}
matrix_to_long <- function( matrix_object, value_name, timepoint ) {
  if ( is.null(matrix_object) || length(matrix_object) == 0L ) {
    return( data.frame( source = character(0), target = character(0), timepoint = character(0), value = numeric(0), stringsAsFactors = FALSE ) )
  }
  output <- as.data.frame( as.table( matrix_object ), stringsAsFactors = FALSE )
  colnames(output) <- c( "source", "target", value_name )
  output$timepoint <- timepoint
  output[ , c( "timepoint", "source", "target", value_name ), drop = FALSE ]
}
pathway_array_to_long <- function( probability_array, timepoint ) {
  if ( is.null(probability_array) || length( dim( probability_array ) ) != 3L ) {
    return( data.frame( timepoint = character(0), source = character(0), target = character(0), pathway = character(0), probability = numeric(0), stringsAsFactors = FALSE ) )
  }
  output <- as.data.frame( as.table( probability_array ), stringsAsFactors = FALSE )
  colnames(output) <- c( "source", "target", "pathway", "probability" )
  output$timepoint <- timepoint
  output[ , c( "timepoint", "source", "target", "pathway", "probability" ), drop = FALSE ]
}
information_flow_table <- function( probability_array, timepoint ) {
  if ( is.null(probability_array) || length( dim( probability_array ) ) != 3L ) {
    return( data.frame( timepoint = character(0), pathway = character(0), information_flow = numeric(0), stringsAsFactors = FALSE ) )
  }
  pathway_names <- dimnames( probability_array )[[3L]]
  information_flow <- vapply( seq_along( pathway_names ), function(pathway_index) {
      sum( probability_array[ , , pathway_index, drop = TRUE ], na.rm = TRUE )
    },
    FUN.VALUE = numeric(1)
  )
  data.frame( timepoint = timepoint, pathway = pathway_names, information_flow = information_flow, stringsAsFactors = FALSE )
}
standardize_interaction_columns <- function( interaction_table, timepoint ) {
  if (nrow(interaction_table) == 0L) {
    return( data.frame( timepoint = character(0), source = character(0), target = character(0), ligand = character(0), receptor = character(0), interaction_name = character(0), pathway_name = character(0), probability = numeric(0), p_value = numeric(0), stringsAsFactors = FALSE ) )
  }
  column_map <- list( source = c( "source" ), target = c( "target" ), ligand = c( "ligand" ), receptor = c( "receptor" ), interaction_name = c( "interaction_name", "interaction_name_2" ), pathway_name = c( "pathway_name", "pathway" ), probability = c( "prob", "probability" ), p_value = c( "pval", "p_value" ) )
  resolved <- list()
  for (canonical_name in names(column_map)) {
    candidates <- column_map[[canonical_name]]
    hit <- candidates[ candidates %in% colnames(interaction_table) ]
    if (length(hit) > 0L) {
      resolved[[canonical_name]] <- hit[[1L]]
    } else {
      resolved[[canonical_name]] <- NA_character_
    }
  }
  required_resolved <- c( "source", "target", "ligand", "receptor", "probability" )
  missing_required <- required_resolved[ is.na( unlist( resolved[ required_resolved ] ) ) ]
  if (length(missing_required) > 0L) {
    stop( paste0( "subsetCommunication output is missing required fields: ", paste( missing_required, collapse = ", " ), ". Current columns: ", paste( colnames(interaction_table), collapse = ", " ) ), call. = FALSE )
  }
  output <- data.frame( timepoint = timepoint, source = as.character( interaction_table[[resolved$source]] ), target = as.character( interaction_table[[resolved$target]] ), ligand = as.character( interaction_table[[resolved$ligand]] ), receptor = as.character( interaction_table[[resolved$receptor]] ), interaction_name = if ( !is.na( resolved$interaction_name ) ) {
      as.character( interaction_table[[resolved$interaction_name]] )
    } else {
      paste0( interaction_table[[resolved$ligand]], "_", interaction_table[[resolved$receptor]] )
    },
    pathway_name = if ( !is.na( resolved$pathway_name ) ) {
      as.character( interaction_table[[resolved$pathway_name]] )
    } else {
      NA_character_
    },
    probability = as.numeric( interaction_table[[resolved$probability]] ), p_value = if ( !is.na( resolved$p_value ) ) {
      as.numeric( interaction_table[[resolved$p_value]] )
    } else {
      NA_real_
    },
    stringsAsFactors = FALSE
  )
  additional_columns <- setdiff( colnames(interaction_table), unname( unlist( resolved ) ) )
  additional_columns <- additional_columns[ !is.na( additional_columns ) ]
  for (column_name in additional_columns) {
    if ( !column_name %in% colnames(output) ) {
      output[[column_name]] <- interaction_table[[column_name]]
    }
  }
  output
}
aggregate_branch_interactions <- function( interaction_table, source_group ) {
  current <- interaction_table[ interaction_table$source == source_group, , drop = FALSE ]
  if (nrow(current) == 0L) {
    return( data.frame( target = character(0), ligand = character(0), receptor = character(0), interaction_name = character(0), pathway_name = character(0), probability = numeric(0), p_value = numeric(0), stringsAsFactors = FALSE ) )
  }
  current$pathway_name[ is.na( current$pathway_name ) ] <- "Unassigned"
  aggregate_probability <- stats::aggregate( probability ~ target + ligand + receptor + interaction_name + pathway_name, data = current, FUN = sum, na.rm = TRUE )
  if ( all( is.na( current$p_value ) ) ) {
    aggregate_p <- aggregate_probability[ , c( "target", "ligand", "receptor", "interaction_name", "pathway_name" ), drop = FALSE ]
    aggregate_p$p_value <- NA_real_
  } else {
    aggregate_p <- stats::aggregate( p_value ~ target + ligand + receptor + interaction_name + pathway_name, data = current, FUN = function(values) {
        finite_values <- values[ is.finite( values ) ]
        if (length(finite_values) == 0L) {
          return(1)
        }
        min( finite_values )
      }
    )
  }
  merge( aggregate_probability, aggregate_p, by = c( "target", "ligand", "receptor", "interaction_name", "pathway_name" ), all = TRUE, sort = FALSE )
}
combine_branch_interactions <- function( arm1_table, arm2_table ) {
  key_columns <- c( "target", "ligand", "receptor", "interaction_name", "pathway_name" )
  arm1 <- arm1_table
  arm2 <- arm2_table
  colnames(arm1)[ colnames(arm1) == "probability" ] <- "probability_Arm_1"
  colnames(arm1)[ colnames(arm1) == "p_value" ] <- "p_value_Arm_1"
  colnames(arm2)[ colnames(arm2) == "probability" ] <- "probability_Arm_2"
  colnames(arm2)[ colnames(arm2) == "p_value" ] <- "p_value_Arm_2"
  output <- merge( arm1, arm2, by = key_columns, all = TRUE, sort = FALSE )
  numeric_columns <- c( "probability_Arm_1", "probability_Arm_2" )
  for (column_name in numeric_columns) {
    output[[column_name]][ is.na( output[[column_name]] ) ] <- 0
  }
  for ( p_column in c( "p_value_Arm_1", "p_value_Arm_2" ) ) {
    output[[p_column]][ is.na( output[[p_column]] ) ] <- 1
  }
  output$delta_Arm2_minus_Arm1 <- output$probability_Arm_2 - output$probability_Arm_1
  epsilon <- max( c( output$probability_Arm_1, output$probability_Arm_2 ), na.rm = TRUE ) * 1e-6
  if ( !is.finite(epsilon) || epsilon <= 0 ) {
    epsilon <- 1e-12
  }
  output$log2_ratio_Arm2_over_Arm1 <- log2( ( output$probability_Arm_2 + epsilon ) / ( output$probability_Arm_1 + epsilon ) )
  output$dominant_branch <- ifelse( output$delta_Arm2_minus_Arm1 > 0, "Arm_2", ifelse( output$delta_Arm2_minus_Arm1 < 0, "Arm_1", "Tie" ) )
  output$significant_Arm_1 <- is.finite( output$p_value_Arm_1 ) & output$p_value_Arm_1 <= CELLCHAT_P_THRESHOLD
  output$significant_Arm_2 <- is.finite( output$p_value_Arm_2 ) & output$p_value_Arm_2 <= CELLCHAT_P_THRESHOLD
  output$significant_either <- output$significant_Arm_1 | output$significant_Arm_2
  output$minimum_fold_change_pass <- abs( output$log2_ratio_Arm2_over_Arm1 ) >= log2( MIN_BRANCH_FOLD_CHANGE )
  output <- output[ order( !output$significant_either, -abs( output$delta_Arm2_minus_Arm1 ), -abs( output$log2_ratio_Arm2_over_Arm1 ) ), , drop = FALSE ]
  output
}
expand_database_entities <- function( entities, complex_table ) {
  entities <- unique( as.character( entities ) )
  entities <- entities[ !is.na( entities ) & nzchar( entities ) ]
  output <- vector( "list", length( entities ) )
  for (entity_index in seq_along(entities)) {
    entity <- entities[ entity_index ]
    components <- entity
    if ( !is.null( complex_table ) && nrow( complex_table ) > 0L && entity %in% rownames( complex_table ) ) {
      components <- as.character( unlist( complex_table[ entity, , drop = FALSE ], use.names = FALSE ) )
      components <- components[ !is.na( components ) & nzchar( components ) ]
      if (length(components) == 0L) {
        components <- entity
      }
    }
    output[[entity_index]] <- data.frame( database_entity = entity, gene = unique( components ), stringsAsFactors = FALSE )
  }
  if (length(output) == 0L) {
    return( data.frame( database_entity = character(0), gene = character(0), stringsAsFactors = FALSE ) )
  }
  unique( do.call( rbind, output ) )
}
entity_components <- function( entity, complex_table, expression_genes ) {
  components <- entity
  if ( entity %in% rownames( complex_table ) ) {
    components <- as.character( unlist( complex_table[ entity, , drop = FALSE ], use.names = FALSE ) )
    components <- components[ !is.na( components ) & nzchar( components ) ]
  }
  intersect( unique( components ), expression_genes )
}
entity_expression_summary <- function( entity, cells, expression_data, complex_table ) {
  cells <- intersect( cells, colnames( expression_data ) )
  components <- entity_components( entity, complex_table, rownames( expression_data ) )
  if ( length(cells) == 0L || length(components) == 0L ) {
    return( c( entity_mean_expression = NA_real_, entity_percent_expressing = NA_real_, n_cells = length(cells), n_components = length(components) ) )
  }
  current_matrix <- expression_data[ components, cells, drop = FALSE ]
  component_means <- as.numeric( Matrix::rowMeans( current_matrix ) )
  component_percent <- as.numeric( Matrix::rowSums( current_matrix > 0 ) / length( cells ) * 100 )
  c( entity_mean_expression = min( component_means ), entity_percent_expressing = min( component_percent ), n_cells = length( cells ), n_components = length( components ) )
}
collapse_unique <- function(x) {
  values <- unique( as.character( x ) )
  values <- values[ !is.na( values ) & nzchar( values ) ]
  if (length(values) == 0L) {
    return(NA_character_)
  }
  paste( sort( values ), collapse = " | " )
}
stratified_sample_for_fast_cellchat <- function( metadata, max_cells_per_stratum, keep_all_groups, seed ) {
  required_columns <- c( "cell", "sample_id", "timepoint_standard", "cellchat_group" )
  missing_columns <- setdiff( required_columns, colnames(metadata) )
  if (length(missing_columns) > 0L) {
    stop( paste0( "Fast sampling metadata is missing: ", paste( missing_columns, collapse = ", " ) ), call. = FALSE )
  }
  metadata <- metadata[ order( metadata$timepoint_standard, metadata$sample_id, metadata$cellchat_group, metadata$cell ), , drop = FALSE ]
  stratum_key <- paste( metadata$timepoint_standard, metadata$sample_id, metadata$cellchat_group, sep = "||" )
  stratum_levels <- unique(stratum_key)
  selected_indices <- integer(0)
  audit_rows <- vector( "list", length(stratum_levels) )
  for (stratum_index in seq_along(stratum_levels)) {
    current_key <- stratum_levels[[stratum_index]]
    current_indices <- which( stratum_key == current_key )
    first_index <- current_indices[[1L]]
    current_group <- metadata$cellchat_group[[first_index]]
    keep_all <- current_group %in% keep_all_groups
    target_n <- if (keep_all) {
      length(current_indices)
    } else {
      min( length(current_indices), max_cells_per_stratum )
    }
    if (target_n < length(current_indices)) {
      set.seed( seed + stratum_index )
      current_selected <- sample( current_indices, size = target_n, replace = FALSE )
    } else {
      current_selected <- current_indices
    }
    selected_indices <- c( selected_indices, current_selected )
    audit_rows[[stratum_index]] <- data.frame( timepoint = metadata$timepoint_standard[[first_index]], sample_id = metadata$sample_id[[first_index]], cellchat_group = current_group, keep_all_group = keep_all, full_n = length(current_indices), sampled_n = length(current_selected), retained_fraction = length(current_selected) / length(current_indices), stringsAsFactors = FALSE )
  }
  sampled_metadata <- metadata[ sort(selected_indices), , drop = FALSE ]
  sampled_metadata <- sampled_metadata[ order(sampled_metadata$cell), , drop = FALSE ]
  audit <- do.call( rbind, audit_rows )
  summary_by_time <- data.table::as.data.table(audit)[ , .( full_n = sum(full_n), sampled_n = sum(sampled_n), retained_fraction = sum(sampled_n) / sum(full_n), n_strata = .N, n_keep_all_strata = sum(keep_all_group) ), by = timepoint ]
  list( metadata = sampled_metadata, audit = audit, summary_by_time = as.data.frame( summary_by_time, stringsAsFactors = FALSE ) )
}
run_cellchat_timepoint <- function( timepoint, expression_data, metadata, eligible_groups, database ) {
  timepoint_start <- Sys.time()
  current_metadata <- metadata[ metadata$timepoint_standard == timepoint & metadata$cellchat_group %in% eligible_groups, , drop = FALSE ]
  current_metadata <- current_metadata[ !duplicated( current_metadata$cell ), , drop = FALSE ]
  current_metadata <- current_metadata[ order( current_metadata$cell ), , drop = FALSE ]
  current_metadata$cellchat_group <- droplevels( factor( current_metadata$cellchat_group ) )
  current_cells <- as.character( current_metadata$cell )
  current_groups <- levels( current_metadata$cellchat_group )
  checkpoint_file <- checkpoint_path( timepoint )
  if ( REUSE_VALID_CHECKPOINTS && file.exists( checkpoint_file ) ) {
    checkpoint <- tryCatch( readRDS( checkpoint_file ), error = function(e) {
        NULL
      }
    )
    if ( !is.null( checkpoint ) && checkpoint_is_valid( checkpoint, timepoint, current_cells, current_groups ) ) {
      return( list( cellchat = checkpoint$cellchat, run_summary = data.frame( timepoint = timepoint, status = "REUSED_CHECKPOINT", n_cells = length( current_cells ), n_groups = length( current_groups ), groups = paste( current_groups, collapse = " | " ), elapsed_seconds = 0, checkpoint = safe_path( checkpoint_file ), stringsAsFactors = FALSE ) ) )
    }
  }
  current_expression <- expression_data[ , current_cells, drop = FALSE ]
  current_meta_for_cellchat <- current_metadata[ , c( "cell", "cellchat_group" ), drop = FALSE ]
  rownames( current_meta_for_cellchat ) <- current_meta_for_cellchat$cell
  cellchat <- CellChat::createCellChat( object = current_expression, meta = current_meta_for_cellchat, group.by = "cellchat_group" )
  cellchat@DB <- database
  cellchat <- CellChat::subsetData( cellchat )
  cellchat <- CellChat::identifyOverExpressedGenes( cellchat )
  cellchat <- CellChat::identifyOverExpressedInteractions( cellchat )
  compute_arguments <- list( object = cellchat, type = CELLCHAT_PARAMETERS$probability_type, trim = CELLCHAT_PARAMETERS$trim, raw.use = CELLCHAT_PARAMETERS$raw.use, population.size = CELLCHAT_PARAMETERS$population.size, distance.use = CELLCHAT_PARAMETERS$distance.use, nboot = CELLCHAT_PARAMETERS$nboot, seed.use = CELLCHAT_PARAMETERS$seed.use )
  cellchat <- call_with_supported_arguments( CellChat::computeCommunProb, compute_arguments )
  cellchat <- CellChat::filterCommunication( cellchat, min.cells = CELLCHAT_PARAMETERS$min.cells )
  cellchat <- CellChat::computeCommunProbPathway( cellchat )
  cellchat <- CellChat::aggregateNet( cellchat )
  cellchat <- CellChat::netAnalysis_computeCentrality( cellchat, slot.name = "netP" )
  elapsed_seconds <- as.numeric( difftime( Sys.time(), timepoint_start, units = "secs" ) )
  checkpoint <- list( step = "Step8B", generated_at = timestamp(), timepoint = timepoint, parameters = CELLCHAT_PARAMETERS, cells = current_cells, groups = current_groups, cellchat = cellchat )
  saveRDS( checkpoint, checkpoint_file, compress = TRUE )
  rm( current_expression )
  gc()
  list( cellchat = cellchat, run_summary = data.frame( timepoint = timepoint, status = "COMPUTED", n_cells = length( current_cells ), n_groups = length( current_groups ), groups = paste( current_groups, collapse = " | " ), elapsed_seconds = elapsed_seconds, checkpoint = safe_path( checkpoint_file ), stringsAsFactors = FALSE ) )
}
run_step8B <- function() {
  overall_start <- Sys.time()
  required_files <- c( ATLAS_RDS_FILE, STEP8A_MANIFEST_FILE, STEP8A_DECISION_FILE, STEP8A_CANDIDATE_LIGAND_FILE )
  missing_files <- required_files[ !file.exists( required_files ) ]
  input_audit <- data.frame( role = c( "Full GSE162610 atlas", "Step8A input manifest", "Step8A final decision", "Step8A Step7G candidate-ligand catalog" ), path = safe_path( required_files ), exists = file.exists( required_files ), size_bytes = ifelse( file.exists( required_files ), file.info( required_files )$size, NA_real_ ), stringsAsFactors = FALSE )
  write_csv_utf8( input_audit, INPUT_AUDIT_FILE )
  if (length(missing_files) > 0L) {
    abort_step( paste0( "Required Step8B input file(s) are missing:\n", paste( safe_path( missing_files ), collapse = "\n" ) ) )
  }
  step8a_decision <- data.table::fread( STEP8A_DECISION_FILE, data.table = FALSE, showProgress = FALSE )
  step8a_pass <- decision_pass( step8a_decision, "Overall Step 8A" )
  if (!step8a_pass) {
    abort_step( "Step8A does not contain an Overall PASS decision." )
  }
  manifest <- readRDS( STEP8A_MANIFEST_FILE )
  required_manifest_fields <- c( "status", "cell_metadata", "eligible_groups_by_time", "formal_7dpi_sender_groups", "core_recipient_groups", "Step8B_parameters" )
  missing_manifest_fields <- setdiff( required_manifest_fields, names( manifest ) )
  if (length(missing_manifest_fields) > 0L) {
    abort_step( paste0( "Step8A manifest is missing: ", paste( missing_manifest_fields, collapse = ", " ) ) )
  }
  if ( !identical( manifest$status, "PASS" ) ) {
    abort_step( "Step8A manifest status is not PASS." )
  }
  manifest_parameters <- manifest$Step8B_parameters
  parameter_names_to_compare <- setdiff( intersect( names( CELLCHAT_PARAMETERS ), names( manifest_parameters ) ), "nboot" )
  parameter_match <- all( vapply( parameter_names_to_compare, function(parameter_name) {
        identical( CELLCHAT_PARAMETERS[[parameter_name]], manifest_parameters[[parameter_name]] )
      },
      FUN.VALUE = logical(1)
    )
  )
  if (!parameter_match) {
    abort_step( paste0( "Fast Step8B structural parameters differ from the frozen Step8A ", "manifest. Only the documented nboot override is allowed." ) )
  }
  manifest_nboot <- as.integer( manifest_parameters$nboot )
  if ( !identical( manifest_nboot, 100L ) || !identical( as.integer( CELLCHAT_PARAMETERS$nboot ), 20L ) || !identical( FAST_NBOOT, 20L ) ) {
    abort_step( paste0( "Unexpected nboot configuration. Step8A should specify 100 and the ", "rapid-preview script should specify ", FAST_NBOOT, "." ) )
  }
  metadata <- as.data.frame( manifest$cell_metadata, stringsAsFactors = FALSE )
  required_metadata_columns <- c( "cell", "sample_id", "timepoint_standard", "cellchat_group", "include_step8_timecourse" )
  missing_metadata_columns <- setdiff( required_metadata_columns, colnames( metadata ) )
  if (length(missing_metadata_columns) > 0L) {
    abort_step( paste0( "Step8A cell metadata is missing: ", paste( missing_metadata_columns, collapse = ", " ) ) )
  }
  metadata$cell <- as.character( metadata$cell )
  metadata$sample_id <- as.character( metadata$sample_id )
  metadata$timepoint_standard <- as.character( metadata$timepoint_standard )
  metadata$cellchat_group <- as.character( metadata$cellchat_group )
  metadata$include_step8_timecourse <- as_logical_safe( metadata$include_step8_timecourse )
  metadata <- metadata[ metadata$include_step8_timecourse %in% TRUE & metadata$timepoint_standard %in% FORMAL_TIMEPOINTS, , drop = FALSE ]
  metadata_full <- metadata
  fast_sampling <- stratified_sample_for_fast_cellchat( metadata = metadata_full, max_cells_per_stratum = FAST_MAX_CELLS_PER_SAMPLE_GROUP, keep_all_groups = FAST_KEEP_ALL_GROUPS, seed = FAST_SAMPLING_SEED )
  metadata_inference <- fast_sampling$metadata
  write_csv_utf8( fast_sampling$audit, FAST_SAMPLING_AUDIT_FILE )
  write_csv_utf8( fast_sampling$summary_by_time, FAST_SAMPLING_SUMMARY_FILE )
  if (anyDuplicated(metadata$cell) > 0L) {
    abort_step( "Duplicated cells were found in the Step8A CellChat metadata." )
  }
  atlas <- readRDS( ATLAS_RDS_FILE )
  if (!inherits(atlas, "Seurat")) {
    abort_step( "The full atlas input is not a Seurat object." )
  }
  Seurat::DefaultAssay( atlas ) <- "RNA"
  expression_data <- get_rna_data( atlas )
  missing_expression_cells <- setdiff( unique( c( metadata_full$cell, metadata_inference$cell ) ), colnames( expression_data ) )
  if (length(missing_expression_cells) > 0L) {
    abort_step( paste0( "Normalized RNA data are missing ", length( missing_expression_cells ), " Step8A cells." ) )
  }
  metadata_full <- metadata_full[ order( metadata_full$cell ), , drop = FALSE ]
  metadata_inference <- metadata_inference[ order( metadata_inference$cell ), , drop = FALSE ]
  rm( atlas )
  gc()
  CellChatDB_mouse <- load_cellchat_mouse_db()
  complex_table <- as.data.frame( CellChatDB_mouse$complex, stringsAsFactors = FALSE )
  if ( is.null( rownames( complex_table ) ) ) {
    rownames( complex_table ) <- character( nrow( complex_table ) )
  }
  if (N_WORKERS > 1L) {
    future::plan( "multisession", workers = N_WORKERS )
  } else {
    future::plan( "sequential" )
  }
  on.exit( future::plan( "sequential" ), add = TRUE )
  cellchat_objects <- list()
  run_summary_parts <- list()
  for (timepoint in FORMAL_TIMEPOINTS) {
    eligible_groups <- as.character( manifest$eligible_groups_by_time[[timepoint]] )
    eligible_groups <- eligible_groups[ !is.na( eligible_groups ) & nzchar( eligible_groups ) ]
    if (length(eligible_groups) < 2L) {
      abort_step( paste0( "Fewer than two eligible groups are available at ", timepoint, "." ) )
    }
    run_result <- run_cellchat_timepoint( timepoint = timepoint, expression_data = expression_data, metadata = metadata_inference, eligible_groups = eligible_groups, database = CellChatDB_mouse )
    cellchat_objects[[timepoint]] <- run_result$cellchat
    run_summary_parts[[timepoint]] <- run_result$run_summary
  }
  timepoint_run_summary <- do.call( rbind, run_summary_parts )
  rownames( timepoint_run_summary ) <- NULL
  write_csv_utf8( timepoint_run_summary, TIMEPOINT_RUN_SUMMARY_FILE )
  network_summary_parts <- list()
  pathway_flow_parts <- list()
  all_interaction_parts <- list()
  pathway_long_parts <- list()
  for (timepoint in FORMAL_TIMEPOINTS) {
    cellchat <- cellchat_objects[[timepoint]]
    count_matrix <- cellchat@net$count
    weight_matrix <- cellchat@net$weight
    pathway_probability <- cellchat@netP$prob
    gene_interactions_raw <- safe_subset_communication( cellchat = cellchat, slot_name = "net", sources_use = NULL, targets_use = NULL, threshold = 1 )
    gene_interactions <- standardize_interaction_columns( gene_interactions_raw, timepoint )
    network_summary_parts[[timepoint]] <- data.frame( timepoint = timepoint, n_cells = timepoint_run_summary$n_cells[ timepoint_run_summary$timepoint == timepoint ], n_groups = length( levels( cellchat@idents ) ), total_interaction_count = sum( count_matrix, na.rm = TRUE ), total_interaction_weight = sum( weight_matrix, na.rm = TRUE ), nonzero_pathways = if ( !is.null( pathway_probability ) && length( dim( pathway_probability ) ) == 3L ) {
        sum( vapply( seq_len( dim( pathway_probability )[[3L]] ), function(pathway_index) {
              sum( pathway_probability[ , , pathway_index, drop = TRUE ], na.rm = TRUE ) >
                0
            },
            FUN.VALUE = logical(1)
          )
        )
      } else {
        0L
      },
      significant_gene_level_rows = sum( gene_interactions$p_value <= CELLCHAT_P_THRESHOLD, na.rm = TRUE ), stringsAsFactors = FALSE
    )
    pathway_flow_parts[[timepoint]] <- information_flow_table( pathway_probability, timepoint )
    all_interaction_parts[[timepoint]] <- gene_interactions
    pathway_long_parts[[timepoint]] <- pathway_array_to_long( pathway_probability, timepoint )
  }
  network_summary <- do.call( rbind, network_summary_parts )
  pathway_flow <- data.table::rbindlist( pathway_flow_parts, use.names = TRUE, fill = TRUE )
  all_interactions <- data.table::rbindlist( all_interaction_parts, use.names = TRUE, fill = TRUE )
  pathway_probability_long <- data.table::rbindlist( pathway_long_parts, use.names = TRUE, fill = TRUE )
  write_csv_utf8( network_summary, NETWORK_SUMMARY_FILE )
  write_tsv( pathway_flow, PATHWAY_FLOW_FILE )
  write_tsv_gz( all_interactions, ALL_INTERACTIONS_FILE )
  merge_arguments <- list( object.list = cellchat_objects, add.names = FORMAL_TIMEPOINTS, cell.prefix = TRUE )
  merged_cellchat <- call_with_supported_arguments( CellChat::mergeCellChat, merge_arguments )
  saveRDS( merged_cellchat, MERGED_CELLCHAT_FILE, compress = TRUE )
  cellchat_7dpi <- cellchat_objects[["7dpi"]]
  seven_dpi_groups <- levels( cellchat_7dpi@idents )
  if ( !all( FORMAL_BRANCH_SENDERS %in% seven_dpi_groups ) ) {
    abort_step( paste0( "The 7 dpi CellChat object is missing formal branch sender(s): ", paste( setdiff( FORMAL_BRANCH_SENDERS, seven_dpi_groups ), collapse = ", " ) ) )
  }
  formal_targets <- setdiff( seven_dpi_groups, EXCLUDED_FORMAL_TARGETS )
  pathway_7dpi_long <- pathway_probability_long[ pathway_probability_long$timepoint == "7dpi" & pathway_probability_long$source %in% FORMAL_BRANCH_SENDERS & pathway_probability_long$target %in% formal_targets, , drop = FALSE ]
  pathway_target_wide <- data.table::dcast( data.table::as.data.table( pathway_7dpi_long ), target + pathway ~ source, value.var = "probability", fun.aggregate = sum, fill = 0 )
  if ( !"Microglia_Arm_1" %in% colnames( pathway_target_wide ) ) {
    pathway_target_wide$Microglia_Arm_1 <- 0
  }
  if ( !"Microglia_Arm_2" %in% colnames( pathway_target_wide ) ) {
    pathway_target_wide$Microglia_Arm_2 <- 0
  }
  pathway_target_wide$delta_Arm2_minus_Arm1 <- pathway_target_wide$Microglia_Arm_2 - pathway_target_wide$Microglia_Arm_1
  pathway_target_wide$dominant_branch <- ifelse( pathway_target_wide$delta_Arm2_minus_Arm1 > 0, "Arm_2", ifelse( pathway_target_wide$delta_Arm2_minus_Arm1 < 0, "Arm_1", "Tie" ) )
  pathway_target_wide <- pathway_target_wide[ order( -abs( pathway_target_wide$delta_Arm2_minus_Arm1 ) ), , drop = FALSE ]
  write_tsv( pathway_target_wide, PATHWAY_TARGET_DIFF_FILE )
  pathway_outgoing <- data.table::dcast( data.table::as.data.table( pathway_7dpi_long ), pathway ~ source, value.var = "probability", fun.aggregate = sum, fill = 0 )
  if ( !"Microglia_Arm_1" %in% colnames( pathway_outgoing ) ) {
    pathway_outgoing$Microglia_Arm_1 <- 0
  }
  if ( !"Microglia_Arm_2" %in% colnames( pathway_outgoing ) ) {
    pathway_outgoing$Microglia_Arm_2 <- 0
  }
  pathway_outgoing$delta_Arm2_minus_Arm1 <- pathway_outgoing$Microglia_Arm_2 - pathway_outgoing$Microglia_Arm_1
  pathway_outgoing$total_branch_outgoing <- pathway_outgoing$Microglia_Arm_1 + pathway_outgoing$Microglia_Arm_2
  pathway_outgoing$Arm2_fraction <- ifelse( pathway_outgoing$total_branch_outgoing > 0, pathway_outgoing$Microglia_Arm_2 / pathway_outgoing$total_branch_outgoing, NA_real_ )
  pathway_outgoing$dominant_branch <- ifelse( pathway_outgoing$delta_Arm2_minus_Arm1 > 0, "Arm_2", ifelse( pathway_outgoing$delta_Arm2_minus_Arm1 < 0, "Arm_1", "Tie" ) )
  pathway_outgoing <- pathway_outgoing[ order( -abs( pathway_outgoing$delta_Arm2_minus_Arm1 ) ), , drop = FALSE ]
  write_tsv( pathway_outgoing, PATHWAY_OUTGOING_DIFF_FILE )
  interactions_7dpi <- as.data.frame( all_interactions[ all_interactions$timepoint == "7dpi" & all_interactions$source %in% FORMAL_BRANCH_SENDERS & all_interactions$target %in% formal_targets, , drop = FALSE ], stringsAsFactors = FALSE )
  arm1_interactions <- aggregate_branch_interactions( interactions_7dpi, "Microglia_Arm_1" )
  arm2_interactions <- aggregate_branch_interactions( interactions_7dpi, "Microglia_Arm_2" )
  branch_lr_differential <- combine_branch_interactions( arm1_interactions, arm2_interactions )
  candidate_catalog <- data.table::fread( STEP8A_CANDIDATE_LIGAND_FILE, data.table = FALSE, showProgress = FALSE )
  if ( !"gene" %in% colnames( candidate_catalog ) ) {
    abort_step( "Step8A candidate-ligand catalog does not contain a gene column." )
  }
  candidate_catalog$gene <- as.character( candidate_catalog$gene )
  ligand_entities <- unique( branch_lr_differential$ligand )
  ligand_component_map <- expand_database_entities( ligand_entities, complex_table )
  candidate_component_hits <- merge( ligand_component_map, candidate_catalog, by = "gene", all = FALSE, sort = FALSE )
  if (nrow(candidate_component_hits) > 0L) {
    candidate_hit_summary <- data.table::as.data.table( candidate_component_hits )[ , .( Step7G_candidate_TFs = collapse_unique( TF ), Step7G_mechanistic_priorities = collapse_unique( mechanistic_priority ), Step7G_candidate_ligand_components = collapse_unique( gene ) ), by = .( ligand = database_entity ) ]
    candidate_hit_summary <- as.data.frame( candidate_hit_summary, stringsAsFactors = FALSE )
  } else {
    candidate_hit_summary <- data.frame( ligand = character(0), Step7G_candidate_TFs = character(0), Step7G_mechanistic_priorities = character(0), Step7G_candidate_ligand_components = character(0), stringsAsFactors = FALSE )
  }
  branch_lr_differential <- merge( branch_lr_differential, candidate_hit_summary, by = "ligand", all.x = TRUE, sort = FALSE )
  branch_lr_differential$Step7G_candidate_ligand <- !is.na( branch_lr_differential$Step7G_candidate_TFs )
  branch_lr_differential <- branch_lr_differential[ order( !branch_lr_differential$significant_either, !branch_lr_differential$Step7G_candidate_ligand, -abs( branch_lr_differential$delta_Arm2_minus_Arm1 ) ), , drop = FALSE ]
  write_tsv( branch_lr_differential, LR_DIFF_FILE )
  seven_dpi_metadata <- metadata_full[ metadata_full$timepoint_standard == "7dpi", , drop = FALSE ]
  seven_dpi_samples <- sort( unique( seven_dpi_metadata$sample_id ) )
  recipient_sample_counts <- as.data.frame( table( seven_dpi_metadata$sample_id, seven_dpi_metadata$cellchat_group ), stringsAsFactors = FALSE )
  colnames(recipient_sample_counts) <- c( "sample_id", "cellchat_group", "n_cells" )
  recipient_sample_counts <- recipient_sample_counts[ recipient_sample_counts$cellchat_group %in% formal_targets, , drop = FALSE ]
  recipient_sample_counts$eligible_target_sample <- recipient_sample_counts$n_cells >= MIN_TARGET_CELLS_PER_SAMPLE
  write_csv_utf8( recipient_sample_counts, RECIPIENT_AUDIT_FILE )
  sample_audit_candidates <- branch_lr_differential[ branch_lr_differential$significant_either & branch_lr_differential$minimum_fold_change_pass, , drop = FALSE ]
  candidate_rows_extra <- branch_lr_differential[ branch_lr_differential$Step7G_candidate_ligand & branch_lr_differential$significant_either, , drop = FALSE ]
  sample_audit_candidates <- unique( rbind( sample_audit_candidates, candidate_rows_extra ) )
  sample_audit_candidates <- sample_audit_candidates[ order( !sample_audit_candidates$Step7G_candidate_ligand, -abs( sample_audit_candidates$delta_Arm2_minus_Arm1 ) ), , drop = FALSE ]
  sample_audit_candidates <- utils::head( sample_audit_candidates, MAX_SAMPLE_AUDIT_INTERACTIONS )
  expression_cache <- new.env( parent = emptyenv() )
  get_cached_entity_summary <- function( entity, sample_id, cellchat_group ) {
    cache_key <- paste( entity, sample_id, cellchat_group, sep = "||" )
    if ( exists( cache_key, envir = expression_cache, inherits = FALSE ) ) {
      return( get( cache_key, envir = expression_cache, inherits = FALSE ) )
    }
    cells <- seven_dpi_metadata$cell[ seven_dpi_metadata$sample_id == sample_id & seven_dpi_metadata$cellchat_group == cellchat_group ]
    summary_value <- entity_expression_summary( entity = entity, cells = cells, expression_data = expression_data, complex_table = complex_table )
    assign( cache_key, summary_value, envir = expression_cache )
    summary_value
  }
  sample_audit_parts <- list()
  if (nrow(sample_audit_candidates) > 0L) {
    for ( interaction_index in seq_len( nrow( sample_audit_candidates ) ) ) {
      current_interaction <- sample_audit_candidates[ interaction_index, , drop = FALSE ]
      for (sample_id in seven_dpi_samples) {
        arm1_ligand <- get_cached_entity_summary( current_interaction$ligand, sample_id, "Microglia_Arm_1" )
        arm2_ligand <- get_cached_entity_summary( current_interaction$ligand, sample_id, "Microglia_Arm_2" )
        target_receptor <- get_cached_entity_summary( current_interaction$receptor, sample_id, current_interaction$target )
        ligand_delta <- arm2_ligand[["entity_mean_expression"]] - arm1_ligand[["entity_mean_expression"]]
        direction_match <- if ( current_interaction$dominant_branch == "Arm_2" ) {
          is.finite( ligand_delta ) && ligand_delta >
              0
        } else if (
          current_interaction$dominant_branch == "Arm_1"
        ) {
          is.finite( ligand_delta ) && ligand_delta <
              0
        } else {
          FALSE
        }
        sender_evaluable <- arm1_ligand[["n_cells"]] >= MIN_SENDER_CELLS_PER_SAMPLE && arm2_ligand[["n_cells"]] >= MIN_SENDER_CELLS_PER_SAMPLE
        target_evaluable <- target_receptor[["n_cells"]] >= MIN_TARGET_CELLS_PER_SAMPLE
        receptor_supported <- target_evaluable && is.finite( target_receptor[["entity_percent_expressing"]] ) && target_receptor[["entity_percent_expressing"]] >= MIN_RECEPTOR_PERCENT_EXPRESSED
        sample_audit_parts[[length( sample_audit_parts ) + 1L]] <- data.frame( sample_id = sample_id, target = current_interaction$target, ligand = current_interaction$ligand, receptor = current_interaction$receptor, interaction_name = current_interaction$interaction_name, pathway_name = current_interaction$pathway_name, dominant_branch = current_interaction$dominant_branch, probability_Arm_1 = current_interaction$probability_Arm_1, probability_Arm_2 = current_interaction$probability_Arm_2, delta_Arm2_minus_Arm1 = current_interaction$delta_Arm2_minus_Arm1, p_value_Arm_1 = current_interaction$p_value_Arm_1, p_value_Arm_2 = current_interaction$p_value_Arm_2, Step7G_candidate_ligand = current_interaction$Step7G_candidate_ligand, Step7G_candidate_TFs = current_interaction$Step7G_candidate_TFs, Arm_1_ligand_mean = arm1_ligand[["entity_mean_expression"]], Arm_1_ligand_percent = arm1_ligand[["entity_percent_expressing"]], Arm_1_sender_cells = arm1_ligand[["n_cells"]], Arm_2_ligand_mean = arm2_ligand[["entity_mean_expression"]], Arm_2_ligand_percent = arm2_ligand[["entity_percent_expressing"]], Arm_2_sender_cells = arm2_ligand[["n_cells"]], ligand_delta_Arm2_minus_Arm1 = ligand_delta, ligand_direction_matches_pooled = direction_match, receptor_mean_in_target = target_receptor[["entity_mean_expression"]], receptor_percent_in_target = target_receptor[["entity_percent_expressing"]], target_cells = target_receptor[["n_cells"]], sender_evaluable = sender_evaluable, target_evaluable = target_evaluable, receptor_supported = receptor_supported, sample_axis_supported = sender_evaluable && target_evaluable && direction_match && receptor_supported, stringsAsFactors = FALSE )
      }
    }
  }
  if (length(sample_audit_parts) > 0L) {
    sample_audit_long <- do.call( rbind, sample_audit_parts )
  } else {
    sample_audit_long <- data.frame( sample_id = character(0), target = character(0), ligand = character(0), receptor = character(0), interaction_name = character(0), pathway_name = character(0), dominant_branch = character(0), probability_Arm_1 = numeric(0), probability_Arm_2 = numeric(0), delta_Arm2_minus_Arm1 = numeric(0), stringsAsFactors = FALSE )
  }
  write_tsv_gz( sample_audit_long, SAMPLE_AUDIT_LONG_FILE )
  if (nrow(sample_audit_long) > 0L) {
    sample_audit_long$Step7G_candidate_TFs[ is.na( sample_audit_long$Step7G_candidate_TFs ) ] <- ""
    sample_audit_long$pathway_name[ is.na( sample_audit_long$pathway_name ) ] <- "Unassigned"
    sample_shortlist <- stats::aggregate( cbind( sender_evaluable = as.integer( sample_audit_long$sender_evaluable ), target_evaluable = as.integer( sample_audit_long$target_evaluable ), direction_match = as.integer( sample_audit_long$ligand_direction_matches_pooled ), receptor_supported = as.integer( sample_audit_long$receptor_supported ), sample_axis_supported = as.integer( sample_audit_long$sample_axis_supported ) ) ~ target + ligand + receptor + interaction_name + pathway_name + dominant_branch + probability_Arm_1 + probability_Arm_2 + delta_Arm2_minus_Arm1 + p_value_Arm_1 + p_value_Arm_2 + Step7G_candidate_ligand + Step7G_candidate_TFs, data = sample_audit_long, FUN = sum )
    colnames(sample_shortlist)[ colnames(sample_shortlist) == "sender_evaluable" ] <- "n_samples_sender_evaluable"
    colnames(sample_shortlist)[ colnames(sample_shortlist) == "target_evaluable" ] <- "n_samples_target_evaluable"
    colnames(sample_shortlist)[ colnames(sample_shortlist) == "direction_match" ] <- "n_samples_ligand_direction_match"
    colnames(sample_shortlist)[ colnames(sample_shortlist) == "receptor_supported" ] <- "n_samples_receptor_supported"
    colnames(sample_shortlist)[ colnames(sample_shortlist) == "sample_axis_supported" ] <- "n_samples_axis_supported"
    sample_shortlist$n_7dpi_samples <- length( seven_dpi_samples )
    sample_shortlist$sample_replication_status <- ifelse( sample_shortlist$n_samples_axis_supported == sample_shortlist$n_7dpi_samples, "REPLICATED_IN_ALL_7DPI_SAMPLES", ifelse( sample_shortlist$n_samples_axis_supported > 0, "PARTIAL_SAMPLE_SUPPORT", "POOLED_ONLY" ) )
    sample_shortlist <- sample_shortlist[ order( -sample_shortlist$n_samples_axis_supported, !sample_shortlist$Step7G_candidate_ligand, -abs( sample_shortlist$delta_Arm2_minus_Arm1 ) ), , drop = FALSE ]
  } else {
    sample_shortlist <- data.frame( target = character(0), ligand = character(0), receptor = character(0), interaction_name = character(0), pathway_name = character(0), dominant_branch = character(0), probability_Arm_1 = numeric(0), probability_Arm_2 = numeric(0), delta_Arm2_minus_Arm1 = numeric(0), p_value_Arm_1 = numeric(0), p_value_Arm_2 = numeric(0), Step7G_candidate_ligand = logical(0), Step7G_candidate_TFs = character(0), n_samples_sender_evaluable = integer(0), n_samples_target_evaluable = integer(0), n_samples_ligand_direction_match = integer(0), n_samples_receptor_supported = integer(0), n_samples_axis_supported = integer(0), n_7dpi_samples = integer(0), sample_replication_status = character(0), stringsAsFactors = FALSE )
  }
  sample_shortlist <- utils::head( sample_shortlist, MAX_SHORTLIST_ROWS )
  write_tsv( sample_shortlist, SAMPLE_SHORTLIST_FILE )
  candidate_integrated <- sample_shortlist[ sample_shortlist$Step7G_candidate_ligand %in% TRUE, , drop = FALSE ]
  candidate_integrated <- candidate_integrated[ order( -candidate_integrated$n_samples_axis_supported, -abs( candidate_integrated$delta_Arm2_minus_Arm1 ) ), , drop = FALSE ]
  write_tsv( candidate_integrated, CANDIDATE_INTEGRATED_FILE )
  network_plot_long <- rbind( data.frame( timepoint = network_summary$timepoint, metric = "Interaction count", value = network_summary$total_interaction_count, stringsAsFactors = FALSE ), data.frame( timepoint = network_summary$timepoint, metric = "Interaction weight", value = network_summary$total_interaction_weight, stringsAsFactors = FALSE ) )
  network_plot_long$timepoint <- factor( network_plot_long$timepoint, levels = FORMAL_TIMEPOINTS )
  p_network <- ggplot( network_plot_long, aes( x = timepoint, y = value, group = metric ) ) + geom_line( linewidth = 0.8 ) + geom_point( size = 2.5 ) + facet_wrap( ~ metric, scales = "free_y", ncol = 1 ) + labs( title = "GSE162610 CellChat network summary", x = NULL, y = NULL ) + theme_bw( base_size = 11 )
  ggsave( FIG_NETWORK_SUMMARY, p_network, width = 7, height = 7 )
  pathway_flow_matrix <- data.table::dcast( data.table::as.data.table( pathway_flow ), pathway ~ timepoint, value.var = "information_flow", fill = 0 )
  pathway_flow_matrix$total_flow <- rowSums( pathway_flow_matrix[ , FORMAL_TIMEPOINTS, with = FALSE ] )
  pathway_flow_matrix <- pathway_flow_matrix[ order( -pathway_flow_matrix$total_flow ) ]
  pathway_flow_top <- utils::head( pathway_flow_matrix, 40L )
  pathway_flow_long <- data.table::melt( pathway_flow_top[ , c( "pathway", FORMAL_TIMEPOINTS ), with = FALSE ], id.vars = "pathway", variable.name = "timepoint", value.name = "information_flow" )
  pathway_flow_long$pathway <- factor( pathway_flow_long$pathway, levels = rev( pathway_flow_top$pathway ) )
  pathway_flow_long$timepoint <- factor( pathway_flow_long$timepoint, levels = FORMAL_TIMEPOINTS )
  p_pathway_flow <- ggplot( pathway_flow_long, aes( x = timepoint, y = pathway, fill = information_flow ) ) + geom_tile( color = "white", linewidth = 0.2 ) + scale_fill_gradient( low = "white", high = "black" ) + labs( title = "Top pathway information flow across 1, 3 and 7 dpi", x = NULL, y = NULL, fill = "Information\nflow" ) + theme_bw( base_size = 9 ) + theme( panel.grid = element_blank(), axis.text.y = element_text( size = 7 ) )
  ggsave( FIG_PATHWAY_FLOW, p_pathway_flow, width = 7, height = 10, limitsize = FALSE )
  pathway_plot <- utils::head( pathway_outgoing, MAX_FIGURE_ROWS )
  pathway_plot$pathway <- factor( pathway_plot$pathway, levels = rev( pathway_plot$pathway ) )
  p_branch_pathway <- ggplot( pathway_plot, aes( x = delta_Arm2_minus_Arm1, y = pathway ) ) + geom_vline( xintercept = 0, linetype = 2 ) + geom_point( aes( size = total_branch_outgoing, shape = dominant_branch ) ) + labs( title = "7 dpi branch-differential outgoing pathways", subtitle = "Positive values indicate stronger Arm_2 outgoing probability", x = "Arm_2 minus Arm_1 outgoing probability", y = NULL, size = "Total branch\noutgoing", shape = "Dominant" ) + theme_bw( base_size = 10 )
  ggsave( FIG_BRANCH_PATHWAY, p_branch_pathway, width = 8, height = 10, limitsize = FALSE )
  lr_plot <- branch_lr_differential[ branch_lr_differential$significant_either, , drop = FALSE ]
  lr_plot <- utils::head( lr_plot, MAX_FIGURE_ROWS )
  lr_plot$axis_label <- paste0( lr_plot$interaction_name, " -> ", lr_plot$target )
  lr_plot$axis_label <- factor( lr_plot$axis_label, levels = rev( lr_plot$axis_label ) )
  p_branch_lr <- ggplot( lr_plot, aes( x = delta_Arm2_minus_Arm1, y = axis_label ) ) + geom_vline( xintercept = 0, linetype = 2 ) + geom_point( aes( size = probability_Arm_1 + probability_Arm_2, shape = Step7G_candidate_ligand ) ) + labs( title = "Top 7 dpi branch-differential ligand-receptor interactions", subtitle = "Positive values indicate stronger Arm_2 sender probability", x = "Arm_2 minus Arm_1 probability", y = NULL, size = "Total\nprobability", shape = "Step7G\ncandidate ligand" ) + theme_bw( base_size = 9 ) + theme( axis.text.y = element_text( size = 6.5 ) )
  ggsave( FIG_BRANCH_LR, p_branch_lr, width = 9, height = 11, limitsize = FALSE )
  candidate_figure_created <- FALSE
  if (nrow(candidate_integrated) > 0L) {
    candidate_plot <- utils::head( candidate_integrated, MAX_FIGURE_ROWS )
    candidate_plot$axis_label <- paste0( candidate_plot$Step7G_candidate_TFs, ": ", candidate_plot$interaction_name, " -> ", candidate_plot$target )
    candidate_plot$axis_label <- factor( candidate_plot$axis_label, levels = rev( candidate_plot$axis_label ) )
    p_candidate <- ggplot( candidate_plot, aes( x = delta_Arm2_minus_Arm1, y = axis_label ) ) + geom_vline( xintercept = 0, linetype = 2 ) + geom_point( aes( size = probability_Arm_1 + probability_Arm_2, shape = sample_replication_status ) ) + labs( title = "Step7G candidate TF-target ligands supported by 7 dpi CellChat", subtitle = "This integrates association layers and does not prove TF-to-ligand causality", x = "Arm_2 minus Arm_1 probability", y = NULL, size = "Total\nprobability", shape = "Sample support" ) + theme_bw( base_size = 9 ) + theme( axis.text.y = element_text( size = 6.2 ) )
    ggsave( FIG_CANDIDATE_AXES, p_candidate, width = 10, height = 11, limitsize = FALSE )
    candidate_figure_created <- TRUE
  } else {
    grDevices::pdf( FIG_CANDIDATE_AXES, width = 8, height = 5 )
    graphics::plot.new()
    graphics::text( 0.5, 0.5, "No Step7G candidate ligand passed the current pooled CellChat and sample-audit filters." )
    grDevices::dev.off()
    candidate_figure_created <- TRUE
  }
  checkpoint_pass <- all( vapply( FORMAL_TIMEPOINTS, function(timepoint) {
        file.exists( checkpoint_path( timepoint ) )
      },
      FUN.VALUE = logical(1)
    )
  )
  interaction_pass <- nrow( all_interactions ) > 0L && all( FORMAL_TIMEPOINTS %in% unique( all_interactions$timepoint ) )
  pathway_pass <- nrow( pathway_flow ) > 0L && all( FORMAL_TIMEPOINTS %in% unique( pathway_flow$timepoint ) )
  branch_pathway_pass <- nrow( pathway_outgoing ) > 0L
  branch_lr_pass <- nrow( branch_lr_differential ) > 0L && any( branch_lr_differential$significant_either )
  sample_audit_pass <- file.exists( SAMPLE_AUDIT_LONG_FILE ) && file.exists( SAMPLE_SHORTLIST_FILE )
  figure_pass <- all( file.exists( c( FIG_NETWORK_SUMMARY, FIG_PATHWAY_FLOW, FIG_BRANCH_PATHWAY, FIG_BRANCH_LR, FIG_CANDIDATE_AXES ) ) )
  decision_rows <- list()
  decision_rows <- add_decision_row( decision_rows, "Required input files", length( missing_files ) == 0L, detail = "4/4 required files present" )
  decision_rows <- add_decision_row( decision_rows, "Step 8A Overall PASS", step8a_pass, detail = safe_path( STEP8A_DECISION_FILE ) )
  fast_sampling_pass <- all( FORMAL_TIMEPOINTS %in% fast_sampling$summary_by_time$timepoint ) && nrow( metadata_inference ) <
      nrow( metadata_full ) && all( FAST_KEEP_ALL_GROUPS[ FAST_KEEP_ALL_GROUPS %in% metadata_full$cellchat_group ] %in% metadata_inference$cellchat_group )
  decision_rows <- add_decision_row( decision_rows, "Fast stratified sampling integrity", fast_sampling_pass, detail = paste( paste0( fast_sampling$summary_by_time$timepoint, "=", fast_sampling$summary_by_time$sampled_n, "/", fast_sampling$summary_by_time$full_n, " (", round( 100 * fast_sampling$summary_by_time$retained_fraction, 1 ), "%)" ), collapse = "; " ) )
  decision_rows <- add_decision_row( decision_rows, "1 dpi, 3 dpi and 7 dpi CellChat checkpoints", checkpoint_pass, detail = paste( paste0( FORMAL_TIMEPOINTS, "=", file.exists( vapply( FORMAL_TIMEPOINTS, checkpoint_path, FUN.VALUE = character(1) ) ) ), collapse = "; " ) )
  decision_rows <- add_decision_row( decision_rows, "Per-timepoint interaction tables", interaction_pass, detail = paste0( nrow( all_interactions ), " rows" ) )
  decision_rows <- add_decision_row( decision_rows, "Per-timepoint pathway information flow", pathway_pass, detail = paste0( nrow( pathway_flow ), " pathway-time rows" ) )
  decision_rows <- add_decision_row( decision_rows, "Merged 1/3/7 dpi CellChat object", file.exists( MERGED_CELLCHAT_FILE ) && inherits( merged_cellchat, "CellChat" ), detail = safe_path( MERGED_CELLCHAT_FILE ) )
  decision_rows <- add_decision_row( decision_rows, "Formal 7 dpi Arm_1 versus Arm_2 pathway comparison", branch_pathway_pass, detail = paste0( nrow( pathway_outgoing ), " pathways" ) )
  decision_rows <- add_decision_row( decision_rows, "Formal 7 dpi Arm_1 versus Arm_2 LR comparison", branch_lr_pass, detail = paste0( nrow( branch_lr_differential ), " LR-target rows; significant_either=", sum( branch_lr_differential$significant_either ) ) )
  decision_rows <- add_decision_row( decision_rows, "Two-sample ligand-receptor expression audit", sample_audit_pass, detail = paste0( nrow( sample_audit_long ), " sample-interaction rows" ) )
  decision_rows <- add_decision_row( decision_rows, "Step7G candidate-ligand integration", file.exists( CANDIDATE_INTEGRATED_FILE ), detail = paste0( nrow( candidate_integrated ), " integrated rows" ) )
  decision_rows <- add_decision_row( decision_rows, "Step8B figures", figure_pass && candidate_figure_created, detail = "5/5 figures produced" )
  decision_rows <- add_decision_row( decision_rows, "Complete Step8B result RDS", FALSE, detail = "pending" )
  mandatory_pre_rds <- all( vapply( decision_rows[ seq_len( length( decision_rows ) - 1L ) ], function(current_row) {
        isTRUE( current_row$pass[[1L]] )
      },
      FUN.VALUE = logical(1)
    )
  )
  preliminary_status <- if ( mandatory_pre_rds ) {
    "PASS"
  } else {
    "FAIL"
  }
  results_object <- list( step = "Step8B", generated_at = timestamp(), preliminary_status = preliminary_status, parameters = CELLCHAT_PARAMETERS, discovery_mode = FAST_DISCOVERY_MODE, analysis_label = "OFFICIAL_RAPID_PREVIEW_NBOOT20", confirmation_required = TRUE, fast_sampling_parameters = list( max_cells_per_sample_group = FAST_MAX_CELLS_PER_SAMPLE_GROUP, keep_all_groups = FAST_KEEP_ALL_GROUPS, sampling_seed = FAST_SAMPLING_SEED ), fast_sampling_audit = fast_sampling$audit, fast_sampling_summary = fast_sampling$summary_by_time, full_input_cells = nrow( metadata_full ), inference_cells = nrow( metadata_inference ), workers = N_WORKERS, timepoint_run_summary = timepoint_run_summary, network_summary = network_summary, pathway_information_flow = pathway_flow, pathway_target_differential_7dpi = pathway_target_wide, pathway_outgoing_differential_7dpi = pathway_outgoing, LR_differential_7dpi = branch_lr_differential, sample_audit_shortlist = sample_shortlist, Step7G_candidate_integrated = candidate_integrated, recipient_sample_coverage = recipient_sample_counts, checkpoint_paths = setNames( vapply( FORMAL_TIMEPOINTS, function(timepoint) {
          safe_path( checkpoint_path( timepoint ) )
        },
        FUN.VALUE = character(1)
      ), FORMAL_TIMEPOINTS
    ), merged_cellchat_path = safe_path(
        MERGED_CELLCHAT_FILE
      ), interpretation_boundary = c(
      "Arm_1 and Arm_2 remain neutral labels.", "CellChat bootstrap P values are not replicate-level P values.", "Sample expression support is descriptive replication evidence.", "Step7G TF-target-ligand integration does not establish causality."
    )
  )
  saveRDS( results_object, RESULT_RDS_FILE, compress = TRUE )
  result_rds_pass <- file.exists( RESULT_RDS_FILE )
  decision_rows[[length( decision_rows )]]$pass <- result_rds_pass
  decision_rows[[length( decision_rows )]]$status <- if ( result_rds_pass ) {
    "PASS"
  } else {
    "FAIL"
  }
  decision_rows[[length( decision_rows )]]$detail <- safe_path( RESULT_RDS_FILE )
  final_decision <- do.call( rbind, decision_rows )
  overall_pass <- all( final_decision$pass )
  final_decision <- rbind( final_decision, data.frame( domain = "Overall Step 8B PREVIEW", pass = overall_pass, status = if ( overall_pass ) {
        "PASS"
      } else {
        "FAIL"
      },
      detail = if ( overall_pass ) {
        "PASS"
      } else {
        "One or more mandatory Step8B domains failed."
      },
      stringsAsFactors = FALSE
    )
  )
  write_csv_utf8( final_decision, FINAL_DECISION_FILE )
  FINAL_STATUS <- if ( overall_pass ) {
    "PASS"
  } else {
    "FAIL"
  }
  results_object$status <- FINAL_STATUS
  results_object$final_decision <- final_decision
  saveRDS( results_object, RESULT_RDS_FILE, compress = TRUE )
  overall_elapsed <- as.numeric( difftime( Sys.time(), overall_start, units = "secs" ) )
  top_arm2_pathways <- pathway_outgoing[ pathway_outgoing$dominant_branch == "Arm_2", , drop = FALSE ]
  top_arm2_pathways <- utils::head( top_arm2_pathways, 20L )
  top_arm1_pathways <- pathway_outgoing[ pathway_outgoing$dominant_branch == "Arm_1", , drop = FALSE ]
  top_arm1_pathways <- top_arm1_pathways[ order( top_arm1_pathways$delta_Arm2_minus_Arm1 ), , drop = FALSE ]
  top_arm1_pathways <- utils::head( top_arm1_pathways, 20L )
  top_replicated_lr <- sample_shortlist[ sample_shortlist$sample_replication_status == "REPLICATED_IN_ALL_7DPI_SAMPLES", , drop = FALSE ]
  top_replicated_lr <- utils::head( top_replicated_lr, 30L )
  capture.output( sessionInfo(), file = SESSION_FILE )
  future::plan( "sequential" )
  gc()
  invisible( results_object )
}
STEP8B_RESULT <- tryCatch( run_step8B(), error = function(e) {
    future::plan( "sequential" )
    error_lines <- c( "======================================================================", "GSE162610 | STEP 8B PREVIEW ERROR", "======================================================================", paste0( "Time: ", timestamp() ), paste0( "Message: ", conditionMessage(e) ), paste0( "Failed call: ", paste( deparse( conditionCall(e) ), collapse = " " ) ), "", paste0( "Valid completed timepoint checkpoints are retained in: ", safe_path( CHECKPOINTDIR ) ), "rerun Step1-Step8A.", "Rerun only this Step8B PREVIEW file after correcting the reported issue." )
    writeLines( error_lines, con = ERROR_FILE, useBytes = TRUE )
    stop( e )
  }
)
