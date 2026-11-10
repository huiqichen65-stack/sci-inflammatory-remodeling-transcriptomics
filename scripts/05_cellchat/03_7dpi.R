rm(list = ls())
gc()
options( stringsAsFactors = FALSE, width = 260, warn = 1, timeout = 86400, future.globals.maxSize = 32 * 1024^3, future.rng.onMisuse = "error" )
set.seed(20260803L)
PROJECT_ROOT <- Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
FORMAL_TIMEPOINT <- "7dpi"
FORMAL_BRANCH_SENDERS <- c( "Microglia_Arm_1", "Microglia_Arm_2" )
EXCLUDED_FORMAL_TARGETS <- c( "Microglia_Arm_1", "Microglia_Arm_2", "Microglia_Shared_prebranch" )
CELLCHAT_P_THRESHOLD <- 0.05
MIN_BRANCH_FOLD_CHANGE <- 1.25
MIN_CELLS_PER_GROUP <- 20L
MIN_SENDER_CELLS_PER_SAMPLE <- 20L
MIN_TARGET_CELLS_PER_SAMPLE <- 20L
MIN_RECEPTOR_PERCENT_EXPRESSED <- 5
NBOOT_CONFIRMATION <- 100L
SEED_CONFIRMATION <- 20260803L
REUSE_VALID_CHECKPOINT <- TRUE
PRIMARY_INTERACTIONS <- c( "TNF_TNFRSF1A", "TNF_TNFRSF1B", "ICAM1_ITGAM_ITGB2", "ICAM1_ITGAX_ITGB2", "ICAM1_ITGAM_ITGB2L", "CCL3_CCR1", "CCL2_CCR2", "GAS6_AXL", "GAS6_MERTK" )
CONTEXT_INTERACTIONS <- c( "TGFB1_TGFBR1_TGFBR2", "SPP1_CD44", "SPP1_ITGAV_ITGB1", "SPP1_ITGAV_ITGB5", "SPP1_ITGA5_ITGB1", "OSM_OSMR_IL6ST", "OSM_LIFR_IL6ST", "PPIA_BSG", "ENTPD1_ADORA1", "ENTPD1_ADORA2B", "GAS6_TYRO3" )
CELLCHAT_PARAMETERS <- list( database = "CellChatDB.mouse_targeted_preview_shortlist", probability_type = "triMean", trim = 0.10, raw.use = TRUE, population.size = FALSE, distance.use = FALSE, min.cells = MIN_CELLS_PER_GROUP, nboot = NBOOT_CONFIRMATION, seed.use = SEED_CONFIRMATION, workers = 1L )
ATLAS_RDS_FILE <- file.path( PROJECT_ROOT, "04_results", "step1F_GSE162610_broad_annotation", "GSE162610_broad_annotation_audit.rds" )
STEP8A_DIR <- file.path( PROJECT_ROOT, "04_results", "step8A_GSE162610_CellChat_input_preflight" )
STEP8A_MANIFEST_FILE <- file.path( STEP8A_DIR, "RDS", "GSE162610_step8A_CellChat_input_manifest.rds" )
STEP8A_DECISION_FILE <- file.path( STEP8A_DIR, "tables", "08_final_decision.csv" )
PREVIEW_DIR <- file.path( PROJECT_ROOT, "04_results", "step8B_PREVIEW_GSE162610_CellChat_nboot20" )
PREVIEW_DECISION_FILE <- file.path( PREVIEW_DIR, "tables", "12_final_decision.csv" )
PREVIEW_LR_FILE <- file.path( PREVIEW_DIR, "tables", "07_7dpi_branch_LR_differential.tsv" )
PREVIEW_SAMPLE_SHORTLIST_FILE <- file.path( PREVIEW_DIR, "tables", "09_7dpi_sample_supported_LR_shortlist.tsv" )
OUTDIR <- file.path( PROJECT_ROOT, "04_results", "step8B_CONFIRM_GSE162610_7dpi_targeted_fullcell_nboot100" )
TABDIR <- file.path(OUTDIR, "tables")
FIGDIR <- file.path(OUTDIR, "figures")
RDSDIR <- file.path(OUTDIR, "RDS")
REPORTDIR <- file.path(OUTDIR, "reports")
for (directory in c( OUTDIR, TABDIR, FIGDIR, RDSDIR, REPORTDIR )) {
  dir.create( directory, recursive = TRUE, showWarnings = FALSE )
}
INPUT_AUDIT_FILE <- file.path( TABDIR, "00_input_file_audit.csv" )
CELL_COUNT_FILE <- file.path( TABDIR, "01_full_7dpi_cell_counts_by_group_sample.csv" )
TARGETED_LR_CATALOG_FILE <- file.path( TABDIR, "02_targeted_LR_catalog.tsv" )
ALL_COMMUNICATION_FILE <- file.path( TABDIR, "03_fullcell_targeted_all_communications.tsv" )
FULL_BRANCH_CONFIRMATION_FILE <- file.path( TABDIR, "04_fullcell_7dpi_branch_LR_confirmation.tsv" )
SAMPLE_AUDIT_LONG_FILE <- file.path( TABDIR, "05_fullcell_7dpi_sample_expression_audit_long.tsv" )
PREVIEW_VS_FULL_FILE <- file.path( TABDIR, "06_preview_vs_full_confirmation.tsv" )
CONFIRMED_CORE_FILE <- file.path( TABDIR, "07_confirmed_core_axes.tsv" )
PATHWAY_SUMMARY_FILE <- file.path( TABDIR, "08_targeted_pathway_confirmation_summary.tsv" )
FINAL_DECISION_FILE <- file.path( TABDIR, "09_final_decision.csv" )
CELLCHAT_CHECKPOINT_FILE <- file.path( RDSDIR, "GSE162610_7dpi_targeted_fullcell_nboot100_CellChat_checkpoint.rds" )
RESULT_RDS_FILE <- file.path( RDSDIR, "GSE162610_step8B_7dpi_targeted_fullcell_nboot100_confirmation_results.rds" )
REPORT_FILE <- file.path( REPORTDIR, "Step8B_7dpi_targeted_fullcell_nboot100_confirmation_report.txt" )
ERROR_FILE <- file.path( REPORTDIR, "Step8B_7dpi_targeted_fullcell_nboot100_confirmation_ERROR.txt" )
FIGURE_DELTA_FILE <- file.path( FIGDIR, "Fig8B_CONFIRM_01_branch_probability_delta.pdf" )
FIGURE_PREVIEW_FULL_FILE <- file.path( FIGDIR, "Fig8B_CONFIRM_02_preview_vs_full_delta.pdf" )
required_packages <- c( "Seurat", "SeuratObject", "Matrix", "CellChat", "data.table", "ggplot2", "future" )
missing_packages <- required_packages[ !vapply( required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1) ) ]
if (length(missing_packages) > 0L) {
  stop( paste0( "Missing required package(s): ", paste( missing_packages, collapse = ", " ) ), call. = FALSE )
}
timestamp <- function() {
  format( Sys.time(), "%Y-%m-%d %H:%M:%S" )
}
safe_path <- function(path) {
  normalizePath( path, winslash = "/", mustWork = FALSE )
}
write_csv_utf8 <- function(data, path) {
  data.table::fwrite( data, file = path, sep = ",", quote = TRUE, na = "NA" )
}
write_tsv <- function(data, path) {
  data.table::fwrite( data, file = path, sep = "\t", quote = FALSE, na = "NA" )
}
abort_step <- function(text) {
  writeLines( c( "======================================================================", "GSE162610 | STEP 8B TARGETED FULL-CELL CONFIRMATION ERROR", "======================================================================", paste0( "Time: ", timestamp() ), "", text ), con = ERROR_FILE, useBytes = TRUE )
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
  if ( !"domain" %in% colnames( decision_table ) ) {
    return(FALSE)
  }
  index <- which( as.character( decision_table$domain ) == domain_name )
  if (length(index) == 0L) {
    return(FALSE)
  }
  if ( "pass" %in% colnames( decision_table ) ) {
    value <- as_logical_safe( decision_table$pass[ index[[1L]] ] )
    return( isTRUE( value[[1L]] ) )
  }
  if ( "status" %in% colnames( decision_table ) ) {
    return( identical( toupper( as.character( decision_table$status[ index[[1L]] ] ) ), "PASS" ) )
  }
  FALSE
}
add_decision_row <- function( decision_list, domain, pass, detail = "", status_if_fail = "FAIL" ) {
  decision_list[[length( decision_list ) + 1L]] <- data.frame( domain = domain, pass = isTRUE(pass), status = if ( isTRUE(pass) ) {
      "PASS"
    } else {
      status_if_fail
    },
    detail = as.character( detail ), stringsAsFactors = FALSE
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
call_with_supported_arguments <- function( function_object, arguments ) {
  supported <- names( formals( function_object ) )
  if ( "..." %in% supported ) {
    return( do.call( function_object, arguments ) )
  }
  arguments <- arguments[ names( arguments ) %in% supported ]
  do.call( function_object, arguments )
}
safe_subset_communication <- function( cellchat, slot_name = "net", sources_use = NULL, targets_use = NULL, threshold = 1 ) {
  arguments <- list( object = cellchat, slot.name = slot_name, sources.use = sources_use, targets.use = targets_use, thresh = threshold )
  output <- call_with_supported_arguments( CellChat::subsetCommunication, arguments )
  as.data.frame( output, stringsAsFactors = FALSE )
}
standardize_interaction_columns <- function( interaction_table, timepoint ) {
  if (nrow(interaction_table) == 0L) {
    return( data.frame( timepoint = character(0), source = character(0), target = character(0), ligand = character(0), receptor = character(0), interaction_name = character(0), pathway_name = character(0), probability = numeric(0), p_value = numeric(0), stringsAsFactors = FALSE ) )
  }
  column_map <- list( source = c( "source" ), target = c( "target" ), ligand = c( "ligand" ), receptor = c( "receptor" ), interaction_name = c( "interaction_name", "interaction_name_2" ), pathway_name = c( "pathway_name", "pathway" ), probability = c( "prob", "probability" ), p_value = c( "pval", "p_value" ) )
  resolved <- list()
  for (canonical_name in names(column_map)) {
    candidates <- column_map[[canonical_name]]
    hit <- candidates[ candidates %in% colnames( interaction_table ) ]
    resolved[[canonical_name]] <- if ( length(hit) > 0L ) {
      hit[[1L]]
    } else {
      NA_character_
    }
  }
  required_resolved <- c( "source", "target", "ligand", "receptor", "probability" )
  missing_required <- required_resolved[ is.na( unlist( resolved[ required_resolved ] ) ) ]
  if (length(missing_required) > 0L) {
    stop( paste0( "subsetCommunication output is missing: ", paste( missing_required, collapse = ", " ), ". Current columns: ", paste( colnames( interaction_table ), collapse = ", " ) ), call. = FALSE )
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
  output
}
aggregate_branch_interactions <- function( interaction_table, source_group ) {
  current <- interaction_table[ interaction_table$source == source_group, , drop = FALSE ]
  if (nrow(current) == 0L) {
    return( data.frame( target = character(0), ligand = character(0), receptor = character(0), interaction_name = character(0), pathway_name = character(0), probability = numeric(0), p_value = numeric(0), stringsAsFactors = FALSE ) )
  }
  current$pathway_name[ is.na( current$pathway_name ) ] <- "Unassigned"
  aggregate_probability <- stats::aggregate( probability ~ target + ligand + receptor + interaction_name + pathway_name, data = current, FUN = sum, na.rm = TRUE )
  aggregate_p <- stats::aggregate( p_value ~ target + ligand + receptor + interaction_name + pathway_name, data = current, FUN = function(values) {
      finite_values <- values[ is.finite( values ) ]
      if (length(finite_values) == 0L) {
        return(1)
      }
      min( finite_values )
    }
  )
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
  for ( column_name in c( "probability_Arm_1", "probability_Arm_2" ) ) {
    output[[column_name]][ is.na( output[[column_name]] ) ] <- 0
  }
  for ( column_name in c( "p_value_Arm_1", "p_value_Arm_2" ) ) {
    output[[column_name]][ is.na( output[[column_name]] ) ] <- 1
  }
  output$delta_Arm2_minus_Arm1 <- output$probability_Arm_2 - output$probability_Arm_1
  epsilon <- max( c( output$probability_Arm_1, output$probability_Arm_2 ), na.rm = TRUE ) * 1e-6
  if ( !is.finite( epsilon ) || epsilon <= 0 ) {
    epsilon <- 1e-12
  }
  output$log2_ratio_Arm2_over_Arm1 <- log2( ( output$probability_Arm_2 + epsilon ) / ( output$probability_Arm_1 + epsilon ) )
  output$dominant_branch <- ifelse( output$delta_Arm2_minus_Arm1 > 0, "Arm_2", ifelse( output$delta_Arm2_minus_Arm1 < 0, "Arm_1", "Tie" ) )
  output$significant_Arm_1 <- is.finite( output$p_value_Arm_1 ) & output$p_value_Arm_1 <= CELLCHAT_P_THRESHOLD
  output$significant_Arm_2 <- is.finite( output$p_value_Arm_2 ) & output$p_value_Arm_2 <= CELLCHAT_P_THRESHOLD
  output$significant_either <- output$significant_Arm_1 | output$significant_Arm_2
  output$minimum_fold_change_pass <- abs( output$log2_ratio_Arm2_over_Arm1 ) >= log2( MIN_BRANCH_FOLD_CHANGE )
  output[ order( !output$significant_either, -abs( output$delta_Arm2_minus_Arm1 ), -abs( output$log2_ratio_Arm2_over_Arm1 ) ), , drop = FALSE ]
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
    return( c( entity_mean_expression = NA_real_, entity_percent_expressing = NA_real_, n_cells = length( cells ), n_components = length( components ) ) )
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
    return( NA_character_ )
  }
  paste( sort( values ), collapse = " | " )
}
checkpoint_is_valid <- function( checkpoint, expected_cells, expected_groups, expected_interactions ) {
  required_names <- c( "cellchat", "cells", "groups", "interactions", "parameters" )
  if ( !is.list( checkpoint ) || !all( required_names %in% names( checkpoint ) ) ) {
    return(FALSE)
  }
  if ( !inherits( checkpoint$cellchat, "CellChat" ) ) {
    return(FALSE)
  }
  identical( checkpoint$cells, expected_cells ) && identical( checkpoint$groups, expected_groups ) && identical( checkpoint$interactions, expected_interactions ) && identical( checkpoint$parameters, CELLCHAT_PARAMETERS )
}
main <- function() {
  start_time <- Sys.time()
  required_files <- c( ATLAS_RDS_FILE, STEP8A_MANIFEST_FILE, STEP8A_DECISION_FILE, PREVIEW_DECISION_FILE, PREVIEW_LR_FILE, PREVIEW_SAMPLE_SHORTLIST_FILE )
  input_audit <- data.frame( role = c( "Full atlas", "Step8A input manifest", "Step8A final decision", "Step8B preview final decision", "Step8B preview LR differential", "Step8B preview sample-supported shortlist" ), path = safe_path( required_files ), exists = file.exists( required_files ), size_bytes = ifelse( file.exists( required_files ), file.info( required_files )$size, NA_real_ ), stringsAsFactors = FALSE )
  write_csv_utf8( input_audit, INPUT_AUDIT_FILE )
  missing_files <- required_files[ !file.exists( required_files ) ]
  if (length(missing_files) > 0L) {
    abort_step( paste0( "Missing required file(s):\n", paste( safe_path( missing_files ), collapse = "\n" ) ) )
  }
  step8a_decision <- data.table::fread( STEP8A_DECISION_FILE, data.table = FALSE, showProgress = FALSE )
  preview_decision <- data.table::fread( PREVIEW_DECISION_FILE, data.table = FALSE, showProgress = FALSE )
  if ( !decision_pass( step8a_decision, "Overall Step 8A" ) ) {
    abort_step( "Step8A Overall PASS was not found." )
  }
  if ( !decision_pass( preview_decision, "Overall Step 8B PREVIEW" ) ) {
    abort_step( "Step8B preview Overall PASS was not found." )
  }
  manifest <- readRDS( STEP8A_MANIFEST_FILE )
  required_manifest_fields <- c( "status", "cell_metadata", "eligible_groups_by_time", "formal_7dpi_sender_groups" )
  missing_manifest_fields <- setdiff( required_manifest_fields, names( manifest ) )
  if (length(missing_manifest_fields) > 0L) {
    abort_step( paste0( "Step8A manifest is missing: ", paste( missing_manifest_fields, collapse = ", " ) ) )
  }
  if ( !identical( manifest$status, "PASS" ) ) {
    abort_step( "Step8A manifest status is not PASS." )
  }
  metadata <- as.data.frame( manifest$cell_metadata, stringsAsFactors = FALSE )
  required_metadata_columns <- c( "cell", "sample_id", "timepoint_standard", "cellchat_group", "include_step8_timecourse" )
  missing_metadata_columns <- setdiff( required_metadata_columns, colnames( metadata ) )
  if (length(missing_metadata_columns) > 0L) {
    abort_step( paste0( "Step8A metadata is missing: ", paste( missing_metadata_columns, collapse = ", " ) ) )
  }
  metadata$cell <- as.character( metadata$cell )
  metadata$sample_id <- as.character( metadata$sample_id )
  metadata$timepoint_standard <- as.character( metadata$timepoint_standard )
  metadata$cellchat_group <- as.character( metadata$cellchat_group )
  metadata$include_step8_timecourse <- as_logical_safe( metadata$include_step8_timecourse )
  eligible_groups <- as.character( manifest$eligible_groups_by_time[[FORMAL_TIMEPOINT]] )
  eligible_groups <- eligible_groups[ !is.na( eligible_groups ) & nzchar( eligible_groups ) ]
  metadata_7dpi <- metadata[ metadata$include_step8_timecourse %in% TRUE & metadata$timepoint_standard == FORMAL_TIMEPOINT & metadata$cellchat_group %in% eligible_groups, , drop = FALSE ]
  metadata_7dpi <- metadata_7dpi[ order( metadata_7dpi$cell ), , drop = FALSE ]
  if (anyDuplicated(metadata_7dpi$cell) > 0L) {
    abort_step( "Duplicated 7 dpi cells were found in the Step8A metadata." )
  }
  group_counts <- as.data.frame( table( metadata_7dpi$sample_id, metadata_7dpi$cellchat_group ), stringsAsFactors = FALSE )
  colnames(group_counts) <- c( "sample_id", "cellchat_group", "n_cells" )
  group_counts <- group_counts[ group_counts$n_cells > 0, , drop = FALSE ]
  write_csv_utf8( group_counts, CELL_COUNT_FILE )
  total_group_counts <- table( metadata_7dpi$cellchat_group )
  insufficient_groups <- names( total_group_counts )[ total_group_counts < MIN_CELLS_PER_GROUP ]
  if (length(insufficient_groups) > 0L) {
    abort_step( paste0( "The following 7 dpi groups have fewer than ", MIN_CELLS_PER_GROUP, " cells: ", paste( insufficient_groups, collapse = ", " ) ) )
  }
  if ( !all( FORMAL_BRANCH_SENDERS %in% unique( metadata_7dpi$cellchat_group ) ) ) {
    abort_step( "Both formal branch sender groups are not present at 7 dpi." )
  }
  preview_lr <- data.table::fread( PREVIEW_LR_FILE, data.table = FALSE, showProgress = FALSE )
  preview_shortlist <- data.table::fread( PREVIEW_SAMPLE_SHORTLIST_FILE, data.table = FALSE, showProgress = FALSE )
  required_preview_columns <- c( "target", "ligand", "receptor", "interaction_name", "pathway_name", "dominant_branch", "probability_Arm_1", "probability_Arm_2", "delta_Arm2_minus_Arm1" )
  missing_preview_columns <- setdiff( required_preview_columns, colnames( preview_lr ) )
  if (length(missing_preview_columns) > 0L) {
    abort_step( paste0( "Preview LR table is missing: ", paste( missing_preview_columns, collapse = ", " ) ) )
  }
  if ( !"interaction_name" %in% colnames( preview_shortlist ) || !"sample_replication_status" %in% colnames( preview_shortlist ) ) {
    abort_step( paste0( "Preview sample shortlist must contain interaction_name and ", "sample_replication_status." ) )
  }
  selected_interactions <- sort( unique( as.character( preview_shortlist$interaction_name ) ) )
  selected_interactions <- selected_interactions[ !is.na( selected_interactions ) & nzchar( selected_interactions ) ]
  required_primary_missing <- setdiff( PRIMARY_INTERACTIONS, selected_interactions )
  if (length(required_primary_missing) > 0L) {
    abort_step( paste0( "Primary preview interaction(s) are missing from the shortlist: ", paste( required_primary_missing, collapse = ", " ) ) )
  }
  atlas <- readRDS( ATLAS_RDS_FILE )
  if (!inherits(atlas, "Seurat")) {
    abort_step( "The full atlas input is not a Seurat object." )
  }
  Seurat::DefaultAssay( atlas ) <- "RNA"
  expression_data <- get_rna_data( atlas )
  missing_cells <- setdiff( metadata_7dpi$cell, colnames( expression_data ) )
  if (length(missing_cells) > 0L) {
    abort_step( paste0( "Normalized RNA data are missing ", length( missing_cells ), " eligible 7 dpi cells." ) )
  }
  rm( atlas )
  gc()
  CellChatDB_mouse <- load_cellchat_mouse_db()
  interaction_db <- as.data.frame( CellChatDB_mouse$interaction, stringsAsFactors = FALSE )
  if ( !"interaction_name" %in% colnames( interaction_db ) ) {
    interaction_db$interaction_name <- rownames( interaction_db )
  }
  interaction_db$interaction_name <- as.character( interaction_db$interaction_name )
  selected_lr_db <- interaction_db[ interaction_db$interaction_name %in% selected_interactions, , drop = FALSE ]
  missing_from_database <- setdiff( selected_interactions, selected_lr_db$interaction_name )
  if (length(missing_from_database) > 0L) {
    abort_step( paste0( "Selected interaction(s) are absent from CellChatDB.mouse: ", paste( missing_from_database, collapse = ", " ) ) )
  }
  selected_lr_db <- selected_lr_db[ match( selected_interactions, selected_lr_db$interaction_name ), , drop = FALSE ]
  rownames( selected_lr_db ) <- selected_lr_db$interaction_name
  targeted_database <- CellChatDB_mouse
  targeted_database$interaction <- selected_lr_db
  preview_support_by_interaction <- data.table::as.data.table( preview_shortlist )[ , .( preview_replication_status = if ( any( sample_replication_status == "REPLICATED_IN_ALL_7DPI_SAMPLES", na.rm = TRUE ) ) {
        "REPLICATED_IN_ALL_7DPI_SAMPLES"
      } else {
        "PARTIAL_SAMPLE_SUPPORT"
      },
      preview_targets = collapse_unique( target ), preview_dominant_branches = collapse_unique( dominant_branch )
    ), by = interaction_name
  ]
  preview_support_by_interaction <- as.data.frame( preview_support_by_interaction, stringsAsFactors = FALSE )
  targeted_lr_catalog <- merge( selected_lr_db, preview_support_by_interaction, by = "interaction_name", all.x = TRUE, sort = FALSE )
  targeted_lr_catalog$confirmation_tier <- ifelse( targeted_lr_catalog$interaction_name %in% PRIMARY_INTERACTIONS, "PRIMARY", ifelse( targeted_lr_catalog$interaction_name %in% CONTEXT_INTERACTIONS, "CONTEXT", "SECONDARY" ) )
  targeted_lr_catalog <- targeted_lr_catalog[ order( match( targeted_lr_catalog$confirmation_tier, c( "PRIMARY", "SECONDARY", "CONTEXT" ) ), targeted_lr_catalog$pathway_name, targeted_lr_catalog$interaction_name ), , drop = FALSE ]
  write_tsv( targeted_lr_catalog, TARGETED_LR_CATALOG_FILE )
  expected_cells <- metadata_7dpi$cell
  expected_groups <- sort( unique( metadata_7dpi$cellchat_group ) )
  expected_interactions <- sort( selected_lr_db$interaction_name )
  cellchat <- NULL
  run_status <- "COMPUTED"
  if ( REUSE_VALID_CHECKPOINT && file.exists( CELLCHAT_CHECKPOINT_FILE ) ) {
    checkpoint <- tryCatch( readRDS( CELLCHAT_CHECKPOINT_FILE ), error = function(e) {
        NULL
      }
    )
    if ( !is.null( checkpoint ) && checkpoint_is_valid( checkpoint, expected_cells, expected_groups, expected_interactions ) ) {
      cellchat <- checkpoint$cellchat
      run_status <- "REUSED_CHECKPOINT"
    }
  }
  if (is.null(cellchat)) {
    current_expression <- expression_data[ , metadata_7dpi$cell, drop = FALSE ]
    current_meta <- data.frame( cellchat_group = factor( metadata_7dpi$cellchat_group, levels = expected_groups ), samples = metadata_7dpi$sample_id, row.names = metadata_7dpi$cell, stringsAsFactors = FALSE )
    future::plan( "sequential" )
    cellchat <- CellChat::createCellChat( object = current_expression, meta = current_meta, group.by = "cellchat_group" )
    cellchat@DB <- targeted_database
    cellchat <- CellChat::subsetData( cellchat )
    if (nrow(cellchat@data.signaling) == 0L) {
      abort_step( "Targeted subsetData produced zero signaling genes." )
    }
    cellchat@LR <- list( LRsig = selected_lr_db )
    cellchat <- call_with_supported_arguments( CellChat::computeCommunProb, list( object = cellchat, type = CELLCHAT_PARAMETERS$probability_type, trim = CELLCHAT_PARAMETERS$trim, LR.use = selected_lr_db, raw.use = CELLCHAT_PARAMETERS$raw.use, population.size = CELLCHAT_PARAMETERS$population.size, distance.use = CELLCHAT_PARAMETERS$distance.use, nboot = CELLCHAT_PARAMETERS$nboot, seed.use = CELLCHAT_PARAMETERS$seed.use ) )
    cellchat <- call_with_supported_arguments( CellChat::filterCommunication, list( object = cellchat, min.cells = CELLCHAT_PARAMETERS$min.cells ) )
    cellchat <- call_with_supported_arguments( CellChat::computeCommunProbPathway, list( object = cellchat, pairLR.use = selected_lr_db, thresh = CELLCHAT_P_THRESHOLD ) )
    cellchat <- call_with_supported_arguments( CellChat::aggregateNet, list( object = cellchat, pairLR.use = selected_lr_db, thresh = CELLCHAT_P_THRESHOLD, return.object = TRUE ) )
    checkpoint <- list( cellchat = cellchat, cells = expected_cells, groups = expected_groups, interactions = expected_interactions, parameters = CELLCHAT_PARAMETERS, created = timestamp() )
    saveRDS( checkpoint, CELLCHAT_CHECKPOINT_FILE, compress = FALSE )
  }
  all_communications_raw <- safe_subset_communication( cellchat = cellchat, slot_name = "net", sources_use = NULL, targets_use = NULL, threshold = 1 )
  all_communications <- standardize_interaction_columns( all_communications_raw, FORMAL_TIMEPOINT )
  write_tsv( all_communications, ALL_COMMUNICATION_FILE )
  formal_targets <- setdiff( expected_groups, EXCLUDED_FORMAL_TARGETS )
  formal_communications <- all_communications[ all_communications$source %in% FORMAL_BRANCH_SENDERS & all_communications$target %in% formal_targets, , drop = FALSE ]
  arm1 <- aggregate_branch_interactions( formal_communications, "Microglia_Arm_1" )
  arm2 <- aggregate_branch_interactions( formal_communications, "Microglia_Arm_2" )
  full_branch_confirmation <- combine_branch_interactions( arm1, arm2 )
  full_branch_confirmation$confirmation_tier <- ifelse( full_branch_confirmation$interaction_name %in% PRIMARY_INTERACTIONS, "PRIMARY", ifelse( full_branch_confirmation$interaction_name %in% CONTEXT_INTERACTIONS, "CONTEXT", "SECONDARY" ) )
  write_tsv( full_branch_confirmation, FULL_BRANCH_CONFIRMATION_FILE )
  complex_table <- as.data.frame( CellChatDB_mouse$complex, stringsAsFactors = FALSE )
  samples_7dpi <- sort( unique( metadata_7dpi$sample_id ) )
  expression_cache <- new.env( parent = emptyenv() )
  get_cached_entity_summary <- function( entity, sample_id, cellchat_group ) {
    cache_key <- paste( entity, sample_id, cellchat_group, sep = "||" )
    if ( exists( cache_key, envir = expression_cache, inherits = FALSE ) ) {
      return( get( cache_key, envir = expression_cache, inherits = FALSE ) )
    }
    cells <- metadata_7dpi$cell[ metadata_7dpi$sample_id == sample_id & metadata_7dpi$cellchat_group == cellchat_group ]
    result <- entity_expression_summary( entity = entity, cells = cells, expression_data = expression_data, complex_table = complex_table )
    assign( cache_key, result, envir = expression_cache )
    result
  }
  sample_audit_parts <- list()
  audit_index <- 0L
  for ( interaction_index in seq_len( nrow( full_branch_confirmation ) ) ) {
    current <- full_branch_confirmation[ interaction_index, , drop = FALSE ]
    for (sample_id in samples_7dpi) {
      arm1_ligand <- get_cached_entity_summary( current$ligand, sample_id, "Microglia_Arm_1" )
      arm2_ligand <- get_cached_entity_summary( current$ligand, sample_id, "Microglia_Arm_2" )
      target_receptor <- get_cached_entity_summary( current$receptor, sample_id, current$target )
      ligand_delta <- arm2_ligand[["entity_mean_expression"]] - arm1_ligand[["entity_mean_expression"]]
      ligand_direction_match <- if ( current$dominant_branch == "Arm_2" ) {
        is.finite( ligand_delta ) && ligand_delta >
            0
      } else if (
        current$dominant_branch == "Arm_1"
      ) {
        is.finite( ligand_delta ) && ligand_delta <
            0
      } else {
        FALSE
      }
      sender_evaluable <- arm1_ligand[["n_cells"]] >= MIN_SENDER_CELLS_PER_SAMPLE && arm2_ligand[["n_cells"]] >= MIN_SENDER_CELLS_PER_SAMPLE
      target_evaluable <- target_receptor[["n_cells"]] >= MIN_TARGET_CELLS_PER_SAMPLE
      receptor_supported <- target_evaluable && is.finite( target_receptor[["entity_percent_expressing"]] ) && target_receptor[["entity_percent_expressing"]] >= MIN_RECEPTOR_PERCENT_EXPRESSED
      axis_supported <- sender_evaluable && target_evaluable && receptor_supported && ligand_direction_match
      audit_index <- audit_index + 1L
      sample_audit_parts[[audit_index]] <- data.frame( sample_id = sample_id, target = current$target, ligand = current$ligand, receptor = current$receptor, interaction_name = current$interaction_name, pathway_name = current$pathway_name, dominant_branch_full = current$dominant_branch, arm1_ligand_mean = arm1_ligand[["entity_mean_expression"]], arm2_ligand_mean = arm2_ligand[["entity_mean_expression"]], ligand_delta_Arm2_minus_Arm1 = ligand_delta, target_receptor_mean = target_receptor[["entity_mean_expression"]], target_receptor_percent_expressing = target_receptor[["entity_percent_expressing"]], arm1_sender_cells = arm1_ligand[["n_cells"]], arm2_sender_cells = arm2_ligand[["n_cells"]], target_cells = target_receptor[["n_cells"]], sender_evaluable = sender_evaluable, target_evaluable = target_evaluable, ligand_direction_match = ligand_direction_match, receptor_supported = receptor_supported, axis_supported = axis_supported, stringsAsFactors = FALSE )
    }
  }
  sample_audit_long <- if ( length( sample_audit_parts ) > 0L ) {
    do.call( rbind, sample_audit_parts )
  } else {
    data.frame()
  }
  write_tsv( sample_audit_long, SAMPLE_AUDIT_LONG_FILE )
  sample_support_summary <- data.table::as.data.table( sample_audit_long )[ , .( n_samples_sender_evaluable = sum( sender_evaluable ), n_samples_target_evaluable = sum( target_evaluable ), n_samples_ligand_direction_match = sum( ligand_direction_match ), n_samples_receptor_supported = sum( receptor_supported ), n_samples_axis_supported = sum( axis_supported ), n_7dpi_samples = data.table::uniqueN( sample_id ) ), by = .( target, ligand, receptor, interaction_name, pathway_name, dominant_branch_full ) ]
  sample_support_summary <- as.data.frame( sample_support_summary, stringsAsFactors = FALSE )
  sample_support_summary$full_sample_replication_status <- ifelse( sample_support_summary$n_samples_axis_supported == sample_support_summary$n_7dpi_samples, "REPLICATED_IN_ALL_7DPI_SAMPLES", ifelse( sample_support_summary$n_samples_axis_supported > 0, "PARTIAL_SAMPLE_SUPPORT", "NOT_SAMPLE_SUPPORTED" ) )
  preview_for_merge <- preview_lr[ preview_lr$interaction_name %in% selected_interactions, , drop = FALSE ]
  preview_keep <- c( "target", "ligand", "receptor", "interaction_name", "pathway_name", "probability_Arm_1", "p_value_Arm_1", "probability_Arm_2", "p_value_Arm_2", "delta_Arm2_minus_Arm1", "log2_ratio_Arm2_over_Arm1", "dominant_branch", "significant_either", "minimum_fold_change_pass", "Step7G_candidate_ligand", "Step7G_candidate_TFs" )
  preview_keep <- preview_keep[ preview_keep %in% colnames( preview_for_merge ) ]
  preview_for_merge <- preview_for_merge[ , preview_keep, drop = FALSE ]
  preview_nonkey <- setdiff( colnames( preview_for_merge ), c( "target", "ligand", "receptor", "interaction_name", "pathway_name" ) )
  colnames( preview_for_merge )[ match( preview_nonkey, colnames( preview_for_merge ) ) ] <- paste0( preview_nonkey, "_preview" )
  full_for_merge <- full_branch_confirmation
  full_nonkey <- setdiff( colnames( full_for_merge ), c( "target", "ligand", "receptor", "interaction_name", "pathway_name" ) )
  colnames( full_for_merge )[ match( full_nonkey, colnames( full_for_merge ) ) ] <- paste0( full_nonkey, "_full" )
  preview_vs_full <- merge( preview_for_merge, full_for_merge, by = c( "target", "ligand", "receptor", "interaction_name", "pathway_name" ), all = TRUE, sort = FALSE )
  preview_vs_full <- merge( preview_vs_full, sample_support_summary, by.x = c( "target", "ligand", "receptor", "interaction_name", "pathway_name", "dominant_branch_full" ), by.y = c( "target", "ligand", "receptor", "interaction_name", "pathway_name", "dominant_branch_full" ), all.x = TRUE, sort = FALSE )
  preview_vs_full$direction_matches_preview <- !is.na( preview_vs_full$dominant_branch_preview ) & !is.na( preview_vs_full$dominant_branch_full ) & preview_vs_full$dominant_branch_preview == preview_vs_full$dominant_branch_full
  preview_vs_full$dominant_branch_significant_full <- ifelse( preview_vs_full$dominant_branch_full == "Arm_2", preview_vs_full$p_value_Arm_2_full <= CELLCHAT_P_THRESHOLD, ifelse( preview_vs_full$dominant_branch_full == "Arm_1", preview_vs_full$p_value_Arm_1_full <= CELLCHAT_P_THRESHOLD, FALSE ) )
  preview_vs_full$confirmation_status <- ifelse( preview_vs_full$direction_matches_preview & preview_vs_full$dominant_branch_significant_full & preview_vs_full$minimum_fold_change_pass_full & preview_vs_full$full_sample_replication_status == "REPLICATED_IN_ALL_7DPI_SAMPLES", "CONFIRMED_CORE", ifelse( preview_vs_full$direction_matches_preview & preview_vs_full$dominant_branch_significant_full & preview_vs_full$full_sample_replication_status == "REPLICATED_IN_ALL_7DPI_SAMPLES", "CONFIRMED_DIRECTION_ONLY", ifelse( preview_vs_full$direction_matches_preview & preview_vs_full$dominant_branch_significant_full & preview_vs_full$full_sample_replication_status == "PARTIAL_SAMPLE_SUPPORT", "PARTIAL_SAMPLE_CONTEXT", "NOT_CONFIRMED" ) ) )
  preview_vs_full$mechanistic_priority <- ifelse( preview_vs_full$confirmation_tier_full == "PRIMARY" & preview_vs_full$confirmation_status == "CONFIRMED_CORE", "TIER_1_CORE", ifelse( preview_vs_full$confirmation_tier_full == "PRIMARY" & preview_vs_full$confirmation_status %in% c( "CONFIRMED_DIRECTION_ONLY", "PARTIAL_SAMPLE_CONTEXT" ), "TIER_2_PRIMARY_CONTEXT", ifelse( preview_vs_full$confirmation_status == "CONFIRMED_CORE", "TIER_2_CONFIRMED_CONTEXT", "TIER_3_NOT_CORE" ) ) )
  preview_vs_full <- preview_vs_full[ order( match( preview_vs_full$mechanistic_priority, c( "TIER_1_CORE", "TIER_2_PRIMARY_CONTEXT", "TIER_2_CONFIRMED_CONTEXT", "TIER_3_NOT_CORE" ) ), -abs( preview_vs_full$delta_Arm2_minus_Arm1_full ) ), , drop = FALSE ]
  write_tsv( preview_vs_full, PREVIEW_VS_FULL_FILE )
  confirmed_core <- preview_vs_full[ preview_vs_full$confirmation_status == "CONFIRMED_CORE", , drop = FALSE ]
  write_tsv( confirmed_core, CONFIRMED_CORE_FILE )
  pathway_summary <- data.table::as.data.table( preview_vs_full )[ , .( n_axes = .N, n_confirmed_core = sum( confirmation_status == "CONFIRMED_CORE", na.rm = TRUE ), n_confirmed_direction_only = sum( confirmation_status == "CONFIRMED_DIRECTION_ONLY", na.rm = TRUE ), n_partial_sample_context = sum( confirmation_status == "PARTIAL_SAMPLE_CONTEXT", na.rm = TRUE ), Arm1_total_probability_full = sum( probability_Arm_1_full, na.rm = TRUE ), Arm2_total_probability_full = sum( probability_Arm_2_full, na.rm = TRUE ), delta_Arm2_minus_Arm1_full = sum( delta_Arm2_minus_Arm1_full, na.rm = TRUE ) ), by = pathway_name ]
  pathway_summary$dominant_branch_full <- ifelse( pathway_summary$delta_Arm2_minus_Arm1_full > 0, "Arm_2", ifelse( pathway_summary$delta_Arm2_minus_Arm1_full < 0, "Arm_1", "Tie" ) )
  pathway_summary <- pathway_summary[ order( -abs( delta_Arm2_minus_Arm1_full ) ) ]
  pathway_summary <- as.data.frame( pathway_summary, stringsAsFactors = FALSE )
  write_tsv( pathway_summary, PATHWAY_SUMMARY_FILE )
  figure_data <- preview_vs_full[ is.finite( preview_vs_full$delta_Arm2_minus_Arm1_full ), , drop = FALSE ]
  if (nrow(figure_data) > 0L) {
    figure_data$axis_label <- paste0( figure_data$interaction_name, " -> ", figure_data$target )
    figure_data <- utils::head( figure_data[ order( -abs( figure_data$delta_Arm2_minus_Arm1_full ) ), , drop = FALSE ], 40L )
    figure_data$axis_label <- factor( figure_data$axis_label, levels = rev( figure_data$axis_label ) )
    plot_delta <- ggplot2::ggplot( figure_data, ggplot2::aes( x = delta_Arm2_minus_Arm1_full, y = axis_label ) ) + ggplot2::geom_col() + ggplot2::geom_vline( xintercept = 0, linetype = 2 ) + ggplot2::labs( x = "Full-cell probability delta: Arm_2 - Arm_1", y = NULL, title = "7 dpi targeted full-cell CellChat confirmation" ) + ggplot2::theme_bw( base_size = 9 )
    ggplot2::ggsave( filename = FIGURE_DELTA_FILE, plot = plot_delta, width = 9, height = 10, units = "in" )
  }
  scatter_data <- preview_vs_full[ is.finite( preview_vs_full$delta_Arm2_minus_Arm1_preview ) & is.finite( preview_vs_full$delta_Arm2_minus_Arm1_full ), , drop = FALSE ]
  if (nrow(scatter_data) > 0L) {
    plot_preview_full <- ggplot2::ggplot( scatter_data, ggplot2::aes( x = delta_Arm2_minus_Arm1_preview, y = delta_Arm2_minus_Arm1_full ) ) + ggplot2::geom_point() + ggplot2::geom_abline( intercept = 0, slope = 1, linetype = 2 ) + ggplot2::geom_hline( yintercept = 0, linetype = 3 ) + ggplot2::geom_vline( xintercept = 0, linetype = 3 ) + ggplot2::labs( x = "Preview delta: Arm_2 - Arm_1", y = "Full-cell nboot=100 delta: Arm_2 - Arm_1", title = "Preview versus formal confirmation" ) + ggplot2::theme_bw( base_size = 10 )
    ggplot2::ggsave( filename = FIGURE_PREVIEW_FULL_FILE, plot = plot_preview_full, width = 7, height = 6, units = "in" )
  }
  decision_rows <- list()
  decision_rows <- add_decision_row( decision_rows, "Required input files", all( input_audit$exists ), detail = paste0( sum( input_audit$exists ), "/", nrow( input_audit ), " required files present" ) )
  decision_rows <- add_decision_row( decision_rows, "Step8A and Step8B preview PASS gates", TRUE, detail = "Both upstream decision gates passed" )
  decision_rows <- add_decision_row( decision_rows, "Full eligible 7 dpi cell set", nrow( metadata_7dpi ) == 14505L, detail = paste0( nrow( metadata_7dpi ), " cells across ", length( expected_groups ), " groups" ), status_if_fail = "REVIEW_REQUIRED" )
  decision_rows <- add_decision_row( decision_rows, "Frozen branch sender counts", sum( metadata_7dpi$cellchat_group == "Microglia_Arm_1" ) == 1399L && sum( metadata_7dpi$cellchat_group == "Microglia_Arm_2" ) == 1142L, detail = paste0( "Arm_1=", sum( metadata_7dpi$cellchat_group == "Microglia_Arm_1" ), "; Arm_2=", sum( metadata_7dpi$cellchat_group == "Microglia_Arm_2" ) ), status_if_fail = "REVIEW_REQUIRED" )
  decision_rows <- add_decision_row( decision_rows, "Targeted LR catalog", nrow( selected_lr_db ) == length( selected_interactions ) && length( selected_interactions ) >= length( PRIMARY_INTERACTIONS ), detail = paste0( nrow( selected_lr_db ), " interactions; ", sum( selected_lr_db$interaction_name %in% PRIMARY_INTERACTIONS ), " primary" ) )
  decision_rows <- add_decision_row( decision_rows, "RNG-safe nboot=100 CellChat inference", inherits( cellchat, "CellChat" ) && identical( as.integer( cellchat@options$parameter$nboot ), NBOOT_CONFIRMATION ), detail = paste0( "status=", run_status, "; workers=1; nboot=", cellchat@options$parameter$nboot ) )
  decision_rows <- add_decision_row( decision_rows, "Formal branch confirmation table", nrow( preview_vs_full ) > 0L, detail = paste0( nrow( preview_vs_full ), " LR-target rows" ) )
  decision_rows <- add_decision_row( decision_rows, "Two-sample full-cell expression audit", nrow( sample_audit_long ) > 0L && length( unique( sample_audit_long$sample_id ) ) == 2L, detail = paste0( nrow( sample_audit_long ), " sample-axis rows; ", sum( sample_support_summary$full_sample_replication_status == "REPLICATED_IN_ALL_7DPI_SAMPLES" ), " axes replicated in both samples" ) )
  decision_rows <- add_decision_row( decision_rows, "Confirmed core axes produced", nrow( confirmed_core ) > 0L, detail = paste0( nrow( confirmed_core ), " CONFIRMED_CORE LR-target axes" ), status_if_fail = "REVIEW_REQUIRED" )
  overall_pass <- all( vapply( decision_rows, function(row) {
        isTRUE( row$pass[[1L]] )
      },
      FUN.VALUE = logical(1)
    )
  )
  decision_rows <- add_decision_row( decision_rows, "Overall Step 8B targeted confirmation", overall_pass, detail = if ( overall_pass ) {
      "PASS"
    } else {
      "At least one required domain needs review"
    },
    status_if_fail = "REVIEW_REQUIRED"
  )
  final_decision <- do.call( rbind, decision_rows )
  rownames( final_decision ) <- NULL
  write_csv_utf8( final_decision, FINAL_DECISION_FILE )
  result_object <- list( status = if ( overall_pass ) {
      "PASS"
    } else {
      "REVIEW_REQUIRED"
    },
    generated = timestamp(), parameters = CELLCHAT_PARAMETERS, input_cells = nrow( metadata_7dpi ), input_groups = expected_groups, selected_interactions = selected_interactions, targeted_lr_catalog = targeted_lr_catalog, cell_counts = group_counts, cellchat_checkpoint = safe_path( CELLCHAT_CHECKPOINT_FILE ), branch_confirmation = full_branch_confirmation, sample_audit_long = sample_audit_long, preview_vs_full = preview_vs_full, confirmed_core = confirmed_core, pathway_summary = pathway_summary, final_decision = final_decision
  )
  saveRDS( result_object, RESULT_RDS_FILE, compress = FALSE )
  elapsed_seconds <- as.numeric( difftime( Sys.time(), start_time, units = "secs" ) )
  top_core_text <- if ( nrow( confirmed_core ) > 0L ) {
    paste( capture.output( print( utils::head( confirmed_core[ , intersect( c( "target", "ligand", "receptor", "interaction_name", "pathway_name", "dominant_branch_full", "probability_Arm_1_full", "probability_Arm_2_full", "delta_Arm2_minus_Arm1_full", "p_value_Arm_1_full", "p_value_Arm_2_full", "full_sample_replication_status", "mechanistic_priority" ), colnames( confirmed_core ) ), drop = FALSE ], 30L ), row.names = FALSE ) ), collapse = "\n" )
  } else {
    "No CONFIRMED_CORE axes were produced."
  }
  invisible( result_object )
}
tryCatch( main(), error = function(e) {
    failed_call <- tryCatch( paste( deparse( sys.call( -1L ) ), collapse = " " ), error = function(inner_error) {
        NA_character_
      }
    )
    error_lines <- c( "======================================================================", "GSE162610 | STEP 8B TARGETED FULL-CELL CONFIRMATION ERROR", "======================================================================", paste0( "Time: ", timestamp() ), paste0( "Message: ", conditionMessage( e ) ), paste0( "Failed call: ", failed_call ), "", paste0( "If the targeted checkpoint was already saved, rerunning this same ", "file will reuse it." ), "rerun Step1-Step8B preview.", "======================================================================" )
    writeLines( error_lines, con = ERROR_FILE, useBytes = TRUE )
    stop( e )
  }
)
