rm(list = ls())
gc()
options( stringsAsFactors = FALSE, width = 240, future.globals.maxSize = 24 * 1024^3, timeout = 1800 )
set.seed(20260802)
PROJECT_ROOT <- Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
EXPECTED_ATLAS_CELLS <- 62842L
EXPECTED_CLEAN_MICROGLIA_CELLS <- 19958L
EXPECTED_FROZEN_BRANCH_CELLS <- 4203L
EXPECTED_BRANCH_COUNTS <- c( Shared_prebranch = 116L, Arm_1 = 2289L, Arm_2 = 1798L )
STEP1F_CLUSTER_TO_BROAD_TYPE <- c( "0" = "Microglia", "1" = "Endothelial", "2" = "Macrophage_Monocyte", "3" = "Ependymal", "4" = "Myeloid_general", "5" = "OPC", "6" = "Neutrophil", "7" = "Pericyte_VSMC", "8" = "Myeloid_general", "9" = "T_cell", "10" = "Oligodendrocyte", "11" = "Fibroblast", "12" = "Astrocyte", "13" = "Neuron", "14" = "Ependymal", "15" = "OPC", "16" = "Neutrophil", "17" = "Microglia", "18" = "Endothelial" )
EXPECTED_STEP1F_CLUSTERS <- names( STEP1F_CLUSTER_TO_BROAD_TYPE )
FORMAL_TIMEPOINTS <- c( "1dpi", "3dpi", "7dpi" )
MIN_CELLS_PER_GROUP <- 20L
MIN_SAMPLES_PER_TIMEPOINT <- 2L
STEP8B_PARAMETERS <- list( database = "CellChatDB.mouse", database_subset = "all_categories", probability_type = "triMean", trim = 0.10, raw.use = TRUE, population.size = FALSE, distance.use = FALSE, min.cells = MIN_CELLS_PER_GROUP, nboot = 100L, seed.use = 20260802L )
ATLAS_RDS_FILE <- file.path( PROJECT_ROOT, "04_results", "step1F_GSE162610_broad_annotation", "GSE162610_broad_annotation_audit.rds" )
CLEAN_MICROGLIA_RDS_FILE <- file.path( PROJECT_ROOT, "04_results", "step3D_GSE162610_time_DEG_pathway_confirmation", "GSE162610_step3D_time_DEG_pathway_confirmed.rds" )
FROZEN_BRANCH_META_FILE <- file.path( PROJECT_ROOT, "04_results", "step7A_GSE162610_pySCENIC_input_environment_audit", "input", "GSE162610_step7A_frozen_4203_cell_metadata.csv" )
STEP7G_DIR <- file.path( PROJECT_ROOT, "04_results", "step7G_GSE162610_regulon_target_program_closure" )
STEP7G_DECISION_FILE <- file.path( STEP7G_DIR, "tables", "08_final_decision.csv" )
STEP7G_SIGNATURE_TABLE_FILE <- file.path( STEP7G_DIR, "tables", "09_candidate_target_signature_summary.tsv" )
STEP7G_CLOSURE_FILE <- file.path( STEP7G_DIR, "tables", "03_timing_candidate_target_closure.tsv" )
OUTDIR <- file.path( PROJECT_ROOT, "04_results", "step8A_GSE162610_CellChat_input_preflight" )
TABDIR <- file.path(OUTDIR, "tables")
FIGDIR <- file.path(OUTDIR, "figures")
RDSDIR <- file.path(OUTDIR, "RDS")
REPORTDIR <- file.path(OUTDIR, "reports")
for (directory in c( OUTDIR, TABDIR, FIGDIR, RDSDIR, REPORTDIR )) {
  dir.create( directory, recursive = TRUE, showWarnings = FALSE )
}
INPUT_AUDIT_FILE <- file.path( TABDIR, "00_input_file_audit.csv" )
LABEL_MAPPING_FILE <- file.path( TABDIR, "01_cellchat_group_mapping.csv" )
CLUSTER_RECONSTRUCTION_FILE <- file.path( TABDIR, "01b_Step1F_cluster_annotation_reconstruction.csv" )
CELL_METADATA_FILE <- file.path( TABDIR, "02_CellChat_cell_metadata.tsv.gz" )
GROUP_TIME_COUNT_FILE <- file.path( TABDIR, "03_cellchat_group_counts_by_timepoint.csv" )
GROUP_SAMPLE_COUNT_FILE <- file.path( TABDIR, "04_cellchat_group_counts_by_sample.csv" )
ELIGIBILITY_FILE <- file.path( TABDIR, "05_timepoint_group_eligibility.csv" )
CANDIDATE_LIGAND_FILE <- file.path( TABDIR, "06_candidate_regulon_ligand_catalog.tsv" )
DATABASE_AUDIT_FILE <- file.path( TABDIR, "07_CellChatDB_mouse_audit.csv" )
FINAL_DECISION_FILE <- file.path( TABDIR, "08_final_decision.csv" )
MANIFEST_RDS_FILE <- file.path( RDSDIR, "GSE162610_step8A_CellChat_input_manifest.rds" )
REPORT_FILE <- file.path( REPORTDIR, "Step8A_report.txt" )
SESSION_FILE <- file.path( REPORTDIR, "Step8A_sessionInfo.txt" )
ERROR_FILE <- file.path( REPORTDIR, "Step8A_ERROR.txt" )
FIG_GROUP_COMPOSITION <- file.path( FIGDIR, "Fig8A_01_CellChat_group_composition_by_timepoint.pdf" )
FIG_BRANCH_TIME <- file.path( FIGDIR, "Fig8A_02_frozen_branch_cells_by_timepoint.pdf" )
FIG_CANDIDATE_LIGANDS <- file.path( FIGDIR, "Fig8A_03_candidate_regulon_ligand_expression.pdf" )
REQUIRED_PACKAGES <- c( "Seurat", "SeuratObject", "Matrix", "data.table", "ggplot2", "CellChat" )
MISSING_PACKAGES <- REQUIRED_PACKAGES[ !vapply( REQUIRED_PACKAGES, requireNamespace, logical(1), quietly = TRUE ) ]
if (length(MISSING_PACKAGES) > 0L) {
  stop( paste0( "Missing required R package(s): ", paste(MISSING_PACKAGES, collapse = ", "), "\nInstall the missing package(s) and rerun this script.", "\nThis script does not automatically update CellChat because an existing ", "working CellChat environment should be preserved." ) )
}
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(data.table)
  library(ggplot2)
  library(CellChat)
})
timestamp <- function() {
  format( Sys.time(), "%Y-%m-%d %H:%M:%S" )
}
safe_path <- function(x) {
  normalizePath( x, winslash = "/", mustWork = FALSE )
}
write_csv_utf8 <- function(x, file) {
  data.table::fwrite( x, file = file, sep = ",", quote = TRUE, na = "NA" )
}
write_tsv <- function(x, file) {
  data.table::fwrite( x, file = file, sep = "\t", quote = FALSE, na = "NA" )
}
write_tsv_gz <- function(x, file) {
  data.table::fwrite( x, file = file, sep = "\t", quote = FALSE, na = "NA", compress = "gzip" )
}
safe_numeric <- function(x) {
  suppressWarnings( as.numeric( as.character(x) ) )
}
collapse_unique <- function(x) {
  x <- unique( as.character(x) )
  x <- x[ !is.na(x) & nzchar(x) ]
  if (length(x) == 0L) {
    return(NA_character_)
  }
  paste( sort(x), collapse = " | " )
}
standardize_timepoint <- function(x) {
  original <- as.character(x)
  clean <- tolower( trimws(original) )
  clean <- gsub( "[[:space:]_.-]+", "", clean )
  out <- rep( NA_character_, length(clean) )
  out[ grepl( "uninj|uninjured|naive|control|sham", clean ) ] <- "Uninjured"
  out[ grepl( "(^1dpi$)|(^dpi1$)|(^day1$)|(^1day$)|(^d1$)", clean ) ] <- "1dpi"
  out[ grepl( "(^3dpi$)|(^dpi3$)|(^day3$)|(^3day$)|(^d3$)", clean ) ] <- "3dpi"
  out[ grepl( "(^7dpi$)|(^dpi7$)|(^day7$)|(^7day$)|(^d7$)", clean ) ] <- "7dpi"
  out
}
detect_column <- function( df, candidates, required = TRUE, description = "column" ) {
  hit <- candidates[ candidates %in% colnames(df) ]
  if (length(hit) > 0L) {
    return(hit[[1L]])
  }
  if (required) {
    stop( paste0( "Could not detect ", description, ". Tried: ", paste( candidates, collapse = ", " ) ) )
  }
  NA_character_
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
    stop( paste0( "Unable to retrieve the normalized RNA data layer from the atlas object. ", "Step8 requires non-negative log-normalized RNA expression." ) )
  }
  output
}
decision_pass <- function( decision_df, domain_pattern ) {
  if ( !"domain" %in% colnames(decision_df) || !"pass" %in% colnames(decision_df) ) {
    return(FALSE)
  }
  index <- grepl( domain_pattern, as.character( decision_df$domain ), ignore.case = TRUE )
  if (!any(index)) {
    return(FALSE)
  }
  pass_values <- tolower( as.character( decision_df$pass[index] ) )
  any( pass_values %in% c( "true", "t", "1", "pass" ) )
}
extract_cell_column <- function(df) {
  candidates <- c( "cell", "barcode", "cell_id", "Cell", "Barcode" )
  detected <- candidates[ candidates %in% colnames(df) ]
  if (length(detected) > 0L) {
    return(detected[[1L]])
  }
  first_column <- colnames(df)[[1L]]
  first_values <- as.character( df[[first_column]] )
  if ( length(unique(first_values)) == nrow(df) ) {
    return(first_column)
  }
  stop( "Could not determine the cell-barcode column in the frozen Step7A metadata." )
}
split_pipe_genes <- function(x) {
  x <- as.character(x)
  if ( is.na(x) || !nzchar(x) ) {
    return(character(0))
  }
  genes <- unlist( strsplit( x, "\\s*\\|\\s*" ), use.names = FALSE )
  genes <- trimws( genes )
  unique( genes[ !is.na(genes) & nzchar(genes) ] )
}
load_cellchat_mouse_db <- function() {
  database <- tryCatch( get( "CellChatDB.mouse", envir = asNamespace("CellChat"), inherits = FALSE ), error = function(e) {
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
  stop( "CellChatDB.mouse could not be loaded from the installed CellChat package." )
}
expand_database_entities <- function( entities, complex_table ) {
  entities <- unique( as.character(entities) )
  entities <- entities[ !is.na(entities) & nzchar(entities) ]
  output <- vector( "list", length(entities) )
  for (entity_index in seq_along(entities)) {
    entity <- entities[[entity_index]]
    components <- entity
    if ( !is.null(complex_table) && nrow(complex_table) > 0L && entity %in% rownames(complex_table) ) {
      components <- as.character( unlist( complex_table[ entity, , drop = FALSE ], use.names = FALSE ) )
      components <- components[ !is.na(components) & nzchar(components) ]
      if (length(components) == 0L) {
        components <- entity
      }
    }
    output[[entity_index]] <- data.frame( database_entity = entity, gene = unique(components), stringsAsFactors = FALSE )
  }
  if (length(output) == 0L) {
    return( data.frame( database_entity = character(0), gene = character(0), stringsAsFactors = FALSE ) )
  }
  unique( do.call( rbind, output ) )
}
summarize_gene_expression <- function( expression_matrix, genes, cells ) {
  genes <- intersect( unique( as.character(genes) ), rownames(expression_matrix) )
  cells <- intersect( unique( as.character(cells) ), colnames(expression_matrix) )
  if ( length(genes) == 0L ) {
    return( data.frame( gene = character(0), mean_expression = numeric(0), percent_expressing = numeric(0), n_cells = integer(0), stringsAsFactors = FALSE ) )
  }
  if ( length(cells) == 0L ) {
    return( data.frame( gene = genes, mean_expression = rep( NA_real_, length(genes) ), percent_expressing = rep( NA_real_, length(genes) ), n_cells = rep( 0L, length(genes) ), stringsAsFactors = FALSE ) )
  }
  current_matrix <- expression_matrix[ genes, cells, drop = FALSE ]
  output <- data.frame( gene = genes, mean_expression = as.numeric( Matrix::rowMeans( current_matrix ) ), percent_expressing = as.numeric( Matrix::rowSums( current_matrix > 0 ) / length(cells) * 100 ), n_cells = length(cells), stringsAsFactors = FALSE )
  output
}
add_decision_row <- function( decision_list, domain, pass, status_if_fail = "REVIEW_REQUIRED", detail = "" ) {
  decision_list[[ length(decision_list) + 1L ]] <- data.frame( domain = domain, pass = isTRUE(pass), status = if ( isTRUE(pass) ) {
      "PASS"
    } else {
      status_if_fail
    },
    detail = as.character(detail), stringsAsFactors = FALSE
  )
  decision_list
}
run_step8A <- function() {
  overall_start <- Sys.time()
  required_files <- c( ATLAS_RDS_FILE, CLEAN_MICROGLIA_RDS_FILE, FROZEN_BRANCH_META_FILE, STEP7G_DECISION_FILE, STEP7G_SIGNATURE_TABLE_FILE, STEP7G_CLOSURE_FILE )
  missing_files <- required_files[ !file.exists( required_files ) ]
  input_audit <- data.frame( role = c( "Full GSE162610 atlas with broad annotation", "Validated clean microglia object", "Frozen Step7A 4203-cell metadata", "Step7G final decision", "Step7G candidate signature table", "Step7G timing-candidate closure table" ), path = safe_path( required_files ), exists = file.exists( required_files ), size_bytes = ifelse( file.exists(required_files), file.info(required_files)$size, NA_real_ ), stringsAsFactors = FALSE )
  write_csv_utf8( input_audit, INPUT_AUDIT_FILE )
  if (length(missing_files) > 0L) {
    stop( paste0( "Required Step8A input file(s) are missing:\n", paste( safe_path(missing_files), collapse = "\n" ) ) )
  }
  step7g_decision <- data.table::fread( STEP7G_DECISION_FILE, data.table = FALSE, showProgress = FALSE )
  step7g_pass <- decision_pass( step7g_decision, "Overall Step 7G" )
  if (!step7g_pass) {
    stop( "Step7G does not contain an Overall PASS decision." )
  }
  atlas <- readRDS( ATLAS_RDS_FILE )
  if (!inherits(atlas, "Seurat")) {
    stop( "The full atlas input is not a Seurat object." )
  }
  atlas_cells <- colnames( atlas )
  atlas_meta <- atlas@meta.data
  atlas_meta$cell <- rownames( atlas_meta )
  if ( !identical( atlas_meta$cell, atlas_cells ) ) {
    atlas_meta <- atlas_meta[ atlas_cells, , drop = FALSE ]
    atlas_meta$cell <- rownames( atlas_meta )
  }
  required_atlas_meta <- c( "sample_id", "timepoint", "cluster_provisional" )
  missing_atlas_meta <- setdiff( required_atlas_meta, colnames(atlas_meta) )
  if (length(missing_atlas_meta) > 0L) {
    stop( paste0( "The atlas object is missing required metadata: ", paste( missing_atlas_meta, collapse = ", " ) ) )
  }
  clean_microglia_object <- readRDS( CLEAN_MICROGLIA_RDS_FILE )
  if (!inherits(clean_microglia_object, "Seurat")) {
    stop( "The clean-microglia input is not a Seurat object." )
  }
  clean_microglia_cells <- colnames( clean_microglia_object )
  rm(clean_microglia_object)
  gc()
  clean_microglia_in_atlas <- intersect( clean_microglia_cells, atlas_cells )
  clean_microglia_missing_from_atlas <- setdiff( clean_microglia_cells, atlas_cells )
  frozen_meta <- data.table::fread( FROZEN_BRANCH_META_FILE, data.table = FALSE, showProgress = FALSE )
  frozen_cell_column <- extract_cell_column( frozen_meta )
  frozen_state_column <- detect_column( frozen_meta, candidates = c( "step5CR_state", "BranchState", "branch_state", "state" ), required = TRUE, description = "frozen branch-state column" )
  frozen_time_column <- detect_column( frozen_meta, candidates = c( "timepoint", "time", "Timepoint" ), required = FALSE, description = "frozen timepoint column" )
  frozen_cells <- as.character( frozen_meta[[frozen_cell_column]] )
  frozen_states <- as.character( frozen_meta[[frozen_state_column]] )
  names(frozen_states) <- frozen_cells
  if (anyDuplicated(frozen_cells) > 0L) {
    stop( "Duplicated cell barcodes were detected in the frozen Step7A metadata." )
  }
  allowed_frozen_states <- names( EXPECTED_BRANCH_COUNTS )
  unexpected_states <- setdiff( unique(frozen_states), allowed_frozen_states )
  if (length(unexpected_states) > 0L) {
    stop( paste0( "Unexpected frozen branch state(s): ", paste( unexpected_states, collapse = ", " ) ) )
  }
  frozen_count_table <- table( factor( frozen_states, levels = allowed_frozen_states ) )
  frozen_count_vector <- as.integer( frozen_count_table )
  names(frozen_count_vector) <- allowed_frozen_states
  frozen_count_exact <- all( unname( frozen_count_vector ) == unname( as.integer( EXPECTED_BRANCH_COUNTS ) ) )
  frozen_cells_in_atlas <- intersect( frozen_cells, atlas_cells )
  frozen_cells_missing_from_atlas <- setdiff( frozen_cells, atlas_cells )
  frozen_cells_in_clean_microglia <- intersect( frozen_cells, clean_microglia_cells )
  frozen_not_in_clean_microglia <- setdiff( frozen_cells, clean_microglia_cells )
  atlas_meta$sample_id <- as.character( atlas_meta$sample_id )
  atlas_meta$timepoint_raw <- as.character( atlas_meta$timepoint )
  atlas_meta$timepoint_standard <- standardize_timepoint( atlas_meta$timepoint_raw )
  atlas_meta$cluster_provisional <- as.character( atlas_meta$cluster_provisional )
  atlas_meta$cluster_key <- tolower( trimws( atlas_meta$cluster_provisional ) )
  atlas_meta$cluster_key <- sub( "^g", "", atlas_meta$cluster_key )
  atlas_meta$cluster_key <- sub( "\\.0+$", "", atlas_meta$cluster_key )
  observed_cluster_keys <- sort( unique( atlas_meta$cluster_key ) )
  unexpected_cluster_keys <- setdiff( observed_cluster_keys, EXPECTED_STEP1F_CLUSTERS )
  missing_expected_cluster_keys <- setdiff( EXPECTED_STEP1F_CLUSTERS, observed_cluster_keys )
  if (length(unexpected_cluster_keys) > 0L) {
    stop( paste0( "Unexpected Step1F provisional cluster ID(s): ", paste( unexpected_cluster_keys, collapse = ", " ) ) )
  }
  if (length(missing_expected_cluster_keys) > 0L) {
    stop( paste0( "Expected Step1F provisional cluster ID(s) are absent: ", paste( missing_expected_cluster_keys, collapse = ", " ) ) )
  }
  if ( "broad_type_suggestion" %in% colnames(atlas_meta) ) {
    atlas_meta$broad_type_suggestion_original <- as.character( atlas_meta$broad_type_suggestion )
  } else {
    atlas_meta$broad_type_suggestion_original <- NA_character_
  }
  atlas_meta$broad_type_reconstructed <- unname( STEP1F_CLUSTER_TO_BROAD_TYPE[ atlas_meta$cluster_key ] )
  if ( any( is.na( atlas_meta$broad_type_reconstructed ) ) ) {
    stop( "The verified Step1F cluster map left one or more atlas cells unresolved." )
  }
  atlas_meta$is_clean_microglia <- atlas_meta$cell %in%
    clean_microglia_cells
  atlas_meta$frozen_branch_state <- unname( frozen_states[ atlas_meta$cell ] )
  broad_group_map <- c( Endothelial = "Endothelial", Oligodendrocyte = "Oligodendrocyte", Fibroblast = "Fibroblast", Astrocyte = "Astrocyte", Neuron = "Neuron", Ependymal = "Ependymal", OPC = "OPC", Neutrophil = "Neutrophil", Pericyte_VSMC = "Pericyte_VSMC", T_cell = "T_cell" )
  atlas_meta$cellchat_group <- unname( broad_group_map[ atlas_meta$broad_type_reconstructed ] )
  atlas_meta$cellchat_group[ atlas_meta$broad_type_reconstructed %in% c( "Microglia", "Myeloid_general", "Macrophage_Monocyte" ) & !atlas_meta$is_clean_microglia ] <- "Other_myeloid"
  atlas_meta$cellchat_group[ atlas_meta$is_clean_microglia ] <- "Microglia_other"
  branch_to_cellchat <- c( Shared_prebranch = "Microglia_Shared_prebranch", Arm_1 = "Microglia_Arm_1", Arm_2 = "Microglia_Arm_2" )
  branch_index <- !is.na( atlas_meta$frozen_branch_state )
  atlas_meta$cellchat_group[ branch_index ] <- unname( branch_to_cellchat[ atlas_meta$frozen_branch_state[ branch_index ] ] )
  atlas_meta$cellchat_group[ is.na( atlas_meta$cellchat_group ) | !nzchar( atlas_meta$cellchat_group ) ] <- "Excluded_unknown"
  atlas_meta$include_step8_timecourse <- atlas_meta$timepoint_standard %in%
      FORMAL_TIMEPOINTS & atlas_meta$cellchat_group != "Excluded_unknown"
  atlas_meta$include_step8_formal_branch <- atlas_meta$timepoint_standard == "7dpi" & atlas_meta$cellchat_group %in%
        c( "Microglia_Arm_1", "Microglia_Arm_2" )
  cluster_reconstruction_rows <- lapply( EXPECTED_STEP1F_CLUSTERS, function(current_cluster) {
      cluster_index <- atlas_meta$cluster_key == current_cluster
      original_values <- atlas_meta$broad_type_suggestion_original[ cluster_index ]
      data.frame( cluster_provisional = current_cluster, cluster_display = paste0( "g", current_cluster ), reconstructed_broad_type = unname( STEP1F_CLUSTER_TO_BROAD_TYPE[ current_cluster ] ), n_cells = sum( cluster_index ), original_broad_type_values = collapse_unique( original_values ), original_nonmissing_fraction = mean( !is.na( original_values ) & nzchar( original_values ) ), stringsAsFactors = FALSE )
    }
  )
  cluster_reconstruction_audit <- do.call( rbind, cluster_reconstruction_rows )
  cluster_reconstruction_audit$reconstruction_complete <- cluster_reconstruction_audit$n_cells > 0L & !is.na( cluster_reconstruction_audit$reconstructed_broad_type )
  write_csv_utf8( cluster_reconstruction_audit, CLUSTER_RECONSTRUCTION_FILE )
  label_mapping <- data.frame( source_rule = c( "Validated clean microglia not in frozen trajectory", "Frozen Step7 Shared_prebranch", "Frozen Step7 Arm_1", "Frozen Step7 Arm_2", names(broad_group_map), "Non-clean broad Microglia/Myeloid_general/Macrophage_Monocyte", "Unresolved broad annotation" ), cellchat_group = c( "Microglia_other", "Microglia_Shared_prebranch", "Microglia_Arm_1", "Microglia_Arm_2", unname(broad_group_map), "Other_myeloid", "Excluded_unknown" ), interpretation = c( "Resident microglia outside the frozen 4203-cell trajectory; neutral background group", "Frozen common prebranch microglia", "Frozen neutral branch label; no repair claim", "Frozen neutral branch label; no scar claim", rep( "Broad recipient-cell class reconstructed from verified Step1F cluster mapping", length(broad_group_map) ), "Myeloid cells not contained in the validated resident-microglia set", "Excluded from formal CellChat input" ), stringsAsFactors = FALSE )
  write_csv_utf8( label_mapping, LABEL_MAPPING_FILE )
  analysis_meta <- atlas_meta[ atlas_meta$include_step8_timecourse, , drop = FALSE ]
  group_time_counts <- as.data.frame( table( analysis_meta$timepoint_standard, analysis_meta$cellchat_group ), stringsAsFactors = FALSE )
  colnames(group_time_counts) <- c( "timepoint", "cellchat_group", "n_cells" )
  group_time_counts <- group_time_counts[ group_time_counts$n_cells > 0L, , drop = FALSE ]
  group_time_counts$eligible_min_cells <- group_time_counts$n_cells >= MIN_CELLS_PER_GROUP
  formal_branch_sender <- group_time_counts$cellchat_group %in%
    c( "Microglia_Arm_1", "Microglia_Arm_2" )
  group_time_counts$analysis_role <- ifelse( group_time_counts$timepoint == "7dpi" & formal_branch_sender & group_time_counts$eligible_min_cells, "Formal_7dpi_branch_sender", ifelse( group_time_counts$eligible_min_cells, "Eligible_timecourse_group", "Below_min_cells_descriptive_only" ) )
  group_time_counts <- group_time_counts[ order( match( group_time_counts$timepoint, FORMAL_TIMEPOINTS ), group_time_counts$cellchat_group ), , drop = FALSE ]
  write_csv_utf8( group_time_counts, GROUP_TIME_COUNT_FILE )
  group_sample_counts <- as.data.frame( table( analysis_meta$timepoint_standard, analysis_meta$sample_id, analysis_meta$cellchat_group ), stringsAsFactors = FALSE )
  colnames(group_sample_counts) <- c( "timepoint", "sample_id", "cellchat_group", "n_cells" )
  group_sample_counts <- group_sample_counts[ group_sample_counts$n_cells > 0L, , drop = FALSE ]
  group_sample_counts <- group_sample_counts[ order( match( group_sample_counts$timepoint, FORMAL_TIMEPOINTS ), group_sample_counts$sample_id, group_sample_counts$cellchat_group ), , drop = FALSE ]
  write_csv_utf8( group_sample_counts, GROUP_SAMPLE_COUNT_FILE )
  sample_time_table <- unique( analysis_meta[ , c( "sample_id", "timepoint_standard" ), drop = FALSE ] )
  samples_per_timepoint <- as.data.frame( table( factor( sample_time_table$timepoint_standard, levels = FORMAL_TIMEPOINTS ) ), stringsAsFactors = FALSE )
  colnames(samples_per_timepoint) <- c( "timepoint", "n_samples" )
  eligibility_rows <- list()
  for (current_timepoint in FORMAL_TIMEPOINTS) {
    current_counts <- group_time_counts[ group_time_counts$timepoint == current_timepoint, , drop = FALSE ]
    eligible_groups <- current_counts$cellchat_group[ current_counts$eligible_min_cells ]
    eligibility_rows[[ length(eligibility_rows) + 1L ]] <- data.frame( timepoint = current_timepoint, n_total_cells = sum( current_counts$n_cells ), n_eligible_groups = length( eligible_groups ), eligible_groups = paste( eligible_groups, collapse = " | " ), n_samples = samples_per_timepoint$n_samples[ samples_per_timepoint$timepoint == current_timepoint ], stringsAsFactors = FALSE )
  }
  eligibility_table <- do.call( rbind, eligibility_rows )
  write_csv_utf8( eligibility_table, ELIGIBILITY_FILE )
  metadata_output <- atlas_meta[ , unique( c( "cell", "sample_id", "timepoint_raw", "timepoint_standard", "cluster_provisional", "cluster_key", "broad_type_suggestion_original", "broad_type_reconstructed", "is_clean_microglia", "frozen_branch_state", "cellchat_group", "include_step8_timecourse", "include_step8_formal_branch", intersect( c( "original_V30", "original_V30_step1F", "final_lineage_label", "microglia_state" ), colnames(atlas_meta) ) ) ), drop = FALSE ]
  write_tsv_gz( metadata_output, CELL_METADATA_FILE )
  DefaultAssay( atlas ) <- "RNA"
  expression_data <- get_rna_data( atlas )
  expression_dimensions <- dim( expression_data )
  expression_cell_order_exact <- identical( colnames(expression_data), atlas_cells )
  if (!expression_cell_order_exact) {
    missing_expression_cells <- setdiff( atlas_cells, colnames(expression_data) )
    if (length(missing_expression_cells) > 0L) {
      stop( paste0( "Normalized RNA data are missing ", length(missing_expression_cells), " atlas cells." ) )
    }
    expression_data <- expression_data[ , atlas_cells, drop = FALSE ]
    expression_cell_order_exact <- identical( colnames(expression_data), atlas_cells )
  }
  audit_cells <- head( analysis_meta$cell, min( 500L, nrow(analysis_meta) ) )
  audit_matrix <- expression_data[ , audit_cells, drop = FALSE ]
  audit_values <- if ( inherits( audit_matrix, "Matrix" ) ) {
    audit_matrix@x
  } else {
    as.numeric( audit_matrix )
  }
  expression_values_finite <- all( is.finite( audit_values ) )
  expression_nonnegative <- if ( length(audit_values) == 0L ) {
    FALSE
  } else {
    min(audit_values) >= 0
  }
  canonical_mouse_markers <- c( "P2ry12", "Tmem119", "Apoe", "Ccl3", "Col1a1", "Aldh1l1", "Pdgfra", "Pecam1" )
  canonical_mouse_markers_present <- intersect( canonical_mouse_markers, rownames(expression_data) )
  cellchat_version <- as.character( packageVersion( "CellChat" ) )
  required_cellchat_functions <- c( "createCellChat", "subsetData", "identifyOverExpressedGenes", "identifyOverExpressedInteractions", "computeCommunProb", "filterCommunication", "computeCommunProbPathway", "aggregateNet", "subsetCommunication", "mergeCellChat", "netAnalysis_computeCentrality" )
  function_available <- vapply( required_cellchat_functions, function(function_name) {
      exists( function_name, envir = asNamespace("CellChat"), inherits = FALSE )
    },
    logical(1)
  )
  cellchat_api_pass <- all( function_available )
  CellChatDB_mouse <- load_cellchat_mouse_db()
  database_components <- names( CellChatDB_mouse )
  database_structure_pass <- all( c( "interaction", "complex" ) %in% database_components )
  interaction_table <- as.data.frame( CellChatDB_mouse$interaction, stringsAsFactors = FALSE )
  if ( !all( c( "ligand", "receptor" ) %in% colnames(interaction_table) ) ) {
    stop( "CellChatDB.mouse interaction table lacks ligand/receptor columns." )
  }
  complex_table <- as.data.frame( CellChatDB_mouse$complex, stringsAsFactors = FALSE )
  ligand_component_map <- expand_database_entities( interaction_table$ligand, complex_table )
  receptor_component_map <- expand_database_entities( interaction_table$receptor, complex_table )
  db_ligand_genes <- unique( ligand_component_map$gene )
  db_receptor_genes <- unique( receptor_component_map$gene )
  db_signaling_genes <- unique( c( db_ligand_genes, db_receptor_genes ) )
  db_genes_present <- intersect( db_signaling_genes, rownames(expression_data) )
  db_overlap_fraction <- length( db_genes_present ) / max( 1L, length( db_signaling_genes ) )
  annotation_column <- if ( "annotation" %in% colnames(interaction_table) ) {
    "annotation"
  } else {
    NA_character_
  }
  annotation_summary <- if ( !is.na(annotation_column) ) {
    as.data.frame( table( interaction_table[[annotation_column]], useNA = "ifany" ), stringsAsFactors = FALSE )
  } else {
    data.frame( Var1 = "annotation_column_unavailable", Freq = nrow( interaction_table ), stringsAsFactors = FALSE )
  }
  colnames(annotation_summary) <- c( "database_category", "n_interactions" )
  database_audit <- rbind( data.frame( metric = c( "CellChat_version", "database_name", "interaction_rows", "complex_rows", "unique_ligand_entities", "unique_receptor_entities", "expanded_ligand_component_genes", "expanded_receptor_component_genes", "database_signaling_genes", "database_signaling_genes_in_atlas", "database_gene_overlap_fraction", "normalized_RNA_genes", "normalized_RNA_cells" ), value = as.character( c( cellchat_version, "CellChatDB.mouse", nrow(interaction_table), nrow(complex_table), length(unique(interaction_table$ligand)), length(unique(interaction_table$receptor)), length(db_ligand_genes), length(db_receptor_genes), length(db_signaling_genes), length(db_genes_present), round(db_overlap_fraction, 6), nrow(expression_data), ncol(expression_data) ) ), stringsAsFactors = FALSE ), data.frame( metric = paste0( "category__", annotation_summary$database_category ), value = as.character( annotation_summary$n_interactions ), stringsAsFactors = FALSE ) )
  write_csv_utf8( database_audit, DATABASE_AUDIT_FILE )
  preferred_smoke_groups <- c( "Microglia_Arm_1", "Microglia_Arm_2", "Astrocyte", "Fibroblast", "Pericyte_VSMC", "Endothelial", "OPC", "Oligodendrocyte", "Neuron", "Other_myeloid" )
  seven_dpi_meta <- analysis_meta[ analysis_meta$timepoint_standard == "7dpi", , drop = FALSE ]
  seven_dpi_group_counts <- sort( table( seven_dpi_meta$cellchat_group ), decreasing = TRUE )
  populated_seven_dpi_groups <- names( seven_dpi_group_counts[ seven_dpi_group_counts >= 10L ] )
  preferred_available <- preferred_smoke_groups[ preferred_smoke_groups %in% populated_seven_dpi_groups ]
  supplemental_available <- setdiff( populated_seven_dpi_groups, preferred_available )
  available_smoke_groups <- unique( c( preferred_available, supplemental_available ) )
  smoke_source_meta <- seven_dpi_meta
  smoke_source_label <- "7dpi"
  if (length(available_smoke_groups) < 2L) {
    all_time_group_counts <- sort( table( analysis_meta$cellchat_group ), decreasing = TRUE )
    populated_all_time_groups <- names( all_time_group_counts[ all_time_group_counts >= 10L ] )
    preferred_available <- preferred_smoke_groups[ preferred_smoke_groups %in% populated_all_time_groups ]
    supplemental_available <- setdiff( populated_all_time_groups, preferred_available )
    available_smoke_groups <- unique( c( preferred_available, supplemental_available ) )
    smoke_source_meta <- analysis_meta
    smoke_source_label <- "combined_1_3_7dpi"
  }
  available_smoke_groups <- head( available_smoke_groups, 4L )
  if (length(available_smoke_groups) < 2L) {
    stop( paste0( "Fewer than two sufficiently populated CellChat groups are available ", "even after the 1/3/7 dpi fallback. Populated 7 dpi groups: ", paste( populated_seven_dpi_groups, collapse = " | " ), "." ) )
  }
  smoke_group_count_table <- data.frame( smoke_source = smoke_source_label, cellchat_group = available_smoke_groups, n_available_cells = vapply( available_smoke_groups, function(current_group) {
        sum( smoke_source_meta$cellchat_group == current_group )
      },
      integer(1)
    ), stringsAsFactors = FALSE
  )
  write_csv_utf8( smoke_group_count_table, file.path( TABDIR, "07b_CellChat_smoke_test_groups.csv" ) )
  smoke_cells <- unlist( lapply( available_smoke_groups, function(current_group) {
        current_cells <- smoke_source_meta$cell[ smoke_source_meta$cellchat_group == current_group ]
        head( current_cells, min( 25L, length(current_cells) ) )
      }
    ), use.names = FALSE
  )
  smoke_meta <- smoke_source_meta[ match( smoke_cells, smoke_source_meta$cell ), c( "cell", "cellchat_group" ), drop = FALSE ]
  rownames(smoke_meta) <- smoke_meta$cell
  smoke_expression <- expression_data[ , smoke_cells, drop = FALSE ]
  smoke_test_message <- NA_character_
  smoke_test_pass <- tryCatch( {
      smoke_cellchat <- CellChat::createCellChat( object = smoke_expression, meta = smoke_meta, group.by = "cellchat_group" )
      smoke_cellchat@DB <- CellChatDB_mouse
      smoke_cellchat <- CellChat::subsetData( smoke_cellchat )
      signaling_rows <- if ( "data.signaling" %in% slotNames(smoke_cellchat) ) {
        nrow( smoke_cellchat@data.signaling )
      } else {
        0L
      }
      smoke_test_message <- paste0( "PASS; source=", smoke_source_label, "; groups=", paste( available_smoke_groups, collapse = "|" ), "; cells=", length(smoke_cells), "; signaling_genes=", signaling_rows )
      signaling_rows > 0L
    },
    error = function(e) {
      smoke_test_message <<- paste0( "FAIL: ", conditionMessage(e) )
      FALSE
    }
  )
  if ( exists( "smoke_cellchat", inherits = FALSE ) ) {
    rm(smoke_cellchat)
  }
  rm( smoke_expression, smoke_meta )
  gc()
  signature_summary <- data.table::fread( STEP7G_SIGNATURE_TABLE_FILE, data.table = FALSE, showProgress = FALSE )
  closure_table <- data.table::fread( STEP7G_CLOSURE_FILE, data.table = FALSE, showProgress = FALSE )
  required_signature_columns <- c( "signature_name", "TF", "regulon", "signature_type", "mechanistic_priority", "genes" )
  missing_signature_columns <- setdiff( required_signature_columns, colnames(signature_summary) )
  if (length(missing_signature_columns) > 0L) {
    stop( paste0( "Step7G signature table is missing columns: ", paste( missing_signature_columns, collapse = ", " ) ) )
  }
  selected_signatures <- signature_summary[ signature_summary$signature_type %in% c( "high_confidence_closure_targets", "direction_consistent_targets" ) & !signature_summary$mechanistic_priority %in% c( "Do_not_prioritize_direction_discordant", "Do_not_prioritize_sample_sensitive" ), , drop = FALSE ]
  candidate_gene_rows <- list()
  if (nrow(selected_signatures) > 0L) {
    for ( signature_index in seq_len( nrow(selected_signatures) ) ) {
      current_genes <- split_pipe_genes( selected_signatures$genes[ signature_index ] )
      if (length(current_genes) == 0L) {
        next
      }
      candidate_gene_rows[[ length(candidate_gene_rows) + 1L ]] <- data.frame( signature_name = selected_signatures$signature_name[ signature_index ], TF = selected_signatures$TF[ signature_index ], regulon = selected_signatures$regulon[ signature_index ], signature_type = selected_signatures$signature_type[ signature_index ], mechanistic_priority = selected_signatures$mechanistic_priority[ signature_index ], gene = current_genes, stringsAsFactors = FALSE )
    }
  }
  if (length(candidate_gene_rows) > 0L) {
    candidate_gene_long <- unique( do.call( rbind, candidate_gene_rows ) )
  } else {
    candidate_gene_long <- data.frame( signature_name = character(0), TF = character(0), regulon = character(0), signature_type = character(0), mechanistic_priority = character(0), gene = character(0), stringsAsFactors = FALSE )
  }
  candidate_ligands <- merge( candidate_gene_long, ligand_component_map, by = "gene", all = FALSE, sort = FALSE )
  candidate_ligands$present_in_atlas <- candidate_ligands$gene %in%
    rownames(expression_data)
  candidate_ligands <- candidate_ligands[ candidate_ligands$present_in_atlas, , drop = FALSE ]
  expression_groups <- c( Microglia_Shared_prebranch = "Microglia_Shared_prebranch", Microglia_Arm_1 = "Microglia_Arm_1", Microglia_Arm_2 = "Microglia_Arm_2" )
  unique_candidate_ligand_genes <- unique( candidate_ligands$gene )
  candidate_expression_tables <- list()
  for ( group_name in names(expression_groups) ) {
    current_cells <- analysis_meta$cell[ analysis_meta$timepoint_standard == "7dpi" & analysis_meta$cellchat_group == expression_groups[[group_name]] ]
    current_summary <- summarize_gene_expression( expression_matrix = expression_data, genes = unique_candidate_ligand_genes, cells = current_cells )
    if (nrow(current_summary) == 0L) {
      next
    }
    colnames(current_summary)[ colnames(current_summary) == "mean_expression" ] <- paste0( "mean_", group_name, "_7dpi" )
    colnames(current_summary)[ colnames(current_summary) == "percent_expressing" ] <- paste0( "pct_", group_name, "_7dpi" )
    colnames(current_summary)[ colnames(current_summary) == "n_cells" ] <- paste0( "n_", group_name, "_7dpi" )
    candidate_expression_tables[[ group_name ]] <- current_summary
  }
  candidate_ligand_catalog <- candidate_ligands
  if (length(candidate_expression_tables) > 0L) {
    for ( current_expression_table in candidate_expression_tables ) {
      candidate_ligand_catalog <- merge( candidate_ligand_catalog, current_expression_table, by = "gene", all.x = TRUE, sort = FALSE )
    }
  }
  if ( all( c( "mean_Microglia_Arm_1_7dpi", "mean_Microglia_Arm_2_7dpi" ) %in% colnames(candidate_ligand_catalog) ) ) {
    candidate_ligand_catalog$delta_mean_Arm2_minus_Arm1_7dpi <- candidate_ligand_catalog$mean_Microglia_Arm_2_7dpi - candidate_ligand_catalog$mean_Microglia_Arm_1_7dpi
    candidate_ligand_catalog$expression_direction_7dpi <- ifelse( candidate_ligand_catalog$delta_mean_Arm2_minus_Arm1_7dpi > 0, "Higher_in_Arm_2", ifelse( candidate_ligand_catalog$delta_mean_Arm2_minus_Arm1_7dpi < 0, "Higher_in_Arm_1", "Equal_or_not_detected" ) )
  } else {
    candidate_ligand_catalog$delta_mean_Arm2_minus_Arm1_7dpi <- NA_real_
    candidate_ligand_catalog$expression_direction_7dpi <- NA_character_
  }
  closure_priority_columns <- intersect( c( "TF", "late_direction", "timing_classification", "onset_progress_qualified", "terminal_effect_z", "target_closure_strength", "target_closure_domains", "seven_dpi_direction_consistent", "CSI_module" ), colnames(closure_table) )
  if ( "TF" %in% closure_priority_columns ) {
    closure_priority <- unique( closure_table[ , closure_priority_columns, drop = FALSE ] )
    candidate_ligand_catalog <- merge( candidate_ligand_catalog, closure_priority, by = "TF", all.x = TRUE, sort = FALSE )
  }
  if (nrow(candidate_ligand_catalog) > 0L) {
    priority_rank <- match( candidate_ligand_catalog$mechanistic_priority, c( "Tier_1_branchpoint_coincident_target_closure", "Tier_2_downstream_Arm2_regulatory_amplifier", "Tier_3_branchpoint_coincident_weak_closure", "Exploratory_target_association" ) )
    priority_rank[ is.na(priority_rank) ] <- 99L
    candidate_ligand_catalog <- candidate_ligand_catalog[ order( priority_rank, candidate_ligand_catalog$TF, -abs( candidate_ligand_catalog$delta_mean_Arm2_minus_Arm1_7dpi ), candidate_ligand_catalog$gene ), , drop = FALSE ]
  }
  write_tsv( candidate_ligand_catalog, CANDIDATE_LIGAND_FILE )
  group_time_plot_data <- group_time_counts
  group_time_plot_data$timepoint <- factor( group_time_plot_data$timepoint, levels = FORMAL_TIMEPOINTS )
  group_time_plot_data$cellchat_group <- factor( group_time_plot_data$cellchat_group, levels = rev( sort( unique( group_time_plot_data$cellchat_group ) ) ) )
  group_time_plot_data$log10_cells <- log10( group_time_plot_data$n_cells + 1 )
  p_group_composition <- ggplot( group_time_plot_data, aes( x = timepoint, y = cellchat_group, fill = log10_cells ) ) + geom_tile( color = "white", linewidth = 0.3 ) + geom_text( aes( label = n_cells ), size = 2.7 ) + scale_fill_gradient( low = "white", high = "black", name = "log10(n+1)" ) + labs( title = "GSE162610 Step8A CellChat input composition", subtitle = paste0( "Numbers are cells; formal minimum per group/timepoint = ", MIN_CELLS_PER_GROUP ), x = NULL, y = NULL ) + theme_bw( base_size = 10 ) + theme( panel.grid = element_blank(), axis.text.y = element_text( size = 8 ) )
  ggsave( FIG_GROUP_COMPOSITION, p_group_composition, width = 9, height = max( 7, 0.35 * length( unique( group_time_plot_data$cellchat_group ) ) + 2 ), limitsize = FALSE )
  branch_time_plot_data <- group_time_counts[ group_time_counts$cellchat_group %in% c( "Microglia_Shared_prebranch", "Microglia_Arm_1", "Microglia_Arm_2" ), , drop = FALSE ]
  branch_time_plot_data$timepoint <- factor( branch_time_plot_data$timepoint, levels = FORMAL_TIMEPOINTS )
  p_branch_time <- ggplot( branch_time_plot_data, aes( x = timepoint, y = n_cells, fill = cellchat_group ) ) + geom_col( position = position_dodge( width = 0.8 ), width = 0.72 ) + geom_text( aes( label = n_cells ), position = position_dodge( width = 0.8 ), vjust = -0.3, size = 3 ) + labs( title = "Frozen neutral microglia branch cells by real timepoint", subtitle = "Counts are descriptive and rename Arm_1 or Arm_2", x = NULL, y = "Cells", fill = NULL ) + theme_bw( base_size = 11 ) + theme( panel.grid.major.x = element_blank() )
  ggsave( FIG_BRANCH_TIME, p_branch_time, width = 8, height = 5 )
  candidate_ligand_figure_created <- FALSE
  if ( nrow(candidate_ligand_catalog) > 0L && all( c( "mean_Microglia_Arm_1_7dpi", "mean_Microglia_Arm_2_7dpi" ) %in% colnames(candidate_ligand_catalog) ) ) {
    tier1_ligands <- candidate_ligand_catalog[ candidate_ligand_catalog$mechanistic_priority == "Tier_1_branchpoint_coincident_target_closure", , drop = FALSE ]
    if (nrow(tier1_ligands) == 0L) {
      tier1_ligands <- head( candidate_ligand_catalog, 30L )
    }
    tier1_ligands <- tier1_ligands[ !duplicated( paste( tier1_ligands$TF, tier1_ligands$gene, sep = "__" ) ), , drop = FALSE ]
    tier1_ligands <- head( tier1_ligands[ order( -abs( tier1_ligands$delta_mean_Arm2_minus_Arm1_7dpi ) ), , drop = FALSE ], 40L )
    if (nrow(tier1_ligands) > 0L) {
      ligand_plot_data <- rbind( data.frame( TF = tier1_ligands$TF, gene = tier1_ligands$gene, branch = "Arm_1", mean_expression = tier1_ligands$mean_Microglia_Arm_1_7dpi, stringsAsFactors = FALSE ), data.frame( TF = tier1_ligands$TF, gene = tier1_ligands$gene, branch = "Arm_2", mean_expression = tier1_ligands$mean_Microglia_Arm_2_7dpi, stringsAsFactors = FALSE ) )
      ligand_plot_data$TF_gene <- paste0( ligand_plot_data$TF, " -> ", ligand_plot_data$gene )
      ligand_order <- unique( paste0( tier1_ligands$TF, " -> ", tier1_ligands$gene ) )
      ligand_plot_data$TF_gene <- factor( ligand_plot_data$TF_gene, levels = rev( ligand_order ) )
      ligand_plot_data$branch <- factor( ligand_plot_data$branch, levels = c( "Arm_1", "Arm_2" ) )
      p_candidate_ligands <- ggplot( ligand_plot_data, aes( x = branch, y = TF_gene, fill = mean_expression ) ) + geom_tile( color = "white", linewidth = 0.3 ) + scale_fill_gradient( low = "white", high = "black", na.value = "grey90", name = "Mean log-normalized\nexpression" ) + labs( title = "Step7G candidate target-ligands in 7 dpi branch cells", subtitle = "Database/target overlap only; not a CellChat interaction result", x = NULL, y = NULL ) + theme_bw( base_size = 10 ) + theme( panel.grid = element_blank(), axis.text.y = element_text( size = 7 ) )
      ggsave( FIG_CANDIDATE_LIGANDS, p_candidate_ligands, width = 7, height = max( 6, 0.25 * length( ligand_order ) + 2 ), limitsize = FALSE )
      candidate_ligand_figure_created <- TRUE
    }
  }
  branch_7dpi_counts <- setNames( rep( 0L, 2L ), c( "Microglia_Arm_1", "Microglia_Arm_2" ) )
  branch_7dpi_observed <- group_time_counts[ group_time_counts$timepoint == "7dpi" & group_time_counts$cellchat_group %in% names(branch_7dpi_counts), , drop = FALSE ]
  branch_7dpi_counts[ branch_7dpi_observed$cellchat_group ] <- branch_7dpi_observed$n_cells
  formal_branch_7dpi_pass <- all( branch_7dpi_counts >= MIN_CELLS_PER_GROUP )
  core_recipient_groups <- c( "Astrocyte", "Fibroblast", "Pericyte_VSMC", "Endothelial", "OPC", "Oligodendrocyte", "Neuron", "Other_myeloid" )
  recipient_7dpi_counts <- setNames( rep( 0L, length( core_recipient_groups ) ), core_recipient_groups )
  recipient_7dpi_observed <- group_time_counts[ group_time_counts$timepoint == "7dpi" & group_time_counts$cellchat_group %in% core_recipient_groups, , drop = FALSE ]
  recipient_7dpi_counts[ recipient_7dpi_observed$cellchat_group ] <- recipient_7dpi_observed$n_cells
  core_recipient_groups_eligible <- names( recipient_7dpi_counts )[ recipient_7dpi_counts >= MIN_CELLS_PER_GROUP ]
  recipient_coverage_pass <- length( core_recipient_groups_eligible ) >= 6L
  timepoint_sample_pass <- all( samples_per_timepoint$n_samples >= MIN_SAMPLES_PER_TIMEPOINT )
  timepoint_group_pass <- all( eligibility_table$n_eligible_groups >= 5L )
  clean_microglia_mapping_pass <- length(clean_microglia_cells) == EXPECTED_CLEAN_MICROGLIA_CELLS && length(clean_microglia_missing_from_atlas) == 0L
  cluster_reconstruction_pass <- identical( sort( observed_cluster_keys ), sort( EXPECTED_STEP1F_CLUSTERS ) ) && all( cluster_reconstruction_audit$reconstruction_complete ) && all( !is.na( atlas_meta$broad_type_reconstructed ) )
  frozen_mapping_pass <- length(frozen_cells) == EXPECTED_FROZEN_BRANCH_CELLS && length(frozen_cells_missing_from_atlas) == 0L && length(frozen_not_in_clean_microglia) == 0L && frozen_count_exact
  timepoint_standardization_pass <- all( FORMAL_TIMEPOINTS %in% unique( atlas_meta$timepoint_standard ) )
  expression_pass <- nrow(expression_data) > 1000L && ncol(expression_data) == ncol(atlas) && expression_cell_order_exact && expression_values_finite && expression_nonnegative && length(canonical_mouse_markers_present) >= 6L
  database_overlap_pass <- database_structure_pass && nrow(interaction_table) > 500L && length(db_genes_present) > 500L && db_overlap_fraction >= 0.50
  decision_rows <- list()
  decision_rows <- add_decision_row( decision_rows, "Required input files", length(missing_files) == 0L, detail = paste0( "6/6 required files present" ) )
  decision_rows <- add_decision_row( decision_rows, "Step 7G Overall PASS", step7g_pass, detail = safe_path( STEP7G_DECISION_FILE ) )
  decision_rows <- add_decision_row( decision_rows, "Full atlas object and metadata", ncol(atlas) == EXPECTED_ATLAS_CELLS && length(missing_atlas_meta) == 0L, detail = paste0( nrow(atlas), " genes x ", ncol(atlas), " cells" ) )
  decision_rows <- add_decision_row( decision_rows, "Validated clean microglia mapped to atlas", clean_microglia_mapping_pass, detail = paste0( "clean=", length(clean_microglia_cells), "; mapped=", length(clean_microglia_in_atlas), "; missing=", length(clean_microglia_missing_from_atlas) ) )
  decision_rows <- add_decision_row( decision_rows, "Verified Step1F cluster annotation reconstructed", cluster_reconstruction_pass, detail = paste0( "observed_clusters=", length( observed_cluster_keys ), "; mapped_cells=", sum( !is.na( atlas_meta$broad_type_reconstructed ) ), "; unexpected_clusters=", length( unexpected_cluster_keys ), "; missing_expected_clusters=", length( missing_expected_cluster_keys ) ) )
  decision_rows <- add_decision_row( decision_rows, "Frozen 4203 branch cells mapped exactly", frozen_mapping_pass, detail = paste0( "Shared_prebranch=", frozen_count_vector[["Shared_prebranch"]], "; Arm_1=", frozen_count_vector[["Arm_1"]], "; Arm_2=", frozen_count_vector[["Arm_2"]], "; missing_from_atlas=", length(frozen_cells_missing_from_atlas), "; outside_clean_microglia=", length(frozen_not_in_clean_microglia) ) )
  decision_rows <- add_decision_row( decision_rows, "1/3/7 dpi timepoint standardization", timepoint_standardization_pass, detail = paste( sort( unique( atlas_meta$timepoint_standard[ !is.na( atlas_meta$timepoint_standard ) ] ) ), collapse = " | " ) )
  decision_rows <- add_decision_row( decision_rows, "Biological sample coverage per timepoint", timepoint_sample_pass, detail = paste( paste0( samples_per_timepoint$timepoint, "=", samples_per_timepoint$n_samples ), collapse = "; " ) )
  decision_rows <- add_decision_row( decision_rows, "At least five eligible groups per timepoint", timepoint_group_pass, detail = paste( paste0( eligibility_table$timepoint, "=", eligibility_table$n_eligible_groups ), collapse = "; " ) )
  decision_rows <- add_decision_row( decision_rows, "Formal 7 dpi Arm_1 and Arm_2 sender coverage", formal_branch_7dpi_pass, detail = paste( paste0( names(branch_7dpi_counts), "=", branch_7dpi_counts ), collapse = "; " ) )
  decision_rows <- add_decision_row( decision_rows, "Core 7 dpi recipient-cell coverage", recipient_coverage_pass, detail = paste0( length(core_recipient_groups_eligible), "/", length(core_recipient_groups), " groups have >=", MIN_CELLS_PER_GROUP, " cells: ", paste( core_recipient_groups_eligible, collapse = " | " ) ) )
  decision_rows <- add_decision_row( decision_rows, "Normalized RNA data suitable for CellChat", expression_pass, detail = paste0( expression_dimensions[[1L]], " genes x ", expression_dimensions[[2L]], " cells; canonical markers=", paste( canonical_mouse_markers_present, collapse = "|" ) ) )
  decision_rows <- add_decision_row( decision_rows, "CellChat required API", cellchat_api_pass, detail = paste( paste0( names(function_available), "=", function_available ), collapse = "; " ) )
  decision_rows <- add_decision_row( decision_rows, "CellChatDB.mouse structure and gene overlap", database_overlap_pass, detail = paste0( "version=", cellchat_version, "; interactions=", nrow(interaction_table), "; signaling_genes_in_atlas=", length(db_genes_present), "; overlap=", round( 100 * db_overlap_fraction, 2 ), "%" ) )
  decision_rows <- add_decision_row( decision_rows, "Non-inferential CellChat smoke test", smoke_test_pass, detail = smoke_test_message )
  decision_rows <- add_decision_row( decision_rows, "Step7G target-ligand hypothesis catalog generated", file.exists( CANDIDATE_LIGAND_FILE ), detail = paste0( nrow(candidate_ligand_catalog), " TF-signature-ligand rows; informational only" ) )
  mandatory_pass <- all( vapply( decision_rows[ seq_len( length(decision_rows) - 1L ) ], function(current_row) {
        isTRUE( current_row$pass[[1L]] )
      },
      logical(1)
    )
  )
  FINAL_STATUS <- if ( mandatory_pass ) {
    "PASS"
  } else {
    "REVIEW_REQUIRED"
  }
  decision_rows <- add_decision_row( decision_rows, "Overall Step 8A", identical( FINAL_STATUS, "PASS" ), detail = FINAL_STATUS )
  final_decision <- do.call( rbind, decision_rows )
  write_csv_utf8( final_decision, FINAL_DECISION_FILE )
  eligible_groups_by_time <- setNames( lapply( FORMAL_TIMEPOINTS, function(current_timepoint) {
        current_rows <- group_time_counts[ group_time_counts$timepoint == current_timepoint & group_time_counts$eligible_min_cells, , drop = FALSE ]
        as.character( current_rows$cellchat_group )
      }
    ), FORMAL_TIMEPOINTS
  )
  timepoint_cells <- setNames( lapply( FORMAL_TIMEPOINTS, function(current_timepoint) {
        analysis_meta$cell[ analysis_meta$timepoint_standard == current_timepoint ]
      }
    ), FORMAL_TIMEPOINTS
  )
  manifest <- list( step = "Step8A", generated_at = timestamp(), status = FINAL_STATUS, project_root = safe_path( PROJECT_ROOT ), source_files = list( atlas_rds = safe_path( ATLAS_RDS_FILE ), clean_microglia_rds = safe_path( CLEAN_MICROGLIA_RDS_FILE ), frozen_branch_metadata = safe_path( FROZEN_BRANCH_META_FILE ), step7g_decision = safe_path( STEP7G_DECISION_FILE ), step7g_signature_table = safe_path( STEP7G_SIGNATURE_TABLE_FILE ) ), formal_timepoints = FORMAL_TIMEPOINTS, minimum_cells_per_group = MIN_CELLS_PER_GROUP, Step8B_parameters = STEP8B_PARAMETERS, neutral_branch_labels = branch_to_cellchat, Step1F_cluster_to_broad_type = STEP1F_CLUSTER_TO_BROAD_TYPE, cluster_reconstruction_audit = cluster_reconstruction_audit, broad_group_mapping = broad_group_map, cell_metadata = metadata_output, group_time_counts = group_time_counts, group_sample_counts = group_sample_counts, eligible_groups_by_time = eligible_groups_by_time, timepoint_cells = timepoint_cells, formal_7dpi_sender_groups = c( "Microglia_Arm_1", "Microglia_Arm_2" ), core_recipient_groups = core_recipient_groups, candidate_regulon_ligand_catalog = candidate_ligand_catalog, database_audit = database_audit, final_decision = final_decision, interpretation_boundary = c( "Arm_1 and Arm_2 remain neutral labels.", "Step8A does not estimate communication.", "Step7G target-ligand overlap is hypothesis generation only.", "CellChat cannot establish TF-to-ligand transcriptional causality." ) )
  saveRDS( manifest, MANIFEST_RDS_FILE, compress = TRUE )
  overall_elapsed <- as.numeric( difftime( Sys.time(), overall_start, units = "mins" ) )
  branch_counts_report <- as.data.frame( frozen_count_table, stringsAsFactors = FALSE )
  colnames(branch_counts_report) <- c( "frozen_state", "n_cells" )
  capture.output( sessionInfo(), file = SESSION_FILE )
  invisible( manifest )
}
STEP8A_RESULT <- tryCatch( run_step8A(), error = function(e) {
    error_lines <- c( "======================================================================", "GSE162610 | STEP 8A ERROR", "======================================================================", paste0( "Time: ", timestamp() ), paste0( "Message: ", conditionMessage(e) ), "", paste0( "Failed call: ", paste( deparse( conditionCall(e) ), collapse = " " ) ), "", "rerun earlier Step7 analyses.", "Fix only the Step8A input/package/metadata issue reported above." )
    writeLines( error_lines, con = ERROR_FILE, useBytes = TRUE )
    stop( e )
  }
)
