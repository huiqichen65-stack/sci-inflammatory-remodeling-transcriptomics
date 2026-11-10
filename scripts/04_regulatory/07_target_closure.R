rm(list = ls())
gc()
options( stringsAsFactors = FALSE, width = 240, warn = 1, timeout = 86400 )
cat("R parser check PASS: the complete Step 7G file was parsed by source().\n")
set.seed(20260802)
PROJECT_ROOT <- Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
EXPECTED_UNIVERSE_GENES <- 7897L
EXPECTED_REGULONS <- 292L
EXPECTED_PRIMARY_REGULONS <- 215L
EXPECTED_REGULON_TARGET_EDGES <- 29609L
EXPECTED_PROGRESSIVE_GENES <- 212L
EXPECTED_STRONG_PROGRESSIVE_GENES <- 172L
EXPECTED_HIGHCONF_MARKERS <- 745L
EXPECTED_STRICT_PB_MARKERS <- 343L
EXPECTED_STEP6B_MODULE_GENES <- 803L
EXPECTED_STEP6B_MODULES <- 5L
PRIMARY_Q_MAX <- 0.05
SENSITIVITY_Q_MAX <- 0.10
MIN_FOLD_ENRICHMENT <- 1.25
MIN_PROGRESSIVE_OVERLAP <- 3L
MIN_STRONG_PROGRESSIVE_OVERLAP <- 2L
MIN_HIGHCONF_MARKER_OVERLAP <- 3L
MIN_FORMAL_MODULE_OVERLAP <- 3L
MIN_STRICT_MARKER_OVERLAP <- 1L
MIN_CORE_EDGE_EVIDENCE_COUNT <- 2L
MAX_NETWORK_TARGETS_PER_TF <- 30L
MAX_HEATMAP_CANDIDATES <- 45L
EARLY_ARM2_CLASSES <- c( "early_Arm2_fate_timing_candidate_requires_Step7G", "early_weak_Arm2_timing_candidate_requires_Step7G" )
LATE_ARM2_CLASS <- "late_Arm2_branch_associated_regulon"
LATE_ARM1_CLASS <- "late_Arm1_branch_associated_regulon"
DIRECTION_DISCORDANT_CLASS <- "RSS_dynamic_direction_discordant"
SAMPLE_SENSITIVE_CLASS <- "sample_sensitive_regulon"
FORMAL_ARM2_MODULES <- c( "M1_Arm2_terminal", "M2_Arm2_terminal" )
EXPLORATORY_ARM2_MODULES <- c( "M4_Arm2_terminal", "M5_Arm2_terminal" )
SHARED_MODULES <- c( "M3_Shared_prebranch" )
KEY_TFS <- unique( c( "Nfkb1", "Nfe2l2", "Batf", "Relb", "Nfkb2", "Fosl2", "Etv4", "Crem", "Mxi1", "Foxo3", "Taf1", "Elf4", "Brf2", "Gabpb1", "Rel", "Irf1", "Maff", "Egr1", "Mafk", "Fosb", "Atf4", "Rela", "Atf3", "Fos", "Junb", "Tfeb", "Cebpb", "Srebf1", "Gabpa", "Cux1", "Borcs8", "Bhlhe40", "Mafb", "Spi1", "Irf8", "Irf7", "Stat1" ) )
STEP7A_DIR <- file.path( PROJECT_ROOT, "04_results", "step7A_GSE162610_pySCENIC_input_environment_audit" )
STEP7A_GENE_FILE <- file.path( STEP7A_DIR, "input", "GSE162610_step7A_frozen_4203_genes.tsv" )
STEP7C_DIR <- file.path( PROJECT_ROOT, "04_results", "step7C_GSE162610_cisTarget_motif_pruning_regulon_inference" )
STEP7C_TARGET_FILE <- file.path( STEP7C_DIR, "tables", "05_regulon_targets_long.tsv" )
STEP7C_DECISION_FILE <- file.path( STEP7C_DIR, "tables", "08_final_decision.csv" )
STEP5CR_DIR <- file.path( PROJECT_ROOT, "04_results", "step5CR_GSE162610_Y29_prebranch_progressive_divergence" )
STEP5CR_RDS_FILE <- file.path( STEP5CR_DIR, "GSE162610_step5CR_results.rds" )
STEP6A_DIR <- file.path( PROJECT_ROOT, "04_results", "step6A_GSE162610_7dpi_branch_specific_DEG" )
STEP6A_HIGHCONF_FILE <- file.path( STEP6A_DIR, "tables", "02_high_confidence_branch_markers.csv" )
STEP6A_ARM1_FILE <- file.path( STEP6A_DIR, "tables", "03_Arm1_biased_high_confidence_genes.csv" )
STEP6A_ARM2_FILE <- file.path( STEP6A_DIR, "tables", "04_Arm2_biased_high_confidence_genes.csv" )
STEP6A_STRICT_FILE <- file.path( STEP6A_DIR, "tables", "05_strict_pseudobulk_supported_genes.csv" )
STEP6B_DIR <- file.path( PROJECT_ROOT, "04_results", "step6B_GSE162610_dynamic_branch_modules" )
STEP6B_GENE_ASSIGNMENT_FILE <- file.path( STEP6B_DIR, "tables", "01_dynamic_module_gene_assignments.csv" )
STEP7F_DIR <- file.path( PROJECT_ROOT, "04_results", "step7F_GSE162610_regulon_fate_timing" )
STEP7F_RDS_FILE <- file.path( STEP7F_DIR, "RDS", "GSE162610_step7F_regulon_fate_timing_results.rds" )
STEP7F_DECISION_FILE <- file.path( STEP7F_DIR, "tables", "08_final_decision.csv" )
OUTDIR <- file.path( PROJECT_ROOT, "04_results", "step7G_GSE162610_regulon_target_program_closure" )
TABDIR <- file.path(OUTDIR, "tables")
FIGDIR <- file.path(OUTDIR, "figures")
RDSDIR <- file.path(OUTDIR, "RDS")
REPORTDIR <- file.path(OUTDIR, "reports")
SIGNATUREDIR <- file.path(OUTDIR, "signatures")
for (directory in c( OUTDIR, TABDIR, FIGDIR, RDSDIR, REPORTDIR, SIGNATUREDIR )) {
  dir.create( directory, recursive = TRUE, showWarnings = FALSE )
}
INPUT_AUDIT_FILE <- file.path( TABDIR, "00_input_file_audit.csv" )
GENE_SET_AUDIT_FILE <- file.path( TABDIR, "01_evidence_gene_set_audit.tsv" )
ALL_ENRICHMENT_FILE <- file.path( TABDIR, "02_all_regulon_target_overlap_enrichment.tsv" )
ALL_CLOSURE_FILE <- file.path( TABDIR, "02b_all_regulon_target_closure_summary.tsv" )
CANDIDATE_CLOSURE_FILE <- file.path( TABDIR, "03_timing_candidate_target_closure.tsv" )
OVERLAP_GENES_FILE <- file.path( TABDIR, "04_candidate_overlap_genes_long.tsv" )
NETWORK_EDGES_FILE <- file.path( TABDIR, "05_priority_TF_target_network_edges.tsv" )
EVIDENCE_MATRIX_FILE <- file.path( TABDIR, "06_candidate_evidence_matrix.tsv" )
GO_BP_FILE <- file.path( TABDIR, "07_optional_candidate_GO_BP.tsv" )
GO_STATUS_FILE <- file.path( TABDIR, "07b_optional_GO_status.csv" )
FINAL_DECISION_FILE <- file.path( TABDIR, "08_final_decision.csv" )
SIGNATURE_TABLE_FILE <- file.path( TABDIR, "09_candidate_target_signature_summary.tsv" )
SIGNATURE_GMT_FILE <- file.path( SIGNATUREDIR, "GSE162610_step7G_candidate_target_signatures.gmt" )
RESULT_RDS_FILE <- file.path( RDSDIR, "GSE162610_step7G_regulon_target_program_closure_results.rds" )
REPORT_FILE <- file.path( REPORTDIR, "Step7G_report.txt" )
SESSION_FILE <- file.path( REPORTDIR, "Step7G_sessionInfo.txt" )
ERROR_FILE <- file.path( REPORTDIR, "Step7G_ERROR.txt" )
FIG_ENRICHMENT_HEATMAP <- file.path( FIGDIR, "Fig7G_01_candidate_enrichment_heatmap.pdf" )
FIG_EVIDENCE_SCORE <- file.path( FIGDIR, "Fig7G_02_candidate_closure_score.pdf" )
FIG_OVERLAP_COUNTS <- file.path( FIGDIR, "Fig7G_03_candidate_overlap_counts.pdf" )
FIG_TARGET_EDGES <- file.path( FIGDIR, "Fig7G_04_priority_target_edges.pdf" )
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
  writeLines( c( "======================================================================", "GSE162610 | STEP 7G ERROR", "======================================================================", paste0("Time: ", timestamp()), "", text ), con = ERROR_FILE, useBytes = TRUE )
  stop( text, call. = FALSE )
}
as_logical_safe <- function(x) {
  if (is.logical(x)) {
    return(x)
  }
  normalized <- toupper( trimws( as.character(x) ) )
  output <- rep( NA, length(normalized) )
  output[ normalized %in% c( "TRUE", "T", "1", "YES", "Y" ) ] <- TRUE
  output[ normalized %in% c( "FALSE", "F", "0", "NO", "N" ) ] <- FALSE
  output
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
    value <- as_logical_safe( decision_table$pass[ index[[1]] ] )
    return( isTRUE( value[[1]] ) )
  }
  if ("status" %in% colnames(decision_table)) {
    return( identical( toupper( as.character( decision_table$status[ index[[1]] ] ) ), "PASS" ) )
  }
  FALSE
}
read_gene_universe <- function(path) {
  lines <- readLines( path, warn = FALSE, encoding = "UTF-8" )
  lines <- trimws( lines )
  lines <- lines[ nzchar(lines) ]
  genes <- vapply( strsplit( lines, "\t", fixed = TRUE ), function(parts) {
      parts[[1]]
    },
    FUN.VALUE = character(1)
  )
  header_like <- tolower( genes ) %in%
    c( "gene", "genes", "gene_id", "feature", "features" )
  genes <- genes[ !header_like ]
  unique( genes[ nzchar(genes) ] )
}
require_columns <- function( data, columns, data_name ) {
  missing_columns <- setdiff( columns, colnames(data) )
  if (length(missing_columns) > 0L) {
    stop( paste0( data_name, " is missing required columns: ", paste( missing_columns, collapse = ", " ), ". Current columns: ", paste( colnames(data), collapse = ", " ) ), call. = FALSE )
  }
}
hypergeometric_overlap <- function( target_genes, evidence_genes, universe_genes ) {
  target_genes <- unique( intersect( target_genes, universe_genes ) )
  evidence_genes <- unique( intersect( evidence_genes, universe_genes ) )
  overlap_genes <- intersect( target_genes, evidence_genes )
  universe_n <- length( universe_genes )
  target_n <- length( target_genes )
  evidence_n <- length( evidence_genes )
  overlap_n <- length( overlap_genes )
  expected_overlap <- if ( universe_n > 0L ) {
    target_n * evidence_n / universe_n
  } else {
    NA_real_
  }
  fold_enrichment <- if ( is.finite(expected_overlap) && expected_overlap > 0 ) {
    overlap_n / expected_overlap
  } else {
    NA_real_
  }
  a <- overlap_n
  b <- target_n - overlap_n
  c_value <- evidence_n - overlap_n
  d <- universe_n - a - b - c_value
  odds_ratio_haldane <- if ( all( c( a, b, c_value, d ) >= 0 ) ) {
    ( ( a + 0.5 ) * ( d + 0.5 ) ) / ( ( b + 0.5 ) * ( c_value + 0.5 ) )
  } else {
    NA_real_
  }
  p_value <- if ( universe_n > 0L && target_n <= universe_n && evidence_n <= universe_n ) {
    stats::phyper( q = overlap_n - 1L, m = evidence_n, n = universe_n - evidence_n, k = target_n, lower.tail = FALSE )
  } else {
    NA_real_
  }
  list( universe_n = universe_n, target_n = target_n, evidence_n = evidence_n, overlap_n = overlap_n, expected_overlap = expected_overlap, fold_enrichment = fold_enrichment, odds_ratio_haldane = odds_ratio_haldane, p_value = p_value, overlap_genes = overlap_genes )
}
get_enrichment_row <- function( enrichment_table, regulon_name, evidence_name ) {
  index <- which( enrichment_table$regulon == regulon_name & enrichment_table$evidence_set == evidence_name )
  if (length(index) == 0L) {
    return(NULL)
  }
  enrichment_table[ index[[1]], , drop = FALSE ]
}
is_supported_overlap <- function( enrichment_row, minimum_overlap, q_max = PRIMARY_Q_MAX ) {
  if ( is.null(enrichment_row) || nrow(enrichment_row) == 0L ) {
    return(FALSE)
  }
  isTRUE( is.finite( enrichment_row$q_within_evidence_set ) && enrichment_row$q_within_evidence_set <= q_max && is.finite( enrichment_row$fold_enrichment ) && enrichment_row$fold_enrichment >= MIN_FOLD_ENRICHMENT && enrichment_row$overlap_n >= minimum_overlap )
}
safe_value <- function( enrichment_row, column_name, default = NA ) {
  if ( is.null(enrichment_row) || nrow(enrichment_row) == 0L || !column_name %in% colnames(enrichment_row) ) {
    return(default)
  }
  enrichment_row[[column_name]][[1]]
}
direction_to_suffix <- function(direction) {
  if ( identical( direction, "Arm_1" ) ) {
    return("Arm1")
  }
  if ( identical( direction, "Arm_2" ) ) {
    return("Arm2")
  }
  "Neutral"
}
sanitize_signature_name <- function(x) {
  gsub( "[^A-Za-z0-9_.-]", "_", x )
}
overall_start <- Sys.time()
msg("==============================================================")
msg("GSE162610 Step 7G | regulon target-program closure")
msg("==============================================================")
msg("Background universe: frozen Step7A 7,897 expressed genes.")
msg("Primary candidates: Step7F timing candidates and key regulons.")
required_files <- c( STEP7A_GENE_FILE, STEP7C_TARGET_FILE, STEP7C_DECISION_FILE, STEP5CR_RDS_FILE, STEP6A_HIGHCONF_FILE, STEP6A_ARM1_FILE, STEP6A_ARM2_FILE, STEP6A_STRICT_FILE, STEP6B_GENE_ASSIGNMENT_FILE, STEP7F_RDS_FILE, STEP7F_DECISION_FILE )
missing_files <- required_files[ !file.exists( required_files ) ]
if (length(missing_files) > 0L) {
  abort_step( paste0( "Required Step7G files are missing:\n", paste( safe_path( missing_files ), collapse = "\n" ) ) )
}
input_audit <- data.frame( role = c( "Step7A_frozen_gene_universe", "Step7C_regulon_targets", "Step7C_final_decision", "Step5CR_results", "Step6A_highconfidence_markers", "Step6A_Arm1_markers", "Step6A_Arm2_markers", "Step6A_strict_markers", "Step6B_gene_assignments", "Step7F_complete_results", "Step7F_final_decision" ), path = safe_path( required_files ), exists = file.exists( required_files ), size_bytes = file.info( required_files )$size, stringsAsFactors = FALSE )
write_csv( input_audit, INPUT_AUDIT_FILE )
step7c_decision <- data.table::fread( STEP7C_DECISION_FILE, data.table = FALSE, showProgress = FALSE )
step7f_decision <- data.table::fread( STEP7F_DECISION_FILE, data.table = FALSE, showProgress = FALSE )
step7c_pass <- decision_pass( step7c_decision, "Overall Step 7C" )
step7f_pass <- decision_pass( step7f_decision, "Overall Step 7F" )
if (!step7c_pass) {
  abort_step( "Step 7C does not contain an Overall PASS decision." )
}
if (!step7f_pass) {
  abort_step( "Step 7F does not contain an Overall PASS decision." )
}
msg("Loading the frozen Step7A gene universe...")
universe_genes <- read_gene_universe( STEP7A_GENE_FILE )
if ( length(universe_genes) != EXPECTED_UNIVERSE_GENES || anyDuplicated( universe_genes ) > 0L ) {
  abort_step( paste0( "Unexpected Step7A gene universe: ", length(universe_genes), " unique genes; expected ", EXPECTED_UNIVERSE_GENES, "." ) )
}
msg("Loading Step7C motif-supported regulon targets...")
regulon_targets <- data.table::fread( STEP7C_TARGET_FILE, data.table = FALSE, showProgress = FALSE )
require_columns( regulon_targets, c( "regulon", "TF", "sign", "target", "weight" ), "Step7C regulon-target table" )
regulon_targets$regulon <- as.character( regulon_targets$regulon )
regulon_targets$TF <- as.character( regulon_targets$TF )
regulon_targets$target <- as.character( regulon_targets$target )
regulon_targets$weight <- suppressWarnings( as.numeric( regulon_targets$weight ) )
if ( nrow(regulon_targets) != EXPECTED_REGULON_TARGET_EDGES || length( unique( regulon_targets$regulon ) ) != EXPECTED_REGULONS || any( !is.finite( regulon_targets$weight ) ) || any( regulon_targets$weight < 0 ) || anyDuplicated( paste( regulon_targets$regulon, regulon_targets$target, sep = "||" ) ) > 0L || any( !regulon_targets$target %in% universe_genes ) ) {
  abort_step( paste0( "Step7C regulon-target audit failed. Rows=", nrow(regulon_targets), "; regulons=", length( unique( regulon_targets$regulon ) ), "." ) )
}
regulon_targets$self_target <- regulon_targets$TF == regulon_targets$target
regulon_targets$weight_rank_within_regulon <- ave( -regulon_targets$weight, regulon_targets$regulon, FUN = function(x) {
      rank( x, ties.method = "min" )
    }
  )
regulon_targets$weight_percentile_within_regulon <- ave( regulon_targets$weight, regulon_targets$regulon, FUN = function(x) {
      rank( x, ties.method = "average" ) / length(x)
    }
  )
target_split <- split( regulon_targets, regulon_targets$regulon )
msg("Loading Step5C-R progressive-gene evidence...")
step5cr <- readRDS( STEP5CR_RDS_FILE )
if ( !"master_table" %in% names(step5cr) ) {
  abort_step( "Step5C-R RDS does not contain master_table." )
}
step5_master <- as.data.frame( step5cr$master_table )
require_columns( step5_master, c( "gene_id", "terminal_direction" ), "Step5C-R master_table" )
progressive_flag <- if ( "progressive_shape_significant" %in% colnames(step5_master) ) {
  as_logical_safe( step5_master$progressive_shape_significant )
} else if (
  "progressive_shape_q_value" %in%
    colnames(step5_master)
) {
  is.finite( step5_master$progressive_shape_q_value ) & step5_master$progressive_shape_q_value <= 0.05
} else {
  abort_step( "Step5C-R master_table has no progressive-shape flag or q-value." )
}
strong_column_candidates <- c( "strong_progressive", "strong_progressive_pattern" )
strong_column_hits <- strong_column_candidates[ strong_column_candidates %in% colnames(step5_master) ]
if (length(strong_column_hits) == 0L) {
  abort_step( "Step5C-R master_table has no strong-progressive column." )
}
strong_column <- strong_column_hits[1]
strong_flag <- as_logical_safe( step5_master[[strong_column]] )
step5_master$gene_id <- as.character( step5_master$gene_id )
step5_master$terminal_direction <- as.character( step5_master$terminal_direction )
progressive_genes <- unique( step5_master$gene_id[ progressive_flag %in% TRUE ] )
strong_progressive_genes <- unique( step5_master$gene_id[ strong_flag %in% TRUE ] )
progressive_arm1 <- unique( step5_master$gene_id[ progressive_flag %in% TRUE & step5_master$terminal_direction == "Arm_1" ] )
progressive_arm2 <- unique( step5_master$gene_id[ progressive_flag %in% TRUE & step5_master$terminal_direction == "Arm_2" ] )
strong_arm1 <- unique( step5_master$gene_id[ strong_flag %in% TRUE & step5_master$terminal_direction == "Arm_1" ] )
strong_arm2 <- unique( step5_master$gene_id[ strong_flag %in% TRUE & step5_master$terminal_direction == "Arm_2" ] )
if ( length(progressive_genes) != EXPECTED_PROGRESSIVE_GENES || length(strong_progressive_genes) != EXPECTED_STRONG_PROGRESSIVE_GENES ) {
  abort_step( paste0( "Unexpected Step5C-R counts. Progressive=", length(progressive_genes), "; strong=", length(strong_progressive_genes), "." ) )
}
msg("Loading Step6A branch-marker evidence...")
step6a_highconf <- data.table::fread( STEP6A_HIGHCONF_FILE, data.table = FALSE, showProgress = FALSE )
step6a_arm1 <- data.table::fread( STEP6A_ARM1_FILE, data.table = FALSE, showProgress = FALSE )
step6a_arm2 <- data.table::fread( STEP6A_ARM2_FILE, data.table = FALSE, showProgress = FALSE )
step6a_strict <- data.table::fread( STEP6A_STRICT_FILE, data.table = FALSE, showProgress = FALSE )
for ( table_name in c( "step6a_highconf", "step6a_arm1", "step6a_arm2", "step6a_strict" ) ) {
  current_table <- get( table_name )
  require_columns( current_table, "gene_id", table_name )
}
step6a_highconf$gene_id <- as.character( step6a_highconf$gene_id )
step6a_arm1$gene_id <- as.character( step6a_arm1$gene_id )
step6a_arm2$gene_id <- as.character( step6a_arm2$gene_id )
step6a_strict$gene_id <- as.character( step6a_strict$gene_id )
highconf_genes <- unique( step6a_highconf$gene_id )
arm1_highconf_genes <- unique( step6a_arm1$gene_id )
arm2_highconf_genes <- unique( step6a_arm2$gene_id )
strict_genes <- unique( step6a_strict$gene_id )
strict_direction_lookup <- NULL
if ( "final_direction" %in% colnames(step6a_strict) ) {
  strict_direction_lookup <- setNames( as.character( step6a_strict$final_direction ), step6a_strict$gene_id )
} else if (
  "final_direction" %in%
    colnames(step6a_highconf)
) {
  strict_direction_lookup <- setNames( as.character( step6a_highconf$final_direction ), step6a_highconf$gene_id )
}
if (is.null(strict_direction_lookup)) {
  strict_arm1 <- intersect( strict_genes, arm1_highconf_genes )
  strict_arm2 <- intersect( strict_genes, arm2_highconf_genes )
} else {
  strict_arm1 <- strict_genes[ strict_direction_lookup[ strict_genes ] == "Arm_1-biased" ]
  strict_arm2 <- strict_genes[ strict_direction_lookup[ strict_genes ] == "Arm_2-biased" ]
}
strict_arm1 <- unique( strict_arm1[ !is.na(strict_arm1) ] )
strict_arm2 <- unique( strict_arm2[ !is.na(strict_arm2) ] )
if ( length(highconf_genes) != EXPECTED_HIGHCONF_MARKERS || length(strict_genes) != EXPECTED_STRICT_PB_MARKERS ) {
  abort_step( paste0( "Unexpected Step6A marker counts. High-confidence=", length(highconf_genes), "; strict=", length(strict_genes), "." ) )
}
msg("Loading Step6B dynamic module assignments...")
step6b_assignments <- data.table::fread( STEP6B_GENE_ASSIGNMENT_FILE, data.table = FALSE, showProgress = FALSE )
require_columns( step6b_assignments, c( "gene_id", "module_name" ), "Step6B gene assignments" )
step6b_assignments$gene_id <- as.character( step6b_assignments$gene_id )
step6b_assignments$module_name <- as.character( step6b_assignments$module_name )
if ( length( unique( step6b_assignments$gene_id ) ) != EXPECTED_STEP6B_MODULE_GENES || length( unique( step6b_assignments$module_name ) ) != EXPECTED_STEP6B_MODULES ) {
  abort_step( paste0( "Unexpected Step6B assignment counts. Genes=", length( unique( step6b_assignments$gene_id ) ), "; modules=", length( unique( step6b_assignments$module_name ) ), "." ) )
}
step6b_module_sets <- split( step6b_assignments$gene_id, step6b_assignments$module_name )
formal_arm2_module_genes <- unique( unlist( step6b_module_sets[ intersect( FORMAL_ARM2_MODULES, names( step6b_module_sets ) ) ], use.names = FALSE ) )
exploratory_arm2_module_genes <- unique( unlist( step6b_module_sets[ intersect( EXPLORATORY_ARM2_MODULES, names( step6b_module_sets ) ) ], use.names = FALSE ) )
shared_module_genes <- unique( unlist( step6b_module_sets[ intersect( SHARED_MODULES, names( step6b_module_sets ) ) ], use.names = FALSE ) )
all_dynamic_module_genes <- unique( step6b_assignments$gene_id )
evidence_sets <- list( Progressive_all = progressive_genes, Progressive_Arm1 = progressive_arm1, Progressive_Arm2 = progressive_arm2, Strong_progressive_all = strong_progressive_genes, Strong_progressive_Arm1 = strong_arm1, Strong_progressive_Arm2 = strong_arm2, Step6A_highconf_all = highconf_genes, Step6A_Arm1_highconf = arm1_highconf_genes, Step6A_Arm2_highconf = arm2_highconf_genes, Step6A_strict_all = strict_genes, Step6A_Arm1_strict = strict_arm1, Step6A_Arm2_strict = strict_arm2, Step6B_M1_Arm2_terminal = step6b_module_sets[["M1_Arm2_terminal"]], Step6B_M2_Arm2_terminal = step6b_module_sets[["M2_Arm2_terminal"]], Step6B_formal_Arm2_M1_M2 = formal_arm2_module_genes, Step6B_M3_Shared_prebranch = shared_module_genes, Step6B_exploratory_Arm2_M4_M5 = exploratory_arm2_module_genes, Step6B_all_dynamic = all_dynamic_module_genes )
evidence_metadata <- data.frame( evidence_set = names( evidence_sets ), evidence_level = c( "Step5C-R", "Step5C-R", "Step5C-R", "Step5C-R", "Step5C-R", "Step5C-R", "Step6A", "Step6A", "Step6A", "Step6A", "Step6A", "Step6A", "Step6B", "Step6B", "Step6B", "Step6B", "Step6B", "Step6B" ), expected_direction = c( "Any", "Arm_1", "Arm_2", "Any", "Arm_1", "Arm_2", "Any", "Arm_1", "Arm_2", "Any", "Arm_1", "Arm_2", "Arm_2", "Arm_2", "Arm_2", "Shared_prebranch", "Arm_2", "Any" ), evidence_role = c( "progressive", "progressive_directional", "progressive_directional", "strong_progressive", "strong_progressive_directional", "strong_progressive_directional", "terminal_marker", "terminal_marker_directional", "terminal_marker_directional", "strict_terminal_marker", "strict_terminal_marker_directional", "strict_terminal_marker_directional", "formal_terminal_module", "formal_terminal_module", "formal_terminal_module_combined", "shared_identity_module", "exploratory_terminal_module", "all_dynamic_modules" ), formal_priority = c( TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE ), stringsAsFactors = FALSE )
gene_set_audit_parts <- vector( "list", length(evidence_sets) )
for ( evidence_index in seq_along(evidence_sets) ) {
  evidence_name <- names( evidence_sets )[ evidence_index ]
  original_genes <- unique( as.character( evidence_sets[[evidence_name]] ) )
  retained_genes <- intersect( original_genes, universe_genes )
  evidence_sets[[evidence_name]] <- retained_genes
  gene_set_audit_parts[[evidence_index]] <- data.frame( evidence_set = evidence_name, original_n = length( original_genes ), in_universe_n = length( retained_genes ), outside_universe_n = length( setdiff( original_genes, universe_genes ) ), universe_fraction = length( retained_genes ) / length( universe_genes ), stringsAsFactors = FALSE )
}
gene_set_audit <- do.call( rbind, gene_set_audit_parts )
gene_set_audit <- merge( evidence_metadata, gene_set_audit, by = "evidence_set", all.x = TRUE, sort = FALSE )
gene_set_audit <- gene_set_audit[ match( evidence_metadata$evidence_set, gene_set_audit$evidence_set ), , drop = FALSE ]
write_tsv( gene_set_audit, GENE_SET_AUDIT_FILE )
msg("Loading Step7F timing results...")
step7f <- readRDS( STEP7F_RDS_FILE )
required_step7f_fields <- c( "regulon_timing_master", "key_candidate_timing", "final_decision" )
missing_step7f_fields <- setdiff( required_step7f_fields, names(step7f) )
if (length(missing_step7f_fields) > 0L) {
  abort_step( paste0( "Step7F RDS is missing: ", paste( missing_step7f_fields, collapse = ", " ) ) )
}
if ( !decision_pass( step7f$final_decision, "Overall Step 7F" ) ) {
  abort_step( "Step7F RDS does not contain an Overall PASS decision." )
}
timing_master <- as.data.frame( step7f$regulon_timing_master )
require_columns( timing_master, c( "regulon", "TF", "n_targets", "late_direction", "onset_progress_qualified", "terminal_effect_z", "seven_dpi_direction_consistent", "sample_direction_consistency_fraction", "CSI_module", "timing_classification", "formal_terminal_reference_progress", "lead_before_formal_terminal_program" ), "Step7F regulon timing master" )
timing_master$regulon <- as.character( timing_master$regulon )
timing_master$TF <- as.character( timing_master$TF )
timing_master$late_direction <- as.character( timing_master$late_direction )
timing_master$timing_classification <- as.character( timing_master$timing_classification )
timing_master$seven_dpi_direction_consistent <- as_logical_safe( timing_master$seven_dpi_direction_consistent )
if ( nrow(timing_master) != EXPECTED_PRIMARY_REGULONS || any( !timing_master$regulon %in% names( target_split ) ) ) {
  abort_step( "Step7F primary regulons do not match Step7C target sets." )
}
enrichment_start <- Sys.time()
msg( "Testing ", nrow(timing_master), " primary regulons against ", length(evidence_sets), " evidence gene sets..." )
enrichment_parts <- vector( "list", nrow(timing_master) * length(evidence_sets) )
enrichment_counter <- 0L
for ( regulon_index in seq_len( nrow(timing_master) ) ) {
  regulon_name <- timing_master$regulon[ regulon_index ]
  target_table <- target_split[[regulon_name]]
  target_genes <- unique( target_table$target )
  total_weight <- sum( target_table$weight )
  for ( evidence_index in seq_along(evidence_sets) ) {
    evidence_name <- names( evidence_sets )[ evidence_index ]
    overlap_result <- hypergeometric_overlap( target_genes = target_genes, evidence_genes = evidence_sets[[evidence_name]], universe_genes = universe_genes )
    overlap_target_rows <- target_table[ target_table$target %in% overlap_result$overlap_genes, , drop = FALSE ]
    overlap_weight_sum <- sum( overlap_target_rows$weight )
    enrichment_counter <- enrichment_counter + 1L
    enrichment_parts[[enrichment_counter]] <- data.frame( regulon = regulon_name, TF = timing_master$TF[ regulon_index ], late_direction = timing_master$late_direction[ regulon_index ], timing_classification = timing_master$timing_classification[ regulon_index ], evidence_set = evidence_name, evidence_level = evidence_metadata$evidence_level[ evidence_index ], expected_direction = evidence_metadata$expected_direction[ evidence_index ], evidence_role = evidence_metadata$evidence_role[ evidence_index ], formal_priority = evidence_metadata$formal_priority[ evidence_index ], universe_n = overlap_result$universe_n, target_n = overlap_result$target_n, evidence_n = overlap_result$evidence_n, overlap_n = overlap_result$overlap_n, expected_overlap = overlap_result$expected_overlap, fold_enrichment = overlap_result$fold_enrichment, odds_ratio_haldane = overlap_result$odds_ratio_haldane, p_value = overlap_result$p_value, target_weight_total = total_weight, overlap_weight_sum = overlap_weight_sum, overlap_weight_fraction = if ( total_weight > 0 ) {
          overlap_weight_sum / total_weight
        } else {
          NA_real_
        },
        overlap_median_weight = if ( nrow(overlap_target_rows) > 0L ) {
          stats::median( overlap_target_rows$weight )
        } else {
          NA_real_
        },
        overlap_top_weight_quartile_n = sum( overlap_target_rows$weight_percentile_within_regulon >= 0.75, na.rm = TRUE ), overlap_genes = paste( sort( overlap_result$overlap_genes ), collapse = " | " ), stringsAsFactors = FALSE
      )
  }
  if ( regulon_index %% 25L == 0L ) {
    msg( "  completed ", regulon_index, "/", nrow(timing_master), " regulons" )
  }
}
all_enrichment <- do.call( rbind, enrichment_parts )
all_enrichment$q_global <- stats::p.adjust( all_enrichment$p_value, method = "BH" )
all_enrichment$q_within_evidence_set <- ave( all_enrichment$p_value, all_enrichment$evidence_set, FUN = function(values) {
      stats::p.adjust( values, method = "BH" )
    }
  )
all_enrichment$direction_match <- all_enrichment$expected_direction == "Any" | all_enrichment$expected_direction == all_enrichment$late_direction | ( all_enrichment$expected_direction == "Shared_prebranch" & all_enrichment$late_direction == "Arm_1" )
all_enrichment <- all_enrichment[ order( all_enrichment$regulon, all_enrichment$q_within_evidence_set, -all_enrichment$fold_enrichment ), , drop = FALSE ]
write_tsv( all_enrichment, ALL_ENRICHMENT_FILE )
enrichment_elapsed <- as.numeric( difftime( Sys.time(), enrichment_start, units = "secs" ) )
msg( "Overlap enrichment complete in ", format_duration( enrichment_elapsed ), "." )
candidate_pool_index <- timing_master$timing_classification %in%
    c( EARLY_ARM2_CLASSES, LATE_ARM2_CLASS, LATE_ARM1_CLASS, DIRECTION_DISCORDANT_CLASS, SAMPLE_SENSITIVE_CLASS ) | timing_master$TF %in%
    KEY_TFS
candidate_pool <- timing_master[ candidate_pool_index, , drop = FALSE ]
candidate_rows <- vector( "list", nrow(candidate_pool) )
for ( candidate_index in seq_len( nrow(candidate_pool) ) ) {
  current <- candidate_pool[ candidate_index, , drop = FALSE ]
  regulon_name <- current$regulon[[1]]
  late_direction <- current$late_direction[[1]]
  direction_suffix <- direction_to_suffix( late_direction )
  progressive_set_name <- paste0( "Progressive_", direction_suffix )
  strong_set_name <- paste0( "Strong_progressive_", direction_suffix )
  highconf_set_name <- paste0( "Step6A_", direction_suffix, "_highconf" )
  strict_set_name <- paste0( "Step6A_", direction_suffix, "_strict" )
  progressive_row <- get_enrichment_row( all_enrichment, regulon_name, progressive_set_name )
  strong_row <- get_enrichment_row( all_enrichment, regulon_name, strong_set_name )
  highconf_row <- get_enrichment_row( all_enrichment, regulon_name, highconf_set_name )
  strict_row <- get_enrichment_row( all_enrichment, regulon_name, strict_set_name )
  formal_module_row <- if ( late_direction == "Arm_2" ) {
    get_enrichment_row( all_enrichment, regulon_name, "Step6B_formal_Arm2_M1_M2" )
  } else {
    get_enrichment_row( all_enrichment, regulon_name, "Step6B_M3_Shared_prebranch" )
  }
  progressive_support <- is_supported_overlap( progressive_row, MIN_PROGRESSIVE_OVERLAP )
  strong_support <- is_supported_overlap( strong_row, MIN_STRONG_PROGRESSIVE_OVERLAP )
  highconf_support <- is_supported_overlap( highconf_row, MIN_HIGHCONF_MARKER_OVERLAP )
  module_support <- is_supported_overlap( formal_module_row, MIN_FORMAL_MODULE_OVERLAP )
  strict_support <- is_supported_overlap( strict_row, MIN_STRICT_MARKER_OVERLAP, q_max = SENSITIVITY_Q_MAX )
  target_closure_domains <- sum( c( progressive_support, strong_support, highconf_support, module_support ) )
  target_closure_strength <- if ( target_closure_domains >= 3L && ( progressive_support || strong_support ) && ( highconf_support || module_support ) ) {
    "strong"
  } else if (
    target_closure_domains >= 2L
  ) {
    "moderate"
  } else if (
    target_closure_domains >= 1L || strict_support
  ) {
    "weak"
  } else {
    "absent"
  }
  early_arm2_timing <- current$timing_classification[[1]] %in%
    EARLY_ARM2_CLASSES
  sample_consistent <- isTRUE( current$seven_dpi_direction_consistent[[1]] )
  earlier_than_terminal_program <- is.finite( current$lead_before_formal_terminal_program[[1]] ) && current$lead_before_formal_terminal_program[[1]] >
      0.02
  internal_fate_driver_criteria_met <- early_arm2_timing && sample_consistent && earlier_than_terminal_program && target_closure_strength == "strong"
  mechanistic_priority <- if ( current$timing_classification[[1]] == SAMPLE_SENSITIVE_CLASS ) {
    "Do_not_prioritize_sample_sensitive"
  } else if (
    current$timing_classification[[1]] == DIRECTION_DISCORDANT_CLASS
  ) {
    "Do_not_prioritize_direction_discordant"
  } else if (
    early_arm2_timing && target_closure_strength == "strong"
  ) {
    "Tier_1_branchpoint_coincident_target_closure"
  } else if (
    early_arm2_timing && target_closure_strength == "moderate"
  ) {
    "Tier_2_branchpoint_coincident_partial_closure"
  } else if (
    early_arm2_timing
  ) {
    "Tier_3_branchpoint_coincident_weak_closure"
  } else if (
    current$timing_classification[[1]] == LATE_ARM2_CLASS && target_closure_strength %in%
      c( "strong", "moderate" )
  ) {
    "Tier_2_downstream_Arm2_regulatory_amplifier"
  } else if (
    current$timing_classification[[1]] == LATE_ARM1_CLASS && target_closure_strength %in%
      c( "strong", "moderate" )
  ) {
    "Arm1_identity_retention_target_association"
  } else {
    "Exploratory_target_association"
  }
  candidate_rows[[candidate_index]] <- data.frame( regulon = regulon_name, TF = current$TF[[1]], n_targets = current$n_targets[[1]], timing_classification = current$timing_classification[[1]], late_direction = late_direction, onset_progress_qualified = current$onset_progress_qualified[[1]], terminal_effect_z = current$terminal_effect_z[[1]], formal_terminal_reference_progress = current$formal_terminal_reference_progress[[1]], lead_before_formal_terminal_program = current$lead_before_formal_terminal_program[[1]], earlier_than_terminal_program = earlier_than_terminal_program, seven_dpi_direction_consistent = sample_consistent, sample_direction_consistency_fraction = current$sample_direction_consistency_fraction[[1]], CSI_module = current$CSI_module[[1]], progressive_set = progressive_set_name, progressive_overlap_n = safe_value( progressive_row, "overlap_n", 0L ), progressive_fold_enrichment = safe_value( progressive_row, "fold_enrichment" ), progressive_q = safe_value( progressive_row, "q_within_evidence_set" ), progressive_support = progressive_support, strong_progressive_set = strong_set_name, strong_progressive_overlap_n = safe_value( strong_row, "overlap_n", 0L ), strong_progressive_fold_enrichment = safe_value( strong_row, "fold_enrichment" ), strong_progressive_q = safe_value( strong_row, "q_within_evidence_set" ), strong_progressive_support = strong_support, highconf_marker_set = highconf_set_name, highconf_marker_overlap_n = safe_value( highconf_row, "overlap_n", 0L ), highconf_marker_fold_enrichment = safe_value( highconf_row, "fold_enrichment" ), highconf_marker_q = safe_value( highconf_row, "q_within_evidence_set" ), highconf_marker_support = highconf_support, strict_marker_set = strict_set_name, strict_marker_overlap_n = safe_value( strict_row, "overlap_n", 0L ), strict_marker_fold_enrichment = safe_value( strict_row, "fold_enrichment" ), strict_marker_q = safe_value( strict_row, "q_within_evidence_set" ), strict_marker_support = strict_support, formal_module_set = if ( late_direction == "Arm_2" ) {
        "Step6B_formal_Arm2_M1_M2"
      } else {
        "Step6B_M3_Shared_prebranch"
      },
      formal_module_overlap_n = safe_value( formal_module_row, "overlap_n", 0L ), formal_module_fold_enrichment = safe_value( formal_module_row, "fold_enrichment" ), formal_module_q = safe_value( formal_module_row, "q_within_evidence_set" ), formal_module_support = module_support, target_closure_domains = target_closure_domains, target_closure_strength = target_closure_strength, internal_fate_driver_criteria_met = internal_fate_driver_criteria_met, mechanistic_priority = mechanistic_priority, stringsAsFactors = FALSE
    )
}
candidate_closure <- do.call( rbind, candidate_rows )
priority_order <- c( "Tier_1_branchpoint_coincident_target_closure", "Tier_2_branchpoint_coincident_partial_closure", "Tier_2_downstream_Arm2_regulatory_amplifier", "Tier_3_branchpoint_coincident_weak_closure", "Arm1_identity_retention_target_association", "Exploratory_target_association", "Do_not_prioritize_direction_discordant", "Do_not_prioritize_sample_sensitive" )
candidate_closure <- candidate_closure[ order( match( candidate_closure$mechanistic_priority, priority_order ), is.na( candidate_closure$onset_progress_qualified ), candidate_closure$onset_progress_qualified, -candidate_closure$target_closure_domains, -abs( candidate_closure$terminal_effect_z ) ), , drop = FALSE ]
write_tsv( candidate_closure, CANDIDATE_CLOSURE_FILE )
all_closure <- timing_master
selected_summary_sets <- c( "Progressive_Arm1", "Progressive_Arm2", "Strong_progressive_Arm1", "Strong_progressive_Arm2", "Step6A_Arm1_highconf", "Step6A_Arm2_highconf", "Step6A_Arm1_strict", "Step6A_Arm2_strict", "Step6B_formal_Arm2_M1_M2", "Step6B_M3_Shared_prebranch" )
for (evidence_name in selected_summary_sets) {
  subset_table <- all_enrichment[ all_enrichment$evidence_set == evidence_name, c( "regulon", "overlap_n", "fold_enrichment", "q_within_evidence_set", "overlap_weight_fraction" ), drop = FALSE ]
  suffix <- gsub( "[^A-Za-z0-9]+", "_", evidence_name )
  colnames(subset_table)[ -1L ] <- paste0( colnames(subset_table)[ -1L ], "__", suffix )
  all_closure <- merge( all_closure, subset_table, by = "regulon", all.x = TRUE, sort = FALSE )
}
all_closure <- all_closure[ match( timing_master$regulon, all_closure$regulon ), , drop = FALSE ]
write_tsv( all_closure, ALL_CLOSURE_FILE )
msg("Annotating candidate target genes with Step5/Step6 evidence...")
step5_direction_lookup <- setNames( step5_master$terminal_direction, step5_master$gene_id )
step6a_direction_lookup <- rep( NA_character_, length( universe_genes ) )
names(step6a_direction_lookup) <- universe_genes
step6a_direction_lookup[ arm1_highconf_genes ] <- "Arm_1"
step6a_direction_lookup[ arm2_highconf_genes ] <- "Arm_2"
step6b_module_lookup <- setNames( step6b_assignments$module_name, step6b_assignments$gene_id )
overlap_gene_parts <- vector( "list", nrow(candidate_closure) )
for ( candidate_index in seq_len( nrow(candidate_closure) ) ) {
  regulon_name <- candidate_closure$regulon[ candidate_index ]
  current_targets <- target_split[[regulon_name]]
  current_targets$progressive <- current_targets$target %in%
    progressive_genes
  current_targets$strong_progressive <- current_targets$target %in%
    strong_progressive_genes
  current_targets$progressive_direction <- unname( step5_direction_lookup[ current_targets$target ] )
  current_targets$highconf_marker <- current_targets$target %in%
    highconf_genes
  current_targets$highconf_direction <- unname( step6a_direction_lookup[ current_targets$target ] )
  current_targets$strict_pseudobulk_marker <- current_targets$target %in%
    strict_genes
  current_targets$Step6B_module <- unname( step6b_module_lookup[ current_targets$target ] )
  current_targets$formal_Arm2_module <- current_targets$Step6B_module %in%
    FORMAL_ARM2_MODULES
  current_targets$exploratory_Arm2_module <- current_targets$Step6B_module %in%
    EXPLORATORY_ARM2_MODULES
  current_targets$shared_prebranch_module <- current_targets$Step6B_module %in%
    SHARED_MODULES
  current_direction <- candidate_closure$late_direction[ candidate_index ]
  current_targets$direction_consistent_progressive <- current_targets$progressive & current_targets$progressive_direction == current_direction
  current_targets$direction_consistent_highconf <- current_targets$highconf_marker & current_targets$highconf_direction == current_direction
  current_targets$direction_consistent_module <- if ( current_direction == "Arm_2" ) {
    current_targets$formal_Arm2_module
  } else {
    current_targets$shared_prebranch_module
  }
  current_targets$evidence_count <- as.integer( current_targets$direction_consistent_progressive ) + as.integer( current_targets$strong_progressive & current_targets$progressive_direction == current_direction ) + as.integer( current_targets$direction_consistent_highconf ) + as.integer( current_targets$strict_pseudobulk_marker & current_targets$highconf_direction == current_direction ) + as.integer( current_targets$direction_consistent_module )
  current_targets$any_relevant_evidence <- current_targets$evidence_count > 0L
  current_targets$regulon_timing_classification <- candidate_closure$timing_classification[ candidate_index ]
  current_targets$regulon_late_direction <- current_direction
  current_targets$target_closure_strength <- candidate_closure$target_closure_strength[ candidate_index ]
  current_targets$mechanistic_priority <- candidate_closure$mechanistic_priority[ candidate_index ]
  overlap_gene_parts[[candidate_index]] <- current_targets[ current_targets$any_relevant_evidence, , drop = FALSE ]
}
candidate_overlap_genes <- data.table::rbindlist( overlap_gene_parts, use.names = TRUE, fill = TRUE )
candidate_overlap_genes <- as.data.frame( candidate_overlap_genes )
candidate_overlap_genes <- candidate_overlap_genes[ order( match( candidate_overlap_genes$mechanistic_priority, priority_order ), candidate_overlap_genes$regulon, -candidate_overlap_genes$evidence_count, -candidate_overlap_genes$weight ), , drop = FALSE ]
write_tsv( candidate_overlap_genes, OVERLAP_GENES_FILE )
priority_classes_for_network <- c( "Tier_1_branchpoint_coincident_target_closure", "Tier_2_branchpoint_coincident_partial_closure", "Tier_2_downstream_Arm2_regulatory_amplifier" )
priority_candidates <- candidate_closure[ candidate_closure$mechanistic_priority %in% priority_classes_for_network, , drop = FALSE ]
network_edge_parts <- list()
for ( candidate_index in seq_len( nrow(priority_candidates) ) ) {
  regulon_name <- priority_candidates$regulon[ candidate_index ]
  current_edges <- candidate_overlap_genes[ candidate_overlap_genes$regulon == regulon_name & candidate_overlap_genes$evidence_count >= MIN_CORE_EDGE_EVIDENCE_COUNT, , drop = FALSE ]
  if (nrow(current_edges) == 0L) {
    current_edges <- candidate_overlap_genes[ candidate_overlap_genes$regulon == regulon_name, , drop = FALSE ]
  }
  current_edges <- current_edges[ order( -current_edges$evidence_count, -current_edges$weight ), , drop = FALSE ]
  current_edges <- utils::head( current_edges, MAX_NETWORK_TARGETS_PER_TF )
  network_edge_parts[[ length(network_edge_parts) + 1L ]] <- current_edges
}
if (length(network_edge_parts) > 0L) {
  priority_network_edges <- data.table::rbindlist( network_edge_parts, use.names = TRUE, fill = TRUE )
  priority_network_edges <- as.data.frame( priority_network_edges )
} else {
  priority_network_edges <- candidate_overlap_genes[ 0, , drop = FALSE ]
}
write_tsv( priority_network_edges, NETWORK_EDGES_FILE )
msg("Exporting candidate target signatures...")
signature_rows <- list()
gmt_lines <- character(0)
signature_candidates <- candidate_closure[ !candidate_closure$mechanistic_priority %in% c( "Do_not_prioritize_direction_discordant", "Do_not_prioritize_sample_sensitive" ), , drop = FALSE ]
for ( candidate_index in seq_len( nrow(signature_candidates) ) ) {
  regulon_name <- signature_candidates$regulon[ candidate_index ]
  tf_name <- signature_candidates$TF[ candidate_index ]
  current_all_targets <- unique( target_split[[regulon_name]]$target )
  current_overlap <- candidate_overlap_genes[ candidate_overlap_genes$regulon == regulon_name, , drop = FALSE ]
  direction_consistent_targets <- unique( current_overlap$target[ current_overlap$direction_consistent_progressive | current_overlap$direction_consistent_highconf | current_overlap$direction_consistent_module ] )
  high_confidence_closure_targets <- unique( current_overlap$target[ current_overlap$evidence_count >= MIN_CORE_EDGE_EVIDENCE_COUNT ] )
  signature_list <- list( all_motif_supported_targets = current_all_targets, direction_consistent_targets = direction_consistent_targets, high_confidence_closure_targets = high_confidence_closure_targets )
  for (signature_type in names(signature_list)) {
    signature_genes <- unique( signature_list[[signature_type]] )
    signature_name <- paste( "GSE162610", "Step7G", sanitize_signature_name(tf_name), signature_type, sep = "__" )
    signature_rows[[ length(signature_rows) + 1L ]] <- data.frame( signature_name = signature_name, TF = tf_name, regulon = regulon_name, signature_type = signature_type, n_genes = length( signature_genes ), mechanistic_priority = signature_candidates$mechanistic_priority[ candidate_index ], genes = paste( sort( signature_genes ), collapse = " | " ), stringsAsFactors = FALSE )
    if (length(signature_genes) > 0L) {
      gmt_lines <- c( gmt_lines, paste( c( signature_name, paste0( "GSE162610 Step7G ", tf_name, " ", signature_type ), sort( signature_genes ) ), collapse = "\t" ) )
    }
  }
}
signature_summary <- do.call( rbind, signature_rows )
write_tsv( signature_summary, SIGNATURE_TABLE_FILE )
writeLines( gmt_lines, con = SIGNATURE_GMT_FILE, useBytes = TRUE )
msg("Checking optional local GO BP packages...")
go_status <- data.frame( packages_available = FALSE, candidate_gene_sets_tested = 0L, significant_rows = 0L, status = "SKIPPED", note = "", stringsAsFactors = FALSE )
go_results <- data.frame( note = character(0), stringsAsFactors = FALSE )
go_packages_available <- requireNamespace( "clusterProfiler", quietly = TRUE ) && requireNamespace( "org.Mm.eg.db", quietly = TRUE )
if (go_packages_available) {
  go_status$packages_available <- TRUE
  orgdb <- get( "org.Mm.eg.db", envir = asNamespace( "org.Mm.eg.db" ) )
  go_parts <- list()
  go_candidates <- signature_summary[ signature_summary$signature_type == "high_confidence_closure_targets" & signature_summary$n_genes >= 5L, , drop = FALSE ]
  for ( candidate_index in seq_len( nrow(go_candidates) ) ) {
    candidate_genes <- strsplit( go_candidates$genes[ candidate_index ], " | ", fixed = TRUE )[[1]]
    candidate_genes <- candidate_genes[ nzchar(candidate_genes) ]
    go_object <- tryCatch( clusterProfiler::enrichGO( gene = candidate_genes, universe = universe_genes, OrgDb = orgdb, keyType = "SYMBOL", ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1, readable = FALSE ), error = function(e) {
        NULL
      }
    )
    if (!is.null(go_object)) {
      go_table <- as.data.frame( go_object )
      if (nrow(go_table) > 0L) {
        go_table$TF <- go_candidates$TF[ candidate_index ]
        go_table$regulon <- go_candidates$regulon[ candidate_index ]
        go_table$signature_name <- go_candidates$signature_name[ candidate_index ]
        go_parts[[ length(go_parts) + 1L ]] <- go_table
      }
    }
  }
  if (length(go_parts) > 0L) {
    go_results <- data.table::rbindlist( go_parts, use.names = TRUE, fill = TRUE )
    go_results <- as.data.frame( go_results )
    go_results <- go_results[ order( go_results$p.adjust, -go_results$Count ), , drop = FALSE ]
  }
  go_status$candidate_gene_sets_tested <- nrow( go_candidates )
  go_status$significant_rows <- if ( nrow(go_results) > 0L ) {
    sum( go_results$p.adjust <= 0.05, na.rm = TRUE )
  } else {
    0L
  }
  go_status$status <- "COMPLETED"
  go_status$note <- "Optional GO BP used local org.Mm.eg.db; it does not determine Step7G PASS."
} else {
  go_status$note <- "clusterProfiler and/or org.Mm.eg.db unavailable; overlap closure remains complete."
}
write_tsv( go_results, GO_BP_FILE )
write_csv( go_status, GO_STATUS_FILE )
heatmap_evidence_sets <- c( "Progressive_Arm2", "Strong_progressive_Arm2", "Step6A_Arm2_highconf", "Step6A_Arm2_strict", "Step6B_M1_Arm2_terminal", "Step6B_M2_Arm2_terminal", "Step6B_formal_Arm2_M1_M2", "Step6B_M3_Shared_prebranch" )
heatmap_candidates <- candidate_closure[ !candidate_closure$mechanistic_priority %in% c( "Do_not_prioritize_sample_sensitive" ), , drop = FALSE ]
heatmap_candidates <- utils::head( heatmap_candidates, MAX_HEATMAP_CANDIDATES )
evidence_matrix <- matrix( 0, nrow = nrow( heatmap_candidates ), ncol = length( heatmap_evidence_sets ), dimnames = list( paste0( heatmap_candidates$TF, " [", heatmap_candidates$late_direction, "]" ), heatmap_evidence_sets ) )
overlap_count_matrix <- matrix( 0, nrow = nrow( heatmap_candidates ), ncol = length( heatmap_evidence_sets ), dimnames = dimnames( evidence_matrix ) )
for ( row_index in seq_len( nrow(heatmap_candidates) ) ) {
  regulon_name <- heatmap_candidates$regulon[ row_index ]
  for ( column_index in seq_along( heatmap_evidence_sets ) ) {
    evidence_name <- heatmap_evidence_sets[ column_index ]
    enrichment_row <- get_enrichment_row( all_enrichment, regulon_name, evidence_name )
    q_value <- safe_value( enrichment_row, "q_within_evidence_set", 1 )
    fold_value <- safe_value( enrichment_row, "fold_enrichment", 1 )
    overlap_value <- safe_value( enrichment_row, "overlap_n", 0 )
    signed_score <- if ( is.finite(q_value) && is.finite(fold_value) ) {
      min( -log10( max( q_value, 1e-300 ) ), 20 ) * sign( log2( max( fold_value, 1e-12 ) ) )
    } else {
      0
    }
    evidence_matrix[ row_index, column_index ] <- signed_score
    overlap_count_matrix[ row_index, column_index ] <- overlap_value
  }
}
evidence_matrix_table <- data.frame( candidate = rownames( evidence_matrix ), evidence_matrix, check.names = FALSE, stringsAsFactors = FALSE )
write_tsv( evidence_matrix_table, EVIDENCE_MATRIX_FILE )
msg("Generating Step7G figures...")
heatmap_limit <- max( abs( evidence_matrix ), na.rm = TRUE )
if ( !is.finite( heatmap_limit ) || heatmap_limit <= 0 ) {
  heatmap_limit <- 1
}
heatmap_breaks <- seq( -heatmap_limit, heatmap_limit, length.out = 101L )
heatmap_palette <- grDevices::colorRampPalette( c( "#3B4CC0", "white", "#B40426" ) )( 100L )
grDevices::pdf( FIG_ENRICHMENT_HEATMAP, width = 12, height = max( 7, nrow( evidence_matrix ) * 0.24 ) )
graphics::par( mar = c( 10, 11, 4, 2 ) )
graphics::image( x = seq_len( ncol( evidence_matrix ) ), y = seq_len( nrow( evidence_matrix ) ), z = t( evidence_matrix[ nrow(evidence_matrix):1, , drop = FALSE ] ), col = heatmap_palette, breaks = heatmap_breaks, xlab = "", ylab = "", xaxt = "n", yaxt = "n", main = "Signed -log10(BH q) target enrichment\nred=enriched; blue=depleted" )
graphics::axis( side = 1, at = seq_len( ncol( evidence_matrix ) ), labels = colnames( evidence_matrix ), las = 2, cex.axis = 0.68 )
graphics::axis( side = 2, at = seq_len( nrow( evidence_matrix ) ), labels = rev( rownames( evidence_matrix ) ), las = 2, cex.axis = 0.68 )
grDevices::dev.off()
score_plot_table <- candidate_closure[ !candidate_closure$mechanistic_priority %in% c( "Do_not_prioritize_sample_sensitive" ), , drop = FALSE ]
score_plot_table <- utils::head( score_plot_table, 40L )
score_colors <- ifelse( score_plot_table$late_direction == "Arm_2", "#B64B45", "#2F6B9A" )
grDevices::pdf( FIG_EVIDENCE_SCORE, width = 10, height = max( 6, nrow( score_plot_table ) * 0.22 ) )
graphics::par( mar = c( 5, 10, 3, 2 ) )
bar_positions <- graphics::barplot( height = score_plot_table$target_closure_domains, names.arg = score_plot_table$TF, horiz = TRUE, las = 1, col = score_colors, xlim = c( 0, 4.5 ), xlab = "Supported target-closure domains (0-4)", main = "Candidate target-program closure score" )
graphics::abline( v = c( 2, 3 ), lty = c( 3, 2 ), col = "grey40" )
graphics::text( x = score_plot_table$target_closure_domains + 0.08, y = bar_positions, labels = score_plot_table$target_closure_strength, pos = 4, cex = 0.65 )
grDevices::dev.off()
count_plot_table <- heatmap_candidates
count_matrix_plot <- cbind( Progressive = count_plot_table$progressive_overlap_n, Strong = count_plot_table$strong_progressive_overlap_n, Highconf = count_plot_table$highconf_marker_overlap_n, Formal_module = count_plot_table$formal_module_overlap_n )
rownames(count_matrix_plot) <- count_plot_table$TF
grDevices::pdf( FIG_OVERLAP_COUNTS, width = 11, height = max( 7, nrow( count_matrix_plot ) * 0.25 ) )
graphics::par( mar = c( 5, 10, 3, 2 ) )
graphics::barplot( t( count_matrix_plot ), beside = FALSE, horiz = TRUE, las = 1, names.arg = rownames( count_matrix_plot ), col = c( "#4C78A8", "#72B7B2", "#F58518", "#E45756" ), xlab = "Direction-consistent target overlaps", main = "Candidate target overlaps across evidence layers" )
graphics::legend( "bottomright", legend = colnames( count_matrix_plot ), fill = c( "#4C78A8", "#72B7B2", "#F58518", "#E45756" ), bty = "n", cex = 0.75 )
grDevices::dev.off()
grDevices::pdf( FIG_TARGET_EDGES, width = 11, height = 8.5, onefile = TRUE )
if (nrow(priority_network_edges) == 0L) {
  graphics::plot.new()
  graphics::text( 0.5, 0.5, "No priority candidate produced target edges under the predefined criteria." )
} else {
  priority_tf_order <- unique( priority_network_edges$TF )
  for (tf_name in priority_tf_order) {
    current_edges <- priority_network_edges[ priority_network_edges$TF == tf_name, , drop = FALSE ]
    current_edges <- current_edges[ order( current_edges$evidence_count, current_edges$weight ), , drop = FALSE ]
    edge_colors <- grDevices::colorRampPalette( c( "#D9E8F5", "#F4A582", "#B2182B" ) )( 6L )
    current_color_index <- pmin( pmax( current_edges$evidence_count + 1L, 1L ), 6L )
    graphics::par( mar = c( 5, 10, 4, 2 ) )
    bar_positions <- graphics::barplot( current_edges$weight, names.arg = current_edges$target, horiz = TRUE, las = 1, col = edge_colors[ current_color_index ], xlab = "pySCENIC target weight", main = paste0( tf_name, " motif-supported target closure" ) )
    graphics::text( current_edges$weight, bar_positions, labels = paste0( "evidence=", current_edges$evidence_count ), pos = 4, cex = 0.62 )
  }
}
grDevices::dev.off()
technical_enrichment_pass <- nrow( all_enrichment ) == EXPECTED_PRIMARY_REGULONS * length( evidence_sets ) && all( is.finite( all_enrichment$p_value ) ) && all( all_enrichment$p_value >= 0 & all_enrichment$p_value <= 1 ) && all( is.finite( all_enrichment$q_within_evidence_set ) )
candidate_annotation_pass <- nrow( candidate_closure ) > 0L && all( !is.na( candidate_closure$target_closure_strength ) ) && all( !is.na( candidate_closure$mechanistic_priority ) )
signature_pass <- file.exists( SIGNATURE_GMT_FILE ) && file.info( SIGNATURE_GMT_FILE )$size >
    0 && nrow( signature_summary ) >
    0L
figure_pass <- all( file.exists( c( FIG_ENRICHMENT_HEATMAP, FIG_EVIDENCE_SCORE, FIG_OVERLAP_COUNTS, FIG_TARGET_EDGES ) ) )
final_decision <- data.frame( domain = c( "Required input files", "Step 7C Overall PASS", "Step 7F Overall PASS", "Frozen 7897-gene universe", "Step7C 292 regulons and 29609 target edges", "Step5C-R 212 progressive genes", "Step5C-R 172 strong progressive genes", "Step6A 745 high-confidence markers", "Step6A 343 strict pseudobulk markers", "Step6B 803 genes in five dynamic modules", "All 215 primary regulons tested against every evidence set", "Hypergeometric p and BH q values valid", "Timing-candidate target closure integrated", "Candidate overlap-gene table produced", "Candidate target signatures produced", "Priority target-edge table produced", "Step7G figures produced", "Complete Step7G RDS produced", "Overall Step 7G" ), pass = c( all( input_audit$exists ), step7c_pass, step7f_pass, length( universe_genes ) == EXPECTED_UNIVERSE_GENES, nrow( regulon_targets ) == EXPECTED_REGULON_TARGET_EDGES && length( unique( regulon_targets$regulon ) ) == EXPECTED_REGULONS, length( progressive_genes ) == EXPECTED_PROGRESSIVE_GENES, length( strong_progressive_genes ) == EXPECTED_STRONG_PROGRESSIVE_GENES, length( highconf_genes ) == EXPECTED_HIGHCONF_MARKERS, length( strict_genes ) == EXPECTED_STRICT_PB_MARKERS, length( unique( step6b_assignments$gene_id ) ) == EXPECTED_STEP6B_MODULE_GENES && length( unique( step6b_assignments$module_name ) ) == EXPECTED_STEP6B_MODULES, nrow( all_enrichment ) == EXPECTED_PRIMARY_REGULONS * length( evidence_sets ), technical_enrichment_pass, candidate_annotation_pass, file.exists( OVERLAP_GENES_FILE ), signature_pass, file.exists( NETWORK_EDGES_FILE ), figure_pass, FALSE, FALSE ), stringsAsFactors = FALSE )
final_decision$status <- ifelse( final_decision$pass, "PASS", "FAIL" )
step7g_results <- list( step = "Step7G", generated_at = timestamp(), parameters = list( universe_n = EXPECTED_UNIVERSE_GENES, primary_q_max = PRIMARY_Q_MAX, sensitivity_q_max = SENSITIVITY_Q_MAX, minimum_fold_enrichment = MIN_FOLD_ENRICHMENT, minimum_core_edge_evidence_count = MIN_CORE_EDGE_EVIDENCE_COUNT ), universe_genes = universe_genes, evidence_metadata = evidence_metadata, evidence_sets = evidence_sets, gene_set_audit = gene_set_audit, regulon_targets = regulon_targets, timing_master = timing_master, all_enrichment = all_enrichment, all_regulon_closure_summary = all_closure, candidate_closure = candidate_closure, candidate_overlap_genes = candidate_overlap_genes, priority_network_edges = priority_network_edges, signature_summary = signature_summary, optional_GO_results = go_results, optional_GO_status = go_status, final_decision = final_decision )
saveRDS( step7g_results, RESULT_RDS_FILE, compress = TRUE )
rds_pass <- file.exists( RESULT_RDS_FILE )
final_decision$pass[ final_decision$domain == "Complete Step7G RDS produced" ] <- rds_pass
non_overall_pass <- final_decision$pass[ final_decision$domain != "Overall Step 7G" ]
final_decision$pass[ final_decision$domain == "Overall Step 7G" ] <- all( non_overall_pass )
final_decision$status <- ifelse( final_decision$pass, "PASS", "FAIL" )
write_csv( final_decision, FINAL_DECISION_FILE )
FINAL_STATUS <- if ( isTRUE( final_decision$pass[ final_decision$domain == "Overall Step 7G" ][[1]] ) ) {
  "PASS"
} else {
  "FAIL"
}
step7g_results$status <- FINAL_STATUS
step7g_results$final_decision <- final_decision
saveRDS( step7g_results, RESULT_RDS_FILE, compress = TRUE )
overall_elapsed <- as.numeric( difftime( Sys.time(), overall_start, units = "secs" ) )
closure_counts <- as.data.frame( table( candidate_closure$target_closure_strength, useNA = "ifany" ), stringsAsFactors = FALSE )
colnames(closure_counts) <- c( "target_closure_strength", "n_candidates" )
priority_counts <- as.data.frame( table( candidate_closure$mechanistic_priority, useNA = "ifany" ), stringsAsFactors = FALSE )
colnames(priority_counts) <- c( "mechanistic_priority", "n_candidates" )
tier1_candidates <- candidate_closure[ candidate_closure$mechanistic_priority == "Tier_1_branchpoint_coincident_target_closure", , drop = FALSE ]
tier2_early_candidates <- candidate_closure[ candidate_closure$mechanistic_priority == "Tier_2_branchpoint_coincident_partial_closure", , drop = FALSE ]
downstream_amplifiers <- candidate_closure[ candidate_closure$mechanistic_priority == "Tier_2_downstream_Arm2_regulatory_amplifier", , drop = FALSE ]
internal_fate_driver_candidates <- candidate_closure[ candidate_closure$internal_fate_driver_criteria_met, , drop = FALSE ]
report_candidate_columns <- c( "TF", "regulon", "timing_classification", "late_direction", "onset_progress_qualified", "terminal_effect_z", "earlier_than_terminal_program", "progressive_overlap_n", "progressive_fold_enrichment", "progressive_q", "strong_progressive_overlap_n", "strong_progressive_fold_enrichment", "strong_progressive_q", "highconf_marker_overlap_n", "highconf_marker_fold_enrichment", "highconf_marker_q", "formal_module_overlap_n", "formal_module_fold_enrichment", "formal_module_q", "strict_marker_overlap_n", "target_closure_domains", "target_closure_strength", "seven_dpi_direction_consistent", "CSI_module", "mechanistic_priority" )
report_lines <- c( "======================================================================", "GSE162610 | STEP 7G REPORT", "Regulon target-program closure", "======================================================================", "", paste0( "Generated: ", timestamp() ), paste0( "Project root: ", PROJECT_ROOT ), paste0( "Overall Step 7G status: ", FINAL_STATUS ), "", "FORMAL SCOPE:", paste0( "Step7G tested motif-supported regulon targets against Step5C-R ", "progressive genes, Step6A branch markers and Step6B dynamic modules." ), paste0( "Hypergeometric tests used the frozen Step7A 7,897-gene universe. ", "BH correction was calculated globally and within each evidence-set family." ), "", "FROZEN INPUT:", paste0( "Gene universe: ", length( universe_genes ) ), paste0( "Regulons: ", length( unique( regulon_targets$regulon ) ), "; target edges: ", nrow( regulon_targets ) ), paste0( "Primary regulons tested: ", nrow( timing_master ) ), paste0( "Progressive genes: ", length( progressive_genes ), "; strong progressive genes: ", length( strong_progressive_genes ) ), paste0( "Step6A high-confidence markers: ", length( highconf_genes ), "; strict pseudobulk markers: ", length( strict_genes ) ), paste0( "Step6B dynamic genes: ", length( all_dynamic_module_genes ), " across ", length( step6b_module_sets ), " modules." ), "", "CANDIDATE CLOSURE DEFINITIONS:", paste0( "A supported domain required within-evidence-set BH q<=", PRIMARY_Q_MAX, ", fold enrichment>=", MIN_FOLD_ENRICHMENT, " and the predefined minimum overlap count." ), paste0( "Strong closure required at least three of four domains and evidence ", "from both a dynamic-gene domain and a terminal marker/module domain." ), paste0( "The four domains were direction-consistent progressive genes, strong ", "progressive genes, Step6A high-confidence markers and the corresponding ", "Step6B module reference." ), "", "TARGET-CLOSURE STRENGTH COUNTS:", paste( capture.output( print( closure_counts, row.names = FALSE ) ), collapse = "\n" ), "", "MECHANISTIC PRIORITY COUNTS:", paste( capture.output( print( priority_counts, row.names = FALSE ) ), collapse = "\n" ), "", "TIER 1 | BRANCHPOINT-COINCIDENT TARGET CLOSURE:", if ( nrow( tier1_candidates ) > 0L ) {
    paste( capture.output( print( tier1_candidates[ , report_candidate_columns, drop = FALSE ], row.names = FALSE ) ), collapse = "\n" )
  } else {
    "None met the predefined Tier-1 target-closure criteria."
  },
  "", "TIER 2 | BRANCHPOINT-COINCIDENT PARTIAL CLOSURE:", if ( nrow( tier2_early_candidates ) > 0L ) {
    paste( capture.output( print( tier2_early_candidates[ , report_candidate_columns, drop = FALSE ], row.names = FALSE ) ), collapse = "\n" )
  } else {
    "None met the predefined Tier-2 early partial-closure criteria."
  },
  "", "TIER 2 | DOWNSTREAM ARM_2 REGULATORY AMPLIFIERS:", if ( nrow( downstream_amplifiers ) > 0L ) {
    paste( capture.output( print( downstream_amplifiers[ , report_candidate_columns, drop = FALSE ], row.names = FALSE ) ), collapse = "\n" )
  } else {
    "None met the downstream-amplifier closure criteria."
  },
  "", "INTERNAL FATE-DRIVER CRITERIA:", paste0( "Candidates earlier than the formal terminal-program reference and with ", "strong target closure: ", nrow( internal_fate_driver_candidates ) ), if ( nrow( internal_fate_driver_candidates ) == 0L ) {
    paste0( "No regulon currently satisfies the internal pre-terminal fate-driver ", "criterion. Candidates with onset=0 remain branchpoint-coincident." )
  } else {
    paste( capture.output( print( internal_fate_driver_candidates[ , report_candidate_columns, drop = FALSE ], row.names = FALSE ) ), collapse = "\n" )
  },
  "", "OPTIONAL GO BP:", paste( capture.output( print( go_status, row.names = FALSE ) ), collapse = "\n" ), "", "FINAL DECISION:", paste( capture.output( print( final_decision, row.names = FALSE ) ), collapse = "\n" ), "", "INTERPRETATION BOUNDARY:", paste0( "Target-program closure supports a regulatory relationship, not causal ", "lineage control. Because the earliest formal terminal programs and the ", "earliest regulon candidates were both detected at Y_29=0, no candidate ", "should be described as demonstrably preceding terminal-program activation." ), paste0( "Tier-1 candidates may be called branchpoint-coincident regulatory ", "candidates with direction-consistent target closure. A fate-driver claim ", "still requires external temporal/spatial replication and functional testing." ), "", "NEXT ACTION:", if ( FINAL_STATUS == "PASS" ) {
    paste0( "Proceed to Step8 CellChat using the frozen neutral branch labels. ", "Candidate TF-target signatures exported here may be connected to ", "branch-specific ligands, but CellChat does not prove transcriptional causality." )
  } else {
    paste0( "Do not proceed to Step8 until the failed Step7G audit domain is resolved." )
  },
  "", "PRIMARY OUTPUTS:", paste0( "All enrichment tests: ", safe_path( ALL_ENRICHMENT_FILE ) ), paste0( "Timing-candidate closure: ", safe_path( CANDIDATE_CLOSURE_FILE ) ), paste0( "Candidate overlap genes: ", safe_path( OVERLAP_GENES_FILE ) ), paste0( "Candidate target signatures: ", safe_path( SIGNATURE_GMT_FILE ) ), paste0( "Complete Step7G RDS: ", safe_path( RESULT_RDS_FILE ) ), "", "KEY OUTPUTS:", "1) reports/Step7G_report.txt", "2) tables/03_timing_candidate_target_closure.tsv", "3) tables/04_candidate_overlap_genes_long.tsv", "4) figures/Fig7G_01_candidate_enrichment_heatmap.pdf", "5) tables/08_final_decision.csv"
)
writeLines( report_lines, con = REPORT_FILE, useBytes = TRUE )
writeLines( capture.output( sessionInfo() ), con = SESSION_FILE, useBytes = TRUE )
gc()
msg("==============================================================")
msg("STEP 7G COMPLETE")
msg("==============================================================")
msg("Step 7G status: ", FINAL_STATUS)
msg( "All enrichment tests: ", nrow( all_enrichment ) )
msg( "Tier-1 branchpoint-coincident closure candidates: ", nrow( tier1_candidates ) )
msg( "Internal pre-terminal fate-driver criteria met: ", nrow( internal_fate_driver_candidates ) )
msg( "Total time: ", format_duration( overall_elapsed ) )
msg( "Report: ", safe_path( REPORT_FILE ) )
msg("")
msg("Key outputs:")
msg("1) reports/Step7G_report.txt")
msg("2) tables/03_timing_candidate_target_closure.tsv")
msg("3) tables/04_candidate_overlap_genes_long.tsv")
msg("4) figures/Fig7G_01_candidate_enrichment_heatmap.pdf")
msg("5) tables/08_final_decision.csv")
