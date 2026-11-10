rm(list = ls())
gc()
options( stringsAsFactors = FALSE, width = 220, warn = 1, timeout = 86400 )
set.seed(20260801)
PROJECT_ROOT <- Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
CONDA_ENV_NAME <- "pyscenic_env"
NUM_WORKERS <- 8L
SMOKE_WORKERS <- 2L
RUN_SMOKE_TEST <- TRUE
FORCE_RERUN_CTX <- FALSE
MIN_GENES <- 10L
RANK_THRESHOLD <- 5000L
AUC_THRESHOLD <- 0.05
NES_THRESHOLD <- 3.0
CHUNK_SIZE <- 100L
CTX_MODE <- "custom_multiprocessing"
MASK_DROPOUTS <- FALSE
ALL_MODULES <- FALSE
EXPECTED_CELLS <- 4203L
EXPECTED_GENES <- 7897L
EXPECTED_GRN_TFS <- 1023L
EXPECTED_GRN_TARGETS <- 7897L
EXPECTED_GRN_EDGES <- 653060L
EXPECTED_STATE_COUNTS <- c( Shared_prebranch = 116L, Arm_1 = 2289L, Arm_2 = 1798L )
MIN_MOTIF_ROWS <- 20L
MIN_REGULONS <- 20L
MIN_REGULON_TFS <- 20L
MIN_KEY_TFS_RETAINED_FOR_AUDIT <- 1L
KEY_TFS <- unique( c( "Spi1", "Junb", "Cebpb", "Egr1", "Fos", "Fosb", "Jun", "Atf3", "Mafb", "Nfkb1", "Nfkb2", "Rela", "Relb", "Stat1", "Stat3", "Irf1", "Irf7", "Tfeb", "Mitf", "Nr1h3", "Srebf1", "Pparg" ) )
SMOKE_TFS_PREFERRED <- c( "Spi1", "Junb", "Cebpb", "Egr1", "Fos", "Mafb", "Nfkb1", "Stat1" )
SMOKE_EDGES_PER_TF <- 300L
STEP7B_DIR <- file.path( PROJECT_ROOT, "04_results", "step7B_GSE162610_GRNBoost2_network_inference" )
ADJ_FILE <- file.path( STEP7B_DIR, "network", "GSE162610_step7B_GRNBoost2_adjacencies.tsv" )
EXPRESSION_LOOM <- file.path( STEP7B_DIR, "input", "GSE162610_step7B_frozen_4203_raw_counts.loom" )
STEP7B_REPORT <- file.path( STEP7B_DIR, "reports", "Step7B_report.txt" )
RESOURCE_DIR <- file.path( PROJECT_ROOT, "01_raw_data", "pySCENIC_resources" )
TF_FILE <- file.path( RESOURCE_DIR, "mm_mgi_tfs.txt" )
MOTIF_ANNOTATION_FILE <- file.path( RESOURCE_DIR, "motifs-v9-nr.mgi-m0.001-o0.0.tbl" )
FEATHER_FILE <- file.path( RESOURCE_DIR, paste0( "mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.", "genes_vs_motifs.rankings.feather" ) )
OUTDIR <- file.path( PROJECT_ROOT, "04_results", "step7C_GSE162610_cisTarget_motif_pruning_regulon_inference" )
INPUTDIR <- file.path(OUTDIR, "input")
REGULONDIR <- file.path(OUTDIR, "regulons")
TABDIR <- file.path(OUTDIR, "tables")
LOGDIR <- file.path(OUTDIR, "logs")
REPORTDIR <- file.path(OUTDIR, "reports")
PYDIR <- file.path(OUTDIR, "python")
TMPDIR_CTX <- file.path(OUTDIR, "tmp")
for (directory in c( OUTDIR, INPUTDIR, REGULONDIR, TABDIR, LOGDIR, REPORTDIR, PYDIR, TMPDIR_CTX )) {
  dir.create( directory, recursive = TRUE, showWarnings = FALSE )
}
SMOKE_ADJ_FILE <- file.path( INPUTDIR, "GSE162610_step7C_smoke_adjacencies.tsv" )
SMOKE_MOTIF_FILE <- file.path( REGULONDIR, "GSE162610_step7C_smoke_enriched_motifs.csv" )
MOTIF_OUTPUT_FILE <- file.path( REGULONDIR, "GSE162610_step7C_enriched_motifs.csv" )
REGULON_DAT_FILE <- file.path( REGULONDIR, "GSE162610_step7C_regulons.dat" )
REGULON_GMT_FILE <- file.path( REGULONDIR, "GSE162610_step7C_regulons.gmt" )
REGULON_JSON_FILE <- file.path( REGULONDIR, "GSE162610_step7C_regulons.json" )
PREFLIGHT_FILE <- file.path( TABDIR, "01_step7C_preflight_audit.tsv" )
AUDIT_METRICS_FILE <- file.path( TABDIR, "03_regulon_pruning_summary.tsv" )
REGULON_SUMMARY_FILE <- file.path( TABDIR, "04_regulon_summary.tsv" )
REGULON_TARGET_FILE <- file.path( TABDIR, "05_regulon_targets_long.tsv" )
KEY_TF_AUDIT_FILE <- file.path( TABDIR, "06_key_TF_regulon_audit.tsv" )
MOTIF_TF_SUMMARY_FILE <- file.path( TABDIR, "07_TF_motif_enrichment_summary.tsv" )
FINAL_DECISION_FILE <- file.path( TABDIR, "08_final_decision.csv" )
REPORT_FILE <- file.path( REPORTDIR, "Step7C_report.txt" )
ERROR_FILE <- file.path( REPORTDIR, "Step7C_ERROR.txt" )
SESSION_FILE <- file.path( REPORTDIR, "Step7C_sessionInfo.txt" )
if (!requireNamespace("data.table", quietly = TRUE)) {
  stop( "The R package 'data.table' is required.", call. = FALSE )
}
timestamp <- function() {
  format( Sys.time(), "%Y-%m-%d %H:%M:%S" )
}
safe_path <- function(path) {
  normalizePath( path, winslash = "/", mustWork = FALSE )
}
quote_cmd <- function(x) {
  shQuote( as.character(x), type = "cmd" )
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
format_bytes <- function(x) {
  x <- as.numeric(x)
  if (!is.finite(x)) {
    return(NA_character_)
  }
  units <- c("B", "KB", "MB", "GB", "TB")
  idx <- 1L
  while (x >= 1024 && idx < length(units)) {
    x <- x / 1024
    idx <- idx + 1L
  }
  paste0( format( round(x, 3), nsmall = 3, trim = TRUE ), " ", units[[idx]] )
}
write_utf8_csv <- function(x, path) {
  utils::write.csv( x, file = path, row.names = FALSE, fileEncoding = "UTF-8", na = "" )
}
write_utf8_tsv <- function(x, path) {
  data.table::fwrite( x, file = path, sep = "\t", quote = FALSE, na = "" )
}
read_log_tail <- function(path, n = 150L) {
  if (!file.exists(path)) {
    return("<log not found>")
  }
  lines <- readLines( path, warn = FALSE, encoding = "UTF-8" )
  if (length(lines) == 0L) {
    return("<empty log>")
  }
  paste( utils::tail(lines, n), collapse = "\n" )
}
abort_step <- function(text) {
  writeLines( c( "======================================================================", "GSE162610 | STEP 7C ERROR", "======================================================================", paste0("Time: ", timestamp()), "", text ), con = ERROR_FILE, useBytes = TRUE )
  stop( text, call. = FALSE )
}
with_ctx_environment <- function(code) {
  variable_names <- c( "TEMP", "TMP", "TMPDIR", "JOBLIB_TEMP_FOLDER", "PYTHONUTF8", "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS" )
  old_values <- Sys.getenv( variable_names, unset = NA_character_ )
  on.exit( {
      for (index in seq_along(variable_names)) {
        variable_name <- variable_names[[index]]
        old_value <- old_values[[index]]
        if (is.na(old_value)) {
          Sys.unsetenv(variable_name)
        } else {
          do.call( Sys.setenv, setNames( list(old_value), variable_name ) )
        }
      }
    },
    add = TRUE
  )
  Sys.setenv( TEMP = safe_path(TMPDIR_CTX), TMP = safe_path(TMPDIR_CTX), TMPDIR = safe_path(TMPDIR_CTX), JOBLIB_TEMP_FOLDER = safe_path(TMPDIR_CTX), PYTHONUTF8 = "1", OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1", NUMEXPR_NUM_THREADS = "1" )
  force(code)
}
run_command <- function( command, args = character(0), log_file ) {
  if (file.exists(log_file)) {
    file.remove(log_file)
  }
  start_time <- Sys.time()
  status <- with_ctx_environment( tryCatch( system2( command = command, args = args, stdout = log_file, stderr = log_file ), error = function(e) {
        999L
      }
    )
  )
  if (is.null(status)) {
    status <- 0L
  }
  elapsed_seconds <- as.numeric( difftime( Sys.time(), start_time, units = "secs" ) )
  list( status = as.integer(status), elapsed_seconds = elapsed_seconds, elapsed_human = format_duration(elapsed_seconds), log_file = log_file )
}
capture_command <- function( command, args = character(0) ) {
  output <- with_ctx_environment( tryCatch( suppressWarnings( system2( command = command, args = args, stdout = TRUE, stderr = TRUE ) ), error = function(e) {
        structure( paste0( "SYSTEM2_ERROR: ", conditionMessage(e) ), status = 999L )
      }
    )
  )
  status <- attr(output, "status")
  if (is.null(status)) {
    status <- 0L
  }
  list( status = as.integer(status), output = as.character(output) )
}
find_environment_python <- function() {
  candidates <- unique( c( file.path( Sys.getenv("USERPROFILE"), "anaconda3", "envs", CONDA_ENV_NAME, "python.exe" ), file.path( Sys.getenv("USERPROFILE"), "miniconda3", "envs", CONDA_ENV_NAME, "python.exe" ), "C:/ProgramData/anaconda3/envs/pyscenic_env/python.exe", "C:/ProgramData/miniconda3/envs/pyscenic_env/python.exe" ) )
  candidates <- candidates[ file.exists(candidates) ]
  if (length(candidates) == 0L) {
    return(NA_character_)
  }
  safe_path(candidates[[1]])
}
valid_motif_output <- function(path) {
  if ( !file.exists(path) || is.na(file.info(path)$size) || file.info(path)$size < 1024 ) {
    return(FALSE)
  }
  first_lines <- tryCatch( readLines( path, n = 3L, warn = FALSE, encoding = "UTF-8" ), error = function(e) {
      character(0)
    }
  )
  length(first_lines) >= 2L && any( grepl( "Enrichment|MotifID|TargetGenes", first_lines, ignore.case = TRUE ) )
}
read_metric_tsv <- function(path) {
  table <- data.table::fread( path, sep = "\t", header = TRUE, data.table = FALSE, showProgress = FALSE )
  if (!all(c("metric", "value") %in% colnames(table))) {
    stop( "Metric table lacks metric/value columns.", call. = FALSE )
  }
  setNames( as.character(table$value), as.character(table$metric) )
}
overall_start <- Sys.time()
required_files <- c( ADJ_FILE, EXPRESSION_LOOM, TF_FILE, MOTIF_ANNOTATION_FILE, FEATHER_FILE )
missing_files <- required_files[ !file.exists(required_files) ]
if (length(missing_files) > 0L) {
  abort_step( paste0( "Required Step 7C files are missing:\n", paste( safe_path(missing_files), collapse = "\n" ) ) )
}
input_file_audit <- data.frame( role = c( "GRNBoost2_adjacencies", "expression_loom", "mouse_TF_list", "mouse_motif_annotation", "mm10_10kb_feather" ), path = safe_path(required_files), size_bytes = file.info(required_files)$size, size_human = vapply( file.info(required_files)$size, format_bytes, FUN.VALUE = character(1) ), exists = file.exists(required_files), stringsAsFactors = FALSE )
write_utf8_csv( input_file_audit, file.path( TABDIR, "00_input_file_audit.csv" ) )
ENV_PYTHON <- find_environment_python()
if ( is.na(ENV_PYTHON) || !file.exists(ENV_PYTHON) ) {
  abort_step( "The pyscenic_env Python executable was not found." )
}
PYSCENIC_EXE <- file.path( dirname(ENV_PYTHON), "Scripts", "pyscenic.exe" )
if (!file.exists(PYSCENIC_EXE)) {
  abort_step( paste0( "pySCENIC CLI was not found at: ", PYSCENIC_EXE ) )
}
PYSCENIC_EXE <- safe_path(PYSCENIC_EXE)
ctx_help <- capture_command( command = PYSCENIC_EXE, args = c( "ctx", "--help" ) )
writeLines( ctx_help$output, con = file.path( LOGDIR, "00_pyscenic_ctx_help.txt" ), useBytes = TRUE )
required_help_tokens <- c( "--annotations_fname", "--expression_mtx_fname", "--mode", "--min_genes", "--num_workers" )
ctx_cli_pass <- ctx_help$status == 0L && all( vapply( required_help_tokens, function(token) {
        any( grepl( token, ctx_help$output, fixed = TRUE ) )
      },
      FUN.VALUE = logical(1)
    )
  )
if (!ctx_cli_pass) {
  abort_step( paste0( "pySCENIC ctx CLI precheck failed.\n\n", paste( ctx_help$output, collapse = "\n" ) ) )
}
adjacency <- data.table::fread( ADJ_FILE, showProgress = TRUE )
required_adjacency_columns <- c( "TF", "target", "importance" )
missing_adjacency_columns <- setdiff( required_adjacency_columns, colnames(adjacency) )
if (length(missing_adjacency_columns) > 0L) {
  abort_step( paste0( "Adjacency file is missing columns: ", paste( missing_adjacency_columns, collapse = ", " ) ) )
}
adjacency <- adjacency[ , ..required_adjacency_columns ]
adjacency[ , TF := as.character(TF) ]
adjacency[ , target := as.character(target) ]
adjacency[ , importance := as.numeric(importance) ]
adjacency_edges <- nrow(adjacency)
adjacency_TFs <- data.table::uniqueN(adjacency$TF)
adjacency_targets <- data.table::uniqueN(adjacency$target)
adjacency_bad_values <- sum( !is.finite(adjacency$importance) | adjacency$importance <= 0 )
adjacency_duplicate_pairs <- adjacency[ , .N, by = .( TF, target ) ][ N > 1L, .N ]
if ( adjacency_edges != EXPECTED_GRN_EDGES || adjacency_TFs != EXPECTED_GRN_TFS || adjacency_targets != EXPECTED_GRN_TARGETS || adjacency_bad_values != 0L || adjacency_duplicate_pairs != 0L ) {
  abort_step( paste0( "The Step 7B adjacency network does not match the frozen PASS result.\n", "Edges: ", adjacency_edges, "; TFs: ", adjacency_TFs, "; targets: ", adjacency_targets, "; bad importance: ", adjacency_bad_values, "; duplicate pairs: ", adjacency_duplicate_pairs ) )
}
PREFLIGHT_PY <- file.path( PYDIR, "36_cisTarget_preflight.py" )
preflight_lines <- c( "import sys", "from pathlib import Path", "", "import numpy as np", "import pandas as pd", "import loompy", "from ctxcore.rnkdb import FeatherRankingDatabase", "", "adjacency_file = Path(sys.argv[1])", "loom_file = Path(sys.argv[2])", "tf_file = Path(sys.argv[3])", "annotation_file = Path(sys.argv[4])", "feather_file = Path(sys.argv[5])", "audit_file = Path(sys.argv[6])", "expected_genes = int(sys.argv[7])", "expected_cells = int(sys.argv[8])", "", "adj = pd.read_csv(adjacency_file, sep='\\t', dtype={'TF': str, 'target': str})", "required = {'TF', 'target', 'importance'}", "if not required.issubset(adj.columns):", " raise RuntimeError(f'Adjacency columns missing: {required - set(adj.columns)}')", "", "with loompy.connect(str(loom_file), mode='r') as ds:", " loom_shape = tuple(ds.shape)", " genes = [str(x) for x in ds.ra.Gene]", " cells = [str(x) for x in ds.ca.CellID]", " branch_states = [str(x) for x in ds.ca.BranchState]", "", "if loom_shape != (expected_genes, expected_cells):", " raise RuntimeError(f'Unexpected loom shape: {loom_shape}')", "if len(set(genes)) != len(genes):", " raise RuntimeError('Duplicated loom genes detected')", "if len(set(cells)) != len(cells):", " raise RuntimeError('Duplicated loom cells detected')", "", "expected_states = {'Shared_prebranch': 116, 'Arm_1': 2289, 'Arm_2': 1798}", "observed_states = {state: branch_states.count(state) for state in expected_states}", "if observed_states != expected_states:", " raise RuntimeError(f'Frozen branch counts differ: {observed_states}')", "", "network_genes = set(adj['TF']).union(set(adj['target']))", "missing_network_genes = sorted(network_genes - set(genes))", "if missing_network_genes:", " raise RuntimeError(f'{len(missing_network_genes)} network genes missing from loom')", "", "with tf_file.open('r', encoding='utf-8-sig') as handle:", " tf_list = [line.strip() for line in handle if line.strip()]", "tf_overlap = sorted(set(adj['TF']).intersection(tf_list))", "", "database = FeatherRankingDatabase(fname=str(feather_file), name='mm10_10kb')", "db_genes = set(database.genes)", "db_overlap = set(genes).intersection(db_genes)", "", "annotation = pd.read_csv(annotation_file, sep='\\t', low_memory=False)", "if annotation.shape[0] < 1000:", " raise RuntimeError(f'Motif annotation unexpectedly small: {annotation.shape}')", "", "metrics = [", " ('loom_genes', len(genes)),", " ('loom_cells', len(cells)),", " ('state_Shared_prebranch', observed_states['Shared_prebranch']),", " ('state_Arm_1', observed_states['Arm_1']),", " ('state_Arm_2', observed_states['Arm_2']),", " ('adjacency_edges', int(adj.shape[0])),", " ('adjacency_TFs', int(adj['TF'].nunique())),", " ('adjacency_targets', int(adj['target'].nunique())),", " ('network_genes_missing_from_loom', len(missing_network_genes)),", " ('GRN_TFs_in_mouse_TF_list', len(tf_overlap)),", " ('feather_database_genes', len(db_genes)),", " ('loom_feather_overlap_genes', len(db_overlap)),", " ('loom_feather_overlap_fraction', len(db_overlap) / len(genes)),", " ('motif_annotation_rows', int(annotation.shape[0])),", " ('motif_annotation_columns', int(annotation.shape[1]))", "]", "", "with audit_file.open('w', encoding='utf-8', newline='\\n') as handle:", " handle.write('metric\\tvalue\\n')", " for metric, value in metrics:", " handle.write(f'{metric}\\t{value}\\n')", "", "print('STEP7C_PREFLIGHT_PASS')", "for metric, value in metrics:", " print(f'{metric}: {value}')" )
writeLines( preflight_lines, con = PREFLIGHT_PY, useBytes = TRUE )
preflight_log <- file.path( LOGDIR, "01_step7C_preflight.log" )
preflight_result <- run_command( command = ENV_PYTHON, args = c( quote_cmd(PREFLIGHT_PY), quote_cmd(ADJ_FILE), quote_cmd(EXPRESSION_LOOM), quote_cmd(TF_FILE), quote_cmd(MOTIF_ANNOTATION_FILE), quote_cmd(FEATHER_FILE), quote_cmd(PREFLIGHT_FILE), as.character(EXPECTED_GENES), as.character(EXPECTED_CELLS) ), log_file = preflight_log )
if ( preflight_result$status != 0L || !file.exists(PREFLIGHT_FILE) ) {
  abort_step( paste0( "Step 7C Python preflight failed.\n\n", read_log_tail( preflight_log, n = 200L ) ) )
}
preflight_metrics <- read_metric_tsv( PREFLIGHT_FILE )
preflight_table <- data.frame( metric = names(preflight_metrics), value = unname(preflight_metrics), stringsAsFactors = FALSE )
write_utf8_csv( preflight_table, file.path( TABDIR, "01_step7C_preflight_audit.csv" ) )
preflight_pass <- as.integer(preflight_metrics[["loom_genes"]]) == EXPECTED_GENES && as.integer(preflight_metrics[["loom_cells"]]) == EXPECTED_CELLS && as.integer(preflight_metrics[["adjacency_edges"]]) == EXPECTED_GRN_EDGES && as.integer(preflight_metrics[["adjacency_TFs"]]) == EXPECTED_GRN_TFS && as.integer(preflight_metrics[["adjacency_targets"]]) == EXPECTED_GRN_TARGETS && as.integer( preflight_metrics[["network_genes_missing_from_loom"]] ) == 0L && as.integer( preflight_metrics[["GRN_TFs_in_mouse_TF_list"]] ) == EXPECTED_GRN_TFS && as.numeric( preflight_metrics[["loom_feather_overlap_fraction"]] ) >= 0.90
if (!preflight_pass) {
  abort_step( paste0( "Step 7C preflight thresholds were not met.\n\n", paste( capture.output( print( preflight_table, row.names = FALSE ) ), collapse = "\n" ) ) )
}
smoke_status <- NA_integer_
smoke_elapsed <- NA_real_
smoke_pass <- TRUE
smoke_motif_size <- NA_real_
smoke_TFs_used <- character(0)
if (RUN_SMOKE_TEST) {
  tf_total_importance <- adjacency[ , .( total_importance = sum( importance, na.rm = TRUE ) ), by = TF ]
  data.table::setorder( tf_total_importance, -total_importance )
  smoke_TFs_used <- unique( c( intersect( SMOKE_TFS_PREFERRED, adjacency$TF ), tf_total_importance$TF ) )
  smoke_TFs_used <- utils::head( smoke_TFs_used, 12L )
  smoke_adjacency <- adjacency[ TF %in% smoke_TFs_used ][ order( TF, -importance ), utils::head( .SD, SMOKE_EDGES_PER_TF ), by = TF ]
  write_utf8_tsv( smoke_adjacency, SMOKE_ADJ_FILE )
  if (file.exists(SMOKE_MOTIF_FILE)) {
    file.remove(SMOKE_MOTIF_FILE)
  }
  smoke_log <- file.path( LOGDIR, "02_ctx_smoke_test.log" )
  smoke_args <- c( "ctx", quote_cmd(SMOKE_ADJ_FILE), quote_cmd(FEATHER_FILE), "--annotations_fname", quote_cmd(MOTIF_ANNOTATION_FILE), "--expression_mtx_fname", quote_cmd(EXPRESSION_LOOM), "--mode", CTX_MODE, "--num_workers", as.character(SMOKE_WORKERS), "--min_genes", as.character(MIN_GENES), "--rank_threshold", as.character(RANK_THRESHOLD), "--auc_threshold", as.character(AUC_THRESHOLD), "--nes_threshold", as.character(NES_THRESHOLD), "--chunk_size", "20", "--cell_id_attribute", "CellID", "--gene_attribute", "Gene", "--output", quote_cmd(SMOKE_MOTIF_FILE) )
  if (MASK_DROPOUTS) {
    smoke_args <- c( smoke_args, "--mask_dropouts" )
  }
  if (ALL_MODULES) {
    smoke_args <- c( smoke_args, "--all_modules" )
  }
  smoke_result <- run_command( command = PYSCENIC_EXE, args = smoke_args, log_file = smoke_log )
  smoke_status <- smoke_result$status
  smoke_elapsed <- smoke_result$elapsed_seconds
  smoke_motif_size <- if ( file.exists(SMOKE_MOTIF_FILE) ) {
    file.info(SMOKE_MOTIF_FILE)$size
  } else {
    NA_real_
  }
  smoke_pass <- smoke_status == 0L && valid_motif_output(SMOKE_MOTIF_FILE)
  smoke_audit <- data.frame( metric = c( "smoke_enabled", "smoke_TFs", "smoke_edges", "smoke_status_code", "smoke_elapsed_seconds", "smoke_elapsed_human", "smoke_output_size_bytes", "smoke_output_size_human", "smoke_pass" ), value = c( TRUE, length(smoke_TFs_used), nrow(smoke_adjacency), smoke_status, smoke_elapsed, format_duration(smoke_elapsed), smoke_motif_size, format_bytes(smoke_motif_size), smoke_pass ), stringsAsFactors = FALSE )
  write_utf8_csv( smoke_audit, file.path( TABDIR, "02_ctx_smoke_test_audit.csv" ) )
  if (!smoke_pass) {
    abort_step( paste0( "cisTarget smoke test failed. The full ctx run was not started.\n\n", read_log_tail( smoke_log, n = 220L ) ) )
  }
  rm( smoke_adjacency, tf_total_importance )
  gc()
} else {
  smoke_audit <- data.frame( metric = "smoke_enabled", value = FALSE, stringsAsFactors = FALSE )
  write_utf8_csv( smoke_audit, file.path( TABDIR, "02_ctx_smoke_test_audit.csv" ) )
}
rm(adjacency)
gc()
ctx_status <- NA_integer_
ctx_elapsed <- NA_real_
ctx_reused_existing <- FALSE
if ( !FORCE_RERUN_CTX && valid_motif_output(MOTIF_OUTPUT_FILE) ) {
  ctx_status <- 0L
  ctx_elapsed <- 0
  ctx_reused_existing <- TRUE
} else {
  if (file.exists(MOTIF_OUTPUT_FILE)) {
    file.remove(MOTIF_OUTPUT_FILE)
  }
  ctx_log <- file.path( LOGDIR, "03_pyscenic_ctx_full_run.log" )
  ctx_args <- c( "ctx", quote_cmd(ADJ_FILE), quote_cmd(FEATHER_FILE), "--annotations_fname", quote_cmd(MOTIF_ANNOTATION_FILE), "--expression_mtx_fname", quote_cmd(EXPRESSION_LOOM), "--mode", CTX_MODE, "--num_workers", as.character(NUM_WORKERS), "--min_genes", as.character(MIN_GENES), "--rank_threshold", as.character(RANK_THRESHOLD), "--auc_threshold", as.character(AUC_THRESHOLD), "--nes_threshold", as.character(NES_THRESHOLD), "--chunk_size", as.character(CHUNK_SIZE), "--cell_id_attribute", "CellID", "--gene_attribute", "Gene", "--output", quote_cmd(MOTIF_OUTPUT_FILE) )
  if (MASK_DROPOUTS) {
    ctx_args <- c( ctx_args, "--mask_dropouts" )
  }
  if (ALL_MODULES) {
    ctx_args <- c( ctx_args, "--all_modules" )
  }
  ctx_result <- run_command( command = PYSCENIC_EXE, args = ctx_args, log_file = ctx_log )
  ctx_status <- ctx_result$status
  ctx_elapsed <- ctx_result$elapsed_seconds
  if ( ctx_status != 0L || !valid_motif_output(MOTIF_OUTPUT_FILE) ) {
    abort_step( paste0( "The full pySCENIC ctx run failed or produced no valid motif table.\n\n", read_log_tail( ctx_log, n = 250L ) ) )
  }
}
AUDIT_PY <- file.path( PYDIR, "36_build_regulons.py" )
audit_lines <- c( "import json", "import pickle", "import sys", "from pathlib import Path", "", "import numpy as np", "import pandas as pd", "import loompy", "", "from ctxcore.genesig import GeneSignature", "from pyscenic.transform import df2regulons", "from pyscenic.utils import load_motifs", "", "motif_file = Path(sys.argv[1])", "adjacency_file = Path(sys.argv[2])", "loom_file = Path(sys.argv[3])", "dat_file = Path(sys.argv[4])", "gmt_file = Path(sys.argv[5])", "json_file = Path(sys.argv[6])", "metrics_file = Path(sys.argv[7])", "regulon_summary_file = Path(sys.argv[8])", "targets_file = Path(sys.argv[9])", "key_tf_file = Path(sys.argv[10])", "motif_tf_file = Path(sys.argv[11])", "key_tfs = [x for x in sys.argv[12].split(',') if x]", "", "motifs = load_motifs(str(motif_file), sep=',')", "if motifs.empty:", " raise RuntimeError('Motif enrichment table is empty')", "", "regulons = list(df2regulons(motifs))", "if len(regulons) == 0:", " raise RuntimeError('No regulons were derived from the motif table')", "", "with loompy.connect(str(loom_file), mode='r') as ds:", " loom_genes = set(str(x) for x in ds.ra.Gene)", "", "adj = pd.read_csv(", " adjacency_file,", " sep='\\t',", " dtype={'TF': str, 'target': str}", ")", "", "grn_summary = (", " adj.groupby('TF', sort=False)", " .agg(", " GRN_n_targets=('target', 'nunique'),", " GRN_total_importance=('importance', 'sum')", " )", " .reset_index()", ")", "", "summary_rows = []", "target_rows = []", "name_to_targets = {}", "", "for regulon in regulons:", " name = str(regulon.name)", " tf = str(getattr(regulon, 'transcription_factor', name.split('(')[0]))", " gene2weight = dict(regulon.gene2weight)", " targets = list(gene2weight.keys())", " weights = [float(gene2weight[g]) for g in targets]", " context = sorted(str(x) for x in getattr(regulon, 'context', []))", " sign = 'negative' if name.endswith('(-)') else 'positive'", "", " name_to_targets[name] = targets", "", " summary_rows.append({", " 'regulon': name,", " 'TF': tf,", " 'sign': sign,", " 'n_targets': len(targets),", " 'total_weight': float(np.sum(weights)) if weights else 0.0,", " 'mean_weight': float(np.mean(weights)) if weights else np.nan,", " 'median_weight': float(np.median(weights)) if weights else np.nan,", " 'maximum_weight': float(np.max(weights)) if weights else np.nan,", " 'context': ' | '.join(context)", " })", "", " for target in targets:", " target_rows.append({", " 'regulon': name,", " 'TF': tf,", " 'sign': sign,", " 'target': str(target),", " 'weight': float(gene2weight[target])", " })", "", "regulon_summary = pd.DataFrame(summary_rows)", "targets_long = pd.DataFrame(target_rows)", "", "if regulon_summary['regulon'].duplicated().any():", " raise RuntimeError('Duplicated regulon names detected')", "", "targets_outside_loom = sorted(set(targets_long['target']) - loom_genes)", "duplicate_regulon_targets = int(", " targets_long.duplicated(['regulon', 'target']).sum()", ")", "", "with dat_file.open('wb') as handle:", " pickle.dump(regulons, handle)", "", "GeneSignature.to_gmt(str(gmt_file), regulons)", "", "with json_file.open('w', encoding='utf-8') as handle:", " json.dump(name_to_targets, handle, ensure_ascii=False, indent=2)", "", "regulon_summary = regulon_summary.sort_values(", " ['n_targets', 'total_weight'],", " ascending=[False, False]", ")", "regulon_summary.to_csv(regulon_summary_file, sep='\\t', index=False)", "targets_long.to_csv(targets_file, sep='\\t', index=False)", "", "regulon_tf_summary = (", " regulon_summary.groupby('TF', sort=False)", " .agg(", " regulon_count=('regulon', 'nunique'),", " pruned_target_count=('n_targets', 'sum'),", " maximum_single_regulon_targets=('n_targets', 'max'),", " regulon_names=('regulon', lambda x: ' | '.join(sorted(set(x))))", " )", " .reset_index()", ")", "", "key_table = pd.DataFrame({'TF': key_tfs})", "key_table = key_table.merge(grn_summary, on='TF', how='left')", "key_table = key_table.merge(regulon_tf_summary, on='TF', how='left')", "key_table['present_in_GRN'] = key_table['GRN_n_targets'].notna()", "key_table['present_in_pruned_regulon'] = key_table['regulon_count'].notna()", "key_table['pruned_target_count'] = key_table['pruned_target_count'].fillna(0).astype(int)", "key_table['regulon_count'] = key_table['regulon_count'].fillna(0).astype(int)", "key_table['retention_fraction_vs_GRN_targets'] = np.where(", " key_table['GRN_n_targets'].fillna(0) > 0,", " key_table['pruned_target_count'] / key_table['GRN_n_targets'],", " np.nan", ")", "key_table.to_csv(key_tf_file, sep='\\t', index=False)", "", "motif_reset = motifs.reset_index()", "flattened_columns = []", "for col in motif_reset.columns:", " if isinstance(col, tuple):", " parts = [str(x) for x in col if str(x) and not str(x).startswith('Unnamed:')]", " flattened_columns.append('__'.join(parts))", " else:", " flattened_columns.append(str(col))", "motif_reset.columns = flattened_columns", "", "tf_column = next(", " (c for c in motif_reset.columns if c == 'TF' or c.endswith('__TF')),", " motif_reset.columns[0]", ")", "nes_column = next(", " (c for c in motif_reset.columns if c == 'NES' or c.endswith('__NES')),", " None", ")", "", "if nes_column is not None:", " motif_tf_summary = (", " motif_reset.groupby(tf_column)", " .agg(", " enriched_motif_rows=(tf_column, 'size'),", " maximum_NES=(nes_column, 'max'),", " median_NES=(nes_column, 'median')", " )", " .reset_index()", " .rename(columns={tf_column: 'TF'})", " )", "else:", " motif_tf_summary = (", " motif_reset.groupby(tf_column)", " .size()", " .reset_index(name='enriched_motif_rows')", " .rename(columns={tf_column: 'TF'})", " )", "", "motif_tf_summary = motif_tf_summary.sort_values(", " 'enriched_motif_rows',", " ascending=False", ")", "motif_tf_summary.to_csv(motif_tf_file, sep='\\t', index=False)", "", "unique_regulon_tfs = int(regulon_summary['TF'].nunique())", "unique_pruned_targets = int(targets_long['target'].nunique())", "positive_regulons = int((regulon_summary['sign'] == 'positive').sum())", "negative_regulons = int((regulon_summary['sign'] == 'negative').sum())", "key_tfs_retained = int(key_table['present_in_pruned_regulon'].sum())", "", "metrics = [", " ('motif_rows', int(motifs.shape[0])),", " ('motif_table_unique_TFs', int(motifs.index.get_level_values(0).nunique())),", " ('regulon_count', int(regulon_summary.shape[0])),", " ('unique_regulon_TFs', unique_regulon_tfs),", " ('positive_regulons', positive_regulons),", " ('negative_regulons', negative_regulons),", " ('pruned_regulon_target_edges', int(targets_long.shape[0])),", " ('unique_pruned_targets', unique_pruned_targets),", " ('minimum_regulon_targets', int(regulon_summary['n_targets'].min())),", " ('median_regulon_targets', float(regulon_summary['n_targets'].median())),", " ('mean_regulon_targets', float(regulon_summary['n_targets'].mean())),", " ('maximum_regulon_targets', int(regulon_summary['n_targets'].max())),", " ('original_GRN_edges', int(adj.shape[0])),", " ('original_GRN_TFs', int(adj['TF'].nunique())),", " ('TF_retention_fraction', unique_regulon_tfs / adj['TF'].nunique()),", " ('edge_retention_fraction', targets_long.shape[0] / adj.shape[0]),", " ('targets_outside_loom', len(targets_outside_loom)),", " ('duplicate_regulon_target_pairs', duplicate_regulon_targets),", " ('key_TFs_retained', key_tfs_retained),", " ('dat_file_size_bytes', dat_file.stat().st_size),", " ('gmt_file_size_bytes', gmt_file.stat().st_size),", " ('json_file_size_bytes', json_file.stat().st_size)", "]", "", "with metrics_file.open('w', encoding='utf-8', newline='\\n') as handle:", " handle.write('metric\\tvalue\\n')", " for metric, value in metrics:", " handle.write(f'{metric}\\t{value}\\n')", "", "print('STEP7C_REGULON_AUDIT_PASS')", "for metric, value in metrics:", " print(f'{metric}: {value}')" )
writeLines( audit_lines, con = AUDIT_PY, useBytes = TRUE )
audit_log <- file.path( LOGDIR, "04_build_regulons_and_audit.log" )
audit_result <- run_command( command = ENV_PYTHON, args = c( quote_cmd(AUDIT_PY), quote_cmd(MOTIF_OUTPUT_FILE), quote_cmd(ADJ_FILE), quote_cmd(EXPRESSION_LOOM), quote_cmd(REGULON_DAT_FILE), quote_cmd(REGULON_GMT_FILE), quote_cmd(REGULON_JSON_FILE), quote_cmd(AUDIT_METRICS_FILE), quote_cmd(REGULON_SUMMARY_FILE), quote_cmd(REGULON_TARGET_FILE), quote_cmd(KEY_TF_AUDIT_FILE), quote_cmd(MOTIF_TF_SUMMARY_FILE), paste( KEY_TFS, collapse = "," ) ), log_file = audit_log )
if ( audit_result$status != 0L || !file.exists(AUDIT_METRICS_FILE) ) {
  abort_step( paste0( "Regulon construction/audit failed.\n\n", read_log_tail( audit_log, n = 250L ) ) )
}
audit_metrics <- read_metric_tsv( AUDIT_METRICS_FILE )
audit_summary <- data.frame( metric = names(audit_metrics), value = unname(audit_metrics), stringsAsFactors = FALSE )
write_utf8_csv( audit_summary, file.path( TABDIR, "03_regulon_pruning_summary.csv" ) )
key_tf_audit <- data.table::fread( KEY_TF_AUDIT_FILE, data.table = FALSE, showProgress = FALSE )
regulon_files_pass <- all( file.exists( c( MOTIF_OUTPUT_FILE, REGULON_DAT_FILE, REGULON_GMT_FILE, REGULON_JSON_FILE, REGULON_SUMMARY_FILE, REGULON_TARGET_FILE, KEY_TF_AUDIT_FILE ) ) )
structural_regulon_pass <- as.integer(audit_metrics[["motif_rows"]]) >= MIN_MOTIF_ROWS && as.integer(audit_metrics[["regulon_count"]]) >= MIN_REGULONS && as.integer(audit_metrics[["unique_regulon_TFs"]]) >= MIN_REGULON_TFS && as.integer(audit_metrics[["minimum_regulon_targets"]]) >= 1L && as.integer(audit_metrics[["targets_outside_loom"]]) == 0L && as.integer( audit_metrics[["duplicate_regulon_target_pairs"]] ) == 0L
key_tf_audit_pass <- as.integer(audit_metrics[["key_TFs_retained"]]) >= MIN_KEY_TFS_RETAINED_FOR_AUDIT
final_decision <- data.frame( domain = c( "Required input/resource files", "pySCENIC ctx CLI", "Frozen Step 7B adjacency network", "Expression loom and branch counts", "Mouse TF-list compatibility", "mm10 feather gene overlap", "Technical cisTarget smoke test", "Full pySCENIC ctx process", "Motif enrichment output", "Regulon DAT/GMT/JSON outputs", "Minimum motif rows", "Minimum regulon count", "Minimum regulon TF count", "All regulon targets present in loom", "No duplicate regulon-target pairs", "Key TF audit produced", "Overall Step 7C" ), pass = c( all(input_file_audit$exists), ctx_cli_pass, adjacency_edges == EXPECTED_GRN_EDGES && adjacency_TFs == EXPECTED_GRN_TFS && adjacency_targets == EXPECTED_GRN_TARGETS, preflight_pass, as.integer( preflight_metrics[["GRN_TFs_in_mouse_TF_list"]] ) == EXPECTED_GRN_TFS, as.numeric( preflight_metrics[["loom_feather_overlap_fraction"]] ) >= 0.90, smoke_pass, ctx_status == 0L && valid_motif_output(MOTIF_OUTPUT_FILE), as.integer(audit_metrics[["motif_rows"]]) >= MIN_MOTIF_ROWS, regulon_files_pass, as.integer(audit_metrics[["motif_rows"]]) >= MIN_MOTIF_ROWS, as.integer(audit_metrics[["regulon_count"]]) >= MIN_REGULONS, as.integer(audit_metrics[["unique_regulon_TFs"]]) >= MIN_REGULON_TFS, as.integer(audit_metrics[["targets_outside_loom"]]) == 0L, as.integer( audit_metrics[["duplicate_regulon_target_pairs"]] ) == 0L, key_tf_audit_pass, preflight_pass && smoke_pass && ctx_status == 0L && regulon_files_pass && structural_regulon_pass && key_tf_audit_pass ), stringsAsFactors = FALSE )
final_decision$status <- ifelse( final_decision$pass, "PASS", "FAIL" )
write_utf8_csv( final_decision, FINAL_DECISION_FILE )
FINAL_STATUS <- if ( isTRUE( final_decision$pass[ final_decision$domain == "Overall Step 7C" ][[1]] ) ) {
  "PASS"
} else {
  "FAIL"
}
overall_elapsed <- as.numeric( difftime( Sys.time(), overall_start, units = "secs" ) )
writeLines( capture.output( sessionInfo() ), con = SESSION_FILE, useBytes = TRUE )
gc()
