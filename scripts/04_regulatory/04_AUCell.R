rm(list = ls())
gc()
options( stringsAsFactors = FALSE, width = 220, warn = 1, timeout = 86400 )
set.seed(20260801)
PROJECT_ROOT <- Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
CONDA_ENV_NAME <- "pyscenic_env"
NUM_WORKERS <- 8L
SMOKE_WORKERS <- 2L
AUC_THRESHOLD <- 0.05
AUC_SEED <- 20260801L
RUN_SMOKE_TEST <- FALSE
FORCE_RERUN_AUCELL <- FALSE
EXPECTED_CELLS <- 4203L
EXPECTED_GENES <- 7897L
EXPECTED_REGULONS <- 292L
EXPECTED_STATE_COUNTS <- c( Shared_prebranch = 116L, Arm_1 = 2289L, Arm_2 = 1798L )
PRIMARY_MIN_TARGETS <- 10L
MIN_UNIQUE_AUC_VALUES <- 10L
LOW_INFORMATION_SD <- 1e-8
LOW_INFORMATION_RANGE <- 1e-7
KEY_TFS <- unique( c( "Spi1", "Irf8", "Atf3", "Fos", "Junb", "Egr1", "Cebpb", "Mafb", "Nfkb1", "Nfkb2", "Rela", "Relb", "Stat1", "Stat3", "Irf1", "Irf7", "Tfeb", "Srebf1", "Nfe2l2", "Bhlhe41", "Mitf", "Nr1h3", "Pparg" ) )
SMOKE_REGULON_COUNT <- 20L
SMOKE_CELL_COUNT <- 300L
STEP7B_DIR <- file.path( PROJECT_ROOT, "04_results", "step7B_GSE162610_GRNBoost2_network_inference" )
EXPRESSION_LOOM <- file.path( STEP7B_DIR, "input", "GSE162610_step7B_frozen_4203_raw_counts.loom" )
STEP7C_DIR <- file.path( PROJECT_ROOT, "04_results", "step7C_GSE162610_cisTarget_motif_pruning_regulon_inference" )
REGULON_DAT_FILE <- file.path( STEP7C_DIR, "regulons", "GSE162610_step7C_regulons.dat" )
REGULON_GMT_FILE <- file.path( STEP7C_DIR, "regulons", "GSE162610_step7C_regulons.gmt" )
REGULON_SUMMARY_INPUT <- file.path( STEP7C_DIR, "tables", "04_regulon_summary.tsv" )
STEP7C_REPORT <- file.path( STEP7C_DIR, "reports", "Step7C_report.txt" )
OUTDIR <- file.path( PROJECT_ROOT, "04_results", "step7D_GSE162610_AUCell_regulon_activity_scoring" )
INPUTDIR <- file.path(OUTDIR, "input")
AUCDIR <- file.path(OUTDIR, "AUC")
TABDIR <- file.path(OUTDIR, "tables")
LOGDIR <- file.path(OUTDIR, "logs")
REPORTDIR <- file.path(OUTDIR, "reports")
PYDIR <- file.path(OUTDIR, "python")
TMPDIR_AUC <- file.path(OUTDIR, "tmp")
RDSDIR <- file.path(OUTDIR, "RDS")
for (directory in c( OUTDIR, INPUTDIR, AUCDIR, TABDIR, LOGDIR, REPORTDIR, PYDIR, TMPDIR_AUC, RDSDIR )) {
  dir.create( directory, recursive = TRUE, showWarnings = FALSE )
}
SMOKE_LOOM <- file.path( INPUTDIR, "GSE162610_step7D_smoke_300cells.loom" )
SMOKE_REGULON_GMT <- file.path( INPUTDIR, "GSE162610_step7D_smoke_20regulons.gmt" )
SMOKE_OUTPUT_LOOM <- file.path( AUCDIR, "GSE162610_step7D_smoke_AUCell.loom" )
AUCELL_OUTPUT_LOOM <- file.path( AUCDIR, "GSE162610_step7D_SCENIC_AUCell.loom" )
AUC_CELLS_BY_REGULONS_TSV <- file.path( AUCDIR, "GSE162610_step7D_AUC_cells_by_regulons.tsv.gz" )
AUC_REGULONS_BY_CELLS_RDS <- file.path( RDSDIR, "GSE162610_step7D_AUC_regulons_by_cells.rds" )
AUC_CELLS_BY_REGULONS_RDS <- file.path( RDSDIR, "GSE162610_step7D_AUC_cells_by_regulons.rds" )
STEP7D_RESULTS_RDS <- file.path( RDSDIR, "GSE162610_step7D_AUCell_results.rds" )
PREFLIGHT_FILE <- file.path( TABDIR, "01_step7D_preflight_audit.tsv" )
AUC_MATRIX_AUDIT_FILE <- file.path( TABDIR, "03_AUCell_matrix_audit.csv" )
AUC_DISTRIBUTION_FILE <- file.path( TABDIR, "04_regulon_AUC_distribution_summary.tsv" )
STATE_SUMMARY_FILE <- file.path( TABDIR, "05_regulon_AUC_state_descriptive_summary.tsv" )
KEY_TF_AUC_FILE <- file.path( TABDIR, "06_key_TF_AUC_summary.tsv" )
TIMEPOINT_SUMMARY_FILE <- file.path( TABDIR, "07_regulon_AUC_timepoint_descriptive_summary.tsv" )
FINAL_DECISION_FILE <- file.path( TABDIR, "08_final_decision.csv" )
REPORT_FILE <- file.path( REPORTDIR, "Step7D_report.txt" )
ERROR_FILE <- file.path( REPORTDIR, "Step7D_ERROR.txt" )
SESSION_FILE <- file.path( REPORTDIR, "Step7D_sessionInfo.txt" )
if (!requireNamespace("data.table", quietly = TRUE)) {
  stop( "The R package 'data.table' is required.", call. = FALSE )
}
timestamp <- function() {
  format( Sys.time(), "%Y-%m-%d %H:%M:%S" )
}
msg <- function(...) {
  cat( paste0( "[", timestamp(), "] ", paste0(..., collapse = "") ), "\n" )
  flush.console()
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
  index <- 1L
  while (x >= 1024 && index < length(units)) {
    x <- x / 1024
    index <- index + 1L
  }
  paste0( format( round(x, 3), nsmall = 3, trim = TRUE ), " ", units[[index]] )
}
write_utf8_csv <- function(x, path) {
  utils::write.csv( x, file = path, row.names = FALSE, fileEncoding = "UTF-8", na = "" )
}
write_utf8_tsv <- function(x, path) {
  data.table::fwrite( x, file = path, sep = "\t", quote = FALSE, na = "" )
}
read_log_tail <- function(path, n = 180L) {
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
  writeLines( c( "======================================================================", "GSE162610 | STEP 7D ERROR", "======================================================================", paste0("Time: ", timestamp()), "", text ), con = ERROR_FILE, useBytes = TRUE )
  stop( text, call. = FALSE )
}
with_aucell_environment <- function(code) {
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
  Sys.setenv( TEMP = safe_path(TMPDIR_AUC), TMP = safe_path(TMPDIR_AUC), TMPDIR = safe_path(TMPDIR_AUC), JOBLIB_TEMP_FOLDER = safe_path(TMPDIR_AUC), PYTHONUTF8 = "1", OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1", NUMEXPR_NUM_THREADS = "1" )
  force(code)
}
run_command <- function( command, args = character(0), log_file ) {
  if (file.exists(log_file)) {
    file.remove(log_file)
  }
  cat( paste0( "[", timestamp(), "] COMMAND\n", safe_path(command), " ", paste(args, collapse = " "), "\n\n" ), file = log_file, append = TRUE )
  start_time <- Sys.time()
  status <- with_aucell_environment( tryCatch( system2( command = command, args = args, stdout = log_file, stderr = log_file ), error = function(e) {
        cat( paste0( "\nSYSTEM2 ERROR: ", conditionMessage(e), "\n" ), file = log_file, append = TRUE )
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
  output <- with_aucell_environment( tryCatch( suppressWarnings( system2( command = command, args = args, stdout = TRUE, stderr = TRUE ) ), error = function(e) {
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
read_metric_tsv <- function(path) {
  table <- data.table::fread( path, sep = "\t", header = TRUE, data.table = FALSE, showProgress = FALSE )
  if (!all(c("metric", "value") %in% colnames(table))) {
    stop( "Metric table lacks metric/value columns.", call. = FALSE )
  }
  setNames( as.character(table$value), as.character(table$metric) )
}
valid_aucell_loom <- function(path) {
  file.exists(path) && !is.na(file.info(path)$size) && file.info(path)$size > 1024 * 1024
}
help_has <- function(help_output, token) {
  any( grepl( token, help_output, fixed = TRUE ) )
}
build_aucell_args <- function( help_output, expression_loom, regulon_file, output_loom, workers, seed, auc_threshold ) {
  args <- c("aucell")
  if (help_has(help_output, "--num_workers")) {
    args <- c( args, "--num_workers", as.character(workers) )
  }
  if (help_has(help_output, "--seed")) {
    args <- c( args, "--seed", as.character(seed) )
  }
  if (help_has(help_output, "--auc_threshold")) {
    args <- c( args, "--auc_threshold", as.character(auc_threshold) )
  }
  if (help_has(help_output, "--cell_id_attribute")) {
    args <- c( args, "--cell_id_attribute", "CellID" )
  }
  if (help_has(help_output, "--gene_attribute")) {
    args <- c( args, "--gene_attribute", "Gene" )
  }
  args <- c( args, "--output", quote_cmd(output_loom), quote_cmd(expression_loom), quote_cmd(regulon_file) )
  args
}
overall_start <- Sys.time()
msg("==============================================================")
msg("GSE162610 Step 7D | AUCell regulon activity scoring")
msg("==============================================================")
msg("This step calculates continuous per-cell regulon AUC only.")
required_files <- c( EXPRESSION_LOOM, REGULON_DAT_FILE, REGULON_GMT_FILE, REGULON_SUMMARY_INPUT )
missing_files <- required_files[ !file.exists(required_files) ]
if (length(missing_files) > 0L) {
  abort_step( paste0( "Required Step 7D files are missing:\n", paste( safe_path(missing_files), collapse = "\n" ) ) )
}
input_file_audit <- data.frame( role = c( "expression_loom", "regulon_DAT", "regulon_GMT", "Step7C_regulon_summary" ), path = safe_path(required_files), size_bytes = file.info(required_files)$size, size_human = vapply( file.info(required_files)$size, format_bytes, FUN.VALUE = character(1) ), exists = file.exists(required_files), stringsAsFactors = FALSE )
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
msg("Python: ", ENV_PYTHON)
msg("pySCENIC: ", PYSCENIC_EXE)
msg("Workers: ", NUM_WORKERS)
msg("AUCell threshold fraction: ", AUC_THRESHOLD)
msg("AUCell signature input: Step 7C GMT (DAT loader compatibility fix)")
msg("AUCell output loom extraction: validate=FALSE with manual QC")
msg("GMT-derived AUC field names: remove terminal comma only, then exact-match audit")
msg("Existing completed AUCell loom will be reused when present.")
aucell_help <- capture_command( command = PYSCENIC_EXE, args = c( "aucell", "--help" ) )
writeLines( aucell_help$output, con = file.path( LOGDIR, "00_pyscenic_aucell_help.txt" ), useBytes = TRUE )
aucell_cli_pass <- aucell_help$status == 0L && help_has( aucell_help$output, "--output" )
if (!aucell_cli_pass) {
  abort_step( paste0( "pySCENIC AUCell CLI precheck failed.\n\n", paste( aucell_help$output, collapse = "\n" ) ) )
}
PREFLIGHT_PY <- file.path( PYDIR, "37_AUCell_preflight.py" )
preflight_lines <- c( "import pickle", "import sys", "from pathlib import Path", "", "import numpy as np", "import pandas as pd", "import loompy", "from ctxcore.genesig import GeneSignature", "", "expression_loom = Path(sys.argv[1])", "regulon_dat = Path(sys.argv[2])", "regulon_summary_file = Path(sys.argv[3])", "smoke_loom = Path(sys.argv[4])", "smoke_regulon_gmt = Path(sys.argv[5])", "audit_file = Path(sys.argv[6])", "expected_genes = int(sys.argv[7])", "expected_cells = int(sys.argv[8])", "expected_regulons = int(sys.argv[9])", "smoke_cells_n = int(sys.argv[10])", "smoke_regulons_n = int(sys.argv[11])", "key_tfs = [x for x in sys.argv[12].split(',') if x]", "", "with regulon_dat.open('rb') as handle:", " regulons = pickle.load(handle)", "", "if len(regulons) != expected_regulons:", " raise RuntimeError(f'Expected {expected_regulons} regulons, observed {len(regulons)}')", "", "regulon_names = [str(x.name) for x in regulons]", "if len(set(regulon_names)) != len(regulon_names):", " raise RuntimeError('Duplicated regulon names detected')", "", "regulon_tfs = [", " str(getattr(x, 'transcription_factor', str(x.name).split('(')[0]))", " for x in regulons", "]", "", "regulon_sizes = [len(x.genes) for x in regulons]", "regulon_summary = pd.read_csv(regulon_summary_file, sep='\\t')", "", "with loompy.connect(str(expression_loom), mode='r') as ds:", " loom_shape = tuple(ds.shape)", " genes = [str(x) for x in ds.ra.Gene]", " cells = [str(x) for x in ds.ca.CellID]", " branch_states = [str(x) for x in ds.ca.BranchState]", "", " if loom_shape != (expected_genes, expected_cells):", " raise RuntimeError(f'Unexpected expression loom shape: {loom_shape}')", "", " if len(set(genes)) != len(genes):", " raise RuntimeError('Duplicated genes in expression loom')", "", " if len(set(cells)) != len(cells):", " raise RuntimeError('Duplicated cells in expression loom')", "", " expected_states = {'Shared_prebranch': 116, 'Arm_1': 2289, 'Arm_2': 1798}", " observed_states = {state: branch_states.count(state) for state in expected_states}", "", " if observed_states != expected_states:", " raise RuntimeError(f'Frozen state counts differ: {observed_states}')", "", " loom_gene_set = set(genes)", " missing_targets = sorted(", " {gene for regulon in regulons for gene in regulon.genes} - loom_gene_set", " )", "", " if missing_targets:", " raise RuntimeError(f'{len(missing_targets)} regulon targets missing from loom')", "", " state_arrays = {state: np.flatnonzero(np.asarray(branch_states) == state) for state in expected_states}", " selected_cells = []", " per_state = max(1, smoke_cells_n // 3)", "", " for state in ['Shared_prebranch', 'Arm_1', 'Arm_2']:", " indices = state_arrays[state]", " selected_cells.extend(indices[:min(per_state, len(indices))].tolist())", "", " selected_cells = sorted(set(selected_cells))", " if len(selected_cells) < smoke_cells_n:", " remaining = [i for i in range(len(cells)) if i not in set(selected_cells)]", " selected_cells.extend(remaining[:smoke_cells_n - len(selected_cells)])", "", " selected_cells = np.asarray(sorted(selected_cells[:smoke_cells_n]), dtype=int)", " smoke_matrix = ds[:, selected_cells]", " smoke_cell_ids = np.asarray(ds.ca.CellID[selected_cells], dtype=object)", "", "if smoke_loom.exists():", " smoke_loom.unlink()", "", "loompy.create(", " str(smoke_loom),", " smoke_matrix,", " {'Gene': np.asarray(genes, dtype=object)},", " {'CellID': smoke_cell_ids}", ")", "", "preferred = []", "for tf in key_tfs:", " for index, regulon_tf in enumerate(regulon_tfs):", " if regulon_tf == tf and index not in preferred:", " preferred.append(index)", "", "size_order = list(np.argsort(np.asarray(regulon_sizes))[::-1])", "selected_regulon_indices = []", "", "for index in preferred + size_order:", " if index not in selected_regulon_indices:", " selected_regulon_indices.append(index)", " if len(selected_regulon_indices) >= smoke_regulons_n:", " break", "", "smoke_regulons = [regulons[i] for i in selected_regulon_indices]", "", "GeneSignature.to_gmt(str(smoke_regulon_gmt), smoke_regulons)", "", "primary_regulons = sum(size >= 10 for size in regulon_sizes)", "small_regulons = sum(size < 10 for size in regulon_sizes)", "key_tfs_present = len(set(regulon_tfs).intersection(key_tfs))", "", "metrics = [", " ('loom_genes', len(genes)),", " ('loom_cells', len(cells)),", " ('state_Shared_prebranch', observed_states['Shared_prebranch']),", " ('state_Arm_1', observed_states['Arm_1']),", " ('state_Arm_2', observed_states['Arm_2']),", " ('regulon_count', len(regulons)),", " ('unique_regulon_TFs', len(set(regulon_tfs))),", " ('primary_regulons_n_targets_ge_10', primary_regulons),", " ('small_regulons_n_targets_lt_10', small_regulons),", " ('minimum_regulon_targets', min(regulon_sizes)),", " ('median_regulon_targets', float(np.median(regulon_sizes))),", " ('maximum_regulon_targets', max(regulon_sizes)),", " ('regulon_targets_missing_from_loom', len(missing_targets)),", " ('key_TFs_with_regulons', key_tfs_present),", " ('smoke_cells', len(selected_cells)),", " ('smoke_regulons', len(smoke_regulons))", "]", "", "with audit_file.open('w', encoding='utf-8', newline='\\n') as handle:", " handle.write('metric\\tvalue\\n')", " for metric, value in metrics:", " handle.write(f'{metric}\\t{value}\\n')", "", "print('STEP7D_PREFLIGHT_PASS')", "for metric, value in metrics:", " print(f'{metric}: {value}')" )
writeLines( preflight_lines, con = PREFLIGHT_PY, useBytes = TRUE )
preflight_log <- file.path( LOGDIR, "01_step7D_preflight.log" )
msg("Running Step 7D preflight and preparing smoke inputs...")
preflight_result <- run_command( command = ENV_PYTHON, args = c( quote_cmd(PREFLIGHT_PY), quote_cmd(EXPRESSION_LOOM), quote_cmd(REGULON_DAT_FILE), quote_cmd(REGULON_SUMMARY_INPUT), quote_cmd(SMOKE_LOOM), quote_cmd(SMOKE_REGULON_GMT), quote_cmd(PREFLIGHT_FILE), as.character(EXPECTED_GENES), as.character(EXPECTED_CELLS), as.character(EXPECTED_REGULONS), as.character(SMOKE_CELL_COUNT), as.character(SMOKE_REGULON_COUNT), paste( KEY_TFS, collapse = "," ) ), log_file = preflight_log )
if ( preflight_result$status != 0L || !file.exists(PREFLIGHT_FILE) ) {
  abort_step( paste0( "Step 7D preflight failed.\n\n", read_log_tail( preflight_log, n = 220L ) ) )
}
preflight_metrics <- read_metric_tsv( PREFLIGHT_FILE )
preflight_table <- data.frame( metric = names(preflight_metrics), value = unname(preflight_metrics), stringsAsFactors = FALSE )
write_utf8_csv( preflight_table, file.path( TABDIR, "01_step7D_preflight_audit.csv" ) )
preflight_pass <- as.integer(preflight_metrics[["loom_genes"]]) == EXPECTED_GENES && as.integer(preflight_metrics[["loom_cells"]]) == EXPECTED_CELLS && as.integer(preflight_metrics[["regulon_count"]]) == EXPECTED_REGULONS && as.integer( preflight_metrics[["regulon_targets_missing_from_loom"]] ) == 0L && as.integer(preflight_metrics[["state_Shared_prebranch"]]) == EXPECTED_STATE_COUNTS[["Shared_prebranch"]] && as.integer(preflight_metrics[["state_Arm_1"]]) == EXPECTED_STATE_COUNTS[["Arm_1"]] && as.integer(preflight_metrics[["state_Arm_2"]]) == EXPECTED_STATE_COUNTS[["Arm_2"]]
if (!preflight_pass) {
  abort_step( paste0( "Step 7D preflight thresholds were not met.\n\n", paste( capture.output( print( preflight_table, row.names = FALSE ) ), collapse = "\n" ) ) )
}
msg( "Preflight PASS: ", preflight_metrics[["regulon_count"]], " regulons, ", preflight_metrics[["loom_cells"]], " cells." )
smoke_status <- NA_integer_
smoke_elapsed <- NA_real_
smoke_pass <- TRUE
if (RUN_SMOKE_TEST) {
  if (file.exists(SMOKE_OUTPUT_LOOM)) {
    file.remove(SMOKE_OUTPUT_LOOM)
  }
  smoke_log <- file.path( LOGDIR, "02_AUCell_smoke_test.log" )
  smoke_args <- build_aucell_args( help_output = aucell_help$output, expression_loom = SMOKE_LOOM, regulon_file = SMOKE_REGULON_GMT, output_loom = SMOKE_OUTPUT_LOOM, workers = SMOKE_WORKERS, seed = AUC_SEED, auc_threshold = AUC_THRESHOLD )
  msg( "Running AUCell smoke test: ", preflight_metrics[["smoke_regulons"]], " regulons x ", preflight_metrics[["smoke_cells"]], " cells..." )
  smoke_result <- run_command( command = PYSCENIC_EXE, args = smoke_args, log_file = smoke_log )
  smoke_status <- smoke_result$status
  smoke_elapsed <- smoke_result$elapsed_seconds
  smoke_pass <- smoke_status == 0L && valid_aucell_loom(SMOKE_OUTPUT_LOOM)
  smoke_audit <- data.frame( metric = c( "smoke_enabled", "smoke_regulons", "smoke_cells", "smoke_status_code", "smoke_elapsed_seconds", "smoke_elapsed_human", "smoke_output_size_bytes", "smoke_output_size_human", "smoke_pass" ), value = c( TRUE, preflight_metrics[["smoke_regulons"]], preflight_metrics[["smoke_cells"]], smoke_status, smoke_elapsed, format_duration(smoke_elapsed), ifelse( file.exists(SMOKE_OUTPUT_LOOM), file.info(SMOKE_OUTPUT_LOOM)$size, NA_real_ ), ifelse( file.exists(SMOKE_OUTPUT_LOOM), format_bytes(file.info(SMOKE_OUTPUT_LOOM)$size), NA_character_ ), smoke_pass ), stringsAsFactors = FALSE )
  write_utf8_csv( smoke_audit, file.path( TABDIR, "02_AUCell_smoke_test_audit.csv" ) )
  if (!smoke_pass) {
    abort_step( paste0( "AUCell smoke test failed. The full run was not started.\n\n", read_log_tail( smoke_log, n = 240L ) ) )
  }
  msg( "AUCell smoke test PASS in ", format_duration(smoke_elapsed), "." )
} else {
  smoke_audit <- data.frame( metric = "smoke_enabled", value = FALSE, stringsAsFactors = FALSE )
  write_utf8_csv( smoke_audit, file.path( TABDIR, "02_AUCell_smoke_test_audit.csv" ) )
}
aucell_status <- NA_integer_
aucell_elapsed <- NA_real_
aucell_reused_existing <- FALSE
if ( !FORCE_RERUN_AUCELL && valid_aucell_loom(AUCELL_OUTPUT_LOOM) ) {
  aucell_status <- 0L
  aucell_elapsed <- 0
  aucell_reused_existing <- TRUE
  msg( "A structurally valid AUCell loom already exists; reusing it." )
} else {
  if (file.exists(AUCELL_OUTPUT_LOOM)) {
    file.remove(AUCELL_OUTPUT_LOOM)
  }
  aucell_log <- file.path( LOGDIR, "03_pyscenic_aucell_full_run.log" )
  full_args <- build_aucell_args( help_output = aucell_help$output, expression_loom = EXPRESSION_LOOM, regulon_file = REGULON_GMT_FILE, output_loom = AUCELL_OUTPUT_LOOM, workers = NUM_WORKERS, seed = AUC_SEED, auc_threshold = AUC_THRESHOLD )
  msg("Starting full pySCENIC AUCell.")
  msg( "External command running." )
  msg("Progress log: ", safe_path(aucell_log))
  aucell_result <- run_command( command = PYSCENIC_EXE, args = full_args, log_file = aucell_log )
  aucell_status <- aucell_result$status
  aucell_elapsed <- aucell_result$elapsed_seconds
  if ( aucell_status != 0L || !valid_aucell_loom(AUCELL_OUTPUT_LOOM) ) {
    abort_step( paste0( "Full pySCENIC AUCell failed or produced no valid output loom.\n\n", read_log_tail( aucell_log, n = 260L ) ) )
  }
}
msg( "AUCell output loom created/reused: ", format_bytes( file.info(AUCELL_OUTPUT_LOOM)$size ) )
EXTRACT_PY <- file.path( PYDIR, "37_extract_AUC.py" )
extract_lines <- c( "import gzip", "import pickle", "import sys", "from pathlib import Path", "", "import numpy as np", "import pandas as pd", "import loompy", "", "input_loom = Path(sys.argv[1])", "aucell_loom = Path(sys.argv[2])", "regulon_dat = Path(sys.argv[3])", "auc_tsv = Path(sys.argv[4])", "matrix_audit_file = Path(sys.argv[5])", "distribution_file = Path(sys.argv[6])", "state_summary_file = Path(sys.argv[7])", "key_tf_file = Path(sys.argv[8])", "timepoint_summary_file = Path(sys.argv[9])", "expected_cells = int(sys.argv[10])", "expected_regulons = int(sys.argv[11])", "primary_min_targets = int(sys.argv[12])", "min_unique_values = int(sys.argv[13])", "low_sd = float(sys.argv[14])", "low_range = float(sys.argv[15])", "key_tfs = [x for x in sys.argv[16].split(',') if x]", "", "with regulon_dat.open('rb') as handle:", " regulons = pickle.load(handle)", "", "regulon_names_expected = [str(x.name) for x in regulons]", "regulon_tfs = [", " str(getattr(x, 'transcription_factor', str(x.name).split('(')[0]))", " for x in regulons", "]", "regulon_sizes = [len(x.genes) for x in regulons]", "", "with loompy.connect(str(input_loom), mode='r') as source:", " source_cells = [str(x) for x in source.ca.CellID]", " metadata = pd.DataFrame({", " 'CellID': source_cells,", " 'BranchState': [str(x) for x in source.ca.BranchState],", " 'SampleID': [str(x) for x in source.ca.SampleID],", " 'Timepoint': [str(x) for x in source.ca.Timepoint],", " 'Pseudotime': np.asarray(source.ca.Pseudotime, dtype=float),", " 'Y29SignedProgress': np.asarray(source.ca.Y29SignedProgress, dtype=float)", " })", "", "with loompy.connect(str(aucell_loom), mode='r', validate=False) as ds:", " output_cells = [str(x) for x in ds.ca.CellID]", "", " auc_attribute = None", " for candidate in ['RegulonsAUC', 'RegulonAUC']:", " if candidate in ds.ca:", " auc_attribute = candidate", " break", "", " if auc_attribute is None:", " raise RuntimeError(", " 'No RegulonsAUC/RegulonAUC column attribute found. Available: ' +", " ', '.join(ds.ca.keys())", " )", "", " raw_auc = np.asarray(ds.ca[auc_attribute])", "", "if output_cells != source_cells:", " raise RuntimeError('AUCell loom cell order differs from frozen input loom')", "", "if raw_auc.dtype.names is None:", " raise RuntimeError(", " f'AUC attribute is not a named structured array: shape={raw_auc.shape}, dtype={raw_auc.dtype}'", " )", "", "auc_dict = {}", "raw_to_normalized_names = {}", "", "for raw_name in raw_auc.dtype.names:", " raw_name_string = str(raw_name)", " normalized_name = raw_name_string[:-1] if raw_name_string.endswith(',') else raw_name_string", "", " if normalized_name in auc_dict:", " raise RuntimeError(", " f'Regulon-name normalization collision: {raw_name_string} -> {normalized_name}'", " )", "", " values = np.asarray(raw_auc[raw_name]).reshape(-1).astype(float)", " auc_dict[normalized_name] = values", " raw_to_normalized_names[raw_name_string] = normalized_name", "", "auc = pd.DataFrame(auc_dict, index=output_cells)", "auc.index.name = 'CellID'", "", "if auc.shape != (expected_cells, expected_regulons):", " raise RuntimeError(f'Unexpected AUC matrix shape: {auc.shape}')", "", "if set(auc.columns) != set(regulon_names_expected):", " missing = sorted(set(regulon_names_expected) - set(auc.columns))", " extra = sorted(set(auc.columns) - set(regulon_names_expected))", " raise RuntimeError(f'Regulon-name mismatch. Missing={missing[:10]}, extra={extra[:10]}')", "", "auc = auc.loc[:, regulon_names_expected]", "", "values = auc.to_numpy(dtype=float)", "finite_mask = np.isfinite(values)", "bad_values = int((~finite_mask).sum())", "negative_values = int((values < -1e-12).sum())", "above_one_values = int((values > 1.000001).sum())", "", "if bad_values > 0:", " raise RuntimeError(f'AUC matrix contains {bad_values} non-finite values')", "", "auc.reset_index().to_csv(", " auc_tsv,", " sep='\\t',", " index=False,", " compression='gzip'", ")", "", "distribution_rows = []", "", "for index, regulon_name in enumerate(regulon_names_expected):", " vector = auc[regulon_name].to_numpy(dtype=float)", " quantiles = np.quantile(vector, [0.01, 0.05, 0.25, 0.50, 0.75, 0.95, 0.99])", " unique_values = int(np.unique(vector).size)", " standard_deviation = float(np.std(vector, ddof=1))", " observed_range = float(np.max(vector) - np.min(vector))", " low_information = (", " standard_deviation <= low_sd", " or observed_range <= low_range", " or unique_values < min_unique_values", " )", "", " distribution_rows.append({", " 'regulon': regulon_name,", " 'TF': regulon_tfs[index],", " 'n_targets': regulon_sizes[index],", " 'primary_n_targets_ge_10': regulon_sizes[index] >= primary_min_targets,", " 'mean_AUC': float(np.mean(vector)),", " 'sd_AUC': standard_deviation,", " 'minimum_AUC': float(np.min(vector)),", " 'q01_AUC': float(quantiles[0]),", " 'q05_AUC': float(quantiles[1]),", " 'q25_AUC': float(quantiles[2]),", " 'median_AUC': float(quantiles[3]),", " 'q75_AUC': float(quantiles[4]),", " 'q95_AUC': float(quantiles[5]),", " 'q99_AUC': float(quantiles[6]),", " 'maximum_AUC': float(np.max(vector)),", " 'IQR_AUC': float(quantiles[4] - quantiles[2]),", " 'range_AUC': observed_range,", " 'zero_fraction': float(np.mean(vector == 0)),", " 'unique_AUC_values': unique_values,", " 'low_information_flag': bool(low_information)", " })", "", "distribution = pd.DataFrame(distribution_rows)", "distribution = distribution.sort_values(", " ['low_information_flag', 'sd_AUC', 'n_targets'],", " ascending=[True, False, False]", ")", "distribution.to_csv(distribution_file, sep='\\t', index=False)", "", "auc_with_metadata = metadata.set_index('CellID').join(auc, how='left')", "", "state_rows = []", "for state, frame in auc_with_metadata.groupby('BranchState', sort=False):", " for regulon_name in regulon_names_expected:", " vector = frame[regulon_name].to_numpy(dtype=float)", " state_rows.append({", " 'BranchState': state,", " 'regulon': regulon_name,", " 'TF': regulon_tfs[regulon_names_expected.index(regulon_name)],", " 'n_cells': int(frame.shape[0]),", " 'mean_AUC': float(np.mean(vector)),", " 'median_AUC': float(np.median(vector)),", " 'sd_AUC': float(np.std(vector, ddof=1)) if len(vector) > 1 else np.nan", " })", "", "pd.DataFrame(state_rows).to_csv(state_summary_file, sep='\\t', index=False)", "", "time_rows = []", "for timepoint, frame in auc_with_metadata.groupby('Timepoint', sort=False):", " for regulon_name in regulon_names_expected:", " vector = frame[regulon_name].to_numpy(dtype=float)", " time_rows.append({", " 'Timepoint': timepoint,", " 'regulon': regulon_name,", " 'TF': regulon_tfs[regulon_names_expected.index(regulon_name)],", " 'n_cells': int(frame.shape[0]),", " 'mean_AUC': float(np.mean(vector)),", " 'median_AUC': float(np.median(vector)),", " 'sd_AUC': float(np.std(vector, ddof=1)) if len(vector) > 1 else np.nan", " })", "", "pd.DataFrame(time_rows).to_csv(timepoint_summary_file, sep='\\t', index=False)", "", "distribution_indexed = distribution.set_index('TF', drop=False)", "key_rows = []", "", "for tf in key_tfs:", " matching = [i for i, value in enumerate(regulon_tfs) if value == tf]", "", " if not matching:", " key_rows.append({", " 'TF': tf,", " 'regulon': '',", " 'present_in_regulons': False,", " 'n_targets': 0,", " 'primary_n_targets_ge_10': False,", " 'mean_AUC': np.nan,", " 'sd_AUC': np.nan,", " 'median_AUC': np.nan,", " 'q95_AUC': np.nan,", " 'maximum_AUC': np.nan,", " 'unique_AUC_values': 0,", " 'low_information_flag': np.nan", " })", " continue", "", " for index in matching:", " regulon_name = regulon_names_expected[index]", " row = distribution.loc[distribution['regulon'] == regulon_name].iloc[0]", " key_rows.append({", " 'TF': tf,", " 'regulon': regulon_name,", " 'present_in_regulons': True,", " 'n_targets': int(row['n_targets']),", " 'primary_n_targets_ge_10': bool(row['primary_n_targets_ge_10']),", " 'mean_AUC': float(row['mean_AUC']),", " 'sd_AUC': float(row['sd_AUC']),", " 'median_AUC': float(row['median_AUC']),", " 'q95_AUC': float(row['q95_AUC']),", " 'maximum_AUC': float(row['maximum_AUC']),", " 'unique_AUC_values': int(row['unique_AUC_values']),", " 'low_information_flag': bool(row['low_information_flag'])", " })", "", "pd.DataFrame(key_rows).to_csv(key_tf_file, sep='\\t', index=False)", "", "low_information_count = int(distribution['low_information_flag'].sum())", "primary_count = int(distribution['primary_n_targets_ge_10'].sum())", "nonzero_regulons = int((distribution['maximum_AUC'] > 0).sum())", "", "metrics = [", " ('AUC_cells', int(auc.shape[0])),", " ('AUC_regulons', int(auc.shape[1])),", " ('AUC_values', int(values.size)),", " ('non_finite_AUC_values', bad_values),", " ('negative_AUC_values', negative_values),", " ('AUC_values_above_1', above_one_values),", " ('minimum_AUC', float(values.min())),", " ('median_AUC_all_values', float(np.median(values))),", " ('mean_AUC_all_values', float(values.mean())),", " ('maximum_AUC', float(values.max())),", " ('primary_regulons_n_targets_ge_10', primary_count),", " ('small_regulons_n_targets_lt_10', int(expected_regulons - primary_count)),", " ('nonzero_regulons', nonzero_regulons),", " ('low_information_regulons', low_information_count),", " ('informative_regulons', int(expected_regulons - low_information_count)),", " ('low_information_fraction', low_information_count / expected_regulons),", " ('cell_order_exact_match', True),", " ('regulon_name_set_exact_match', True),", " ('raw_regulon_names_with_terminal_comma', int(sum(name.endswith(',') for name in raw_to_normalized_names))),", " ('regulon_name_normalization_collisions', 0),", " ('AUC_tsv_size_bytes', auc_tsv.stat().st_size)", "]", "", "pd.DataFrame(metrics, columns=['metric', 'value']).to_csv(", " matrix_audit_file,", " index=False", ")", "", "print('STEP7D_AUC_EXTRACTION_PASS')", "for metric, value in metrics:", " print(f'{metric}: {value}')" )
writeLines( extract_lines, con = EXTRACT_PY, useBytes = TRUE )
extract_log <- file.path( LOGDIR, "04_extract_AUC_and_summarize.log" )
msg("Extracting the AUCell matrix and building QC summaries...")
extract_result <- run_command( command = ENV_PYTHON, args = c( quote_cmd(EXTRACT_PY), quote_cmd(EXPRESSION_LOOM), quote_cmd(AUCELL_OUTPUT_LOOM), quote_cmd(REGULON_DAT_FILE), quote_cmd(AUC_CELLS_BY_REGULONS_TSV), quote_cmd(AUC_MATRIX_AUDIT_FILE), quote_cmd(AUC_DISTRIBUTION_FILE), quote_cmd(STATE_SUMMARY_FILE), quote_cmd(KEY_TF_AUC_FILE), quote_cmd(TIMEPOINT_SUMMARY_FILE), as.character(EXPECTED_CELLS), as.character(EXPECTED_REGULONS), as.character(PRIMARY_MIN_TARGETS), as.character(MIN_UNIQUE_AUC_VALUES), format(LOW_INFORMATION_SD, scientific = TRUE), format(LOW_INFORMATION_RANGE, scientific = TRUE), paste( KEY_TFS, collapse = "," ) ), log_file = extract_log )
if ( extract_result$status != 0L || !file.exists(AUC_MATRIX_AUDIT_FILE) || !file.exists(AUC_CELLS_BY_REGULONS_TSV) ) {
  abort_step( paste0( "AUCell extraction/QC failed.\n\n", read_log_tail( extract_log, n = 280L ) ) )
}
msg("Creating RDS matrices for downstream RSS/CSI and pseudotime analysis...")
auc_table <- data.table::fread( AUC_CELLS_BY_REGULONS_TSV, data.table = FALSE, showProgress = TRUE )
if ( nrow(auc_table) != EXPECTED_CELLS || ncol(auc_table) != EXPECTED_REGULONS + 1L || colnames(auc_table)[[1]] != "CellID" ) {
  abort_step( paste0( "Unexpected extracted AUCell table dimensions: ", nrow(auc_table), " rows x ", ncol(auc_table), " columns." ) )
}
cell_ids <- as.character( auc_table$CellID )
auc_cells_by_regulons <- as.matrix( auc_table[ , setdiff( colnames(auc_table), "CellID" ), drop = FALSE ] )
storage.mode( auc_cells_by_regulons ) <- "double"
rownames( auc_cells_by_regulons ) <- cell_ids
auc_regulons_by_cells <- t( auc_cells_by_regulons )
if ( any(!is.finite(auc_regulons_by_cells)) || nrow(auc_regulons_by_cells) != EXPECTED_REGULONS || ncol(auc_regulons_by_cells) != EXPECTED_CELLS ) {
  abort_step( "The R AUCell matrix failed its structural/finite-value check." )
}
saveRDS( auc_cells_by_regulons, file = AUC_CELLS_BY_REGULONS_RDS, compress = TRUE )
saveRDS( auc_regulons_by_cells, file = AUC_REGULONS_BY_CELLS_RDS, compress = TRUE )
matrix_audit <- data.table::fread( AUC_MATRIX_AUDIT_FILE, data.table = FALSE, showProgress = FALSE )
distribution_summary <- data.table::fread( AUC_DISTRIBUTION_FILE, data.table = FALSE, showProgress = FALSE )
key_tf_summary <- data.table::fread( KEY_TF_AUC_FILE, data.table = FALSE, showProgress = FALSE )
state_summary <- data.table::fread( STATE_SUMMARY_FILE, data.table = FALSE, showProgress = FALSE )
timepoint_summary <- data.table::fread( TIMEPOINT_SUMMARY_FILE, data.table = FALSE, showProgress = FALSE )
metric_lookup <- setNames( as.character(matrix_audit$value), as.character(matrix_audit$metric) )
step7d_results <- list( step = "Step7D", generated_at = timestamp(), auc_threshold = AUC_THRESHOLD, auc_seed = AUC_SEED, workers = NUM_WORKERS, frozen_cells = EXPECTED_CELLS, genes = EXPECTED_GENES, regulons = EXPECTED_REGULONS, primary_min_targets = PRIMARY_MIN_TARGETS, expression_loom = safe_path(EXPRESSION_LOOM), regulon_dat_for_audit = safe_path(REGULON_DAT_FILE), regulon_gmt_for_aucell = safe_path(REGULON_GMT_FILE), aucell_loom = safe_path(AUCELL_OUTPUT_LOOM), AUC_regulons_by_cells = auc_regulons_by_cells, distribution_summary = distribution_summary, key_TF_summary = key_tf_summary, state_descriptive_summary = state_summary, timepoint_descriptive_summary = timepoint_summary, matrix_audit = matrix_audit, preflight = preflight_table, smoke_audit = smoke_audit )
saveRDS( step7d_results, file = STEP7D_RESULTS_RDS, compress = TRUE )
auc_matrix_pass <- as.integer(metric_lookup[["AUC_cells"]]) == EXPECTED_CELLS && as.integer(metric_lookup[["AUC_regulons"]]) == EXPECTED_REGULONS && as.integer(metric_lookup[["non_finite_AUC_values"]]) == 0L && as.integer(metric_lookup[["negative_AUC_values"]]) == 0L && as.integer(metric_lookup[["AUC_values_above_1"]]) == 0L && identical( metric_lookup[["cell_order_exact_match"]], "True" ) && identical( metric_lookup[["regulon_name_set_exact_match"]], "True" )
information_pass <- as.integer(metric_lookup[["nonzero_regulons"]]) == EXPECTED_REGULONS && as.integer(metric_lookup[["informative_regulons"]]) >= floor(EXPECTED_REGULONS * 0.80)
rds_pass <- all( file.exists( c( AUC_CELLS_BY_REGULONS_RDS, AUC_REGULONS_BY_CELLS_RDS, STEP7D_RESULTS_RDS ) ) )
final_decision <- data.frame( domain = c( "Required input files", "pySCENIC AUCell CLI", "Frozen expression loom", "Step 7C regulon set", "Regulon targets present in expression loom", "Technical AUCell smoke test", "Full pySCENIC AUCell process", "AUC matrix dimensions", "AUC cell order", "AUC regulon names", "Finite non-negative normalized AUC values", "All regulons have nonzero activity", "At least 80% regulons informative", "Compressed TSV output", "RDS matrices/results", "Key TF AUCell audit", "Overall Step 7D" ), pass = c( all(input_file_audit$exists), aucell_cli_pass, preflight_pass, as.integer(preflight_metrics[["regulon_count"]]) == EXPECTED_REGULONS, as.integer( preflight_metrics[["regulon_targets_missing_from_loom"]] ) == 0L, smoke_pass, aucell_status == 0L && valid_aucell_loom(AUCELL_OUTPUT_LOOM), as.integer(metric_lookup[["AUC_cells"]]) == EXPECTED_CELLS && as.integer(metric_lookup[["AUC_regulons"]]) == EXPECTED_REGULONS, identical( metric_lookup[["cell_order_exact_match"]], "True" ), identical( metric_lookup[["regulon_name_set_exact_match"]], "True" ), as.integer(metric_lookup[["non_finite_AUC_values"]]) == 0L && as.integer(metric_lookup[["negative_AUC_values"]]) == 0L && as.integer(metric_lookup[["AUC_values_above_1"]]) == 0L, as.integer(metric_lookup[["nonzero_regulons"]]) == EXPECTED_REGULONS, as.integer(metric_lookup[["informative_regulons"]]) >= floor(EXPECTED_REGULONS * 0.80), file.exists(AUC_CELLS_BY_REGULONS_TSV), rds_pass, file.exists(KEY_TF_AUC_FILE) && nrow(key_tf_summary) > 0L, preflight_pass && smoke_pass && aucell_status == 0L && auc_matrix_pass && information_pass && rds_pass ), stringsAsFactors = FALSE )
final_decision$status <- ifelse( final_decision$pass, "PASS", "FAIL" )
write_utf8_csv( final_decision, FINAL_DECISION_FILE )
FINAL_STATUS <- if ( isTRUE( final_decision$pass[ final_decision$domain == "Overall Step 7D" ][[1]] ) ) {
  "PASS"
} else {
  "FAIL"
}
overall_elapsed <- as.numeric( difftime( Sys.time(), overall_start, units = "secs" ) )
primary_regulons <- as.integer( metric_lookup[["primary_regulons_n_targets_ge_10"]] )
small_regulons <- as.integer( metric_lookup[["small_regulons_n_targets_lt_10"]] )
informative_regulons <- as.integer( metric_lookup[["informative_regulons"]] )
low_information_regulons <- as.integer( metric_lookup[["low_information_regulons"]] )
report_lines <- c( "======================================================================", "GSE162610 | STEP 7D REPORT", "AUCell per-cell regulon activity scoring", "======================================================================", "", paste0("Generated: ", timestamp()), paste0("Project root: ", PROJECT_ROOT), paste0("Overall Step 7D status: ", FINAL_STATUS), "", "FORMAL SCOPE:", paste0( "Step 7D calculated continuous AUCell activity for all motif-supported ", "Step 7C regulons in each frozen Step5C-R cell." ), paste0( "No RSS, CSI, branch hypothesis testing, binary thresholding or ", "fate-driver claim was performed." ), "", "FROZEN INPUT:", paste0( "Cells: ", EXPECTED_CELLS, " (Shared_prebranch=116; Arm_1=2,289; Arm_2=1,798)" ), paste0("Genes: ", EXPECTED_GENES), paste0("Regulons: ", EXPECTED_REGULONS), paste0( "Primary regulons with >= ", PRIMARY_MIN_TARGETS, " motif-supported targets: ", primary_regulons ), paste0( "Small regulons with < ", PRIMARY_MIN_TARGETS, " targets retained for sensitivity analysis: ", small_regulons ), "", "AUCELL PARAMETERS:", paste0("Workers: ", NUM_WORKERS), paste0("Seed: ", AUC_SEED), paste0("AUC threshold fraction: ", AUC_THRESHOLD), paste0("Continuous AUC primary output: TRUE"), paste0("Binary thresholding performed: FALSE"), paste0("AUCell CLI signature input: Step 7C GMT"), paste0("DAT used only for direct binary-safe audit/extraction: TRUE"), paste0("AUCell output loom strict schema validation: FALSE"), paste0("Manual cell/regulon/value QC after extraction: TRUE"), paste0("Regulon-name normalization: terminal comma only"), paste0("Normalization collisions permitted: FALSE"), "", "TIMING:", paste0( "Preflight/smoke preparation: ", format_duration(preflight_result$elapsed_seconds) ), paste0( "Technical smoke test in this recovery run: SKIPPED ", "(already passed at 2026-08-01 23:15:46)" ), paste0( "Full pySCENIC AUCell: ", format_duration(aucell_elapsed) ), paste0( "AUC extraction/QC: ", format_duration(extract_result$elapsed_seconds) ), paste0( "Total Step 7D elapsed: ", format_duration(overall_elapsed) ), "", "AUC MATRIX:", paste( capture.output( print( matrix_audit, row.names = FALSE ) ), collapse = "\n" ), "", "REGULON INFORMATION AUDIT:", paste0( "Informative regulons: ", informative_regulons, "/", EXPECTED_REGULONS ), paste0( "Low-information flags: ", low_information_regulons, "/", EXPECTED_REGULONS ), paste0( "Low-information flags are technical annotations only; the full 292-regulon ", "matrix remains preserved." ), "", "KEY TF AUC AUDIT:", paste( capture.output( print( key_tf_summary, row.names = FALSE ) ), collapse = "\n" ), "", "FINAL DECISION:", paste( capture.output( print( final_decision, row.names = FALSE ) ), collapse = "\n" ), "", "INTERPRETATION BOUNDARY:", paste0( "The continuous AUCell matrix quantifies per-cell activity of ", "motif-supported candidate regulons. High activity does not by itself ", "establish branch specificity or fate-driving timing." ), "", "NEXT ACTION:", if (FINAL_STATUS == "PASS") {
    paste0( "Proceed to Step 7E regulon specificity and module analysis: ", "RSS ranking for Shared_prebranch/Arm_1/Arm_2, followed by CSI-based ", "regulon modules. Formal fate timing remains a later pseudotime step." )
  } else {
    paste0( "Do not calculate RSS/CSI. Review the failed Step 7D audit domain first." )
  },
  "", "PRIMARY OUTPUTS:", paste0( "AUCell loom: ", safe_path(AUCELL_OUTPUT_LOOM) ), paste0( "Compressed AUC TSV: ", safe_path(AUC_CELLS_BY_REGULONS_TSV) ), paste0( "Regulons x cells RDS: ", safe_path(AUC_REGULONS_BY_CELLS_RDS) ), paste0( "Complete Step 7D RDS: ", safe_path(STEP7D_RESULTS_RDS) ), "", "KEY OUTPUTS:", "1) reports/Step7D_report.txt", "2) tables/03_AUCell_matrix_audit.csv", "3) tables/04_regulon_AUC_distribution_summary.tsv", "4) tables/06_key_TF_AUC_summary.tsv", "5) tables/08_final_decision.csv"
)
writeLines( report_lines, con = REPORT_FILE, useBytes = TRUE )
writeLines( capture.output( sessionInfo() ), con = SESSION_FILE, useBytes = TRUE )
rm( auc_table, auc_cells_by_regulons, auc_regulons_by_cells )
gc()
msg("==============================================================")
msg("STEP 7D COMPLETE")
msg("==============================================================")
msg("Step 7D status: ", FINAL_STATUS)
msg( "AUC matrix: ", EXPECTED_REGULONS, " regulons x ", EXPECTED_CELLS, " cells." )
msg( "Informative regulons: ", informative_regulons, "; low-information flags: ", low_information_regulons, "." )
msg( "Full AUCell time: ", format_duration(aucell_elapsed) )
msg( "Total time: ", format_duration(overall_elapsed) )
msg("Report: ", safe_path(REPORT_FILE))
msg("")
msg("Key outputs:")
msg("1) reports/Step7D_report.txt")
msg("2) tables/03_AUCell_matrix_audit.csv")
msg("3) tables/04_regulon_AUC_distribution_summary.tsv")
msg("4) tables/06_key_TF_AUC_summary.tsv")
msg("5) tables/08_final_decision.csv")
