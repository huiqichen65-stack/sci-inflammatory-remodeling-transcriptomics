# ==============================================================================
# 34_step7A3_GSE162610_pySCENIC_ctxcore_resource_repair.R
# ==============================================================================
# Project:
#   GSE162610 adult SCI activated/DAM-like microglia fate-bifurcation project
#
# STEP 7A-3
#   Targeted repair of:
#     1) ctxcore / pkg_resources import failure
#     2) mouse motif-annotation download failure
#     3) mm10 10 kb feather download failure
#     4) final Step 7A re-audit
#
# This script REUSES the existing conda environment:
#   <USERPROFILE>/anaconda3/envs/pyscenic_env
#
# It DOES NOT:
#   - delete/recreate the conda environment
#   - run GRNBoost2
#   - run motif pruning
#   - run AUCell
#   - alter the frozen 4,203-cell set
#
# Run in RStudio:
#   source("<PROJECT_ROOT>/34_step7A3_GSE162610_pySCENIC_ctxcore_resource_repair.R")
# ==============================================================================

rm(list = ls())
gc()

options(
  stringsAsFactors = FALSE,
  width = 220,
  warn = 1,
  timeout = 86400
)

# ==============================================================================
# 0. SETTINGS
# ==============================================================================

PROJECT_ROOT <- Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
CONDA_ENV_NAME <- "pyscenic_env"

RESOURCE_DIR <- file.path(
  PROJECT_ROOT,
  "01_raw_data",
  "pySCENIC_resources"
)

OUTDIR <- file.path(
  PROJECT_ROOT,
  "04_results",
  "step7A3_GSE162610_pySCENIC_ctxcore_resource_repair"
)

LOGDIR <- file.path(OUTDIR, "logs")
TABDIR <- file.path(OUTDIR, "tables")
REPORTDIR <- file.path(OUTDIR, "reports")

for (dd in c(RESOURCE_DIR, OUTDIR, LOGDIR, TABDIR, REPORTDIR)) {
  dir.create(
    dd,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

REPORT_FILE <- file.path(
  REPORTDIR,
  "Step7A3_report.txt"
)

ERROR_FILE <- file.path(
  REPORTDIR,
  "Step7A3_ERROR.txt"
)

SESSION_FILE <- file.path(
  REPORTDIR,
  "Step7A3_sessionInfo.txt"
)

STEP7A_SCRIPT <- file.path(
  PROJECT_ROOT,
  "32_step7A_GSE162610_pySCENIC_input_environment_audit.R"
)

STEP7A_REPORT <- file.path(
  PROJECT_ROOT,
  "04_results",
  "step7A_GSE162610_pySCENIC_input_environment_audit",
  "reports",
  "Step7A_report.txt"
)

STEP7A_GENE_FILE <- file.path(
  PROJECT_ROOT,
  "04_results",
  "step7A_GSE162610_pySCENIC_input_environment_audit",
  "input",
  "GSE162610_step7A_frozen_4203_genes.tsv"
)

# One 10 kb database is sufficient for the current discovery analysis.
MOTIF_URL <- paste0(
  "https://resources.aertslab.org/cistarget/motif2tf/",
  "motifs-v9-nr.mgi-m0.001-o0.0.tbl"
)

FEATHER_URL <- paste0(
  "https://resources.aertslab.org/cistarget/databases/",
  "mus_musculus/mm10/refseq_r80/mc9nr/gene_based/",
  "mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.",
  "genes_vs_motifs.rankings.feather"
)

MOTIF_FILE <- file.path(
  RESOURCE_DIR,
  basename(MOTIF_URL)
)

FEATHER_FILE <- file.path(
  RESOURCE_DIR,
  basename(FEATHER_URL)
)

TF_FILE <- file.path(
  RESOURCE_DIR,
  "mm_mgi_tfs.txt"
)

# Conservative minimum sizes used only to reject empty/partial/HTML files.
MOTIF_MIN_BYTES <- 80 * 1024^2
FEATHER_MIN_BYTES <- 800 * 1024^2
TF_MIN_BYTES <- 5000

# Download controls.
MAX_DOWNLOAD_ROUNDS <- 4L
CURL_CONNECT_TIMEOUT <- 180L
CURL_RETRY_COUNT <- 10L
CURL_RETRY_DELAY <- 15L
CURL_SPEED_TIME <- 900L
CURL_SPEED_LIMIT <- 512L

# ==============================================================================
# 1. UTILITIES
# ==============================================================================

timestamp <- function() {
  format(
    Sys.time(),
    "%Y-%m-%d %H:%M:%S"
  )
}

msg <- function(...) {
  cat(
    paste0(
      "[",
      timestamp(),
      "] ",
      paste0(..., collapse = "")
    ),
    "\n"
  )
  flush.console()
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

  paste0(
    format(
      round(x, 3),
      nsmall = 3,
      trim = TRUE
    ),
    " ",
    units[[idx]]
  )
}

write_utf8_csv <- function(x, path) {
  utils::write.csv(
    x,
    file = path,
    row.names = FALSE,
    fileEncoding = "UTF-8",
    na = ""
  )
}

safe_path <- function(path) {
  normalizePath(
    path,
    winslash = "/",
    mustWork = FALSE
  )
}

quote_cmd <- function(x) {
  shQuote(
    as.character(x),
    type = "cmd"
  )
}

read_log_tail <- function(path, n = 80L) {
  if (!file.exists(path)) {
    return("<log not found>")
  }

  lines <- readLines(
    path,
    warn = FALSE,
    encoding = "UTF-8"
  )

  if (length(lines) == 0) {
    return("<empty log>")
  }

  paste(
    utils::tail(lines, n),
    collapse = "\n"
  )
}

abort_step <- function(text) {
  writeLines(
    c(
      "======================================================================",
      "GSE162610 | STEP 7A-3 ERROR",
      "======================================================================",
      paste0("Time: ", timestamp()),
      "",
      text
    ),
    con = ERROR_FILE,
    useBytes = TRUE
  )

  stop(
    text,
    call. = FALSE
  )
}

run_command <- function(
  command,
  args = character(0),
  log_file,
  append = FALSE,
  echo_command = TRUE
) {
  if (!append && file.exists(log_file)) {
    file.remove(log_file)
  }

  dir.create(
    dirname(log_file),
    recursive = TRUE,
    showWarnings = FALSE
  )

  if (echo_command) {
    cat(
      paste0(
        "[",
        timestamp(),
        "] COMMAND\n",
        safe_path(command),
        " ",
        paste(args, collapse = " "),
        "\n\n"
      ),
      file = log_file,
      append = TRUE
    )
  }

  start_time <- Sys.time()

  status <- tryCatch(
    system2(
      command = command,
      args = args,
      stdout = log_file,
      stderr = log_file
    ),
    error = function(e) {
      cat(
        paste0(
          "\nSYSTEM2 ERROR: ",
          conditionMessage(e),
          "\n"
        ),
        file = log_file,
        append = TRUE
      )

      999L
    }
  )

  if (is.null(status)) {
    status <- 0L
  }

  elapsed <- as.numeric(
    difftime(
      Sys.time(),
      start_time,
      units = "secs"
    )
  )

  list(
    status = as.integer(status),
    elapsed_seconds = elapsed,
    log_file = log_file
  )
}

capture_command <- function(
  command,
  args = character(0)
) {
  output <- tryCatch(
    suppressWarnings(
      system2(
        command = command,
        args = args,
        stdout = TRUE,
        stderr = TRUE
      )
    ),
    error = function(e) {
      structure(
        paste0(
          "SYSTEM2_ERROR: ",
          conditionMessage(e)
        ),
        status = 999L
      )
    }
  )

  status <- attr(output, "status")

  if (is.null(status)) {
    status <- 0L
  }

  list(
    status = as.integer(status),
    output = as.character(output)
  )
}

find_conda <- function() {
  candidates <- unique(
    c(
      Sys.getenv("CONDA_EXE"),
      unname(Sys.which("conda")),
      file.path(
        Sys.getenv("USERPROFILE"),
        "anaconda3",
        "Scripts",
        "conda.exe"
      ),
      file.path(
        Sys.getenv("USERPROFILE"),
        "miniconda3",
        "Scripts",
        "conda.exe"
      ),
      "C:/ProgramData/anaconda3/Scripts/conda.exe",
      "C:/ProgramData/miniconda3/Scripts/conda.exe"
    )
  )

  candidates <- candidates[
    nzchar(candidates) &
      file.exists(candidates)
  ]

  if (length(candidates) == 0) {
    return(NA_character_)
  }

  candidates[[1]]
}

conda_run_args <- function(...) {
  c(
    "run",
    "-n",
    CONDA_ENV_NAME,
    ...
  )
}

get_env_python <- function(conda_exe) {
  probe <- capture_command(
    conda_exe,
    conda_run_args(
      "python",
      "-c",
      quote_cmd(
        "import sys; print(sys.executable)"
      )
    )
  )

  if (probe$status != 0L) {
    return(NA_character_)
  }

  paths <- trimws(probe$output)
  paths <- paths[nzchar(paths)]
  paths <- paths[file.exists(paths)]

  if (length(paths) == 0) {
    return(NA_character_)
  }

  safe_path(paths[[length(paths)]])
}

is_html_or_error_text <- function(path) {
  if (
    !file.exists(path) ||
    is.na(file.info(path)$size) ||
    file.info(path)$size == 0
  ) {
    return(TRUE)
  }

  con <- file(
    path,
    open = "rb"
  )

  on.exit(
    close(con),
    add = TRUE
  )

  raw_head <- readBin(
    con,
    what = "raw",
    n = 8192
  )

  if (length(raw_head) == 0) {
    return(TRUE)
  }

  text_head <- tryCatch(
    rawToChar(raw_head),
    error = function(e) ""
  )

  grepl(
    "<html|<!doctype|access denied|forbidden|not found|bad gateway|service unavailable",
    text_head,
    ignore.case = TRUE
  )
}

file_basic_valid <- function(path, minimum_bytes) {
  file.exists(path) &&
    is.finite(file.info(path)$size) &&
    file.info(path)$size >= minimum_bytes &&
    !is_html_or_error_text(path)
}

download_once_with_curl <- function(
  curl_exe,
  url,
  destination,
  role,
  round_number
) {
  log_file <- file.path(
    LOGDIR,
    paste0(
      "download_",
      role,
      "_round",
      round_number,
      ".log"
    )
  )

  run_command(
    command = curl_exe,
    args = c(
      "--location",
      "--fail",
      "--show-error",
      "--retry",
      as.character(CURL_RETRY_COUNT),
      "--retry-all-errors",
      "--retry-delay",
      as.character(CURL_RETRY_DELAY),
      "--connect-timeout",
      as.character(CURL_CONNECT_TIMEOUT),
      "--speed-time",
      as.character(CURL_SPEED_TIME),
      "--speed-limit",
      as.character(CURL_SPEED_LIMIT),
      "--http1.1",
      "--continue-at",
      "-",
      "--output",
      quote_cmd(destination),
      quote_cmd(url)
    ),
    log_file = log_file
  )
}

download_once_with_r <- function(
  url,
  destination,
  role,
  round_number
) {
  log_file <- file.path(
    LOGDIR,
    paste0(
      "download_",
      role,
      "_Rfallback_round",
      round_number,
      ".log"
    )
  )

  start_time <- Sys.time()

  status <- tryCatch(
    {
      sink(
        log_file,
        append = FALSE,
        split = FALSE
      )

      on.exit(
        sink(),
        add = TRUE
      )

      utils::download.file(
        url = url,
        destfile = destination,
        method = "libcurl",
        mode = "wb",
        quiet = FALSE
      )
    },
    error = function(e) {
      cat(
        paste0(
          "download.file ERROR: ",
          conditionMessage(e),
          "\n"
        )
      )

      999L
    }
  )

  while (sink.number() > 0) {
    sink()
  }

  if (is.null(status)) {
    status <- 0L
  }

  elapsed <- as.numeric(
    difftime(
      Sys.time(),
      start_time,
      units = "secs"
    )
  )

  list(
    status = as.integer(status),
    elapsed_seconds = elapsed,
    log_file = log_file
  )
}

download_resumable <- function(
  url,
  destination,
  minimum_bytes,
  role
) {
  if (file_basic_valid(destination, minimum_bytes)) {
    msg(
      role,
      " already exists and passes the basic file check: ",
      format_bytes(file.info(destination)$size)
    )

    return(
      data.frame(
        role = role,
        method = "existing",
        round = 0L,
        status_code = 0L,
        file_size_bytes = file.info(destination)$size,
        file_size_human = format_bytes(file.info(destination)$size),
        basic_valid = TRUE,
        elapsed_seconds = 0,
        stringsAsFactors = FALSE
      )
    )
  }

  curl_exe <- unname(
    Sys.which("curl")
  )

  audit_rows <- list()

  for (round_number in seq_len(MAX_DOWNLOAD_ROUNDS)) {
    before_size <- if (file.exists(destination)) {
      file.info(destination)$size
    } else {
      0
    }

    msg(
      "Downloading ",
      role,
      " | round ",
      round_number,
      "/",
      MAX_DOWNLOAD_ROUNDS,
      " | current size ",
      format_bytes(before_size)
    )

    if (
      nzchar(curl_exe) &&
      file.exists(curl_exe)
    ) {
      result <- download_once_with_curl(
        curl_exe = curl_exe,
        url = url,
        destination = destination,
        role = role,
        round_number = round_number
      )

      method <- "curl_resume"
    } else {
      result <- download_once_with_r(
        url = url,
        destination = destination,
        role = role,
        round_number = round_number
      )

      method <- "R_libcurl"
    }

    after_size <- if (file.exists(destination)) {
      file.info(destination)$size
    } else {
      0
    }

    valid_now <- file_basic_valid(
      destination,
      minimum_bytes
    )

    audit_rows[[length(audit_rows) + 1L]] <- data.frame(
      role = role,
      method = method,
      round = round_number,
      status_code = result$status,
      file_size_bytes = after_size,
      file_size_human = format_bytes(after_size),
      basic_valid = valid_now,
      elapsed_seconds = result$elapsed_seconds,
      stringsAsFactors = FALSE
    )

    if (valid_now) {
      msg(
        role,
        " download completed: ",
        format_bytes(after_size)
      )

      break
    }

    msg(
      role,
      " is not complete yet. The next round will resume the partial file."
    )

    Sys.sleep(
      min(
        30L * round_number,
        120L
      )
    )
  }

  do.call(
    rbind,
    audit_rows
  )
}

# ==============================================================================
# 2. START AND INPUT CHECKS
# ==============================================================================

overall_start <- Sys.time()

msg("==============================================================")
msg("GSE162610 Step 7A-3 | ctxcore and resource targeted repair")
msg("==============================================================")
msg("The existing pyscenic_env will be reused.")
msg("This script will not run GRNBoost2.")

if (!dir.exists(PROJECT_ROOT)) {
  abort_step(
    paste0(
      "Project root does not exist: ",
      PROJECT_ROOT
    )
  )
}

for (required_path in c(
  STEP7A_SCRIPT,
  STEP7A_GENE_FILE,
  TF_FILE
)) {
  if (!file.exists(required_path)) {
    abort_step(
      paste0(
        "Required input is missing: ",
        required_path
      )
    )
  }
}

if (!file_basic_valid(TF_FILE, TF_MIN_BYTES)) {
  abort_step(
    paste0(
      "The existing mouse TF list is missing or invalid: ",
      TF_FILE
    )
  )
}

CONDA_EXE <- find_conda()

if (is.na(CONDA_EXE)) {
  abort_step(
    "Conda could not be found."
  )
}

ENV_PYTHON <- get_env_python(
  CONDA_EXE
)

if (
  is.na(ENV_PYTHON) ||
  !file.exists(ENV_PYTHON)
) {
  abort_step(
    paste0(
      "The existing conda environment '",
      CONDA_ENV_NAME,
      "' is not callable."
    )
  )
}

msg("Conda: ", safe_path(CONDA_EXE))
msg("Environment Python: ", ENV_PYTHON)

# ==============================================================================
# 3. REPAIR pkg_resources / ctxcore
# ==============================================================================

repair_start <- Sys.time()

msg("Installing a setuptools version that still provides pkg_resources...")

setuptools_result <- run_command(
  command = CONDA_EXE,
  args = c(
    "install",
    "-y",
    "-n",
    CONDA_ENV_NAME,
    "--override-channels",
    "--strict-channel-priority",
    "-c",
    "conda-forge",
    quote_cmd("setuptools<82")
  ),
  log_file = file.path(
    LOGDIR,
    "01_conda_install_setuptools_lt82.log"
  )
)

if (setuptools_result$status != 0L) {
  abort_step(
    paste0(
      "setuptools<82 installation failed.\n\n",
      read_log_tail(
        setuptools_result$log_file,
        n = 120L
      )
    )
  )
}

# Import probe.
IMPORT_PROBE_FILE <- file.path(
  OUTDIR,
  "probe_ctxcore_pyscenic.py"
)

writeLines(
  c(
    "import sys",
    "import importlib",
    "import importlib.metadata as metadata",
    "",
    "packages = ['setuptools', 'pkg_resources', 'ctxcore', 'pyscenic', 'arboreto', 'loompy']",
    "",
    "print('python_executable\\t' + sys.executable)",
    "print('python_version\\t' + sys.version.replace('\\n', ' '))",
    "",
    "for package in packages:",
    "    try:",
    "        module = importlib.import_module(package)",
    "        distribution_name = 'setuptools' if package == 'pkg_resources' else package",
    "        try:",
    "            version = metadata.version(distribution_name)",
    "        except Exception:",
    "            version = getattr(module, '__version__', 'UNKNOWN')",
    "        print(package + '\\tPASS\\t' + str(version) + '\\t')",
    "    except Exception as exc:",
    "        text = str(exc).replace('\\t', ' ').replace('\\n', ' ')",
    "        print(package + '\\tFAIL\\tNA\\t' + text)"
  ),
  con = IMPORT_PROBE_FILE,
  useBytes = TRUE
)

import_probe <- capture_command(
  CONDA_EXE,
  conda_run_args(
    "python",
    quote_cmd(IMPORT_PROBE_FILE)
  )
)

import_lines <- import_probe$output[
  grepl(
    "^(setuptools|pkg_resources|ctxcore|pyscenic|arboreto|loompy)\\t",
    import_probe$output
  )
]

if (length(import_lines) > 0) {
  split_imports <- strsplit(
    import_lines,
    "\t",
    fixed = TRUE
  )

  import_audit <- do.call(
    rbind,
    lapply(
      split_imports,
      function(xx) {
        length(xx) <- 4L

        data.frame(
          package = xx[[1]],
          import_status = xx[[2]],
          version = xx[[3]],
          error = xx[[4]],
          stringsAsFactors = FALSE
        )
      }
    )
  )

  rownames(import_audit) <- NULL
} else {
  import_audit <- data.frame(
    package = c(
      "setuptools",
      "pkg_resources",
      "ctxcore",
      "pyscenic",
      "arboreto",
      "loompy"
    ),
    import_status = "FAIL",
    version = NA_character_,
    error = paste(
      import_probe$output,
      collapse = " | "
    ),
    stringsAsFactors = FALSE
  )
}

write_utf8_csv(
  import_audit,
  file.path(
    TABDIR,
    "01_ctxcore_python_import_audit.csv"
  )
)

required_imports <- c(
  "pkg_resources",
  "ctxcore",
  "pyscenic",
  "arboreto",
  "loompy"
)

import_pass <- all(
  required_imports %in%
    import_audit$package
) &&
  all(
    import_audit$import_status[
      match(
        required_imports,
        import_audit$package
      )
    ] == "PASS"
  )

if (!import_pass) {
  abort_step(
    paste0(
      "ctxcore/pkg_resources repair did not pass.\n\n",
      paste(
        capture.output(
          print(
            import_audit,
            row.names = FALSE
          )
        ),
        collapse = "\n"
      )
    )
  )
}

msg("pkg_resources and ctxcore imports: PASS")

# pySCENIC CLI probe.
cli_probe <- capture_command(
  CONDA_EXE,
  conda_run_args(
    "pyscenic",
    "--help"
  )
)

cli_pass <- cli_probe$status == 0L &&
  any(
    grepl(
      "pyscenic|Single-Cell Regulatory Network Inference|grn|aucell",
      cli_probe$output,
      ignore.case = TRUE
    )
  )

cli_audit <- data.frame(
  command = "conda run -n pyscenic_env pyscenic --help",
  status_code = cli_probe$status,
  callable = cli_pass,
  output_preview = paste(
    utils::head(
      cli_probe$output,
      20L
    ),
    collapse = " | "
  ),
  stringsAsFactors = FALSE
)

write_utf8_csv(
  cli_audit,
  file.path(
    TABDIR,
    "02_pyscenic_CLI_audit.csv"
  )
)

if (!cli_pass) {
  abort_step(
    paste0(
      "pySCENIC CLI is still not callable.\n\n",
      paste(
        utils::tail(
          cli_probe$output,
          50L
        ),
        collapse = "\n"
      )
    )
  )
}

msg("pySCENIC CLI: PASS")

repair_elapsed <- as.numeric(
  difftime(
    Sys.time(),
    repair_start,
    units = "secs"
  )
)

# ==============================================================================
# 4. DOWNLOAD MOUSE MOTIF ANNOTATION AND ONE 10 KB FEATHER
# ==============================================================================

download_start <- Sys.time()

motif_download_audit <- download_resumable(
  url = MOTIF_URL,
  destination = MOTIF_FILE,
  minimum_bytes = MOTIF_MIN_BYTES,
  role = "mouse_motif_annotation"
)

feather_download_audit <- download_resumable(
  url = FEATHER_URL,
  destination = FEATHER_FILE,
  minimum_bytes = FEATHER_MIN_BYTES,
  role = "mm10_10kb_feather"
)

download_audit <- rbind(
  motif_download_audit,
  feather_download_audit
)

download_audit$elapsed_human <- vapply(
  download_audit$elapsed_seconds,
  format_duration,
  FUN.VALUE = character(1)
)

write_utf8_csv(
  download_audit,
  file.path(
    TABDIR,
    "03_resource_download_retry_audit.csv"
  )
)

motif_basic_pass <- file_basic_valid(
  MOTIF_FILE,
  MOTIF_MIN_BYTES
)

feather_basic_pass <- file_basic_valid(
  FEATHER_FILE,
  FEATHER_MIN_BYTES
)

if (!motif_basic_pass || !feather_basic_pass) {
  abort_step(
    paste0(
      "One or more official mouse resources remain incomplete.\n\n",
      paste(
        capture.output(
          print(
            download_audit,
            row.names = FALSE
          )
        ),
        collapse = "\n"
      ),
      "\n\nDo not delete partial files; the next repair run can resume them."
    )
  )
}

download_elapsed <- as.numeric(
  difftime(
    Sys.time(),
    download_start,
    units = "secs"
  )
)

# ==============================================================================
# 5. CONTENT AND DATABASE READABILITY VALIDATION
# ==============================================================================

msg("Validating motif annotation and feather database...")

motif_head <- readLines(
  MOTIF_FILE,
  n = 10L,
  warn = FALSE,
  encoding = "UTF-8"
)

motif_content_pass <- length(motif_head) >= 2L &&
  any(
    grepl(
      "\t",
      motif_head,
      fixed = TRUE
    )
  ) &&
  !any(
    grepl(
      "<html|<!doctype",
      motif_head,
      ignore.case = TRUE
    )
  )

FEATHER_PROBE_FILE <- file.path(
  OUTDIR,
  "probe_feather_database.py"
)

writeLines(
  c(
    "import sys",
    "from pathlib import Path",
    "from ctxcore.rnkdb import FeatherRankingDatabase",
    "",
    "gene_file = Path(sys.argv[1])",
    "db_file = Path(sys.argv[2])",
    "",
    "with gene_file.open('r', encoding='utf-8') as handle:",
    "    expression_genes = {line.strip().split('\\t')[0] for line in handle if line.strip()}",
    "",
    "try:",
    "    db = FeatherRankingDatabase(fname=str(db_file), name='mm10_10kb')",
    "    db_genes = set(db.genes)",
    "    overlap = expression_genes.intersection(db_genes)",
    "    fraction = len(overlap) / max(1, len(expression_genes))",
    "    print(",
    "        str(db_file) + '\\tPASS\\t' + str(db.total_genes) + '\\t' +",
    "        str(len(expression_genes)) + '\\t' + str(len(overlap)) + '\\t' +",
    "        str(fraction) + '\\t'",
    "    )",
    "except Exception as exc:",
    "    text = str(exc).replace('\\t', ' ').replace('\\n', ' ')",
    "    print(str(db_file) + '\\tFAIL\\tNA\\t' + str(len(expression_genes)) +",
    "          '\\tNA\\tNA\\t' + text)"
  ),
  con = FEATHER_PROBE_FILE,
  useBytes = TRUE
)

feather_probe <- capture_command(
  CONDA_EXE,
  conda_run_args(
    "python",
    quote_cmd(FEATHER_PROBE_FILE),
    quote_cmd(STEP7A_GENE_FILE),
    quote_cmd(FEATHER_FILE)
  )
)

feather_line <- feather_probe$output[
  grepl(
    "\\.feather\\t",
    feather_probe$output
  )
]

if (length(feather_line) >= 1L) {
  xx <- strsplit(
    feather_line[[1]],
    "\t",
    fixed = TRUE
  )[[1]]

  length(xx) <- 7L

  feather_audit <- data.frame(
    path = xx[[1]],
    read_status = xx[[2]],
    database_gene_count = suppressWarnings(
      as.numeric(xx[[3]])
    ),
    expression_gene_count = suppressWarnings(
      as.numeric(xx[[4]])
    ),
    overlap_gene_count = suppressWarnings(
      as.numeric(xx[[5]])
    ),
    overlap_fraction = suppressWarnings(
      as.numeric(xx[[6]])
    ),
    error = xx[[7]],
    stringsAsFactors = FALSE
  )
} else {
  feather_audit <- data.frame(
    path = safe_path(FEATHER_FILE),
    read_status = "FAIL",
    database_gene_count = NA_real_,
    expression_gene_count = NA_real_,
    overlap_gene_count = NA_real_,
    overlap_fraction = NA_real_,
    error = paste(
      feather_probe$output,
      collapse = " | "
    ),
    stringsAsFactors = FALSE
  )
}

write_utf8_csv(
  feather_audit,
  file.path(
    TABDIR,
    "04_feather_readability_gene_overlap_audit.csv"
  )
)

feather_read_pass <- nrow(feather_audit) == 1L &&
  feather_audit$read_status[[1]] == "PASS" &&
  is.finite(
    feather_audit$overlap_fraction[[1]]
  ) &&
  feather_audit$overlap_fraction[[1]] >= 0.70

resource_validation <- data.frame(
  resource = c(
    "mouse_TF_list",
    "mouse_motif_annotation",
    "mm10_10kb_feather"
  ),
  path = safe_path(
    c(
      TF_FILE,
      MOTIF_FILE,
      FEATHER_FILE
    )
  ),
  size_bytes = file.info(
    c(
      TF_FILE,
      MOTIF_FILE,
      FEATHER_FILE
    )
  )$size,
  size_human = vapply(
    file.info(
      c(
        TF_FILE,
        MOTIF_FILE,
        FEATHER_FILE
      )
    )$size,
    format_bytes,
    FUN.VALUE = character(1)
  ),
  basic_pass = c(
    file_basic_valid(
      TF_FILE,
      TF_MIN_BYTES
    ),
    motif_basic_pass,
    feather_basic_pass
  ),
  content_or_readability_pass = c(
    TRUE,
    motif_content_pass,
    feather_read_pass
  ),
  final_pass = c(
    file_basic_valid(
      TF_FILE,
      TF_MIN_BYTES
    ),
    motif_basic_pass &&
      motif_content_pass,
    feather_basic_pass &&
      feather_read_pass
  ),
  stringsAsFactors = FALSE
)

write_utf8_csv(
  resource_validation,
  file.path(
    TABDIR,
    "05_final_resource_validation.csv"
  )
)

if (!all(resource_validation$final_pass)) {
  abort_step(
    paste0(
      "Resource validation failed.\n\n",
      paste(
        capture.output(
          print(
            resource_validation,
            row.names = FALSE
          )
        ),
        collapse = "\n"
      )
    )
  )
}

msg("Mouse TF, motif annotation and 10 kb feather validation: PASS")

# ==============================================================================
# 6. AUTOMATICALLY RE-RUN STEP 7A
# ==============================================================================

msg("Re-running Step 7A in a clean R process...")

rscript_exe <- file.path(
  R.home("bin"),
  "Rscript.exe"
)

if (!file.exists(rscript_exe)) {
  rscript_exe <- file.path(
    R.home("bin"),
    "Rscript"
  )
}

if (!file.exists(rscript_exe)) {
  abort_step(
    "Rscript executable could not be located."
  )
}

rerun_result <- run_command(
  command = rscript_exe,
  args = quote_cmd(STEP7A_SCRIPT),
  log_file = file.path(
    LOGDIR,
    "05_rerun_step7A.log"
  )
)

if (rerun_result$status != 0L) {
  abort_step(
    paste0(
      "Step 7A re-run failed.\n\n",
      read_log_tail(
        rerun_result$log_file,
        n = 120L
      )
    )
  )
}

step7a_status <- NA_character_

if (file.exists(STEP7A_REPORT)) {
  step7a_lines <- readLines(
    STEP7A_REPORT,
    warn = FALSE,
    encoding = "UTF-8"
  )

  status_line <- step7a_lines[
    grepl(
      "^Overall status:",
      step7a_lines
    )
  ]

  if (length(status_line) > 0) {
    step7a_status <- trimws(
      sub(
        "^Overall status:",
        "",
        status_line[[1]]
      )
    )
  }
}

step7a_pass <- identical(
  step7a_status,
  "PASS"
)

# ==============================================================================
# 7. FINAL DECISION AND REPORT
# ==============================================================================

overall_elapsed <- as.numeric(
  difftime(
    Sys.time(),
    overall_start,
    units = "secs"
  )
)

final_decision <- data.frame(
  domain = c(
    "Existing conda environment callable",
    "setuptools provides pkg_resources",
    "ctxcore import",
    "pySCENIC import",
    "pySCENIC CLI",
    "Mouse TF list",
    "Mouse motif annotation",
    "mm10 10 kb feather basic file check",
    "mm10 feather readable by ctxcore",
    "mm10 feather expression-gene overlap >=70%",
    "Step 7A re-run",
    "Updated Step 7A overall status",
    "Overall Step 7A-3"
  ),
  pass = c(
    file.exists(ENV_PYTHON),
    import_audit$import_status[
      match(
        "pkg_resources",
        import_audit$package
      )
    ] == "PASS",
    import_audit$import_status[
      match(
        "ctxcore",
        import_audit$package
      )
    ] == "PASS",
    import_audit$import_status[
      match(
        "pyscenic",
        import_audit$package
      )
    ] == "PASS",
    cli_pass,
    resource_validation$final_pass[
      resource_validation$resource ==
        "mouse_TF_list"
    ],
    resource_validation$final_pass[
      resource_validation$resource ==
        "mouse_motif_annotation"
    ],
    feather_basic_pass,
    feather_audit$read_status[[1]] == "PASS",
    feather_audit$overlap_fraction[[1]] >= 0.70,
    rerun_result$status == 0L,
    step7a_pass,
    step7a_pass
  ),
  stringsAsFactors = FALSE
)

final_decision$status <- ifelse(
  final_decision$pass,
  "PASS",
  "FAIL"
)

write_utf8_csv(
  final_decision,
  file.path(
    TABDIR,
    "06_final_decision.csv"
  )
)

FINAL_STATUS <- if (
  all(
    final_decision$pass[
      final_decision$domain !=
        "Overall Step 7A-3"
    ]
  )
) {
  "PASS"
} else {
  "FAIL"
}

report_lines <- c(
  "======================================================================",
  "GSE162610 | STEP 7A-3 REPORT",
  "ctxcore/pkg_resources and mouse-resource targeted repair",
  "======================================================================",
  "",
  paste0("Generated: ", timestamp()),
  paste0("Project root: ", PROJECT_ROOT),
  paste0("Existing environment: ", CONDA_ENV_NAME),
  paste0("Environment Python: ", ENV_PYTHON),
  paste0("Overall Step 7A-3 status: ", FINAL_STATUS),
  "",
  "SCOPE:",
  "This script reused the existing Bioconda pyscenic_env.",
  "It did not delete/recreate the environment and did not run GRNBoost2.",
  "",
  "TIMING:",
  paste0(
    "setuptools/ctxcore repair: ",
    format_duration(repair_elapsed)
  ),
  paste0(
    "resource download/retry: ",
    format_duration(download_elapsed)
  ),
  paste0(
    "Step7A re-run: ",
    format_duration(rerun_result$elapsed_seconds)
  ),
  paste0(
    "total elapsed: ",
    format_duration(overall_elapsed)
  ),
  "",
  "PYTHON IMPORT AUDIT:",
  paste(
    capture.output(
      print(
        import_audit,
        row.names = FALSE
      )
    ),
    collapse = "\n"
  ),
  "",
  "PYSCENIC CLI AUDIT:",
  paste(
    capture.output(
      print(
        cli_audit,
        row.names = FALSE
      )
    ),
    collapse = "\n"
  ),
  "",
  "RESOURCE DOWNLOAD AUDIT:",
  paste(
    capture.output(
      print(
        download_audit,
        row.names = FALSE
      )
    ),
    collapse = "\n"
  ),
  "",
  "FEATHER READABILITY / GENE OVERLAP:",
  paste(
    capture.output(
      print(
        feather_audit,
        row.names = FALSE
      )
    ),
    collapse = "\n"
  ),
  "",
  "FINAL RESOURCE VALIDATION:",
  paste(
    capture.output(
      print(
        resource_validation,
        row.names = FALSE
      )
    ),
    collapse = "\n"
  ),
  "",
  "FINAL DECISION:",
  paste(
    capture.output(
      print(
        final_decision,
        row.names = FALSE
      )
    ),
    collapse = "\n"
  ),
  "",
  paste0(
    "Updated Step 7A status: ",
    ifelse(
      is.na(step7a_status),
      "<not found>",
      step7a_status
    )
  ),
  "",
  "NEXT ACTION:",
  if (FINAL_STATUS == "PASS") {
    paste0(
      "Step 7A environment/input audit now passes. ",
      "The next formal step is Step 7B GRNBoost2."
    )
  } else {
    paste0(
      "Do not run GRNBoost2. Upload the Step7A3 report and the failed audit table."
    )
  },
  "",
  "UPLOAD THESE FILES (maximum 5):",
  "1) reports/Step7A3_report.txt",
  "2) tables/01_ctxcore_python_import_audit.csv",
  "3) tables/03_resource_download_retry_audit.csv",
  "4) tables/04_feather_readability_gene_overlap_audit.csv",
  "5) the updated Step7A report",
  "",
  "Do not upload the feather database."
)

writeLines(
  report_lines,
  con = REPORT_FILE,
  useBytes = TRUE
)

writeLines(
  capture.output(
    sessionInfo()
  ),
  con = SESSION_FILE,
  useBytes = TRUE
)

msg("==============================================================")
msg("STEP 7A-3 COMPLETE")
msg("==============================================================")
msg("Step 7A-3 status: ", FINAL_STATUS)
msg(
  "Updated Step 7A status: ",
  ifelse(
    is.na(step7a_status),
    "<not found>",
    step7a_status
  )
)
msg("Total elapsed: ", format_duration(overall_elapsed))
msg("Report: ", REPORT_FILE)
msg("")
msg("请上传：")
msg("1) reports/Step7A3_report.txt")
msg("2) tables/01_ctxcore_python_import_audit.csv")
msg("3) tables/03_resource_download_retry_audit.csv")
msg("4) tables/04_feather_readability_gene_overlap_audit.csv")
msg("5) 更新后的 Step7A_report.txt")
