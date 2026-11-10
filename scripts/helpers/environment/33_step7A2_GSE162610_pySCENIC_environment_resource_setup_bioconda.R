# ==============================================================================
# 33_step7A2_GSE162610_pySCENIC_environment_resource_setup_v3_bioconda.R
# ==============================================================================
# Purpose:
#   Rebuild an isolated Windows conda environment using the Bioconda pySCENIC 0.12.1 package,
#   download the official Mus musculus mm10/mc9nr resources, validate them,
#   and then re-run Step 7A in a separate R process.
#
# Run in RStudio:
#   source("<PROJECT_ROOT>/33_step7A2_GSE162610_pySCENIC_environment_resource_setup_v3_bioconda.R")
#
# This script DOES NOT run GRNBoost2, motif pruning, AUCell, RSS, or CSI.
# Large downloads are resumable when Windows curl.exe is available.
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
# 0. USER-EDITABLE SETTINGS
# ==============================================================================

PROJECT_ROOT <- Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
CONDA_ENV_NAME <- "pyscenic_env"

# FALSE means: reuse/repair an existing environment if present.
# TRUE means: delete and rebuild the environment from scratch.
FORCE_RECREATE_ENV <- TRUE

# Two official mm10/mc9nr gene-based databases are downloaded:
#   10 kb around TSS
#   500 bp upstream / 100 bp downstream
DOWNLOAD_SECOND_FEATHER <- FALSE

# Automatically re-run Step 7A after setup passes.
AUTO_RERUN_STEP7A <- TRUE

# Official package index.
PIP_INDEX_URL <- "https://pypi.org/simple"

# ==============================================================================
# 1. PATHS
# ==============================================================================

OUTDIR <- file.path(
  PROJECT_ROOT,
  "04_results",
  "step7A2_GSE162610_pySCENIC_environment_resource_setup"
)

LOGDIR <- file.path(OUTDIR, "logs")
TABDIR <- file.path(OUTDIR, "tables")
REPORTDIR <- file.path(OUTDIR, "reports")

RESOURCE_DIR <- file.path(
  PROJECT_ROOT,
  "01_raw_data",
  "pySCENIC_resources"
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

for (dd in c(OUTDIR, LOGDIR, TABDIR, REPORTDIR, RESOURCE_DIR)) {
  dir.create(
    dd,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

REPORT_FILE <- file.path(
  REPORTDIR,
  "Step7A2_report.txt"
)

ERROR_FILE <- file.path(
  REPORTDIR,
  "Step7A2_ERROR.txt"
)

SESSION_FILE <- file.path(
  REPORTDIR,
  "Step7A2_sessionInfo.txt"
)

# ==============================================================================
# 2. OFFICIAL RESOURCES
# ==============================================================================

URLS <- c(
  mouse_TF_list =
    "https://raw.githubusercontent.com/aertslab/pySCENIC/master/resources/mm_mgi_tfs.txt",

  mouse_motif_annotation =
    "https://resources.aertslab.org/cistarget/motif2tf/motifs-v9-nr.mgi-m0.001-o0.0.tbl",

  mm10_10kb_feather =
    paste0(
      "https://resources.aertslab.org/cistarget/databases/",
      "mus_musculus/mm10/refseq_r80/mc9nr/gene_based/",
      "mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.",
      "genes_vs_motifs.rankings.feather"
    ),

  mm10_500bp_feather =
    paste0(
      "https://resources.aertslab.org/cistarget/databases/",
      "mus_musculus/mm10/refseq_r80/mc9nr/gene_based/",
      "mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.",
      "genes_vs_motifs.rankings.feather"
    ),

  database_sha256 =
    "https://resources.aertslab.org/cistarget/databases/sha256sum.txt"
)

if (!DOWNLOAD_SECOND_FEATHER) {
  URLS <- URLS[
    names(URLS) != "mm10_500bp_feather"
  ]
}

RESOURCE_PATHS <- file.path(
  RESOURCE_DIR,
  basename(URLS)
)

names(RESOURCE_PATHS) <- names(URLS)

# Minimum plausible sizes. These are deliberately below the official sizes,
# so partial files and HTML error pages are rejected without requiring an
# exact byte count.
MIN_SIZE_BYTES <- c(
  mouse_TF_list = 5000,
  mouse_motif_annotation = 80 * 1024^2,
  mm10_10kb_feather = 800 * 1024^2,
  mm10_500bp_feather = 800 * 1024^2,
  database_sha256 = 10000
)

MIN_SIZE_BYTES <- MIN_SIZE_BYTES[
  names(RESOURCE_PATHS)
]

# ==============================================================================
# 3. REPRODUCIBLE BIOCONDA INSTALL SPECIFICATION
# ==============================================================================

INSTALL_SPEC_FILE <- file.path(
  OUTDIR,
  "pyscenic_bioconda_install_spec.txt"
)

install_spec <- c(
  "Channels, in strict priority order:",
  "  1. conda-forge",
  "  2. bioconda",
  "",
  "Environment:",
  paste0("  name: ", CONDA_ENV_NAME),
  "  python=3.10",
  "  numpy<1.24",
  "  pyscenic=0.12.1",
  "",
  "Reason:",
  "  Use the Bioconda pySCENIC package to avoid PyPI TLS/SSL and dependency",
  "  backtracking failures observed in the previous pip-based installation."
)

writeLines(
  install_spec,
  con = INSTALL_SPEC_FILE,
  useBytes = TRUE
)

# ==============================================================================
# 4. UTILITIES
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

quote_arg <- function(x) {
  shQuote(
    as.character(x),
    type = "cmd"
  )
}

read_log_tail <- function(path, n = 40L) {
  if (!file.exists(path)) {
    return("<log not created>")
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
      "GSE162610 | STEP 7A-2 ERROR",
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
  dir.create(
    dirname(log_file),
    recursive = TRUE,
    showWarnings = FALSE
  )

  if (!append && file.exists(log_file)) {
    file.remove(log_file)
  }

  if (echo_command) {
    command_text <- paste(
      c(
        safe_path(command),
        args
      ),
      collapse = " "
    )

    cat(
      paste0(
        "\n[",
        timestamp(),
        "] COMMAND\n",
        command_text,
        "\n"
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
    n = 4096
  )

  if (length(raw_head) == 0) {
    return(TRUE)
  }

  text_head <- tryCatch(
    rawToChar(raw_head),
    error = function(e) ""
  )

  grepl(
    "<html|<!doctype|access denied|forbidden|not found|bad gateway",
    text_head,
    ignore.case = TRUE
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
      file.path(
        Sys.getenv("LOCALAPPDATA"),
        "anaconda3",
        "Scripts",
        "conda.exe"
      ),
      file.path(
        Sys.getenv("LOCALAPPDATA"),
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

conda_args <- function(...) {
  c(...)
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
      quote_arg(
        "import sys; print(sys.executable)"
      )
    )
  )

  if (probe$status != 0L) {
    return(NA_character_)
  }

  lines <- trimws(probe$output)
  lines <- lines[nzchar(lines)]
  lines <- lines[file.exists(lines)]

  if (length(lines) == 0) {
    return(NA_character_)
  }

  safe_path(lines[[length(lines)]])
}

download_with_resume <- function(
  url,
  dest,
  minimum_size,
  role
) {
  start_time <- Sys.time()
  before_size <- if (file.exists(dest)) {
    file.info(dest)$size
  } else {
    0
  }

  if (
    file.exists(dest) &&
    is.finite(file.info(dest)$size) &&
    file.info(dest)$size >= minimum_size &&
    !is_html_or_error_text(dest)
  ) {
    msg(
      "Resource already present; skipping download: ",
      basename(dest)
    )

    return(
      data.frame(
        role = role,
        url = url,
        path = safe_path(dest),
        action = "SKIPPED_VALID_EXISTING",
        download_status = 0L,
        before_size_bytes = before_size,
        after_size_bytes = file.info(dest)$size,
        elapsed_seconds = 0,
        stringsAsFactors = FALSE
      )
    )
  }

  msg("Downloading: ", basename(dest))
  msg("Large files are resumable; do not close RStudio.")

  curl_exe <- unname(
    Sys.which("curl")
  )

  log_file <- file.path(
    LOGDIR,
    paste0(
      "download_",
      role,
      ".log"
    )
  )

  if (
    nzchar(curl_exe) &&
    file.exists(curl_exe)
  ) {
    result <- run_command(
      command = curl_exe,
      args = c(
        "-L",
        "--fail",
        "--retry",
        "6",
        "--retry-all-errors",
        "--retry-delay",
        "10",
        "--connect-timeout",
        "60",
        "--continue-at",
        "-",
        "--output",
        quote_arg(dest),
        quote_arg(url)
      ),
      log_file = log_file
    )

    status <- result$status
  } else {
    status <- tryCatch(
      {
        utils::download.file(
          url = url,
          destfile = dest,
          mode = "wb",
          method = "libcurl",
          quiet = FALSE
        )
      },
      error = function(e) {
        cat(
          paste0(
            "download.file ERROR: ",
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
  }

  elapsed <- as.numeric(
    difftime(
      Sys.time(),
      start_time,
      units = "secs"
    )
  )

  after_size <- if (file.exists(dest)) {
    file.info(dest)$size
  } else {
    0
  }

  data.frame(
    role = role,
    url = url,
    path = safe_path(dest),
    action = ifelse(
      before_size > 0,
      "RESUME_OR_REDOWNLOAD",
      "NEW_DOWNLOAD"
    ),
    download_status = as.integer(status),
    before_size_bytes = before_size,
    after_size_bytes = after_size,
    elapsed_seconds = elapsed,
    stringsAsFactors = FALSE
  )
}

# ==============================================================================
# 5. START
# ==============================================================================

overall_start <- Sys.time()

msg("==============================================================")
msg("GSE162610 Step 7A-2 | pySCENIC environment/resource setup")
msg("==============================================================")
msg("This step will not run GRNBoost2.")

if (!dir.exists(PROJECT_ROOT)) {
  abort_step(
    paste0(
      "Project root does not exist: ",
      PROJECT_ROOT
    )
  )
}

if (!file.exists(STEP7A_SCRIPT)) {
  abort_step(
    paste0(
      "Step 7A script is missing: ",
      STEP7A_SCRIPT
    )
  )
}

# ==============================================================================
# 6. DETECT CONDA
# ==============================================================================

CONDA_EXE <- find_conda()

if (is.na(CONDA_EXE)) {
  abort_step(
    paste0(
      "Conda was not found. Expected locations included ",
      "<USERPROFILE>/anaconda3/Scripts/conda.exe. ",
      "Upload Step7A2_ERROR.txt; do not install packages into R manually."
    )
  )
}

msg("Conda detected: ", safe_path(CONDA_EXE))

conda_version_probe <- capture_command(
  CONDA_EXE,
  "--version"
)

conda_audit <- data.frame(
  conda_executable = safe_path(CONDA_EXE),
  conda_version_status = conda_version_probe$status,
  conda_version_output = paste(
    conda_version_probe$output,
    collapse = " | "
  ),
  environment_name = CONDA_ENV_NAME,
  force_recreate = FORCE_RECREATE_ENV,
  stringsAsFactors = FALSE
)

write_utf8_csv(
  conda_audit,
  file.path(
    TABDIR,
    "01_conda_detection_audit.csv"
  )
)

# ==============================================================================
# 7. CREATE OR REPAIR ISOLATED ENVIRONMENT
# ==============================================================================

env_stage_start <- Sys.time()

if (FORCE_RECREATE_ENV) {
  msg("FORCE_RECREATE_ENV=TRUE: removing existing environment if present.")

  run_command(
    command = CONDA_EXE,
    args = conda_args(
      "env",
      "remove",
      "-y",
      "-n",
      CONDA_ENV_NAME
    ),
    log_file = file.path(
      LOGDIR,
      "01_conda_remove_environment.log"
    )
  )
}

ENV_PYTHON <- get_env_python(
  CONDA_EXE
)

if (is.na(ENV_PYTHON)) {
  msg(
    "Creating isolated conda environment '",
    CONDA_ENV_NAME,
    "' with Python 3.10..."
  )

  create_result <- run_command(
    command = CONDA_EXE,
    args = conda_args(
      "create",
      "-y",
      "-n",
      CONDA_ENV_NAME,
      "--override-channels",
      "-c",
      "conda-forge",
      "python=3.10",
      "pip"
    ),
    log_file = file.path(
      LOGDIR,
      "02_conda_create_environment.log"
    )
  )

  if (create_result$status != 0L) {
    abort_step(
      paste0(
        "Conda environment creation failed.\n\n",
        read_log_tail(
          create_result$log_file
        )
      )
    )
  }

  ENV_PYTHON <- get_env_python(
    CONDA_EXE
  )
}

if (
  is.na(ENV_PYTHON) ||
  !file.exists(ENV_PYTHON)
) {
  abort_step(
    paste0(
      "The conda environment exists or was created, but its Python ",
      "executable could not be resolved."
    )
  )
}

msg("Environment Python: ", ENV_PYTHON)

# Install the complete tested package stack from conda-forge + Bioconda.
# This avoids PyPI package-index TLS failures and pip dependency backtracking.
msg("Installing pySCENIC 0.12.1 from Bioconda...")

bioconda_install_result <- run_command(
  command = CONDA_EXE,
  args = conda_args(
    "install",
    "-y",
    "-n",
    CONDA_ENV_NAME,
    "--override-channels",
    "--strict-channel-priority",
    "-c",
    "conda-forge",
    "-c",
    "bioconda",
    "python=3.10",
    quote_arg("numpy<1.24"),
    "pyscenic=0.12.1"
  ),
  log_file = file.path(
    LOGDIR,
    "03_conda_install_pyscenic_bioconda.log"
  )
)

if (bioconda_install_result$status != 0L) {
  abort_step(
    paste0(
      "Bioconda pySCENIC installation failed.\n\n",
      read_log_tail(
        bioconda_install_result$log_file,
        n = 120L
      )
    )
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
      "Bioconda installation completed, but the environment Python ",
      "executable could not be resolved."
    )
  )
}

msg("Environment Python after Bioconda install: ", ENV_PYTHON)

# Record the exact solved environment for reproducibility.
explicit_result <- run_command(
  command = CONDA_EXE,
  args = conda_args(
    "list",
    "-n",
    CONDA_ENV_NAME,
    "--explicit"
  ),
  log_file = file.path(
    LOGDIR,
    "04_conda_environment_explicit.txt"
  ),
  echo_command = FALSE
)

list_result <- run_command(
  command = CONDA_EXE,
  args = conda_args(
    "list",
    "-n",
    CONDA_ENV_NAME
  ),
  log_file = file.path(
    LOGDIR,
    "04b_conda_environment_package_list.txt"
  ),
  echo_command = FALSE
)

env_elapsed <- as.numeric(
  difftime(
    Sys.time(),
    env_stage_start,
    units = "secs"
  )
)

# ==============================================================================
# 8. PYTHON PACKAGE AND CLI AUDIT
# ==============================================================================

PYTHON_PROBE <- file.path(
  OUTDIR,
  "probe_pyscenic_environment.py"
)

writeLines(
  c(
    "import sys",
    "import importlib",
    "import importlib.metadata as metadata",
    "",
    "packages = [",
    "    'pyscenic',",
    "    'arboreto',",
    "    'ctxcore',",
    "    'numpy',",
    "    'pandas',",
    "    'scipy',",
    "    'sklearn',",
    "    'numba',",
    "    'dask',",
    "    'distributed',",
    "    'pyarrow',",
    "    'loompy',",
    "    'networkx'",
    "]",
    "",
    "print('PYTHON_EXECUTABLE\\t' + sys.executable)",
    "print('PYTHON_VERSION\\t' + sys.version.replace('\\n', ' '))",
    "",
    "for package in packages:",
    "    try:",
    "        module = importlib.import_module(package)",
    "        distribution_name = 'scikit-learn' if package == 'sklearn' else package",
    "        try:",
    "            version = metadata.version(distribution_name)",
    "        except Exception:",
    "            version = getattr(module, '__version__', 'UNKNOWN')",
    "        print(package + '\\tPASS\\t' + str(version) + '\\t')",
    "    except Exception as exc:",
    "        error = str(exc).replace('\\t', ' ').replace('\\n', ' ')",
    "        print(package + '\\tFAIL\\tNA\\t' + error)"
  ),
  con = PYTHON_PROBE,
  useBytes = TRUE
)

package_probe <- capture_command(
  CONDA_EXE,
  conda_run_args(
    "python",
    quote_arg(PYTHON_PROBE)
  )
)

probe_lines <- package_probe$output[
  grepl(
    "^(pyscenic|arboreto|ctxcore|numpy|pandas|scipy|sklearn|numba|dask|distributed|pyarrow|loompy|networkx)\\t",
    package_probe$output
  )
]

if (length(probe_lines) > 0) {
  split_probe <- strsplit(
    probe_lines,
    "\t",
    fixed = TRUE
  )

  python_package_audit <- do.call(
    rbind,
    lapply(
      split_probe,
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

  rownames(python_package_audit) <- NULL
} else {
  python_package_audit <- data.frame(
    package = requirements,
    import_status = "FAIL",
    version = NA_character_,
    error = paste(
      package_probe$output,
      collapse = " | "
    ),
    stringsAsFactors = FALSE
  )
}

python_package_audit$environment_name <- CONDA_ENV_NAME
python_package_audit$python_executable <- ENV_PYTHON

write_utf8_csv(
  python_package_audit,
  file.path(
    TABDIR,
    "02_python_package_audit.csv"
  )
)

required_imports <- c(
  "pyscenic",
  "arboreto",
  "ctxcore",
  "numpy",
  "pandas",
  "scipy",
  "sklearn",
  "numba",
  "dask",
  "distributed",
  "pyarrow",
  "loompy",
  "networkx"
)

package_pass <- all(
  required_imports %in%
    python_package_audit$package
) &&
  all(
    python_package_audit$import_status[
      match(
        required_imports,
        python_package_audit$package
      )
    ] == "PASS"
  )

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
      "Single-Cell Regulatory Network Inference|pyscenic",
      cli_probe$output,
      ignore.case = TRUE
    )
  )

cli_audit <- data.frame(
  environment_name = CONDA_ENV_NAME,
  python_executable = ENV_PYTHON,
  command = "conda run -n pyscenic_env pyscenic --help",
  status_code = cli_probe$status,
  callable = cli_pass,
  output_first_lines = paste(
    utils::head(
      cli_probe$output,
      15L
    ),
    collapse = " | "
  ),
  stringsAsFactors = FALSE
)

write_utf8_csv(
  cli_audit,
  file.path(
    TABDIR,
    "03_pyscenic_CLI_audit.csv"
  )
)

# ==============================================================================
# 9. DOWNLOAD OFFICIAL MOUSE RESOURCES
# ==============================================================================

resource_stage_start <- Sys.time()

download_rows <- vector(
  mode = "list",
  length = length(URLS)
)

for (ii in seq_along(URLS)) {
  role <- names(URLS)[[ii]]

  download_rows[[ii]] <- download_with_resume(
    url = URLS[[ii]],
    dest = RESOURCE_PATHS[[ii]],
    minimum_size = MIN_SIZE_BYTES[[role]],
    role = role
  )
}

download_audit <- do.call(
  rbind,
  download_rows
)

rownames(download_audit) <- NULL

download_audit$after_size_human <- vapply(
  download_audit$after_size_bytes,
  format_bytes,
  FUN.VALUE = character(1)
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
    "04_resource_download_audit.csv"
  )
)

resource_elapsed <- as.numeric(
  difftime(
    Sys.time(),
    resource_stage_start,
    units = "secs"
  )
)

# ==============================================================================
# 10. BASIC RESOURCE VALIDATION
# ==============================================================================

resource_validation <- data.frame(
  role = names(RESOURCE_PATHS),
  path = safe_path(unname(RESOURCE_PATHS)),
  exists = file.exists(unname(RESOURCE_PATHS)),
  size_bytes = ifelse(
    file.exists(unname(RESOURCE_PATHS)),
    file.info(unname(RESOURCE_PATHS))$size,
    NA_real_
  ),
  minimum_size_bytes = as.numeric(
    MIN_SIZE_BYTES[
      names(RESOURCE_PATHS)
    ]
  ),
  html_or_error_text = vapply(
    unname(RESOURCE_PATHS),
    function(path) {
      if (!file.exists(path)) {
        return(TRUE)
      }

      is_html_or_error_text(path)
    },
    FUN.VALUE = logical(1)
  ),
  stringsAsFactors = FALSE
)

resource_validation$size_human <- vapply(
  resource_validation$size_bytes,
  format_bytes,
  FUN.VALUE = character(1)
)

resource_validation$basic_pass <-
  resource_validation$exists &
  resource_validation$size_bytes >=
    resource_validation$minimum_size_bytes &
  !resource_validation$html_or_error_text

# Text-resource checks.
TF_FILE <- RESOURCE_PATHS[["mouse_TF_list"]]
MOTIF_FILE <- RESOURCE_PATHS[["mouse_motif_annotation"]]

tf_lines <- if (
  file.exists(TF_FILE) &&
  !is_html_or_error_text(TF_FILE)
) {
  trimws(
    readLines(
      TF_FILE,
      warn = FALSE,
      encoding = "UTF-8"
    )
  )
} else {
  character(0)
}

tf_lines <- unique(
  tf_lines[
    nzchar(tf_lines) &
      !grepl("^#", tf_lines)
  ]
)

tf_text_pass <- length(tf_lines) >= 1000L &&
  all(
    c(
      "Egr1",
      "Nfkb1",
      "Rela",
      "Stat1",
      "Junb",
      "Mafb"
    ) %in% tf_lines
  )

motif_head <- if (
  file.exists(MOTIF_FILE) &&
  !is_html_or_error_text(MOTIF_FILE)
) {
  readLines(
    MOTIF_FILE,
    n = 5L,
    warn = FALSE,
    encoding = "UTF-8"
  )
} else {
  character(0)
}

motif_text_pass <- length(motif_head) >= 2L &&
  any(
    grepl(
      "\t",
      motif_head,
      fixed = TRUE
    )
  )

resource_validation$text_content_pass <- NA

resource_validation$text_content_pass[
  resource_validation$role ==
    "mouse_TF_list"
] <- tf_text_pass

resource_validation$text_content_pass[
  resource_validation$role ==
    "mouse_motif_annotation"
] <- motif_text_pass

# ==============================================================================
# 11. FEATHER READABILITY AND GENE-OVERLAP AUDIT
# ==============================================================================

GENE_FILE <- file.path(
  PROJECT_ROOT,
  "04_results",
  "step7A_GSE162610_pySCENIC_input_environment_audit",
  "input",
  "GSE162610_step7A_frozen_4203_genes.tsv"
)

if (!file.exists(GENE_FILE)) {
  abort_step(
    paste0(
      "The Step 7A gene file is missing: ",
      GENE_FILE
    )
  )
}

FEATHER_FILES <- unname(
  RESOURCE_PATHS[
    grepl(
      "_feather$",
      names(RESOURCE_PATHS)
    )
  ]
)

FEATHER_PROBE <- file.path(
  OUTDIR,
  "probe_mouse_feather_databases.py"
)

python_code <- c(
  "import sys",
  "from pathlib import Path",
  "from ctxcore.rnkdb import FeatherRankingDatabase",
  "",
  "gene_file = Path(sys.argv[1])",
  "db_files = [Path(x) for x in sys.argv[2:]]",
  "",
  "with gene_file.open('r', encoding='utf-8') as handle:",
  "    expression_genes = {line.strip().split('\\t')[0] for line in handle if line.strip()}",
  "",
  "for idx, db_file in enumerate(db_files, start=1):",
  "    try:",
  "        db = FeatherRankingDatabase(fname=str(db_file), name=f'db_{idx}')",
  "        db_genes = set(db.genes)",
  "        overlap = expression_genes.intersection(db_genes)",
  "        fraction = len(overlap) / max(1, len(expression_genes))",
  "        print(",
  "            str(db_file) + '\\tPASS\\t' + str(db.total_genes) + '\\t' +",
  "            str(len(expression_genes)) + '\\t' + str(len(overlap)) + '\\t' +",
  "            str(fraction) + '\\t'",
  "        )",
  "    except Exception as exc:",
  "        text = str(exc).replace('\\t', ' ').replace('\\n', ' ')",
  "        print(str(db_file) + '\\tFAIL\\tNA\\t' + str(len(expression_genes)) +",
  "              '\\tNA\\tNA\\t' + text)"
)

writeLines(
  python_code,
  con = FEATHER_PROBE,
  useBytes = TRUE
)

feather_probe <- capture_command(
  CONDA_EXE,
  c(
    conda_run_args(
      "python",
      quote_arg(FEATHER_PROBE),
      quote_arg(GENE_FILE)
    ),
    vapply(
      FEATHER_FILES,
      quote_arg,
      FUN.VALUE = character(1)
    )
  )
)

feather_lines <- feather_probe$output[
  grepl(
    "\\.feather\\t",
    feather_probe$output
  )
]

if (length(feather_lines) > 0) {
  split_feather <- strsplit(
    feather_lines,
    "\t",
    fixed = TRUE
  )

  feather_audit <- do.call(
    rbind,
    lapply(
      split_feather,
      function(xx) {
        length(xx) <- 7L

        data.frame(
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
      }
    )
  )

  rownames(feather_audit) <- NULL
} else {
  feather_audit <- data.frame(
    path = safe_path(FEATHER_FILES),
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
    "05_feather_readability_gene_overlap_audit.csv"
  )
)

feather_pass <- nrow(feather_audit) ==
  length(FEATHER_FILES) &&
  all(
    feather_audit$read_status == "PASS"
  ) &&
  all(
    feather_audit$overlap_fraction >= 0.70
  )

# ==============================================================================
# 12. TF OVERLAP WITH THE 7,491-GENE INPUT
# ==============================================================================

expression_genes <- trimws(
  readLines(
    GENE_FILE,
    warn = FALSE,
    encoding = "UTF-8"
  )
)

expression_genes <- unique(
  expression_genes[
    nzchar(expression_genes)
  ]
)

tf_overlap <- intersect(
  expression_genes,
  tf_lines
)

tf_overlap_audit <- data.frame(
  metric = c(
    "Mouse TF-list entries",
    "Step7A retained expression genes",
    "TFs present in Step7A retained genes",
    "TF overlap fraction of TF list",
    "Egr1 retained",
    "Junb retained",
    "Nfkb1 retained",
    "Rela retained",
    "Stat1 retained",
    "Mafb retained"
  ),
  value = c(
    length(tf_lines),
    length(expression_genes),
    length(tf_overlap),
    ifelse(
      length(tf_lines) > 0,
      length(tf_overlap) /
        length(tf_lines),
      NA_real_
    ),
    "Egr1" %in% tf_overlap,
    "Junb" %in% tf_overlap,
    "Nfkb1" %in% tf_overlap,
    "Rela" %in% tf_overlap,
    "Stat1" %in% tf_overlap,
    "Mafb" %in% tf_overlap
  ),
  stringsAsFactors = FALSE
)

write_utf8_csv(
  tf_overlap_audit,
  file.path(
    TABDIR,
    "06_mouse_TF_expression_overlap_audit.csv"
  )
)

tf_overlap_pass <- tf_text_pass &&
  length(tf_overlap) >= 100L

# ==============================================================================
# 13. OPTIONAL SHA256 VERIFICATION FOR FEATHER DATABASES
# ==============================================================================

SHA_FILE <- RESOURCE_PATHS[["database_sha256"]]

sha256_audit <- data.frame(
  file = basename(FEATHER_FILES),
  expected_sha256 = NA_character_,
  observed_sha256 = NA_character_,
  checksum_available = FALSE,
  checksum_match = NA,
  stringsAsFactors = FALSE
)

if (
  file.exists(SHA_FILE) &&
  !is_html_or_error_text(SHA_FILE)
) {
  sha_lines <- readLines(
    SHA_FILE,
    warn = FALSE,
    encoding = "UTF-8"
  )

  SHA_PROBE <- file.path(
    OUTDIR,
    "calculate_sha256.py"
  )

  writeLines(
    c(
      "import hashlib",
      "import sys",
      "from pathlib import Path",
      "",
      "for item in sys.argv[1:]:",
      "    path = Path(item)",
      "    digest = hashlib.sha256()",
      "    with path.open('rb') as handle:",
      "        while True:",
      "            block = handle.read(16 * 1024 * 1024)",
      "            if not block:",
      "                break",
      "            digest.update(block)",
      "    print(path.name + '\\t' + digest.hexdigest())"
    ),
    con = SHA_PROBE,
    useBytes = TRUE
  )

  msg("Computing SHA256 checksums for the large feather files...")

  sha_probe <- capture_command(
    CONDA_EXE,
    c(
      conda_run_args(
        "python",
        quote_arg(SHA_PROBE)
      ),
      vapply(
        FEATHER_FILES,
        quote_arg,
        FUN.VALUE = character(1)
      )
    )
  )

  observed_lines <- sha_probe$output[
    grepl(
      "\\.feather\\t[0-9a-fA-F]{64}$",
      sha_probe$output
    )
  ]

  observed_map <- setNames(
    sub(
      "^.*\\t",
      "",
      observed_lines
    ),
    sub(
      "\\t.*$",
      "",
      observed_lines
    )
  )

  for (ii in seq_len(nrow(sha256_audit))) {
    bn <- sha256_audit$file[[ii]]

    matching <- sha_lines[
      grepl(
        paste0(
          "(^|/)",
          gsub(
            "\\.",
            "\\\\.",
            bn
          ),
          "$"
        ),
        sha_lines
      )
    ]

    if (length(matching) > 0) {
      expected <- sub(
        "^([0-9a-fA-F]{64}).*$",
        "\\1",
        matching[[1]]
      )

      sha256_audit$expected_sha256[[ii]] <- tolower(expected)
      sha256_audit$checksum_available[[ii]] <- TRUE
    }

    if (bn %in% names(observed_map)) {
      sha256_audit$observed_sha256[[ii]] <- tolower(
        observed_map[[bn]]
      )
    }

    if (
      sha256_audit$checksum_available[[ii]] &&
      !is.na(
        sha256_audit$observed_sha256[[ii]]
      )
    ) {
      sha256_audit$checksum_match[[ii]] <-
        sha256_audit$expected_sha256[[ii]] ==
        sha256_audit$observed_sha256[[ii]]
    }
  }
}

write_utf8_csv(
  sha256_audit,
  file.path(
    TABDIR,
    "07_feather_SHA256_audit.csv"
  )
)

# Checksum is a strengthening audit. If the official list contains the file,
# a mismatch is a failure. If no matching hash is listed, readability + size
# checks remain decisive.
checksum_pass <- !any(
  sha256_audit$checksum_available &
    !is.na(
      sha256_audit$checksum_match
    ) &
    !sha256_audit$checksum_match
)

resource_validation$final_pass <-
  resource_validation$basic_pass

resource_validation$final_pass[
  resource_validation$role ==
    "mouse_TF_list"
] <- resource_validation$basic_pass[
  resource_validation$role ==
    "mouse_TF_list"
] &&
  tf_text_pass

resource_validation$final_pass[
  resource_validation$role ==
    "mouse_motif_annotation"
] <- resource_validation$basic_pass[
  resource_validation$role ==
    "mouse_motif_annotation"
] &&
  motif_text_pass

resource_validation$final_pass[
  grepl(
    "_feather$",
    resource_validation$role
  )
] <- resource_validation$basic_pass[
  grepl(
    "_feather$",
    resource_validation$role
  )
] &&
  feather_pass &&
  checksum_pass

write_utf8_csv(
  resource_validation,
  file.path(
    TABDIR,
    "08_final_resource_validation.csv"
  )
)

# The checksum index is optional strengthening metadata. The required resources
# are the TF list, motif annotation, and the selected feather databases.
resources_pass <- all(
  resource_validation$final_pass[
    resource_validation$role !=
      "database_sha256"
  ]
)

# ==============================================================================
# 14. FINAL SETUP DECISION
# ==============================================================================

setup_decision <- data.frame(
  domain = c(
    "Conda detected",
    "Isolated Python 3.10 environment resolved",
    "Required Python package imports",
    "pySCENIC CLI callable",
    "Mouse TF-list content",
    "Mouse TF overlap with Step7A genes",
    "Mouse motif annotation",
    "mm10 feather readability and >=70% gene overlap",
    "Feather checksum audit",
    "All required resources",
    "Overall Step7A-2 setup"
  ),
  pass = c(
    file.exists(CONDA_EXE),
    file.exists(ENV_PYTHON),
    package_pass,
    cli_pass,
    tf_text_pass,
    tf_overlap_pass,
    motif_text_pass,
    feather_pass,
    checksum_pass,
    resources_pass,
    package_pass &&
      cli_pass &&
      tf_overlap_pass &&
      resources_pass
  ),
  stringsAsFactors = FALSE
)

setup_decision$status <- ifelse(
  setup_decision$pass,
  "PASS",
  "FAIL"
)

SETUP_PASS <- all(
  setup_decision$pass[
    setup_decision$domain !=
      "Overall Step7A-2 setup"
  ]
)

setup_decision$pass[
  setup_decision$domain ==
    "Overall Step7A-2 setup"
] <- SETUP_PASS

setup_decision$status[
  setup_decision$domain ==
    "Overall Step7A-2 setup"
] <- ifelse(
  SETUP_PASS,
  "PASS",
  "FAIL"
)

write_utf8_csv(
  setup_decision,
  file.path(
    TABDIR,
    "09_final_setup_decision.csv"
  )
)

# ==============================================================================
# 15. RE-RUN STEP 7A IN A SEPARATE R PROCESS
# ==============================================================================

rerun_status <- NA_integer_
rerun_elapsed <- NA_real_
rerun_log <- file.path(
  LOGDIR,
  "05_rerun_step7A.log"
)

if (
  SETUP_PASS &&
  AUTO_RERUN_STEP7A
) {
  msg("Setup passed. Re-running Step 7A in a separate R process...")

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

  rerun_result <- run_command(
    command = rscript_exe,
    args = quote_arg(STEP7A_SCRIPT),
    log_file = rerun_log
  )

  rerun_status <- rerun_result$status
  rerun_elapsed <- rerun_result$elapsed_seconds
}

step7a_final_status <- NA_character_

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
    step7a_final_status <- trimws(
      sub(
        "^Overall status:",
        "",
        status_line[[1]]
      )
    )
  }
}

# ==============================================================================
# 16. REPORT
# ==============================================================================

overall_elapsed <- as.numeric(
  difftime(
    Sys.time(),
    overall_start,
    units = "secs"
  )
)

report_lines <- c(
  "======================================================================",
  "GSE162610 | STEP 7A-2 REPORT",
  "Bioconda pySCENIC environment and mouse-resource setup",
  "======================================================================",
  "",
  paste0("Generated: ", timestamp()),
  paste0("Project root: ", PROJECT_ROOT),
  paste0("Conda environment: ", CONDA_ENV_NAME),
  paste0("Environment Python: ", ENV_PYTHON),
  "Installation method: conda-forge + Bioconda (strict channel priority)",
  "",
  "SCOPE:",
  "This script rebuilt the environment from conda-forge + Bioconda and configured resources only.",
  "It did not run GRNBoost2, motif pruning, AUCell, RSS, or CSI.",
  "",
  "TIMING:",
  paste0(
    "Bioconda environment creation/installation: ",
    format_duration(env_elapsed)
  ),
  paste0(
    "Resource downloads: ",
    format_duration(resource_elapsed)
  ),
  paste0(
    "Step7A re-run: ",
    format_duration(rerun_elapsed)
  ),
  paste0(
    "Total elapsed time: ",
    format_duration(overall_elapsed)
  ),
  "",
  "PYTHON PACKAGE AUDIT:",
  paste(
    capture.output(
      print(
        python_package_audit,
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
        download_audit[
          ,
          c(
            "role",
            "action",
            "download_status",
            "after_size_human",
            "elapsed_human"
          ),
          drop = FALSE
        ],
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
  "TF OVERLAP:",
  paste(
    capture.output(
      print(
        tf_overlap_audit,
        row.names = FALSE
      )
    ),
    collapse = "\n"
  ),
  "",
  "FINAL SETUP DECISION:",
  paste(
    capture.output(
      print(
        setup_decision,
        row.names = FALSE
      )
    ),
    collapse = "\n"
  ),
  "",
  paste0(
    "Step7A automatic re-run requested: ",
    AUTO_RERUN_STEP7A
  ),
  paste0(
    "Step7A automatic re-run status code: ",
    rerun_status
  ),
  paste0(
    "Updated Step7A overall status: ",
    ifelse(
      is.na(step7a_final_status),
      "<not available>",
      step7a_final_status
    )
  ),
  "",
  "NEXT ACTION:",
  if (
    SETUP_PASS &&
    identical(
      step7a_final_status,
      "PASS"
    )
  ) {
    paste0(
      "Environment/resources and the repeated Step7A audit passed. ",
      "The next formal step can be Step 7B GRNBoost2."
    )
  } else if (SETUP_PASS) {
    paste0(
      "Environment/resources passed, but the repeated Step7A report did not ",
      "show PASS. Upload the requested files before starting GRNBoost2."
    )
  } else {
    paste0(
      "Do not start GRNBoost2. Upload the Step7A2 report and failed audit ",
      "tables so the failed setup component can be repaired."
    )
  },
  "",
  "UPLOAD THESE FILES (maximum 5):",
  "1) reports/Step7A2_report.txt",
  "2) tables/02_python_package_audit.csv",
  "3) tables/05_feather_readability_gene_overlap_audit.csv",
  "4) tables/06_mouse_TF_expression_overlap_audit.csv",
  "5) the updated Step7A report: step7A.../reports/Step7A_report.txt",
  "",
  "Do not upload the large feather databases."
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
msg("STEP 7A-2 COMPLETE")
msg("==============================================================")
msg(
  "Environment/resource setup status: ",
  ifelse(
    SETUP_PASS,
    "PASS",
    "FAIL"
  )
)
msg(
  "Updated Step 7A status: ",
  ifelse(
    is.na(step7a_final_status),
    "<not available>",
    step7a_final_status
  )
)
msg(
  "Total elapsed time: ",
  format_duration(overall_elapsed)
)
msg("Report: ", REPORT_FILE)
msg("")
msg("请上传以下文件（最多5个）：")
msg("1) reports/Step7A2_report.txt")
msg("2) tables/02_python_package_audit.csv")
msg("3) tables/05_feather_readability_gene_overlap_audit.csv")
msg("4) tables/06_mouse_TF_expression_overlap_audit.csv")
msg("5) 更新后的 Step7A_report.txt")
