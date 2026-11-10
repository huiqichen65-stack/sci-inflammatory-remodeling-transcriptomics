# ==============================================================================
# 35_step7B_GSE162610_GRNBoost2_network_inference_v6_arboreto_core_fix.R
# ==============================================================================
# Stable targeted repair for:
#   AttributeError: 'csc_matrix' object has no attribute 'A'
#
# The previous v5 script wrote a .pth file to the wrong Conda directory.
# This script first removes every previous compatibility .pth/module, then
# patches the exact obsolete line in arboreto 0.1.6:
#
#   target_gene_expression.A.flatten()
#       ->
#   target_gene_expression.toarray().flatten()
#
# A backup of arboreto/core.py is created before modification.
# No Conda package is reinstalled or downgraded.
# After verification, the existing Step 7B v2 script is rerun automatically.
# ==============================================================================

rm(list = ls())
gc()

options(
  stringsAsFactors = FALSE,
  width = 220,
  warn = 1,
  timeout = 86400
)

PROJECT_ROOT <- Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
CONDA_ENV_NAME <- "pyscenic_env"

timestamp <- function() {
  format(Sys.time(), "%Y-%m-%d %H:%M:%S")
}

msg <- function(...) {
  cat(
    sprintf(
      "[%s] %s\n",
      timestamp(),
      paste0(..., collapse = "")
    )
  )
  flush.console()
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

run_capture <- function(
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

msg("==============================================================")
msg("GSE162610 Step 7B | arboreto core compatibility repair v6")
msg("==============================================================")

python_candidates <- unique(
  c(
    file.path(
      Sys.getenv("USERPROFILE"),
      "anaconda3",
      "envs",
      CONDA_ENV_NAME,
      "python.exe"
    ),
    file.path(
      Sys.getenv("USERPROFILE"),
      "miniconda3",
      "envs",
      CONDA_ENV_NAME,
      "python.exe"
    ),
    "C:/ProgramData/anaconda3/envs/pyscenic_env/python.exe",
    "C:/ProgramData/miniconda3/envs/pyscenic_env/python.exe"
  )
)

python_candidates <- python_candidates[
  file.exists(python_candidates)
]

if (length(python_candidates) == 0L) {
  stop(
    "Cannot find the pyscenic_env Python executable.",
    call. = FALSE
  )
}

ENV_PYTHON <- safe_path(
  python_candidates[[1]]
)

ENV_ROOT <- safe_path(
  dirname(ENV_PYTHON)
)

msg("Python: ", ENV_PYTHON)
msg("Conda environment root: ", ENV_ROOT)

# ------------------------------------------------------------------------------
# 1. Remove every failed compatibility module/.pth created by v4/v5.
# ------------------------------------------------------------------------------

compatibility_files <- list.files(
  path = ENV_ROOT,
  pattern = "^pyscenic_scipy_compat(\\.py|\\.pth|\\.pyc)$",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = FALSE
)

compatibility_cache_dirs <- list.files(
  path = ENV_ROOT,
  pattern = "^__pycache__$",
  recursive = TRUE,
  full.names = TRUE,
  include.dirs = TRUE
)

compatibility_cache_files <- character(0)

if (length(compatibility_cache_dirs) > 0L) {
  compatibility_cache_files <- unlist(
    lapply(
      compatibility_cache_dirs,
      function(cache_dir) {
        list.files(
          path = cache_dir,
          pattern = "^pyscenic_scipy_compat.*\\.pyc$",
          full.names = TRUE
        )
      }
    ),
    use.names = FALSE
  )
}

files_to_remove <- unique(
  c(
    compatibility_files,
    compatibility_cache_files
  )
)

if (length(files_to_remove) > 0L) {
  for (target in files_to_remove) {
    unlink(
      target,
      force = TRUE
    )
  }
}

remaining_compatibility_files <- list.files(
  path = ENV_ROOT,
  pattern = "^pyscenic_scipy_compat(\\.py|\\.pth|\\.pyc)$",
  recursive = TRUE,
  full.names = TRUE
)

if (length(remaining_compatibility_files) > 0L) {
  stop(
    paste0(
      "Failed to remove previous compatibility files:\n",
      paste(
        safe_path(remaining_compatibility_files),
        collapse = "\n"
      )
    ),
    call. = FALSE
  )
}

msg(
  "Removed previous failed compatibility files: ",
  length(files_to_remove)
)

# ------------------------------------------------------------------------------
# 2. Confirm Python starts cleanly and locate the installed arboreto/core.py.
# ------------------------------------------------------------------------------

clean_probe <- run_capture(
  command = ENV_PYTHON,
  args = c(
    "-c",
    quote_cmd(
      paste0(
        "import scipy, arboreto, pathlib; ",
        "print('scipy=' + scipy.__version__); ",
        "print(pathlib.Path(arboreto.__file__).resolve().parent / 'core.py')"
      )
    )
  )
)

if (clean_probe$status != 0L) {
  stop(
    paste(
      c(
        "Python still does not start cleanly after removing the bad .pth file:",
        clean_probe$output
      ),
      collapse = "\n"
    ),
    call. = FALSE
  )
}

probe_lines <- clean_probe$output[
  nzchar(clean_probe$output)
]

CORE_FILE <- safe_path(
  tail(
    probe_lines,
    1L
  )
)

if (!file.exists(CORE_FILE)) {
  stop(
    paste0(
      "arboreto/core.py was not found at: ",
      CORE_FILE
    ),
    call. = FALSE
  )
}

msg(
  "Python clean-start PASS: ",
  paste(
    head(probe_lines, -1L),
    collapse = " | "
  )
)

msg("arboreto core: ", CORE_FILE)

# ------------------------------------------------------------------------------
# 3. Back up and patch the one obsolete SciPy sparse-matrix expression.
# ------------------------------------------------------------------------------

OLD_CODE <- "target_gene_expression = target_gene_expression.A.flatten()"
NEW_CODE <- "target_gene_expression = target_gene_expression.toarray().flatten()"

core_text <- paste(
  readLines(
    CORE_FILE,
    warn = FALSE,
    encoding = "UTF-8"
  ),
  collapse = "\n"
)

old_occurrences <- lengths(
  gregexpr(
    OLD_CODE,
    core_text,
    fixed = TRUE
  )
)

if (
  length(old_occurrences) == 1L &&
  old_occurrences[[1]] == 1L
) {
  BACKUP_FILE <- paste0(
    CORE_FILE,
    ".backup_before_step7B_scipy_fix"
  )

  if (!file.exists(BACKUP_FILE)) {
    copied <- file.copy(
      from = CORE_FILE,
      to = BACKUP_FILE,
      overwrite = FALSE,
      copy.mode = TRUE,
      copy.date = TRUE
    )

    if (!isTRUE(copied)) {
      stop(
        "Could not create the arboreto/core.py backup.",
        call. = FALSE
      )
    }
  }

  patched_text <- sub(
    pattern = OLD_CODE,
    replacement = NEW_CODE,
    x = core_text,
    fixed = TRUE
  )

  writeLines(
    strsplit(
      patched_text,
      split = "\n",
      fixed = TRUE
    )[[1]],
    con = CORE_FILE,
    useBytes = TRUE
  )

  patch_action <- "PATCHED"
} else if (
  grepl(
    NEW_CODE,
    core_text,
    fixed = TRUE
  )
) {
  BACKUP_FILE <- paste0(
    CORE_FILE,
    ".backup_before_step7B_scipy_fix"
  )

  patch_action <- "ALREADY_PATCHED"
} else {
  stop(
    paste0(
      "The expected arboreto 0.1.6 code line was not found, so no automatic ",
      "modification was made. Core file: ",
      CORE_FILE
    ),
    call. = FALSE
  )
}

msg("arboreto patch action: ", patch_action)

# ------------------------------------------------------------------------------
# 4. Verify syntax, installed source, and an actual sparse target conversion.
# ------------------------------------------------------------------------------

VERIFY_PY <- file.path(
  tempdir(),
  "verify_step7B_arboreto_core_fix_v6.py"
)

writeLines(
  c(
    "import inspect",
    "import numpy as np",
    "from scipy.sparse import csc_matrix",
    "import arboreto.core as core",
    "",
    "source = inspect.getsource(core.fit_model)",
    "if '.A.flatten()' in source:",
    "    raise RuntimeError('Obsolete .A call remains in arboreto.core.fit_model')",
    "if '.toarray().flatten()' not in source:",
    "    raise RuntimeError('Expected toarray replacement not found')",
    "",
    "target = csc_matrix([[1.0], [0.0], [2.0]])",
    "dense_target = target.toarray().flatten()",
    "expected = np.array([1.0, 0.0, 2.0])",
    "if not np.array_equal(dense_target, expected):",
    "    raise RuntimeError('Sparse target conversion is not numerically exact')",
    "",
    "compile(open(core.__file__, 'rb').read(), core.__file__, 'exec')",
    "",
    "print('arboreto_core=' + core.__file__)",
    "print('sparse_conversion_exact=True')",
    "print('ARBORETO_CORE_COMPATIBILITY_REPAIR_PASS')"
  ),
  con = VERIFY_PY,
  useBytes = TRUE
)

verify_result <- run_capture(
  command = ENV_PYTHON,
  args = c(
    quote_cmd(VERIFY_PY)
  )
)

REPORT_DIR <- file.path(
  PROJECT_ROOT,
  "04_results",
  "step7B_GSE162610_GRNBoost2_network_inference",
  "reports"
)

dir.create(
  REPORT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

REPAIR_REPORT <- file.path(
  REPORT_DIR,
  "Step7B_arboreto_core_compatibility_repair_v6.txt"
)

repair_pass <- verify_result$status == 0L &&
  any(
    grepl(
      "ARBORETO_CORE_COMPATIBILITY_REPAIR_PASS",
      verify_result$output,
      fixed = TRUE
    )
  )

writeLines(
  c(
    "======================================================================",
    "GSE162610 | STEP 7B ARBORETO CORE COMPATIBILITY REPAIR V6",
    "======================================================================",
    paste0("Generated: ", timestamp()),
    paste0("Python: ", ENV_PYTHON),
    paste0("Environment root: ", ENV_ROOT),
    paste0("arboreto core: ", CORE_FILE),
    paste0("Backup: ", BACKUP_FILE),
    paste0("Patch action: ", patch_action),
    paste0(
      "Previous compatibility files removed: ",
      length(files_to_remove)
    ),
    "",
    "EXACT REPLACEMENT:",
    paste0("OLD: ", OLD_CODE),
    paste0("NEW: ", NEW_CODE),
    "",
    "VERIFICATION:",
    verify_result$output,
    "",
    paste0(
      "Repair status: ",
      ifelse(
        repair_pass,
        "PASS",
        "FAIL"
      )
    ),
    "",
    paste0(
      "The repair changes only the removed SciPy sparse-matrix alias .A to ",
      "the equivalent supported method .toarray(). No GRNBoost2 parameters ",
      "or biological input are changed."
    )
  ),
  con = REPAIR_REPORT,
  useBytes = TRUE
)

if (!repair_pass) {
  stop(
    paste(
      c(
        "arboreto core compatibility repair verification failed:",
        verify_result$output,
        paste0(
          "Report: ",
          safe_path(REPAIR_REPORT)
        )
      ),
      collapse = "\n"
    ),
    call. = FALSE
  )
}

msg("arboreto compatibility verification PASS.")
msg("Repair report: ", safe_path(REPAIR_REPORT))

# ------------------------------------------------------------------------------
# 5. Re-run the already validated Step 7B Windows script.
# ------------------------------------------------------------------------------

STEP7B_CANDIDATES <- unique(
  c(
    file.path(
      Sys.getenv("USERPROFILE"),
      "Downloads",
      "35_step7B_GSE162610_GRNBoost2_network_inference_v2_windows_fix.R"
    ),
    file.path(
      PROJECT_ROOT,
      "35_step7B_GSE162610_GRNBoost2_network_inference_v2_windows_fix.R"
    )
  )
)

STEP7B_CANDIDATES <- STEP7B_CANDIDATES[
  file.exists(STEP7B_CANDIDATES)
]

if (length(STEP7B_CANDIDATES) == 0L) {
  stop(
    paste0(
      "The compatibility repair passed, but the Step 7B v2 script was not found.\n",
      "Keep the v2 script in Downloads or the project root, then run it manually."
    ),
    call. = FALSE
  )
}

STEP7B_SCRIPT <- safe_path(
  STEP7B_CANDIDATES[[1]]
)

msg("The existing sparse loom will be reused.")
msg(
  "Automatically re-running Step 7B: ",
  STEP7B_SCRIPT
)
msg("==============================================================")

source(
  STEP7B_SCRIPT,
  chdir = FALSE
)
