# Central path configuration for reproducible execution.
# Run scripts from the repository root, or set the environment variables below.

PROJECT_ROOT <- Sys.getenv(
  "GLIA_PROJECT_ROOT",
  unset = normalizePath(".", winslash = "/", mustWork = FALSE)
)
PROJECT_ROOT <- normalizePath(PROJECT_ROOT, winslash = "/", mustWork = FALSE)

RAW_DATA_DIR <- Sys.getenv(
  "GLIA_RAW_DATA_DIR",
  unset = file.path(PROJECT_ROOT, "01_raw_data")
)
RAW_DATA_DIR <- normalizePath(RAW_DATA_DIR, winslash = "/", mustWork = FALSE)

RESULTS_DIR <- Sys.getenv(
  "GLIA_RESULTS_DIR",
  unset = file.path(PROJECT_ROOT, "04_results")
)
RESULTS_DIR <- normalizePath(RESULTS_DIR, winslash = "/", mustWork = FALSE)

FIGURE_DIR <- Sys.getenv(
  "GLIA_FIGURE_DIR",
  unset = file.path(PROJECT_ROOT, "09_final_figures")
)
FIGURE_DIR <- normalizePath(FIGURE_DIR, winslash = "/", mustWork = FALSE)

REFERENCE_DIR <- Sys.getenv(
  "GLIA_REFERENCE_DIR",
  unset = file.path(RAW_DATA_DIR, "pySCENIC_resources")
)
REFERENCE_DIR <- normalizePath(REFERENCE_DIR, winslash = "/", mustWork = FALSE)

PYSCENIC_PYTHON <- Sys.getenv(
  "PYSCENIC_PYTHON",
  unset = Sys.which("python")
)

project_path <- function(path = "") {
  if (!nzchar(path)) return(PROJECT_ROOT)
  file.path(PROJECT_ROOT, path)
}

raw_data_path <- function(path = "") {
  if (!nzchar(path)) return(RAW_DATA_DIR)
  file.path(RAW_DATA_DIR, path)
}

results_path <- function(path = "") {
  if (!nzchar(path)) return(RESULTS_DIR)
  file.path(RESULTS_DIR, path)
}
