# Repository-level release validation.
# This script performs syntax/consistency checks only. It does NOT rerun analyses.

options(stringsAsFactors = FALSE)

repo_root <- normalizePath(".", winslash = "/", mustWork = TRUE)

required_files <- c(
  "README.md",
  "config.R",
  "public_datasets.csv",
  "analysis_parameters.csv",
  "random_seeds.csv",
  "software_versions.csv",
  "script_manifest.csv",
  "metadata/GSE162610_cell_annotations.csv.gz",
  "metadata/GSE234774_myeloid_cell_annotations.csv.gz"
)
missing_required <- required_files[!file.exists(required_files)]
stopifnot(length(missing_required) == 0L)

# Parse every released R script without executing it.
r_scripts <- list.files("scripts", pattern = "\\.R$", recursive = TRUE, full.names = TRUE)
stopifnot(length(r_scripts) == 114L)
invisible(lapply(r_scripts, parse, keep.source = FALSE, encoding = "UTF-8"))

# Check that known private workstation paths do not occur in released text/code.
text_files <- list.files(
  ".",
  pattern = "\\.(R|py|md|csv|txt)$",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)
forbidden <- c(
  paste0("C:/Users/", "Jiang Wu"),
  paste0("C:\\Users\\", "JIANGW"),
  paste0("E:/", "4scRNA"),
  paste0("Figure5_", "handoff_staging_")
)
bad_hits <- list()
for (f in text_files) {
  x <- tryCatch(readLines(f, warn = FALSE, encoding = "UTF-8"),
                error = function(e) character())
  if (!length(x)) next
  hit <- which(vapply(forbidden, function(p) any(grepl(p, x, fixed = TRUE)),
                       logical(1)))
  if (length(hit)) bad_hits[[f]] <- forbidden[hit]
}
stopifnot(length(bad_hits) == 0L)

# No large intermediate analysis-object formats should be tracked.
all_files <- list.files(".", recursive = TRUE, full.names = TRUE, all.files = TRUE)
blocked <- grepl(
  "\\.(rds|rdata|rda|loom|mtx|h5|h5ad|hdf5|feather|tar|tgz)$",
  all_files,
  ignore.case = TRUE
)
stopifnot(!any(blocked))

# GitHub hard file-size safety margin.
fi <- file.info(all_files[file.exists(all_files) & !dir.exists(all_files)])
stopifnot(all(fi$size < 50 * 1024^2, na.rm = TRUE))

# GSE162610 public annotation audit.
g162 <- read.csv(
  gzfile("metadata/GSE162610_cell_annotations.csv.gz"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
stopifnot(
  nrow(g162) == 62842L,
  length(unique(g162$cell_id)) == 62842L,
  sum(g162$myeloid_included %in% TRUE, na.rm = TRUE) == 32405L,
  sum(g162$clean_resident_microglia %in% TRUE, na.rm = TRUE) == 19958L,
  sum(g162$trajectory_included %in% TRUE, na.rm = TRUE) == 14005L,
  sum(g162$formal_branch_included %in% TRUE, na.rm = TRUE) == 4203L,
  sum(g162$formal_branch_included %in% TRUE & is.na(g162$branch_label)) == 0L
)

# GSE234774 public myeloid annotation audit.
g234 <- read.csv(
  gzfile("metadata/GSE234774_myeloid_cell_annotations.csv.gz"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
stopifnot(
  nrow(g234) == 25205L,
  length(unique(g234$cell_barcode)) == 25205L,
  sum(g234$resident_lineage_status == "RESIDENT_LINEAGE_INCLUDED", na.rm = TRUE) == 20829L,
  sum(g234$resident_lineage_status == "NON_MICROGLIAL_MYELOID_EXCLUDED", na.rm = TRUE) == 3997L,
  sum(g234$resident_lineage_status == "CONTAMINANT_EXCLUDED", na.rm = TRUE) == 379L
)

# Optional Python syntax check when Python is available.
python <- Sys.which("python")
if (!nzchar(python)) python <- Sys.which("python3")
if (nzchar(python)) {
  py_files <- list.files("scripts", pattern = "\\.py$", recursive = TRUE, full.names = TRUE)
  stopifnot(length(py_files) == 10L)
  for (f in py_files) {
    status <- system2(python, c("-m", "py_compile", shQuote(f)),
                      stdout = TRUE, stderr = TRUE)
    code <- attr(status, "status")
    if (is.null(code)) code <- 0L
    stopifnot(code == 0L)
  }
}

cat(
  "PASS: repository structure, R parsing, private-path audit, file-size/object checks, ",
  "GSE162610 annotations, GSE234774 annotations, and optional Python syntax checks.\\n",
  sep = ""
)
