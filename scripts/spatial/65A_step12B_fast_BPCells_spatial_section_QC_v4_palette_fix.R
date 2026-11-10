# ==============================================================================
# 65A_step12B_fast_BPCells_spatial_section_QC_v4_palette_fix.R
#
#
# v2 repair
#   The first version called base::colSums() explicitly. That bypassed the
#   BPCells/Matrix S4 method and caused:
#   "'x' must be an array of at least two dimensions".
#   This version uses Matrix::colSums() for all BPCells column-sum operations.
#   A valid BPCells cache from the interrupted first run is reused automatically.
#
# v3 repair
#   Recent data.table versions can reject DT[!logical_column, ...] when the
#   logical column is supplied as a single symbol in i. This version uses the
#   explicit row condition qc_review_flag == FALSE.
#   The valid BPCells cache is reused automatically; do not delete it.
#
# v4 repair
#   grDevices::hcl.colors() requires a published palette name. The prior script
#   used lowercase "viridis" and the unsupported name "magma", producing:
#   "'palette' does not match any given palette".
#   This version resolves common aliases case-insensitively and falls back to
#   a valid installed HCL palette. The BPCells cache is reused automatically.
#
# Project
#   GSE162610/GSE234774 adult SCI microglia project
#
# Step
#   Step 12B-FAST: GSE234774 spatial 2D BPCells cache, section-aware QC,
#   and image-coordinate review
#
# Why this fast version exists
#   The original Step 12B script materialized the entire Matrix Market matrix
#   in RAM and saved 12 full raw Seurat objects. This version imports the matrix
#   once into a disk-backed BPCells directory, computes only the technical QC
#   required at this stage, and saves a lightweight handoff object rather than
#   12 large count-containing RDS files.
#
# Scope
#   1. Verify the repaired Step 12A2 PASS/PASS_WITH_REVIEW result.
#   2. Import/reuse the full 25,966 x 33,941 spatial 2D count matrix in BPCells.
#   3. Compute nCount, nFeature, and percent.mt without dense conversion.
#   4. Preserve the audited 12 capture areas and 36 condition-specific sections.
#   5. Verify image, position, scalefactor, and processed-barcode linkage.
#   6. Produce one compact six-panel QC page per capture area.
#   7. Save spot-level and section-level QC tables plus a lightweight handoff RDS.
#
# This script does NOT
#   - create or save 12 full raw Seurat objects;
#   - delete any spot;
#   - normalize, select variable genes, scale, run PCA/UMAP, or cluster;
#   - run RCTD or another deconvolution method;
#   - calculate the target inflammatory-remodeling program score;
#   - perform lesion/core/border spatial statistics;
#   - infer age, injury severity, or treatment from filenames.
#
# Required packages
#   BPCells, data.table, jsonlite, png
#
# Run in RStudio
#   source(
#     "<PROJECT_ROOT>/65A_step12B_fast_BPCells_spatial_section_QC_v4_palette_fix.R",
#     encoding = "UTF-8"
#   )
# ==============================================================================

options(
  stringsAsFactors = FALSE,
  warn = 1,
  width = 180
)

# ------------------------------------------------------------------------------
# 0. User controls
# ------------------------------------------------------------------------------

PROJECT_ROOT <- Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
# Keep FALSE in normal use. Set TRUE only when the BPCells matrix directory is
# incomplete/corrupted and must be rebuilt from the original mtx.gz file.
FORCE_REIMPORT_BPCELLS <- FALSE

# A rerun may safely overwrite Step 12B-fast tables/figures while preserving a
# valid BPCells cache. Keep TRUE for convenient checkpoint/resume behavior.
OVERWRITE_SMALL_OUTPUTS <- TRUE

# BPCells import settings. These match the fixed-memory strategy already used in
# the GSE234774 snRNA workflow. The sort buffer is 4 GiB; the import itself is
# disk-backed and does not materialize the complete matrix in RAM.
BPCELLS_SORT_BYTES <- 4294967296
BPCELLS_LOAD_BYTES <- 67108864

# Conservative section-aware review thresholds. These create review labels only;
# no spot is removed in Step 12B-fast.
ABSOLUTE_MIN_NCOUNT <- 100
ABSOLUTE_MIN_NFEATURE <- 100
ABSOLUTE_MIN_SECTION_SPOTS_STOP <- 50
ABSOLUTE_HIGH_MT_FLOOR <- 25
ABSOLUTE_HIGH_MT_CEILING <- 50
ROBUST_MAD_MULTIPLIER <- 5
SECTION_REVIEW_FLAG_FRACTION <- 0.20

EXPECTED_CAPTURE_AREAS <- c(
  "078_A", "078_B", "078_C", "078_D",
  "105_A", "105_B", "105_C", "105_D",
  "106_A", "106_B", "106_C", "106_D"
)

EXPECTED_TIMEPOINTS <- c(
  "uninjured",
  "7days",
  "2months"
)

# ------------------------------------------------------------------------------
# 1. Paths
# ------------------------------------------------------------------------------

INPUT_DIR <- file.path(
  PROJECT_ROOT,
  "01_raw_data",
  "GSE234774"
)

STEP12A2_DIR <- file.path(
  PROJECT_ROOT,
  "04_results",
  "step12A2_GSE234774_spatial_metadata_asset_linkage_repair"
)

STEP12A2_REPAIR_DIR <- file.path(
  PROJECT_ROOT,
  "04_results",
  "step12A2_false_STOP_repair_only"
)

OUTPUT_DIR <- file.path(
  PROJECT_ROOT,
  "04_results",
  "step12B_fast_BPCells_spatial_section_QC"
)

TABLE_DIR <- file.path(OUTPUT_DIR, "tables")
FIGURE_DIR <- file.path(OUTPUT_DIR, "figures")
BPCELLS_DIR <- file.path(OUTPUT_DIR, "BPCells")
TMP_DIR <- file.path(OUTPUT_DIR, "BPCells_tmp")
ASSET_PREVIEW_DIR <- file.path(OUTPUT_DIR, "lowres_asset_preview")
RDS_DIR <- file.path(OUTPUT_DIR, "RDS")

MATRIX_FILE <- file.path(
  INPUT_DIR,
  "GSE234774_spatial_2d_filtered_spatial_2d.mtx.gz"
)

BARCODE_FILE <- file.path(
  INPUT_DIR,
  "GSE234774_spatial_2d_barcodes.txt.gz"
)

FEATURE_FILE <- file.path(
  INPUT_DIR,
  "GSE234774_spatial_2d_features.txt.gz"
)

METADATA_FILE <- file.path(
  INPUT_DIR,
  "GSE234774_spatial_2d_meta.txt.gz"
)

REPAIRED_ASSET_MANIFEST_FILE <- file.path(
  STEP12A2_REPAIR_DIR,
  "01_repaired_2d_capture_area_asset_manifest.csv"
)

REPAIRED_STEP12A2_REPORT <- file.path(
  STEP12A2_REPAIR_DIR,
  "10_STEP12A2_REPAIRED_FINAL_REPORT.txt"
)

PRIOR_SECTION_MANIFEST_FILE <- file.path(
  STEP12A2_DIR,
  "03_condition_specific_section_manifest.csv"
)

AUDIT_EXTRACTED_ASSET_DIR <- file.path(
  STEP12A2_DIR,
  "audit_extracted_2d_assets"
)

FULL_BPCELLS_MATRIX_DIR <- file.path(
  BPCELLS_DIR,
  "GSE234774_spatial_2d_full_counts"
)

RUN_LOG_FILE <- file.path(OUTPUT_DIR, "00_run_log.txt")
SESSION_INFO_FILE <- file.path(OUTPUT_DIR, "00_sessionInfo.txt")

INPUT_AUDIT_OUTPUT <- file.path(
  TABLE_DIR,
  "01_input_package_BPCells_audit.csv"
)

CAPTURE_ASSET_AUDIT_OUTPUT <- file.path(
  TABLE_DIR,
  "02_capture_area_image_coordinate_audit.csv"
)

SECTION_QC_OUTPUT <- file.path(
  TABLE_DIR,
  "03_section_QC_summary.csv"
)

SPOT_QC_OUTPUT <- file.path(
  TABLE_DIR,
  "04_spot_QC_flags.csv.gz"
)

BPCELLS_HANDOFF_OUTPUT <- file.path(
  TABLE_DIR,
  "05_BPCells_spatial_handoff_index.csv"
)

AUDIT_CHECKS_OUTPUT <- file.path(
  TABLE_DIR,
  "09_STEP12B_FAST_audit_checks.csv"
)

QC_PDF_OUTPUT <- file.path(
  FIGURE_DIR,
  "01_spatial_section_and_QC_maps_fast.pdf"
)

HANDOFF_RDS_OUTPUT <- file.path(
  RDS_DIR,
  "GSE234774_step12B_FAST_lightweight_handoff.rds"
)

FINAL_REPORT_OUTPUT <- file.path(
  OUTPUT_DIR,
  "10_STEP12B_FAST_FINAL_REPORT.txt"
)

for (directory in c(
  OUTPUT_DIR,
  TABLE_DIR,
  FIGURE_DIR,
  BPCELLS_DIR,
  TMP_DIR,
  ASSET_PREVIEW_DIR,
  RDS_DIR
)) {
  dir.create(
    directory,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

if (!OVERWRITE_SMALL_OUTPUTS) {
  small_outputs <- c(
    INPUT_AUDIT_OUTPUT,
    CAPTURE_ASSET_AUDIT_OUTPUT,
    SECTION_QC_OUTPUT,
    SPOT_QC_OUTPUT,
    BPCELLS_HANDOFF_OUTPUT,
    AUDIT_CHECKS_OUTPUT,
    QC_PDF_OUTPUT,
    HANDOFF_RDS_OUTPUT,
    FINAL_REPORT_OUTPUT
  )

  if (any(file.exists(small_outputs))) {
    stop(
      "Step 12B-fast output already exists. Set ",
      "OVERWRITE_SMALL_OUTPUTS <- TRUE to replace the small tables/figures."
    )
  }
}

# ------------------------------------------------------------------------------
# 2. Package preflight
# ------------------------------------------------------------------------------

cran_required <- c(
  "data.table",
  "jsonlite",
  "png",
  "Matrix"
)

missing_cran <- cran_required[
  !vapply(
    cran_required,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_cran) > 0L) {
  stop(
    "Missing required package(s): ",
    paste(missing_cran, collapse = ", "),
    ". Install once with:\n",
    "install.packages(c(",
    paste(
      sprintf('"%s"', missing_cran),
      collapse = ", "
    ),
    "))"
  )
}

if (!requireNamespace("BPCells", quietly = TRUE)) {
  stop(
    "BPCells is not installed. Install it once, restart RStudio, and rerun:\n",
    'install.packages("BPCells", repos=c(',
    '"https://bnprks.r-universe.dev",',
    '"https://cloud.r-project.org"))'
  )
}

if (!exists(
  "import_matrix_market",
  where = asNamespace("BPCells"),
  inherits = FALSE
)) {
  stop(
    "The installed BPCells version lacks import_matrix_market(). ",
    "Update BPCells before rerunning Step 12B-fast."
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(BPCells)
})

# ------------------------------------------------------------------------------
# 3. Logging and general helpers
# ------------------------------------------------------------------------------

run_log <- character()

log_message <- function(...) {
  message_text <- paste0(...)
  timestamp <- format(
    Sys.time(),
    "%Y-%m-%d %H:%M:%S"
  )
  line <- paste0(
    "[",
    timestamp,
    "] ",
    message_text
  )
  cat(line, "\n")
  run_log <<- c(run_log, line)
}

write_run_log <- function() {
  writeLines(
    run_log,
    RUN_LOG_FILE,
    useBytes = TRUE
  )
}

elapsed_minutes <- function(start_time) {
  round(
    as.numeric(
      difftime(
        Sys.time(),
        start_time,
        units = "mins"
      )
    ),
    3
  )
}

format_bytes <- function(number_bytes) {
  if (!is.finite(number_bytes)) {
    return(NA_character_)
  }

  units <- c(
    "B",
    "KiB",
    "MiB",
    "GiB",
    "TiB"
  )

  value <- number_bytes
  unit_index <- 1L

  while (
    value >= 1024 &&
      unit_index < length(units)
  ) {
    value <- value / 1024
    unit_index <- unit_index + 1L
  }

  paste0(
    round(value, 2),
    " ",
    units[[unit_index]]
  )
}

trim_na <- function(x) {
  x <- trimws(as.character(x))
  x[x == ""] <- NA_character_
  x
}

collapse_unique <- function(
    x,
    separator = " | ",
    max_values = 30L) {

  x <- unique(trim_na(x))
  x <- x[!is.na(x)]

  if (length(x) == 0L) {
    return(NA_character_)
  }

  if (length(x) > max_values) {
    x <- c(
      x[seq_len(max_values)],
      paste0(
        "... +",
        length(x) - max_values,
        " more"
      )
    )
  }

  paste(
    x,
    collapse = separator
  )
}

as_logical_strict <- function(
    x,
    field_name) {

  if (is.logical(x)) {
    return(x)
  }

  value <- toupper(
    trimws(
      as.character(x)
    )
  )

  output <- rep(
    NA,
    length(value)
  )

  output[value %in% c(
    "TRUE",
    "T",
    "1"
  )] <- TRUE

  output[value %in% c(
    "FALSE",
    "F",
    "0"
  )] <- FALSE

  bad <- is.na(output) &
    !is.na(value) &
    nzchar(value)

  if (any(bad)) {
    stop(
      "Column '",
      field_name,
      "' has non-logical value(s): ",
      paste(
        unique(value[bad]),
        collapse = " | "
      )
    )
  }

  output
}

required_columns_present <- function(
    table,
    required_columns,
    table_name) {

  missing_columns <- setdiff(
    required_columns,
    names(table)
  )

  if (length(missing_columns) > 0L) {
    stop(
      table_name,
      " is missing required column(s): ",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  }

  invisible(TRUE)
}

regex_escape <- function(x) {
  gsub(
    "([][{}()+*^$|\\\\?.])",
    "\\\\\\1",
    x
  )
}

safe_median <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]

  if (length(x) == 0L) {
    return(NA_real_)
  }

  stats::median(
    x,
    na.rm = TRUE
  )
}

safe_quantile <- function(
    x,
    probability) {

  x <- as.numeric(x)
  x <- x[is.finite(x)]

  if (length(x) == 0L) {
    return(NA_real_)
  }

  as.numeric(
    stats::quantile(
      x,
      probs = probability,
      names = FALSE,
      na.rm = TRUE,
      type = 8
    )
  )
}

safe_mad <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]

  if (length(x) == 0L) {
    return(NA_real_)
  }

  stats::mad(
    x,
    center = stats::median(x),
    constant = 1.4826,
    na.rm = TRUE
  )
}

robust_lower_threshold <- function(
    values,
    absolute_floor,
    mad_multiplier = ROBUST_MAD_MULTIPLIER) {

  center <- safe_median(values)
  dispersion <- safe_mad(values)

  if (!is.finite(center)) {
    return(as.numeric(absolute_floor))
  }

  if (!is.finite(dispersion) || dispersion == 0) {
    dispersion <- stats::IQR(
      values,
      na.rm = TRUE
    ) / 1.349
  }

  if (!is.finite(dispersion) || dispersion == 0) {
    dispersion <- 0
  }

  max(
    as.numeric(absolute_floor),
    floor(
      center -
        mad_multiplier * dispersion
    )
  )
}

robust_upper_mt_threshold <- function(
    values,
    floor_value = ABSOLUTE_HIGH_MT_FLOOR,
    ceiling_value = ABSOLUTE_HIGH_MT_CEILING,
    mad_multiplier = ROBUST_MAD_MULTIPLIER) {

  center <- safe_median(values)
  dispersion <- safe_mad(values)

  if (!is.finite(center)) {
    return(as.numeric(ceiling_value))
  }

  if (!is.finite(dispersion) || dispersion == 0) {
    dispersion <- stats::IQR(
      values,
      na.rm = TRUE
    ) / 1.349
  }

  if (!is.finite(dispersion) || dispersion == 0) {
    dispersion <- 0
  }

  proposed <- center +
    mad_multiplier * dispersion

  min(
    as.numeric(ceiling_value),
    max(
      as.numeric(floor_value),
      proposed
    )
  )
}

append_reason <- function(
    current,
    condition,
    label) {

  output <- current
  indices <- which(condition)

  if (length(indices) > 0L) {
    output[indices] <- ifelse(
      nzchar(output[indices]),
      paste0(
        output[indices],
        ";",
        label
      ),
      label
    )
  }

  output
}

open_bpcells_safe <- function(path) {
  tryCatch(
    BPCells::open_matrix_dir(path),
    error = function(e) e
  )
}

directory_size_bytes <- function(path) {
  if (!dir.exists(path)) {
    return(NA_real_)
  }

  files <- list.files(
    path,
    recursive = TRUE,
    full.names = TRUE,
    all.files = TRUE,
    no.. = TRUE
  )

  if (length(files) == 0L) {
    return(0)
  }

  info <- file.info(files)
  sum(
    info$size[!info$isdir],
    na.rm = TRUE
  )
}

decompress_gzip_binary <- function(
    source_file,
    destination_file) {

  if (!file.exists(source_file)) {
    stop(
      "Compressed file not found: ",
      source_file
    )
  }

  dir.create(
    dirname(destination_file),
    recursive = TRUE,
    showWarnings = FALSE
  )

  input_connection <- gzfile(
    source_file,
    open = "rb"
  )

  output_connection <- file(
    destination_file,
    open = "wb"
  )

  on.exit(
    {
      try(
        close(input_connection),
        silent = TRUE
      )
      try(
        close(output_connection),
        silent = TRUE
      )
    },
    add = TRUE
  )

  repeat {
    buffer <- readBin(
      input_connection,
      what = "raw",
      n = 1024L * 1024L
    )

    if (length(buffer) == 0L) {
      break
    }

    writeBin(
      buffer,
      output_connection
    )
  }

  close(input_connection)
  close(output_connection)

  if (
    !file.exists(destination_file) ||
      file.info(destination_file)$size <= 0
  ) {
    stop(
      "Decompression produced an empty file: ",
      destination_file
    )
  }

  invisible(destination_file)
}

read_gzip_text <- function(path) {
  connection <- gzfile(
    path,
    open = "rt"
  )

  on.exit(
    close(connection),
    add = TRUE
  )

  paste(
    readLines(
      connection,
      warn = FALSE,
      encoding = "UTF-8"
    ),
    collapse = "\n"
  )
}

read_position_file <- function(path) {
  position_table <- data.table::fread(
    path,
    header = FALSE,
    data.table = TRUE,
    showProgress = FALSE
  )

  if (ncol(position_table) < 6L) {
    stop(
      "Position file has fewer than six columns: ",
      basename(path)
    )
  }

  first_value <- tolower(
    trimws(
      as.character(
        position_table[[1L]][1L]
      )
    )
  )

  if (first_value %in% c(
    "barcode",
    "barcodes"
  )) {
    position_table <- data.table::fread(
      path,
      header = TRUE,
      data.table = TRUE,
      showProgress = FALSE
    )

    normalized_names <- tolower(
      gsub(
        "[^a-z0-9]+",
        "_",
        names(position_table)
      )
    )

    setnames(
      position_table,
      normalized_names
    )

    barcode_candidates <- c(
      "barcode",
      "barcodes"
    )

    barcode_column <- barcode_candidates[
      barcode_candidates %in%
        names(position_table)
    ][1L]

    if (
      length(barcode_column) == 0L ||
        is.na(barcode_column)
    ) {
      stop(
        "No barcode column found in position file: ",
        basename(path)
      )
    }

    requested_columns <- c(
      barcode_column,
      "in_tissue",
      "array_row",
      "array_col",
      "pxl_row_in_fullres",
      "pxl_col_in_fullres"
    )

    missing_columns <- setdiff(
      requested_columns,
      names(position_table)
    )

    if (length(missing_columns) > 0L) {
      stop(
        "Position file is missing: ",
        paste(
          missing_columns,
          collapse = ", "
        )
      )
    }

    position_table <- position_table[
      ,
      ..requested_columns
    ]

    setnames(
      position_table,
      c(
        "barcode",
        "in_tissue",
        "array_row",
        "array_col",
        "pxl_row_in_fullres",
        "pxl_col_in_fullres"
      )
    )
  } else {
    position_table <- position_table[
      ,
      1:6
    ]

    setnames(
      position_table,
      c(
        "barcode",
        "in_tissue",
        "array_row",
        "array_col",
        "pxl_row_in_fullres",
        "pxl_col_in_fullres"
      )
    )
  }

  position_table[
    ,
    barcode := trimws(
      as.character(barcode)
    )
  ]

  position_table[
    ,
    in_tissue := suppressWarnings(
      as.integer(in_tissue)
    )
  ]

  position_table[
    ,
    array_row := suppressWarnings(
      as.integer(array_row)
    )
  ]

  position_table[
    ,
    array_col := suppressWarnings(
      as.integer(array_col)
    )
  ]

  position_table[
    ,
    pxl_row_in_fullres := suppressWarnings(
      as.numeric(pxl_row_in_fullres)
    )
  ]

  position_table[
    ,
    pxl_col_in_fullres := suppressWarnings(
      as.numeric(pxl_col_in_fullres)
    )
  ]

  position_table
}

read_features <- function(path) {
  feature_table <- data.table::fread(
    path,
    header = FALSE,
    sep = "\t",
    data.table = FALSE,
    fill = TRUE,
    showProgress = FALSE
  )

  if (ncol(feature_table) == 1L) {
    names(feature_table) <- "feature_id"
    feature_table$feature_name <-
      feature_table$feature_id
    feature_table$feature_type <-
      "Gene Expression"
  } else if (ncol(feature_table) == 2L) {
    names(feature_table) <- c(
      "feature_id",
      "feature_name"
    )
    feature_table$feature_type <-
      "Gene Expression"
  } else {
    names(feature_table)[1:3] <- c(
      "feature_id",
      "feature_name",
      "feature_type"
    )
  }

  feature_table$feature_id <- trim_na(
    feature_table$feature_id
  )

  feature_table$feature_name <- trim_na(
    feature_table$feature_name
  )

  missing_name <- is.na(
    feature_table$feature_name
  )

  feature_table$feature_name[missing_name] <-
    feature_table$feature_id[missing_name]

  feature_table$feature_name_unique <- make.unique(
    feature_table$feature_name,
    sep = "_dup"
  )

  feature_table
}

# ------------------------------------------------------------------------------
# 4. Audit-check helpers
# ------------------------------------------------------------------------------

audit_checks <- data.table(
  check_id = character(),
  category = character(),
  description = character(),
  status = character(),
  observed = character(),
  required = character(),
  critical = logical(),
  action_if_not_pass = character()
)

add_check <- function(
    check_id,
    category,
    description,
    passed,
    observed,
    required,
    critical = TRUE,
    action_if_not_pass = "") {

  audit_checks <<- rbind(
    audit_checks,
    data.table(
      check_id = as.character(check_id),
      category = as.character(category),
      description = as.character(description),
      status = if (
        isTRUE(passed)
      ) {
        "PASS"
      } else {
        "FAIL"
      },
      observed = as.character(observed),
      required = as.character(required),
      critical = isTRUE(critical),
      action_if_not_pass = as.character(
        action_if_not_pass
      )
    ),
    fill = TRUE
  )
}

add_review <- function(
    check_id,
    category,
    description,
    observed,
    required,
    action_if_not_pass) {

  audit_checks <<- rbind(
    audit_checks,
    data.table(
      check_id = as.character(check_id),
      category = as.character(category),
      description = as.character(description),
      status = "REVIEW",
      observed = as.character(observed),
      required = as.character(required),
      critical = FALSE,
      action_if_not_pass = as.character(
        action_if_not_pass
      )
    ),
    fill = TRUE
  )
}

write_fatal_report <- function(error_message) {
  report_lines <- c(
    "============================================================",
    "GSE234774 Step 12B-FAST BPCells spatial section QC",
    "============================================================",
    paste0(
      "Generated: ",
      format(
        Sys.time(),
        "%Y-%m-%d %H:%M:%S"
      )
    ),
    paste0(
      "Output directory: ",
      OUTPUT_DIR
    ),
    "",
    "FINAL_STATUS: STOP",
    "",
    "FATAL ERROR",
    paste0(
      "- ",
      error_message
    ),
    "",
    "BOUNDARIES",
    "- No RCTD was run.",
    "- No target-program score was calculated.",
    "- No normalization, PCA, clustering, or spatial statistics were run.",
    "",
    "NEXT ACTION",
    "- Upload this report and 00_run_log.txt.",
    "- Do not delete a valid BPCells cache.",
    "- Do not proceed to RCTD or target scoring.",
    "============================================================"
  )

  writeLines(
    report_lines,
    FINAL_REPORT_OUTPUT,
    useBytes = TRUE
  )

  log_message(
    "FATAL ERROR: ",
    error_message
  )

  write_run_log()
}

# ------------------------------------------------------------------------------
# 5. Plotting helpers
# ------------------------------------------------------------------------------

resolve_hcl_palette <- function(
    palette_name,
    fallback = "Viridis") {

  available_palettes <- grDevices::hcl.pals()

  requested <- as.character(
    palette_name
  )[1L]

  aliases <- c(
    "viridis" = "Viridis",
    "magma" = "Inferno",
    "inferno" = "Inferno",
    "plasma" = "Plasma",
    "rocket" = "Rocket",
    "mako" = "Mako"
  )

  alias_match <- aliases[
    tolower(requested)
  ]

  candidate <- if (
    length(alias_match) == 1L &&
      !is.na(alias_match)
  ) {
    unname(alias_match)
  } else {
    requested
  }

  exact_match <- available_palettes[
    tolower(available_palettes) ==
      tolower(candidate)
  ]

  if (length(exact_match) >= 1L) {
    return(exact_match[[1L]])
  }

  fallback_match <- available_palettes[
    tolower(available_palettes) ==
      tolower(fallback)
  ]

  if (length(fallback_match) >= 1L) {
    return(fallback_match[[1L]])
  }

  available_palettes[[1L]]
}

continuous_point_colors <- function(
    values,
    palette_name = "Viridis") {

  values <- as.numeric(values)
  finite_values <- values[is.finite(values)]

  if (length(finite_values) == 0L) {
    return(
      rep(
        "#808080AA",
        length(values)
      )
    )
  }

  lower_limit <- safe_quantile(
    finite_values,
    0.01
  )

  upper_limit <- safe_quantile(
    finite_values,
    0.99
  )

  if (
    !is.finite(lower_limit) ||
      !is.finite(upper_limit) ||
      upper_limit <= lower_limit
  ) {
    lower_limit <- min(
      finite_values
    )
    upper_limit <- max(
      finite_values
    )
  }

  if (
    !is.finite(lower_limit) ||
      !is.finite(upper_limit) ||
      upper_limit <= lower_limit
  ) {
    return(
      rep(
        "#3B528BAA",
        length(values)
      )
    )
  }

  clipped <- pmin(
    upper_limit,
    pmax(
      lower_limit,
      values
    )
  )

  scaled <- (
    clipped -
      lower_limit
  ) / (
    upper_limit -
      lower_limit
  )

  resolved_palette <- resolve_hcl_palette(
    palette_name
  )

  palette_values <- grDevices::hcl.colors(
    100L,
    palette = resolved_palette
  )

  indices <- pmax(
    1L,
    pmin(
      100L,
      floor(
        scaled * 99
      ) + 1L
    )
  )

  output <- rep(
    "#808080AA",
    length(values)
  )

  valid <- is.finite(indices)

  output[valid] <- grDevices::adjustcolor(
    palette_values[indices[valid]],
    alpha.f = 0.80
  )

  output
}

draw_background_image <- function(
    image_array,
    image_width,
    image_height) {

  graphics::plot.new()

  graphics::plot.window(
    xlim = c(
      0,
      image_width
    ),
    ylim = c(
      image_height,
      0
    ),
    asp = 1
  )

  graphics::rasterImage(
    image_array,
    xleft = 0,
    ybottom = image_height,
    xright = image_width,
    ytop = 0,
    interpolate = FALSE
  )

  graphics::box(
    col = "#666666"
  )
}

draw_continuous_spatial_panel <- function(
    image_array,
    image_width,
    image_height,
    x,
    y,
    values,
    title_text,
    palette_name = "viridis") {

  draw_background_image(
    image_array,
    image_width,
    image_height
  )

  point_colors <- continuous_point_colors(
    values,
    palette_name = palette_name
  )

  graphics::points(
    x,
    y,
    pch = 16,
    cex = 0.42,
    col = point_colors
  )

  graphics::title(
    main = title_text,
    cex.main = 0.92
  )

  finite_values <- values[
    is.finite(values)
  ]

  if (length(finite_values) > 0L) {
    graphics::mtext(
      paste0(
        "min/median/max: ",
        round(
          min(finite_values),
          2
        ),
        " / ",
        round(
          stats::median(finite_values),
          2
        ),
        " / ",
        round(
          max(finite_values),
          2
        )
      ),
      side = 1,
      line = -1.2,
      cex = 0.65
    )
  }
}

draw_section_panel <- function(
    image_array,
    image_width,
    image_height,
    x,
    y,
    timepoint,
    section_id,
    title_text) {

  draw_background_image(
    image_array,
    image_width,
    image_height
  )

  timepoint_levels <- c(
    "uninjured",
    "7days",
    "2months"
  )

  timepoint_colors <- c(
    "uninjured" = "#2C7BB6CC",
    "7days" = "#FDAE61CC",
    "2months" = "#7B3294CC"
  )

  color_vector <- timepoint_colors[
    as.character(timepoint)
  ]

  color_vector[
    is.na(color_vector)
  ] <- "#808080CC"

  graphics::points(
    x,
    y,
    pch = 16,
    cex = 0.47,
    col = color_vector
  )

  graphics::title(
    main = title_text,
    cex.main = 0.92
  )

  present_levels <- timepoint_levels[
    timepoint_levels %in%
      unique(
        as.character(timepoint)
      )
  ]

  graphics::legend(
    "topright",
    legend = present_levels,
    col = timepoint_colors[
      present_levels
    ],
    pch = 16,
    pt.cex = 1.0,
    cex = 0.68,
    bty = "n"
  )

  section_centers <- data.table(
    x = x,
    y = y,
    section_id = section_id
  )[
    ,
    .(
      center_x = stats::median(
        x,
        na.rm = TRUE
      ),
      center_y = stats::median(
        y,
        na.rm = TRUE
      )
    ),
    by = section_id
  ]

  graphics::text(
    section_centers$center_x,
    section_centers$center_y,
    labels = section_centers$section_id,
    cex = 0.54,
    font = 2,
    col = "#111111"
  )
}

draw_qc_flag_panel <- function(
    image_array,
    image_width,
    image_height,
    x,
    y,
    qc_review_flag,
    title_text) {

  draw_background_image(
    image_array,
    image_width,
    image_height
  )

  colors <- ifelse(
    qc_review_flag,
    "#D73027DD",
    "#1A9850AA"
  )

  graphics::points(
    x,
    y,
    pch = 16,
    cex = 0.47,
    col = colors
  )

  graphics::title(
    main = title_text,
    cex.main = 0.92
  )

  graphics::legend(
    "topright",
    legend = c(
      "provisional pass",
      "QC review"
    ),
    col = c(
      "#1A9850AA",
      "#D73027DD"
    ),
    pch = 16,
    pt.cex = 1.0,
    cex = 0.68,
    bty = "n"
  )
}

draw_summary_panel <- function(
    capture_area_id,
    capture_spot_count,
    section_ids,
    section_spot_counts,
    qc_review_count,
    qc_review_fraction,
    zero_count_spots,
    position_overlap_fraction,
    in_tissue_overlap_fraction,
    bpcells_mode) {

  graphics::plot.new()

  summary_lines <- c(
    paste0(
      "Capture area: ",
      capture_area_id
    ),
    paste0(
      "Processed spots: ",
      capture_spot_count
    ),
    paste0(
      "Sections: ",
      paste(
        paste0(
          section_ids,
          " (n=",
          section_spot_counts,
          ")"
        ),
        collapse = "; "
      )
    ),
    paste0(
      "QC review spots: ",
      qc_review_count,
      " (",
      round(
        100 * qc_review_fraction,
        2
      ),
      "%)"
    ),
    paste0(
      "Zero-count spots: ",
      zero_count_spots
    ),
    paste0(
      "Processed-position overlap: ",
      round(
        100 * position_overlap_fraction,
        2
      ),
      "%"
    ),
    paste0(
      "Processed-in-tissue overlap: ",
      round(
        100 * in_tissue_overlap_fraction,
        2
      ),
      "%"
    ),
    paste0(
      "BPCells mode: ",
      bpcells_mode
    ),
    "",
    "No spot was removed.",
    "No RCTD or target score was run."
  )

  graphics::text(
    x = 0.02,
    y = seq(
      0.93,
      0.08,
      length.out = length(summary_lines)
    ),
    labels = summary_lines,
    adj = c(
      0,
      0.5
    ),
    cex = 0.78
  )

  graphics::title(
    main = "Capture-area audit summary",
    cex.main = 0.92
  )
}

# ------------------------------------------------------------------------------
# 6. Main workflow
# ------------------------------------------------------------------------------

main <- function() {
  workflow_start <- Sys.time()

  log_message(
    "============================================================"
  )

  log_message(
    "GSE234774 Step 12B-FAST started."
  )

  log_message(
    "============================================================"
  )

  # --------------------------------------------------------------------------
  # 6A. Validate required inputs and Step 12A2 release
  # --------------------------------------------------------------------------

  required_files <- c(
    matrix = MATRIX_FILE,
    barcodes = BARCODE_FILE,
    features = FEATURE_FILE,
    metadata = METADATA_FILE,
    repaired_asset_manifest =
      REPAIRED_ASSET_MANIFEST_FILE,
    repaired_report =
      REPAIRED_STEP12A2_REPORT,
    prior_section_manifest =
      PRIOR_SECTION_MANIFEST_FILE
  )

  input_exists <- file.exists(
    required_files
  )

  add_check(
    "C01",
    "input",
    "All required Step 12B-fast inputs exist.",
    all(input_exists),
    paste(
      paste0(
        names(required_files),
        "=",
        input_exists
      ),
      collapse = "; "
    ),
    "All required files TRUE",
    TRUE,
    paste(
      "Restore missing local inputs.",
      "Do not silently download or substitute files."
    )
  )

  if (!all(input_exists)) {
    stop(
      "Missing input file(s): ",
      paste(
        names(required_files)[
          !input_exists
        ],
        collapse = ", "
      )
    )
  }

  prior_report <- readLines(
    REPAIRED_STEP12A2_REPORT,
    warn = FALSE,
    encoding = "UTF-8"
  )

  prior_status_line <- prior_report[
    grepl(
      "^FINAL_STATUS:",
      prior_report
    )
  ]

  prior_status_ok <- length(
    prior_status_line
  ) == 1L &&
    grepl(
      "^FINAL_STATUS: (PASS|PASS_WITH_REVIEW)$",
      prior_status_line
    )

  add_check(
    "C02",
    "input",
    "Repaired Step 12A2 status releases Step 12B.",
    prior_status_ok,
    if (
      length(prior_status_line) > 0L
    ) {
      prior_status_line[[1L]]
    } else {
      "No FINAL_STATUS line"
    },
    "PASS or PASS_WITH_REVIEW",
    TRUE,
    "Do not continue until Step 12A2 is technically complete."
  )

  if (!prior_status_ok) {
    stop(
      "Repaired Step 12A2 status does not permit Step 12B."
    )
  }

  asset_manifest <- data.table::fread(
    REPAIRED_ASSET_MANIFEST_FILE,
    data.table = TRUE,
    showProgress = FALSE
  )

  required_columns_present(
    asset_manifest,
    c(
      "capture_area_id",
      "n_spots",
      "n_section_ids",
      "image_lowres_member",
      "position_member",
      "scalefactor_member",
      "linkage_complete"
    ),
    "Repaired asset manifest"
  )

  asset_manifest[
    ,
    capture_area_id := toupper(
      trimws(
        as.character(capture_area_id)
      )
    )
  ]

  asset_manifest[
    ,
    linkage_complete := as_logical_strict(
      linkage_complete,
      "linkage_complete"
    )
  ]

  asset_manifest_ok <- nrow(
    asset_manifest
  ) == 12L &&
    setequal(
      asset_manifest$capture_area_id,
      EXPECTED_CAPTURE_AREAS
    ) &&
    all(
      asset_manifest$linkage_complete
    )

  add_check(
    "C03",
    "input",
    "The repaired asset manifest contains 12 complete spatial 2D capture areas.",
    asset_manifest_ok,
    paste0(
      "rows=",
      nrow(asset_manifest),
      "; complete=",
      sum(asset_manifest$linkage_complete)
    ),
    "12 rows and 12 complete linkages",
    TRUE,
    "Return to Step 12A2 if an asset linkage is incomplete."
  )

  if (!asset_manifest_ok) {
    stop(
      "Repaired capture-area asset manifest is not complete."
    )
  }

  # --------------------------------------------------------------------------
  # 6B. Read metadata, barcodes, and features
  # --------------------------------------------------------------------------

  log_message(
    "[1/6] Reading metadata, barcodes, and features."
  )

  metadata <- data.table::fread(
    METADATA_FILE,
    header = TRUE,
    data.table = TRUE,
    showProgress = FALSE
  )

  setnames(
    metadata,
    trimws(
      names(metadata)
    )
  )

  required_columns_present(
    metadata,
    c(
      "barcode",
      "id",
      "slide",
      "spot",
      "group",
      "replicate",
      "sections",
      "x",
      "y",
      "z",
      "nCount_Spatial",
      "nFeature_Spatial"
    ),
    "Spatial 2D metadata"
  )

  metadata[
    ,
    barcode := trimws(
      as.character(barcode)
    )
  ]

  metadata[
    ,
    id := trimws(
      as.character(id)
    )
  ]

  metadata[
    ,
    slide := trimws(
      as.character(slide)
    )
  ]

  metadata[
    ,
    spot := toupper(
      trimws(
        as.character(spot)
      )
    )
  ]

  metadata[
    ,
    group := trimws(
      as.character(group)
    )
  ]

  metadata[
    ,
    replicate := trimws(
      as.character(replicate)
    )
  ]

  metadata[
    ,
    sections := trimws(
      as.character(sections)
    )
  ]

  metadata[
    ,
    x := suppressWarnings(
      as.numeric(x)
    )
  ]

  metadata[
    ,
    y := suppressWarnings(
      as.numeric(y)
    )
  ]

  metadata[
    ,
    z := suppressWarnings(
      as.numeric(z)
    )
  ]

  metadata[
    ,
    nCount_Spatial_source := suppressWarnings(
      as.numeric(nCount_Spatial)
    )
  ]

  metadata[
    ,
    nFeature_Spatial_source := suppressWarnings(
      as.numeric(nFeature_Spatial)
    )
  ]

  metadata[
    ,
    capture_area_id := spot
  ]

  metadata[
    ,
    section_id := id
  ]

  metadata[
    ,
    timepoint := group
  ]

  metadata[
    ,
    replicate_index := replicate
  ]

  metadata[
    ,
    biological_replicate_id := paste(
      group,
      replicate,
      sep = "_"
    )
  ]

  metadata[
    ,
    injury_condition := ifelse(
      tolower(group) == "uninjured",
      "Uninjured",
      "SCI"
    )
  ]

  metadata[
    ,
    section_order := sections
  ]

  barcodes <- data.table::fread(
    BARCODE_FILE,
    header = FALSE,
    sep = "\t",
    data.table = FALSE,
    showProgress = FALSE
  )[[1L]]

  barcodes <- trimws(
    as.character(barcodes)
  )

  features <- read_features(
    FEATURE_FILE
  )

  metadata_barcode_order_ok <-
    nrow(metadata) == length(barcodes) &&
    !anyDuplicated(metadata$barcode) &&
    !anyDuplicated(barcodes) &&
    identical(
      metadata$barcode,
      barcodes
    )

  add_check(
    "C04",
    "metadata",
    "Metadata barcode order exactly matches the matrix barcode file.",
    metadata_barcode_order_ok,
    paste0(
      "metadata rows=",
      nrow(metadata),
      "; barcode rows=",
      length(barcodes),
      "; exact order=",
      identical(
        metadata$barcode,
        barcodes
      ),
      "; duplicates=",
      anyDuplicated(barcodes)
    ),
    "Exact order and zero duplicated barcodes",
    TRUE,
    "Do not import the matrix with misaligned column names."
  )

  if (!metadata_barcode_order_ok) {
    stop(
      "Metadata and barcode order are not identical."
    )
  }

  metadata_structure_ok <-
    setequal(
      sort(
        unique(metadata$capture_area_id)
      ),
      sort(EXPECTED_CAPTURE_AREAS)
    ) &&
    setequal(
      sort(
        unique(metadata$timepoint)
      ),
      sort(EXPECTED_TIMEPOINTS)
    ) &&
    uniqueN(
      metadata$section_id
    ) == 36L

  add_check(
    "C05",
    "metadata",
    "Metadata preserves the audited 12 capture areas, 36 sections, and three timepoints.",
    metadata_structure_ok,
    paste0(
      "capture areas=",
      uniqueN(
        metadata$capture_area_id
      ),
      "; sections=",
      uniqueN(
        metadata$section_id
      ),
      "; timepoints=",
      paste(
        sort(
          unique(metadata$timepoint)
        ),
        collapse = " | "
      )
    ),
    "12 capture areas; 36 sections; uninjured, 7days, 2months",
    TRUE,
    "Return to the repaired Step 12A2 metadata mapping."
  )

  if (!metadata_structure_ok) {
    stop(
      "Spatial metadata structure differs from the repaired audit."
    )
  }

  # --------------------------------------------------------------------------
  # 6C. Import or reopen the disk-backed BPCells matrix
  # --------------------------------------------------------------------------

  log_message(
    "[2/6] Preparing the disk-backed BPCells spatial matrix."
  )

  if (
    isTRUE(
      FORCE_REIMPORT_BPCELLS
    ) &&
      dir.exists(
        FULL_BPCELLS_MATRIX_DIR
      )
  ) {
    log_message(
      "FORCE_REIMPORT_BPCELLS=TRUE; deleting the existing BPCells directory."
    )

    unlink(
      FULL_BPCELLS_MATRIX_DIR,
      recursive = TRUE,
      force = TRUE
    )
  }

  bpcells_import_start <- Sys.time()

  cache_reused <- FALSE

  if (
    dir.exists(
      FULL_BPCELLS_MATRIX_DIR
    )
  ) {
    existing_matrix <- open_bpcells_safe(
      FULL_BPCELLS_MATRIX_DIR
    )

    if (
      inherits(
        existing_matrix,
        "error"
      )
    ) {
      stop(
        "The existing BPCells directory cannot be opened: ",
        conditionMessage(existing_matrix),
        ". Set FORCE_REIMPORT_BPCELLS <- TRUE and rerun."
      )
    }

    expected_dimensions <- c(
      nrow(features),
      length(barcodes)
    )

    if (
      !identical(
        as.integer(
          dim(existing_matrix)
        ),
        as.integer(
          expected_dimensions
        )
      ) ||
        !identical(
          colnames(existing_matrix),
          barcodes
        )
    ) {
      stop(
        "The existing BPCells cache has the wrong dimensions or barcode order. ",
        "Set FORCE_REIMPORT_BPCELLS <- TRUE and rerun."
      )
    }

    full_counts <- existing_matrix
    rm(existing_matrix)
    cache_reused <- TRUE

    log_message(
      "Reusing the existing valid BPCells matrix."
    )
  } else {
    log_message(
      "No valid spatial BPCells cache exists. Starting the one-time import."
    )

    log_message(
      "The compressed mtx.gz must be scanned once. Do not interrupt RStudio."
    )

    imported_matrix <- BPCells::import_matrix_market(
      mtx_path = MATRIX_FILE,
      outdir = FULL_BPCELLS_MATRIX_DIR,
      row_names =
        features$feature_name_unique,
      col_names = barcodes,
      row_major = FALSE,
      tmpdir = TMP_DIR,
      load_bytes = BPCELLS_LOAD_BYTES,
      sort_bytes = BPCELLS_SORT_BYTES
    )

    rm(imported_matrix)
    invisible(gc())

    full_counts <- open_bpcells_safe(
      FULL_BPCELLS_MATRIX_DIR
    )

    if (
      inherits(
        full_counts,
        "error"
      )
    ) {
      stop(
        "The newly imported BPCells directory cannot be opened: ",
        conditionMessage(full_counts)
      )
    }

    log_message(
      "Completed the one-time BPCells import."
    )
  }

  bpcells_import_minutes <- elapsed_minutes(
    bpcells_import_start
  )

  bpcells_dimensions_ok <-
    identical(
      as.integer(
        dim(full_counts)
      ),
      as.integer(
        c(
          nrow(features),
          length(barcodes)
        )
      )
    ) &&
    identical(
      rownames(full_counts),
      features$feature_name_unique
    ) &&
    identical(
      colnames(full_counts),
      barcodes
    )

  add_check(
    "C06",
    "BPCells",
    "The BPCells matrix has exact feature and barcode dimensions/order.",
    bpcells_dimensions_ok,
    paste0(
      "dimensions=",
      nrow(full_counts),
      " x ",
      ncol(full_counts),
      "; cache mode=",
      if (
        cache_reused
      ) {
        "reused"
      } else {
        "new import"
      },
      "; elapsed minutes=",
      bpcells_import_minutes
    ),
    paste0(
      nrow(features),
      " x ",
      length(barcodes),
      " with exact row/column names"
    ),
    TRUE,
    "Rebuild only if the BPCells cache is incomplete or mismatched."
  )

  if (!bpcells_dimensions_ok) {
    stop(
      "BPCells matrix integrity check failed."
    )
  }

  # --------------------------------------------------------------------------
  # 6D. Compute QC metrics with BPCells
  # --------------------------------------------------------------------------

  log_message(
    "[3/6] Computing nCount, nFeature, and percent.mt from BPCells."
  )

  qc_metric_start <- Sys.time()

  ncount_bpcells <- as.numeric(
    Matrix::colSums(
      full_counts
    )
  )

  nfeature_bpcells <- as.numeric(
    Matrix::colSums(
      full_counts > 0
    )
  )

  mitochondrial_index <- grepl(
    "^mt-",
    rownames(full_counts),
    ignore.case = TRUE
  )

  mitochondrial_gene_count <- sum(
    mitochondrial_index
  )

  if (mitochondrial_gene_count > 0L) {
    mitochondrial_counts <- as.numeric(
      Matrix::colSums(
        full_counts[
          mitochondrial_index,
          ,
          drop = FALSE
        ]
      )
    )

    percent_mt <- ifelse(
      ncount_bpcells > 0,
      100 *
        mitochondrial_counts /
        ncount_bpcells,
      NA_real_
    )
  } else {
    mitochondrial_counts <- rep(
      NA_real_,
      length(ncount_bpcells)
    )

    percent_mt <- rep(
      NA_real_,
      length(ncount_bpcells)
    )
  }

  qc_metric_minutes <- elapsed_minutes(
    qc_metric_start
  )

  metadata[
    ,
    nCount_Spatial_BPCells :=
      ncount_bpcells
  ]

  metadata[
    ,
    nFeature_Spatial_BPCells :=
      nfeature_bpcells
  ]

  metadata[
    ,
    mitochondrial_counts_BPCells :=
      mitochondrial_counts
  ]

  metadata[
    ,
    percent.mt := percent_mt
  ]

  ncount_match <- all(
    metadata$nCount_Spatial_source ==
      metadata$nCount_Spatial_BPCells
  )

  nfeature_match <- all(
    metadata$nFeature_Spatial_source ==
      metadata$nFeature_Spatial_BPCells
  )

  add_check(
    "C07",
    "QC",
    "BPCells-derived nCount exactly matches the processed metadata.",
    ncount_match,
    paste0(
      "matching spots=",
      sum(
        metadata$nCount_Spatial_source ==
          metadata$nCount_Spatial_BPCells
      ),
      "/",
      nrow(metadata)
    ),
    "33,941/33,941 exact matches",
    TRUE,
    "Do not use QC metrics if matrix and metadata counts disagree."
  )

  add_check(
    "C08",
    "QC",
    "BPCells-derived nFeature exactly matches the processed metadata.",
    nfeature_match,
    paste0(
      "matching spots=",
      sum(
        metadata$nFeature_Spatial_source ==
          metadata$nFeature_Spatial_BPCells
      ),
      "/",
      nrow(metadata)
    ),
    "33,941/33,941 exact matches",
    TRUE,
    "Do not use QC metrics if matrix and metadata feature counts disagree."
  )

  if (
    !ncount_match ||
      !nfeature_match
  ) {
    stop(
      "BPCells-derived QC metrics differ from the processed metadata."
    )
  }

  zero_count_spots <- sum(
    !is.finite(
      metadata$nCount_Spatial_BPCells
    ) |
      metadata$nCount_Spatial_BPCells <= 0
  )

  add_check(
    "C09",
    "QC",
    "No processed spot has zero total counts.",
    zero_count_spots == 0L,
    paste0(
      "zero-count spots=",
      zero_count_spots
    ),
    "Zero",
    TRUE,
    "Review matrix integrity before any downstream spatial analysis."
  )

  if (mitochondrial_gene_count == 0L) {
    add_review(
      "R_MT",
      "QC",
      "No mitochondrial gene was detected from the feature names.",
      "percent.mt is unavailable.",
      "Review feature naming before defining mitochondrial QC.",
      "Do not infer mitochondrial genes from another dataset."
    )
  }

  # --------------------------------------------------------------------------
  # 6E. Derive section-aware review flags
  # --------------------------------------------------------------------------

  log_message(
    "[4/6] Deriving conservative section-aware QC review flags."
  )

  section_thresholds <- metadata[
    ,
    .(
      low_count_threshold =
        robust_lower_threshold(
          nCount_Spatial_BPCells,
          ABSOLUTE_MIN_NCOUNT
        ),
      low_feature_threshold =
        robust_lower_threshold(
          nFeature_Spatial_BPCells,
          ABSOLUTE_MIN_NFEATURE
        ),
      high_mt_threshold =
        if (
          all(
            is.na(percent.mt)
          )
        ) {
          NA_real_
        } else {
          robust_upper_mt_threshold(
            percent.mt
          )
        }
    ),
    by = section_id
  ]

  metadata <- merge(
    metadata,
    section_thresholds,
    by = "section_id",
    all.x = TRUE,
    sort = FALSE
  )

  metadata <- metadata[
    match(
      barcodes,
      barcode
    )
  ]

  if (
    anyNA(
      metadata$barcode
    ) ||
      !identical(
        metadata$barcode,
        barcodes
      )
  ) {
    stop(
      "Metadata order changed unexpectedly after attaching section thresholds."
    )
  }

  metadata[
    ,
    zero_count_flag := (
      !is.finite(
        nCount_Spatial_BPCells
      ) |
        nCount_Spatial_BPCells <= 0
    )
  ]

  metadata[
    ,
    low_count_review_flag := (
      !zero_count_flag &
        nCount_Spatial_BPCells <
        low_count_threshold
    )
  ]

  metadata[
    ,
    low_feature_review_flag := (
      !is.finite(
        nFeature_Spatial_BPCells
      ) |
        nFeature_Spatial_BPCells <
        low_feature_threshold
    )
  ]

  metadata[
    ,
    high_mt_review_flag := if (
      all(
        is.na(percent.mt)
      )
    ) {
      FALSE
    } else {
      is.finite(
        percent.mt
      ) &
        percent.mt >
        high_mt_threshold
    }
  ]

  metadata[
    ,
    qc_review_flag := (
      zero_count_flag |
        low_count_review_flag |
        low_feature_review_flag |
        high_mt_review_flag
    )
  ]

  metadata[
    ,
    qc_provisional_pass :=
      !qc_review_flag
  ]

  metadata[
    ,
    qc_review_reason := ""
  ]

  metadata[
    ,
    qc_review_reason := append_reason(
      qc_review_reason,
      zero_count_flag,
      "zero_count"
    )
  ]

  metadata[
    ,
    qc_review_reason := append_reason(
      qc_review_reason,
      low_count_review_flag,
      "low_count"
    )
  ]

  metadata[
    ,
    qc_review_reason := append_reason(
      qc_review_reason,
      low_feature_review_flag,
      "low_feature"
    )
  ]

  metadata[
    ,
    qc_review_reason := append_reason(
      qc_review_reason,
      high_mt_review_flag,
      "high_percent_mt"
    )
  ]

  # data.table compatibility:
  # Do not use DT[!logical_column, ...] here. In recent data.table versions,
  # a single logical-column symbol in i may be interpreted as column selection.
  metadata[
    qc_review_flag == FALSE,
    qc_review_reason := "none"
  ]

  section_qc <- metadata[
    ,
    .(
      capture_area_id =
        collapse_unique(
          capture_area_id
        ),
      slide_id =
        collapse_unique(
          slide
        ),
      timepoint =
        collapse_unique(
          timepoint
        ),
      replicate_index =
        collapse_unique(
          replicate_index
        ),
      biological_replicate_id =
        collapse_unique(
          biological_replicate_id
        ),
      injury_condition =
        collapse_unique(
          injury_condition
        ),
      section_order =
        collapse_unique(
          section_order
        ),
      z_coordinate =
        collapse_unique(
          z
        ),
      n_spots = .N,
      n_zero_count =
        sum(
          zero_count_flag
        ),
      n_qc_review =
        sum(
          qc_review_flag
        ),
      qc_review_fraction =
        mean(
          qc_review_flag
        ),
      n_provisional_pass =
        sum(
          qc_provisional_pass
        ),
      provisional_pass_fraction =
        mean(
          qc_provisional_pass
        ),
      nCount_min =
        min(
          nCount_Spatial_BPCells
        ),
      nCount_q01 =
        safe_quantile(
          nCount_Spatial_BPCells,
          0.01
        ),
      nCount_median =
        safe_median(
          nCount_Spatial_BPCells
        ),
      nCount_q99 =
        safe_quantile(
          nCount_Spatial_BPCells,
          0.99
        ),
      nCount_max =
        max(
          nCount_Spatial_BPCells
        ),
      nFeature_min =
        min(
          nFeature_Spatial_BPCells
        ),
      nFeature_q01 =
        safe_quantile(
          nFeature_Spatial_BPCells,
          0.01
        ),
      nFeature_median =
        safe_median(
          nFeature_Spatial_BPCells
        ),
      nFeature_q99 =
        safe_quantile(
          nFeature_Spatial_BPCells,
          0.99
        ),
      nFeature_max =
        max(
          nFeature_Spatial_BPCells
        ),
      percent_mt_median =
        safe_median(
          percent.mt
        ),
      percent_mt_q99 =
        safe_quantile(
          percent.mt,
          0.99
        ),
      low_count_threshold =
        unique(
          low_count_threshold
        ),
      low_feature_threshold =
        unique(
          low_feature_threshold
        ),
      high_mt_threshold =
        unique(
          high_mt_threshold
        )
    ),
    by = section_id
  ]

  section_qc[
    ,
    section_stop_flag :=
      n_spots <
      ABSOLUTE_MIN_SECTION_SPOTS_STOP
  ]

  section_qc[
    ,
    high_review_fraction_flag :=
      qc_review_fraction >
      SECTION_REVIEW_FLAG_FRACTION
  ]

  # --------------------------------------------------------------------------
  # 6F. Validate position/image linkage and make compact QC pages
  # --------------------------------------------------------------------------

  log_message(
    "[5/6] Validating capture-area assets and drawing compact QC pages."
  )

  capture_audit_list <- vector(
    mode = "list",
    length = length(
      EXPECTED_CAPTURE_AREAS
    )
  )

  grDevices::pdf(
    QC_PDF_OUTPUT,
    width = 12,
    height = 8.5,
    onefile = TRUE,
    useDingbats = FALSE
  )

  pdf_open <- TRUE

  on.exit(
    {
      if (isTRUE(pdf_open)) {
        try(
          grDevices::dev.off(),
          silent = TRUE
        )
      }
    },
    add = TRUE
  )

  for (
    capture_index in
      seq_along(
        EXPECTED_CAPTURE_AREAS
      )
  ) {
    capture_area <- EXPECTED_CAPTURE_AREAS[
      capture_index
    ]

    log_message(
      "[",
      capture_index,
      "/12] Reviewing ",
      capture_area,
      "."
    )

    current_manifest <- asset_manifest[
      capture_area_id == capture_area
    ]

    if (nrow(current_manifest) != 1L) {
      stop(
        "Expected one asset-manifest row for ",
        capture_area,
        "; observed ",
        nrow(current_manifest)
      )
    }

    source_lowres_image <- file.path(
      AUDIT_EXTRACTED_ASSET_DIR,
      current_manifest$image_lowres_member
    )

    source_position_file <- file.path(
      AUDIT_EXTRACTED_ASSET_DIR,
      current_manifest$position_member
    )

    source_scalefactor_file <- file.path(
      AUDIT_EXTRACTED_ASSET_DIR,
      current_manifest$scalefactor_member
    )

    capture_sources <- c(
      lowres_image =
        source_lowres_image,
      position =
        source_position_file,
      scalefactor =
        source_scalefactor_file
    )

    if (!all(file.exists(capture_sources))) {
      stop(
        "Missing Step 12A2 controlled asset(s) for ",
        capture_area,
        ": ",
        paste(
          names(capture_sources)[
            !file.exists(capture_sources)
          ],
          collapse = ", "
        )
      )
    }

    preview_capture_dir <- file.path(
      ASSET_PREVIEW_DIR,
      capture_area
    )

    dir.create(
      preview_capture_dir,
      recursive = TRUE,
      showWarnings = FALSE
    )

    preview_lowres_image <- file.path(
      preview_capture_dir,
      "tissue_lowres_image.png"
    )

    if (
      !file.exists(
        preview_lowres_image
      ) ||
        file.info(
          preview_lowres_image
        )$size <= 0
    ) {
      decompress_gzip_binary(
        source_lowres_image,
        preview_lowres_image
      )
    }

    lowres_image <- png::readPNG(
      preview_lowres_image,
      native = FALSE,
      info = TRUE
    )

    image_dimensions <- dim(
      lowres_image
    )

    if (
      length(image_dimensions) < 2L ||
        any(
          image_dimensions[1:2] <= 0L
        )
    ) {
      stop(
        "Unreadable low-resolution image for ",
        capture_area
      )
    }

    image_height <- image_dimensions[1L]
    image_width <- image_dimensions[2L]

    position_table <- read_position_file(
      source_position_file
    )

    scalefactor_text <- read_gzip_text(
      source_scalefactor_file
    )

    scalefactors <- jsonlite::fromJSON(
      scalefactor_text
    )

    required_scalefactor_fields <- c(
      "tissue_lowres_scalef",
      "spot_diameter_fullres"
    )

    missing_scalefactor_fields <- setdiff(
      required_scalefactor_fields,
      names(scalefactors)
    )

    if (
      length(
        missing_scalefactor_fields
      ) > 0L
    ) {
      stop(
        "Scalefactor JSON for ",
        capture_area,
        " is missing: ",
        paste(
          missing_scalefactor_fields,
          collapse = ", "
        )
      )
    }

    current_metadata <- metadata[
      capture_area_id ==
        capture_area
    ]

    processed_prefix <- paste0(
      capture_area,
      "-"
    )

    prefix_ok <- all(
      startsWith(
        current_metadata$barcode,
        processed_prefix
      )
    )

    if (!prefix_ok) {
      stop(
        "Processed barcode prefix mismatch for ",
        capture_area
      )
    }

    current_metadata[
      ,
      raw_spot_barcode := sub(
        paste0(
          "^",
          regex_escape(
            processed_prefix
          )
        ),
        "",
        barcode,
        perl = TRUE
      )
    ]

    if (
      anyDuplicated(
        current_metadata$raw_spot_barcode
      )
    ) {
      stop(
        "Duplicate raw barcodes after prefix removal for ",
        capture_area
      )
    }

    current_metadata <- merge(
      current_metadata,
      position_table[
        ,
        .(
          raw_spot_barcode = barcode,
          position_in_tissue = in_tissue,
          array_row,
          array_col,
          pxl_row_in_fullres,
          pxl_col_in_fullres
        )
      ],
      by = "raw_spot_barcode",
      all.x = TRUE,
      sort = FALSE
    )

    current_metadata <- current_metadata[
      match(
        metadata[
          capture_area_id ==
            capture_area,
          barcode
        ],
        barcode
      )
    ]

    position_overlap_fraction <- mean(
      !is.na(
        current_metadata$pxl_row_in_fullres
      ) &
        !is.na(
          current_metadata$pxl_col_in_fullres
        )
    )

    in_tissue_overlap_fraction <- mean(
      current_metadata$position_in_tissue ==
        1L,
      na.rm = TRUE
    )

    current_metadata[
      ,
      image_x :=
        pxl_col_in_fullres *
        as.numeric(
          scalefactors$tissue_lowres_scalef
        )
    ]

    current_metadata[
      ,
      image_y :=
        pxl_row_in_fullres *
        as.numeric(
          scalefactors$tissue_lowres_scalef
        )
    ]

    coordinate_bounds_ok <- all(
      is.finite(
        current_metadata$image_x
      ) &
        is.finite(
          current_metadata$image_y
      ) &
        current_metadata$image_x >= 0 &
        current_metadata$image_x <=
        image_width + 2 &
        current_metadata$image_y >= 0 &
        current_metadata$image_y <=
        image_height + 2
    )

    section_rows <- section_qc[
      capture_area_id ==
        capture_area
    ][
      order(
        timepoint,
        section_id
      )
    ]

    graphics::par(
      mfrow = c(
        2,
        3
      ),
      mar = c(
        0.5,
        0.5,
        2.5,
        0.5
      ),
      oma = c(
        0.5,
        0.5,
        2.2,
        0.5
      ),
      xaxs = "i",
      yaxs = "i"
    )

    draw_section_panel(
      image_array = lowres_image,
      image_width = image_width,
      image_height = image_height,
      x = current_metadata$image_x,
      y = current_metadata$image_y,
      timepoint =
        current_metadata$timepoint,
      section_id =
        current_metadata$section_id,
      title_text =
        "Condition-specific sections"
    )

    draw_continuous_spatial_panel(
      image_array = lowres_image,
      image_width = image_width,
      image_height = image_height,
      x = current_metadata$image_x,
      y = current_metadata$image_y,
      values =
        current_metadata$nCount_Spatial_BPCells,
      title_text =
        "nCount_Spatial",
      palette_name =
        "Viridis"
    )

    draw_continuous_spatial_panel(
      image_array = lowres_image,
      image_width = image_width,
      image_height = image_height,
      x = current_metadata$image_x,
      y = current_metadata$image_y,
      values =
        current_metadata$nFeature_Spatial_BPCells,
      title_text =
        "nFeature_Spatial",
      palette_name =
        "Viridis"
    )

    draw_continuous_spatial_panel(
      image_array = lowres_image,
      image_width = image_width,
      image_height = image_height,
      x = current_metadata$image_x,
      y = current_metadata$image_y,
      values =
        current_metadata$percent.mt,
      title_text =
        "percent.mt",
      palette_name =
        "Inferno"
    )

    draw_qc_flag_panel(
      image_array = lowres_image,
      image_width = image_width,
      image_height = image_height,
      x = current_metadata$image_x,
      y = current_metadata$image_y,
      qc_review_flag =
        current_metadata$qc_review_flag,
      title_text =
        "Conservative QC review flags"
    )

    draw_summary_panel(
      capture_area_id = capture_area,
      capture_spot_count =
        nrow(current_metadata),
      section_ids =
        section_rows$section_id,
      section_spot_counts =
        section_rows$n_spots,
      qc_review_count =
        sum(
          current_metadata$qc_review_flag
        ),
      qc_review_fraction =
        mean(
          current_metadata$qc_review_flag
        ),
      zero_count_spots =
        sum(
          current_metadata$zero_count_flag
        ),
      position_overlap_fraction =
        position_overlap_fraction,
      in_tissue_overlap_fraction =
        in_tissue_overlap_fraction,
      bpcells_mode = if (
        cache_reused
      ) {
        "reused cache"
      } else {
        "new cache"
      }
    )

    graphics::mtext(
      paste0(
        "GSE234774 spatial 2D | ",
        capture_area,
        " | Step 12B-FAST"
      ),
      side = 3,
      outer = TRUE,
      line = 0.4,
      cex = 1.05,
      font = 2
    )

    capture_audit_list[[capture_index]] <- data.table(
      capture_area_id = capture_area,
      processed_spots =
        nrow(current_metadata),
      expected_spots_from_manifest =
        as.integer(
          current_manifest$n_spots
        ),
      section_count =
        uniqueN(
          current_metadata$section_id
        ),
      timepoints =
        collapse_unique(
          current_metadata$timepoint
        ),
      lowres_image_file =
        source_lowres_image,
      lowres_image_width =
        image_width,
      lowres_image_height =
        image_height,
      position_file =
        source_position_file,
      position_rows =
        nrow(position_table),
      position_in_tissue_rows =
        sum(
          position_table$in_tissue == 1L,
          na.rm = TRUE
        ),
      processed_position_overlap_fraction =
        position_overlap_fraction,
      processed_in_tissue_overlap_fraction =
        in_tissue_overlap_fraction,
      scalefactor_file =
        source_scalefactor_file,
      tissue_lowres_scalef =
        as.numeric(
          scalefactors$tissue_lowres_scalef
        ),
      spot_diameter_fullres =
        as.numeric(
          scalefactors$spot_diameter_fullres
        ),
      coordinate_bounds_ok =
        coordinate_bounds_ok,
      qc_review_spots =
        sum(
          current_metadata$qc_review_flag
        ),
      qc_review_fraction =
        mean(
          current_metadata$qc_review_flag
        ),
      asset_linkage_status =
        if (
          position_overlap_fraction == 1 &&
            in_tissue_overlap_fraction == 1 &&
            coordinate_bounds_ok
        ) {
          "PASS"
        } else {
          "FAIL"
        }
    )

    rm(
      current_metadata,
      lowres_image,
      position_table,
      scalefactors
    )

    invisible(gc())
  }

  grDevices::dev.off()
  pdf_open <- FALSE

  capture_audit <- rbindlist(
    capture_audit_list,
    fill = TRUE
  )

  # --------------------------------------------------------------------------
  # 6G. Final checks and outputs
  # --------------------------------------------------------------------------

  log_message(
    "[6/6] Writing tables, handoff, checks, and final report."
  )

  capture_linkage_ok <-
    nrow(capture_audit) == 12L &&
    all(
      capture_audit$processed_spots ==
        capture_audit$expected_spots_from_manifest
    ) &&
    all(
      capture_audit$section_count == 3L
    ) &&
    all(
      capture_audit$processed_position_overlap_fraction ==
        1
    ) &&
    all(
      capture_audit$processed_in_tissue_overlap_fraction ==
        1
    ) &&
    all(
      capture_audit$coordinate_bounds_ok
    ) &&
    all(
      capture_audit$asset_linkage_status ==
        "PASS"
    )

  add_check(
    "C10",
    "spatial_asset",
    "All 12 capture areas have exact image-position-processed-spot linkage.",
    capture_linkage_ok,
    paste0(
      "complete=",
      sum(
        capture_audit$asset_linkage_status ==
          "PASS"
      ),
      "/",
      nrow(capture_audit)
    ),
    "12/12 complete capture areas",
    TRUE,
    "Do not construct final spatial objects until all capture areas link exactly."
  )

  all_spots_preserved <-
    nrow(metadata) == length(barcodes) &&
    !anyDuplicated(metadata$barcode) &&
    sum(
      capture_audit$processed_spots
    ) == length(barcodes)

  add_check(
    "C11",
    "spot",
    "All processed spatial barcodes are preserved exactly once.",
    all_spots_preserved,
    paste0(
      "spot rows=",
      nrow(metadata),
      "; capture-area sum=",
      sum(
        capture_audit$processed_spots
      ),
      "; expected=",
      length(barcodes),
      "; duplicate index=",
      anyDuplicated(
        metadata$barcode
      )
    ),
    "33,941 unique spots",
    TRUE,
    "Resolve missing or duplicated spots before any downstream analysis."
  )

  section_structure_ok <-
    nrow(section_qc) == 36L &&
    uniqueN(
      section_qc$section_id
    ) == 36L &&
    !any(
      section_qc$section_stop_flag
    )

  add_check(
    "C12",
    "section",
    "All 36 condition-specific sections are retained with sufficient spots.",
    section_structure_ok,
    paste0(
      "sections=",
      uniqueN(
        section_qc$section_id
      ),
      "; below ",
      ABSOLUTE_MIN_SECTION_SPOTS_STOP,
      " spots=",
      sum(
        section_qc$section_stop_flag
      )
    ),
    paste0(
      "36 sections; each >= ",
      ABSOLUTE_MIN_SECTION_SPOTS_STOP,
      " spots"
    ),
    TRUE,
    "Do not perform section-level inference for a missing or extremely small section."
  )

  high_review_sections <- section_qc[
    high_review_fraction_flag == TRUE
  ]

  if (nrow(high_review_sections) > 0L) {
    add_review(
      "R_QC_SECTION",
      "QC",
      "One or more sections exceed the conservative QC-review fraction.",
      paste0(
        "sections=",
        paste(
          high_review_sections$section_id,
          collapse = " | "
        ),
        "; threshold=",
        SECTION_REVIEW_FLAG_FRACTION
      ),
      "Review the maps and section-level distributions before final filtering.",
      "Do not automatically delete all review-flagged spots."
    )
  }

  add_review(
    "R_AGE",
    "metadata",
    "Age is not a dedicated processed spatial 2D metadata field.",
    "Age was not added or inferred.",
    "Confirm adult cohort status from authoritative source metadata.",
    "Do not infer age from filenames or from the snRNA cohort."
  )

  add_review(
    "R_SEVERITY_TREATMENT",
    "metadata",
    "Severity and treatment are not dedicated processed spatial 2D metadata fields.",
    "Severity and treatment were not added or inferred.",
    "Document source-level cohort selection separately.",
    "Do not manufacture severity or treatment labels."
  )

  package_versions <- data.table(
    item_type = "package",
    item = c(
      "R",
      "BPCells",
      "data.table",
      "jsonlite",
      "png"
    ),
    value = c(
      R.version.string,
      as.character(
        utils::packageVersion(
          "BPCells"
        )
      ),
      as.character(
        utils::packageVersion(
          "data.table"
        )
      ),
      as.character(
        utils::packageVersion(
          "jsonlite"
        )
      ),
      as.character(
        utils::packageVersion(
          "png"
        )
      )
    ),
    note = "Version used in Step 12B-fast"
  )

  input_file_audit <- data.table(
    item_type = "input_file",
    item = names(required_files),
    value = unname(required_files),
    note = paste0(
      "exists=",
      file.exists(
        unname(required_files)
      )
    )
  )

  bpcells_audit <- data.table(
    item_type = "BPCells",
    item = c(
      "matrix_directory",
      "cache_mode",
      "dimensions",
      "storage_order",
      "matrix_type",
      "import_or_open_minutes",
      "QC_metric_minutes",
      "directory_size_bytes",
      "directory_size_readable",
      "sort_bytes",
      "load_bytes",
      "mitochondrial_gene_count"
    ),
    value = c(
      FULL_BPCELLS_MATRIX_DIR,
      if (
        cache_reused
      ) {
        "reused_existing_cache"
      } else {
        "new_cache_created"
      },
      paste(
        dim(full_counts),
        collapse = " x "
      ),
      as.character(
        BPCells::storage_order(
          full_counts
        )
      ),
      as.character(
        BPCells::matrix_type(
          full_counts
        )
      ),
      as.character(
        bpcells_import_minutes
      ),
      as.character(
        qc_metric_minutes
      ),
      as.character(
        directory_size_bytes(
          FULL_BPCELLS_MATRIX_DIR
        )
      ),
      format_bytes(
        directory_size_bytes(
          FULL_BPCELLS_MATRIX_DIR
        )
      ),
      as.character(
        BPCELLS_SORT_BYTES
      ),
      as.character(
        BPCELLS_LOAD_BYTES
      ),
      as.character(
        mitochondrial_gene_count
      )
    ),
    note = c(
      "Persistent disk-backed count matrix",
      "First run imports; later runs reuse",
      "Features x spots",
      "BPCells storage orientation",
      "BPCells data type",
      "Elapsed minutes",
      "Elapsed minutes",
      "Persistent cache size",
      "Persistent cache size",
      "Import sorting memory",
      "Import contiguous load size",
      "Detected by ^mt- ignoring case"
    )
  )

  input_audit <- rbindlist(
    list(
      package_versions,
      input_file_audit,
      bpcells_audit
    ),
    fill = TRUE
  )

  spot_qc_output <- metadata[
    ,
    .(
      barcode,
      capture_area_id,
      section_id,
      slide_id = slide,
      timepoint,
      replicate_index,
      biological_replicate_id,
      injury_condition,
      section_order,
      x,
      y,
      z,
      nCount_Spatial_source,
      nFeature_Spatial_source,
      nCount_Spatial_BPCells,
      nFeature_Spatial_BPCells,
      mitochondrial_counts_BPCells,
      percent.mt,
      low_count_threshold,
      low_feature_threshold,
      high_mt_threshold,
      zero_count_flag,
      low_count_review_flag,
      low_feature_review_flag,
      high_mt_review_flag,
      qc_review_flag,
      qc_provisional_pass,
      qc_review_reason
    )
  ]

  bpcells_handoff_index <- data.table(
    resource = c(
      "BPCells_full_counts",
      "processed_metadata",
      "repaired_asset_manifest",
      "capture_area_QC",
      "section_QC",
      "spot_QC",
      "QC_maps",
      "lightweight_handoff_RDS"
    ),
    path = c(
      FULL_BPCELLS_MATRIX_DIR,
      METADATA_FILE,
      REPAIRED_ASSET_MANIFEST_FILE,
      CAPTURE_ASSET_AUDIT_OUTPUT,
      SECTION_QC_OUTPUT,
      SPOT_QC_OUTPUT,
      QC_PDF_OUTPUT,
      HANDOFF_RDS_OUTPUT
    ),
    description = c(
      "Persistent disk-backed raw spatial 2D counts",
      "Original processed spatial metadata",
      "Step 12A2 repaired asset linkage",
      "Capture-area image/coordinate audit",
      "36-section QC summary",
      "33,941-spot QC table",
      "One compact six-panel page per capture area",
      "Paths, metadata, QC summaries, and provenance; no embedded count matrix"
    )
  )

  lightweight_handoff <- list(
    step = "Step12B_FAST",
    completed_at = as.character(
      Sys.time()
    ),
    project_root = PROJECT_ROOT,
    matrix = list(
      type = "BPCells_disk_matrix",
      path = FULL_BPCELLS_MATRIX_DIR,
      dimensions = dim(
        full_counts
      ),
      rownames_identical_to_features = identical(
        rownames(full_counts),
        features$feature_name_unique
      ),
      colnames_identical_to_barcodes = identical(
        colnames(full_counts),
        barcodes
      ),
      cache_mode = if (
        cache_reused
      ) {
        "reused_existing_cache"
      } else {
        "new_cache_created"
      }
    ),
    metadata_roles = list(
      barcode = "barcode",
      capture_area_id = "spot",
      section_id = "id",
      timepoint = "group",
      replicate_index = "replicate",
      biological_replicate_id =
        "paste(group, replicate, sep='_')",
      section_order = "sections",
      z_coordinate = "z"
    ),
    spot_metadata_QC = spot_qc_output,
    section_QC = section_qc,
    capture_area_audit = capture_audit,
    repaired_asset_manifest =
      asset_manifest,
    boundaries = c(
      "No spot removed.",
      "No normalization or clustering.",
      "No RCTD.",
      "No target-program score.",
      "Age/severity/treatment not inferred."
    )
  )

  data.table::fwrite(
    input_audit,
    INPUT_AUDIT_OUTPUT,
    bom = TRUE
  )

  data.table::fwrite(
    capture_audit,
    CAPTURE_ASSET_AUDIT_OUTPUT,
    bom = TRUE
  )

  data.table::fwrite(
    section_qc[
      order(
        timepoint,
        biological_replicate_id,
        capture_area_id,
        section_id
      )
    ],
    SECTION_QC_OUTPUT,
    bom = TRUE
  )

  data.table::fwrite(
    spot_qc_output,
    SPOT_QC_OUTPUT,
    bom = TRUE,
    compress = "gzip"
  )

  data.table::fwrite(
    bpcells_handoff_index,
    BPCELLS_HANDOFF_OUTPUT,
    bom = TRUE
  )

  saveRDS(
    lightweight_handoff,
    HANDOFF_RDS_OUTPUT,
    compress = FALSE
  )

  critical_failures <- audit_checks[
    critical == TRUE &
      status == "FAIL"
  ]

  review_items <- audit_checks[
    status == "REVIEW"
  ]

  final_status <- if (
    nrow(
      critical_failures
    ) > 0L
  ) {
    "STOP"
  } else if (
    nrow(
      review_items
    ) > 0L
  ) {
    "PASS_WITH_REVIEW"
  } else {
    "PASS"
  }

  data.table::fwrite(
    audit_checks,
    AUDIT_CHECKS_OUTPUT,
    bom = TRUE
  )

  total_elapsed_minutes <- elapsed_minutes(
    workflow_start
  )

  report_lines <- c(
    "============================================================",
    "GSE234774 Step 12B-FAST BPCells spatial section QC",
    "============================================================",
    paste0(
      "Generated: ",
      format(
        Sys.time(),
        "%Y-%m-%d %H:%M:%S"
      )
    ),
    paste0(
      "Project root: ",
      PROJECT_ROOT
    ),
    paste0(
      "Output directory: ",
      OUTPUT_DIR
    ),
    "",
    paste0(
      "FINAL_STATUS: ",
      final_status
    ),
    "",
    "SCOPE",
    "- The spatial 2D count matrix was imported/reused as a persistent BPCells disk matrix.",
    "- Technical QC was calculated without dense conversion.",
    "- Image, coordinate, scalefactor, capture-area, and section linkage was reviewed.",
    "- One compact six-panel QC page was generated per capture area.",
    "- No spot was removed.",
    "- No 12-object raw Seurat RDS set was created.",
    "- No normalization, PCA, clustering, RCTD, target scoring, or spatial statistics were run.",
    "",
    "BPCells MATRIX",
    paste0(
      "- Cache mode: ",
      if (
        cache_reused
      ) {
        "reused_existing_cache"
      } else {
        "new_cache_created"
      }
    ),
    paste0(
      "- Matrix dimensions: ",
      nrow(full_counts),
      " features x ",
      ncol(full_counts),
      " spots"
    ),
    paste0(
      "- Import/open time: ",
      bpcells_import_minutes,
      " minutes"
    ),
    paste0(
      "- QC metric time: ",
      qc_metric_minutes,
      " minutes"
    ),
    paste0(
      "- Cache size: ",
      format_bytes(
        directory_size_bytes(
          FULL_BPCELLS_MATRIX_DIR
        )
      )
    ),
    paste0(
      "- Cache path: ",
      FULL_BPCELLS_MATRIX_DIR
    ),
    "",
    "QC CONSISTENCY",
    paste0(
      "- nCount exact match with processed metadata: ",
      sum(
        metadata$nCount_Spatial_source ==
          metadata$nCount_Spatial_BPCells
      ),
      "/",
      nrow(metadata)
    ),
    paste0(
      "- nFeature exact match with processed metadata: ",
      sum(
        metadata$nFeature_Spatial_source ==
          metadata$nFeature_Spatial_BPCells
      ),
      "/",
      nrow(metadata)
    ),
    paste0(
      "- Mitochondrial genes detected: ",
      mitochondrial_gene_count
    ),
    paste0(
      "- Zero-count spots: ",
      zero_count_spots
    ),
    "",
    "SPATIAL STRUCTURE",
    paste0(
      "- Capture areas with complete linkage: ",
      sum(
        capture_audit$asset_linkage_status ==
          "PASS"
      ),
      "/12"
    ),
    paste0(
      "- Processed spots preserved: ",
      nrow(metadata),
      "/",
      length(barcodes)
    ),
    paste0(
      "- Condition-specific sections retained: ",
      uniqueN(
        section_qc$section_id
      ),
      "/36"
    ),
    paste0(
      "- Sections below ",
      ABSOLUTE_MIN_SECTION_SPOTS_STOP,
      " spots: ",
      sum(
        section_qc$section_stop_flag
      )
    ),
    "",
    "CONSERVATIVE QC REVIEW",
    paste0(
      "- QC-review spots: ",
      sum(
        metadata$qc_review_flag
      ),
      "/",
      nrow(metadata),
      " (",
      round(
        100 *
          mean(
            metadata$qc_review_flag
          ),
        2
      ),
      "%)"
    ),
    paste0(
      "- Sections with QC-review fraction > ",
      SECTION_REVIEW_FLAG_FRACTION,
      ": ",
      nrow(
        high_review_sections
      )
    ),
    "- Review flags are not final exclusion calls.",
    "",
    "CRITICAL CHECKS",
    paste0(
      "- PASS: ",
      nrow(
        audit_checks[
          critical == TRUE &
            status == "PASS"
        ]
      ),
      "; FAIL: ",
      nrow(
        critical_failures
      )
    )
  )

  if (
    nrow(
      critical_failures
    ) > 0L
  ) {
    report_lines <- c(
      report_lines,
      "",
      "STOP REASONS",
      paste0(
        "- ",
        critical_failures$check_id,
        ": ",
        critical_failures$description,
        " Observed: ",
        critical_failures$observed,
        " Action: ",
        critical_failures$action_if_not_pass
      )
    )
  }

  if (
    nrow(
      review_items
    ) > 0L
  ) {
    report_lines <- c(
      report_lines,
      "",
      "REVIEW ITEMS",
      paste0(
        "- ",
        review_items$check_id,
        ": ",
        review_items$description,
        " Observed: ",
        review_items$observed,
        " Action: ",
        review_items$action_if_not_pass
      )
    )
  }

  report_lines <- c(
    report_lines,
    "",
    "INTERPRETATION",
    if (
      final_status == "STOP"
    ) {
      paste(
        "- Step 12B-fast is not released.",
        "Do not run RCTD or target-program scoring."
      )
    } else {
      paste(
        "- Step 12B-fast technical reconstruction is complete.",
        "The persistent BPCells matrix, 12 capture-area linkages,",
        "36 section assignments, and QC tables are available for review."
      )
    },
    "- The lightweight handoff RDS contains paths and QC metadata but does not embed the count matrix.",
    "- Age, severity, and treatment remain source-metadata review items.",
    "",
    "NEXT ACTION",
    if (
      final_status == "STOP"
    ) {
      "- Upload the seven key outputs for diagnosis."
    } else {
      paste(
        "- Upload the seven key outputs.",
        "After reviewing the section QC summary and spatial QC PDF,",
        "freeze the final spot-retention rule before constructing",
        "the filtered spatial query used for RCTD."
      )
    },
    "",
    "KEY OUTPUTS TO UPLOAD",
    "- tables/01_input_package_BPCells_audit.csv",
    "- tables/02_capture_area_image_coordinate_audit.csv",
    "- tables/03_section_QC_summary.csv",
    "- tables/04_spot_QC_flags.csv.gz",
    "- figures/01_spatial_section_and_QC_maps_fast.pdf",
    "- tables/09_STEP12B_FAST_audit_checks.csv",
    "- 10_STEP12B_FAST_FINAL_REPORT.txt",
    "",
    "DO NOT UPLOAD UNLESS REQUESTED",
    paste0(
      "- BPCells directory: ",
      FULL_BPCELLS_MATRIX_DIR
    ),
    paste0(
      "- Lightweight handoff RDS: ",
      HANDOFF_RDS_OUTPUT
    ),
    "",
    paste0(
      "TOTAL ELAPSED TIME: ",
      total_elapsed_minutes,
      " minutes"
    ),
    "============================================================"
  )

  writeLines(
    report_lines,
    FINAL_REPORT_OUTPUT,
    useBytes = TRUE
  )

  writeLines(
    capture.output(
      sessionInfo()
    ),
    SESSION_INFO_FILE,
    useBytes = TRUE
  )

  log_message(
    "FINAL_STATUS: ",
    final_status
  )

  log_message(
    "BPCells cache mode: ",
    if (
      cache_reused
    ) {
      "reused"
    } else {
      "new"
    }
  )

  log_message(
    "Capture areas: ",
    nrow(capture_audit),
    "/12"
  )

  log_message(
    "Sections: ",
    uniqueN(section_qc$section_id),
    "/36"
  )

  log_message(
    "Spots: ",
    nrow(metadata),
    "/",
    length(barcodes)
  )

  log_message(
    "Total elapsed time: ",
    total_elapsed_minutes,
    " minutes"
  )

  log_message(
    "Open first: ",
    FINAL_REPORT_OUTPUT
  )

  write_run_log()

  cat("\n")
  cat(
    "============================================================\n"
  )
  cat(
    "Step 12B-FAST completed: ",
    final_status,
    "\n",
    sep = ""
  )
  cat(
    "Capture areas: ",
    nrow(capture_audit),
    "/12\n",
    sep = ""
  )
  cat(
    "Sections: ",
    uniqueN(section_qc$section_id),
    "/36\n",
    sep = ""
  )
  cat(
    "Spots: ",
    nrow(metadata),
    "/",
    length(barcodes),
    "\n",
    sep = ""
  )
  cat(
    "Open first: ",
    FINAL_REPORT_OUTPUT,
    "\n",
    sep = ""
  )
  cat(
    "============================================================\n"
  )

  invisible(
    list(
      final_status = final_status,
      audit_checks = audit_checks,
      capture_audit = capture_audit,
      section_qc = section_qc,
      bpcells_path =
        FULL_BPCELLS_MATRIX_DIR,
      handoff_rds =
        HANDOFF_RDS_OUTPUT
    )
  )
}

tryCatch(
  main(),
  error = function(e) {
    write_fatal_report(
      conditionMessage(e)
    )

    message(
      "Step 12B-FAST stopped: ",
      conditionMessage(e)
    )
  }
)
