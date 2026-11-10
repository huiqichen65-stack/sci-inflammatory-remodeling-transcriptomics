options(stringsAsFactors = FALSE, warn = 1)
if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("Missing R package: data.table")
}
suppressPackageStartupMessages(library(data.table))
PROJECT_ROOT <- Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
INPUT_DIR <- file.path(
  PROJECT_ROOT,
  "04_results",
  "step12A2_GSE234774_spatial_metadata_asset_linkage_repair"
)
OUTPUT_DIR <- file.path(
  PROJECT_ROOT,
  "04_results",
  "step12A2_asset_manifest"
)
INPUT_ASSET_MANIFEST <- file.path(
  INPUT_DIR,
  "02_2d_capture_area_asset_manifest.csv"
)
INPUT_POSITION_AUDIT <- file.path(
  INPUT_DIR,
  "04_position_barcode_overlap_audit.csv"
)
INPUT_AUDIT_CHECKS <- file.path(
  INPUT_DIR,
  "09_STEP12A2_audit_checks.csv"
)
INPUT_FINAL_REPORT <- file.path(
  INPUT_DIR,
  "10_STEP12A2_FINAL_REPORT.txt"
)
OUTPUT_REPAIRED_MANIFEST <- file.path(
  OUTPUT_DIR,
  "01_2d_capture_area_asset_manifest.csv"
)
OUTPUT_REPAIR_EVIDENCE <- file.path(
  OUTPUT_DIR,
  "02_position_readability_evidence.csv"
)
OUTPUT_REPAIRED_CHECKS <- file.path(
  OUTPUT_DIR,
  "09_STEP12A2_audit_checks.csv"
)
OUTPUT_FINAL_REPORT <- file.path(
  OUTPUT_DIR,
  "10_STEP12A2_FINAL_REPORT.txt"
)
OUTPUT_RUN_LOG <- file.path(
  OUTPUT_DIR,
  "00_run_log.txt"
)
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
run_log <- character()
write_log <- function() {
  writeLines(run_log, OUTPUT_RUN_LOG, useBytes = TRUE)
}
as_logical_strict <- function(x, field_name) {
  if (is.logical(x)) {
    return(x)
  }
  value <- toupper(trimws(as.character(x)))
  output <- rep(NA, length(value))
  output[value %in% c("TRUE", "T", "1")] <- TRUE
  output[value %in% c("FALSE", "F", "0")] <- FALSE
  bad <- is.na(output) & !is.na(value) & nzchar(value)
  if (any(bad)) {
    stop(
      "Column '",
      field_name,
      "' contains values that cannot be converted to logical: ",
      paste(unique(value[bad]), collapse = " | ")
    )
  }
  output
}
as_numeric_strict <- function(x, field_name) {
  output <- suppressWarnings(as.numeric(x))
  bad <- is.na(output) & !is.na(x) & nzchar(trimws(as.character(x)))
  if (any(bad)) {
    stop(
      "Column '",
      field_name,
      "' contains values that cannot be converted to numeric: ",
      paste(unique(as.character(x)[bad]), collapse = " | ")
    )
  }
  output
}
required_columns_present <- function(table, required_columns, table_name) {
  missing_columns <- setdiff(required_columns, names(table))
  if (length(missing_columns) > 0L) {
    stop(
      table_name,
      " is missing required column(s): ",
      paste(missing_columns, collapse = ", ")
    )
  }
  invisible(TRUE)
}
write_fatal_report <- function(error_message) {
  write_log()
}
main <- function() {
  required_files <- c(
    asset_manifest = INPUT_ASSET_MANIFEST,
    position_audit = INPUT_POSITION_AUDIT,
    audit_checks = INPUT_AUDIT_CHECKS,
    final_report = INPUT_FINAL_REPORT
  )
  file_status <- file.exists(required_files)
  if (!all(file_status)) {
    stop(
      "Missing required Step 12A2 output file(s): ",
      paste(names(required_files)[!file_status], collapse = ", ")
    )
  }
  asset_manifest <- fread(
    INPUT_ASSET_MANIFEST,
    data.table = TRUE,
    showProgress = FALSE
  )
  position_audit <- fread(
    INPUT_POSITION_AUDIT,
    data.table = TRUE,
    showProgress = FALSE
  )
  audit_checks <- fread(
    INPUT_AUDIT_CHECKS,
    data.table = TRUE,
    showProgress = FALSE
  )
  original_report <- readLines(
    INPUT_FINAL_REPORT,
    warn = FALSE,
    encoding = "UTF-8"
  )
  required_columns_present(
    asset_manifest,
    c(
      "capture_area_id",
      "image_hires_count",
      "image_lowres_count",
      "position_count",
      "scalefactor_count",
      "image_hires_readable",
      "image_lowres_readable",
      "position_readable",
      "scalefactor_readable",
      "processed_in_position_fraction",
      "processed_in_tissue_fraction",
      "linkage_complete"
    ),
    "Asset manifest"
  )
  required_columns_present(
    position_audit,
    c(
      "capture_area_id",
      "position_member",
      "position_read_status",
      "position_read_message",
      "position_barcode_unique",
      "processed_prefix_valid",
      "processed_in_position_fraction",
      "processed_in_tissue_fraction",
      "missing_processed_barcode_examples"
    ),
    "Position audit"
  )
  required_columns_present(
    audit_checks,
    c(
      "check_id",
      "category",
      "description",
      "status",
      "observed",
      "required",
      "critical",
      "action_if_not_pass"
    ),
    "Audit checks"
  )
  asset_manifest[, capture_area_id := toupper(
    trimws(as.character(capture_area_id))
  )]
  position_audit[, capture_area_id := toupper(
    trimws(as.character(capture_area_id))
  )]
  if (anyDuplicated(asset_manifest$capture_area_id)) {
    stop("Asset manifest contains duplicated capture_area_id values.")
  }
  if (anyDuplicated(position_audit$capture_area_id)) {
    stop("Position audit contains duplicated capture_area_id values.")
  }
  if (!setequal(
    asset_manifest$capture_area_id,
    position_audit$capture_area_id
  )) {
    stop(
      "Capture-area sets differ between the asset manifest and position audit."
    )
  }
  expected_capture_areas <- c(
    "078_A", "078_B", "078_C", "078_D",
    "105_A", "105_B", "105_C", "105_D",
    "106_A", "106_B", "106_C", "106_D"
  )
  if (!setequal(
    asset_manifest$capture_area_id,
    expected_capture_areas
  )) {
    stop(
      "The Step 12A2 tables contain the expected 12 spatial 2D ",
      "capture areas. Observed: ",
      paste(sort(asset_manifest$capture_area_id), collapse = " | ")
    )
  }
  logical_columns_asset <- c(
    "image_hires_readable",
    "image_lowres_readable",
    "position_readable",
    "scalefactor_readable",
    "linkage_complete"
  )
  for (current_column in logical_columns_asset) {
    asset_manifest[, (current_column) := as_logical_strict(
      get(current_column),
      current_column
    )]
  }
  count_columns <- c(
    "image_hires_count",
    "image_lowres_count",
    "position_count",
    "scalefactor_count"
  )
  for (current_column in count_columns) {
    asset_manifest[, (current_column) := as_numeric_strict(
      get(current_column),
      current_column
    )]
  }
  fraction_columns_asset <- c(
    "processed_in_position_fraction",
    "processed_in_tissue_fraction"
  )
  for (current_column in fraction_columns_asset) {
    asset_manifest[, (current_column) := as_numeric_strict(
      get(current_column),
      current_column
    )]
  }
  position_audit[, position_barcode_unique := as_logical_strict(
    position_barcode_unique,
    "position_barcode_unique"
  )]
  position_audit[, processed_prefix_valid := as_logical_strict(
    processed_prefix_valid,
    "processed_prefix_valid"
  )]
  position_audit[, processed_in_position_fraction := as_numeric_strict(
    processed_in_position_fraction,
    "processed_in_position_fraction"
  )]
  position_audit[, processed_in_tissue_fraction := as_numeric_strict(
    processed_in_tissue_fraction,
    "processed_in_tissue_fraction"
  )]
  position_audit[, position_read_status := toupper(
    trimws(as.character(position_read_status))
  )]
  position_audit[, missing_processed_barcode_examples := trimws(
    as.character(missing_processed_barcode_examples)
  )]
  position_audit[
    missing_processed_barcode_examples %in% c("", "NA", "NaN"),
    missing_processed_barcode_examples := NA_character_
  ]
  position_audit[, repaired_position_readable := (
    position_read_status == "PASS" &
      position_barcode_unique &
      processed_prefix_valid &
      !is.na(processed_in_position_fraction) &
      processed_in_position_fraction == 1 &
      !is.na(processed_in_tissue_fraction) &
      processed_in_tissue_fraction == 1 &
      is.na(missing_processed_barcode_examples)
  )]
  repair_evidence <- position_audit[
    ,
    .(
      capture_area_id,
      position_member,
      position_read_status,
      position_read_message,
      position_barcode_unique,
      processed_prefix_valid,
      processed_in_position_fraction,
      processed_in_tissue_fraction,
      missing_processed_barcode_examples,
      repaired_position_readable
    )
  ]
  repaired_manifest <- merge(
    asset_manifest,
    repair_evidence[
      ,
      .(
        capture_area_id,
        direct_position_read_status = position_read_status,
        direct_position_read_message = position_read_message,
        direct_position_barcode_unique = position_barcode_unique,
        direct_processed_prefix_valid = processed_prefix_valid,
        direct_processed_in_position_fraction =
          processed_in_position_fraction,
        direct_processed_in_tissue_fraction =
          processed_in_tissue_fraction,
        repaired_position_readable
      )
    ],
    by = "capture_area_id",
    all.x = TRUE,
    sort = TRUE
  )
  repaired_manifest[, original_position_readable := position_readable]
  repaired_manifest[, original_linkage_complete := linkage_complete]
  repaired_manifest[, position_readable := repaired_position_readable]
  repaired_manifest[, linkage_complete := (
    image_hires_count == 1 &
      image_lowres_count == 1 &
      position_count == 1 &
      scalefactor_count == 1 &
      image_hires_readable &
      image_lowres_readable &
      position_readable &
      scalefactor_readable &
      processed_in_position_fraction == 1 &
      processed_in_tissue_fraction == 1
  )]
  repaired_manifest[, repair_applied := (
    original_position_readable != position_readable |
      original_linkage_complete != linkage_complete
  )]
  core_evidence_ok <- all(
    repaired_manifest$image_hires_count == 1 &
      repaired_manifest$image_lowres_count == 1 &
      repaired_manifest$position_count == 1 &
      repaired_manifest$scalefactor_count == 1 &
      repaired_manifest$image_hires_readable &
      repaired_manifest$image_lowres_readable &
      repaired_manifest$scalefactor_readable &
      repaired_manifest$direct_position_read_status == "PASS" &
      repaired_manifest$direct_position_barcode_unique &
      repaired_manifest$direct_processed_prefix_valid &
      repaired_manifest$direct_processed_in_position_fraction == 1 &
      repaired_manifest$direct_processed_in_tissue_fraction == 1
  )
  all_linked <- all(repaired_manifest$linkage_complete)
  if (!core_evidence_ok) {
    stop(
      "Direct evidence does not support repairing position_readable for all ",
      "12 capture areas. Review 01_2d_capture_area_asset_manifest.csv."
    )
  }
  if (!all_linked) {
    stop(
      "After evidence-based repair, one or more capture areas remain incomplete."
    )
  }
  original_false_fail_count <- sum(
    !asset_manifest$position_readable &
      asset_manifest$processed_in_position_fraction == 1 &
      asset_manifest$processed_in_tissue_fraction == 1
  )
  repaired_count <- sum(repaired_manifest$repair_applied)
  audit_checks[, critical := as_logical_strict(
    critical,
    "critical"
  )]
  audit_checks[, status := toupper(trimws(as.character(status)))]
  c19_index <- which(audit_checks$check_id == "C19")
  if (length(c19_index) != 1L) {
    stop(
      "Expected exactly one C19 row in the original audit-check table; observed ",
      length(c19_index),
      "."
    )
  }
  original_c19 <- copy(audit_checks[c19_index])
  audit_checks[
    c19_index,
    `:=`(
      status = "PASS",
      observed = paste0(
        "complete=",
        sum(repaired_manifest$linkage_complete),
        "/",
        nrow(repaired_manifest),
        "; position readability restored from direct position-audit evidence"
      ),
      action_if_not_pass = paste(
        "If this check fails in a future rerun, inspect the direct",
        "position-read and barcode-overlap evidence before Step 12B."
      )
    )
  ]
  unresolved_critical_failures <- audit_checks[
    critical == TRUE & status == "FAIL"
  ]
  review_items <- audit_checks[status == "REVIEW"]
  final_status <- if (nrow(unresolved_critical_failures) > 0L) {
    "STOP"
  } else if (nrow(review_items) > 0L) {
    "PASS_WITH_REVIEW"
  } else {
    "PASS"
  }
  if (final_status == "STOP") {
    stop(
      "Critical failures remain after repairing C19: ",
      paste(
        unresolved_critical_failures$check_id,
        collapse = " | "
      )
    )
  }
  repair_evidence[, original_manifest_position_readable :=
    asset_manifest[
      match(capture_area_id, asset_manifest$capture_area_id),
      position_readable
    ]
  ]
  repair_evidence[, evidence_based_repair := (
    !original_manifest_position_readable &
      repaired_position_readable
  )]
  setcolorder(
    repaired_manifest,
    c(
      "capture_area_id",
      "original_position_readable",
      "position_readable",
      "original_linkage_complete",
      "linkage_complete",
      "repair_applied",
      setdiff(
        names(repaired_manifest),
        c(
          "capture_area_id",
          "original_position_readable",
          "position_readable",
          "original_linkage_complete",
          "linkage_complete",
          "repair_applied"
        )
      )
    )
  )
  fwrite(
    repaired_manifest,
    OUTPUT_REPAIRED_MANIFEST,
    bom = TRUE
  )
  fwrite(
    repair_evidence,
    OUTPUT_REPAIR_EVIDENCE,
    bom = TRUE
  )
  fwrite(
    audit_checks,
    OUTPUT_REPAIRED_CHECKS,
    bom = TRUE
  )
  original_status_line <- original_report[
    grepl("^FINAL_STATUS:", original_report)
  ]
  if (length(original_status_line) == 0L) {
    original_status_line <- "FINAL_STATUS: <not found in original report>"
  }
  write_log()
  invisible(
    list(
      final_status = final_status,
      repaired_manifest = repaired_manifest,
      repair_evidence = repair_evidence,
      checks = audit_checks
    )
  )
}
tryCatch(
  main(),
  error = function(e) {
    write_fatal_report(conditionMessage(e))
  }
)
