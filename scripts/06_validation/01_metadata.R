options(stringsAsFactors = FALSE, warn = 1)
project_root <- Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
input_dir <- file.path(project_root, "01_raw_data", "GSE234774")
output_dir <- file.path( project_root, "04_results", "step9A2_GSE234774_targeted_metadata_resolution" )
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
run_log <- file.path(output_dir, "00_run_log.txt")
if (file.exists(run_log)) {
  file.remove(run_log)
}
write_csv_safe <- function(x, filename) {
  out <- file.path(output_dir, filename)
  utils::write.csv(x, out, row.names = FALSE, na = "")
  invisible(out)
}
trim_na <- function(x) {
  x <- trimws(as.character(x))
  x[x == ""] <- NA_character_
  x
}
first_nonempty <- function(x) {
  x <- trim_na(x)
  x <- x[!is.na(x)]
  if (length(x) == 0L) {
    return(NA_character_)
  }
  x[[1L]]
}
collapse_unique <- function(x, sep = " | ", max_values = 30L) {
  x <- unique(trim_na(x))
  x <- x[!is.na(x)]
  if (length(x) == 0L) {
    return(NA_character_)
  }
  if (length(x) > max_values) {
    x <- c( x[seq_len(max_values)], paste0("... +", length(x) - max_values, " more") )
  }
  paste(x, collapse = sep)
}
normalise_name <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  x
}
read_delimited_gz <- function(path) {
  con <- gzfile(path, open = "rt")
  on.exit(try(close(con), silent = TRUE), add = TRUE)
  x <- utils::read.delim( con, header = TRUE, sep = "\t", quote = "", comment.char = "", stringsAsFactors = FALSE, check.names = FALSE )
  names(x) <- make.unique(names(x), sep = "_dup")
  x
}
read_lines_gz <- function(path) {
  con <- gzfile(path, open = "rt")
  on.exit(try(close(con), silent = TRUE), add = TRUE)
  readLines(con, warn = FALSE, encoding = "UTF-8")
}
find_column <- function(df, aliases) {
  if (ncol(df) == 0L) {
    return(NA_character_)
  }
  actual <- names(df)
  normalised_actual <- normalise_name(actual)
  normalised_aliases <- normalise_name(aliases)
  exact_idx <- match(normalised_aliases, normalised_actual)
  exact_idx <- exact_idx[!is.na(exact_idx)]
  if (length(exact_idx) > 0L) {
    return(actual[[exact_idx[[1L]]]])
  }
  NA_character_
}
get_column_or_na <- function(df, column_name) {
  if (is.na(column_name) || !(column_name %in% names(df))) {
    return(rep(NA_character_, nrow(df)))
  }
  trim_na(df[[column_name]])
}
safe_match_fraction <- function(x, y) {
  x <- trim_na(x)
  y <- trim_na(y)
  keep <- !is.na(x) & !is.na(y)
  if (!any(keep)) {
    return(NA_real_)
  }
  mean(tolower(x[keep]) == tolower(y[keep]))
}
standardise_timepoint <- function(experiment, condition_label) {
  experiment <- tolower(trim_na(experiment))
  condition_label <- tolower(trim_na(condition_label))
  out <- rep(NA_character_, length(experiment))
  uninjured_idx <- experiment == "uninjured" | condition_label == "uninjured"
  out[uninjured_idx] <- "Uninjured"
  timecourse_idx <- experiment == "timecourse"
  out[timecourse_idx & condition_label == "1d"] <- "1 dpi"
  out[timecourse_idx & condition_label == "4d"] <- "4 dpi"
  out[timecourse_idx & condition_label == "7d"] <- "7 dpi"
  out[timecourse_idx & condition_label == "14d"] <- "14 dpi"
  out[timecourse_idx & condition_label == "1m"] <- "1 month"
  out[timecourse_idx & condition_label == "2m"] <- "2 months"
  out[timecourse_idx & condition_label == "old"] <- "OLD_UNRESOLVED"
  out
}
timepoint_order <- function(x) {
  levels <- c( "Uninjured", "1 dpi", "4 dpi", "7 dpi", "14 dpi", "1 month", "2 months", "OLD_UNRESOLVED" )
  match(x, levels)
}
meta_path <- file.path(input_dir, "GSE234774_rnaseq_meta.txt.gz")
barcode_path <- file.path(input_dir, "GSE234774_rnaseq_barcodes.txt.gz")
required_paths <- c(meta_path, barcode_path)
if (!all(file.exists(required_paths))) {
  missing_paths <- required_paths[!file.exists(required_paths)]
  stop("Input validation failed.", call. = FALSE)
}
meta <- read_delimited_gz(meta_path)
expression_barcodes <- read_lines_gz(barcode_path)
expression_barcodes <- trimws(expression_barcodes)
expression_barcodes <- expression_barcodes[expression_barcodes != ""]
metadata_columns <- data.frame( column_index = seq_along(meta), column_name = names(meta), normalised_name = normalise_name(names(meta)), class = vapply(meta, function(x) paste(class(x), collapse = "/"), character(1)), n_nonmissing = vapply( meta, function(x) sum(!is.na(trim_na(x))), integer(1) ), n_unique_nonmissing = vapply( meta, function(x) length(unique(trim_na(x)[!is.na(trim_na(x))])), integer(1) ), example_values = vapply( meta, function(x) collapse_unique(head(trim_na(x)[!is.na(trim_na(x))], 12L)), character(1) ), stringsAsFactors = FALSE )
write_csv_safe(metadata_columns, "01_metadata_columns.csv")
key_columns <- data.frame( role = c( "barcode", "base_cell_barcode", "experimental_block", "replicate", "experiment", "condition_label", "library", "cluster", "annotation" ), detected_column = c( find_column(meta, c("barcode", "cell_barcode")), find_column(meta, c("bc", "base_barcode", "cell_bc")), find_column(meta, c("orig.ident", "orig_ident")), find_column(meta, c("replicate", "rep", "replicate_id")), find_column(meta, c("experiment", "experiment_clean")), find_column(meta, c("label_clean", "label", "condition_clean")), find_column(meta, c("library_clean", "library", "sample_id", "sample")), find_column(meta, c("seurat_clusters", "cluster")), find_column( meta, c( "celltype", "cell_type", "annotation", "celltype_clean", "cell_type_clean", "cellclass", "cell_class" ) ) ), stringsAsFactors = FALSE )
write_csv_safe(key_columns, "02_detected_key_columns.csv")
barcode_col <- first_nonempty( key_columns$detected_column[key_columns$role == "barcode"] )
base_barcode_col <- first_nonempty( key_columns$detected_column[key_columns$role == "base_cell_barcode"] )
orig_ident_col <- first_nonempty( key_columns$detected_column[key_columns$role == "experimental_block"] )
replicate_col <- first_nonempty( key_columns$detected_column[key_columns$role == "replicate"] )
experiment_col <- first_nonempty( key_columns$detected_column[key_columns$role == "experiment"] )
label_col <- first_nonempty( key_columns$detected_column[key_columns$role == "condition_label"] )
library_col <- first_nonempty( key_columns$detected_column[key_columns$role == "library"] )
cluster_col <- first_nonempty( key_columns$detected_column[key_columns$role == "cluster"] )
annotation_col <- first_nonempty( key_columns$detected_column[key_columns$role == "annotation"] )
if (is.na(barcode_col)) {
  stop("Input validation failed.", call. = FALSE)
}
barcode_meta <- get_column_or_na(meta, barcode_col)
library_from_barcode <- sub( "-[ACGTN]+(?:-[0-9]+)?$", "", barcode_meta, perl = TRUE, ignore.case = TRUE )
base_barcode_from_barcode <- sub( "^.*-([ACGTN]+(?:-[0-9]+)?)$", "\\1", barcode_meta, perl = TRUE, ignore.case = TRUE )
experiment_from_barcode <- sub( "_.*$", "", library_from_barcode, perl = TRUE )
replicate_from_barcode <- sub( "^.*_([0-9]+)$", "\\1", library_from_barcode, perl = TRUE )
replicate_pattern_valid <- grepl( "_[0-9]+$", library_from_barcode, perl = TRUE )
replicate_from_barcode[!replicate_pattern_valid] <- NA_character_
condition_from_barcode <- library_from_barcode
condition_from_barcode <- sub( "^[^_]+_", "", condition_from_barcode, perl = TRUE )
condition_from_barcode <- sub( "_[0-9]+$", "", condition_from_barcode, perl = TRUE )
condition_from_barcode[ tolower(experiment_from_barcode) == "uninjured" & is.na(trim_na(condition_from_barcode)) ] <- "uninjured"
metadata_barcode_overlap <- mean(barcode_meta %in% expression_barcodes)
expression_barcode_overlap <- mean(expression_barcodes %in% barcode_meta)
barcode_order_identical <- identical(barcode_meta, expression_barcodes)
base_barcode_meta <- get_column_or_na(meta, base_barcode_col)
orig_ident_meta <- get_column_or_na(meta, orig_ident_col)
replicate_meta <- get_column_or_na(meta, replicate_col)
experiment_meta <- get_column_or_na(meta, experiment_col)
label_meta <- get_column_or_na(meta, label_col)
library_meta <- get_column_or_na(meta, library_col)
comparison_table <- data.frame( comparison = c( "metadata barcode vs expression barcode", "parsed base barcode vs metadata bc", "parsed experimental block vs orig.ident", "parsed experiment vs experiment column", "parsed condition vs label_clean", "parsed replicate vs replicate column", "parsed library vs library_clean" ), metadata_column = c( barcode_col, base_barcode_col, orig_ident_col, experiment_col, label_col, replicate_col, library_col ), comparison_available = c( TRUE, !is.na(base_barcode_col), !is.na(orig_ident_col), !is.na(experiment_col), !is.na(label_col), !is.na(replicate_col), !is.na(library_col) ), match_fraction = c( mean(c(metadata_barcode_overlap, expression_barcode_overlap)), safe_match_fraction(base_barcode_from_barcode, base_barcode_meta), safe_match_fraction(experiment_from_barcode, orig_ident_meta), safe_match_fraction(experiment_from_barcode, experiment_meta), safe_match_fraction(condition_from_barcode, label_meta), safe_match_fraction(replicate_from_barcode, replicate_meta), safe_match_fraction(library_from_barcode, library_meta) ), detail = c( paste0( "metadata->expression=", round(100 * metadata_barcode_overlap, 3), "%; expression->metadata=", round(100 * expression_barcode_overlap, 3), "%; order_identical=", barcode_order_identical ), "The base barcode excludes the library prefix.", "orig.ident is expected to represent the experimental block.", "Compared only when an explicit experiment column exists.", "Compared only when label_clean or an equivalent column exists.", "Compared only when a replicate column exists.", "Compared only when library_clean or an equivalent column exists." ), stringsAsFactors = FALSE )
write_csv_safe( comparison_table, "03_barcode_and_metadata_consistency.csv" )
library_meta_match <- safe_match_fraction( library_from_barcode, library_meta )
if ( !is.na(library_col) && !is.na(library_meta_match) && library_meta_match >= 0.95 ) {
  biological_library <- library_meta
  biological_library_source <- paste0( library_col, " (confirmed against barcode prefix)" )
} else {
  biological_library <- library_from_barcode
  biological_library_source <- "parsed barcode prefix"
}
experiment_final <- experiment_from_barcode
condition_final <- condition_from_barcode
replicate_final <- replicate_from_barcode
experiment_meta_match <- safe_match_fraction( experiment_from_barcode, experiment_meta )
label_meta_match <- safe_match_fraction( condition_from_barcode, label_meta )
replicate_meta_match <- safe_match_fraction( replicate_from_barcode, replicate_meta )
if ( !is.na(experiment_col) && !is.na(experiment_meta_match) && experiment_meta_match >= 0.95 ) {
  experiment_final <- experiment_meta
}
if ( !is.na(label_col) && !is.na(label_meta_match) && label_meta_match >= 0.95 ) {
  condition_final <- label_meta
}
if ( !is.na(replicate_col) && !is.na(replicate_meta_match) && replicate_meta_match >= 0.95 ) {
  replicate_final <- replicate_meta
}
experiment_final <- tolower(trim_na(experiment_final))
condition_final <- tolower(trim_na(condition_final))
replicate_final <- trim_na(replicate_final)
timepoint_standardised <- standardise_timepoint( experiment = experiment_final, condition_label = condition_final )
candidate_main_cohort <- ( experiment_final == "uninjured" | condition_final == "uninjured" | ( experiment_final == "timecourse" & condition_final %in% c("1d", "4d", "7d", "14d", "1m", "2m") ) )
candidate_exclusion_reason <- ifelse( candidate_main_cohort, NA_character_, ifelse( experiment_final == "timecourse" & condition_final == "old", "timecourse_OLD requires manual interpretation and is not auto-included", paste0( "separate experimental block: ", experiment_final ) ) )
cell_manifest <- data.frame( barcode = barcode_meta, biological_library = biological_library, biological_library_source = biological_library_source, experiment = experiment_final, condition_label = condition_final, replicate = replicate_final, timepoint_standardised = timepoint_standardised, candidate_main_timecourse_cohort = candidate_main_cohort, candidate_exclusion_reason = candidate_exclusion_reason, stringsAsFactors = FALSE )
preview_idx <- unique(c( seq_len(min(200L, nrow(cell_manifest))), which( is.na(cell_manifest$biological_library) | is.na(cell_manifest$experiment) | is.na(cell_manifest$replicate) ) ))
cell_manifest_preview <- cell_manifest[ preview_idx, , drop = FALSE ]
write_csv_safe( cell_manifest_preview, "04_cell_manifest_preview_and_unresolved.csv" )
split_idx <- split( seq_len(nrow(cell_manifest)), cell_manifest$biological_library )
library_manifest_rows <- lapply(names(split_idx), function(lib) {
  idx <- split_idx[[lib]]
  data.frame( biological_library = lib, n_nuclei = length(idx), experiment = collapse_unique(cell_manifest$experiment[idx]), condition_label = collapse_unique(cell_manifest$condition_label[idx]), replicate = collapse_unique(cell_manifest$replicate[idx]), timepoint_standardised = collapse_unique( cell_manifest$timepoint_standardised[idx] ), candidate_main_timecourse_cohort = all( cell_manifest$candidate_main_timecourse_cohort[idx] ), candidate_exclusion_reason = collapse_unique( cell_manifest$candidate_exclusion_reason[idx] ), n_unique_cell_barcodes = length(unique(cell_manifest$barcode[idx])), duplicated_cell_barcodes_within_library = sum( duplicated(cell_manifest$barcode[idx]) ), stringsAsFactors = FALSE )
})
library_manifest <- do.call(rbind, library_manifest_rows)
rownames(library_manifest) <- NULL
library_manifest <- library_manifest[ order( timepoint_order(library_manifest$timepoint_standardised), library_manifest$experiment, library_manifest$condition_label, suppressWarnings(as.numeric(library_manifest$replicate)), library_manifest$biological_library, na.last = TRUE ), , drop = FALSE ]
write_csv_safe( library_manifest, "05_library_manifest_resolved.csv" )
candidate_library_manifest <- library_manifest[ library_manifest$candidate_main_timecourse_cohort, , drop = FALSE ]
write_csv_safe( candidate_library_manifest, "06_candidate_timecourse_cohort_manifest.csv" )
if (nrow(candidate_library_manifest) > 0L) {
  timepoint_split <- split( seq_len(nrow(candidate_library_manifest)), candidate_library_manifest$timepoint_standardised )
  candidate_timepoint_summary <- do.call( rbind, lapply(names(timepoint_split), function(tp) {
      idx <- timepoint_split[[tp]]
      data.frame( timepoint = tp, n_biological_libraries = length(idx), n_nuclei = sum(candidate_library_manifest$n_nuclei[idx]), library_ids = paste( candidate_library_manifest$biological_library[idx], collapse = ";" ), replicate_labels = collapse_unique( candidate_library_manifest$replicate[idx] ), stringsAsFactors = FALSE )
    })
  )
  candidate_timepoint_summary <- candidate_timepoint_summary[ order(timepoint_order(candidate_timepoint_summary$timepoint)), , drop = FALSE ]
} else {
  candidate_timepoint_summary <- data.frame( timepoint = character(0), n_biological_libraries = integer(0), n_nuclei = integer(0), library_ids = character(0), replicate_labels = character(0), stringsAsFactors = FALSE )
}
write_csv_safe( candidate_timepoint_summary, "06B_candidate_timepoint_summary.csv" )
experiment_split <- split( seq_len(nrow(library_manifest)), library_manifest$experiment )
experimental_block_summary <- do.call( rbind, lapply(names(experiment_split), function(exp_name) {
    idx <- experiment_split[[exp_name]]
    data.frame( experimental_block = exp_name, n_biological_libraries = length(idx), n_nuclei = sum(library_manifest$n_nuclei[idx]), condition_labels = collapse_unique( library_manifest$condition_label[idx] ), replicate_labels = collapse_unique( library_manifest$replicate[idx] ), included_in_candidate_main_cohort = any( library_manifest$candidate_main_timecourse_cohort[idx] ), stringsAsFactors = FALSE )
  })
)
experimental_block_summary <- experimental_block_summary[ order(experimental_block_summary$experimental_block), , drop = FALSE ]
write_csv_safe( experimental_block_summary, "06C_experimental_block_summary.csv" )
barcode_count_match <- nrow(meta) == length(expression_barcodes)
barcode_correspondence_pass <- metadata_barcode_overlap >= 0.99 && expression_barcode_overlap >= 0.99
library_resolved_fraction <- mean( !is.na(trim_na(cell_manifest$biological_library)) )
experiment_resolved_fraction <- mean( !is.na(trim_na(cell_manifest$experiment)) )
replicate_resolved_fraction <- mean( !is.na(trim_na(cell_manifest$replicate)) )
n_libraries <- nrow(library_manifest)
n_candidate_libraries <- nrow(candidate_library_manifest)
expected_main_timepoints <- c( "Uninjured", "1 dpi", "4 dpi", "7 dpi", "14 dpi", "1 month", "2 months" )
observed_main_timepoints <- unique( candidate_timepoint_summary$timepoint )
all_main_timepoints_present <- all( expected_main_timepoints %in% observed_main_timepoints )
minimum_two_libraries_each <- if ( all_main_timepoints_present ) {
  all( candidate_timepoint_summary$n_biological_libraries[ match( expected_main_timepoints, candidate_timepoint_summary$timepoint ) ] >= 2L )
} else {
  FALSE
}
late_timepoints_present <- all( c("14 dpi", "1 month", "2 months") %in% observed_main_timepoints )
old_libraries <- library_manifest[ tolower(library_manifest$experiment) == "timecourse" & tolower(library_manifest$condition_label) == "old", , drop = FALSE ]
audit_checks <- data.frame( check = c( "Metadata row count equals expression barcode count", "Metadata barcodes and expression barcodes correspond bidirectionally at >=99%", "Biological library ID is resolved for >=99% of nuclei", "Experimental block is resolved for >=99% of nuclei", "Replicate label is resolved for >=99% of nuclei", "At least two biological libraries are identified", "Candidate uninjured plus time-course cohort is identified", "Expected main timepoints are present", "Each expected main timepoint contains at least two libraries", "Late validation timepoints 14 dpi, 1 month, and 2 months are present", "timecourse_OLD libraries are kept separate pending interpretation", "Age of the candidate time-course cohort is directly encoded in cell metadata" ), required_for_schema_pass = c( TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE ), status = c( ifelse(barcode_count_match, "PASS", "STOP"), ifelse(barcode_correspondence_pass, "PASS", "STOP"), ifelse(library_resolved_fraction >= 0.99, "PASS", "STOP"), ifelse(experiment_resolved_fraction >= 0.99, "PASS", "STOP"), ifelse(replicate_resolved_fraction >= 0.99, "PASS", "STOP"), ifelse(n_libraries >= 2L, "PASS", "STOP"), ifelse(n_candidate_libraries >= 2L, "PASS", "STOP"), ifelse(all_main_timepoints_present, "PASS", "STOP"), ifelse(minimum_two_libraries_each, "PASS", "STOP"), ifelse(late_timepoints_present, "PASS", "STOP"), ifelse(nrow(old_libraries) > 0L, "REVIEW", "NOT_APPLICABLE"), "REVIEW" ), detail = c( paste0( "metadata rows=", nrow(meta), "; expression barcodes=", length(expression_barcodes) ), paste0( "metadata->expression=", round(100 * metadata_barcode_overlap, 3), "%; expression->metadata=", round(100 * expression_barcode_overlap, 3), "%; order_identical=", barcode_order_identical ), paste0( round(100 * library_resolved_fraction, 3), "%; source=", biological_library_source ), paste0( round(100 * experiment_resolved_fraction, 3), "%" ), paste0( round(100 * replicate_resolved_fraction, 3), "%" ), paste0( "resolved biological libraries=", n_libraries ), paste0( "candidate libraries=", n_candidate_libraries, "; candidate nuclei=", sum(candidate_library_manifest$n_nuclei) ), paste0( "observed=", paste( candidate_timepoint_summary$timepoint, collapse = ";" ) ), paste0( "library counts=", paste( paste0( candidate_timepoint_summary$timepoint, ":", candidate_timepoint_summary$n_biological_libraries ), collapse = ";" ) ), paste0( "late timepoints present=", late_timepoints_present ), paste0( "OLD libraries=", nrow(old_libraries), "; not auto-included" ), paste( "Age is not accepted from barcode text or experimental-block labels.", "Confirm the time-course cohort age from the study methods before final freeze." ) ), stringsAsFactors = FALSE )
write_csv_safe( audit_checks, "07_STEP9A2_audit_checks.csv" )
schema_fail <- audit_checks$required_for_schema_pass & audit_checks$status != "PASS"
if (any(schema_fail)) {
  final_decision <- "STOP"
} else {
  final_decision <- "PASS_WITH_REVIEW"
}
failed_required <- audit_checks[ audit_checks$required_for_schema_pass & audit_checks$status != "PASS", , drop = FALSE ]
timepoint_report_lines <- if (nrow(candidate_timepoint_summary) > 0L) {
  paste0( candidate_timepoint_summary$timepoint, ": ", candidate_timepoint_summary$n_biological_libraries, " libraries; ", candidate_timepoint_summary$n_nuclei, " nuclei" )
} else {
  "No candidate time-course libraries were resolved."
}
if (final_decision == "STOP") {
  for (i in seq_len(nrow(failed_required))) {
  }
} else {
}
if (final_decision == "STOP") {
  stop( "GSE234774 Step 9A2 decision is STOP. Review the final report.", call. = FALSE )
}
