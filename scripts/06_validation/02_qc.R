options( stringsAsFactors = FALSE, warn = 1, future.globals.maxSize = 12 * 1024^3 )
SEED <- 20260803L
FORCE_REIMPORT_FULL_MATRIX <- FALSE
FORCE_REWRITE_CANDIDATE_MATRIX <- FALSE
FORCE_REWRITE_FILTERED_MATRIX <- FALSE
BPCELLS_SORT_BYTES <- 4294967296
BPCELLS_LOAD_BYTES <- 67108864
MIN_FEATURES_ABSOLUTE <- 200L
MIN_COUNTS_ABSOLUTE <- 500L
LOWER_MAD_K <- 5
UPPER_MAD_K <- 6
UPPER_QUANTILE <- 0.995
MITO_MIN_CUTOFF <- 5
MITO_MAX_CUTOFF <- 20
HB_MIN_CUTOFF <- 1
HB_MAX_CUTOFF <- 10
MIN_GLOBAL_QC_RETENTION <- 0.50
MIN_LIBRARY_QC_RETENTION <- 0.30
MIN_QC_CELLS_PER_LIBRARY <- 100L
MIN_LIBRARIES_PER_TIMEPOINT <- 2L
expected_timepoints <- c( "Uninjured", "1 dpi", "4 dpi", "7 dpi", "14 dpi", "1 month", "2 months" )
project_root <- Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
input_dir <- file.path( project_root, "01_raw_data", "GSE234774" )
step9A2_dir <- file.path( project_root, "04_results", "step9A2_GSE234774_targeted_metadata_resolution" )
output_dir <- file.path( project_root, "04_results", "step9B_GSE234774_timecourse_QC_only" )
tables_dir <- file.path(output_dir, "tables")
figures_dir <- file.path(output_dir, "figures")
rds_dir <- file.path(output_dir, "RDS")
bpcells_dir <- file.path(output_dir, "BPCells")
tmp_dir <- file.path(output_dir, "BPCells_tmp")
for (d in c( output_dir, tables_dir, figures_dir, rds_dir, bpcells_dir, tmp_dir )) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}
run_log <- file.path(output_dir, "00_run_log.txt")
if (file.exists(run_log)) {
  file.remove(run_log)
}
mtx_path <- file.path( input_dir, "filtered_scRNA.mtx" )
barcode_path <- file.path( input_dir, "GSE234774_rnaseq_barcodes.txt.gz" )
feature_path <- file.path( input_dir, "GSE234774_rnaseq_features.txt.gz" )
metadata_path <- file.path( input_dir, "GSE234774_rnaseq_meta.txt.gz" )
candidate_manifest_path <- file.path( step9A2_dir, "06_candidate_timecourse_cohort_manifest.csv" )
log_message <- function(...) {
  msg <- paste0(...)
  stamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  line <- paste0("[", stamp, "] ", msg)
  cat(line, "\n")
  cat(line, "\n", file = run_log, append = TRUE)
}
write_csv_safe <- function(x, filename) {
  path <- file.path(tables_dir, filename)
  utils::write.csv(x, path, row.names = FALSE, na = "")
  log_message("Wrote: ", path)
  invisible(path)
}
stop_with_report <- function(reason) {
  report_path <- file.path( output_dir, "10_STEP9B_QC_FINAL_REPORT.txt" )
  report <- c( "GSE234774 STEP 9B QC-ONLY", "============================================================", paste0("Run time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), paste0("Output directory: ", output_dir), "", "FINAL DECISION: STOP", "", paste0("Reason: ", reason) )
  writeLines(report, report_path, useBytes = TRUE)
  stop(reason, call. = FALSE)
}
trim_na <- function(x) {
  x <- trimws(as.character(x))
  x[x == ""] <- NA_character_
  x
}
collapse_unique <- function(x, sep = " | ", max_values = 20L) {
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
safe_median <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) {
    return(NA_real_)
  }
  stats::median(x, na.rm = TRUE)
}
safe_mad <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) {
    return(NA_real_)
  }
  stats::mad( x, center = stats::median(x), constant = 1.4826, na.rm = TRUE )
}
safe_quantile <- function(x, prob) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) {
    return(NA_real_)
  }
  as.numeric( stats::quantile( x, probs = prob, names = FALSE, na.rm = TRUE, type = 8 ) )
}
safe_package_version <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    return(NA_character_)
  }
  as.character(utils::packageVersion(pkg))
}
parse_library_from_barcode <- function(barcode) {
  sub( "-[ACGTN]+(?:-[0-9]+)?$", "", as.character(barcode), perl = TRUE, ignore.case = TRUE )
}
read_gz_lines <- function(path) {
  con <- gzfile(path, open = "rt")
  on.exit(try(close(con), silent = TRUE), add = TRUE)
  readLines(con, warn = FALSE, encoding = "UTF-8")
}
read_features <- function(path) {
  x <- data.table::fread( path, header = FALSE, sep = "\t", data.table = FALSE, showProgress = FALSE )
  if (ncol(x) == 1L) {
    names(x) <- "feature_id"
    x$feature_name <- x$feature_id
    x$feature_type <- "Gene Expression"
  } else if (ncol(x) == 2L) {
    names(x) <- c("feature_id", "feature_name")
    x$feature_type <- "Gene Expression"
  } else {
    names(x)[1:3] <- c( "feature_id", "feature_name", "feature_type" )
  }
  x$feature_id <- as.character(x$feature_id)
  x$feature_name <- as.character(x$feature_name)
  missing_name <- is.na(trim_na(x$feature_name))
  x$feature_name[missing_name] <- x$feature_id[missing_name]
  x$feature_name_unique <- make.unique( x$feature_name, sep = "_dup" )
  x
}
find_exact_column <- function(df, candidates) {
  actual <- names(df)
  actual_norm <- tolower( gsub("[^a-z0-9]+", "_", actual) )
  candidate_norm <- tolower( gsub("[^a-z0-9]+", "_", candidates) )
  idx <- match(candidate_norm, actual_norm)
  idx <- idx[!is.na(idx)]
  if (length(idx) == 0L) {
    return(NA_character_)
  }
  actual[[idx[[1L]]]]
}
open_bpcells_safe <- function(path) {
  tryCatch( BPCells::open_matrix_dir(path), error = function(e) e )
}
log_message("Starting package preflight.")
cran_required <- c( "Matrix", "data.table", "ggplot2", "patchwork", "scales", "Seurat", "SeuratObject" )
missing_cran <- cran_required[ !vapply( cran_required, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1) ) ]
if (length(missing_cran) > 0L) {
  stop_with_report( paste0( "Missing required package(s): ", paste(missing_cran, collapse = ", "), ". Install them before rerunning." ) )
}
if (utils::packageVersion("Seurat") < "5.0.0") {
  stop_with_report( paste0( "Seurat >=5.0.0 is required. Detected version: ", utils::packageVersion("Seurat") ) )
}
if (!requireNamespace("BPCells", quietly = TRUE)) {
  log_message( "BPCells is not installed. Attempting installation." )
  try( utils::install.packages( "BPCells", repos = c( "https://bnprks.r-universe.dev", "https://cloud.r-project.org" ), dependencies = TRUE ), silent = TRUE )
}
if (!requireNamespace("BPCells", quietly = TRUE)) {
  stop_with_report( paste( "BPCells could not be installed automatically.", "Run the following command, restart RStudio, and rerun:", 'install.packages("BPCells", repos=c("https://bnprks.r-universe.dev","https://cloud.r-project.org"))' ) )
}
if (!exists( "import_matrix_market", where = asNamespace("BPCells"), inherits = FALSE )) {
  stop_with_report( "The installed BPCells version lacks import_matrix_market(). Update BPCells." )
}
suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(Seurat)
  library(SeuratObject)
  library(BPCells)
})
package_versions <- data.frame( package = c( "R", "Matrix", "data.table", "ggplot2", "patchwork", "scales", "Seurat", "SeuratObject", "BPCells" ), version = c( paste0(R.version$major, ".", R.version$minor), safe_package_version("Matrix"), safe_package_version("data.table"), safe_package_version("ggplot2"), safe_package_version("patchwork"), safe_package_version("scales"), safe_package_version("Seurat"), safe_package_version("SeuratObject"), safe_package_version("BPCells") ), stringsAsFactors = FALSE )
write_csv_safe( package_versions, "00_package_versions.csv" )
required_paths <- c( mtx_path, barcode_path, feature_path, metadata_path, candidate_manifest_path )
missing_paths <- required_paths[!file.exists(required_paths)]
if (length(missing_paths) > 0L) {
  stop_with_report( paste0( "Missing required input file(s): ", paste(missing_paths, collapse = " ; ") ) )
}
input_inventory <- data.frame( role = c( "Uncompressed Matrix Market counts", "Expression barcodes", "Expression features", "Author cell metadata", "Frozen Step 9A2 cohort manifest" ), path = normalizePath( required_paths, winslash = "/", mustWork = TRUE ), size_bytes = as.numeric( file.info(required_paths)$size ), stringsAsFactors = FALSE )
write_csv_safe( input_inventory, "01_input_inventory.csv" )
log_message( "Using uncompressed MTX: ", mtx_path )
log_message( "BPCells sort memory: ", round(BPCELLS_SORT_BYTES / 1024^3, 2), " GiB" )
log_message( "BPCells contiguous load size: ", round(BPCELLS_LOAD_BYTES / 1024^2, 2), " MiB" )
log_message("Reading frozen Step 9A2 cohort manifest.")
candidate_manifest <- utils::read.csv( candidate_manifest_path, stringsAsFactors = FALSE, check.names = FALSE )
required_manifest_columns <- c( "biological_library", "n_nuclei", "experiment", "condition_label", "replicate", "timepoint_standardised", "candidate_main_timecourse_cohort" )
missing_manifest_columns <- setdiff( required_manifest_columns, names(candidate_manifest) )
if (length(missing_manifest_columns) > 0L) {
  stop_with_report( paste0( "Candidate manifest is missing required column(s): ", paste(missing_manifest_columns, collapse = ", ") ) )
}
candidate_manifest <- candidate_manifest[ as.logical( candidate_manifest$candidate_main_timecourse_cohort ), , drop = FALSE ]
candidate_manifest$biological_library <- trim_na( candidate_manifest$biological_library )
candidate_manifest$timepoint_standardised <- trim_na( candidate_manifest$timepoint_standardised )
candidate_manifest$replicate <- trim_na( candidate_manifest$replicate )
candidate_libraries <- candidate_manifest$biological_library
expected_candidate_cells <- sum(candidate_manifest$n_nuclei)
if (length(candidate_libraries) != 20L) {
  stop_with_report( paste0( "Expected 20 candidate libraries, detected ", length(candidate_libraries), "." ) )
}
if (expected_candidate_cells != 151712L) {
  stop_with_report( paste0( "Expected 151,712 candidate nuclei, detected ", expected_candidate_cells, "." ) )
}
if (!all( expected_timepoints %in% candidate_manifest$timepoint_standardised )) {
  stop_with_report( "The frozen manifest does not contain all expected timepoints." )
}
log_message("Reading expression barcodes.")
all_barcodes <- read_gz_lines(barcode_path)
all_barcodes <- trimws(all_barcodes)
all_barcodes <- all_barcodes[all_barcodes != ""]
log_message("Reading feature table.")
features <- read_features(feature_path)
log_message("Reading complete author metadata.")
cell_metadata_all <- data.table::fread( metadata_path, header = TRUE, sep = "\t", data.table = FALSE, showProgress = TRUE, check.names = FALSE )
names(cell_metadata_all) <- make.unique( names(cell_metadata_all), sep = "_dup" )
barcode_column <- find_exact_column( cell_metadata_all, c("barcode", "cell_barcode") )
if (is.na(barcode_column)) {
  stop_with_report( "Could not identify the barcode column in metadata." )
}
metadata_barcodes <- trim_na( cell_metadata_all[[barcode_column]] )
if (length(all_barcodes) != nrow(cell_metadata_all)) {
  stop_with_report( paste0( "Barcode count and metadata row count differ: ", length(all_barcodes), " versus ", nrow(cell_metadata_all), "." ) )
}
if (!identical(all_barcodes, metadata_barcodes)) {
  metadata_match <- match( all_barcodes, metadata_barcodes )
  if (anyNA(metadata_match)) {
    stop_with_report( "Author metadata cannot be aligned to every expression barcode." )
  }
  cell_metadata_all <- cell_metadata_all[ metadata_match, , drop = FALSE ]
  metadata_barcodes <- trim_na( cell_metadata_all[[barcode_column]] )
}
if (!identical(all_barcodes, metadata_barcodes)) {
  stop_with_report( "Barcode order still differs after metadata alignment." )
}
if (anyDuplicated(all_barcodes) > 0L) {
  stop_with_report( "Duplicated expression barcodes were detected." )
}
parsed_library_all <- parse_library_from_barcode( all_barcodes )
candidate_indices <- which( parsed_library_all %in% candidate_libraries )
candidate_barcodes <- all_barcodes[ candidate_indices ]
candidate_library_from_barcode <- parsed_library_all[ candidate_indices ]
if (length(candidate_indices) != expected_candidate_cells) {
  stop_with_report( paste0( "Frozen manifest expects ", expected_candidate_cells, " nuclei, but barcode parsing selected ", length(candidate_indices), "." ) )
}
library_count_check <- as.data.frame( table( biological_library = candidate_library_from_barcode ), stringsAsFactors = FALSE )
names(library_count_check)[ names(library_count_check) == "Freq" ] <- "n_cells_from_barcode"
library_count_check <- merge( candidate_manifest[ , c( "biological_library", "n_nuclei", "timepoint_standardised", "replicate" ), drop = FALSE ], library_count_check, by = "biological_library", all = TRUE, sort = FALSE )
library_count_check$count_match <- with( library_count_check, n_nuclei == n_cells_from_barcode )
write_csv_safe( library_count_check, "02_frozen_library_count_check.csv" )
if ( anyNA(library_count_check$count_match) || !all(library_count_check$count_match) ) {
  stop_with_report( "At least one frozen library count does not match." )
}
full_matrix_dir <- file.path( bpcells_dir, "GSE234774_full_435099" )
candidate_matrix_dir <- file.path( bpcells_dir, "GSE234774_candidate_timecourse_151712" )
filtered_matrix_dir <- file.path( bpcells_dir, "GSE234774_candidate_timecourse_QC_filtered" )
if ( isTRUE(FORCE_REIMPORT_FULL_MATRIX) && dir.exists(full_matrix_dir) ) {
  log_message( "Deleting the existing full BPCells directory before re-import." )
  unlink( full_matrix_dir, recursive = TRUE, force = TRUE )
}
if (!dir.exists(full_matrix_dir)) {
  log_message( "Starting BPCells import from the uncompressed MTX." )
  log_message( "Import in progress." )
  imported_full <- BPCells::import_matrix_market( mtx_path = mtx_path, outdir = full_matrix_dir, row_names = features$feature_name_unique, col_names = all_barcodes, row_major = FALSE, tmpdir = tmp_dir, load_bytes = BPCELLS_LOAD_BYTES, sort_bytes = BPCELLS_SORT_BYTES )
  rm(imported_full)
  invisible(gc())
  log_message("Completed full BPCells import.")
} else {
  log_message( "Reusing existing full BPCells matrix: ", full_matrix_dir )
}
full_counts <- open_bpcells_safe( full_matrix_dir )
if (inherits(full_counts, "error")) {
  stop_with_report( paste0( "The full BPCells directory cannot be opened and may be incomplete: ", conditionMessage(full_counts), ". Set FORCE_REIMPORT_FULL_MATRIX <- TRUE and rerun." ) )
}
expected_full_dim <- c( nrow(features), length(all_barcodes) )
if (!identical( as.integer(dim(full_counts)), as.integer(expected_full_dim) )) {
  stop_with_report( paste0( "Full BPCells dimensions are ", paste(dim(full_counts), collapse = " x "), "; expected ", paste(expected_full_dim, collapse = " x "), "." ) )
}
if (!identical( colnames(full_counts), all_barcodes )) {
  stop_with_report( "Full BPCells barcode order differs from the barcode file." )
}
if ( isTRUE(FORCE_REWRITE_CANDIDATE_MATRIX) && dir.exists(candidate_matrix_dir) ) {
  log_message( "Deleting the existing candidate BPCells directory." )
  unlink( candidate_matrix_dir, recursive = TRUE, force = TRUE )
}
if (!dir.exists(candidate_matrix_dir)) {
  log_message( "Writing the frozen 151,712-nucleus cohort to BPCells." )
  candidate_view <- full_counts[ , candidate_indices, drop = FALSE ]
  BPCells::write_matrix_dir( mat = candidate_view, dir = candidate_matrix_dir, compress = TRUE, overwrite = FALSE )
  rm(candidate_view)
  invisible(gc())
  log_message("Completed candidate BPCells matrix.")
} else {
  log_message( "Reusing existing candidate BPCells matrix: ", candidate_matrix_dir )
}
candidate_counts <- open_bpcells_safe( candidate_matrix_dir )
if (inherits(candidate_counts, "error")) {
  stop_with_report( paste0( "Candidate BPCells directory cannot be opened: ", conditionMessage(candidate_counts), ". Set FORCE_REWRITE_CANDIDATE_MATRIX <- TRUE and rerun." ) )
}
if (!identical( as.integer(dim(candidate_counts)), as.integer(c( nrow(features), expected_candidate_cells )) )) {
  stop_with_report( "Candidate BPCells dimensions do not match the frozen cohort." )
}
if (!identical( colnames(candidate_counts), candidate_barcodes )) {
  stop_with_report( "Candidate BPCells barcode order differs from the frozen cohort." )
}
rm(full_counts)
invisible(gc())
candidate_metadata <- cell_metadata_all[ candidate_indices, , drop = FALSE ]
rownames(candidate_metadata) <- candidate_barcodes
if ("nCount_RNA" %in% names(candidate_metadata)) {
  candidate_metadata$author_nCount_RNA <- suppressWarnings( as.numeric(candidate_metadata$nCount_RNA) )
  candidate_metadata$nCount_RNA <- NULL
}
if ("nFeature_RNA" %in% names(candidate_metadata)) {
  candidate_metadata$author_nFeature_RNA <- suppressWarnings( as.numeric(candidate_metadata$nFeature_RNA) )
  candidate_metadata$nFeature_RNA <- NULL
}
manifest_match <- match( candidate_library_from_barcode, candidate_manifest$biological_library )
if (anyNA(manifest_match)) {
  stop_with_report( "At least one candidate nucleus does not map to the frozen manifest." )
}
candidate_metadata$biological_library <- candidate_library_from_barcode
candidate_metadata$timepoint <- candidate_manifest$timepoint_standardised[ manifest_match ]
candidate_metadata$replicate_id <- candidate_manifest$replicate[ manifest_match ]
candidate_metadata$experiment_block <- candidate_manifest$experiment[ manifest_match ]
candidate_metadata$condition_label <- candidate_manifest$condition_label[ manifest_match ]
candidate_metadata$validation_cohort <- "GSE234774_frozen_timecourse"
candidate_metadata$age_status <- "METHODS_REVIEW_REQUIRED"
candidate_metadata$statistical_unit <- "biological_library"
candidate_metadata$timepoint <- factor( candidate_metadata$timepoint, levels = expected_timepoints, ordered = TRUE )
if ( anyNA(candidate_metadata$biological_library) || anyNA(candidate_metadata$timepoint) || anyNA(candidate_metadata$replicate_id) ) {
  stop_with_report( "Candidate metadata contains unresolved library, timepoint, or replicate." )
}
rm( cell_metadata_all, metadata_barcodes, parsed_library_all )
invisible(gc())
log_message( "Creating the BPCells-backed pre-QC Seurat object." )
set.seed(SEED)
obj_preQC <- Seurat::CreateSeuratObject( counts = candidate_counts, assay = "RNA", project = "GSE234774_timecourse_preQC", meta.data = candidate_metadata, min.cells = 0, min.features = 0 )
Seurat::DefaultAssay(obj_preQC) <- "RNA"
if ( nrow(obj_preQC) != nrow(features) || ncol(obj_preQC) != expected_candidate_cells ) {
  stop_with_report( "Pre-QC Seurat object dimensions are incorrect." )
}
metadata_alignment_check <- data.frame( check = c( "Metadata rows", "Matrix columns", "Barcode order identical", "Biological libraries", "Timepoints" ), observed = c( nrow(obj_preQC@meta.data), ncol(obj_preQC), identical( rownames(obj_preQC@meta.data), colnames(obj_preQC) ), length(unique( obj_preQC$biological_library )), length(unique( obj_preQC$timepoint )) ), expected = c( expected_candidate_cells, expected_candidate_cells, TRUE, 20L, 7L ), stringsAsFactors = FALSE )
metadata_alignment_check$pass <- metadata_alignment_check$observed == metadata_alignment_check$expected
write_csv_safe( metadata_alignment_check, "03_candidate_metadata_alignment_check.csv" )
if (!all(metadata_alignment_check$pass)) {
  stop_with_report( "Candidate metadata alignment check failed." )
}
qc_metric_consistency <- data.frame( metric = c( "nCount_RNA", "nFeature_RNA" ), author_column_available = c( "author_nCount_RNA" %in% names(obj_preQC@meta.data), "author_nFeature_RNA" %in% names(obj_preQC@meta.data) ), exact_match_fraction = c( if ( "author_nCount_RNA" %in% names(obj_preQC@meta.data) ) {
      mean( obj_preQC$nCount_RNA == obj_preQC$author_nCount_RNA, na.rm = TRUE )
    } else {
      NA_real_
    },
    if ( "author_nFeature_RNA" %in% names(obj_preQC@meta.data) ) {
      mean( obj_preQC$nFeature_RNA == obj_preQC$author_nFeature_RNA, na.rm = TRUE )
    } else {
      NA_real_
    }
  ), maximum_absolute_difference = c(
    if ( "author_nCount_RNA" %in% names(obj_preQC@meta.data) ) {
      max( abs( obj_preQC$nCount_RNA - obj_preQC$author_nCount_RNA ), na.rm = TRUE )
    } else {
      NA_real_
    },
    if ( "author_nFeature_RNA" %in% names(obj_preQC@meta.data) ) {
      max( abs( obj_preQC$nFeature_RNA - obj_preQC$author_nFeature_RNA ), na.rm = TRUE )
    } else {
      NA_real_
    }
  ), stringsAsFactors = FALSE
)
write_csv_safe( qc_metric_consistency, "04_author_vs_recalculated_QC_metrics.csv" )
log_message( "Calculating percent.mt and percent.HB." )
gene_names <- rownames(obj_preQC)
mt_genes <- gene_names[ grepl( "^mt-", gene_names, ignore.case = TRUE ) ]
hb_genes <- gene_names[ grepl( "^hb[ab]", gene_names, ignore.case = TRUE ) ]
if (length(mt_genes) > 0L) {
  obj_preQC[["percent.mt"]] <- Seurat::PercentageFeatureSet( obj_preQC, assay = "RNA", features = mt_genes )
} else {
  obj_preQC$percent.mt <- 0
}
if (length(hb_genes) > 0L) {
  obj_preQC[["percent.HB"]] <- Seurat::PercentageFeatureSet( obj_preQC, assay = "RNA", features = hb_genes )
} else {
  obj_preQC$percent.HB <- 0
}
obj_preQC$log10GenesPerUMI <- with( obj_preQC@meta.data, log10(pmax(nFeature_RNA, 1)) / log10(pmax(nCount_RNA, 2)) )
feature_set_audit <- data.frame( feature_set = c( "Mitochondrial", "Hemoglobin" ), n_genes = c( length(mt_genes), length(hb_genes) ), examples = c( collapse_unique(head(mt_genes, 20L)), collapse_unique(head(hb_genes, 20L)) ), stringsAsFactors = FALSE )
write_csv_safe( feature_set_audit, "04B_QC_feature_set_audit.csv" )
qc_meta <- obj_preQC@meta.data
qc_meta$barcode <- rownames(qc_meta)
preQC_summary_by_library <- data.table::as.data.table( qc_meta )[ , .( n_cells = .N, median_nFeature_RNA = median(nFeature_RNA), q01_nFeature_RNA = quantile(nFeature_RNA, 0.01), q99_nFeature_RNA = quantile(nFeature_RNA, 0.99), median_nCount_RNA = median(nCount_RNA), q01_nCount_RNA = quantile(nCount_RNA, 0.01), q99_nCount_RNA = quantile(nCount_RNA, 0.99), median_percent_mt = median(percent.mt), q99_percent_mt = quantile(percent.mt, 0.99), median_percent_HB = median(percent.HB), q99_percent_HB = quantile(percent.HB, 0.99), median_log10GenesPerUMI = median(log10GenesPerUMI) ), by = .( biological_library, timepoint, replicate_id ) ]
write_csv_safe( as.data.frame(preQC_summary_by_library), "05_preQC_summary_by_library.csv" )
derive_thresholds <- function(df) {
  feature_med <- safe_median( df$nFeature_RNA )
  feature_mad <- safe_mad( df$nFeature_RNA )
  count_med <- safe_median( df$nCount_RNA )
  count_mad <- safe_mad( df$nCount_RNA )
  min_features <- max( MIN_FEATURES_ABSOLUTE, floor( min( safe_quantile( df$nFeature_RNA, 0.01 ), feature_med - LOWER_MAD_K * feature_mad, na.rm = TRUE ) ) )
  max_features <- ceiling( max( safe_quantile( df$nFeature_RNA, UPPER_QUANTILE ), feature_med + UPPER_MAD_K * feature_mad, na.rm = TRUE ) )
  min_counts <- max( MIN_COUNTS_ABSOLUTE, floor( min( safe_quantile( df$nCount_RNA, 0.01 ), count_med - LOWER_MAD_K * count_mad, na.rm = TRUE ) ) )
  max_counts <- ceiling( max( safe_quantile( df$nCount_RNA, UPPER_QUANTILE ), count_med + UPPER_MAD_K * count_mad, na.rm = TRUE ) )
  mt_med <- safe_median( df$percent.mt )
  mt_mad <- safe_mad( df$percent.mt )
  max_percent_mt <- min( MITO_MAX_CUTOFF, max( MITO_MIN_CUTOFF, safe_quantile( df$percent.mt, UPPER_QUANTILE ), mt_med + UPPER_MAD_K * mt_mad, na.rm = TRUE ) )
  hb_med <- safe_median( df$percent.HB )
  hb_mad <- safe_mad( df$percent.HB )
  max_percent_HB <- min( HB_MAX_CUTOFF, max( HB_MIN_CUTOFF, safe_quantile( df$percent.HB, UPPER_QUANTILE ), hb_med + UPPER_MAD_K * hb_mad, na.rm = TRUE ) )
  data.frame( n_cells_before = nrow(df), min_nFeature_RNA = as.integer( min_features ), max_nFeature_RNA = as.integer( max_features ), min_nCount_RNA = as.integer( min_counts ), max_nCount_RNA = as.integer( max_counts ), max_percent_mt = max_percent_mt, max_percent_HB = max_percent_HB, median_nFeature_RNA = feature_med, MAD_nFeature_RNA = feature_mad, median_nCount_RNA = count_med, MAD_nCount_RNA = count_mad, median_percent_mt = mt_med, MAD_percent_mt = mt_mad, median_percent_HB = hb_med, MAD_percent_HB = hb_mad, stringsAsFactors = FALSE )
}
library_split <- split( qc_meta, qc_meta$biological_library )
threshold_rows <- lapply( names(library_split), function(lib) {
    out <- derive_thresholds( library_split[[lib]] )
    out$biological_library <- lib
    out
  }
)
qc_thresholds <- do.call( rbind, threshold_rows )
qc_thresholds <- merge( candidate_manifest[ , c( "biological_library", "timepoint_standardised", "replicate" ), drop = FALSE ], qc_thresholds, by = "biological_library", all.x = TRUE, sort = FALSE )
qc_thresholds <- qc_thresholds[ match( candidate_libraries, qc_thresholds$biological_library ), , drop = FALSE ]
write_csv_safe( qc_thresholds, "06_QC_thresholds_by_library.csv" )
threshold_match <- match( qc_meta$biological_library, qc_thresholds$biological_library )
if (anyNA(threshold_match)) {
  stop_with_report( "At least one nucleus lacks its library-specific QC thresholds." )
}
qc_meta$qc_fail_low_features <- qc_meta$nFeature_RNA <
  qc_thresholds$min_nFeature_RNA[ threshold_match ]
qc_meta$qc_fail_high_features <- qc_meta$nFeature_RNA >
  qc_thresholds$max_nFeature_RNA[ threshold_match ]
qc_meta$qc_fail_low_counts <- qc_meta$nCount_RNA <
  qc_thresholds$min_nCount_RNA[ threshold_match ]
qc_meta$qc_fail_high_counts <- qc_meta$nCount_RNA >
  qc_thresholds$max_nCount_RNA[ threshold_match ]
qc_meta$qc_fail_high_mt <- qc_meta$percent.mt >
  qc_thresholds$max_percent_mt[ threshold_match ]
qc_meta$qc_fail_high_HB <- qc_meta$percent.HB >
  qc_thresholds$max_percent_HB[ threshold_match ]
failure_columns <- c( "qc_fail_low_features", "qc_fail_high_features", "qc_fail_low_counts", "qc_fail_high_counts", "qc_fail_high_mt", "qc_fail_high_HB" )
qc_meta$qc_pass <- rowSums( qc_meta[ , failure_columns, drop = FALSE ] ) == 0L
qc_meta$qc_failure_reason <- apply( qc_meta[ , failure_columns, drop = FALSE ], 1L, function(z) {
    failed <- failure_columns[as.logical(z)]
    if (length(failed) == 0L) {
      return("PASS")
    }
    paste( sub("^qc_fail_", "", failed), collapse = ";" )
  }
)
obj_preQC <- Seurat::AddMetaData( obj_preQC, metadata = qc_meta[ , c( failure_columns, "qc_pass", "qc_failure_reason" ), drop = FALSE ] )
qc_retention_by_library <- data.table::as.data.table( qc_meta )[ , .( n_before = .N, n_QC_pass = sum(qc_pass), n_QC_fail = sum(!qc_pass), QC_retention = mean(qc_pass), fail_low_features = sum(qc_fail_low_features), fail_high_features = sum(qc_fail_high_features), fail_low_counts = sum(qc_fail_low_counts), fail_high_counts = sum(qc_fail_high_counts), fail_high_mt = sum(qc_fail_high_mt), fail_high_HB = sum(qc_fail_high_HB) ), by = .( biological_library, timepoint, replicate_id ) ]
write_csv_safe( as.data.frame(qc_retention_by_library), "07_QC_retention_by_library.csv" )
qc_retention_by_timepoint <- data.table::as.data.table( qc_meta )[ , .( n_before = .N, n_QC_pass = sum(qc_pass), n_QC_fail = sum(!qc_pass), QC_retention = mean(qc_pass), n_libraries = uniqueN( biological_library ) ), by = timepoint ]
write_csv_safe( as.data.frame(qc_retention_by_timepoint), "07B_QC_retention_by_timepoint.csv" )
failure_reason_summary <- data.table::as.data.table( qc_meta )[ , .( n_cells = .N ), by = .( biological_library, timepoint, qc_failure_reason ) ]
write_csv_safe( as.data.frame(failure_reason_summary), "07C_QC_failure_reasons_by_library.csv" )
cell_qc_path <- file.path( tables_dir, "07D_cell_level_QC_status.tsv.gz" )
data.table::fwrite( qc_meta, file = cell_qc_path, sep = "\t", quote = FALSE, na = "NA", compress = "gzip" )
log_message("Wrote: ", cell_qc_path)
set.seed(SEED)
plot_indices <- unlist( lapply( split( seq_len(nrow(qc_meta)), qc_meta$timepoint ), function(idx) {
      sample( idx, size = min(length(idx), 10000L), replace = FALSE )
    }
  ), use.names = FALSE
)
plot_meta <- qc_meta[ plot_indices, , drop = FALSE ]
p_feature <- ggplot( plot_meta, aes( x = timepoint, y = nFeature_RNA ) ) + geom_violin( scale = "width", trim = TRUE ) + scale_y_log10() + labs( x = NULL, y = "nFeature_RNA (log10)", title = "Features before QC" ) + theme_bw() + theme( axis.text.x = element_text( angle = 45, hjust = 1 ) )
p_count <- ggplot( plot_meta, aes( x = timepoint, y = nCount_RNA ) ) + geom_violin( scale = "width", trim = TRUE ) + scale_y_log10() + labs( x = NULL, y = "nCount_RNA (log10)", title = "Counts before QC" ) + theme_bw() + theme( axis.text.x = element_text( angle = 45, hjust = 1 ) )
p_mt <- ggplot( plot_meta, aes( x = timepoint, y = percent.mt ) ) + geom_violin( scale = "width", trim = TRUE ) + labs( x = NULL, y = "percent.mt", title = "Mitochondrial fraction before QC" ) + theme_bw() + theme( axis.text.x = element_text( angle = 45, hjust = 1 ) )
p_hb <- ggplot( plot_meta, aes( x = timepoint, y = percent.HB ) ) + geom_violin( scale = "width", trim = TRUE ) + labs( x = NULL, y = "percent.HB", title = "Hemoglobin fraction before QC" ) + theme_bw() + theme( axis.text.x = element_text( angle = 45, hjust = 1 ) )
qc_violin <- ( p_feature + p_count ) / ( p_mt + p_hb )
ggsave( filename = file.path( figures_dir, "01_QC_distributions_before_filter.pdf" ), plot = qc_violin, width = 14, height = 10, limitsize = FALSE )
p_retention <- ggplot( as.data.frame(qc_retention_by_library), aes( x = reorder( biological_library, QC_retention ), y = QC_retention ) ) + geom_col() + coord_flip() + scale_y_continuous( limits = c(0, 1), labels = scales::percent ) + labs( x = NULL, y = "QC retention", title = "Library-level QC retention" ) + theme_bw()
ggsave( filename = file.path( figures_dir, "02_QC_retention_by_library.pdf" ), plot = p_retention, width = 10, height = 8, limitsize = FALSE )
qc_pass_cells <- qc_meta$barcode[ qc_meta$qc_pass ]
qc_pass_index <- match( qc_pass_cells, colnames(candidate_counts) )
if (anyNA(qc_pass_index)) {
  stop_with_report( "At least one QC-pass barcode is absent from the candidate matrix." )
}
if ( isTRUE(FORCE_REWRITE_FILTERED_MATRIX) && dir.exists(filtered_matrix_dir) ) {
  log_message( "Deleting the existing QC-filtered BPCells directory." )
  unlink( filtered_matrix_dir, recursive = TRUE, force = TRUE )
}
if (!dir.exists(filtered_matrix_dir)) {
  log_message( "Writing QC-pass nuclei to the filtered BPCells matrix." )
  filtered_view <- candidate_counts[ , qc_pass_index, drop = FALSE ]
  BPCells::write_matrix_dir( mat = filtered_view, dir = filtered_matrix_dir, compress = TRUE, overwrite = FALSE )
  rm(filtered_view)
  invisible(gc())
  log_message( "Completed QC-filtered BPCells matrix." )
} else {
  log_message( "Reusing existing QC-filtered BPCells matrix: ", filtered_matrix_dir )
}
filtered_counts <- open_bpcells_safe( filtered_matrix_dir )
if (inherits(filtered_counts, "error")) {
  stop_with_report( paste0( "QC-filtered BPCells directory cannot be opened: ", conditionMessage(filtered_counts), ". Set FORCE_REWRITE_FILTERED_MATRIX <- TRUE and rerun." ) )
}
if ( ncol(filtered_counts) != length(qc_pass_cells) || !identical( colnames(filtered_counts), qc_pass_cells ) ) {
  stop_with_report( paste( "The existing filtered matrix does not match the current QC-pass cells.", "Set FORCE_REWRITE_FILTERED_MATRIX <- TRUE and rerun." ) )
}
filtered_metadata <- qc_meta[ match( qc_pass_cells, qc_meta$barcode ), , drop = FALSE ]
rownames(filtered_metadata) <- qc_pass_cells
log_message( "Creating the QC-filtered Seurat object." )
obj_QC <- Seurat::CreateSeuratObject( counts = filtered_counts, assay = "RNA", project = "GSE234774_timecourse_QC_filtered", meta.data = filtered_metadata, min.cells = 0, min.features = 0 )
Seurat::DefaultAssay(obj_QC) <- "RNA"
obj_QC@misc$GSE234774_step9B_QC <- list( input_mtx = normalizePath( mtx_path, winslash = "/", mustWork = TRUE ), biological_unit = "biological_library", cohort_age_status = "METHODS_REVIEW_REQUIRED", timepoints = expected_timepoints, seed = SEED, bpcells_sort_bytes = BPCELLS_SORT_BYTES, bpcells_load_bytes = BPCELLS_LOAD_BYTES, completed_at = as.character(Sys.time()), note = paste( "QC only.", "Doublet detection, normalization, integration, clustering,", "annotation, signature transfer, and spatial analysis were not run." ) )
preQC_object_path <- file.path( rds_dir, "GSE234774_step9B_timecourse_preQC_BPCells.rds" )
qc_object_path <- file.path( rds_dir, "GSE234774_step9B_timecourse_QC_filtered_BPCells.rds" )
saveRDS( obj_preQC, preQC_object_path, compress = FALSE )
saveRDS( obj_QC, qc_object_path, compress = FALSE )
log_message( "Saved pre-QC object: ", preQC_object_path )
log_message( "Saved QC-filtered object: ", qc_object_path )
results_manifest <- list( package_versions = package_versions, input_inventory = input_inventory, candidate_manifest = candidate_manifest, feature_set_audit = feature_set_audit, qc_thresholds = qc_thresholds, qc_retention_by_library = as.data.frame( qc_retention_by_library ), qc_retention_by_timepoint = as.data.frame( qc_retention_by_timepoint ), paths = list( preQC_object = preQC_object_path, QC_filtered_object = qc_object_path, full_BPCell_matrix = full_matrix_dir, candidate_BPCell_matrix = candidate_matrix_dir, QC_filtered_BPCell_matrix = filtered_matrix_dir ) )
saveRDS( results_manifest, file.path( rds_dir, "GSE234774_step9B_QC_results_manifest.rds" ), compress = TRUE )
global_QC_retention <- mean( qc_meta$qc_pass )
minimum_library_QC_retention <- min( qc_retention_by_library$QC_retention )
minimum_QC_cells_per_library <- min( qc_retention_by_library$n_QC_pass )
timepoint_library_check <- data.table::as.data.table( qc_retention_by_library )[ n_QC_pass >= MIN_QC_CELLS_PER_LIBRARY, .( n_qualified_libraries = uniqueN( biological_library ) ), by = timepoint ]
timepoint_library_check <- merge( data.frame( timepoint = factor( expected_timepoints, levels = expected_timepoints, ordered = TRUE ) ), as.data.frame(timepoint_library_check), by = "timepoint", all.x = TRUE, sort = FALSE )
timepoint_library_check$n_qualified_libraries[ is.na( timepoint_library_check$n_qualified_libraries ) ] <- 0L
write_csv_safe( timepoint_library_check, "08_timepoint_library_retention_check.csv" )
reload_test <- tryCatch( readRDS(qc_object_path), error = function(e) e )
reload_pass <- !inherits( reload_test, "error" ) && nrow(reload_test) == nrow(features) && ncol(reload_test) == length(qc_pass_cells)
rm(reload_test)
invisible(gc())
audit_checks <- data.frame( check = c( "Uncompressed MTX was used", "Full BPCells matrix dimensions match the audited input", "Frozen cohort contains 20 libraries and 151,712 nuclei", "Candidate metadata and matrix barcodes are identical", "Global QC retention is at least 50%", "Every library QC retention is at least 30%", "Every library retains at least 100 QC-pass nuclei", "Every expected timepoint retains at least two qualified libraries", "QC-filtered BPCells matrix matches QC-pass barcodes", "QC-filtered Seurat object saves and reloads", "Doublet detection was not run", "Normalization/integration/clustering/annotation were not run", "Cohort age has been confirmed from study methods" ), required_for_pass = c( TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE ), status = c( ifelse( identical( basename(mtx_path), "filtered_scRNA.mtx" ), "PASS", "STOP" ), ifelse( identical( as.integer(dim(candidate_counts)), as.integer(c( nrow(features), expected_candidate_cells )) ), "PASS", "STOP" ), ifelse( length(candidate_libraries) == 20L && expected_candidate_cells == 151712L, "PASS", "STOP" ), ifelse( identical( rownames(candidate_metadata), colnames(candidate_counts) ), "PASS", "STOP" ), ifelse( global_QC_retention >= MIN_GLOBAL_QC_RETENTION, "PASS", "STOP" ), ifelse( minimum_library_QC_retention >= MIN_LIBRARY_QC_RETENTION, "PASS", "STOP" ), ifelse( minimum_QC_cells_per_library >= MIN_QC_CELLS_PER_LIBRARY, "PASS", "STOP" ), ifelse( all( timepoint_library_check$n_qualified_libraries >= MIN_LIBRARIES_PER_TIMEPOINT ), "PASS", "STOP" ), ifelse( ncol(filtered_counts) == length(qc_pass_cells) && identical( colnames(filtered_counts), qc_pass_cells ), "PASS", "STOP" ), ifelse( reload_pass, "PASS", "STOP" ), "CONFIRMED_NOT_RUN", "CONFIRMED_NOT_RUN", "REVIEW" ), detail = c( mtx_path, paste0( nrow(features), " features x ", length(all_barcodes), " total barcodes" ), paste0( length(candidate_libraries), " libraries; ", expected_candidate_cells, " candidate nuclei" ), paste0( "metadata rows=", nrow(candidate_metadata), "; matrix columns=", ncol(candidate_counts) ), paste0( "global QC retention=", round( 100 * global_QC_retention, 2 ), "%" ), paste0( "minimum library QC retention=", round( 100 * minimum_library_QC_retention, 2 ), "%" ), paste0( "minimum QC-pass nuclei per library=", minimum_QC_cells_per_library ), paste( paste0( timepoint_library_check$timepoint, ":", timepoint_library_check$n_qualified_libraries ), collapse = ";" ), paste0( nrow(filtered_counts), " features x ", ncol(filtered_counts), " QC-pass nuclei" ), ifelse( reload_pass, "RDS reload PASS", "RDS reload failed or dimensions differ" ), "Deferred to the next independent step.", "Deferred to later steps.", "Age remains METHODS_REVIEW_REQUIRED." ), stringsAsFactors = FALSE )
write_csv_safe( audit_checks, "09_STEP9B_QC_audit_checks.csv" )
required_fail <- audit_checks$required_for_pass & audit_checks$status != "PASS"
final_decision <- if (any(required_fail)) {
  "STOP"
} else {
  "PASS_WITH_METHODS_REVIEW"
}
library_lines <- paste0( qc_retention_by_library$biological_library, " [", qc_retention_by_library$timepoint, "]: ", qc_retention_by_library$n_before, " -> ", qc_retention_by_library$n_QC_pass, " QC-pass nuclei; retention=", round( 100 * qc_retention_by_library$QC_retention, 2 ), "%" )
timepoint_lines <- paste0( qc_retention_by_timepoint$timepoint, ": ", qc_retention_by_timepoint$n_before, " -> ", qc_retention_by_timepoint$n_QC_pass, " QC-pass nuclei; retention=", round( 100 * qc_retention_by_timepoint$QC_retention, 2 ), "%" )
report_lines <- c( "GSE234774 STEP 9B QC-ONLY", "============================================================", paste0( "Run time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S") ), paste0("Project root: ", project_root), paste0("Output directory: ", output_dir), "", paste0("FINAL DECISION: ", final_decision), "", "SCOPE", "The frozen 20-library long-term cohort was reconstructed and QC-filtered.", "The uncompressed filtered_scRNA.mtx file was used.", "Raw counts remain stored in BPCells format.", "No doublet detection, normalization, integration, clustering, annotation, signature transfer, pseudotime, or spatial analysis was run.", "", "INPUT", paste0( "Full matrix: ", nrow(features), " features x ", length(all_barcodes), " nuclei" ), paste0( "Frozen candidate cohort: ", nrow(candidate_counts), " features x ", ncol(candidate_counts), " nuclei" ), paste0( "Biological libraries: ", length(candidate_libraries) ), paste0( "BPCells sort memory: ", round(BPCELLS_SORT_BYTES / 1024^3, 2), " GiB" ), "", "QC RESULT", paste0( "QC-pass object: ", nrow(obj_QC), " features x ", ncol(obj_QC), " nuclei" ), paste0( "Global QC retention: ", round( 100 * global_QC_retention, 2 ), "%" ), paste0( "Minimum library QC retention: ", round( 100 * minimum_library_QC_retention, 2 ), "%" ), paste0( "Minimum QC-pass nuclei in one library: ", minimum_QC_cells_per_library ), "", "TIMEPOINT RETENTION", timepoint_lines, "", "LIBRARY RETENTION", library_lines, "", "BOUNDARIES", "1. biological_library is the statistical replicate.", "2. Individual nuclei are not independent biological replicates.", "3. orig.ident is an experimental block, not a sample identifier.", "4. Cohort age remains METHODS_REVIEW_REQUIRED.", "5. GSE234774 validates frozen GSE162610 programs and does not redefine Arm_1/Arm_2.", "", "PRIMARY QC-FILTERED OBJECT", qc_object_path )
if (final_decision == "STOP") {
  failed <- audit_checks[ audit_checks$required_for_pass & audit_checks$status != "PASS", , drop = FALSE ]
  report_lines <- c( report_lines, "", "STOP REASONS" )
  for (i in seq_len(nrow(failed))) {
    report_lines <- c( report_lines, paste0( i, ". ", failed$check[[i]], " -- ", failed$detail[[i]] ) )
  }
  report_lines <- c( report_lines, "", "Do not start doublet detection or downstream analysis." )
} else {
  report_lines <- c( report_lines, "", "STEP 9B QC DECISION", "Object reconstruction and library-aware QC passed.", "Stop here. Review the QC outputs before running doublet detection." )
}
report_path <- file.path( output_dir, "10_STEP9B_QC_FINAL_REPORT.txt" )
writeLines( report_lines, report_path, useBytes = TRUE )
log_message("Final decision: ", final_decision)
log_message("Final report: ", report_path)
cat("\n")
cat("============================================================\n")
cat("GSE234774 STEP 9B QC-ONLY COMPLETED\n")
cat("FINAL DECISION:", final_decision, "\n")
cat("Output:", output_dir, "\n")
cat("Primary QC-filtered object:\n")
cat(" ", qc_object_path, "\n", sep = "")
cat("Key outputs:\n")
cat(" 1) tables/05_preQC_summary_by_library.csv\n")
cat(" 2) tables/06_QC_thresholds_by_library.csv\n")
cat(" 3) tables/07_QC_retention_by_library.csv\n")
cat(" 4) tables/09_STEP9B_QC_audit_checks.csv\n")
cat(" 5) 10_STEP9B_QC_FINAL_REPORT.txt\n")
cat("============================================================\n")
if (final_decision == "STOP") {
  stop( "GSE234774 Step 9B QC decision is STOP. Review the final report.", call. = FALSE )
}
