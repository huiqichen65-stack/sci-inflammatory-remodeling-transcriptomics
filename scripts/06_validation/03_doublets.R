options( stringsAsFactors = FALSE, warn = 1, future.globals.maxSize = 12 * 1024^3 )
SEED <- 20260803L
DBR_PER_1K <- 0.008
USE_CLUSTERS <- FALSE
MIN_CELLS_PER_GENE <- 3L
FORCE_RERUN_ALL_LIBRARIES <- FALSE
FORCE_REWRITE_SINGLET_MATRIX <- FALSE
MAX_ACCEPTABLE_DOUBLET_RATE <- 0.25
MIN_SINGLET_CELLS_PER_LIBRARY <- 100L
MIN_QUALIFIED_LIBRARIES_PER_TIMEPOINT <- 2L
expected_timepoints <- c( "Uninjured", "1 dpi", "4 dpi", "7 dpi", "14 dpi", "1 month", "2 months" )
project_root <- Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
input_object_path <- file.path( project_root, "04_results", "step9B_GSE234774_timecourse_QC_only", "RDS", "GSE234774_step9B_timecourse_QC_filtered_BPCells.rds" )
output_dir <- file.path( project_root, "04_results", "step9C_GSE234774_samplewise_scDblFinder" )
tables_dir <- file.path(output_dir, "tables")
figures_dir <- file.path(output_dir, "figures")
rds_dir <- file.path(output_dir, "RDS")
checkpoint_dir <- file.path(output_dir, "checkpoints_by_library")
bpcells_dir <- file.path(output_dir, "BPCells")
for (d in c( output_dir, tables_dir, figures_dir, rds_dir, checkpoint_dir, bpcells_dir )) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}
run_log <- file.path(output_dir, "00_run_log.txt")
if (file.exists(run_log)) {
  file.remove(run_log)
}
annotated_object_path <- file.path( rds_dir, "GSE234774_step9C_QC_doublet_annotated_BPCells.rds" )
singlet_object_path <- file.path( rds_dir, "GSE234774_step9C_QC_singlets_BPCells.rds" )
singlet_matrix_dir <- file.path( bpcells_dir, "GSE234774_step9C_QC_singlets" )
write_csv_safe <- function(x, filename) {
  out <- file.path(tables_dir, filename)
  utils::write.csv(x, out, row.names = FALSE, na = "")
  invisible(out)
}
stop_with_report <- function(reason) {
  report_path <- file.path( output_dir, "10_STEP9C_FINAL_REPORT.txt" )
  stop(reason, call. = FALSE)
}
trim_na <- function(x) {
  x <- trimws(as.character(x))
  x[x == ""] <- NA_character_
  x
}
safe_package_version <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    return(NA_character_)
  }
  as.character(utils::packageVersion(pkg))
}
open_bpcells_safe <- function(path) {
  tryCatch( BPCells::open_matrix_dir(path), error = function(e) e )
}
sanitize_filename <- function(x) {
  gsub("[^A-Za-z0-9_.-]+", "_", as.character(x))
}
expected_doublet_rate <- function(n_cells) {
  DBR_PER_1K * (n_cells / 1000)
}
checkpoint_is_usable <- function(path, expected_cells, expected_library) {
  if (!file.exists(path)) {
    return(FALSE)
  }
  x <- tryCatch( data.table::fread( path, data.table = FALSE, showProgress = FALSE ), error = function(e) NULL )
  if (is.null(x)) {
    return(FALSE)
  }
  required <- c( "barcode", "biological_library", "scDblFinder.score", "scDblFinder.class" )
  if (!all(required %in% names(x))) {
    return(FALSE)
  }
  if (nrow(x) != length(expected_cells)) {
    return(FALSE)
  }
  if (!setequal(x$barcode, expected_cells)) {
    return(FALSE)
  }
  if ( length(unique(x$biological_library)) != 1L || unique(x$biological_library) != expected_library ) {
    return(FALSE)
  }
  if ( anyNA(x$scDblFinder.score) || anyNA(x$scDblFinder.class) ) {
    return(FALSE)
  }
  TRUE
}
cran_required <- c( "Matrix", "data.table", "ggplot2", "scales", "Seurat", "SeuratObject", "BPCells" )
missing_cran <- cran_required[ !vapply( cran_required, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1) ) ]
if (length(missing_cran) > 0L) {
  stop_with_report( paste0( "Missing required package(s): ", paste(missing_cran, collapse = ", "), ". Install them before rerunning." ) )
}
bioc_required <- c( "BiocParallel", "SingleCellExperiment", "scDblFinder" )
missing_bioc <- bioc_required[ !vapply( bioc_required, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1) ) ]
if (length(missing_bioc) > 0L) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    utils::install.packages( "BiocManager", repos = "https://cloud.r-project.org" )
  }
  try( BiocManager::install( missing_bioc, ask = FALSE, update = FALSE ), silent = TRUE )
}
missing_bioc_after <- bioc_required[ !vapply( bioc_required, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1) ) ]
if (length(missing_bioc_after) > 0L) {
  stop_with_report( paste0( "Missing Bioconductor package(s): ", paste(missing_bioc_after, collapse = ", "), ". Install them and rerun." ) )
}
suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(ggplot2)
  library(scales)
  library(Seurat)
  library(SeuratObject)
  library(BPCells)
  library(BiocParallel)
  library(SingleCellExperiment)
  library(scDblFinder)
})
package_versions <- data.frame( package = c( "R", "Matrix", "data.table", "ggplot2", "scales", "Seurat", "SeuratObject", "BPCells", "BiocParallel", "SingleCellExperiment", "scDblFinder" ), version = c( paste0(R.version$major, ".", R.version$minor), safe_package_version("Matrix"), safe_package_version("data.table"), safe_package_version("ggplot2"), safe_package_version("scales"), safe_package_version("Seurat"), safe_package_version("SeuratObject"), safe_package_version("BPCells"), safe_package_version("BiocParallel"), safe_package_version("SingleCellExperiment"), safe_package_version("scDblFinder") ), stringsAsFactors = FALSE )
write_csv_safe( package_versions, "00_package_versions.csv" )
if (!file.exists(input_object_path)) {
  stop_with_report( paste0( "QC-filtered object not found: ", input_object_path ) )
}
obj_QC <- readRDS(input_object_path)
if (!inherits(obj_QC, "Seurat")) {
  stop_with_report( "The Step 9B input RDS is not a Seurat object." )
}
if (!"RNA" %in% names(obj_QC@assays)) {
  stop_with_report( "RNA assay is absent from the Step 9B object." )
}
Seurat::DefaultAssay(obj_QC) <- "RNA"
required_metadata <- c( "biological_library", "timepoint", "replicate_id", "nCount_RNA", "nFeature_RNA", "percent.mt", "percent.HB", "qc_pass" )
missing_metadata <- setdiff( required_metadata, names(obj_QC@meta.data) )
if (length(missing_metadata) > 0L) {
  stop_with_report( paste0( "Step 9B object is missing required metadata: ", paste(missing_metadata, collapse = ", ") ) )
}
if (!all(as.logical(obj_QC$qc_pass))) {
  stop_with_report( "The Step 9B input object contains nuclei that did not pass QC." )
}
if (ncol(obj_QC) != 148567L) {
  stop_with_report( paste0( "Expected 148,567 QC-pass nuclei, detected ", ncol(obj_QC), "." ) )
}
libraries <- unique( as.character(obj_QC$biological_library) )
if (length(libraries) != 20L) {
  stop_with_report( paste0( "Expected 20 biological libraries, detected ", length(libraries), "." ) )
}
observed_timepoints <- unique( as.character(obj_QC$timepoint) )
if (!all(expected_timepoints %in% observed_timepoints)) {
  stop_with_report( "The QC object does not contain all seven expected timepoints." )
}
layer_names <- SeuratObject::Layers( obj_QC, assay = "RNA", search = NA )
counts_layers <- layer_names[ grepl("^counts($|\\.)", layer_names) ]
if (length(counts_layers) != 1L) {
  stop_with_report( paste0( "Expected exactly one counts layer, detected: ", paste(counts_layers, collapse = ", ") ) )
}
counts_layer <- counts_layers[[1L]]
counts_all <- SeuratObject::LayerData( object = obj_QC, assay = "RNA", layer = counts_layer, fast = FALSE )
if ( nrow(counts_all) != nrow(obj_QC) || ncol(counts_all) != ncol(obj_QC) ) {
  stop_with_report( "Counts-layer dimensions differ from the Seurat object." )
}
if (!identical(colnames(counts_all), colnames(obj_QC))) {
  stop_with_report( "Counts-layer barcodes differ from the Seurat object." )
}
input_audit <- data.frame( item = c( "Input object", "Features", "QC-pass nuclei", "Biological libraries", "Timepoints", "Counts layer", "All qc_pass" ), value = c( normalizePath( input_object_path, winslash = "/", mustWork = TRUE ), nrow(obj_QC), ncol(obj_QC), length(libraries), length(observed_timepoints), counts_layer, all(as.logical(obj_QC$qc_pass)) ), stringsAsFactors = FALSE )
write_csv_safe( input_audit, "01_input_object_audit.csv" )
library_manifest <- unique( obj_QC@meta.data[ , c( "biological_library", "timepoint", "replicate_id" ), drop = FALSE ] )
library_manifest$timepoint <- factor( as.character(library_manifest$timepoint), levels = expected_timepoints, ordered = TRUE )
library_manifest <- library_manifest[ order( library_manifest$timepoint, suppressWarnings( as.numeric( as.character( library_manifest$replicate_id ) ) ), library_manifest$biological_library ), , drop = FALSE ]
libraries <- as.character( library_manifest$biological_library )
pre_doublet_library_counts <- as.data.frame( table( biological_library = as.character( obj_QC$biological_library ) ), stringsAsFactors = FALSE )
names(pre_doublet_library_counts)[ names(pre_doublet_library_counts) == "Freq" ] <- "n_QC_cells"
pre_doublet_library_counts <- merge( library_manifest, pre_doublet_library_counts, by = "biological_library", all.x = TRUE, sort = FALSE )
pre_doublet_library_counts <- pre_doublet_library_counts[ match( libraries, pre_doublet_library_counts$biological_library ), , drop = FALSE ]
pre_doublet_library_counts$expected_doublet_rate <- vapply( pre_doublet_library_counts$n_QC_cells, expected_doublet_rate, numeric(1) )
write_csv_safe( pre_doublet_library_counts, "02_pre_doublet_library_manifest.csv" )
result_list <- vector( mode = "list", length = length(libraries) )
names(result_list) <- libraries
error_list <- list()
for (i in seq_along(libraries)) {
  lib <- libraries[[i]]
  lib_cells <- rownames(obj_QC@meta.data)[ as.character(obj_QC$biological_library) == lib ]
  checkpoint_path <- file.path( checkpoint_dir, paste0( sprintf("%02d", i), "_", sanitize_filename(lib), "_scDblFinder.tsv.gz" ) )
  if ( !isTRUE(FORCE_RERUN_ALL_LIBRARIES) && checkpoint_is_usable( checkpoint_path, expected_cells = lib_cells, expected_library = lib ) ) {
    result_i <- data.table::fread( checkpoint_path, data.table = FALSE, showProgress = FALSE )
    result_i <- result_i[ match(lib_cells, result_i$barcode), , drop = FALSE ]
    result_list[[lib]] <- result_i
    next
  }
  result_i <- tryCatch( {
      cell_idx <- match( lib_cells, colnames(counts_all) )
      if (anyNA(cell_idx)) {
        stop( "At least one library barcode was not found in the counts layer." )
      }
      lib_view <- counts_all[ , cell_idx, drop = FALSE ]
      lib_counts <- methods::as( lib_view, "dgCMatrix" )
      rownames(lib_counts) <- rownames(lib_view)
      colnames(lib_counts) <- colnames(lib_view)
      if ( ncol(lib_counts) != length(lib_cells) || !identical( colnames(lib_counts), lib_cells ) ) {
        stop( "In-memory sparse matrix does not match the expected library cells." )
      }
      keep_genes <- Matrix::rowSums( lib_counts > 0 ) >= MIN_CELLS_PER_GENE
      if (sum(keep_genes) < 500L) {
        stop( paste0( "Only ", sum(keep_genes), " genes are detected in at least ", MIN_CELLS_PER_GENE, " nuclei." ) )
      }
      lib_counts <- lib_counts[ keep_genes, , drop = FALSE ]
      if ( any(Matrix::colSums(lib_counts) <= 0) ) {
        stop( "At least one nucleus has zero counts after gene filtering." )
      }
      sce <- SingleCellExperiment::SingleCellExperiment( assays = list( counts = lib_counts ) )
      colnames(sce) <- lib_cells
      colData(sce)$biological_library <- lib
      set.seed(SEED + i)
      sce <- scDblFinder::scDblFinder( sce, clusters = USE_CLUSTERS, dbr = NULL, dbr.per1k = DBR_PER_1K, BPPARAM = BiocParallel::SerialParam( progressbar = FALSE ), verbose = TRUE )
      score <- as.numeric( colData(sce)$scDblFinder.score )
      class <- as.character( colData(sce)$scDblFinder.class )
      if ( length(score) != length(lib_cells) || length(class) != length(lib_cells) || anyNA(score) || anyNA(class) ) {
        stop( "scDblFinder returned incomplete scores/classes." )
      }
      result <- data.frame( barcode = colnames(sce), biological_library = lib, scDblFinder.score = score, scDblFinder.class = class, n_genes_used = sum(keep_genes), n_library_cells = length(lib_cells), expected_doublet_rate = expected_doublet_rate( length(lib_cells) ), dbr_per_1k = DBR_PER_1K, scDblFinder_version = safe_package_version( "scDblFinder" ), seed = SEED + i, stringsAsFactors = FALSE )
      result <- result[ match(lib_cells, result$barcode), , drop = FALSE ]
      data.table::fwrite( result, file = checkpoint_path, sep = "\t", quote = FALSE, na = "NA", compress = "gzip" )
      rm( lib_view, lib_counts, sce )
      invisible(gc())
      result
    },
    error = function(e) {
      error_list[[length(error_list) + 1L]] <<- data.frame( biological_library = lib, n_cells = length(lib_cells), error = conditionMessage(e), stringsAsFactors = FALSE )
      NULL
    }
  )
  result_list[[lib]] <- result_i
}
valid_results <- result_list[ !vapply( result_list, is.null, logical(1) ) ]
if (length(error_list) > 0L) {
  errors <- do.call( rbind, error_list )
} else {
  errors <- data.frame( biological_library = character(0), n_cells = integer(0), error = character(0), stringsAsFactors = FALSE )
}
write_csv_safe( errors, "03_scDblFinder_errors.csv" )
if (length(valid_results) != length(libraries)) {
  stop_with_report( paste0( "scDblFinder completed for ", length(valid_results), "/", length(libraries), " libraries. Review tables/03_scDblFinder_errors.csv." ) )
}
doublet_results <- do.call( rbind, valid_results )
rownames(doublet_results) <- NULL
if ( nrow(doublet_results) != ncol(obj_QC) || anyDuplicated(doublet_results$barcode) > 0L || !setequal( doublet_results$barcode, colnames(obj_QC) ) ) {
  stop_with_report( "Combined scDblFinder results map one-to-one to all QC nuclei." )
}
doublet_results <- doublet_results[ match( colnames(obj_QC), doublet_results$barcode ), , drop = FALSE ]
if (!identical( doublet_results$barcode, colnames(obj_QC) )) {
  stop_with_report( "Combined scDblFinder result order does not match the Seurat object." )
}
doublet_metadata <- data.frame( scDblFinder.score = doublet_results$scDblFinder.score, scDblFinder.class = doublet_results$scDblFinder.class, stringsAsFactors = FALSE, row.names = doublet_results$barcode )
obj_annotated <- Seurat::AddMetaData( obj_QC, metadata = doublet_metadata )
obj_annotated$doublet_method <- "scDblFinder_samplewise_random"
obj_annotated$doublet_unit <- "biological_library"
obj_annotated$doublet_dbr_per_1k <- DBR_PER_1K
if ( anyNA(obj_annotated$scDblFinder.score) || anyNA(obj_annotated$scDblFinder.class) ) {
  stop_with_report( "The doublet metadata was not added completely to the Seurat object." )
}
doublet_summary_by_library <- data.table::as.data.table( obj_annotated@meta.data )[ , .( n_QC_input = .N, n_singlet = sum( scDblFinder.class == "singlet" ), n_doublet = sum( scDblFinder.class == "doublet" ), doublet_rate = mean( scDblFinder.class == "doublet" ), median_doublet_score = median( scDblFinder.score, na.rm = TRUE ), q95_doublet_score = quantile( scDblFinder.score, 0.95, na.rm = TRUE ), maximum_doublet_score = max( scDblFinder.score, na.rm = TRUE ) ), by = .( biological_library, timepoint, replicate_id ) ]
doublet_summary_by_library$expected_doublet_rate <- expected_doublet_rate( doublet_summary_by_library$n_QC_input )
doublet_summary_by_library$observed_expected_ratio <- doublet_summary_by_library$doublet_rate / doublet_summary_by_library$expected_doublet_rate
doublet_summary_by_library$rate_review_flag <- with( doublet_summary_by_library, doublet_rate > expected_doublet_rate * 2 + 0.02 | doublet_rate < pmax(expected_doublet_rate * 0.25, 0.001) )
doublet_summary_by_library <- as.data.frame( doublet_summary_by_library )
doublet_summary_by_library <- doublet_summary_by_library[ match( libraries, doublet_summary_by_library$biological_library ), , drop = FALSE ]
write_csv_safe( doublet_summary_by_library, "04_scDblFinder_summary_by_library.csv" )
doublet_summary_by_timepoint <- data.table::as.data.table( obj_annotated@meta.data )[ , .( n_QC_input = .N, n_singlet = sum( scDblFinder.class == "singlet" ), n_doublet = sum( scDblFinder.class == "doublet" ), doublet_rate = mean( scDblFinder.class == "doublet" ), n_libraries = uniqueN( biological_library ) ), by = timepoint ]
write_csv_safe( as.data.frame( doublet_summary_by_timepoint ), "04B_scDblFinder_summary_by_timepoint.csv" )
complete_cell_status_path <- file.path( tables_dir, "05_cell_level_scDblFinder_status.tsv.gz" )
cell_status <- obj_annotated@meta.data
cell_status$barcode <- rownames(cell_status)
data.table::fwrite( cell_status, file = complete_cell_status_path, sep = "\t", quote = FALSE, na = "NA", compress = "gzip" )
p_rate <- ggplot( doublet_summary_by_library, aes( x = reorder( biological_library, doublet_rate ), y = doublet_rate ) ) + geom_col() + geom_point( aes( y = expected_doublet_rate ), size = 2 ) + coord_flip() + scale_y_continuous( labels = scales::percent ) + labs( x = NULL, y = "Doublet rate", title = "Observed scDblFinder doublet rates by library", subtitle = "Points indicate the 0.8% per 1000-cell expectation" ) + theme_bw()
ggsave( filename = file.path( figures_dir, "01_doublet_rate_by_library.pdf" ), plot = p_rate, width = 10, height = 8, limitsize = FALSE )
set.seed(SEED)
plot_idx <- unlist( lapply( split( seq_len(nrow(cell_status)), cell_status$biological_library ), function(idx) {
      sample( idx, min(length(idx), 3000L), replace = FALSE )
    }
  ), use.names = FALSE
)
plot_meta <- cell_status[ plot_idx, , drop = FALSE ]
p_score <- ggplot( plot_meta, aes( x = biological_library, y = scDblFinder.score, fill = scDblFinder.class ) ) + geom_boxplot( outlier.shape = NA ) + coord_flip() + labs( x = NULL, y = "scDblFinder score", fill = NULL, title = "Doublet score distributions by library" ) + theme_bw()
ggsave( filename = file.path( figures_dir, "02_doublet_score_by_library.pdf" ), plot = p_score, width = 11, height = 9, limitsize = FALSE )
p_complexity <- ggplot( plot_meta, aes( x = scDblFinder.class, y = nFeature_RNA ) ) + geom_violin( scale = "width", trim = TRUE ) + scale_y_log10() + facet_wrap( ~timepoint, scales = "free_y" ) + labs( x = NULL, y = "nFeature_RNA (log10)", title = "Transcriptomic complexity by doublet class" ) + theme_bw()
ggsave( filename = file.path( figures_dir, "03_nFeature_by_doublet_class.pdf" ), plot = p_complexity, width = 12, height = 8, limitsize = FALSE )
obj_annotated@misc$GSE234774_step9C <- list( method = "scDblFinder", mode = "sample-wise; clusters=FALSE", library_field = "biological_library", dbr_per_1k = DBR_PER_1K, package_version = safe_package_version( "scDblFinder" ), seed = SEED, completed_at = as.character(Sys.time()), note = paste( "Doublet detection only.", "No normalization, integration, clustering, annotation,", "signature transfer, pseudotime, or spatial analysis was run." ) )
saveRDS( obj_annotated, annotated_object_path, compress = FALSE )
singlet_cells <- colnames(obj_annotated)[ obj_annotated$scDblFinder.class == "singlet" ]
singlet_idx <- match( singlet_cells, colnames(counts_all) )
if (anyNA(singlet_idx)) {
  stop_with_report( "At least one singlet barcode is absent from the counts layer." )
}
if ( isTRUE(FORCE_REWRITE_SINGLET_MATRIX) && dir.exists(singlet_matrix_dir) ) {
  unlink( singlet_matrix_dir, recursive = TRUE, force = TRUE )
}
if (!dir.exists(singlet_matrix_dir)) {
  singlet_view <- counts_all[ , singlet_idx, drop = FALSE ]
  BPCells::write_matrix_dir( mat = singlet_view, dir = singlet_matrix_dir, compress = TRUE, overwrite = FALSE )
  rm(singlet_view)
  invisible(gc())
} else {
}
singlet_counts <- open_bpcells_safe( singlet_matrix_dir )
if (inherits(singlet_counts, "error")) {
  stop_with_report( paste0( "Singlet BPCells matrix cannot be opened: ", conditionMessage(singlet_counts), ". Set FORCE_REWRITE_SINGLET_MATRIX <- TRUE and rerun." ) )
}
if ( ncol(singlet_counts) != length(singlet_cells) || !identical( colnames(singlet_counts), singlet_cells ) ) {
  stop_with_report( paste( "The singlet BPCells matrix does not match the current singlet calls.", "Set FORCE_REWRITE_SINGLET_MATRIX <- TRUE and rerun." ) )
}
singlet_metadata <- obj_annotated@meta.data[ match( singlet_cells, rownames(obj_annotated@meta.data) ), , drop = FALSE ]
rownames(singlet_metadata) <- singlet_cells
obj_singlet <- Seurat::CreateSeuratObject( counts = singlet_counts, assay = "RNA", project = "GSE234774_timecourse_QC_singlets", meta.data = singlet_metadata, min.cells = 0, min.features = 0 )
Seurat::DefaultAssay(obj_singlet) <- "RNA"
obj_singlet@misc$GSE234774_step9C <- obj_annotated@misc$GSE234774_step9C
obj_singlet@misc$GSE234774_step9C$singlet_matrix_dir <- normalizePath( singlet_matrix_dir, winslash = "/", mustWork = TRUE )
saveRDS( obj_singlet, singlet_object_path, compress = FALSE )
minimum_singlets <- min( doublet_summary_by_library$n_singlet )
maximum_doublet_rate <- max( doublet_summary_by_library$doublet_rate )
global_doublet_rate <- mean( obj_annotated$scDblFinder.class == "doublet" )
qualified_timepoint_libraries <- data.table::as.data.table( doublet_summary_by_library )[ n_singlet >= MIN_SINGLET_CELLS_PER_LIBRARY, .( n_qualified_libraries = uniqueN( biological_library ) ), by = timepoint ]
qualified_timepoint_libraries <- merge( data.frame( timepoint = factor( expected_timepoints, levels = expected_timepoints, ordered = TRUE ) ), as.data.frame( qualified_timepoint_libraries ), by = "timepoint", all.x = TRUE, sort = FALSE )
qualified_timepoint_libraries$n_qualified_libraries[ is.na( qualified_timepoint_libraries$n_qualified_libraries ) ] <- 0L
write_csv_safe( qualified_timepoint_libraries, "06_timepoint_singlet_library_check.csv" )
annotated_reload <- tryCatch( readRDS(annotated_object_path), error = function(e) e )
singlet_reload <- tryCatch( readRDS(singlet_object_path), error = function(e) e )
annotated_reload_pass <- !inherits(annotated_reload, "error") && ncol(annotated_reload) == ncol(obj_QC) && all( c( "scDblFinder.score", "scDblFinder.class" ) %in% names(annotated_reload@meta.data) )
singlet_reload_pass <- !inherits(singlet_reload, "error") && ncol(singlet_reload) == length(singlet_cells) && nrow(singlet_reload) == nrow(obj_QC)
rm( annotated_reload, singlet_reload )
invisible(gc())
audit_checks <- data.frame( check = c( "Step 9B QC object contains 148,567 nuclei", "Twenty biological libraries were processed separately", "All libraries completed scDblFinder", "Every QC nucleus has one score and one class", "No library doublet rate exceeds 25%", "Every library retains at least 100 singlets", "Every timepoint retains at least two qualified libraries", "Singlet BPCells matrix matches singlet barcodes", "Doublet-annotated object saves and reloads", "Singlet-only object saves and reloads", "Normalization/PCA/UMAP/clustering were not run", "Cohort age has been confirmed from study methods" ), required_for_pass = c( TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE ), status = c( ifelse( ncol(obj_QC) == 148567L, "PASS", "STOP" ), ifelse( length(libraries) == 20L, "PASS", "STOP" ), ifelse( length(valid_results) == 20L, "PASS", "STOP" ), ifelse( nrow(doublet_results) == ncol(obj_QC) && !anyNA( doublet_results$scDblFinder.score ) && !anyNA( doublet_results$scDblFinder.class ), "PASS", "STOP" ), ifelse( maximum_doublet_rate <= MAX_ACCEPTABLE_DOUBLET_RATE, "PASS", "STOP" ), ifelse( minimum_singlets >= MIN_SINGLET_CELLS_PER_LIBRARY, "PASS", "STOP" ), ifelse( all( qualified_timepoint_libraries$n_qualified_libraries >= MIN_QUALIFIED_LIBRARIES_PER_TIMEPOINT ), "PASS", "STOP" ), ifelse( ncol(singlet_counts) == length(singlet_cells) && identical( colnames(singlet_counts), singlet_cells ), "PASS", "STOP" ), ifelse( annotated_reload_pass, "PASS", "STOP" ), ifelse( singlet_reload_pass, "PASS", "STOP" ), "CONFIRMED_NOT_RUN", "REVIEW" ), detail = c( paste0( "detected=", ncol(obj_QC) ), paste0( "libraries=", length(libraries) ), paste0( length(valid_results), "/20 completed" ), paste0( "scored nuclei=", nrow(doublet_results) ), paste0( "maximum library doublet rate=", round( 100 * maximum_doublet_rate, 2 ), "%" ), paste0( "minimum singlets in one library=", minimum_singlets ), paste( paste0( qualified_timepoint_libraries$timepoint, ":", qualified_timepoint_libraries$n_qualified_libraries ), collapse = ";" ), paste0( nrow(singlet_counts), " features x ", ncol(singlet_counts), " singlets" ), ifelse( annotated_reload_pass, "RDS reload PASS", "RDS reload failed" ), ifelse( singlet_reload_pass, "RDS reload PASS", "RDS reload failed" ), "No NormalizeData, SCTransform, ScaleData, PCA, UMAP, tSNE, Harmony, neighbors, or clustering calls are present.", "Age remains METHODS_REVIEW_REQUIRED." ), stringsAsFactors = FALSE )
write_csv_safe( audit_checks, "09_STEP9C_audit_checks.csv" )
required_fail <- audit_checks$required_for_pass & audit_checks$status != "PASS"
final_decision <- if (any(required_fail)) {
  "STOP"
} else {
  "PASS_WITH_RATE_REVIEW"
}
library_lines <- paste0( doublet_summary_by_library$biological_library, " [", doublet_summary_by_library$timepoint, "]: ", doublet_summary_by_library$n_QC_input, " QC nuclei -> ", doublet_summary_by_library$n_singlet, " singlets + ", doublet_summary_by_library$n_doublet, " doublets; rate=", round( 100 * doublet_summary_by_library$doublet_rate, 2 ), "%" )
timepoint_lines <- paste0( doublet_summary_by_timepoint$timepoint, ": ", doublet_summary_by_timepoint$n_QC_input, " QC nuclei -> ", doublet_summary_by_timepoint$n_singlet, " singlets + ", doublet_summary_by_timepoint$n_doublet, " doublets; rate=", round( 100 * doublet_summary_by_timepoint$doublet_rate, 2 ), "%" )
rate_review_libraries <- doublet_summary_by_library$biological_library[ doublet_summary_by_library$rate_review_flag ]
if (final_decision == "STOP") {
  failed <- audit_checks[ audit_checks$required_for_pass & audit_checks$status != "PASS", , drop = FALSE ]
  for (i in seq_len(nrow(failed))) {
  }
} else {
}
report_path <- file.path( output_dir, "10_STEP9C_FINAL_REPORT.txt" )
if (final_decision == "STOP") {
  stop( "GSE234774 Step 9C decision is STOP. Review the final report.", call. = FALSE )
}
