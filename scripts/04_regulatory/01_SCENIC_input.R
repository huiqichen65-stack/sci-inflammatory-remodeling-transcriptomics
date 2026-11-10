rm(list = ls())
gc()
options( stringsAsFactors = FALSE, width = 220, warn = 1 )
set.seed(20260801)
PROJECT_ROOT <- Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
OUTDIR <- file.path( PROJECT_ROOT, "04_results", "step7A_GSE162610_pySCENIC_input_environment_audit" )
INPUTDIR <- file.path(OUTDIR, "input")
TABDIR <- file.path(OUTDIR, "tables")
REPORTDIR <- file.path(OUTDIR, "reports")
PROBEDIR <- file.path(OUTDIR, "environment_probe")
for (dd in c(OUTDIR, INPUTDIR, TABDIR, REPORTDIR, PROBEDIR)) {
  dir.create( dd, recursive = TRUE, showWarnings = FALSE )
}
OUT_REPORT <- file.path( REPORTDIR, "Step7A_report.txt" )
OUT_SESSION <- file.path( REPORTDIR, "Step7A_sessionInfo.txt" )
OUT_ERROR <- file.path( REPORTDIR, "Step7A_ERROR.txt" )
INPUT_FILES <- c( full_clean_microglia = file.path( PROJECT_ROOT, "04_results", "step3D_GSE162610_time_DEG_pathway_confirmation", "GSE162610_step3D_time_DEG_pathway_confirmed.rds" ), step4A_scored_object = file.path( PROJECT_ROOT, "04_results", "step4A_GSE162610_repair_scar_candidate_program_discovery", "GSE162610_step4A_candidate_program_scored.rds" ), step5BR_cds = file.path( PROJECT_ROOT, "04_results", "step5BR_GSE162610_injury_response_trajectory_refinement", "GSE162610_step5BR_injury_monocle3_cds.rds" ), step5BR_seurat = file.path( PROJECT_ROOT, "04_results", "step5BR_GSE162610_injury_response_trajectory_refinement", "GSE162610_step5BR_injury_trajectory_seurat.rds" ), step5CR_results = file.path( PROJECT_ROOT, "04_results", "step5CR_GSE162610_Y29_prebranch_progressive_divergence", "GSE162610_step5CR_results.rds" ), step6A_results = file.path( PROJECT_ROOT, "04_results", "step6A_GSE162610_7dpi_branch_specific_DEG", "GSE162610_step6A_7dpi_branch_DEG_results.rds" ), step6B_results = file.path( PROJECT_ROOT, "04_results", "step6B_GSE162610_dynamic_branch_modules", "GSE162610_step6B_dynamic_modules_results.rds" ), step6CR_results = file.path( PROJECT_ROOT, "04_results", "step6CR_GSE162610_branch_DEG_GO_KEGG", "GSE162610_step6CR_branch_DEG_GO_KEGG_results.rds" ), step6E_results = file.path( PROJECT_ROOT, "04_results", "step6E_GSE162610_pseudobulk_ssGSEA", "GSE162610_step6E_pseudobulk_ssGSEA_results.rds" ) )
EXPECTED_STATE_COUNTS <- c( Shared_prebranch = 116L, Arm_1 = 2289L, Arm_2 = 1798L )
EXPECTED_TOTAL_CELLS <- sum(EXPECTED_STATE_COUNTS)
MIN_DETECTION_FRACTION <- 0.01
MIN_COUNTS_PER_DETECTED_FRACTION <- 3
MIN_TF_DETECTION_FRACTION <- 0.001
MIN_TF_TOTAL_COUNTS <- 10
STATE_COLORS <- c( Shared_prebranch = "#756BB1", Arm_1 = "#2F6B9A", Arm_2 = "#B64B45" )
RESOURCE_DIR_CANDIDATES <- unique( c( file.path(PROJECT_ROOT, "01_raw_data", "pySCENIC_resources"), file.path(PROJECT_ROOT, "02_resources", "pySCENIC"), file.path(PROJECT_ROOT, "resources", "pySCENIC"), PROJECT_ROOT ) )
REQUIRED_R_PACKAGES <- c( "Seurat", "SeuratObject", "Matrix" )
missing_r_packages <- REQUIRED_R_PACKAGES[ !vapply( REQUIRED_R_PACKAGES, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1) ) ]
if (length(missing_r_packages) > 0) {
  stop( paste0( "Missing required R packages: ", paste(missing_r_packages, collapse = ", "), ". Install them first; this script does not auto-install packages." ), call. = FALSE )
}
suppressPackageStartupMessages( library(Seurat) )
suppressPackageStartupMessages( library(Matrix) )
timestamp <- function() {
  format( Sys.time(), "%Y-%m-%d %H:%M:%S" )
}
write_utf8_csv <- function(x, path) {
  utils::write.csv( x, file = path, row.names = FALSE, fileEncoding = "UTF-8", na = "" )
}
capture_table <- function(x, n = 100L) {
  if (is.null(x) || nrow(x) == 0) {
    return("<empty>")
  }
  paste( capture.output( print( utils::head(x, n), row.names = FALSE ) ), collapse = "\n" )
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
  paste0( format( round(x, 3), nsmall = 3, trim = TRUE ), " ", units[[idx]] )
}
safe_normalize_path <- function(path) {
  normalizePath( path, winslash = "/", mustWork = FALSE )
}
abort_step <- function(message_text) {
  lines <- c( "======================================================================", "GSE162610 | STEP 7A ERROR", "======================================================================", paste0("Time: ", timestamp()), "", message_text )
  writeLines( lines, con = OUT_ERROR, useBytes = TRUE )
  stop( message_text, call. = FALSE )
}
get_rna_layer <- function( object, layer_name = c("counts", "data") ) {
  layer_name <- match.arg(layer_name)
  if (!inherits(object, "Seurat")) {
    stop( "Input object is not a Seurat object.", call. = FALSE )
  }
  if (!"RNA" %in% names(object@assays)) {
    stop( "Seurat object does not contain an RNA assay.", call. = FALSE )
  }
  result <- tryCatch( SeuratObject::LayerData( object = object, assay = "RNA", layer = layer_name ), error = function(e) {
      NULL
    }
  )
  if (!is.null(result) && nrow(result) > 0 && ncol(result) > 0) {
    return(result)
  }
  result <- tryCatch( Seurat::GetAssayData( object = object, assay = "RNA", slot = layer_name ), error = function(e) {
      NULL
    }
  )
  if (!is.null(result) && nrow(result) > 0 && ncol(result) > 0) {
    return(result)
  }
  joined_object <- tryCatch( SeuratObject::JoinLayers( object = object, assay = "RNA" ), error = function(e) {
      NULL
    }
  )
  if (!is.null(joined_object)) {
    result <- tryCatch( SeuratObject::LayerData( object = joined_object, assay = "RNA", layer = layer_name ), error = function(e) {
        NULL
      }
    )
  }
  if (is.null(result) || nrow(result) == 0 || ncol(result) == 0) {
    stop( paste0( "Unable to retrieve RNA ", layer_name, " layer from the Seurat object." ), call. = FALSE )
  }
  result
}
sparse_log_normalize <- function(count_matrix, scale_factor = 10000) {
  library_sizes <- Matrix::colSums(count_matrix)
  if (any(!is.finite(library_sizes)) || any(library_sizes <= 0)) {
    stop( "At least one frozen cell has zero or invalid RNA library size.", call. = FALSE )
  }
  normalized <- count_matrix %*%
    Matrix::Diagonal( x = scale_factor / library_sizes )
  normalized@x <- log1p(normalized@x)
  normalized
}
read_tf_file <- function(path) {
  if ( is.na(path) || !nzchar(path) || !file.exists(path) ) {
    return(character(0))
  }
  lines <- readLines( path, warn = FALSE, encoding = "UTF-8" )
  lines <- trimws(lines)
  lines <- lines[nzchar(lines)]
  lines <- lines[!grepl("^#", lines)]
  if (length(lines) == 0) {
    return(character(0))
  }
  if ( any( grepl( "<html|<!doctype|github", lines[seq_len(min(20L, length(lines)))], ignore.case = TRUE ) ) ) {
    return(character(0))
  }
  tf <- sub( "[,\t ].*$", "", lines )
  unique(tf[nzchar(tf)])
}
is_probable_html <- function(path) {
  if ( !file.exists(path) || file.info(path)$size == 0 ) {
    return(FALSE)
  }
  lines <- tryCatch( readLines( path, n = 20L, warn = FALSE, encoding = "UTF-8" ), error = function(e) {
      character(0)
    }
  )
  any( grepl( "<html|<!doctype|github", lines, ignore.case = TRUE ) )
}
find_resource_files <- function(directories) {
  existing_dirs <- directories[ dir.exists(directories) ]
  if (length(existing_dirs) == 0) {
    return(character(0))
  }
  all_files <- unique( unlist( lapply( existing_dirs, function(dd) {
          list.files( dd, recursive = TRUE, full.names = TRUE, include.dirs = FALSE )
        }
      ), use.names = FALSE
    )
  )
  all_files[file.exists(all_files)]
}
clean_python_output <- function(x) {
  x <- as.character(x)
  x <- x[nzchar(trimws(x))]
  x
}
system2_capture <- function(command, args = character(0)) {
  out <- tryCatch( suppressWarnings( system2( command = command, args = args, stdout = TRUE, stderr = TRUE ) ), error = function(e) {
      structure( paste0("SYSTEM2_ERROR: ", conditionMessage(e)), status = 999L )
    }
  )
  status <- attr(out, "status")
  if (is.null(status)) {
    status <- 0L
  }
  list( status = as.integer(status), output = clean_python_output(out) )
}
input_file_audit <- data.frame( input_name = names(INPUT_FILES), path = safe_normalize_path(unname(INPUT_FILES)), exists = file.exists(unname(INPUT_FILES)), size_bytes = ifelse( file.exists(unname(INPUT_FILES)), file.info(unname(INPUT_FILES))$size, NA_real_ ), stringsAsFactors = FALSE )
input_file_audit$size_human <- vapply( input_file_audit$size_bytes, format_bytes, FUN.VALUE = character(1) )
write_utf8_csv( input_file_audit, file.path( TABDIR, "01_input_file_audit.csv" ) )
if (any(!input_file_audit$exists)) {
  missing_paths <- input_file_audit$path[ !input_file_audit$exists ]
  abort_step( paste0( "Formal Step 7 input files are missing:\n", paste(missing_paths, collapse = "\n") ) )
}
full_obj <- readRDS( INPUT_FILES[["full_clean_microglia"]] )
if (!inherits(full_obj, "Seurat")) {
  abort_step( "The Step3D full clean microglia RDS is not a Seurat object." )
}
step5cr <- readRDS( INPUT_FILES[["step5CR_results"]] )
step6b <- readRDS( INPUT_FILES[["step6B_results"]] )
required_step5cr <- c( "prebranch_cells", "arm1_cells", "arm2_cells", "analysis_metadata", "progressive_shape_genes", "strong_progressive_genes" )
missing_step5cr <- setdiff( required_step5cr, names(step5cr) )
if (length(missing_step5cr) > 0) {
  abort_step( paste0( "Step5C-R result object is missing components: ", paste(missing_step5cr, collapse = ", ") ) )
}
required_step6b <- c( "module_annotation", "gene_modules", "selected_k" )
missing_step6b <- setdiff( required_step6b, names(step6b) )
if (length(missing_step6b) > 0) {
  abort_step( paste0( "Step6B result object is missing components: ", paste(missing_step6b, collapse = ", ") ) )
}
object_component_audit <- data.frame( object = c( rep("Step5C-R", length(names(step5cr))), rep("Step6B", length(names(step6b))) ), component = c( names(step5cr), names(step6b) ), class = c( vapply( step5cr, function(x) paste(class(x), collapse = ";"), FUN.VALUE = character(1) ), vapply( step6b, function(x) paste(class(x), collapse = ";"), FUN.VALUE = character(1) ) ), stringsAsFactors = FALSE )
write_utf8_csv( object_component_audit, file.path( TABDIR, "02_object_component_audit.csv" ) )
analysis_meta <- as.data.frame( step5cr[["analysis_metadata"]] )
required_analysis_meta <- c( "cell", "state", "timepoint", "sample", "retained_step5BR_pseudotime", "graph_distance_from_Y29", "signed_progress", "post_progress", "branch_contrast" )
missing_analysis_meta <- setdiff( required_analysis_meta, colnames(analysis_meta) )
if (length(missing_analysis_meta) > 0) {
  abort_step( paste0( "Step5C-R analysis_metadata is missing columns: ", paste(missing_analysis_meta, collapse = ", ") ) )
}
analysis_meta <- analysis_meta[ analysis_meta$state %in% names(EXPECTED_STATE_COUNTS), , drop = FALSE ]
if (anyDuplicated(analysis_meta$cell) > 0) {
  abort_step( "Step5C-R analysis_metadata contains duplicated cell IDs." )
}
state_counts <- table( factor( analysis_meta$state, levels = names(EXPECTED_STATE_COUNTS) ) )
state_count_audit <- data.frame( state = names(EXPECTED_STATE_COUNTS), expected_n = as.integer(EXPECTED_STATE_COUNTS), observed_n = as.integer(state_counts), exact_match = as.integer(state_counts) == as.integer(EXPECTED_STATE_COUNTS), stringsAsFactors = FALSE )
write_utf8_csv( state_count_audit, file.path( TABDIR, "03_frozen_state_count_audit.csv" ) )
if ( nrow(analysis_meta) != EXPECTED_TOTAL_CELLS || any(!state_count_audit$exact_match) ) {
  abort_step( paste0( "Frozen Step5C-R cell counts exactly match the formal set.\n", capture_table(state_count_audit) ) )
}
frozen_cells <- as.character( analysis_meta$cell )
if ( !setequal( frozen_cells, c( as.character(step5cr[["prebranch_cells"]]), as.character(step5cr[["arm1_cells"]]), as.character(step5cr[["arm2_cells"]]) ) ) ) {
  abort_step( paste0( "analysis_metadata cells exactly match ", "prebranch_cells + arm1_cells + arm2_cells." ) )
}
missing_in_seurat <- setdiff( frozen_cells, colnames(full_obj) )
if (length(missing_in_seurat) > 0) {
  abort_step( paste0( length(missing_in_seurat), " frozen cells are absent from the Step3D Seurat object." ) )
}
counts_all <- get_rna_layer( full_obj, layer_name = "counts" )
if (!inherits(counts_all, "sparseMatrix")) {
  counts_all <- methods::as( counts_all, "dgCMatrix" )
}
if (anyDuplicated(rownames(counts_all)) > 0) {
  abort_step( paste0( "RNA count matrix contains duplicated gene names. ", "use make.unique() because it would break motif/TF matching." ) )
}
counts_frozen <- counts_all[ , frozen_cells, drop = FALSE ]
counts_frozen <- methods::as( counts_frozen, "dgCMatrix" )
if ( !identical( colnames(counts_frozen), frozen_cells ) ) {
  abort_step( "Raw-count cell order does not match the frozen metadata order." )
}
if ( length(counts_frozen@x) > 0 && ( any(!is.finite(counts_frozen@x)) || any(counts_frozen@x < 0) || any(abs(counts_frozen@x - round(counts_frozen@x)) > 1e-8) ) ) {
  abort_step( "RNA counts are not non-negative integer-like raw counts." )
}
rm(counts_all)
gc()
seurat_meta <- full_obj@meta.data
seurat_meta$cell <- rownames(seurat_meta)
seurat_meta_frozen <- seurat_meta[ match( frozen_cells, seurat_meta$cell ), , drop = FALSE ]
cell_metadata <- analysis_meta[ match( frozen_cells, analysis_meta$cell ), required_analysis_meta, drop = FALSE ]
colnames(cell_metadata)[ colnames(cell_metadata) == "state" ] <- "step5CR_state"
colnames(cell_metadata)[ colnames(cell_metadata) == "sample" ] <- "sample_id"
colnames(cell_metadata)[ colnames(cell_metadata) == "timepoint" ] <- "timepoint"
cell_metadata$sample_7dpi <- ifelse( cell_metadata$timepoint == "7 dpi", cell_metadata$sample_id, NA_character_ )
cell_metadata$Y29_signed_progress <- cell_metadata$signed_progress
cell_metadata$Y29_postbranch_progress <- cell_metadata$post_progress
seurat_keep_patterns <- c( "^orig.ident$", "^sample_id$", "^sample$", "^timepoint$", "^time$", "^dpi$", "^microglia_cluster", "^state", "^shared_transition", "^step3", "^step5" )
seurat_keep <- unique( unlist( lapply( seurat_keep_patterns, function(pp) {
        grep( pp, colnames(seurat_meta_frozen), value = TRUE, ignore.case = TRUE )
      }
    ), use.names = FALSE
  )
)
seurat_keep <- setdiff( seurat_keep, "cell" )
for (cc in seurat_keep) {
  new_name <- paste0( "seurat__", cc )
  if (!new_name %in% colnames(cell_metadata)) {
    cell_metadata[[new_name]] <- seurat_meta_frozen[[cc]]
  }
}
gene_modules <- as.data.frame( step6b[["gene_modules"]] )
required_gene_module_cols <- c( "gene_id", "module_name" )
missing_gene_module_cols <- setdiff( required_gene_module_cols, colnames(gene_modules) )
if (length(missing_gene_module_cols) > 0) {
  abort_step( paste0( "Step6B gene_modules is missing columns: ", paste(missing_gene_module_cols, collapse = ", ") ) )
}
norm_frozen <- tryCatch( {
    norm_all <- get_rna_layer( full_obj, layer_name = "data" )
    norm_use <- norm_all[ , frozen_cells, drop = FALSE ]
    rm(norm_all)
    gc()
    methods::as( norm_use, "dgCMatrix" )
  },
  error = function(e) {
    sparse_log_normalize( counts_frozen )
  }
)
module_names <- unique( as.character( gene_modules$module_name ) )
module_names <- module_names[ !is.na(module_names) & nzchar(module_names) ]
module_score_audit <- list()
for (module_name_current in module_names) {
  module_genes <- unique( as.character( gene_modules$gene_id[ gene_modules$module_name == module_name_current ] ) )
  available_genes <- intersect( module_genes, rownames(norm_frozen) )
  score_col <- paste0( "step6B_mean_norm_score__", make.names(module_name_current) )
  if (length(available_genes) == 0) {
    cell_metadata[[score_col]] <- NA_real_
  } else {
    cell_metadata[[score_col]] <- as.numeric( Matrix::colMeans( norm_frozen[ available_genes, frozen_cells, drop = FALSE ] ) )
  }
  module_score_audit[[module_name_current]] <- data.frame( module_name = module_name_current, genes_in_step6B = length(module_genes), genes_available = length(available_genes), availability_fraction = ifelse( length(module_genes) > 0, length(available_genes) / length(module_genes), NA_real_ ), score_column = score_col, score_method = paste0( "Per-cell mean of RNA log-normalized expression; ", "for metadata/audit only; not used to define branches" ), stringsAsFactors = FALSE )
}
module_score_audit <- do.call( rbind, module_score_audit )
rownames(module_score_audit) <- NULL
write_utf8_csv( module_score_audit, file.path( TABDIR, "04_step6B_module_score_audit.csv" ) )
resource_files <- find_resource_files( RESOURCE_DIR_CANDIDATES )
resource_basename <- basename(resource_files)
mouse_feather_candidates <- resource_files[ grepl( "\\.feather$", resource_basename, ignore.case = TRUE ) & grepl( "mm10|mm9|mus_musculus", resource_basename, ignore.case = TRUE ) ]
mouse_annotation_candidates <- resource_files[ grepl( "\\.tbl$", resource_basename, ignore.case = TRUE ) & grepl( "mgi|mus_musculus|mouse", resource_basename, ignore.case = TRUE ) ]
mouse_tf_candidates <- resource_files[ grepl( "(mm.*tfs|mgi.*tfs|mouse.*tfs).*\\.txt$", resource_basename, ignore.case = TRUE ) ]
resource_category <- function(path) {
  bn <- basename(path)
  if ( grepl("\\.feather$", bn, ignore.case = TRUE) ) {
    return("mouse_cisTarget_feather")
  }
  if ( grepl("\\.tbl$", bn, ignore.case = TRUE) ) {
    return("mouse_motif_annotation")
  }
  "mouse_TF_list"
}
resource_audit_paths <- unique( c( mouse_feather_candidates, mouse_annotation_candidates, mouse_tf_candidates ) )
if (length(resource_audit_paths) > 0) {
  cistarget_audit <- data.frame( category = vapply( resource_audit_paths, resource_category, FUN.VALUE = character(1) ), path = safe_normalize_path(resource_audit_paths), exists = file.exists(resource_audit_paths), readable = file.access(resource_audit_paths, 4L) == 0L, size_bytes = file.info(resource_audit_paths)$size, stringsAsFactors = FALSE )
  cistarget_audit$size_human <- vapply( cistarget_audit$size_bytes, format_bytes, FUN.VALUE = character(1) )
  cistarget_audit$probable_html <- vapply( resource_audit_paths, function(pp) {
      if ( grepl( "\\.(txt|tbl)$", pp, ignore.case = TRUE ) ) {
        is_probable_html(pp)
      } else {
        FALSE
      }
    },
    FUN.VALUE = logical(1)
  )
  cistarget_audit$mouse_name_match <- grepl( "mm10|mm9|mus_musculus|mgi|mouse", basename(resource_audit_paths), ignore.case = TRUE )
} else {
  cistarget_audit <- data.frame( category = c( "mouse_cisTarget_feather", "mouse_motif_annotation", "mouse_TF_list" ), path = NA_character_, exists = FALSE, readable = FALSE, size_bytes = NA_real_, size_human = NA_character_, probable_html = FALSE, mouse_name_match = FALSE, stringsAsFactors = FALSE )
}
write_utf8_csv( cistarget_audit, file.path( TABDIR, "05_cistarget_database_audit.csv" ) )
valid_tf_candidates <- mouse_tf_candidates[ file.exists(mouse_tf_candidates) & file.access(mouse_tf_candidates, 4L) == 0L & !vapply( mouse_tf_candidates, is_probable_html, FUN.VALUE = logical(1) ) ]
TF_FILE_SELECTED <- if ( length(valid_tf_candidates) > 0 ) {
  valid_tf_candidates[ which.max( file.info(valid_tf_candidates)$size ) ][[1]]
} else {
  NA_character_
}
mouse_tfs <- read_tf_file( TF_FILE_SELECTED )
if (length(mouse_tfs) > 0) {
} else {
}
detected_cells <- Matrix::rowSums( counts_frozen > 0 )
total_counts <- Matrix::rowSums( counts_frozen )
n_frozen <- ncol(counts_frozen)
standard_min_cells <- max( 10L, ceiling( MIN_DETECTION_FRACTION * n_frozen ) )
standard_min_counts <- max( 30L, ceiling( MIN_COUNTS_PER_DETECTED_FRACTION * MIN_DETECTION_FRACTION * n_frozen ) )
tf_min_cells <- max( 3L, ceiling( MIN_TF_DETECTION_FRACTION * n_frozen ) )
standard_keep <- detected_cells >= standard_min_cells & total_counts >= standard_min_counts
is_mouse_tf <- rownames(counts_frozen) %in%
  mouse_tfs
tf_rescue <- is_mouse_tf & detected_cells >= tf_min_cells & total_counts >= MIN_TF_TOTAL_COUNTS
final_keep <- standard_keep | tf_rescue
gene_names <- rownames(counts_frozen)
gene_category <- rep( "other", length(gene_names) )
gene_category[ grepl( "^mt-", gene_names, ignore.case = TRUE ) ] <- "mitochondrial"
gene_category[ grepl( "^Rpl|^Rps", gene_names ) ] <- "ribosomal"
gene_category[ grepl( "^Hba|^Hbb", gene_names ) ] <- "hemoglobin"
immediate_early_candidates <- c( "Jun", "Junb", "Jund", "Fos", "Fosb", "Fosl1", "Fosl2", "Egr1", "Egr2", "Egr3", "Atf3", "Nfkb1", "Nfkb2", "Rela", "Relb", "Nfkbia", "Nfkbiz" )
gene_filter_audit <- data.frame( gene = gene_names, detected_cells = as.integer(detected_cells), detection_fraction = as.numeric( detected_cells / n_frozen ), total_counts = as.numeric(total_counts), gene_category = gene_category, in_mouse_TF_list = is_mouse_tf, standard_keep = standard_keep, TF_rescue = tf_rescue, final_keep = final_keep, immediate_early_candidate = gene_names %in% immediate_early_candidates, stringsAsFactors = FALSE )
write_utf8_csv( gene_filter_audit, file.path( TABDIR, "06_gene_filter_and_TF_audit.csv" ) )
gene_filter_summary <- data.frame( metric = c( "Frozen cells", "Genes before filter", "Standard minimum detected cells", "Standard minimum total counts", "TF rescue minimum detected cells", "TF rescue minimum total counts", "Genes passing standard filter", "TF genes in raw matrix", "TF genes passing standard filter", "TF genes rescued", "Genes retained final", "Mitochondrial genes retained", "Ribosomal genes retained", "Immediate-early candidates retained" ), value = c( n_frozen, nrow(counts_frozen), standard_min_cells, standard_min_counts, tf_min_cells, MIN_TF_TOTAL_COUNTS, sum(standard_keep), sum(is_mouse_tf), sum(is_mouse_tf & standard_keep), sum(tf_rescue & !standard_keep), sum(final_keep), sum(final_keep & gene_category == "mitochondrial"), sum(final_keep & gene_category == "ribosomal"), sum( final_keep & gene_names %in% immediate_early_candidates ) ), stringsAsFactors = FALSE )
write_utf8_csv( gene_filter_summary, file.path( TABDIR, "07_gene_filter_summary.csv" ) )
if (sum(final_keep) < 1000) {
  abort_step( paste0( "Fewer than 1,000 genes remain after filtering. ", "proceed to GRN." ) )
}
counts_filtered <- counts_frozen[ final_keep, , drop = FALSE ]
counts_filtered <- methods::as( counts_filtered, "dgCMatrix" )
cell_metadata <- cell_metadata[ match( colnames(counts_filtered), cell_metadata$cell ), , drop = FALSE ]
if ( !identical( as.character(cell_metadata$cell), colnames(counts_filtered) ) ) {
  abort_step( "Final metadata order does not match filtered count-matrix columns." )
}
module_score_cols <- grep( "^step6B_mean_norm_score__", colnames(cell_metadata), value = TRUE )
metadata_required_fields <- c( "cell", "sample_id", "timepoint", "step5CR_state", "graph_distance_from_Y29", "Y29_signed_progress", "Y29_postbranch_progress", "retained_step5BR_pseudotime", "sample_7dpi" )
metadata_audit <- data.frame( field = c( metadata_required_fields, module_score_cols ), exists = c( metadata_required_fields, module_score_cols ) %in% colnames(cell_metadata), n_missing = vapply( c( metadata_required_fields, module_score_cols ), function(cc) {
      if (!cc %in% colnames(cell_metadata)) {
        return(NA_integer_)
      }
      sum( is.na(cell_metadata[[cc]]) | ( is.character(cell_metadata[[cc]]) & !nzchar(cell_metadata[[cc]]) ) )
    },
    FUN.VALUE = integer(1)
  ), stringsAsFactors = FALSE
)
write_utf8_csv( metadata_audit, file.path( TABDIR, "08_cell_metadata_field_audit.csv" ) )
critical_metadata_fields <- c( "cell", "sample_id", "timepoint", "step5CR_state", "Y29_signed_progress", "retained_step5BR_pseudotime" )
critical_metadata_pass <- all( critical_metadata_fields %in% colnames(cell_metadata) ) && all( vapply( critical_metadata_fields, function(cc) {
        all( !is.na( cell_metadata[[cc]] ) )
      },
      FUN.VALUE = logical(1)
    )
  )
module_metadata_pass <- length(module_score_cols) == length(module_names) && all( vapply( module_score_cols, function(cc) {
        all( is.finite( cell_metadata[[cc]] ) )
      },
      FUN.VALUE = logical(1)
    )
  )
write_utf8_csv( cell_metadata, file.path( INPUTDIR, "GSE162610_step7A_frozen_4203_cell_metadata.csv" ) )
MTX_FILE <- file.path( INPUTDIR, "GSE162610_step7A_frozen_4203_raw_counts_genes_by_cells.mtx" )
GENE_FILE <- file.path( INPUTDIR, "GSE162610_step7A_frozen_4203_genes.tsv" )
BARCODE_FILE <- file.path( INPUTDIR, "GSE162610_step7A_frozen_4203_barcodes.tsv" )
Matrix::writeMM( obj = counts_filtered, file = MTX_FILE )
utils::write.table( data.frame(gene = rownames(counts_filtered)), file = GENE_FILE, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE )
utils::write.table( data.frame(cell = colnames(counts_filtered)), file = BARCODE_FILE, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE )
INPUT_README <- file.path( INPUTDIR, "README_step7A_pySCENIC_input.txt" )
writeLines( c( "GSE162610 Step 7A pySCENIC input", "", paste0( "Matrix orientation: genes x cells" ), paste0( "Genes: ", nrow(counts_filtered) ), paste0( "Cells: ", ncol(counts_filtered) ), paste0( "Frozen states: Shared_prebranch=116; Arm_1=2289; Arm_2=1798" ), "", "Files:", basename(MTX_FILE), basename(GENE_FILE), basename(BARCODE_FILE), "GSE162610_step7A_frozen_4203_cell_metadata.csv", "", "The matrix contains sparse raw RNA counts.", "It is not log-normalized and was never converted to a dense CSV.", "", "Mitochondrial/ribosomal genes were not mechanically removed.", "Immediate-early TF candidates were retained if they passed expression rules.", "", "This input has not yet been converted to loom.", "Loom conversion and GRNBoost2 are deferred until environment/database PASS." ), con = INPUT_README, useBytes = TRUE )
input_manifest <- list( project = "GSE162610 adult SCI activated/DAM-like microglia", step = "Step 7A", generated_at = timestamp(), matrix_orientation = "genes x cells", frozen_state_counts = EXPECTED_STATE_COUNTS, matrix_file = MTX_FILE, gene_file = GENE_FILE, barcode_file = BARCODE_FILE, metadata_file = file.path( INPUTDIR, "GSE162610_step7A_frozen_4203_cell_metadata.csv" ), genes = rownames(counts_filtered), cells = colnames(counts_filtered), filtering_parameters = list( minimum_detection_fraction = MIN_DETECTION_FRACTION, minimum_total_counts = standard_min_counts, TF_rescue_minimum_detection_fraction = MIN_TF_DETECTION_FRACTION, TF_rescue_minimum_total_counts = MIN_TF_TOTAL_COUNTS ), selected_mouse_TF_file = TF_FILE_SELECTED, state_colors = STATE_COLORS )
saveRDS( input_manifest, file.path( INPUTDIR, "GSE162610_step7A_input_manifest.rds" ), compress = TRUE )
export_audit <- data.frame( file_role = c( "raw_counts_matrix", "genes", "barcodes", "cell_metadata", "input_readme", "input_manifest" ), path = safe_normalize_path( c( MTX_FILE, GENE_FILE, BARCODE_FILE, file.path( INPUTDIR, "GSE162610_step7A_frozen_4203_cell_metadata.csv" ), INPUT_README, file.path( INPUTDIR, "GSE162610_step7A_input_manifest.rds" ) ) ), exists = file.exists( c( MTX_FILE, GENE_FILE, BARCODE_FILE, file.path( INPUTDIR, "GSE162610_step7A_frozen_4203_cell_metadata.csv" ), INPUT_README, file.path( INPUTDIR, "GSE162610_step7A_input_manifest.rds" ) ) ), size_bytes = file.info( c( MTX_FILE, GENE_FILE, BARCODE_FILE, file.path( INPUTDIR, "GSE162610_step7A_frozen_4203_cell_metadata.csv" ), INPUT_README, file.path( INPUTDIR, "GSE162610_step7A_input_manifest.rds" ) ) )$size, stringsAsFactors = FALSE )
export_audit$size_human <- vapply( export_audit$size_bytes, format_bytes, FUN.VALUE = character(1) )
write_utf8_csv( export_audit, file.path( TABDIR, "09_export_file_audit.csv" ) )
export_pass <- all(export_audit$exists)
PY_PROBE_FILE <- file.path( PROBEDIR, "python_environment_probe.py" )
writeLines( c( "import sys", "import importlib", "import importlib.metadata as metadata", "", "packages = [", " 'pyscenic',", " 'arboreto',", " 'ctxcore',", " 'pandas',", " 'numpy',", " 'loompy',", " 'dask',", " 'numba',", " 'scanpy',", " 'pyarrow'", "]", "", "print('python_executable\\t' + str(sys.executable))", "print('python_version\\t' + str(sys.version).replace('\\n', ' '))", "", "for package in packages:", " try:", " module = importlib.import_module(package)", " try:", " version = metadata.version(package)", " except Exception:", " version = getattr(module, '__version__', 'installed_version_unknown')", " print(package + '\\tPASS\\t' + str(version) + '\\t')", " except Exception as exc:", " text = str(exc).replace('\\t', ' ').replace('\\n', ' ')", " print(package + '\\tFAIL\\tNA\\t' + text)" ), con = PY_PROBE_FILE, useBytes = TRUE )
known_python_paths <- c( Sys.getenv("PYSCENIC_PYTHON"), unname(Sys.which("python")), unname(Sys.which("python3")), file.path( Sys.getenv("USERPROFILE"), "miniconda3", "envs", "pyscenic_env", "python.exe" ), file.path( Sys.getenv("USERPROFILE"), "anaconda3", "envs", "pyscenic_env", "python.exe" ), file.path( Sys.getenv("LOCALAPPDATA"), "miniconda3", "envs", "pyscenic_env", "python.exe" ) )
known_python_paths <- unique( known_python_paths[ nzchar(known_python_paths) ] )
known_python_paths <- known_python_paths[ file.exists(known_python_paths) ]
python_audit_list <- list()
if (length(known_python_paths) > 0) {
  for (python_path in known_python_paths) {
    probe <- system2_capture( command = python_path, args = shQuote(PY_PROBE_FILE) )
    probe_lines <- probe$output
    package_lines <- probe_lines[ grepl( "^(pyscenic|arboreto|ctxcore|pandas|numpy|loompy|dask|numba|scanpy|pyarrow)\\t", probe_lines ) ]
    if (length(package_lines) > 0) {
      split_rows <- strsplit( package_lines, "\t", fixed = TRUE )
      package_df <- do.call( rbind, lapply( split_rows, function(xx) {
            length(xx) <- 4L
            data.frame( package = xx[[1]], import_status = xx[[2]], version = xx[[3]], error = xx[[4]], stringsAsFactors = FALSE )
          }
        )
      )
    } else {
      package_df <- data.frame( package = c( "pyscenic", "arboreto", "ctxcore", "pandas", "numpy", "loompy", "dask", "numba", "scanpy", "pyarrow" ), import_status = "FAIL", version = NA_character_, error = paste( probe_lines, collapse = " | " ), stringsAsFactors = FALSE )
    }
    package_df$python_path <- safe_normalize_path( python_path )
    package_df$system2_status <- probe$status
    python_audit_list[[ safe_normalize_path(python_path) ]] <- package_df
  }
  python_environment_audit <- do.call( rbind, python_audit_list )
  rownames(python_environment_audit) <- NULL
} else {
  python_environment_audit <- data.frame( package = c( "pyscenic", "arboreto", "ctxcore", "pandas", "numpy", "loompy", "dask", "numba", "scanpy", "pyarrow" ), import_status = "FAIL", version = NA_character_, error = "No callable Python executable was found.", python_path = NA_character_, system2_status = 127L, stringsAsFactors = FALSE )
}
write_utf8_csv( python_environment_audit, file.path( TABDIR, "10_python_environment_audit.csv" ) )
required_python_packages <- c( "pyscenic", "arboreto", "ctxcore", "pandas", "numpy", "loompy" )
python_paths_tested <- unique( python_environment_audit$python_path[ !is.na( python_environment_audit$python_path ) ] )
python_path_pass <- setNames( logical(length(python_paths_tested)), python_paths_tested )
for (python_path in python_paths_tested) {
  sub_audit <- python_environment_audit[ python_environment_audit$python_path == python_path & python_environment_audit$package %in% required_python_packages, , drop = FALSE ]
  python_path_pass[[python_path]] <- nrow(sub_audit) == length(required_python_packages) && all( sub_audit$import_status == "PASS" ) && all( sub_audit$system2_status == 0L )
}
PYTHON_SELECTED <- if ( length(python_path_pass) > 0 && any(python_path_pass) ) {
  names(python_path_pass)[ which(python_path_pass)[[1]] ]
} else {
  NA_character_
}
pyscenic_cli_audit <- data.frame( python_path = PYTHON_SELECTED, cli_path = NA_character_, exists = FALSE, system2_status = NA_integer_, callable = FALSE, first_output_line = NA_character_, stringsAsFactors = FALSE )
if (!is.na(PYTHON_SELECTED)) {
  cli_candidates <- unique( c( file.path( dirname(PYTHON_SELECTED), "Scripts", "pyscenic.exe" ), file.path( dirname(PYTHON_SELECTED), "pyscenic.exe" ), unname( Sys.which("pyscenic") ) ) )
  cli_candidates <- cli_candidates[ nzchar(cli_candidates) & file.exists(cli_candidates) ]
  if (length(cli_candidates) > 0) {
    cli_path <- cli_candidates[[1]]
    cli_probe <- system2_capture( command = cli_path, args = "--help" )
    pyscenic_cli_audit <- data.frame( python_path = PYTHON_SELECTED, cli_path = safe_normalize_path( cli_path ), exists = TRUE, system2_status = cli_probe$status, callable = cli_probe$status == 0L, first_output_line = if ( length(cli_probe$output) > 0 ) {
        cli_probe$output[[1]]
      } else {
        NA_character_
      },
      stringsAsFactors = FALSE
    )
  }
}
write_utf8_csv( pyscenic_cli_audit, file.path( TABDIR, "11_pyscenic_CLI_audit.csv" ) )
python_import_pass <- !is.na( PYTHON_SELECTED )
python_cli_pass <- isTRUE( pyscenic_cli_audit$callable[[1]] )
valid_feather <- cistarget_audit[ cistarget_audit$category == "mouse_cisTarget_feather" & cistarget_audit$exists & cistarget_audit$readable & cistarget_audit$mouse_name_match & !cistarget_audit$probable_html, , drop = FALSE ]
valid_annotation <- cistarget_audit[ cistarget_audit$category == "mouse_motif_annotation" & cistarget_audit$exists & cistarget_audit$readable & cistarget_audit$mouse_name_match & !cistarget_audit$probable_html, , drop = FALSE ]
valid_tf_audit <- cistarget_audit[ cistarget_audit$category == "mouse_TF_list" & cistarget_audit$exists & cistarget_audit$readable & cistarget_audit$mouse_name_match & !cistarget_audit$probable_html, , drop = FALSE ]
database_pass <- nrow(valid_feather) >= 1L && nrow(valid_annotation) >= 1L && nrow(valid_tf_audit) >= 1L && length(mouse_tfs) >= 100L
dense_equivalent_bytes <- as.numeric(nrow(counts_filtered)) * as.numeric(ncol(counts_filtered)) * 8
sparse_object_bytes <- as.numeric( object.size(counts_filtered) )
resource_estimate <- data.frame( metric = c( "Cells", "Genes", "Non-zero entries", "Matrix density", "Sparse R object bytes", "Sparse R object human", "Dense double-matrix equivalent bytes", "Dense double-matrix equivalent human", "Approximate 3x dense working-memory lower guide", "Approximate 8x dense working-memory upper guide", "Logical CPU cores reported by R", "Project path contains spaces", "Maximum exported path length", "GRN runtime benchmark performed" ), value = c( ncol(counts_filtered), nrow(counts_filtered), length(counts_filtered@x), signif( length(counts_filtered@x) / ( as.numeric(nrow(counts_filtered)) * as.numeric(ncol(counts_filtered)) ), 6 ), sparse_object_bytes, format_bytes(sparse_object_bytes), dense_equivalent_bytes, format_bytes(dense_equivalent_bytes), format_bytes( 3 * dense_equivalent_bytes ), format_bytes( 8 * dense_equivalent_bytes ), parallel::detectCores( logical = TRUE ), grepl( "\\s", PROJECT_ROOT ), max( nchar( export_audit$path ) ), "No; Step 7A deliberately does not run GRNBoost2" ), stringsAsFactors = FALSE )
write_utf8_csv( resource_estimate, file.path( TABDIR, "12_resource_estimate.csv" ) )
input_matrix_pass <- ncol(counts_filtered) == EXPECTED_TOTAL_CELLS && nrow(counts_filtered) >= 1000L && export_pass
input_metadata_pass <- critical_metadata_pass && module_metadata_pass
tf_filter_pass <- length(mouse_tfs) >= 100L && sum( gene_filter_audit$in_mouse_TF_list & gene_filter_audit$final_keep ) >= 50L
input_pass <- input_matrix_pass && input_metadata_pass && tf_filter_pass
environment_pass <- python_import_pass && python_cli_pass
if ( input_pass && environment_pass && database_pass ) {
  FINAL_STATUS <- "PASS"
  NEXT_ACTION <- paste0( "Step 7A passed. Proceed to Step 7B GRNBoost2 planning, ", "using the frozen 4,203-cell input. change branch labels." )
} else if (
  input_matrix_pass && input_metadata_pass
) {
  FINAL_STATUS <- "PARTIAL PASS"
  NEXT_ACTION <- paste0( "The frozen sparse input is prepared, but one or more of the ", "mouse TF list, Python/pySCENIC CLI, or cisTarget resources is incomplete. ", "Repair the environment/resources before running GRNBoost2." )
} else {
  FINAL_STATUS <- "FAIL"
  NEXT_ACTION <- paste0( "start GRNBoost2. Resolve the failed input/metadata audit first." )
}
decision_table <- data.frame( domain = c( "Formal input files", "Frozen state counts", "Sparse raw-count matrix", "Required cell metadata", "Step6B module metadata", "Mouse TF filtering", "Python required-package imports", "pySCENIC CLI via system2", "Mouse cisTarget feather", "Mouse motif annotation", "Mouse TF list", "Overall Step 7A" ), pass = c( all(input_file_audit$exists), all(state_count_audit$exact_match), input_matrix_pass, critical_metadata_pass, module_metadata_pass, tf_filter_pass, python_import_pass, python_cli_pass, nrow(valid_feather) >= 1L, nrow(valid_annotation) >= 1L, nrow(valid_tf_audit) >= 1L && length(mouse_tfs) >= 100L, FINAL_STATUS == "PASS" ), status = c( ifelse( all(input_file_audit$exists), "PASS", "FAIL" ), ifelse( all(state_count_audit$exact_match), "PASS", "FAIL" ), ifelse( input_matrix_pass, "PASS", "FAIL" ), ifelse( critical_metadata_pass, "PASS", "FAIL" ), ifelse( module_metadata_pass, "PASS", "FAIL" ), ifelse( tf_filter_pass, "PASS", "PARTIAL/FAIL" ), ifelse( python_import_pass, "PASS", "FAIL" ), ifelse( python_cli_pass, "PASS", "FAIL" ), ifelse( nrow(valid_feather) >= 1L, "PASS", "FAIL" ), ifelse( nrow(valid_annotation) >= 1L, "PASS", "FAIL" ), ifelse( nrow(valid_tf_audit) >= 1L && length(mouse_tfs) >= 100L, "PASS", "FAIL" ), FINAL_STATUS ), stringsAsFactors = FALSE )
write_utf8_csv( decision_table, file.path( TABDIR, "13_final_decision_table.csv" ) )
selected_resource_lines <- c( paste0( "Selected Python: ", ifelse( is.na(PYTHON_SELECTED), "<none>", PYTHON_SELECTED ) ), paste0( "Selected mouse TF file: ", ifelse( is.na(TF_FILE_SELECTED), "<none>", safe_normalize_path( TF_FILE_SELECTED ) ) ), paste0( "Valid mouse feather files: ", nrow(valid_feather) ), paste0( "Valid mouse motif annotation files: ", nrow(valid_annotation) ), paste0( "Valid mouse TF-list files: ", nrow(valid_tf_audit) ) )
writeLines( capture.output( sessionInfo() ), con = OUT_SESSION, useBytes = TRUE )
rm( norm_frozen, counts_frozen, full_obj, seurat_meta, seurat_meta_frozen )
gc()
