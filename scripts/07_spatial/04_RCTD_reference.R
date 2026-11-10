options( stringsAsFactors = FALSE, warn = 1, width = 180, Seurat.object.assay.version = "v5", future.globals.maxSize = 16 * 1024^3 )
PROJECT_ROOT <- Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
SEED <- 20260804L
MAX_NUCLEI_PER_RCTD_TYPE <- 2000L
MIN_NUCLEI_PER_RCTD_TYPE <- 25L
MIN_COMMON_GENES <- 10000L
MIN_REFERENCE_OVERLAP_FRACTION <- 0.80
MIN_SPATIAL_OVERLAP_FRACTION <- 0.80
EXPECTED_ATLAS_NUCLEI <- 144359L
EXPECTED_MYELOID_NUCLEI <- 25205L
EXPECTED_SPATIAL_SPOTS <- 33941L
EXPECTED_CAPTURE_AREAS <- 12L
EXPECTED_SECTIONS <- 36L
OVERWRITE_OUTPUT <- TRUE
STEP10D_OBJECT <- file.path( PROJECT_ROOT, "04_results", "step10D_GSE234774_marker_first_major_lineage_annotation", "RDS", "GSE234774_step10D_major_lineage_annotation_evidence_BPCells.rds" )
STEP10G_OBJECT <- file.path( PROJECT_ROOT, "04_results", "step10G_GSE234774_frozen_myeloid_annotations_resident_objects", "RDS", "GSE234774_step10G_full_myeloid_final_annotations_BPCells.rds" )
STEP10G_REPORT <- file.path( PROJECT_ROOT, "04_results", "step10G_GSE234774_frozen_myeloid_annotations_resident_objects", "10_STEP10G_FINAL_REPORT.txt" )
STEP12B_DIR <- file.path( PROJECT_ROOT, "04_results", "step12B_fast_BPCells_spatial_section_QC" )
STEP12B_REPORT <- file.path( STEP12B_DIR, "10_STEP12B_FAST_FINAL_REPORT.txt" )
STEP12B_SPOT_QC <- file.path( STEP12B_DIR, "tables", "04_spot_QC_flags.csv.gz" )
STEP12B_CAPTURE_AUDIT <- file.path( STEP12B_DIR, "tables", "02_capture_area_image_coordinate_audit.csv" )
SPATIAL_BPCELLS_DIR <- file.path( STEP12B_DIR, "BPCells", "GSE234774_spatial_2d_full_counts" )
STEP12A2_ASSET_MANIFEST <- file.path( PROJECT_ROOT, "04_results", "step12A2_asset_manifest", "01_2d_capture_area_asset_manifest.csv" )
OUTPUT_DIR <- file.path( PROJECT_ROOT, "04_results", "step12C_GSE234774_snRNA_reference_RCTD_input_audit" )
TABLE_DIR <- file.path(OUTPUT_DIR, "tables")
RDS_DIR <- file.path(OUTPUT_DIR, "RDS")
RUN_LOG <- file.path(OUTPUT_DIR, "00_run_log.txt")
SESSION_INFO <- file.path(OUTPUT_DIR, "00_sessionInfo.txt")
OUT_INPUT_AUDIT <- file.path( TABLE_DIR, "01_input_and_field_audit.csv" )
OUT_LABEL_MAP <- file.path( TABLE_DIR, "02_RCTD_reference_label_mapping.csv" )
OUT_SAMPLING <- file.path( TABLE_DIR, "03_RCTD_reference_celltype_sampling_summary.csv" )
OUT_GENE_OVERLAP <- file.path( TABLE_DIR, "04_reference_spatial_gene_overlap.csv" )
OUT_QUERY_INDEX <- file.path( TABLE_DIR, "05_RCTD_query_capture_area_index.csv" )
OUT_CHECKS <- file.path( TABLE_DIR, "09_STEP12C_audit_checks.csv" )
OUT_REPORT <- file.path( OUTPUT_DIR, "10_STEP12C_FINAL_REPORT.txt" )
OUT_REFERENCE_RDS <- file.path( RDS_DIR, "GSE234774_step12C_RCTD_reference.rds" )
OUT_REFERENCE_METADATA_RDS <- file.path( RDS_DIR, "GSE234774_step12C_RCTD_reference_cell_metadata.rds" )
OUT_HANDOFF_RDS <- file.path( RDS_DIR, "GSE234774_step12C_RCTD_reference_query_handoff.rds" )
if (OVERWRITE_OUTPUT && dir.exists(OUTPUT_DIR)) {
  unlink(OUTPUT_DIR, recursive = TRUE, force = TRUE)
}
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(RDS_DIR, recursive = TRUE, showWarnings = FALSE)
required_packages <- c( "data.table", "Matrix", "Seurat", "SeuratObject", "BPCells", "spacexr" )
missing_packages <- required_packages[ !vapply( required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1) ) ]
if (length(missing_packages) > 0L) {
  stop( "Missing required package(s): ", paste(missing_packages, collapse = ", "), call. = FALSE )
}
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(Seurat)
  library(SeuratObject)
  library(BPCells)
  library(spacexr)
})
log_lines <- character()
write_log <- function() {
  writeLines(log_lines, RUN_LOG, useBytes = TRUE)
}
collapse_unique <- function(x, max_values = 30L) {
  x <- trimws(as.character(x))
  x <- unique(x[!is.na(x) & nzchar(x)])
  if (length(x) == 0L) {
    return(NA_character_)
  }
  if (length(x) > max_values) {
    x <- c( x[seq_len(max_values)], paste0("... +", length(x) - max_values, " more") )
  }
  paste(x, collapse = " | ")
}
required_columns <- function(x, columns, label) {
  missing <- setdiff(columns, names(x))
  if (length(missing) > 0L) {
    stop( label, " is missing required column(s): ", paste(missing, collapse = ", "), call. = FALSE )
  }
}
as_logical_strict <- function(x, label) {
  if (is.logical(x)) {
    return(x)
  }
  value <- toupper(trimws(as.character(x)))
  out <- rep(NA, length(value))
  out[value %in% c("TRUE", "T", "1")] <- TRUE
  out[value %in% c("FALSE", "F", "0")] <- FALSE
  bad <- is.na(out) & !is.na(value) & nzchar(value)
  if (any(bad)) {
    stop( label, " contains invalid logical values: ", paste(unique(value[bad]), collapse = " | "), call. = FALSE )
  }
  out
}
object_valid <- function(x) {
  tryCatch( {
      methods::validObject(x)
      TRUE
    },
    error = function(e) FALSE
  )
}
single_counts_layer <- function(object, label) {
  layers <- SeuratObject::Layers( object = object, assay = "RNA", search = NA )
  counts_layers <- layers[ grepl("^counts($|\\.)", layers) ]
  if (length(counts_layers) != 1L) {
    stop( label, " must have exactly one RNA counts layer; detected ", length(counts_layers), ": ", paste(counts_layers, collapse = " | "), call. = FALSE )
  }
  counts_layers[[1L]]
}
choose_field <- function(metadata, candidates, role) {
  present <- candidates[candidates %in% names(metadata)]
  if (length(present) == 0L) {
    stop( "No metadata field found for ", role, ". Candidates: ", paste(candidates, collapse = " | "), call. = FALSE )
  }
  informative <- present[ vapply( present, function(name) {
        data.table::uniqueN(metadata[[name]], na.rm = TRUE) > 1L
      },
      FUN.VALUE = logical(1)
    )
  ]
  if (length(informative) > 0L) {
    informative[[1L]]
  } else {
    present[[1L]]
  }
}
to_dgCMatrix <- function(x, label) {
  if (inherits(x, "dgCMatrix")) {
    return(x)
  }
  converted <- tryCatch( methods::as(x, "dgCMatrix"), error = function(e) e )
  if (inherits(converted, "error")) {
    intermediate <- tryCatch( methods::as(x, "CsparseMatrix"), error = function(e) e )
    if (!inherits(intermediate, "error")) {
      converted <- tryCatch( methods::as(intermediate, "dgCMatrix"), error = function(e) e )
    }
  }
  if (inherits(converted, "error")) {
    stop( "Could not convert ", label, " to dgCMatrix without dense conversion: ", conditionMessage(converted), call. = FALSE )
  }
  converted
}
balanced_sample <- function( metadata, current_type, library_column, max_cells, seed) {
  current <- metadata[ metadata$rctd_cell_type == current_type, , drop = FALSE ]
  if (nrow(current) <= max_cells) {
    return(rownames(current))
  }
  library_values <- as.character( current[[library_column]] )
  library_values[ is.na(library_values) | !nzchar(library_values) ] <- "MISSING_LIBRARY"
  libraries <- sort(unique(library_values))
  set.seed(seed)
  base_quota <- floor(max_cells / length(libraries))
  selected <- character()
  for (library_id in libraries) {
    candidates <- rownames(current)[ library_values == library_id ]
    take_n <- min(length(candidates), base_quota)
    if (take_n > 0L) {
      selected <- c( selected, sample(candidates, take_n, replace = FALSE) )
    }
  }
  remaining_n <- max_cells - length(selected)
  if (remaining_n > 0L) {
    remaining <- setdiff(rownames(current), selected)
    selected <- c( selected, sample(remaining, remaining_n, replace = FALSE) )
  }
  selected
}
file_size_safe <- function(path) {
  if (!file.exists(path)) {
    return(NA_real_)
  }
  as.numeric(file.info(path)$size)
}
checks <- data.table( check_id = character(), category = character(), description = character(), status = character(), observed = character(), required = character(), critical = logical(), action_if_not_pass = character() )
add_check <- function( check_id, category, description, passed, observed, required, critical = TRUE, action_if_not_pass = "") {
  checks <<- rbind( checks, data.table( check_id = as.character(check_id), category = as.character(category), description = as.character(description), status = if (isTRUE(passed)) "PASS" else "FAIL", observed = as.character(observed), required = as.character(required), critical = isTRUE(critical), action_if_not_pass = as.character(action_if_not_pass) ), fill = TRUE )
}
add_review <- function( check_id, category, description, observed, required, action_if_not_pass) {
  checks <<- rbind( checks, data.table( check_id = as.character(check_id), category = as.character(category), description = as.character(description), status = "REVIEW", observed = as.character(observed), required = as.character(required), critical = FALSE, action_if_not_pass = as.character(action_if_not_pass) ), fill = TRUE )
}
write_fatal_report <- function(reason) {
  write_log()
}
main <- function() {
  start_time <- Sys.time()
  required_inputs <- c( step10D_object = STEP10D_OBJECT, step10G_object = STEP10G_OBJECT, step10G_report = STEP10G_REPORT, step12B_report = STEP12B_REPORT, step12B_spot_QC = STEP12B_SPOT_QC, step12B_capture_audit = STEP12B_CAPTURE_AUDIT, step12A2_asset_manifest = STEP12A2_ASSET_MANIFEST )
  input_exists <- file.exists(required_inputs)
  spatial_cache_exists <- dir.exists(SPATIAL_BPCELLS_DIR)
  add_check( "C01", "input", "All required frozen inputs and the spatial BPCells cache exist.", all(input_exists) && spatial_cache_exists, paste0( paste( paste0(names(required_inputs), "=", input_exists), collapse = "; " ), "; spatial_BPCells=", spatial_cache_exists ), "All inputs TRUE", TRUE, "Restore frozen outputs; silently substitute files." )
  if (!all(input_exists) || !spatial_cache_exists) {
    stop("At least one required Step 12C input is missing.", call. = FALSE)
  }
  step12B_text <- readLines( STEP12B_REPORT, warn = FALSE, encoding = "UTF-8" )
  step12B_status <- step12B_text[ grepl("^FINAL_STATUS:", step12B_text) ]
  step12B_pass <- length(step12B_status) == 1L && grepl( "^FINAL_STATUS: (PASS|PASS_WITH_REVIEW)$", step12B_status )
  add_check( "C02", "upstream", "Step 12B-FAST releases the spatial inputs.", step12B_pass, if (length(step12B_status) > 0L) {
      step12B_status[[1L]]
    } else {
      "No status line"
    },
    "PASS or PASS_WITH_REVIEW", TRUE, "prepare RCTD inputs until Step 12B passes."
  )
  if (!step12B_pass) {
    stop("Step 12B-FAST does not release Step 12C.", call. = FALSE)
  }
  step10G_text <- readLines( STEP10G_REPORT, warn = FALSE, encoding = "UTF-8" )
  step10G_status <- step10G_text[ grepl("^FINAL DECISION:", step10G_text) ]
  step10G_pass <- length(step10G_status) == 1L && grepl("PASS", step10G_status)
  add_check( "C03", "upstream", "Step 10G releases the refined myeloid annotations.", step10G_pass, if (length(step10G_status) > 0L) {
      step10G_status[[1L]]
    } else {
      "No decision line"
    },
    "A Step 10G PASS decision", TRUE, "backfill provisional myeloid labels."
  )
  if (!step10G_pass) {
    stop("Step 10G does not release its annotations.", call. = FALSE)
  }
  spot_qc <- fread( STEP12B_SPOT_QC, data.table = TRUE, showProgress = FALSE )
  capture_audit <- fread( STEP12B_CAPTURE_AUDIT, data.table = TRUE, showProgress = FALSE )
  asset_manifest <- fread( STEP12A2_ASSET_MANIFEST, data.table = TRUE, showProgress = FALSE )
  required_columns( spot_qc, c( "barcode", "capture_area_id", "section_id", "timepoint", "nCount_Spatial_BPCells", "qc_review_flag" ), "Step 12B spot-QC table" )
  required_columns( capture_audit, c( "capture_area_id", "processed_spots", "section_count", "asset_linkage_status", "processed_position_overlap_fraction", "processed_in_tissue_overlap_fraction", "coordinate_bounds_ok" ), "Step 12B capture audit" )
  required_columns( asset_manifest, c( "capture_area_id", "position_member", "image_lowres_member", "scalefactor_member" ), "Step 12A2 asset manifest" )
  spot_qc[, qc_review_flag := as_logical_strict( qc_review_flag, "qc_review_flag" )]
  all_spot_rule_ok <- nrow(spot_qc) == EXPECTED_SPATIAL_SPOTS && !anyDuplicated(spot_qc$barcode) && uniqueN(spot_qc$capture_area_id) == EXPECTED_CAPTURE_AREAS && uniqueN(spot_qc$section_id) == EXPECTED_SECTIONS
  add_check( "C04", "spatial_rule", "The frozen primary spatial query retains all 33,941 spots.", all_spot_rule_ok, paste0( "spots=", nrow(spot_qc), "; capture areas=", uniqueN(spot_qc$capture_area_id), "; sections=", uniqueN(spot_qc$section_id), "; retained review spots=", sum(spot_qc$qc_review_flag) ), "33,941 unique spots, 12 capture areas, 36 sections", TRUE, "silently remove the five conservative review spots." )
  if (!all_spot_rule_ok) {
    stop("The frozen all-spot rule was not reproduced.", call. = FALSE)
  }
  atlas <- readRDS(STEP10D_OBJECT)
  myeloid <- readRDS(STEP10G_OBJECT)
  object_check <- object_valid(atlas) && object_valid(myeloid) && ncol(atlas) == EXPECTED_ATLAS_NUCLEI && ncol(myeloid) == EXPECTED_MYELOID_NUCLEI
  add_check( "C05", "snRNA_object", "Step 10D and Step 10G objects reload with expected nucleus counts.", object_check, paste0( "atlas=", nrow(atlas), " x ", ncol(atlas), "; myeloid=", nrow(myeloid), " x ", ncol(myeloid) ), paste0( "atlas nuclei=", EXPECTED_ATLAS_NUCLEI, "; myeloid nuclei=", EXPECTED_MYELOID_NUCLEI ), TRUE, "rebuild or reinterpret mismatched frozen objects." )
  if (!object_check) {
    stop("Frozen snRNA object check failed.", call. = FALSE)
  }
  required_columns( atlas@meta.data, c( "major_lineage_candidate", "major_lineage_candidate_confidence", "hcon_res_0_2" ), "Step 10D metadata" )
  required_columns( myeloid@meta.data, c( "myeloid_final_label", "resident_lineage_status", "step10G_annotation_status" ), "Step 10G metadata" )
  library_column <- choose_field( atlas@meta.data, c( "biological_library", "library_id", "sample_id", "orig.ident" ), "biological library" )
  timepoint_column <- choose_field( atlas@meta.data, c( "timepoint", "time", "group" ), "timepoint" )
  counts_layer <- single_counts_layer( atlas, "Step 10D atlas" )
  atlas_counts <- SeuratObject::LayerData( object = atlas, assay = "RNA", layer = counts_layer, fast = FALSE )
  counts_aligned <- identical(colnames(atlas_counts), colnames(atlas)) && identical(rownames(atlas_counts), rownames(atlas))
  add_check( "C06", "snRNA_counts", "Step 10D counts are aligned to the atlas object.", counts_aligned, paste0( "layer=", counts_layer, "; dimensions=", nrow(atlas_counts), " x ", ncol(atlas_counts), "; class=", paste(class(atlas_counts), collapse = " | ") ), "Exact feature and nucleus order", TRUE, "materialize a misaligned reference." )
  if (!counts_aligned) {
    stop("Step 10D counts are not aligned.", call. = FALSE)
  }
  atlas_meta <- atlas@meta.data
  atlas_meta$atlas_nucleus <- rownames(atlas_meta)
  major_map <- c( "Oligodendrocyte" = "Oligodendrocyte", "Neuron" = "Neuron", "Resident_microglia" = "Microglia", "Fibroblast" = "Fibroblast_VLMC", "Astrocyte" = "Astrocyte", "OPC" = "OPC", "Endothelial" = "Endothelial", "Ependymal" = "Ependymal", "NK_cell" = "NK_cell", "Pericyte" = "Pericyte" )
  observed_major <- sort(unique( as.character(atlas_meta$major_lineage_candidate) ))
  unknown_major <- setdiff(observed_major, names(major_map))
  if (length(unknown_major) > 0L) {
    stop( "Unmapped Step 10D major-lineage label(s): ", paste(unknown_major, collapse = " | "), call. = FALSE )
  }
  atlas_meta$rctd_cell_type <- unname( major_map[ as.character(atlas_meta$major_lineage_candidate) ] )
  atlas_meta$rctd_label_source <- "Step10D_major_lineage_candidate"
  myeloid_meta <- myeloid@meta.data
  myeloid_meta$myeloid_nucleus <- rownames(myeloid_meta)
  overlap_index <- match( myeloid_meta$myeloid_nucleus, atlas_meta$atlas_nucleus )
  exact_backfill <- !anyNA(overlap_index) && length(unique(overlap_index)) == nrow(myeloid_meta)
  add_check( "C07", "myeloid_backfill", "Every Step 10G myeloid nucleus maps uniquely to Step 10D.", exact_backfill, paste0( "mapped=", sum(!is.na(overlap_index)), "/", nrow(myeloid_meta), "; unique indices=", length(unique(overlap_index[!is.na(overlap_index)])) ), paste0( EXPECTED_MYELOID_NUCLEI, "/", EXPECTED_MYELOID_NUCLEI, " exact matches" ), TRUE, "Transfer labels by exact nucleus barcode only." )
  if (!exact_backfill) {
    stop("Step 10G-to-Step 10D barcode mapping failed.", call. = FALSE)
  }
  myeloid_labels <- as.character( myeloid_meta$myeloid_final_label )
  resident_status <- as.character( myeloid_meta$resident_lineage_status )
  refined_label <- rep(NA_character_, nrow(myeloid_meta))
  refined_label[ resident_status == "RESIDENT_LINEAGE_INCLUDED" ] <- "Microglia"
  macrophage_like <- grepl( "BAM_like|macrophage_like_non_microglial_myeloid", myeloid_labels, ignore.case = TRUE )
  refined_label[macrophage_like] <- "BAM_macrophage"
  contaminant <- resident_status == "CONTAMINANT_EXCLUDED"
  unresolved <- is.na(refined_label) & !contaminant
  if (any(unresolved)) {
    stop( "Unresolved non-contaminant Step 10G labels: ", paste(unique(myeloid_labels[unresolved]), collapse = " | "), call. = FALSE )
  }
  atlas_meta$rctd_cell_type[overlap_index] <- refined_label
  atlas_meta$rctd_label_source[overlap_index] <- ifelse( contaminant, "Step10G_contaminant_excluded", "Step10G_refined_myeloid_backfill" )
  atlas_meta$rctd_reference_eligible <- !is.na(atlas_meta$rctd_cell_type)
  remaining_step10D_microglia <- atlas_meta[ atlas_meta$rctd_reference_eligible, "major_lineage_candidate" ] == "Resident_microglia" & atlas_meta[ atlas_meta$rctd_reference_eligible, "rctd_label_source" ] == "Step10D_major_lineage_candidate"
  forbidden_labels <- grepl( "target|arm[_ -]?[12]|dam|inflammatory.remodel", atlas_meta$rctd_cell_type, ignore.case = TRUE )
  broad_label_boundary <- !any(remaining_step10D_microglia) && !any(forbidden_labels, na.rm = TRUE)
  add_check( "C08", "reference_label", "All microglia are Step 10G-refined and the target program is not a cell type.", broad_label_boundary, paste0( "unbackfilled candidate microglia=", sum(remaining_step10D_microglia), "; forbidden target/Arm/DAM labels=", sum(forbidden_labels, na.rm = TRUE) ), "Zero for both", TRUE, "Keep target program as a later score, not a deconvolution lineage." )
  if (!broad_label_boundary) {
    stop("Broad reference-label boundary failed.", call. = FALSE)
  }
  atlas_map_table <- as.data.table( atlas_meta[ atlas_meta$rctd_label_source == "Step10D_major_lineage_candidate", , drop = FALSE ] )[ , .( source_layer = "Step10D_major_lineage_candidate", source_label = as.character( major_lineage_candidate[1L] ), rctd_cell_type = as.character( rctd_cell_type[1L] ), n_nuclei = .N, annotation_confidence = collapse_unique( major_lineage_candidate_confidence ), reference_action = "INCLUDE" ), by = major_lineage_candidate ]
  atlas_map_table[, major_lineage_candidate := NULL]
  myeloid_map_data <- data.table( source_label = myeloid_labels, refined_rctd_label = refined_label )
  myeloid_map_table <- myeloid_map_data[ , .( source_layer = "Step10G_myeloid_final_label", rctd_cell_type = collapse_unique( refined_rctd_label ), n_nuclei = .N, annotation_confidence = "FROZEN_STEP10G", reference_action = if ( all(is.na(refined_rctd_label)) ) {
        "EXCLUDE_CONTAMINANT"
      } else {
        "INCLUDE_BACKFILL"
      }
    ), by = source_label
  ]
  label_mapping <- rbindlist( list(atlas_map_table, myeloid_map_table), fill = TRUE )
  setorder(label_mapping, source_layer, source_label)
  ambiguous_n <- sum( grepl( "AMBIG|LOW", as.character( atlas_meta$major_lineage_candidate_confidence ), ignore.case = TRUE ) )
  if (ambiguous_n > 0L) {
    add_review( "R01", "annotation", "Some Step 10D nuclei have ambiguous/low-confidence broad-lineage evidence.", paste0("nuclei=", ambiguous_n), "Use only broad lineage interpretation.", "turn these labels into new fine biological claims." )
  }
  eligible_meta <- atlas_meta[ atlas_meta$rctd_reference_eligible, , drop = FALSE ]
  rownames(eligible_meta) <- eligible_meta$atlas_nucleus
  cell_types <- sort(unique(eligible_meta$rctd_cell_type))
  eligible_counts <- table(eligible_meta$rctd_cell_type)
  small_types <- names( eligible_counts[ eligible_counts < MIN_NUCLEI_PER_RCTD_TYPE ] )
  add_check( "C09", "reference_sampling", "Every broad RCTD cell type has at least 25 eligible nuclei.", length(small_types) == 0L, if (length(small_types) == 0L) {
      paste0( "types=", length(eligible_counts), "; minimum=", min(eligible_counts) )
    } else {
      paste( paste0( small_types, "=", eligible_counts[small_types] ), collapse = " | " )
    },
    paste0( "At least ", MIN_NUCLEI_PER_RCTD_TYPE, " nuclei/type" ), TRUE, "run RCTD with an unsupported rare type."
  )
  if (length(small_types) > 0L) {
    stop("A reference type has fewer than 25 nuclei.", call. = FALSE)
  }
  sampled_by_type <- vector("list", length(cell_types))
  names(sampled_by_type) <- cell_types
  for (i in seq_along(cell_types)) {
    sampled_by_type[[i]] <- balanced_sample( metadata = eligible_meta, current_type = cell_types[[i]], library_column = library_column, max_cells = MAX_NUCLEI_PER_RCTD_TYPE, seed = SEED + i )
  }
  sampled_cells <- unlist(sampled_by_type, use.names = FALSE)
  if (anyDuplicated(sampled_cells)) {
    stop("Reference sampling produced duplicate nuclei.", call. = FALSE)
  }
  sampled_meta <- eligible_meta[ sampled_cells, , drop = FALSE ]
  sampling_summary <- rbindlist( lapply( cell_types, function(cell_type) {
        full_rows <- eligible_meta[ eligible_meta$rctd_cell_type == cell_type, , drop = FALSE ]
        sample_rows <- sampled_meta[ sampled_meta$rctd_cell_type == cell_type, , drop = FALSE ]
        data.table( rctd_cell_type = cell_type, eligible_nuclei = nrow(full_rows), sampled_nuclei = nrow(sample_rows), sampling_fraction = nrow(sample_rows) / nrow(full_rows), eligible_libraries = uniqueN( full_rows[[library_column]], na.rm = TRUE ), sampled_libraries = uniqueN( sample_rows[[library_column]], na.rm = TRUE ), eligible_timepoints = uniqueN( full_rows[[timepoint_column]], na.rm = TRUE ), sampled_timepoints = uniqueN( sample_rows[[timepoint_column]], na.rm = TRUE ), sampled_library_values = collapse_unique( sample_rows[[library_column]] ), sampled_timepoint_values = collapse_unique( sample_rows[[timepoint_column]] ), minimum_cell_rule_pass = nrow(sample_rows) >= MIN_NUCLEI_PER_RCTD_TYPE, maximum_cell_cap_pass = nrow(sample_rows) <= MAX_NUCLEI_PER_RCTD_TYPE )
      }
    )
  )
  sampling_ok <- all(sampling_summary$minimum_cell_rule_pass) && all(sampling_summary$maximum_cell_cap_pass) && sum(sampling_summary$sampled_nuclei) == length(sampled_cells)
  add_check( "C10", "reference_sampling", "The sampled reference obeys the minimum and per-type cap.", sampling_ok, paste0( "types=", nrow(sampling_summary), "; sampled nuclei=", length(sampled_cells), "; min/type=", min(sampling_summary$sampled_nuclei), "; max/type=", max(sampling_summary$sampled_nuclei) ), paste0( MIN_NUCLEI_PER_RCTD_TYPE, "-", MAX_NUCLEI_PER_RCTD_TYPE, " nuclei/type" ), TRUE, "build an undersized or unbounded reference." )
  if (!sampling_ok) {
    stop("Reference sampling check failed.", call. = FALSE)
  }
  low_library_types <- sampling_summary[ sampled_libraries < 2L, rctd_cell_type ]
  if (length(low_library_types) > 0L) {
    add_review( "R02", "reference_sampling", "Some sampled broad types use fewer than two biological libraries.", paste(low_library_types, collapse = " | "), "Report limited reference diversity.", "invent replicate labels." )
  }
  spatial_counts <- BPCells::open_matrix_dir( SPATIAL_BPCELLS_DIR )
  if ( ncol(spatial_counts) != EXPECTED_SPATIAL_SPOTS || nrow(spatial_counts) <= 0L ) {
    stop("Spatial BPCells dimensions are unexpected.", call. = FALSE)
  }
  reference_genes <- rownames(atlas_counts)
  spatial_genes <- rownames(spatial_counts)
  if ( anyDuplicated(reference_genes) || anyDuplicated(spatial_genes) ) {
    stop( "Reference or spatial BPCells feature names are duplicated.", call. = FALSE )
  }
  common_genes <- intersect(reference_genes, spatial_genes)
  overlap_reference <- length(common_genes) / length(reference_genes)
  overlap_spatial <- length(common_genes) / length(spatial_genes)
  overlap_ok <- length(common_genes) >= MIN_COMMON_GENES && overlap_reference >= MIN_REFERENCE_OVERLAP_FRACTION && overlap_spatial >= MIN_SPATIAL_OVERLAP_FRACTION
  add_check( "C11", "gene_overlap", "The snRNA and spatial matrices have sufficient exact shared genes.", overlap_ok, paste0( "common=", length(common_genes), "; reference fraction=", round(overlap_reference, 4), "; spatial fraction=", round(overlap_spatial, 4) ), paste0( ">=", MIN_COMMON_GENES, " genes and both fractions >=0.80" ), TRUE, "run RCTD across mismatched identifiers." )
  if (!overlap_ok) {
    stop("Reference-spatial gene overlap is insufficient.", call. = FALSE)
  }
  gene_overlap <- data.table( metric = c( "reference_genes", "spatial_genes", "common_genes", "common_fraction_of_reference", "common_fraction_of_spatial", "reference_only_genes", "spatial_only_genes", "reference_only_examples", "spatial_only_examples" ), value = c( length(reference_genes), length(spatial_genes), length(common_genes), overlap_reference, overlap_spatial, length(setdiff(reference_genes, spatial_genes)), length(setdiff(spatial_genes, reference_genes)), collapse_unique( head(setdiff(reference_genes, spatial_genes), 30L) ), collapse_unique( head(setdiff(spatial_genes, reference_genes), 30L) ) ), interpretation = c( "Step 10D RNA count features", "Step 12B spatial BPCells features", "Exact overlap passed to RCTD", "Shared fraction of reference", "Shared fraction of spatial", "Not used in Reference", "Unavailable in reference", "Examples only", "Examples only" ) )
  gene_index <- match(common_genes, reference_genes)
  cell_index <- match(sampled_cells, colnames(atlas_counts))
  if (anyNA(gene_index) || anyNA(cell_index)) {
    stop( "Common genes or sampled nuclei cannot be indexed.", call. = FALSE )
  }
  reference_view <- atlas_counts[ gene_index, cell_index, drop = FALSE ]
  reference_counts <- to_dgCMatrix( reference_view, "sampled RCTD reference" )
  rm(reference_view)
  invisible(gc())
  if ( !identical(rownames(reference_counts), common_genes) || !identical(colnames(reference_counts), sampled_cells) ) {
    stop("Materialized reference names are misaligned.", call. = FALSE)
  }
  expressed_genes <- Matrix::rowSums(reference_counts) > 0
  reference_counts <- reference_counts[ expressed_genes, , drop = FALSE ]
  reference_nUMI <- as.numeric( Matrix::colSums(reference_counts) )
  names(reference_nUMI) <- colnames(reference_counts)
  if (any(!is.finite(reference_nUMI) | reference_nUMI <= 0)) {
    stop( "The sampled reference contains zero/non-finite nUMI.", call. = FALSE )
  }
  reference_labels <- factor( sampled_meta[ colnames(reference_counts), "rctd_cell_type" ], levels = sort(unique(sampled_meta$rctd_cell_type)) )
  names(reference_labels) <- colnames(reference_counts)
  reference_input_ok <- nrow(reference_counts) >= MIN_COMMON_GENES && ncol(reference_counts) == length(reference_labels) && ncol(reference_counts) == length(reference_nUMI) && identical(names(reference_labels), colnames(reference_counts)) && identical(names(reference_nUMI), colnames(reference_counts))
  add_check( "C12", "reference_matrix", "Reference counts, cell-type labels, and nUMI are exactly aligned.", reference_input_ok, paste0( "genes=", nrow(reference_counts), "; nuclei=", ncol(reference_counts), "; cell types=", nlevels(reference_labels), "; nonzero=", length(reference_counts@x) ), "Exact column alignment and >=10,000 genes", TRUE, "construct a misaligned spacexr Reference." )
  if (!reference_input_ok) {
    stop("Reference input alignment failed.", call. = FALSE)
  }
  Reference_constructor <- get( "Reference", envir = asNamespace("spacexr") )
  rctd_reference <- Reference_constructor( reference_counts, reference_labels, reference_nUMI )
  sampled_meta_save <- sampled_meta[ colnames(reference_counts), , drop = FALSE ]
  sampled_meta_save$reference_nUMI <- reference_nUMI[ rownames(sampled_meta_save) ]
  saveRDS( rctd_reference, OUT_REFERENCE_RDS, compress = FALSE )
  saveRDS( sampled_meta_save, OUT_REFERENCE_METADATA_RDS, compress = FALSE )
  reference_saved <- file.exists(OUT_REFERENCE_RDS) && file_size_safe(OUT_REFERENCE_RDS) > 0 && file.exists(OUT_REFERENCE_METADATA_RDS) && file_size_safe(OUT_REFERENCE_METADATA_RDS) > 0
  add_check( "C13", "reference_output", "The spacexr Reference and sampled metadata were saved.", reference_saved, paste0( "reference bytes=", file_size_safe(OUT_REFERENCE_RDS), "; metadata bytes=", file_size_safe(OUT_REFERENCE_METADATA_RDS) ), "Both nonempty", TRUE, "proceed to Step 12D without a saved Reference." )
  if (!reference_saved) {
    stop("Reference outputs could not be saved.", call. = FALSE)
  }
  capture_audit[ , capture_area_id := toupper( trimws(as.character(capture_area_id)) ) ]
  capture_audit[ , asset_linkage_status := toupper( trimws(as.character(asset_linkage_status)) ) ]
  capture_audit[ , coordinate_bounds_ok := as_logical_strict( coordinate_bounds_ok, "coordinate_bounds_ok" ) ]
  asset_manifest[ , capture_area_id := toupper( trimws(as.character(capture_area_id)) ) ]
  query_index <- spot_qc[ , .( n_spots = .N, n_sections = uniqueN(section_id), section_ids = collapse_unique(section_id), timepoints = collapse_unique(timepoint), n_QC_review_spots_retained = sum(qc_review_flag), total_nUMI = sum( as.numeric(nCount_Spatial_BPCells) ), min_spot_nUMI = min( as.numeric(nCount_Spatial_BPCells) ), median_spot_nUMI = median( as.numeric(nCount_Spatial_BPCells) ) ), by = capture_area_id ]
  query_index <- merge( query_index, capture_audit[ , .( capture_area_id, image_coordinate_linkage = asset_linkage_status, processed_position_overlap_fraction, processed_in_tissue_overlap_fraction, coordinate_bounds_ok ) ], by = "capture_area_id", all.x = TRUE, sort = TRUE )
  query_index <- merge( query_index, asset_manifest[ , .( capture_area_id, position_member, image_lowres_member, scalefactor_member ) ], by = "capture_area_id", all.x = TRUE, sort = TRUE )
  query_index[, spatial_counts_BPCells_path := SPATIAL_BPCELLS_DIR ]
  query_index[, common_reference_genes := nrow(reference_counts) ]
  query_index[, spot_retention_rule := "PRIMARY_RETAIN_ALL_33941_SPOTS" ]
  query_index[, planned_RCTD_mode := "full"]
  query_index[, Step12D_action := paste( "Create one SpatialRNA query per capture area;", "then create.RCTD and run.RCTD(full)." ) ]
  query_ready <- nrow(query_index) == EXPECTED_CAPTURE_AREAS && sum(query_index$n_spots) == EXPECTED_SPATIAL_SPOTS && all(query_index$n_sections == 3L) && all(query_index$image_coordinate_linkage == "PASS") && all( query_index$processed_position_overlap_fraction == 1 ) && all( query_index$processed_in_tissue_overlap_fraction == 1 ) && all(query_index$coordinate_bounds_ok)
  add_check( "C14", "query_audit", "All 12 all-spot queries are ready for later full-mode RCTD.", query_ready, paste0( "capture areas=", nrow(query_index), "; spots=", sum(query_index$n_spots), "; sections=", sum(query_index$n_sections), "; linkage PASS=", sum(query_index$image_coordinate_linkage == "PASS") ), "12 areas, 33,941 spots, 36 sections, complete linkage", TRUE, "create SpatialRNA queries until all areas are complete." )
  if (!query_ready) {
    stop("One or more future RCTD queries are not ready.", call. = FALSE)
  }
  handoff <- list( step = "Step12C", completed_at = as.character(Sys.time()), frozen_spot_rule = "retain_all_33941_spots_no_sensitivity_branch", reference = list( rds_path = OUT_REFERENCE_RDS, metadata_rds_path = OUT_REFERENCE_METADATA_RDS, n_genes = nrow(reference_counts), n_nuclei = ncol(reference_counts), cell_types = levels(reference_labels), max_nuclei_per_type = MAX_NUCLEI_PER_RCTD_TYPE, target_program_is_cell_type = FALSE ), spatial = list( BPCells_path = SPATIAL_BPCELLS_DIR, n_spots = EXPECTED_SPATIAL_SPOTS, n_capture_areas = EXPECTED_CAPTURE_AREAS, n_sections = EXPECTED_SECTIONS, query_index = query_index, planned_mode = "full" ), next_step = paste( "Step 12D creates one SpatialRNA query per capture area", "and runs full-mode RCTD. Target scoring remains deferred." ) )
  saveRDS(handoff, OUT_HANDOFF_RDS, compress = FALSE)
  input_audit <- rbindlist( list( data.table( item_type = "package", item = required_packages, value = vapply( required_packages, function(package_name) {
            as.character( utils::packageVersion(package_name) )
          },
          FUN.VALUE = character(1)
        ), note = "Version used in Step 12C"
      ), data.table(
        item_type = "input", item = names(required_inputs), value = unname(required_inputs), note = paste0( "exists=", file.exists(unname(required_inputs)) )
      ), data.table(
        item_type = "field_or_rule", item = c( "atlas_dimensions", "myeloid_dimensions", "RNA_counts_layer", "biological_library_field", "timepoint_field", "eligible_reference_nuclei", "sampled_reference_nuclei", "spatial_spots_retained", "QC_review_spots_retained", "target_program_as_cell_type" ), value = c( paste(dim(atlas), collapse = " x "), paste(dim(myeloid), collapse = " x "), counts_layer, library_column, timepoint_column, nrow(eligible_meta), length(sampled_cells), nrow(spot_qc), sum(spot_qc$qc_review_flag), "FALSE" ), note = c( "features x nuclei", "features x nuclei", "raw counts", "balanced sampling", "coverage audit", "contaminants excluded", "max 2,000/type", "all retained", "retained", "program remains a score/state" )
      )
    ), fill = TRUE
  )
  fwrite(input_audit, OUT_INPUT_AUDIT, bom = TRUE)
  fwrite(label_mapping, OUT_LABEL_MAP, bom = TRUE)
  fwrite(sampling_summary, OUT_SAMPLING, bom = TRUE)
  fwrite(gene_overlap, OUT_GENE_OVERLAP, bom = TRUE)
  fwrite(query_index, OUT_QUERY_INDEX, bom = TRUE)
  critical_failures <- checks[ critical == TRUE & status == "FAIL" ]
  reviews <- checks[status == "REVIEW"]
  final_status <- if (nrow(critical_failures) > 0L) {
    "STOP"
  } else if (nrow(reviews) > 0L) {
    "PASS_WITH_REVIEW"
  } else {
    "PASS"
  }
  fwrite(checks, OUT_CHECKS, bom = TRUE)
  elapsed <- round( as.numeric( difftime( Sys.time(), start_time, units = "mins" ) ), 3 )
  if (nrow(critical_failures) > 0L) {
  }
  if (nrow(reviews) > 0L) {
  }
  writeLines( capture.output(sessionInfo()), SESSION_INFO, useBytes = TRUE )
  write_log()
  invisible( list( final_status = final_status, reference_path = OUT_REFERENCE_RDS, metadata_path = OUT_REFERENCE_METADATA_RDS, handoff_path = OUT_HANDOFF_RDS, sampling_summary = sampling_summary, query_index = query_index, checks = checks ) )
}
tryCatch( main(), error = function(e) {
    write_fatal_report(conditionMessage(e))
  }
)
