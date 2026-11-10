options(stringsAsFactors = FALSE, warn = 1)
if (!requireNamespace("data.table", quietly = TRUE)) {
  stop( "Required package 'data.table' is not installed. Run ", "install.packages('data.table') once, then source this script again." )
}
suppressPackageStartupMessages(library(data.table))
PROJECT_ROOT <- Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
INPUT_DIR <- file.path(PROJECT_ROOT, "01_raw_data", "GSE234774")
OUTPUT_DIR <- file.path( PROJECT_ROOT, "04_results", "step12A2_GSE234774_spatial_metadata_asset_linkage_repair" )
EXTRACT_DIR <- file.path(OUTPUT_DIR, "audit_extracted_2d_assets")
MATRIX_FILE <- file.path( INPUT_DIR, "GSE234774_spatial_2d_filtered_spatial_2d.mtx.gz" )
BARCODE_FILE <- file.path( INPUT_DIR, "GSE234774_spatial_2d_barcodes.txt.gz" )
FEATURE_FILE <- file.path( INPUT_DIR, "GSE234774_spatial_2d_features.txt.gz" )
METADATA_FILE <- file.path( INPUT_DIR, "GSE234774_spatial_2d_meta.txt.gz" )
RAW_TAR <- file.path(INPUT_DIR, "GSE234774_RAW.tar")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
RUN_LOG_FILE <- file.path(OUTPUT_DIR, "00_run_log.txt")
REPORT_FILE <- file.path(OUTPUT_DIR, "10_STEP12A2_FINAL_REPORT.txt")
log_lines <- character()
log_message <- function(...) {
  text <- paste0(...)
  cat(text, "\n")
  log_lines <<- c(log_lines, text)
}
write_run_log <- function() {
  writeLines(log_lines, RUN_LOG_FILE, useBytes = TRUE)
}
collapse_unique <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x) == 0L) {
    return(NA_character_)
  }
  paste(sort(unique(x)), collapse = " | ")
}
n_unique_nonmissing <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  uniqueN(x[!is.na(x) & nzchar(x)])
}
safe_character <- function(x) {
  out <- as.character(x)
  out[is.na(out)] <- ""
  trimws(out)
}
regex_escape <- function(x) {
  gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", x)
}
extract_capture_area <- function(x) {
  x <- basename(x)
  match_start <- regexpr("[0-9]{3}_[A-Za-z]", x, perl = TRUE)
  out <- rep(NA_character_, length(x))
  has_match <- match_start > 0L
  out[has_match] <- regmatches(x, match_start)[has_match]
  toupper(out)
}
extract_gsm <- function(x) {
  x <- basename(x)
  out <- ifelse( grepl("^GSM[0-9]+", x, perl = TRUE), sub("^(GSM[0-9]+).*", "\\1", x, perl = TRUE), NA_character_ )
  out
}
classify_tar_asset <- function(x) {
  b <- tolower(basename(x))
  out <- rep("other", length(b))
  out[grepl("_tissue_hires_image\\.png\\.gz$", b, perl = TRUE)] <- "image_hires"
  out[grepl("_tissue_lowres_image\\.png\\.gz$", b, perl = TRUE)] <- "image_lowres"
  out[grepl("_tissue_positions(_list)?\\.csv\\.gz$", b, perl = TRUE)] <- "position"
  out[grepl("_scalefactors_json\\.json\\.gz$", b, perl = TRUE)] <- "scalefactor"
  out
}
is_unsafe_tar_member <- function(x) {
  grepl("^(/|\\\\)", x) | grepl("^[A-Za-z]:[/\\\\]", x) | grepl("(^|[/\\\\])\\.\\.([/\\\\]|$)", x)
}
read_matrix_market_dimensions <- function(path) {
  connection <- gzfile(path, open = "rt")
  on.exit(close(connection), add = TRUE)
  first_line <- readLines(connection, n = 1L, warn = FALSE)
  if (length(first_line) != 1L || !grepl("^%%MatrixMarket", first_line)) {
    stop("The Matrix Market header is missing or unreadable: ", path)
  }
  dimension_line <- character()
  repeat {
    current <- readLines(connection, n = 1L, warn = FALSE)
    if (length(current) == 0L) {
      break
    }
    if (!grepl("^%", current)) {
      dimension_line <- current
      break
    }
  }
  if (length(dimension_line) != 1L) {
    stop("The Matrix Market dimension line could not be found: ", path)
  }
  values <- scan( text = dimension_line, what = numeric(), quiet = TRUE )
  if (length(values) < 3L) {
    stop("The Matrix Market dimension line is malformed: ", dimension_line)
  }
  data.table( n_features = as.integer(values[1]), n_barcodes = as.integer(values[2]), n_nonzero = as.numeric(values[3]) )
}
read_single_column <- function(path, output_name) {
  result <- fread( path, header = FALSE, select = 1L, col.names = output_name, data.table = TRUE, showProgress = FALSE )
  result[, (output_name) := trimws(as.character(get(output_name)))]
  result
}
read_spatial_metadata <- function(path) {
  result <- fread( path, header = TRUE, data.table = TRUE, showProgress = FALSE )
  setnames(result, trimws(names(result)))
  result
}
read_position_file <- function(path) {
  result <- fread( path, header = FALSE, data.table = TRUE, showProgress = FALSE )
  if (ncol(result) < 6L) {
    stop( "Position file has fewer than six columns: ", basename(path), " (observed ", ncol(result), ")" )
  }
  first_value <- tolower(trimws(as.character(result[[1L]][1L])))
  if (first_value %in% c("barcode", "barcodes")) {
    result <- fread( path, header = TRUE, data.table = TRUE, showProgress = FALSE )
    normalized <- tolower(gsub("[^a-z0-9]+", "_", names(result)))
    setnames(result, normalized)
    barcode_candidates <- c("barcode", "barcodes")
    barcode_column <- barcode_candidates[barcode_candidates %in% names(result)][1L]
    if (length(barcode_column) == 0L || is.na(barcode_column)) {
      stop("Position file header is present, but no barcode column was found.")
    }
    required_order <- c( barcode_column, "in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres" )
    missing_columns <- setdiff(required_order, names(result))
    if (length(missing_columns) > 0L) {
      stop( "Position file is missing required columns: ", paste(missing_columns, collapse = ", ") )
    }
    result <- result[, ..required_order]
    setnames( result, c( "barcode", "in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres" ) )
  } else {
    result <- result[, 1:6]
    setnames( result, c( "barcode", "in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres" ) )
  }
  result[, barcode := trimws(as.character(barcode))]
  result[, in_tissue := suppressWarnings(as.integer(in_tissue))]
  result
}
check_png_gz_signature <- function(path) {
  connection <- gzfile(path, open = "rb")
  on.exit(close(connection), add = TRUE)
  signature <- readBin(connection, what = "raw", n = 8L)
  expected <- as.raw(c(137, 80, 78, 71, 13, 10, 26, 10))
  length(signature) == 8L && identical(signature, expected)
}
check_scalefactor_gz <- function(path) {
  connection <- gzfile(path, open = "rt")
  on.exit(close(connection), add = TRUE)
  text <- paste(readLines(connection, warn = FALSE), collapse = "\n")
  required_keys <- c( "tissue_hires_scalef", "tissue_lowres_scalef", "spot_diameter_fullres" )
  present <- vapply( required_keys, function(key) grepl(paste0('"', key, '"'), text, fixed = TRUE), logical(1) )
  list( readable = nzchar(trimws(text)), all_required_keys = all(present), present_keys = paste(required_keys[present], collapse = " | "), missing_keys = paste(required_keys[!present], collapse = " | ") )
}
checks <- data.table( check_id = character(), category = character(), description = character(), status = character(), observed = character(), required = character(), critical = logical(), action_if_not_pass = character() )
add_check <- function( check_id, category, description, passed, observed, required, critical = TRUE, action_if_not_pass = "") {
  checks <<- rbind( checks, data.table( check_id = as.character(check_id), category = as.character(category), description = as.character(description), status = if (isTRUE(passed)) "PASS" else "FAIL", observed = as.character(observed), required = as.character(required), critical = isTRUE(critical), action_if_not_pass = as.character(action_if_not_pass) ), fill = TRUE )
}
add_review <- function( check_id, description, observed, required, action_if_not_pass) {
  checks <<- rbind( checks, data.table( check_id = as.character(check_id), category = "metadata_review", description = as.character(description), status = "REVIEW", observed = as.character(observed), required = as.character(required), critical = FALSE, action_if_not_pass = as.character(action_if_not_pass) ), fill = TRUE )
}
write_failure_report <- function(error_message) {
  report <- c( "============================================================", "GSE234774 Step 12A2 metadata / asset linkage repair", "============================================================", paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), paste0("Project root: ", PROJECT_ROOT), paste0("Input directory: ", INPUT_DIR), paste0("Output directory: ", OUTPUT_DIR), "", "FINAL_STATUS: STOP", "", "FATAL ERROR", paste0("- ", error_message), "", "SCOPE", "- No Seurat object was created.", "- No RCTD, target scoring, normalization, clustering, or spatial statistics were run.", "- Do not proceed to Step 12B until this script completes." )
  writeLines(report, REPORT_FILE, useBytes = TRUE)
  log_message("FATAL ERROR: ", error_message)
  write_run_log()
}
main <- function() {
  log_message("============================================================")
  log_message("GSE234774 Step 12A2 linkage repair started")
  log_message("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  log_message("============================================================")
  required_files <- c( matrix = MATRIX_FILE, barcodes = BARCODE_FILE, features = FEATURE_FILE, metadata = METADATA_FILE, raw_tar = RAW_TAR )
  required_exists <- file.exists(required_files)
  add_check( "C01", "input", "All five required spatial 2D input files exist.", all(required_exists), paste( paste(names(required_files), required_exists, sep = "="), collapse = "; " ), "matrix/barcodes/features/metadata/RAW.tar all TRUE", TRUE, "Restore missing local files. Do not download or replace files silently." )
  if (!all(required_exists)) {
    stop( "Missing required input file(s): ", paste(names(required_files)[!required_exists], collapse = ", ") )
  }
  matrix_dimensions <- read_matrix_market_dimensions(MATRIX_FILE)
  barcode_table <- read_single_column(BARCODE_FILE, "barcode")
  feature_table <- read_single_column(FEATURE_FILE, "feature")
  metadata <- read_spatial_metadata(METADATA_FILE)
  add_check( "C02", "matrix", "Matrix dimensions equal feature and barcode file lengths.", matrix_dimensions$n_features == nrow(feature_table) && matrix_dimensions$n_barcodes == nrow(barcode_table), paste0( "matrix=", matrix_dimensions$n_features, " x ", matrix_dimensions$n_barcodes, "; feature rows=", nrow(feature_table), "; barcode rows=", nrow(barcode_table) ), "Exact equality", TRUE, "Do not create a spatial object until all dimensions match." )
  required_metadata_columns <- c( "barcode", "id", "slide", "spot", "group", "replicate", "sections", "x", "y", "z" )
  missing_metadata_columns <- setdiff( required_metadata_columns, names(metadata) )
  add_check( "C03", "metadata", "Observed Step 12A2 metadata fields are present.", length(missing_metadata_columns) == 0L, if (length(missing_metadata_columns) == 0L) {
      paste(required_metadata_columns, collapse = " | ")
    } else {
      paste0( "missing=", paste(missing_metadata_columns, collapse = " | ") )
    },
    paste(required_metadata_columns, collapse = " | "), TRUE, "Do not substitute fields by memory; inspect the processed metadata header."
  )
  if (length(missing_metadata_columns) > 0L) {
    stop( "Required metadata columns are missing: ", paste(missing_metadata_columns, collapse = ", ") )
  }
  metadata[, barcode := trimws(as.character(barcode))]
  metadata[, id := trimws(as.character(id))]
  metadata[, slide := trimws(as.character(slide))]
  metadata[, spot := toupper(trimws(as.character(spot)))]
  metadata[, group := trimws(as.character(group))]
  metadata[, replicate := trimws(as.character(replicate))]
  metadata[, sections := trimws(as.character(sections))]
  metadata[, z := trimws(as.character(z))]
  barcode_table[, barcode := trimws(as.character(barcode))]
  barcode_exact_match <- setequal( barcode_table$barcode, metadata$barcode )
  barcode_unique <- !anyDuplicated(metadata$barcode) && !anyDuplicated(barcode_table$barcode)
  add_check( "C04", "metadata", "Processed barcode file and metadata barcode column match uniquely.", barcode_exact_match && barcode_unique && nrow(metadata) == nrow(barcode_table), paste0( "barcode file rows=", nrow(barcode_table), "; metadata rows=", nrow(metadata), "; setequal=", barcode_exact_match, "; unique=", barcode_unique ), "Equal row count, exact set match, no duplicates", TRUE, "Resolve barcode formatting or duplicate metadata rows before linkage." )
  metadata[, capture_area_id := spot]
  metadata[, section_id := id]
  metadata[, timepoint := group]
  metadata[, replicate_index := replicate]
  metadata[, biological_replicate_id := paste(group, replicate, sep = "_")]
  metadata[, injury_condition := ifelse( tolower(group) == "uninjured", "Uninjured", "SCI" )]
  metadata[, section_order := sections]
  metadata[, expected_section_id := paste( spot, group, replicate, sep = "_" )]
  id_consistency <- all( metadata$section_id == metadata$expected_section_id )
  add_check( "C05", "metadata", "The observed id field equals spot + group + replicate.", id_consistency, paste0( "matching rows=", sum(metadata$section_id == metadata$expected_section_id), "/", nrow(metadata) ), "All rows match", TRUE, "STOP: section_id cannot be reconstructed consistently from observed fields." )
  expected_timepoints <- c("uninjured", "7days", "2months")
  observed_timepoints <- sort(unique(metadata$timepoint))
  timepoint_set_ok <- setequal(observed_timepoints, expected_timepoints)
  add_check( "C06", "metadata", "The timepoint field is recovered from group with the observed three levels.", timepoint_set_ok, paste(observed_timepoints, collapse = " | "), paste(expected_timepoints, collapse = " | "), TRUE, "Review the group field before defining the spatial time course." )
  observed_capture_areas <- sort(unique(metadata$capture_area_id))
  observed_sections <- sort(unique(metadata$section_id))
  add_check( "C07", "metadata", "The processed spatial 2D object contains the audited 12 capture areas and 36 condition-specific sections.", length(observed_capture_areas) == 12L && length(observed_sections) == 36L, paste0( "capture areas=", length(observed_capture_areas), "; section_ids=", length(observed_sections) ), "12 capture areas; 36 section_ids", TRUE, "Verify that the spatial 2D metadata file, not spatial 3D metadata, was loaded." )
  barcode_prefix_ok <- vapply( observed_capture_areas, function(current_capture) {
      current_barcodes <- metadata[ capture_area_id == current_capture, barcode ]
      all(startsWith(current_barcodes, paste0(current_capture, "-")))
    },
    logical(1)
  )
  add_check( "C08", "metadata", "Every processed barcode prefix matches its capture_area_id.", all(barcode_prefix_ok), paste0( sum(barcode_prefix_ok), "/", length(barcode_prefix_ok), " capture areas pass" ), "All capture areas pass", TRUE, "Resolve spot-prefix parsing before raw position linkage." )
  section_manifest <- metadata[, .( n_spots = .N, n_unique_barcodes = uniqueN(barcode), capture_area_id = collapse_unique(capture_area_id), capture_area_n_unique = n_unique_nonmissing(capture_area_id), slide_id = collapse_unique(slide), slide_n_unique = n_unique_nonmissing(slide), timepoint = collapse_unique(timepoint), timepoint_n_unique = n_unique_nonmissing(timepoint), replicate_index = collapse_unique(replicate_index), replicate_index_n_unique = n_unique_nonmissing(replicate_index), biological_replicate_id = collapse_unique(biological_replicate_id), biological_replicate_n_unique = n_unique_nonmissing( biological_replicate_id ), injury_condition = collapse_unique(injury_condition), injury_condition_n_unique = n_unique_nonmissing(injury_condition), section_order = collapse_unique(section_order), section_order_n_unique = n_unique_nonmissing(section_order), z_coordinate = collapse_unique(z), z_n_unique = n_unique_nonmissing(z), x_min = suppressWarnings(min(as.numeric(x), na.rm = TRUE)), x_max = suppressWarnings(max(as.numeric(x), na.rm = TRUE)), y_min = suppressWarnings(min(as.numeric(y), na.rm = TRUE)), y_max = suppressWarnings(max(as.numeric(y), na.rm = TRUE)) ), by = section_id]
  section_manifest[, section_conflict := ( capture_area_n_unique != 1L | slide_n_unique != 1L | timepoint_n_unique != 1L | replicate_index_n_unique != 1L | biological_replicate_n_unique != 1L | injury_condition_n_unique != 1L | section_order_n_unique != 1L | z_n_unique != 1L )]
  add_check( "C09", "section_linkage", "Each condition-specific section maps to one capture area, timepoint, replicate, condition, section order, and z coordinate.", !any(section_manifest$section_conflict), paste0( "conflicting sections=", sum(section_manifest$section_conflict), "/", nrow(section_manifest) ), "Zero conflicting sections", TRUE, "Do not use section-level inference until all section conflicts are resolved." )
  capture_metadata <- metadata[, .( n_spots = .N, n_section_ids = uniqueN(section_id), section_ids = collapse_unique(section_id), n_timepoints = uniqueN(timepoint), timepoints = collapse_unique(timepoint), replicate_indices = collapse_unique(replicate_index), biological_replicate_ids = collapse_unique(biological_replicate_id), slide_ids = collapse_unique(slide), section_orders = collapse_unique(section_order), z_coordinates = collapse_unique(z) ), by = capture_area_id]
  capture_metadata[, complete_timepoint_set := vapply( strsplit(timepoints, " \\| "), function(x) setequal(x, expected_timepoints), logical(1) )]
  add_check( "C10", "section_linkage", "Every capture area contains one section for each audited timepoint.", all(capture_metadata$complete_timepoint_set) && all(capture_metadata$n_section_ids == length(expected_timepoints)), paste0( sum( capture_metadata$complete_timepoint_set & capture_metadata$n_section_ids == length(expected_timepoints) ), "/", nrow(capture_metadata), " capture areas complete" ), "All capture areas contain three timepoint-specific sections", TRUE, "Review section construction and missing processed spots." )
  mapping_table <- data.table( analysis_role = c( "barcode", "capture_area_id", "condition_specific_section_id", "slide_id", "timepoint", "replicate_index", "biological_replicate_id", "injury_condition", "section_order", "x_coordinate", "y_coordinate", "z_coordinate" ), source_column = c( "barcode", "spot", "id", "slide", "group", "replicate", "group + replicate", "group", "sections", "x", "y", "z" ), derivation = c( "direct", "direct", "direct; validated as spot_group_replicate", "direct", "direct", "direct", "paste(group, replicate, sep='_')", "uninjured -> Uninjured; all other observed timepoints -> SCI", "direct", "direct", "direct", "direct" ), confidence = c( "exact barcode match", "exact barcode-prefix match", "exact row-wise identity reconstruction", "observed metadata field", "observed metadata field and audited values", "observed metadata field", "derived condition-specific biological replicate identifier", "derived from observed timepoint labels", "observed metadata field", "observed metadata field", "observed metadata field", "observed metadata field" ), notes = c( NA, "Used to link RAW.tar image, position, and scalefactor assets.", "Primary section-level analysis identifier.", NA, "Observed levels: uninjured, 7days, 2months.", "Not globally unique when used alone.", "Use this identifier for biological-replicate grouping; preserve derivation provenance.", "This is a broad injury/control label, not a treatment label.", "Serial-section/order field; not used as the unique section_id.", NA, NA, "Spatial z coordinate/order field; not used as the unique section_id." ) )
  tar_members <- utils::untar(RAW_TAR, list = TRUE)
  if (length(tar_members) == 0L) {
    stop("RAW.tar member listing returned zero files.")
  }
  tar_inventory <- data.table(member_path = tar_members)
  tar_inventory[, member_basename := basename(member_path)]
  tar_inventory[, gsm_id := extract_gsm(member_path)]
  tar_inventory[, capture_area_id := extract_capture_area(member_path)]
  tar_inventory[, asset_category := classify_tar_asset(member_path)]
  tar_inventory[, unsafe_member_path := is_unsafe_tar_member(member_path)]
  tar_inventory[, belongs_to_processed_2d := ( capture_area_id %in% observed_capture_areas )]
  add_check( "C11", "archive", "RAW.tar lists successfully and contains no unsafe member path.", nrow(tar_inventory) > 0L && !any(tar_inventory$unsafe_member_path), paste0( "members=", nrow(tar_inventory), "; unsafe=", sum(tar_inventory$unsafe_member_path) ), "Non-empty listing; zero unsafe paths", TRUE, "Do not extract an archive containing absolute or parent-traversal paths." )
  relevant_tar <- tar_inventory[ belongs_to_processed_2d == TRUE & asset_category %in% c( "position", "image_hires", "image_lowres", "scalefactor" ) ]
  raw_capture_areas <- sort(unique( tar_inventory[ !is.na(capture_area_id), capture_area_id ] ))
  raw_capture_match <- setequal( observed_capture_areas, sort(unique(relevant_tar$capture_area_id)) )
  add_check( "C12", "archive", "All processed 2D capture areas are represented in RAW.tar.", raw_capture_match, paste0( "processed=", paste(observed_capture_areas, collapse = " | "), "; matched RAW=", paste(sort(unique(relevant_tar$capture_area_id)), collapse = " | ") ), "Exact processed capture-area set match", TRUE, "Do not link assets by GSM alone; inspect capture-area tokens." )
  asset_counts <- dcast( relevant_tar[ , .N, by = .(capture_area_id, asset_category) ], capture_area_id ~ asset_category, value.var = "N", fill = 0L )
  for (required_category in c( "position", "image_hires", "image_lowres", "scalefactor" )) {
    if (!required_category %in% names(asset_counts)) {
      asset_counts[, (required_category) := 0L]
    }
  }
  exact_asset_counts <- all( asset_counts$position == 1L & asset_counts$image_hires == 1L & asset_counts$image_lowres == 1L & asset_counts$scalefactor == 1L ) && nrow(asset_counts) == length(observed_capture_areas)
  add_check( "C13", "spatial_asset", "Every processed 2D capture area has one position file, two compressed PNG images, and one scalefactor JSON.", exact_asset_counts, paste0( "complete capture areas=", sum( asset_counts$position == 1L & asset_counts$image_hires == 1L & asset_counts$image_lowres == 1L & asset_counts$scalefactor == 1L ), "/", length(observed_capture_areas) ), "12/12 complete capture areas", TRUE, "Review duplicate or missing RAW.tar assets before spatial loading." )
  if (any(tar_inventory$unsafe_member_path)) {
    stop("Unsafe TAR member path detected. Extraction has been blocked.")
  }
  if (dir.exists(EXTRACT_DIR)) {
    unlink(EXTRACT_DIR, recursive = TRUE, force = TRUE)
  }
  dir.create(EXTRACT_DIR, recursive = TRUE, showWarnings = FALSE)
  selected_members <- sort(unique(relevant_tar$member_path))
  utils::untar( RAW_TAR, files = selected_members, exdir = EXTRACT_DIR )
  relevant_tar[, extracted_path := file.path( EXTRACT_DIR, member_path )]
  relevant_tar[, extracted_exists := file.exists(extracted_path)]
  relevant_tar[, extracted_size_bytes := ifelse( extracted_exists, file.info(extracted_path)$size, NA_real_ )]
  extraction_ok <- all(relevant_tar$extracted_exists) && all(relevant_tar$extracted_size_bytes > 0, na.rm = TRUE)
  add_check( "C14", "controlled_extraction", "All processed 2D assets were extracted to the controlled audit directory.", extraction_ok && nrow(relevant_tar) == 4L * length(observed_capture_areas), paste0( "extracted=", sum(relevant_tar$extracted_exists), "/", nrow(relevant_tar), "; directory=", EXTRACT_DIR ), "48/48 linked compressed assets extracted", TRUE, "Do not manually unpack or rename files; review TAR extraction errors." )
  relevant_tar[, asset_readable := FALSE]
  relevant_tar[, asset_check_detail := NA_character_]
  image_rows <- which( relevant_tar$asset_category %in% c( "image_hires", "image_lowres" ) )
  for (row_index in image_rows) {
    current_path <- relevant_tar$extracted_path[row_index]
    current_ok <- FALSE
    current_detail <- ""
    try_result <- tryCatch( {
        current_ok <- check_png_gz_signature(current_path)
        current_detail <- if (current_ok) {
          "Decompressed PNG signature is valid."
        } else {
          "Decompressed bytes do not match the PNG signature."
        }
        TRUE
      },
      error = function(e) {
        current_detail <<- paste0("PNG read error: ", conditionMessage(e))
        FALSE
      }
    )
    relevant_tar$asset_readable[row_index] <- isTRUE(try_result) && isTRUE(current_ok)
    relevant_tar$asset_check_detail[row_index] <- current_detail
  }
  png_ok <- nrow(relevant_tar[ asset_category %in% c("image_hires", "image_lowres") ]) == 2L * length(observed_capture_areas) && all(relevant_tar[ asset_category %in% c("image_hires", "image_lowres"), asset_readable ])
  add_check( "C15", "spatial_asset", "All linked .png.gz files decompress to a valid PNG signature.", png_ok, paste0( "valid PNG files=", sum(relevant_tar[ asset_category %in% c("image_hires", "image_lowres"), asset_readable ]), "/", nrow(relevant_tar[ asset_category %in% c("image_hires", "image_lowres") ]) ), "24/24 valid compressed PNG files", TRUE, "Do not create a Seurat image object until every linked image is readable." )
  scalefactor_rows <- which( relevant_tar$asset_category == "scalefactor" )
  for (row_index in scalefactor_rows) {
    current_path <- relevant_tar$extracted_path[row_index]
    current_detail <- ""
    result <- tryCatch( check_scalefactor_gz(current_path), error = function(e) {
        list( readable = FALSE, all_required_keys = FALSE, present_keys = "", missing_keys = paste0("read error: ", conditionMessage(e)) )
      }
    )
    relevant_tar$asset_readable[row_index] <- isTRUE(result$readable) && isTRUE(result$all_required_keys)
    current_detail <- paste0( "present=", result$present_keys, "; missing=", ifelse(nzchar(result$missing_keys), result$missing_keys, "none") )
    relevant_tar$asset_check_detail[row_index] <- current_detail
  }
  scalefactors_ok <- nrow(relevant_tar[ asset_category == "scalefactor" ]) == length(observed_capture_areas) && all(relevant_tar[ asset_category == "scalefactor", asset_readable ])
  add_check( "C16", "spatial_asset", "Every linked scalefactor JSON is readable and contains the required scale keys.", scalefactors_ok, paste0( "valid scalefactors=", sum(relevant_tar[ asset_category == "scalefactor", asset_readable ]), "/", nrow(relevant_tar[ asset_category == "scalefactor" ]) ), "12/12 readable scalefactor files", TRUE, "Review corrupted or nonstandard scalefactor files before image loading." )
  position_audit_list <- vector( mode = "list", length = length(observed_capture_areas) )
  for (capture_index in seq_along(observed_capture_areas)) {
    current_capture <- observed_capture_areas[capture_index]
    position_row <- relevant_tar[ capture_area_id == current_capture & asset_category == "position" ]
    if (nrow(position_row) != 1L) {
      position_audit_list[[capture_index]] <- data.table( capture_area_id = current_capture, position_member = collapse_unique(position_row$member_path), position_read_status = "FAIL", position_read_message = paste0( "Expected one position file; observed ", nrow(position_row), "." ), processed_barcode_n = metadata[ capture_area_id == current_capture, uniqueN(barcode) ], position_barcode_n = NA_integer_, position_barcode_unique = FALSE, processed_prefix_valid = FALSE, matching_mode = NA_character_, processed_in_position_n = NA_integer_, processed_in_position_fraction = NA_real_, processed_in_tissue_n = NA_integer_, processed_in_tissue_fraction = NA_real_, position_in_tissue_n = NA_integer_, missing_processed_barcode_examples = NA_character_ )
      next
    }
    current_path <- position_row$extracted_path[1L]
    current_processed_full <- metadata[ capture_area_id == current_capture, unique(barcode) ]
    current_prefix <- paste0(current_capture, "-")
    prefix_valid <- all(startsWith( current_processed_full, current_prefix ))
    current_processed_unprefixed <- sub( paste0("^", regex_escape(current_prefix)), "", current_processed_full, perl = TRUE )
    position_result <- tryCatch( {
        position_table <- read_position_file(current_path)
        raw_position_barcodes <- unique(position_table$barcode)
        overlap_unprefixed <- sum( current_processed_unprefixed %in% raw_position_barcodes )
        overlap_full <- sum( current_processed_full %in% raw_position_barcodes )
        if (overlap_unprefixed >= overlap_full) {
          matching_mode <- "processed barcode after capture-area prefix removal"
          processed_for_match <- current_processed_unprefixed
        } else {
          matching_mode <- "processed full barcode"
          processed_for_match <- current_processed_full
        }
        processed_in_position <- processed_for_match %in%
          raw_position_barcodes
        in_tissue_barcodes <- unique(position_table[ !is.na(in_tissue) & in_tissue == 1L, barcode ])
        processed_in_tissue <- processed_for_match %in%
          in_tissue_barcodes
        missing_examples <- head( processed_for_match[!processed_in_position], 10L )
        data.table( capture_area_id = current_capture, position_member = position_row$member_path[1L], position_read_status = "PASS", position_read_message = "Position file read successfully.", processed_barcode_n = length(processed_for_match), position_barcode_n = length(raw_position_barcodes), position_barcode_unique = !anyDuplicated( position_table$barcode ), processed_prefix_valid = prefix_valid, matching_mode = matching_mode, processed_in_position_n = sum(processed_in_position), processed_in_position_fraction = mean( processed_in_position ), processed_in_tissue_n = sum(processed_in_tissue), processed_in_tissue_fraction = mean( processed_in_tissue ), position_in_tissue_n = length(in_tissue_barcodes), missing_processed_barcode_examples = if ( length(missing_examples) == 0L ) {
            NA_character_
          } else {
            paste(missing_examples, collapse = " | ")
          }
        )
      },
      error = function(e) {
        data.table( capture_area_id = current_capture, position_member = position_row$member_path[1L], position_read_status = "FAIL", position_read_message = conditionMessage(e), processed_barcode_n = length(current_processed_full), position_barcode_n = NA_integer_, position_barcode_unique = FALSE, processed_prefix_valid = prefix_valid, matching_mode = NA_character_, processed_in_position_n = NA_integer_, processed_in_position_fraction = NA_real_, processed_in_tissue_n = NA_integer_, processed_in_tissue_fraction = NA_real_, position_in_tissue_n = NA_integer_, missing_processed_barcode_examples = NA_character_ )
      }
    )
    position_audit_list[[capture_index]] <- position_result
  }
  position_audit <- rbindlist( position_audit_list, fill = TRUE )
  position_read_ok <- all( position_audit$position_read_status == "PASS" )
  position_overlap_ok <- position_read_ok && all( !is.na(position_audit$processed_in_position_fraction) & position_audit$processed_in_position_fraction == 1 )
  add_check( "C17", "position_barcode", "Every processed barcode is present in the matched RAW position file.", position_overlap_ok, paste0( "capture areas with 100% overlap=", sum( position_audit$processed_in_position_fraction == 1, na.rm = TRUE ), "/", nrow(position_audit) ), "12/12 capture areas with 100% processed-to-position overlap", TRUE, "Do not create a spatial image object until barcode-position overlap is exact." )
  in_tissue_overlap_ok <- position_read_ok && all( !is.na(position_audit$processed_in_tissue_fraction) & position_audit$processed_in_tissue_fraction == 1 )
  add_check( "C18", "position_barcode", "Every processed filtered barcode is marked in_tissue in the matched position file.", in_tissue_overlap_ok, paste0( "capture areas with 100% in-tissue overlap=", sum( position_audit$processed_in_tissue_fraction == 1, na.rm = TRUE ), "/", nrow(position_audit) ), "12/12 capture areas with 100% processed-to-in-tissue overlap", TRUE, "Review whether the processed matrix and position file derive from the same capture area." )
  relevant_asset_wide <- dcast( relevant_tar, capture_area_id ~ asset_category, value.var = "member_path", fun.aggregate = function(x) collapse_unique(x) )
  relevant_readable_wide <- dcast( relevant_tar, capture_area_id ~ asset_category, value.var = "asset_readable", fun.aggregate = function(x) all(x), fill = FALSE )
  readable_names <- setdiff( names(relevant_readable_wide), "capture_area_id" )
  setnames( relevant_readable_wide, readable_names, paste0(readable_names, "_readable") )
  capture_asset_manifest <- merge( capture_metadata, asset_counts, by = "capture_area_id", all.x = TRUE, sort = TRUE )
  capture_asset_manifest <- merge( capture_asset_manifest, relevant_asset_wide, by = "capture_area_id", all.x = TRUE, sort = TRUE, suffixes = c("_count", "_member") )
  capture_asset_manifest <- merge( capture_asset_manifest, relevant_readable_wide, by = "capture_area_id", all.x = TRUE, sort = TRUE )
  capture_asset_manifest <- merge( capture_asset_manifest, position_audit[ , .( capture_area_id, processed_in_position_fraction, processed_in_tissue_fraction ) ], by = "capture_area_id", all.x = TRUE, sort = TRUE )
  capture_asset_manifest[, linkage_complete := ( position_count == 1L & image_hires_count == 1L & image_lowres_count == 1L & scalefactor_count == 1L & position_readable & image_hires_readable & image_lowres_readable & scalefactor_readable & processed_in_position_fraction == 1 & processed_in_tissue_fraction == 1 )]
  add_check( "C19", "final_linkage", "Every processed capture area has complete and readable position/image/scalefactor linkage.", all(capture_asset_manifest$linkage_complete), paste0( "complete=", sum(capture_asset_manifest$linkage_complete), "/", nrow(capture_asset_manifest) ), "12/12 complete capture areas", TRUE, "STOP before Step 12B and inspect the failed capture-area rows." )
  non2d_raw_inventory <- tar_inventory[ !is.na(capture_area_id) & !capture_area_id %in% observed_capture_areas, .( gsm_id = collapse_unique(gsm_id), n_members = .N, asset_categories = collapse_unique(asset_category), member_examples = paste( head(member_basename, 4L), collapse = " | " ) ), by = capture_area_id ][order(capture_area_id)]
  add_check( "C20", "scope", "RAW capture areas outside the processed spatial 2D set are documented and excluded from 2D linkage.", length(setdiff(raw_capture_areas, observed_capture_areas)) > 0L, paste0( "excluded RAW capture areas=", paste( setdiff(raw_capture_areas, observed_capture_areas), collapse = " | " ) ), "Document, do not merge into the processed spatial 2D object", FALSE, "Keep the non-2D RAW assets separate from the current spatial 2D workflow." )
  add_review( "R_AGE", "Age is not represented by a dedicated field in the processed spatial 2D metadata.", paste0( "available fields=", paste(names(metadata), collapse = " | ") ), "Confirm adult cohort status from authoritative source metadata before manuscript wording.", "Do not infer age from filenames or from the snRNA cohort." )
  add_review( "R_SEVERITY_TREATMENT", "Severity and treatment are not represented by dedicated fields in the processed spatial 2D metadata.", "No dedicated severity or treatment columns were observed.", "Document source-level cohort selection separately.", "Do not manufacture these fields in the spatial object." )
  fwrite( mapping_table, file.path(OUTPUT_DIR, "01_metadata_role_mapping.csv"), bom = TRUE )
  fwrite( capture_asset_manifest, file.path(OUTPUT_DIR, "02_2d_capture_area_asset_manifest.csv"), bom = TRUE )
  fwrite( section_manifest[order( timepoint, biological_replicate_id, capture_area_id, section_id )], file.path( OUTPUT_DIR, "03_condition_specific_section_manifest.csv" ), bom = TRUE )
  fwrite( position_audit[order(capture_area_id)], file.path( OUTPUT_DIR, "04_position_barcode_overlap_audit.csv" ), bom = TRUE )
  fwrite( non2d_raw_inventory, file.path( OUTPUT_DIR, "05_non2d_RAW_capture_area_inventory.csv" ), bom = TRUE )
  critical_failures <- checks[ critical == TRUE & status == "FAIL" ]
  review_items <- checks[status == "REVIEW"]
  final_status <- if (nrow(critical_failures) > 0L) {
    "STOP"
  } else if (nrow(review_items) > 0L) {
    "PASS_WITH_REVIEW"
  } else {
    "PASS"
  }
  fwrite( checks, file.path(OUTPUT_DIR, "09_STEP12A2_audit_checks.csv"), bom = TRUE )
  report <- c( "============================================================", "GSE234774 Step 12A2 spatial metadata / asset linkage repair", "============================================================", paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), paste0("Project root: ", PROJECT_ROOT), paste0("Input directory: ", INPUT_DIR), paste0("Output directory: ", OUTPUT_DIR), paste0("Controlled asset directory: ", EXTRACT_DIR), "", paste0("FINAL_STATUS: ", final_status), "", "SCOPE", "- Metadata and raw-asset linkage repair only.", "- No Seurat object was created.", "- No normalization, clustering, RCTD, target scoring, or spatial statistics were run.", "- Only RAW.tar assets linked to processed spatial 2D capture areas were extracted.", "", "PROCESSED SPATIAL 2D INPUT", paste0( "- Matrix dimensions: features=", matrix_dimensions$n_features, "; barcodes=", matrix_dimensions$n_barcodes, "; nonzero=", format( matrix_dimensions$n_nonzero, scientific = FALSE, trim = TRUE ) ), paste0("- Metadata rows: ", nrow(metadata)), paste0( "- Capture areas: ", length(observed_capture_areas), " [", paste(observed_capture_areas, collapse = ", "), "]" ), paste0("- Condition-specific sections: ", nrow(section_manifest)), paste0( "- Timepoints: ", paste(observed_timepoints, collapse = ", ") ), "", "REPAIRED METADATA ROLES", "- barcode = barcode", "- capture_area_id = spot", "- section_id = id", "- timepoint = group", "- replicate_index = replicate", "- biological_replicate_id = paste(group, replicate, sep='_')", "- injury_condition = Uninjured for uninjured; SCI for injured timepoints", "- section_order = sections", "- z_coordinate = z", "", "RAW ASSET LINKAGE", paste0( "- Linked processed 2D assets: ", nrow(relevant_tar), " compressed files." ), paste0( "- Fully linked capture areas: ", sum(capture_asset_manifest$linkage_complete), "/", nrow(capture_asset_manifest) ), paste0( "- Valid compressed PNG files: ", sum(relevant_tar[ asset_category %in% c("image_hires", "image_lowres"), asset_readable ]), "/", nrow(relevant_tar[ asset_category %in% c("image_hires", "image_lowres") ]) ), paste0( "- Valid scalefactor files: ", sum(relevant_tar[ asset_category == "scalefactor", asset_readable ]), "/", nrow(relevant_tar[ asset_category == "scalefactor" ]) ), paste0( "- Capture areas with 100% processed-to-position overlap: ", sum( position_audit$processed_in_position_fraction == 1, na.rm = TRUE ), "/", nrow(position_audit) ), paste0( "- Capture areas with 100% processed-to-in-tissue overlap: ", sum( position_audit$processed_in_tissue_fraction == 1, na.rm = TRUE ), "/", nrow(position_audit) ), "", "CRITICAL CHECKS", paste0( "- PASS: ", nrow(checks[critical == TRUE & status == "PASS"]), "; FAIL: ", nrow(critical_failures) ) )
  if (nrow(critical_failures) > 0L) {
    report <- c( report, "", "STOP REASONS", paste0( "- ", critical_failures$check_id, ": ", critical_failures$description, " Observed: ", critical_failures$observed, " Action: ", critical_failures$action_if_not_pass ) )
  }
  if (nrow(review_items) > 0L) {
    report <- c( report, "", "REVIEW ITEMS", paste0( "- ", review_items$check_id, ": ", review_items$description, " Observed: ", review_items$observed, " Action: ", review_items$action_if_not_pass ) )
  }
  report <- c( report, "", "NEXT ACTION", if (final_status == "STOP") {
      paste( "Step 12A2 STOP.", "review the seven key outputs.", "Do not create a Seurat spatial object or run RCTD." )
    } else {
      paste( "Step 12A2 linkage PASS.", "The audited 12 spatial 2D capture areas may proceed to Step 12B", "for section-wise spatial object construction and QC.", "Age/severity/treatment remain source-metadata review items." )
    },
    "", "KEY OUTPUTS", "- 01_metadata_role_mapping.csv", "- 02_2d_capture_area_asset_manifest.csv", "- 03_condition_specific_section_manifest.csv", "- 04_position_barcode_overlap_audit.csv", "- 05_non2d_RAW_capture_area_inventory.csv", "- 09_STEP12A2_audit_checks.csv", "- 10_STEP12A2_FINAL_REPORT.txt", "============================================================"
  )
  writeLines(report, REPORT_FILE, useBytes = TRUE)
  log_message("FINAL_STATUS: ", final_status)
  log_message("Output directory: ", OUTPUT_DIR)
  write_run_log()
  cat("\n")
  cat("============================================================\n")
  cat("Step 12A2 completed: ", final_status, "\n", sep = "")
  cat("Open first: ", REPORT_FILE, "\n", sep = "")
  cat("============================================================\n")
  invisible( list( final_status = final_status, checks = checks, mapping = mapping_table, capture_asset_manifest = capture_asset_manifest, section_manifest = section_manifest, position_audit = position_audit ) )
}
tryCatch( main(), error = function(e) {
    write_failure_report(conditionMessage(e))
    message( "Step 12A2 stopped with a fatal error: ", conditionMessage(e) )
  }
)
