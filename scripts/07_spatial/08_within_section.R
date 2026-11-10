options( stringsAsFactors = FALSE, warn = 1, width = 180, future.globals.maxSize = 24 * 1024^3 )
PROJECT_ROOT <- Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
INJURED_TIMEPOINTS <- c( "7days", "2months" )
UNINJURED_TIMEPOINT <- "uninjured"
STRUCTURAL_POSITIVE_TYPES <- c( "BAM_macrophage", "Fibroblast_VLMC", "Endothelial", "Astrocyte", "Pericyte" )
STRUCTURAL_NEGATIVE_TYPES <- c( "Neuron", "Oligodendrocyte", "OPC" )
K_SMOOTH <- 12L
CORE_SEED_QUANTILE <- 0.90
SEED_CONNECTIVITY_RADIUS_NN <- 1.80
MIN_CORE_SPOTS <- 12L
MIN_ZONE_SPOTS <- 20L
BORDER_MAX_DISTANCE_NN <- 2.50
PERILESIONAL_MAX_DISTANCE_NN <- 6.00
FALLBACK_BORDER_DISTANCE_QUANTILE <- 0.25
FALLBACK_PERILESIONAL_DISTANCE_QUANTILE <- 0.60
DISTANCE_PROFILE_BREAKS_NN <- c( 0, 1, 2, 3, 4, 6, 8, 12, Inf )
POINT_CEX <- 0.45
PLOT_QUANTILE_LOW <- 0.01
PLOT_QUANTILE_HIGH <- 0.99
EXPECTED_SPOTS <- 33941L
EXPECTED_CAPTURE_AREAS <- 12L
EXPECTED_SECTIONS <- 36L
EXPECTED_INJURED_SECTIONS <- 24L
EXPECTED_ANIMALS_PER_INJURED_TIMEPOINT <- 3L
STEP12E_DIR <- file.path( PROJECT_ROOT, "04_results", "step12E_GSE234774_target_program_spatial_scoring_microglia_adjustment" )
STEP12E_REPORT <- file.path( STEP12E_DIR, "10_STEP12E_FINAL_REPORT.txt" )
STEP12E_SPOT_SCORES_RDS <- file.path( STEP12E_DIR, "RDS", "GSE234774_step12E_spot_target_shared_microglia_scores.rds" )
STEP12D_DIR <- file.path( PROJECT_ROOT, "04_results", "step12D_GSE234774_full_mode_RCTD_per_capture_area" )
STEP12D_REPORT <- file.path( STEP12D_DIR, "10_STEP12D_FINAL_REPORT.txt" )
STEP12D_NORMALIZED_WEIGHTS_RDS <- file.path( STEP12D_DIR, "RDS", "GSE234774_step12D_all_spots_normalized_RCTD_weights.rds" )
STEP12D_CAPTURE_ROOT <- file.path( STEP12D_DIR, "capture_areas" )
STEP12B_DIR <- file.path( PROJECT_ROOT, "04_results", "step12B_fast_BPCells_spatial_section_QC" )
STEP12B_REPORT <- file.path( STEP12B_DIR, "10_STEP12B_FAST_FINAL_REPORT.txt" )
STEP12B_CAPTURE_AUDIT <- file.path( STEP12B_DIR, "tables", "02_capture_area_image_coordinate_audit.csv" )
LOWRES_PREVIEW_ROOT <- file.path( STEP12B_DIR, "lowres_asset_preview" )
OUTPUT_DIR <- file.path( PROJECT_ROOT, "04_results", "step12F_GSE234774_R2MC4_shared_within_section_revision" )
TABLE_DIR <- file.path( OUTPUT_DIR, "tables" )
FIGURE_DIR <- file.path( OUTPUT_DIR, "figures" )
RDS_DIR <- file.path( OUTPUT_DIR, "RDS" )
RUN_LOG_FILE <- file.path( OUTPUT_DIR, "00_run_log.txt" )
SESSION_INFO_FILE <- file.path( OUTPUT_DIR, "00_sessionInfo.txt" )
OUT_INPUT_AUDIT <- file.path( TABLE_DIR, "01_step12F_input_and_zone_definition_audit.csv" )
OUT_SPOT_ZONES <- file.path( TABLE_DIR, "02_spot_independent_lesion_zones_and_scores.csv.gz" )
OUT_SECTION_ZONE_SUMMARY <- file.path( TABLE_DIR, "03_section_zone_score_summary.csv" )
OUT_SECTION_CONTRASTS <- file.path( TABLE_DIR, "04_section_zone_vs_remote_contrasts.csv" )
OUT_ANIMAL_CONTRASTS <- file.path( TABLE_DIR, "05_animal_zone_vs_remote_contrasts.csv" )
OUT_REVISED_TABLE_S10 <- file.path( TABLE_DIR, "08_Table_S10_REVISED.xlsx" )
OUT_TABLE_S10_AUDIT <- file.path( TABLE_DIR, "08A_Table_S10_provenance_audit.csv" )
ORIGINAL_STEP12F_ANIMAL_CSV <- file.path( PROJECT_ROOT, "04_results", "step12F_GSE234774_independent_lesion_zones_target_program_validation", "tables", "05_animal_zone_vs_remote_contrasts.csv" )
STEP12H_ANIMAL_CSV <- file.path( PROJECT_ROOT, "04_results", "step12H_GSE234774_spatial_cell_source_sensitivity", "tables", "04_animal_zone_vs_remote_contrasts.csv" )
OUT_ANIMAL_INFERENCE <- file.path( TABLE_DIR, "06_animal_level_inference_summary.csv" )
OUT_DISTANCE_PROFILE <- file.path( TABLE_DIR, "07_section_distance_to_core_profiles.csv" )
OUT_CHECKS <- file.path( TABLE_DIR, "09_STEP12F_audit_checks.csv" )
OUT_ZONE_MAPS <- file.path( FIGURE_DIR, "01_independent_lesion_zone_and_target_validation_maps.pdf" )
OUT_ANIMAL_CONTRAST_FIGURE <- file.path( FIGURE_DIR, "02_animal_level_target_zone_contrasts.pdf" )
OUT_DISTANCE_PROFILE_FIGURE <- file.path( FIGURE_DIR, "03_target_score_distance_to_core_profiles.pdf" )
OUT_FINAL_REPORT <- file.path( OUTPUT_DIR, "10_STEP12F_FINAL_REPORT.txt" )
OUT_SPOT_ZONES_RDS <- file.path( RDS_DIR, "GSE234774_step12F_spot_independent_lesion_zones_and_scores.rds" )
OUT_HANDOFF_RDS <- file.path( RDS_DIR, "GSE234774_step12F_independent_lesion_zone_validation_handoff.rds" )
for (directory in c( OUTPUT_DIR, TABLE_DIR, FIGURE_DIR, RDS_DIR )) {
  dir.create( directory, recursive = TRUE, showWarnings = FALSE )
}
required_packages <- c( "data.table", "Matrix", "png" )
missing_packages <- required_packages[ !vapply( required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1) ) ]
if (length(missing_packages) > 0L) {
  stop( "Missing required package(s): ", paste( missing_packages, collapse = ", " ), ". Install them before rerunning Step 12F.", call. = FALSE )
}
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(png)
})
run_log <- character()
write_run_log <- function() {
  writeLines( run_log, RUN_LOG_FILE, useBytes = TRUE )
}
elapsed_minutes <- function(start_time) {
  round( as.numeric( difftime( Sys.time(), start_time, units = "mins" ) ), 3 )
}
file_size_safe <- function(path) {
  if (!file.exists(path)) {
    return(NA_real_)
  }
  as.numeric( file.info(path)$size )
}
required_columns <- function( table, columns, label) {
  missing <- setdiff( columns, names(table) )
  if (length(missing) > 0L) {
    stop( label, " is missing required column(s): ", paste( missing, collapse = ", " ), call. = FALSE )
  }
  invisible(TRUE)
}
collapse_unique <- function( x, max_values = 50L) {
  x <- trimws( as.character(x) )
  x <- unique( x[ !is.na(x) & nzchar(x) ] )
  if (length(x) == 0L) {
    return(NA_character_)
  }
  if (length(x) > max_values) {
    x <- c( x[seq_len(max_values)], paste0( "... +", length(x) - max_values, " more" ) )
  }
  paste( x, collapse = " | " )
}
robust_z <- function(x) {
  x <- as.numeric(x)
  output <- rep( NA_real_, length(x) )
  valid <- is.finite(x)
  if (sum(valid) < 2L) {
    output[valid] <- 0
    return(output)
  }
  center <- stats::median( x[valid] )
  scale_value <- stats::mad( x[valid], center = center, constant = 1.4826 )
  if ( !is.finite(scale_value) || scale_value <= 0 ) {
    scale_value <- stats::sd( x[valid] )
  }
  if ( !is.finite(scale_value) || scale_value <= 0 ) {
    output[valid] <- 0
    return(output)
  }
  output[valid] <- ( x[valid] - center ) / scale_value
  output[valid] <- pmax( -6, pmin( 6, output[valid] ) )
  output
}
safe_mean <- function(x) {
  x <- x[ is.finite(x) ]
  if (length(x) == 0L) {
    return(NA_real_)
  }
  mean(x)
}
safe_median <- function(x) {
  x <- x[ is.finite(x) ]
  if (length(x) == 0L) {
    return(NA_real_)
  }
  stats::median(x)
}
safe_sd <- function(x) {
  x <- x[ is.finite(x) ]
  if (length(x) < 2L) {
    return(NA_real_)
  }
  stats::sd(x)
}
safe_quantile <- function( x, probability) {
  x <- x[ is.finite(x) ]
  if (length(x) == 0L) {
    return(NA_real_)
  }
  as.numeric( stats::quantile( x, probs = probability, names = FALSE, na.rm = TRUE, type = 8 ) )
}
safe_spearman <- function( x, y) {
  valid <- is.finite(x) & is.finite(y)
  if ( sum(valid) < 4L || stats::sd(x[valid]) == 0 || stats::sd(y[valid]) == 0 ) {
    return(NA_real_)
  }
  suppressWarnings( stats::cor( x[valid], y[valid], method = "spearman" ) )
}
exact_sign_test_p <- function(x) {
  x <- x[ is.finite(x) & x != 0 ]
  n <- length(x)
  if (n == 0L) {
    return(NA_real_)
  }
  positives <- sum( x > 0 )
  lower_tail <- stats::pbinom( min( positives, n - positives ), size = n, prob = 0.5 )
  min( 1, 2 * lower_tail )
}
safe_wilcoxon_signed_p <- function(x) {
  x <- x[ is.finite(x) ]
  if ( length(x) < 2L || all(x == 0) ) {
    return(NA_real_)
  }
  tryCatch( {
      suppressWarnings( stats::wilcox.test( x, mu = 0, alternative = "two.sided", exact = TRUE, paired = FALSE, conf.int = FALSE )$p.value )
    },
    error = function(e) {
      NA_real_
    }
  )
}
safe_wilcoxon_rank_p <- function( x, y) {
  x <- x[ is.finite(x) ]
  y <- y[ is.finite(y) ]
  if ( length(x) < 2L || length(y) < 2L ) {
    return(NA_real_)
  }
  tryCatch( {
      suppressWarnings( stats::wilcox.test( x, y, alternative = "two.sided", exact = TRUE, paired = FALSE, conf.int = FALSE )$p.value )
    },
    error = function(e) {
      NA_real_
    }
  )
}
pairwise_distance_matrix <- function( x, y) {
  coordinates <- cbind( as.numeric(x), as.numeric(y) )
  if ( nrow(coordinates) < 2L || any( !is.finite( coordinates ) ) ) {
    stop( "Invalid coordinates supplied to the section distance matrix.", call. = FALSE )
  }
  as.matrix( stats::dist( coordinates, method = "euclidean", diag = TRUE, upper = TRUE ) )
}
local_nn_spacing <- function( distance_matrix) {
  diag(distance_matrix) <- Inf
  nearest <- apply( distance_matrix, 1L, min, na.rm = TRUE )
  spacing <- stats::median( nearest[ is.finite(nearest) & nearest > 0 ] )
  if ( !is.finite(spacing) || spacing <= 0 ) {
    stop( "Could not estimate positive nearest-neighbor spacing.", call. = FALSE )
  }
  spacing
}
smooth_index_knn <- function( structural_index, distance_matrix, k) {
  n <- length( structural_index )
  k_use <- min( as.integer(k), n - 1L )
  if (k_use < 1L) {
    return( structural_index )
  }
  output <- numeric(n)
  for (i in seq_len(n)) {
    ordered <- order( distance_matrix[i, ], na.last = NA )
    ordered <- ordered[ ordered != i ]
    neighbors <- head( ordered, k_use )
    neighbor_distances <- distance_matrix[ i, neighbors ]
    positive_distances <- neighbor_distances[ is.finite(neighbor_distances) & neighbor_distances > 0 ]
    local_distance_scale <- if ( length(positive_distances) > 0L ) {
      stats::median( positive_distances )
    } else {
      1
    }
    if ( !is.finite(local_distance_scale) || local_distance_scale <= 0 ) {
      local_distance_scale <- 1
    }
    weights <- 1 / ( 1 + neighbor_distances / local_distance_scale )
    values <- structural_index[ neighbors ]
    valid <- is.finite(values) & is.finite(weights) & weights > 0
    self_value <- structural_index[i]
    if (!is.finite(self_value)) {
      self_value <- 0
    }
    output[i] <- ( self_value + sum( values[valid] * weights[valid] ) ) / ( 1 + sum( weights[valid] ) )
  }
  output
}
connected_components_from_distance <- function( candidate_indices, distance_matrix, maximum_distance) {
  candidate_indices <- unique( as.integer( candidate_indices ) )
  if (length(candidate_indices) == 0L) {
    return(list())
  }
  unvisited <- candidate_indices
  components <- list()
  while (length(unvisited) > 0L) {
    seed <- unvisited[[1L]]
    queue <- seed
    component <- integer()
    while (length(queue) > 0L) {
      current <- queue[[1L]]
      queue <- queue[-1L]
      if (current %in% component) {
        next
      }
      component <- c( component, current )
      neighbors <- candidate_indices[ distance_matrix[ current, candidate_indices ] <= maximum_distance ]
      new_neighbors <- setdiff( neighbors, component )
      queue <- unique( c( queue, new_neighbors ) )
    }
    components[[ length( components ) + 1L ]] <- sort( unique( component ) )
    unvisited <- setdiff( unvisited, component )
  }
  components
}
select_core_component <- function( smoothed_index, distance_matrix, nn_spacing) {
  threshold <- safe_quantile( smoothed_index, CORE_SEED_QUANTILE )
  candidates <- which( is.finite( smoothed_index ) & smoothed_index >= threshold )
  components <- connected_components_from_distance( candidate_indices = candidates, distance_matrix = distance_matrix, maximum_distance = SEED_CONNECTIVITY_RADIUS_NN * nn_spacing )
  if (length(components) == 0L) {
    components <- list( order( smoothed_index, decreasing = TRUE, na.last = NA )[ seq_len( min( MIN_CORE_SPOTS, sum( is.finite( smoothed_index ) ) ) ) ] )
  }
  component_table <- data.table( component_index = seq_along( components ), n_spots = vapply( components, length, integer(1) ), mean_index = vapply( components, function(indices) {
          safe_mean( smoothed_index[ indices ] )
        },
        numeric(1)
      ), sum_positive_index = vapply(
        components, function(indices) {
          sum( pmax( 0, smoothed_index[ indices ] ), na.rm = TRUE )
        },
        numeric(1)
      )
  )
  data.table::setorder( component_table, -n_spots, -sum_positive_index, -mean_index )
  selected <- components[[ component_table$component_index[ 1L ] ]]
  if (length(selected) < MIN_CORE_SPOTS) {
    maximum_index_spot <- which.max( smoothed_index )
    ordered_near_max <- order( distance_matrix[ maximum_index_spot, ], na.last = NA )
    selected <- head( ordered_near_max, min( MIN_CORE_SPOTS, length( ordered_near_max ) ) )
  }
  sort( unique( selected ) )
}
assign_lesion_zones <- function( core_indices, distance_matrix, nn_spacing) {
  n <- nrow( distance_matrix )
  distance_to_core <- apply( distance_matrix[ , core_indices, drop = FALSE ], 1L, min, na.rm = TRUE )
  distance_to_core_nn <- distance_to_core / nn_spacing
  zone <- rep( NA_character_, n )
  zone[ core_indices ] <- "core"
  non_core <- setdiff( seq_len(n), core_indices )
  zone[ non_core[ distance_to_core_nn[ non_core ] <= BORDER_MAX_DISTANCE_NN ] ] <- "border"
  zone[ non_core[ distance_to_core_nn[ non_core ] > BORDER_MAX_DISTANCE_NN & distance_to_core_nn[ non_core ] <= PERILESIONAL_MAX_DISTANCE_NN ] ] <- "perilesional"
  zone[ non_core[ distance_to_core_nn[ non_core ] > PERILESIONAL_MAX_DISTANCE_NN ] ] <- "remote"
  zone_counts <- table( factor( zone, levels = c( "core", "border", "perilesional", "remote" ) ) )
  fixed_thresholds_pass <- all( zone_counts >= MIN_ZONE_SPOTS )
  zone_method <- "fixed_NN_thresholds"
  border_threshold <- BORDER_MAX_DISTANCE_NN
  perilesional_threshold <- PERILESIONAL_MAX_DISTANCE_NN
  if (!fixed_thresholds_pass) {
    non_core_distance <- distance_to_core_nn[ non_core ]
    border_threshold <- safe_quantile( non_core_distance, FALLBACK_BORDER_DISTANCE_QUANTILE )
    perilesional_threshold <- safe_quantile( non_core_distance, FALLBACK_PERILESIONAL_DISTANCE_QUANTILE )
    zone[non_core] <- ifelse( distance_to_core_nn[non_core] <= border_threshold, "border", ifelse( distance_to_core_nn[non_core] <= perilesional_threshold, "perilesional", "remote" ) )
    zone_method <- "prespecified_distance_quantile_fallback"
  }
  list( zone = factor( zone, levels = c( "core", "border", "perilesional", "remote" ) ), distance_to_core = distance_to_core, distance_to_core_nn = distance_to_core_nn, zone_method = zone_method, border_threshold_nn = border_threshold, perilesional_threshold_nn = perilesional_threshold )
}
build_section_zone_assignment <- function( section_table, normalized_weights) {
  section_id <- as.character( section_table$section_id[ 1L ] )
  timepoint <- as.character( section_table$timepoint[ 1L ] )
  barcodes <- as.character( section_table$barcode )
  if (!all( barcodes %in% rownames( normalized_weights ) )) {
    stop( "RCTD weights are missing section barcodes for ", section_id, ".", call. = FALSE )
  }
  section_weights <- normalized_weights[ barcodes, , drop = FALSE ]
  positive_z <- vapply( STRUCTURAL_POSITIVE_TYPES, function(cell_type) {
      robust_z( section_weights[ , cell_type ] )
    },
    numeric( nrow( section_weights ) )
  )
  negative_z <- vapply( STRUCTURAL_NEGATIVE_TYPES, function(cell_type) {
      robust_z( section_weights[ , cell_type ] )
    },
    numeric( nrow( section_weights ) )
  )
  if (is.null( dim( positive_z ) )) {
    positive_z <- matrix( positive_z, ncol = 1L )
  }
  if (is.null( dim( negative_z ) )) {
    negative_z <- matrix( negative_z, ncol = 1L )
  }
  structural_index_raw <- rowMeans( positive_z ) - rowMeans( negative_z )
  distance_matrix <- pairwise_distance_matrix( x = section_table$coordinate_x_fullres, y = section_table$coordinate_y_fullres )
  nn_spacing <- local_nn_spacing( distance_matrix )
  structural_index_smoothed <- smooth_index_knn( structural_index = structural_index_raw, distance_matrix = distance_matrix, k = K_SMOOTH )
  if (timepoint == UNINJURED_TIMEPOINT) {
    return( data.table( barcode = barcodes, structural_lesion_index_raw = structural_index_raw, structural_lesion_index_smoothed = structural_index_smoothed, local_nn_spacing_fullres = nn_spacing, lesion_core_seed = FALSE, distance_to_core_fullres = NA_real_, distance_to_core_NN = NA_real_, lesion_zone = factor( rep( "uninjured_reference", length( barcodes ) ), levels = c( "core", "border", "perilesional", "remote", "uninjured_reference" ) ), lesion_zone_method = "uninjured_no_artificial_lesion_zone", border_threshold_NN = NA_real_, perilesional_threshold_NN = NA_real_ ) )
  }
  core_indices <- select_core_component( smoothed_index = structural_index_smoothed, distance_matrix = distance_matrix, nn_spacing = nn_spacing )
  zone_result <- assign_lesion_zones( core_indices = core_indices, distance_matrix = distance_matrix, nn_spacing = nn_spacing )
  output <- data.table( barcode = barcodes, structural_lesion_index_raw = structural_index_raw, structural_lesion_index_smoothed = structural_index_smoothed, local_nn_spacing_fullres = nn_spacing, lesion_core_seed = seq_along( barcodes ) %in% core_indices, distance_to_core_fullres = zone_result$distance_to_core, distance_to_core_NN = zone_result$distance_to_core_nn, lesion_zone = factor( as.character( zone_result$zone ), levels = c( "core", "border", "perilesional", "remote", "uninjured_reference" ) ), lesion_zone_method = zone_result$zone_method, border_threshold_NN = zone_result$border_threshold_nn, perilesional_threshold_NN = zone_result$perilesional_threshold_nn )
  output
}
resolve_hcl_palette <- function( requested, fallback = "Viridis") {
  palettes <- grDevices::hcl.pals()
  matched <- palettes[ tolower(palettes) == tolower(requested) ]
  if (length(matched) > 0L) {
    return( matched[[1L]] )
  }
  fallback_match <- palettes[ tolower(palettes) == tolower(fallback) ]
  if (length(fallback_match) > 0L) {
    return( fallback_match[[1L]] )
  }
  palettes[[1L]]
}
continuous_colors <- function( values, palette = "Viridis") {
  values <- as.numeric( values )
  finite_values <- values[ is.finite( values ) ]
  output <- rep( "#808080AA", length( values ) )
  if (length( finite_values ) == 0L) {
    return( output )
  }
  lower <- safe_quantile( finite_values, PLOT_QUANTILE_LOW )
  upper <- safe_quantile( finite_values, PLOT_QUANTILE_HIGH )
  if ( !is.finite(lower) || !is.finite(upper) || upper <= lower ) {
    lower <- min( finite_values )
    upper <- max( finite_values )
  }
  if ( !is.finite(lower) || !is.finite(upper) || upper <= lower ) {
    output[ is.finite( values ) ] <- "#3B528BAA"
    return( output )
  }
  scaled <- ( pmin( upper, pmax( lower, values ) ) - lower ) / ( upper - lower )
  palette_values <- grDevices::hcl.colors( 100L, palette = resolve_hcl_palette( palette ) )
  indices <- floor( scaled * 99 ) + 1L
  valid <- is.finite( indices )
  output[valid] <- grDevices::adjustcolor( palette_values[ pmax( 1L, pmin( 100L, indices[valid] ) ) ], alpha.f = 0.82 )
  output
}
zone_colors <- function(zone) {
  zone <- as.character( zone )
  palette <- c( core = "#B2182B", border = "#EF8A62", perilesional = "#67A9CF", remote = "#2166AC", uninjured_reference = "#7F7F7F" )
  output <- unname( palette[ zone ] )
  output[ is.na( output ) ] <- "#BDBDBD"
  grDevices::adjustcolor( output, alpha.f = 0.82 )
}
draw_image_panel <- function( image_array, image_width, image_height, x, y, colors, title, subtitle = NULL) {
  graphics::plot.new()
  graphics::plot.window( xlim = c( 0, image_width ), ylim = c( image_height, 0 ), asp = 1 )
  graphics::rasterImage( image_array, xleft = 0, ybottom = image_height, xright = image_width, ytop = 0, interpolate = FALSE )
  graphics::points( x, y, pch = 16, cex = POINT_CEX, col = colors )
  graphics::box( col = "#666666" )
  graphics::title( main = title, cex.main = 0.90 )
  if ( !is.null( subtitle ) && nzchar( subtitle ) ) {
    graphics::mtext( subtitle, side = 1, line = -1.1, cex = 0.60 )
  }
}
checks <- data.table( check_id = character(), category = character(), description = character(), status = character(), observed = character(), required = character(), critical = logical(), action_if_not_pass = character() )
add_check <- function( check_id, category, description, passed, observed, required, critical = TRUE, action_if_not_pass = "") {
  checks <<- rbind( checks, data.table( check_id = as.character( check_id ), category = as.character( category ), description = as.character( description ), status = if ( isTRUE( passed ) ) {
        "PASS"
      } else {
        "FAIL"
      },
      observed = as.character( observed ), required = as.character( required ), critical = isTRUE( critical ), action_if_not_pass = as.character( action_if_not_pass )
    ), fill = TRUE
  )
}
add_review <- function( check_id, category, description, observed, required, action_if_not_pass) {
  checks <<- rbind( checks, data.table( check_id = as.character( check_id ), category = as.character( category ), description = as.character( description ), status = "REVIEW", observed = as.character( observed ), required = as.character( required ), critical = FALSE, action_if_not_pass = as.character( action_if_not_pass ) ), fill = TRUE )
}
write_fatal_report <- function(reason) {
  write_run_log()
}
main <- function() {
  workflow_start <- Sys.time()
  required_inputs <- c( Step12E_report = STEP12E_REPORT, Step12E_spot_scores_RDS = STEP12E_SPOT_SCORES_RDS, Step12D_report = STEP12D_REPORT, Step12D_normalized_weights_RDS = STEP12D_NORMALIZED_WEIGHTS_RDS, Step12B_report = STEP12B_REPORT, Step12B_capture_audit = STEP12B_CAPTURE_AUDIT )
  input_exists <- file.exists( required_inputs )
  add_check( "C01", "input", "All frozen Step 12B/12D/12E inputs exist.", all( input_exists ), paste( paste0( names( required_inputs ), "=", input_exists ), collapse = "; " ), "All TRUE", TRUE, "Restore frozen upstream outputs; reconstruct scores or weights." )
  if (!all( input_exists )) {
    stop( "At least one required Step 12F input is missing.", call. = FALSE )
  }
  read_status <- function( path, pattern = "^FINAL_STATUS:") {
    lines <- readLines( path, warn = FALSE, encoding = "UTF-8" )
    lines[ grepl( pattern, lines ) ]
  }
  step12B_status <- read_status( STEP12B_REPORT )
  step12D_status <- read_status( STEP12D_REPORT )
  step12E_status <- read_status( STEP12E_REPORT )
  upstream_pass <- length( step12B_status ) == 1L && length( step12D_status ) == 1L && length( step12E_status ) == 1L && grepl( "^FINAL_STATUS: (PASS|PASS_WITH_REVIEW)$", step12B_status ) && grepl( "^FINAL_STATUS: (PASS|PASS_WITH_REVIEW)$", step12D_status ) && grepl( "^FINAL_STATUS: (PASS|PASS_WITH_REVIEW)$", step12E_status )
  add_check( "C02", "upstream", "Step 12B, Step 12D, and Step 12E release lesion-zone validation.", upstream_pass, paste0( "Step12B=", collapse_unique( step12B_status ), "; Step12D=", collapse_unique( step12D_status ), "; Step12E=", collapse_unique( step12E_status ) ), "All PASS/PASS_WITH_REVIEW", TRUE, "define lesion zones before upstream spatial outputs pass." )
  if (!upstream_pass) {
    stop( "An upstream frozen stage does not release Step 12F.", call. = FALSE )
  }
  spot_scores <- readRDS( STEP12E_SPOT_SCORES_RDS )
  spot_scores <- as.data.table( spot_scores )
  required_columns( spot_scores, c( "barcode", "capture_area_id", "section_id", "timepoint", "replicate_index", "biological_replicate_id", "injury_condition", "qc_review_flag", "Microglia_RCTD_weight", "target_program_AUCell_raw", "shared_prebranch_AUCell_raw", "target_program_microglia_adjusted_z", "target_program_adjusted_z_within_section", "shared_prebranch_microglia_adjusted_z", "shared_prebranch_adjusted_z_within_section" ), "Step 12E spot-score table" )
  normalized_weights <- readRDS( STEP12D_NORMALIZED_WEIGHTS_RDS )
  normalized_weights <- as.matrix( normalized_weights )
  required_structural_types <- c( STRUCTURAL_POSITIVE_TYPES, STRUCTURAL_NEGATIVE_TYPES )
  structural_types_present <- all( required_structural_types %in% colnames( normalized_weights ) )
  alignment_ok <- nrow( spot_scores ) == EXPECTED_SPOTS && nrow( normalized_weights ) == EXPECTED_SPOTS && !anyDuplicated( spot_scores$barcode ) && !anyDuplicated( rownames( normalized_weights ) ) && setequal( spot_scores$barcode, rownames( normalized_weights ) ) && structural_types_present
  add_check( "C03", "alignment", "Spot scores and normalized RCTD weights align and contain all structural cell types.", alignment_ok, paste0( "score rows=", nrow( spot_scores ), "; weight rows=", nrow( normalized_weights ), "; barcode sets equal=", setequal( spot_scores$barcode, rownames( normalized_weights ) ), "; structural types present=", structural_types_present ), "33,941 aligned spots and all prespecified structural cell types", TRUE, "define lesions from incomplete or misaligned weights." )
  if (!alignment_ok) {
    stop( "Step 12E score / Step 12D weight alignment failed.", call. = FALSE )
  }
  spot_scores <- spot_scores[ match( rownames( normalized_weights ), barcode ) ]
  if ( anyNA( spot_scores$barcode ) || !identical( spot_scores$barcode, rownames( normalized_weights ) ) ) {
    stop( "Spot scores could not be reordered to normalized weights.", call. = FALSE )
  }
  capture_areas <- sort( unique( as.character( spot_scores$capture_area_id ) ) )
  coordinate_list <- lapply( capture_areas, function(capture_area) {
      path <- file.path( STEP12D_CAPTURE_ROOT, capture_area, paste0( capture_area, "_spot_metadata.csv.gz" ) )
      if (!file.exists( path )) {
        stop( "Missing Step 12D coordinate metadata for ", capture_area, ".", call. = FALSE )
      }
      table <- data.table::fread( path, data.table = TRUE, showProgress = FALSE )
      required_columns( table, c( "barcode", "coordinate_x_fullres", "coordinate_y_fullres" ), paste0( capture_area, " coordinate table" ) )
      table[ , .( barcode, coordinate_x_fullres, coordinate_y_fullres ) ]
    }
  )
  coordinates <- rbindlist( coordinate_list, fill = TRUE )
  coordinates <- coordinates[ match( spot_scores$barcode, barcode ) ]
  coordinate_alignment_ok <- nrow( coordinates ) == EXPECTED_SPOTS && !anyNA( coordinates$barcode ) && identical( coordinates$barcode, spot_scores$barcode ) && all( is.finite( coordinates$coordinate_x_fullres ) ) && all( is.finite( coordinates$coordinate_y_fullres ) )
  add_check( "C04", "coordinates", "All spot scores have aligned finite full-resolution image coordinates.", coordinate_alignment_ok, paste0( "rows=", nrow( coordinates ), "; missing barcode=", sum( is.na( coordinates$barcode ) ), "; finite x=", sum( is.finite( coordinates$coordinate_x_fullres ) ), "; finite y=", sum( is.finite( coordinates$coordinate_y_fullres ) ) ), "33,941 aligned finite coordinates", TRUE, "calculate lesion distance with missing image coordinates." )
  if (!coordinate_alignment_ok) {
    stop( "Coordinate alignment failed.", call. = FALSE )
  }
  spot_scores[ , coordinate_x_fullres := coordinates$coordinate_x_fullres ]
  spot_scores[ , coordinate_y_fullres := coordinates$coordinate_y_fullres ]
  spatial_structure_ok <- uniqueN( spot_scores$capture_area_id ) == EXPECTED_CAPTURE_AREAS && uniqueN( spot_scores$section_id ) == EXPECTED_SECTIONS && uniqueN( spot_scores[ timepoint %in% INJURED_TIMEPOINTS, section_id ] ) == EXPECTED_INJURED_SECTIONS && all( vapply( INJURED_TIMEPOINTS, function(current_timepoint) {
          uniqueN( spot_scores[ timepoint == current_timepoint, biological_replicate_id ] ) == EXPECTED_ANIMALS_PER_INJURED_TIMEPOINT
        },
        logical(1)
      )
    )
  add_check( "C05", "design", "Spatial metadata contains 36 sections, 24 injured sections, and three animals per injured timepoint.", spatial_structure_ok, paste0( "capture areas=", uniqueN( spot_scores$capture_area_id ), "; sections=", uniqueN( spot_scores$section_id ), "; injured sections=", uniqueN( spot_scores[ timepoint %in% INJURED_TIMEPOINTS, section_id ] ), "; 7days animals=", uniqueN( spot_scores[ timepoint == "7days", biological_replicate_id ] ), "; 2months animals=", uniqueN( spot_scores[ timepoint == "2months", biological_replicate_id ] ) ), "12 capture areas, 36 sections, 24 injured sections, 3 animals/timepoint", TRUE, "use spots as independent replicates." )
  if (!spatial_structure_ok) {
    stop( "Spatial replicate structure is incomplete.", call. = FALSE )
  }
  forbidden_zone_variables <- c( "Microglia", "target_program_AUCell_raw", "target_program_microglia_adjusted_z", "target_program_adjusted_z_within_section", "shared_prebranch_AUCell_raw", "shared_prebranch_microglia_adjusted_z", "shared_prebranch_adjusted_z_within_section", "nCount_Spatial_BPCells", "nFeature_Spatial_BPCells", "percent.mt" )
  zone_definition_variables <- c( STRUCTURAL_POSITIVE_TYPES, STRUCTURAL_NEGATIVE_TYPES, "coordinate_x_fullres", "coordinate_y_fullres" )
  independence_ok <- length( intersect( forbidden_zone_variables, zone_definition_variables ) ) == 0L
  add_check( "C06", "independence", "Lesion-zone definition excludes target scores, Shared_prebranch, Microglia, and QC metrics.", independence_ok, paste0( "zone variables=", paste( zone_definition_variables, collapse = " | " ), "; forbidden overlap=", paste( intersect( forbidden_zone_variables, zone_definition_variables ), collapse = " | " ) ), "Forbidden overlap = none", TRUE, "use the tested target signal to choose lesion locations." )
  if (!independence_ok) {
    stop( "Lesion-zone independence boundary was violated.", call. = FALSE )
  }
  section_ids <- unique( as.character( spot_scores$section_id ) )
  section_ids <- section_ids[ order( match( section_ids, unique( as.character( spot_scores$section_id ) ) ) ) ]
  zone_assignment_list <- vector( "list", length( section_ids ) )
  zone_audit_list <- vector( "list", length( section_ids ) )
  for ( section_index in seq_along( section_ids ) ) {
    current_section_id <- section_ids[ section_index ]
    section_table <- spot_scores[ section_id == current_section_id ]
    assignment <- build_section_zone_assignment( section_table = section_table, normalized_weights = normalized_weights )
    assignment[ , section_id := current_section_id ]
    zone_assignment_list[[ section_index ]] <- assignment
    counts <- table( factor( as.character( assignment$lesion_zone ), levels = c( "core", "border", "perilesional", "remote", "uninjured_reference" ) ) )
    zone_audit_list[[ section_index ]] <- data.table( section_id = current_section_id, capture_area_id = as.character( section_table$capture_area_id[ 1L ] ), timepoint = as.character( section_table$timepoint[ 1L ] ), biological_replicate_id = as.character( section_table$biological_replicate_id[ 1L ] ), n_spots = nrow( assignment ), core_spots = unname( counts[ "core" ] ), border_spots = unname( counts[ "border" ] ), perilesional_spots = unname( counts[ "perilesional" ] ), remote_spots = unname( counts[ "remote" ] ), uninjured_reference_spots = unname( counts[ "uninjured_reference" ] ), zone_method = collapse_unique( assignment$lesion_zone_method ), local_nn_spacing_fullres = safe_median( assignment$local_nn_spacing_fullres ), border_threshold_NN = safe_median( assignment$border_threshold_NN ), perilesional_threshold_NN = safe_median( assignment$perilesional_threshold_NN ), structural_index_min = min( assignment$structural_lesion_index_smoothed, na.rm = TRUE ), structural_index_median = stats::median( assignment$structural_lesion_index_smoothed, na.rm = TRUE ), structural_index_max = max( assignment$structural_lesion_index_smoothed, na.rm = TRUE ) )
    rm( assignment, section_table )
    invisible( gc() )
  }
  zone_assignment <- rbindlist( zone_assignment_list, fill = TRUE )
  zone_audit <- rbindlist( zone_audit_list, fill = TRUE )
  zone_assignment <- zone_assignment[ match( spot_scores$barcode, barcode ) ]
  zone_alignment_ok <- nrow( zone_assignment ) == EXPECTED_SPOTS && !anyNA( zone_assignment$barcode ) && !anyDuplicated( zone_assignment$barcode ) && identical( zone_assignment$barcode, spot_scores$barcode )
  injured_zone_counts_ok <- all( zone_audit[ timepoint %in% INJURED_TIMEPOINTS, core_spots >= MIN_CORE_SPOTS & border_spots >= MIN_ZONE_SPOTS & perilesional_spots >= MIN_ZONE_SPOTS & remote_spots >= MIN_ZONE_SPOTS ] )
  uninjured_zone_ok <- all( zone_audit[ timepoint == UNINJURED_TIMEPOINT, uninjured_reference_spots == n_spots & core_spots == 0L & border_spots == 0L & perilesional_spots == 0L & remote_spots == 0L ] )
  add_check( "C07", "zone_assignment", "All spots receive exactly one section-aware zone, with adequate injured-zone sizes and no artificial uninjured lesion.", zone_alignment_ok && injured_zone_counts_ok && uninjured_zone_ok, paste0( "aligned spots=", nrow( zone_assignment ), "; injured sections with adequate zones=", sum( zone_audit[ timepoint %in% INJURED_TIMEPOINTS, core_spots >= MIN_CORE_SPOTS & border_spots >= MIN_ZONE_SPOTS & perilesional_spots >= MIN_ZONE_SPOTS & remote_spots >= MIN_ZONE_SPOTS ] ), "/", EXPECTED_INJURED_SECTIONS, "; uninjured sections valid=", sum( zone_audit[ timepoint == UNINJURED_TIMEPOINT, uninjured_reference_spots == n_spots ] ), "/12" ), "33,941 aligned spots; 24 injured sections with four adequate zones; 12 uninjured reference sections", TRUE, "Review section maps before any lesion-proximity interpretation." )
  if ( !zone_alignment_ok || !injured_zone_counts_ok || !uninjured_zone_ok ) {
    stop( "Independent lesion-zone assignment failed.", call. = FALSE )
  }
  spot_zones <- cbind( spot_scores, zone_assignment[ , !c( "barcode", "section_id" ) ] )
  injured_spots <- spot_zones[ timepoint %in% INJURED_TIMEPOINTS ]
  zone_levels <- c( "core", "border", "perilesional", "remote" )
  section_zone_summary <- injured_spots[ , .( n_spots = .N, n_QC_review_spots = sum( qc_review_flag ), structural_index_mean = safe_mean( structural_lesion_index_smoothed ), structural_index_median = safe_median( structural_lesion_index_smoothed ), distance_to_core_NN_mean = safe_mean( distance_to_core_NN ), distance_to_core_NN_median = safe_median( distance_to_core_NN ), Microglia_weight_mean = safe_mean( Microglia_RCTD_weight ), Microglia_weight_median = safe_median( Microglia_RCTD_weight ), target_raw_mean = safe_mean( target_program_AUCell_raw ), target_raw_median = safe_median( target_program_AUCell_raw ), target_adjusted_global_z_mean = safe_mean( target_program_microglia_adjusted_z ), target_adjusted_global_z_median = safe_median( target_program_microglia_adjusted_z ), target_adjusted_within_section_z_mean = safe_mean( target_program_adjusted_z_within_section ), target_adjusted_within_section_z_median = safe_median( target_program_adjusted_z_within_section ), shared_raw_mean = safe_mean( shared_prebranch_AUCell_raw ), shared_raw_median = safe_median( shared_prebranch_AUCell_raw ), shared_adjusted_within_section_z_mean = safe_mean( shared_prebranch_adjusted_z_within_section ), shared_adjusted_within_section_z_median = safe_median( shared_prebranch_adjusted_z_within_section ) ), by = .( capture_area_id, section_id, timepoint, replicate_index, biological_replicate_id, injury_condition, lesion_zone ) ]
  section_zone_summary[ , lesion_zone := factor( as.character( lesion_zone ), levels = zone_levels ) ]
  metrics_to_contrast <- c( target_adjusted_within_section_z = "target_program_adjusted_z_within_section", target_adjusted_global_z = "target_program_microglia_adjusted_z", target_raw_AUCell = "target_program_AUCell_raw", shared_adjusted_within_section_z = "shared_prebranch_adjusted_z_within_section", shared_raw_AUCell = "shared_prebranch_AUCell_raw", Microglia_RCTD_weight = "Microglia_RCTD_weight", structural_lesion_index = "structural_lesion_index_smoothed" )
  contrast_zones <- c( "core", "border", "perilesional" )
  section_contrast_rows <- list()
  row_index <- 0L
  injured_section_ids <- unique( as.character( injured_spots$section_id ) )
  for ( current_section in injured_section_ids ) {
    section_data <- injured_spots[ section_id == current_section ]
    remote_data <- section_data[ lesion_zone == "remote" ]
    if (nrow( remote_data ) < MIN_ZONE_SPOTS) {
      stop( "Remote zone is undersized in ", current_section, ".", call. = FALSE )
    }
    for ( metric_name in names( metrics_to_contrast ) ) {
      metric_column <- unname( metrics_to_contrast[ metric_name ] )
      remote_value <- safe_mean( remote_data[[ metric_column ]] )
      for ( contrast_zone in contrast_zones ) {
        zone_data <- section_data[ lesion_zone == contrast_zone ]
        zone_value <- safe_mean( zone_data[[ metric_column ]] )
        row_index <- row_index + 1L
        section_contrast_rows[[ row_index ]] <- data.table( capture_area_id = as.character( section_data$capture_area_id[ 1L ] ), section_id = current_section, timepoint = as.character( section_data$timepoint[ 1L ] ), replicate_index = section_data$replicate_index[ 1L ], biological_replicate_id = as.character( section_data$biological_replicate_id[ 1L ] ), metric = metric_name, contrast = paste0( contrast_zone, "_minus_remote" ), zone = contrast_zone, reference_zone = "remote", n_zone_spots = nrow( zone_data ), n_remote_spots = nrow( remote_data ), zone_mean = zone_value, remote_mean = remote_value, zone_minus_remote = zone_value - remote_value )
      }
    }
  }
  section_contrasts <- rbindlist( section_contrast_rows, fill = TRUE )
  expected_section_contrast_rows <- EXPECTED_INJURED_SECTIONS * length( metrics_to_contrast ) * length( contrast_zones )
  section_outputs_ok <- nrow( section_zone_summary ) == EXPECTED_INJURED_SECTIONS * 4L && nrow( section_contrasts ) == expected_section_contrast_rows && all( is.finite( section_contrasts$zone_minus_remote ) )
  add_check( "C08", "section_contrasts", "Section summaries and prespecified zone-vs-remote contrasts are complete.", section_outputs_ok, paste0( "zone-summary rows=", nrow( section_zone_summary ), "; contrast rows=", nrow( section_contrasts ), "/", expected_section_contrast_rows, "; finite contrasts=", sum( is.finite( section_contrasts$zone_minus_remote ) ) ), "96 section-zone rows and all expected finite contrasts", TRUE, "aggregate incomplete section contrasts to animals." )
  if (!section_outputs_ok) {
    stop( "Section-level zone summary or contrast validation failed.", call. = FALSE )
  }
  animal_contrasts <- section_contrasts[ , .( n_sections = uniqueN( section_id ), animal_mean_contrast = safe_mean( zone_minus_remote ), animal_median_contrast = safe_median( zone_minus_remote ), animal_min_section_contrast = min( zone_minus_remote ), animal_max_section_contrast = max( zone_minus_remote ), positive_sections = sum( zone_minus_remote > 0 ), negative_sections = sum( zone_minus_remote < 0 ), section_contrasts = paste( round( zone_minus_remote, 6 ), collapse = " | " ) ), by = .( timepoint, biological_replicate_id, metric, contrast, zone, reference_zone ) ]
  animals_per_group_ok <- all( animal_contrasts[ , uniqueN( biological_replicate_id ), by = .( timepoint, metric, contrast ) ]$V1 == EXPECTED_ANIMALS_PER_INJURED_TIMEPOINT )
  add_check( "C09", "animal_unit", "Every timepoint/metric/contrast has three independent animal summaries.", animals_per_group_ok, paste0( "minimum animals/group=", min( animal_contrasts[ , uniqueN( biological_replicate_id ), by = .( timepoint, metric, contrast ) ]$V1 ), "; maximum animals/group=", max( animal_contrasts[ , uniqueN( biological_replicate_id ), by = .( timepoint, metric, contrast ) ]$V1 ) ), "Exactly three animals per injured timepoint", TRUE, "replace animal-level inference with spot-level testing." )
  if (!animals_per_group_ok) {
    stop( "Animal-level replicate audit failed.", call. = FALSE )
  }
  animal_inference <- animal_contrasts[ , .( n_animals = .N, animal_mean_contrast_mean = safe_mean( animal_mean_contrast ), animal_mean_contrast_median = safe_median( animal_mean_contrast ), animal_mean_contrast_sd = safe_sd( animal_mean_contrast ), animal_mean_contrast_min = min( animal_mean_contrast ), animal_mean_contrast_max = max( animal_mean_contrast ), positive_animals = sum( animal_mean_contrast > 0 ), negative_animals = sum( animal_mean_contrast < 0 ), zero_animals = sum( animal_mean_contrast == 0 ), direction_consistency_fraction = max( sum( animal_mean_contrast > 0 ), sum( animal_mean_contrast < 0 ) ) / .N, exact_sign_test_p = exact_sign_test_p( animal_mean_contrast ), wilcoxon_signed_rank_p = safe_wilcoxon_signed_p( animal_mean_contrast ), inference_unit = "biological_animal", spot_level_p_values_used = FALSE ), by = .( timepoint, metric, contrast, zone, reference_zone ) ]
  between_timepoint_rows <- animal_contrasts[ , {
      seven_day <- animal_mean_contrast[ timepoint == "7days" ]
      two_month <- animal_mean_contrast[ timepoint == "2months" ]
      list( timepoint = "7days_vs_2months", n_animals = length( seven_day ) + length( two_month ), animal_mean_contrast_mean = safe_mean( two_month ) - safe_mean( seven_day ), animal_mean_contrast_median = safe_median( two_month ) - safe_median( seven_day ), animal_mean_contrast_sd = NA_real_, animal_mean_contrast_min = NA_real_, animal_mean_contrast_max = NA_real_, positive_animals = NA_integer_, negative_animals = NA_integer_, zero_animals = NA_integer_, direction_consistency_fraction = NA_real_, exact_sign_test_p = NA_real_, wilcoxon_signed_rank_p = safe_wilcoxon_rank_p( seven_day, two_month ), inference_unit = "biological_animal_between_timepoints", spot_level_p_values_used = FALSE )
    },
    by = .( metric, contrast, zone, reference_zone )
  ]
  animal_inference <- rbind( animal_inference, between_timepoint_rows, fill = TRUE )
  primary_animal_results <- animal_inference[ metric == "target_adjusted_within_section_z" & timepoint %in% INJURED_TIMEPOINTS ]
  primary_complete <- nrow( primary_animal_results ) == length( INJURED_TIMEPOINTS ) * length( contrast_zones ) && all( primary_animal_results$n_animals == EXPECTED_ANIMALS_PER_INJURED_TIMEPOINT )
  add_check( "C10", "primary_inference", "Primary adjusted-target animal summaries are complete for both injured timepoints and all lesion-proximity contrasts.", primary_complete, paste0( "rows=", nrow( primary_animal_results ), "; expected=", length( INJURED_TIMEPOINTS ) * length( contrast_zones ), "; animals/row=", paste( primary_animal_results$n_animals, collapse = " | " ) ), "Six primary rows, each with three animals", TRUE, "make target-localization claims from incomplete animal summaries." )
  if (!primary_complete) {
    stop( "Primary animal-level inference summary is incomplete.", call. = FALSE )
  }
  injured_spots[ , distance_bin := cut( distance_to_core_NN, breaks = DISTANCE_PROFILE_BREAKS_NN, include.lowest = TRUE, right = FALSE ) ]
  distance_profile <- injured_spots[ , .( n_spots = .N, distance_to_core_NN_mean = safe_mean( distance_to_core_NN ), distance_to_core_NN_median = safe_median( distance_to_core_NN ), target_adjusted_within_section_z_mean = safe_mean( target_program_adjusted_z_within_section ), target_adjusted_within_section_z_median = safe_median( target_program_adjusted_z_within_section ), target_raw_mean = safe_mean( target_program_AUCell_raw ), shared_adjusted_within_section_z_mean = safe_mean( shared_prebranch_adjusted_z_within_section ), Microglia_weight_mean = safe_mean( Microglia_RCTD_weight ), structural_index_mean = safe_mean( structural_lesion_index_smoothed ) ), by = .( capture_area_id, section_id, timepoint, biological_replicate_id, distance_bin ) ]
  distance_profile_ok <- nrow( distance_profile ) > 0L && all( distance_profile$n_spots > 0 ) && all( is.finite( distance_profile$ target_adjusted_within_section_z_mean ) )
  add_check( "C11", "distance_profile", "Section-level distance-to-core profiles are nonempty and finite.", distance_profile_ok, paste0( "rows=", nrow( distance_profile ), "; sections=", uniqueN( distance_profile$section_id ), "; finite target means=", sum( is.finite( distance_profile$ target_adjusted_within_section_z_mean ) ) ), "All injured sections represented with finite target profiles", TRUE, "interpret an incomplete radial profile." )
  if (!distance_profile_ok) {
    stop( "Distance-to-core profile validation failed.", call. = FALSE )
  }
  input_audit <- rbindlist( list( data.table( item_type = "package", item = required_packages, value = vapply( required_packages, function(package_name) {
              as.character( utils::packageVersion( package_name ) )
            },
            FUN.VALUE = character(1)
          ), note = "Version used in Step 12F"
      ), data.table(
        item_type = "zone_definition", item = c( "positive_structural_cell_types", "negative_structural_cell_types", "explicitly_excluded_variables", "K_SMOOTH", "CORE_SEED_QUANTILE", "SEED_CONNECTIVITY_RADIUS_NN", "MIN_CORE_SPOTS", "MIN_ZONE_SPOTS", "BORDER_MAX_DISTANCE_NN", "PERILESIONAL_MAX_DISTANCE_NN", "fallback_border_distance_quantile", "fallback_perilesional_distance_quantile", "primary_outcome", "inferential_unit" ), value = c( paste( STRUCTURAL_POSITIVE_TYPES, collapse = " | " ), paste( STRUCTURAL_NEGATIVE_TYPES, collapse = " | " ), paste( forbidden_zone_variables, collapse = " | " ), K_SMOOTH, CORE_SEED_QUANTILE, SEED_CONNECTIVITY_RADIUS_NN, MIN_CORE_SPOTS, MIN_ZONE_SPOTS, BORDER_MAX_DISTANCE_NN, PERILESIONAL_MAX_DISTANCE_NN, FALLBACK_BORDER_DISTANCE_QUANTILE, FALLBACK_PERILESIONAL_DISTANCE_QUANTILE, "target_program_adjusted_z_within_section", "biological_animal" ), note = c( "Within-section robust Z scores", "Subtracted from positive structural mean", "Never used to locate lesion zones", "Spatial smoothing neighbors", "Top smoothed-index seed fraction", "Seed graph radius in local NN units", "Minimum connected lesion core", "Minimum spot count per injured zone", "Fixed border threshold", "Fixed perilesional threshold", "Prespecified fallback only", "Prespecified fallback only", "Microglia-adjusted target localization", "Section contrasts averaged within animal" )
      ), zone_audit[
        , .( item_type = "section_zone_audit", item = section_id, value = paste0( "core=", core_spots, "; border=", border_spots, "; perilesional=", perilesional_spots, "; remote=", remote_spots, "; uninjured=", uninjured_reference_spots ), note = paste0( "timepoint=", timepoint, "; method=", zone_method, "; NN_spacing=", round( local_nn_spacing_fullres, 4 ) ) )
      ]
    ), fill = TRUE
  )
  data.table::fwrite( input_audit, OUT_INPUT_AUDIT, bom = TRUE )
  data.table::fwrite( spot_zones, OUT_SPOT_ZONES, bom = TRUE, compress = "gzip" )
  data.table::fwrite( section_zone_summary, OUT_SECTION_ZONE_SUMMARY, bom = TRUE )
  data.table::fwrite( section_contrasts, OUT_SECTION_CONTRASTS, bom = TRUE )
  data.table::fwrite( animal_contrasts, OUT_ANIMAL_CONTRASTS, bom = TRUE )
  if (!file.exists(ORIGINAL_STEP12F_ANIMAL_CSV)) {
    stop("Frozen original Step12F animal CSV missing: ", ORIGINAL_STEP12F_ANIMAL_CSV, call. = FALSE)
  }
  if (!file.exists(STEP12H_ANIMAL_CSV)) {
    stop("Independent Step12H check table missing: ", STEP12H_ANIMAL_CSV, call. = FALSE)
  }
  old_animal <- data.table::fread(ORIGINAL_STEP12F_ANIMAL_CSV, showProgress=FALSE)
  revised_animal <- data.table::copy(animal_contrasts)
  stopifnot(nrow(old_animal)==126L, nrow(revised_animal)==126L, ncol(old_animal)==14L, identical(names(old_animal),names(revised_animal)), sum(old_animal$metric=="shared_adjusted_z")==18L, sum(revised_animal$metric=="shared_adjusted_within_section_z")==18L, sum(revised_animal$metric=="shared_adjusted_z")==0L)
  key_cols <- c("timepoint", "biological_replicate_id", "metric", "contrast", "zone", "reference_zone")
  old_nonshared <- old_animal[metric != "shared_adjusted_z"]
  new_nonshared <- revised_animal[metric != "shared_adjusted_within_section_z"]
  stopifnot(nrow(old_nonshared)==108L, nrow(new_nonshared)==108L)
  data.table::setorderv(old_nonshared,key_cols)
  data.table::setorderv(new_nonshared,key_cols)
  if (!identical(old_nonshared[, ..key_cols],new_nonshared[, ..key_cols])) {
    stop("Frozen NON-shared animal contrast keys changed unexpectedly.", call. = FALSE)
  }
  numeric_compare <- c("n_sections", "animal_mean_contrast", "animal_median_contrast", "animal_min_section_contrast", "animal_max_section_contrast", "positive_sections", "negative_sections")
  max_nonshared_difference <- max(vapply(numeric_compare, function(nm) {
    max(abs(as.numeric(old_nonshared[[nm]])-as.numeric(new_nonshared[[nm]])))
  }, numeric(1)))
  if (!is.finite(max_nonshared_difference) || max_nonshared_difference > 1e-8) {
    stop("NON-shared outcome differs from original Step12F beyond tolerance: ", max_nonshared_difference, call. = FALSE)
  }
  h <- data.table::fread(STEP12H_ANIMAL_CSV,showProgress=FALSE)
  h <- h[analysis=="unrestricted_original_shared" & metric=="shared_prebranch_adjusted_z_within_section"]
  new_shared <- revised_animal[metric=="shared_adjusted_within_section_z"]
  shared_keys <- c("timepoint", "biological_replicate_id", "contrast")
  if (nrow(h)!=18L || nrow(new_shared)!=18L) {
    stop("Incorrect shared metric coverage in Step12H or revised Step12F.",call.=FALSE)
  }
  data.table::setorderv(h,shared_keys)
  data.table::setorderv(new_shared,shared_keys)
  if (!identical(h[, ..shared_keys],new_shared[, ..shared_keys])) {
    stop("Revised shared metric and Step12H cross-check keys differ.",call.=FALSE)
  }
  max_shared_difference <- max(abs(h$animal_mean_contrast-new_shared$animal_mean_contrast))
  if (!is.finite(max_shared_difference) || max_shared_difference > 1e-8) {
    stop("New shared score does not reproduce independent Step12H result: ", max_shared_difference,call.=FALSE)
  }
  provenance_audit <- data.table::data.table( check=c("original_step12F_animal_csv_exists","revised_126_rows", "7_metrics_times_6_animals_times_3_zones", "108_nonshared_rows_reproduce_original","18_shared_rows_updated_from_step12E", "18_new_shared_rows_agree_with_step12H","original_step12F_CSV_md5","Step12H_CSV_md5"), value=c(ORIGINAL_STEP12F_ANIMAL_CSV,as.character(nrow(revised_animal)), paste(length(unique(revised_animal$metric)), length(unique(revised_animal$biological_replicate_id)), length(unique(revised_animal$contrast)),sep=" x "), format(max_nonshared_difference,digits=12), as.character(nrow(new_shared)), format(max_shared_difference,digits=12), as.character(tools::md5sum(ORIGINAL_STEP12F_ANIMAL_CSV)), as.character(tools::md5sum(STEP12H_ANIMAL_CSV))), pass=c(TRUE,nrow(revised_animal)==126L, length(unique(revised_animal$metric))==7L, max_nonshared_difference<=1e-8, nrow(new_shared)==18L, max_shared_difference<=1e-8,TRUE,TRUE) )
  data.table::fwrite(provenance_audit,OUT_TABLE_S10_AUDIT,bom=TRUE)
  if (!all(provenance_audit$pass)) stop("Table S10 provenance audit failed.",call.=FALSE)
  submitted_metric_labels <- c( target_adjusted_within_section_z = "target_program_adjusted_within_section_z", target_adjusted_global_z = "target_program_adjusted_global_z", target_raw_AUCell = "target_program_raw_AUCell", shared_adjusted_within_section_z = "shared_prebranch_adjusted_within_section_z", shared_raw_AUCell = "shared_prebranch_raw_AUCell", Microglia_RCTD_weight = "Microglia_RCTD_weight", structural_lesion_index = "structural_lesion_index" )
  submitted <- data.table::copy(animal_contrasts)
  stopifnot(all(submitted$metric %in% names(submitted_metric_labels)))
  submitted[, metric := unname(submitted_metric_labels[metric])]
  stopifnot(nrow(submitted)==126L,ncol(submitted)==14L, length(unique(submitted$metric))==7L)
  if (!requireNamespace("openxlsx",quietly=TRUE)) {
    stop("Missing package openxlsx. Install with install.packages('openxlsx').",call.=FALSE)
  }
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb,"Table S10",gridLines=FALSE)
  openxlsx::writeData(wb,"Table S10", "Animal-level lesion-zone versus remote spatial contrasts in GSE234774", startRow=1,startCol=1)
  openxlsx::writeData(wb,"Table S10",as.data.frame(submitted), startRow=2,startCol=1,withFilter=FALSE)
  openxlsx::mergeCells(wb,"Table S10",cols=1:14,rows=1)
  openxlsx::addStyle(wb,"Table S10",openxlsx::createStyle( textDecoration="bold",fgFill="#193851",fontColour="#FFFFFF",wrapText=TRUE), rows=2,cols=1:14,gridExpand=TRUE)
  openxlsx::setColWidths(wb,"Table S10",cols=1:14, widths=c(12,22,43,24,14,19,12,24,24,27,27,18,18,70))
  openxlsx::freezePane(wb,"Table S10",firstActiveRow=3)
  openxlsx::saveWorkbook(wb,OUT_REVISED_TABLE_S10,overwrite=TRUE)
  if (!file.exists(OUT_REVISED_TABLE_S10)) stop("Table S10 XLSX export failed.",call.=FALSE)
  submission_dir <- file.path(PROJECT_ROOT,"GLIAv01")
  dir.create(submission_dir,recursive=TRUE,showWarnings=FALSE)
  submission_copy <- file.path(submission_dir,"Table_S10_R2MC4_REVISED.xlsx")
  if (!file.copy(OUT_REVISED_TABLE_S10,submission_copy,overwrite=TRUE)) {
    stop("Failed to copy revised Table S10 into GLIAv01.",call.=FALSE)
  }
  data.table::fwrite( animal_inference, OUT_ANIMAL_INFERENCE, bom = TRUE )
  data.table::fwrite( distance_profile, OUT_DISTANCE_PROFILE, bom = TRUE )
  saveRDS( spot_zones, OUT_SPOT_ZONES_RDS, compress = FALSE )
  capture_audit <- data.table::fread( STEP12B_CAPTURE_AUDIT, data.table = TRUE, showProgress = FALSE )
  required_columns( capture_audit, c( "capture_area_id", "tissue_lowres_scalef", "asset_linkage_status" ), "Step 12B capture-area image audit" )
  capture_audit[ , capture_area_id := toupper( trimws( as.character( capture_area_id ) ) ) ]
  grDevices::pdf( OUT_ZONE_MAPS, width = 11, height = 8.5, onefile = TRUE, useDingbats = FALSE )
  pdf_open <- TRUE
  on.exit( {
      if (isTRUE( pdf_open )) {
        try( grDevices::dev.off(), silent = TRUE )
      }
    },
    add = TRUE
  )
  map_pass <- logical( length( capture_areas ) )
  for ( area_index in seq_along( capture_areas ) ) {
    capture_area <- capture_areas[ area_index ]
    current <- spot_zones[ capture_area_id == capture_area ]
    image_file <- file.path( LOWRES_PREVIEW_ROOT, capture_area, "tissue_lowres_image.png" )
    audit_row <- capture_audit[ capture_area_id == capture_area ]
    if ( nrow( audit_row ) != 1L || !file.exists( image_file ) ) {
      stop( "Missing image/audit asset for ", capture_area, ".", call. = FALSE )
    }
    image_array <- png::readPNG( image_file, native = FALSE, info = TRUE )
    image_height <- dim( image_array )[1L]
    image_width <- dim( image_array )[2L]
    lowres_scale <- as.numeric( audit_row$tissue_lowres_scalef )
    image_x <- current$coordinate_x_fullres * lowres_scale
    image_y <- current$coordinate_y_fullres * lowres_scale
    alignment_pass <- all( is.finite( image_x ) ) && all( is.finite( image_y ) ) && all( image_x >= 0 & image_x <= image_width + 2 ) && all( image_y >= 0 & image_y <= image_height + 2 )
    map_pass[ area_index ] <- alignment_pass
    if (!alignment_pass) {
      stop( "Image-coordinate alignment failed for ", capture_area, ".", call. = FALSE )
    }
    graphics::par( mfrow = c( 2, 2 ), mar = c( 0.5, 0.5, 2.4, 0.5 ), oma = c( 0.4, 0.4, 2.0, 0.4 ), xaxs = "i", yaxs = "i" )
    draw_image_panel( image_array = image_array, image_width = image_width, image_height = image_height, x = image_x, y = image_y, colors = continuous_colors( current$ structural_lesion_index_smoothed, palette = "Inferno" ), title = "Independent structural lesion index", subtitle = "Target, Shared and Microglia excluded" )
    draw_image_panel( image_array = image_array, image_width = image_width, image_height = image_height, x = image_x, y = image_y, colors = zone_colors( current$lesion_zone ), title = "Independent lesion zones", subtitle = "core / border / perilesional / remote" )
    draw_image_panel( image_array = image_array, image_width = image_width, image_height = image_height, x = image_x, y = image_y, colors = continuous_colors( current$ target_program_adjusted_z_within_section, palette = "Blue-Red 3" ), title = "Primary adjusted target score", subtitle = "Not used for zone definition" )
    draw_image_panel( image_array = image_array, image_width = image_width, image_height = image_height, x = image_x, y = image_y, colors = continuous_colors( current$ Microglia_RCTD_weight, palette = "Viridis" ), title = "Microglia RCTD weight", subtitle = "Not used for zone definition" )
    graphics::mtext( paste0( "GSE234774 Step 12F | ", capture_area, " | independent zone definition" ), side = 3, outer = TRUE, line = 0.35, cex = 1.02, font = 2 )
    rm( image_array )
    invisible( gc() )
  }
  grDevices::dev.off()
  pdf_open <- FALSE
  map_output_ok <- all( map_pass ) && file.exists( OUT_ZONE_MAPS ) && file_size_safe( OUT_ZONE_MAPS ) > 0
  add_check( "C12", "spatial_maps", "All 12 capture areas produced aligned independent-zone validation maps.", map_output_ok, paste0( "aligned capture areas=", sum( map_pass ), "/12; PDF bytes=", file_size_safe( OUT_ZONE_MAPS ) ), "12/12 aligned and nonempty PDF", TRUE, "Visually review the structural core and zone continuity before biological interpretation." )
  if (!map_output_ok) {
    stop( "Independent lesion-zone map generation failed.", call. = FALSE )
  }
  primary_animals <- animal_contrasts[ metric == "target_adjusted_within_section_z" ]
  zone_order <- c( "core", "border", "perilesional" )
  grDevices::pdf( OUT_ANIMAL_CONTRAST_FIGURE, width = 10, height = 6.5, useDingbats = FALSE )
  graphics::par( mar = c( 6, 4.5, 2.5, 1 ) )
  plot_grid <- data.table::CJ( timepoint = INJURED_TIMEPOINTS, zone = zone_order, sorted = FALSE )
  plot_labels <- paste( plot_grid$timepoint, plot_grid$zone, sep = "\n" )
  plot_data <- lapply( seq_len( nrow( plot_grid ) ), function(index) {
      current_timepoint <- plot_grid$timepoint[ index ]
      current_zone <- plot_grid$zone[ index ]
      primary_animals[ timepoint == current_timepoint & zone == current_zone, animal_mean_contrast ]
    }
  )
  graphics::boxplot( plot_data, names = plot_labels, ylab = "Animal mean: adjusted target zone − remote", xlab = "", main = "Independent lesion-zone validation (animal is the replicate)", las = 2, outline = FALSE )
  graphics::abline( h = 0, lty = 2, lwd = 1.2 )
  for ( index in seq_along( plot_data ) ) {
    graphics::points( rep( index, length( plot_data[[ index ]] ) ), plot_data[[ index ]], pch = 16, cex = 1.0 )
  }
  grDevices::dev.off()
  animal_distance_profile <- distance_profile[ , .( distance_midpoint = safe_mean( distance_to_core_NN_mean ), target_adjusted_mean = safe_mean( target_adjusted_within_section_z_mean ) ), by = .( timepoint, biological_replicate_id, distance_bin ) ]
  grDevices::pdf( OUT_DISTANCE_PROFILE_FIGURE, width = 9, height = 6.5, useDingbats = FALSE )
  graphics::plot( NA, xlim = range( animal_distance_profile$distance_midpoint, finite = TRUE ), ylim = range( animal_distance_profile$target_adjusted_mean, finite = TRUE ), xlab = "Distance from independent lesion core (local NN units)", ylab = "Animal-averaged adjusted target score", main = "Target-program profile along distance from independent lesion core" )
  line_types <- c( "7days" = 1, "2months" = 2 )
  point_types <- c( "7days" = 16, "2months" = 17 )
  for ( current_timepoint in INJURED_TIMEPOINTS ) {
    current <- animal_distance_profile[ timepoint == current_timepoint ]
    animals <- unique( as.character( current$biological_replicate_id ) )
    for ( animal in animals ) {
      animal_data <- current[ biological_replicate_id == animal ]
      data.table::setorder( animal_data, distance_midpoint )
      graphics::lines( animal_data$distance_midpoint, animal_data$target_adjusted_mean, lty = line_types[ current_timepoint ], lwd = 1 )
      graphics::points( animal_data$distance_midpoint, animal_data$target_adjusted_mean, pch = point_types[ current_timepoint ], cex = 0.75 )
    }
  }
  graphics::abline( h = 0, lty = 3 )
  graphics::legend( "topright", legend = INJURED_TIMEPOINTS, lty = line_types[ INJURED_TIMEPOINTS ], pch = point_types[ INJURED_TIMEPOINTS ], bty = "n" )
  grDevices::dev.off()
  figure_outputs_ok <- all( file.exists( c( OUT_ANIMAL_CONTRAST_FIGURE, OUT_DISTANCE_PROFILE_FIGURE ) ) ) && all( file.info( c( OUT_ANIMAL_CONTRAST_FIGURE, OUT_DISTANCE_PROFILE_FIGURE ) )$size > 0 )
  add_check( "C13", "summary_figures", "Animal-level contrast and distance-profile figures were written.", figure_outputs_ok, paste0( "contrast PDF bytes=", file_size_safe( OUT_ANIMAL_CONTRAST_FIGURE ), "; distance PDF bytes=", file_size_safe( OUT_DISTANCE_PROFILE_FIGURE ) ), "Both nonempty", TRUE, "report spatial localization without animal-level visual summaries." )
  if (!figure_outputs_ok) {
    stop( "Animal-level summary figure generation failed.", call. = FALSE )
  }
  fallback_sections <- zone_audit[ timepoint %in% INJURED_TIMEPOINTS & zone_method == "prespecified_distance_quantile_fallback", section_id ]
  if (length( fallback_sections ) > 0L) {
    add_review( "R01", "zone_threshold", "Some injured sections required the prespecified distance-quantile fallback to maintain adequate zone sizes.", paste( fallback_sections, collapse = " | " ), "Inspect the 12-page zone-map PDF before interpreting exact border widths.", "Keep the fallback fixed; tune thresholds using target-score results." )
  }
  primary_direction <- animal_inference[ metric == "target_adjusted_within_section_z" & timepoint %in% INJURED_TIMEPOINTS, .( timepoint, contrast, n_animals, animal_mean_contrast_mean, positive_animals, negative_animals, exact_sign_test_p, wilcoxon_signed_rank_p ) ]
  all_primary_positive <- all( primary_direction$ animal_mean_contrast_mean > 0 )
  if (!all_primary_positive) {
    add_review( "R02", "target_localization", "Not every prespecified timepoint/zone contrast has a positive mean adjusted-target contrast.", paste( paste0( primary_direction$timepoint, ":", primary_direction$contrast, "=", round( primary_direction$ animal_mean_contrast_mean, 4 ) ), collapse = " | " ), "Interpret the program as spatially heterogeneous and retain the persistent inflammatory-remodeling name.", "selectively report only positive zones or animals." )
  }
  critical_failures <- checks[ critical == TRUE & status == "FAIL" ]
  review_items <- checks[ status == "REVIEW" ]
  final_status <- if ( nrow( critical_failures ) > 0L ) {
    "STOP"
  } else if (
    nrow( review_items ) > 0L
  ) {
    "PASS_WITH_REVIEW"
  } else {
    "PASS"
  }
  handoff <- list( step = "Step12F", completed_at = as.character( Sys.time() ), final_status = final_status, lesion_zone_definition = list( positive_structural_types = STRUCTURAL_POSITIVE_TYPES, negative_structural_types = STRUCTURAL_NEGATIVE_TYPES, explicitly_excluded = forbidden_zone_variables, K_smooth = K_SMOOTH, core_seed_quantile = CORE_SEED_QUANTILE, connectivity_radius_NN = SEED_CONNECTIVITY_RADIUS_NN, border_max_distance_NN = BORDER_MAX_DISTANCE_NN, perilesional_max_distance_NN = PERILESIONAL_MAX_DISTANCE_NN, fallback_border_quantile = FALLBACK_BORDER_DISTANCE_QUANTILE, fallback_perilesional_quantile = FALLBACK_PERILESIONAL_DISTANCE_QUANTILE ), primary_outcome = "target_program_adjusted_z_within_section", inferential_unit = "biological_animal", n_spots = nrow( spot_zones ), n_sections = uniqueN( spot_zones$section_id ), n_injured_sections = uniqueN( injured_spots$section_id ), spot_zone_RDS = OUT_SPOT_ZONES_RDS, animal_contrast_CSV = OUT_ANIMAL_CONTRASTS, animal_inference_CSV = OUT_ANIMAL_INFERENCE, zone_maps_PDF = OUT_ZONE_MAPS, program_name_upgraded_to_scar_associated = FALSE, recipient_neighborhood_analysis_run = FALSE, next_step = paste( "Step 12G should be chosen only after reviewing Step 12F.", "If lesion-proximal enrichment is consistent across animals, test", "recipient-cell neighborhoods and ligand-receptor spatial closure;", "otherwise retain the descriptive spatial-validation conclusion." ) )
  saveRDS( handoff, OUT_HANDOFF_RDS, compress = FALSE )
  data.table::fwrite( checks, OUT_CHECKS, bom = TRUE )
  total_elapsed <- elapsed_minutes( workflow_start )
  primary_lines <- paste0( "- ", primary_direction$timepoint, " | ", primary_direction$contrast, ": mean animal contrast=", round( primary_direction$ animal_mean_contrast_mean, 4 ), "; positive animals=", primary_direction$positive_animals, "/", primary_direction$n_animals, "; exact sign P=", format( primary_direction$ exact_sign_test_p, digits = 4 ), "; Wilcoxon P=", format( primary_direction$ wilcoxon_signed_rank_p, digits = 4 ) )
  if (nrow( critical_failures ) > 0L) {
  }
  if (nrow( review_items ) > 0L) {
  }
  writeLines( capture.output( sessionInfo() ), SESSION_INFO_FILE, useBytes = TRUE )
  write_run_log()
  invisible( list( final_status = final_status, spot_zone_RDS = OUT_SPOT_ZONES_RDS, handoff_RDS = OUT_HANDOFF_RDS, primary_animal_results = primary_direction, checks = checks ) )
}
tryCatch( main(), error = function(e) {
    write_fatal_report( conditionMessage( e ) )
  }
)
