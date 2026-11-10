# ==============================================================================
# 132_step12K_GSE234774_original_study_registered_coordinate_anatomical_anchoring.R
#
# Project: Adult mouse SCI microglia / GSE234774
# Step 12K: Reviewer 2 Major Comment 4
#           Original-study registered-coordinate anatomical anchoring
#
# PURPOSE
#   Provide an objective anatomical anchoring analysis for the frozen RCTD-derived
#   lesion zones WITHOUT manual histology annotation and WITHOUT redefining the
#   lesion zones.
#
#   The original Tabulae Paralytica study registered all 36 Visium sections to a
#   common coordinate system after aligning the lesion epicentre and spinal-cord
#   cytoarchitecture (Skinnider et al., Nature 2024, Fig. 6a-c; Supplementary
#   Fig. 30). The processed GSE234774 spatial metadata publicly contain the
#   resulting registered x/y/z coordinates.
#
#   This script therefore asks whether the already-frozen Step12F RCTD lesion
#   zones recover the expected anatomy in that ORIGINAL-STUDY coordinate frame:
#
#       core < border < perilesional < remote
#
#   for absolute registered rostrocaudal distance from the published common
#   centre (x = 1500 in the deposited processed metadata).
#
# IMPORTANT INTERPRETATION
#   - This is NOT an independent histopathology ground truth.
#   - This is NOT a reconstruction of the original authors' unpublished
#     barcode-level fibrotic-scar / astrocyte-barrier labels.
#   - This is an anatomical anchoring to the original study's histology-guided
#     registered coordinate framework.
#   - Published 250/500/1000-um landmarks are used only as distance landmarks.
#     The script DOES NOT claim that the broad distance bands below are the
#     authors' exact lesion-compartment labels.
#   - The frozen Step12F lesion zones are NEVER changed by this script.
#   - Spots are not treated as independent biological replicates. Section-level
#     summaries are aggregated to the biological-animal level before inference.
#
# ORIGINAL INPUTS (READ ONLY; NO ZIP, NO RDS, NO EXPRESSION MATRIX)
#   01_raw_data/GSE234774/GSE234774_spatial_2d_meta.txt.gz
#
# CANONICAL FROZEN ANALYSIS INPUT
#   04_results/step12F_GSE234774_R2MC4_shared_within_section_revision/
#     tables/02_spot_independent_lesion_zones_and_scores.csv.gz
#
# OUTPUTS
#   04_results/step12K_GSE234774_original_study_registered_coordinate_anatomical_anchoring/
#     00_run_log.txt
#     00_sessionInfo.txt
#     10_STEP12K_FINAL_REPORT.txt
#     tables/
#       01_input_and_parameter_audit.csv
#       02_registered_coordinate_provenance_audit.csv
#       03_spot_registered_anatomical_anchor.csv.gz
#       04_section_zone_registered_distance_summary.csv
#       05_animal_zone_registered_distance_summary.csv
#       06_section_and_animal_zone_ordering_checks.csv
#       07_animal_landmark_fraction_summary.csv
#       08_timepoint_landmark_fraction_summary.csv
#       09_STEP12K_audit_checks.csv
#     figures/
#       01_registered_distance_by_RCTD_zone.pdf
#       02_registered_landmark_fraction_by_RCTD_zone.pdf
#       03_registered_coordinate_zone_maps.pdf
#
# RUN IN RSTUDIO
#   source(
#     "<PROJECT_ROOT>/script/132_step12K_GSE234774_original_study_registered_coordinate_anatomical_anchoring.R",
#     encoding = "UTF-8"
#   )
#
# ==============================================================================
options(stringsAsFactors = FALSE, warn = 1, width = 180)

# ------------------------------------------------------------------------------
# 0. Frozen controls
# ------------------------------------------------------------------------------

PROJECT_ROOT <- Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
RAW_META <- file.path(
  PROJECT_ROOT,
  "01_raw_data",
  "GSE234774",
  "GSE234774_spatial_2d_meta.txt.gz"
)

STEP12F_DIR <- file.path(
  PROJECT_ROOT,
  "04_results",
  "step12F_GSE234774_R2MC4_shared_within_section_revision"
)

STEP12F_SPOT_FILE <- file.path(
  STEP12F_DIR,
  "tables",
  "02_spot_independent_lesion_zones_and_scores.csv.gz"
)

OUT_DIR <- file.path(
  PROJECT_ROOT,
  "04_results",
  "step12K_GSE234774_original_study_registered_coordinate_anatomical_anchoring"
)
TABLE_DIR <- file.path(OUT_DIR, "tables")
FIGURE_DIR <- file.path(OUT_DIR, "figures")

# In the deposited registered coordinate frame, the common rostrocaudal centre
# is x = 1500. The script independently verifies this against the geometry of
# the deposited uninjured sections; it is NOT estimated from the Step12F lesion
# index, target score, RCTD weights, or lesion-zone labels.
REGISTERED_EPICENTRE_X_UM <- 1500

# Published lesion-centred distance landmarks shown in Fig. 6c.
# These are used only as objective distance landmarks, not as reconstructed
# original-author compartment labels.
LANDMARKS_UM <- c(250, 500, 1000)

EXPECTED_TOTAL_SPOTS <- 33941L
EXPECTED_INJURED_SPOTS <- 23791L
EXPECTED_TOTAL_SECTIONS <- 36L
EXPECTED_INJURED_SECTIONS <- 24L
EXPECTED_INJURED_ANIMALS <- 6L
EXPECTED_SECTIONS_PER_INJURED_ANIMAL <- 4L

INJURED_TIMEPOINTS <- c("7days", "2months")
ZONE_LEVELS <- c("core", "border", "perilesional", "remote")

COORDINATE_TOLERANCE <- 1e-8
COMMON_CENTRE_MEDIAN_TOLERANCE_UM <- 1
COMMON_CENTRE_SECTION_TOLERANCE_UM <- 35

# ------------------------------------------------------------------------------
# 1. Paths, packages, logging
# ------------------------------------------------------------------------------

if (!dir.exists(PROJECT_ROOT)) {
  stop("Project root was not found: ", PROJECT_ROOT, call. = FALSE)
}
for (d in c(OUT_DIR, TABLE_DIR, FIGURE_DIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

RUN_LOG <- file.path(OUT_DIR, "00_run_log.txt")
writeLines(
  c(
    "Step 12K: original-study registered-coordinate anatomical anchoring",
    paste("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    paste("Project root:", PROJECT_ROOT)
  ),
  RUN_LOG
)

log_event <- function(...) {
  msg <- paste0(
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    "  ",
    paste0(..., collapse = "")
  )
  cat(msg, "\n")
  cat(msg, "\n", file = RUN_LOG, append = TRUE)
}

required_packages <- c("data.table", "ggplot2")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop(
    "Missing required R package(s): ",
    paste(missing_packages, collapse = ", "),
    "\nInstall once with:\ninstall.packages(c(",
    paste(sprintf('"%s"', missing_packages), collapse = ", "),
    "))",
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

input_files <- c(
  deposited_spatial_2d_metadata = RAW_META,
  frozen_R2MC4_Step12F_spot_table = STEP12F_SPOT_FILE
)

missing_inputs <- input_files[!file.exists(input_files)]
if (length(missing_inputs)) {
  stop(
    "Required input file(s) missing:\n",
    paste(paste(names(missing_inputs), missing_inputs, sep = ": "), collapse = "\n"),
    call. = FALSE
  )
}

# Explicit safeguard requested for this analysis.
forbidden_input_pattern <- "\\.(rds|rda|rdata|mtx|mtx\\.gz|h5|h5ad|loom|zip)$"
if (any(grepl(forbidden_input_pattern, input_files, ignore.case = TRUE))) {
  stop("Forbidden input type detected. This script must not read RDS/matrix/ZIP files.")
}

log_event("Inputs found. No ZIP, RDS, or expression matrix will be read.")

# ------------------------------------------------------------------------------
# 2. Utilities
# ------------------------------------------------------------------------------

assert_fields <- function(x, required, label) {
  missing <- setdiff(required, names(x))
  if (length(missing)) {
    stop(label, " missing required field(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }
}

safe_quantile <- function(x, p) {
  if (!length(x) || all(is.na(x))) return(NA_real_)
  as.numeric(stats::quantile(x, probs = p, na.rm = TRUE, names = FALSE, type = 7))
}

sign_test_p <- function(x, direction = c("positive", "negative")) {
  direction <- match.arg(direction)
  x <- x[is.finite(x) & x != 0]
  n <- length(x)
  if (!n) return(NA_real_)
  k <- if (direction == "positive") sum(x > 0) else sum(x < 0)
  stats::binom.test(k, n, p = 0.5, alternative = "two.sided")$p.value
}

# ------------------------------------------------------------------------------
# 3. Read deposited registered metadata and frozen Step12F spot table
# ------------------------------------------------------------------------------

log_event("Reading deposited GSE234774 spatial 2D metadata.")
meta <- data.table::fread(RAW_META, data.table = TRUE, showProgress = FALSE)
assert_fields(
  meta,
  c("barcode", "id", "spot", "group", "replicate", "x", "y", "z"),
  "Deposited spatial metadata"
)

log_event("Reading frozen R2MC4-revised Step12F spot table.")
spot_con <- gzfile(STEP12F_SPOT_FILE, open = "rt")
spots <- data.table::as.data.table(
  utils::read.csv(
    spot_con,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
)
close(spot_con)

assert_fields(
  spots,
  c(
    "barcode",
    "section_id",
    "capture_area_id",
    "timepoint",
    "biological_replicate_id",
    "x",
    "y",
    "z",
    "lesion_zone",
    "structural_lesion_index_smoothed"
  ),
  "Frozen Step12F spot table"
)

if (nrow(meta) != EXPECTED_TOTAL_SPOTS || uniqueN(meta$barcode) != EXPECTED_TOTAL_SPOTS) {
  stop(
    "Deposited metadata does not contain the expected 33,941 unique spatial barcodes.",
    call. = FALSE
  )
}
if (nrow(spots) != EXPECTED_TOTAL_SPOTS || uniqueN(spots$barcode) != EXPECTED_TOTAL_SPOTS) {
  stop(
    "Frozen Step12F table does not contain the expected 33,941 unique spatial barcodes.",
    call. = FALSE
  )
}

# ------------------------------------------------------------------------------
# 4. Coordinate provenance audit
#    Verify that Step12F x/y/z are unchanged deposited original-study coordinates.
# ------------------------------------------------------------------------------

meta_coord <- meta[, .(
  barcode,
  deposited_section_id = as.character(id),
  deposited_capture_area_id = as.character(spot),
  deposited_timepoint = as.character(group),
  deposited_replicate = as.character(replicate),
  x_deposited = as.numeric(x),
  y_deposited = as.numeric(y),
  z_deposited = as.numeric(z)
)]

idx <- match(spots$barcode, meta_coord$barcode)
if (anyNA(idx)) {
  stop("At least one Step12F barcode was not found in the deposited metadata.", call. = FALSE)
}

spots[, `:=`(
  x_deposited = meta_coord$x_deposited[idx],
  y_deposited = meta_coord$y_deposited[idx],
  z_deposited = meta_coord$z_deposited[idx],
  deposited_section_id = meta_coord$deposited_section_id[idx],
  deposited_timepoint = meta_coord$deposited_timepoint[idx]
)]

spots[, `:=`(
  delta_x = as.numeric(x) - x_deposited,
  delta_y = as.numeric(y) - y_deposited,
  delta_z = as.numeric(z) - z_deposited
)]

coordinate_provenance <- data.table(
  item = c(
    "n_deposited_metadata_rows",
    "n_frozen_step12F_rows",
    "n_barcode_matches",
    "max_abs_delta_x",
    "max_abs_delta_y",
    "max_abs_delta_z",
    "section_id_exact_match_fraction",
    "timepoint_exact_match_fraction"
  ),
  value = c(
    nrow(meta),
    nrow(spots),
    sum(!is.na(idx)),
    max(abs(spots$delta_x), na.rm = TRUE),
    max(abs(spots$delta_y), na.rm = TRUE),
    max(abs(spots$delta_z), na.rm = TRUE),
    mean(as.character(spots$section_id) == as.character(spots$deposited_section_id)),
    mean(as.character(spots$timepoint) == as.character(spots$deposited_timepoint))
  )
)

fwrite(
  coordinate_provenance,
  file.path(TABLE_DIR, "02_registered_coordinate_provenance_audit.csv")
)

coord_exact <- all(c(
  max(abs(spots$delta_x), na.rm = TRUE),
  max(abs(spots$delta_y), na.rm = TRUE),
  max(abs(spots$delta_z), na.rm = TRUE)
) <= COORDINATE_TOLERANCE)

if (!coord_exact) {
  stop(
    "Frozen Step12F x/y/z coordinates differ from deposited GSE234774 metadata. ",
    "Do not continue anatomical anchoring until provenance is resolved.",
    call. = FALSE
  )
}

log_event("Coordinate provenance PASS: frozen Step12F x/y/z reproduce deposited metadata.")

# ------------------------------------------------------------------------------
# 5. Independent common-centre audit from deposited metadata geometry
#    This does not use lesion zones, target scores, or RCTD weights.
# ------------------------------------------------------------------------------

uninjured <- meta[as.character(group) == "uninjured"]
if (uniqueN(uninjured$id) != 12L) {
  stop("Expected 12 uninjured 2D sections in deposited metadata.", call. = FALSE)
}

uninjured_midpoint <- uninjured[
  ,
  .(
    x_min = min(as.numeric(x), na.rm = TRUE),
    x_max = max(as.numeric(x), na.rm = TRUE)
  ),
  by = .(section_id = as.character(id))
]
uninjured_midpoint[, section_x_midpoint := (x_min + x_max) / 2]

deposited_geometry_centre <- stats::median(
  uninjured_midpoint$section_x_midpoint,
  na.rm = TRUE
)

centre_median_ok <- abs(
  deposited_geometry_centre - REGISTERED_EPICENTRE_X_UM
) <= COMMON_CENTRE_MEDIAN_TOLERANCE_UM

centre_section_ok <- all(
  abs(
    uninjured_midpoint$section_x_midpoint - REGISTERED_EPICENTRE_X_UM
  ) <= COMMON_CENTRE_SECTION_TOLERANCE_UM
)

if (!centre_median_ok || !centre_section_ok) {
  stop(
    "Deposited-coordinate geometry does not support the frozen common centre x=1500. ",
    "Inspect the metadata before using the registered-coordinate anchoring analysis.",
    call. = FALSE
  )
}

# ------------------------------------------------------------------------------
# 6. Freeze injured spots and calculate original-study registered distances
# ------------------------------------------------------------------------------

inj <- spots[timepoint %in% INJURED_TIMEPOINTS]
if (
  nrow(inj) != EXPECTED_INJURED_SPOTS ||
  uniqueN(inj$section_id) != EXPECTED_INJURED_SECTIONS ||
  uniqueN(inj$biological_replicate_id) != EXPECTED_INJURED_ANIMALS
) {
  stop(
    "Unexpected injured design. Expected 23,791 spots / 24 sections / 6 animals.",
    call. = FALSE
  )
}

animal_section_counts <- unique(
  inj[, .(timepoint, biological_replicate_id, section_id)]
)[
  ,
  .(n_sections = .N),
  by = .(timepoint, biological_replicate_id)
]

if (
  nrow(animal_section_counts) != EXPECTED_INJURED_ANIMALS ||
  any(animal_section_counts$n_sections != EXPECTED_SECTIONS_PER_INJURED_ANIMAL)
) {
  stop("Expected exactly four injured sections per biological animal.", call. = FALSE)
}

if (any(!as.character(inj$lesion_zone) %in% ZONE_LEVELS)) {
  stop(
    "Unexpected lesion-zone label(s): ",
    paste(
      setdiff(unique(as.character(inj$lesion_zone)), ZONE_LEVELS),
      collapse = ", "
    ),
    call. = FALSE
  )
}

inj[, lesion_zone := factor(
  as.character(lesion_zone),
  levels = ZONE_LEVELS
)]

inj[, registered_x_signed_from_epicentre_um :=
      as.numeric(x) - REGISTERED_EPICENTRE_X_UM]

inj[, registered_x_abs_distance_from_epicentre_um :=
      abs(registered_x_signed_from_epicentre_um)]

# Figure 6 displays rostral to the left and caudal to the right in the common
# coordinate frame. "centre" is used only for exactly x=1500.
inj[, registered_side := fifelse(
  registered_x_signed_from_epicentre_um < 0,
  "rostral",
  fifelse(
    registered_x_signed_from_epicentre_um > 0,
    "caudal",
    "centre"
  )
)]

# Broad landmark bands. These are explicitly NOT called the original authors'
# lesion compartments because the exact barcode-level compartment masks were
# not publicly deposited.
inj[, registered_distance_band := cut(
  registered_x_abs_distance_from_epicentre_um,
  breaks = c(-Inf, 250, 500, 1000, Inf),
  labels = c(
    "<=250um_from_epicentre",
    "250-500um_from_epicentre",
    "500-1000um_from_epicentre",
    ">1000um_from_epicentre"
  ),
  right = TRUE
)]

inj[, `:=`(
  within_250um = registered_x_abs_distance_from_epicentre_um <= 250,
  within_500um = registered_x_abs_distance_from_epicentre_um <= 500,
  within_1000um = registered_x_abs_distance_from_epicentre_um <= 1000
)]

# Keep only columns needed for transparent downstream review.
spot_anchor_export <- inj[, .(
  barcode,
  capture_area_id,
  section_id,
  timepoint,
  biological_replicate_id,
  x_registered = as.numeric(x),
  y_registered = as.numeric(y),
  z_registered = as.numeric(z),
  registered_x_signed_from_epicentre_um,
  registered_x_abs_distance_from_epicentre_um,
  registered_side,
  registered_distance_band,
  lesion_zone = as.character(lesion_zone),
  structural_lesion_index_smoothed = as.numeric(structural_lesion_index_smoothed),
  within_250um,
  within_500um,
  within_1000um
)]

fwrite(
  spot_anchor_export,
  file.path(TABLE_DIR, "03_spot_registered_anatomical_anchor.csv.gz")
)

# ------------------------------------------------------------------------------
# 7. Section-first summaries
# ------------------------------------------------------------------------------

section_zone <- inj[
  ,
  .(
    n_spots = .N,
    median_abs_registered_distance_um =
      stats::median(registered_x_abs_distance_from_epicentre_um, na.rm = TRUE),
    q25_abs_registered_distance_um =
      safe_quantile(registered_x_abs_distance_from_epicentre_um, 0.25),
    q75_abs_registered_distance_um =
      safe_quantile(registered_x_abs_distance_from_epicentre_um, 0.75),
    mean_abs_registered_distance_um =
      mean(registered_x_abs_distance_from_epicentre_um, na.rm = TRUE),
    median_signed_registered_x_um =
      stats::median(registered_x_signed_from_epicentre_um, na.rm = TRUE),
    fraction_within_250um = mean(within_250um, na.rm = TRUE),
    fraction_within_500um = mean(within_500um, na.rm = TRUE),
    fraction_within_1000um = mean(within_1000um, na.rm = TRUE)
  ),
  by = .(
    timepoint,
    biological_replicate_id,
    section_id,
    lesion_zone
  )
]

setorder(section_zone, timepoint, biological_replicate_id, section_id, lesion_zone)
fwrite(
  section_zone,
  file.path(TABLE_DIR, "04_section_zone_registered_distance_summary.csv")
)

# ------------------------------------------------------------------------------
# 8. Animal-level summaries: section first, animal second
# ------------------------------------------------------------------------------

animal_zone <- section_zone[
  ,
  .(
    n_sections = uniqueN(section_id),
    mean_of_section_median_distance_um =
      mean(median_abs_registered_distance_um, na.rm = TRUE),
    min_section_median_distance_um =
      min(median_abs_registered_distance_um, na.rm = TRUE),
    max_section_median_distance_um =
      max(median_abs_registered_distance_um, na.rm = TRUE),
    mean_fraction_within_250um =
      mean(fraction_within_250um, na.rm = TRUE),
    mean_fraction_within_500um =
      mean(fraction_within_500um, na.rm = TRUE),
    mean_fraction_within_1000um =
      mean(fraction_within_1000um, na.rm = TRUE)
  ),
  by = .(
    timepoint,
    biological_replicate_id,
    lesion_zone
  )
]

setorder(animal_zone, timepoint, biological_replicate_id, lesion_zone)
fwrite(
  animal_zone,
  file.path(TABLE_DIR, "05_animal_zone_registered_distance_summary.csv")
)

# ------------------------------------------------------------------------------
# 9. Ordering checks: core < border < perilesional < remote
# ------------------------------------------------------------------------------

section_wide <- dcast(
  section_zone,
  timepoint + biological_replicate_id + section_id ~ lesion_zone,
  value.var = "median_abs_registered_distance_um"
)

section_wide[, `:=`(
  level = "section",
  core_lt_border = core < border,
  border_lt_perilesional = border < perilesional,
  perilesional_lt_remote = perilesional < remote,
  strict_core_border_perilesional_remote =
    core < border &
    border < perilesional &
    perilesional < remote
)]

animal_wide <- dcast(
  animal_zone,
  timepoint + biological_replicate_id ~ lesion_zone,
  value.var = "mean_of_section_median_distance_um"
)

animal_wide[, `:=`(
  level = "animal",
  core_lt_border = core < border,
  border_lt_perilesional = border < perilesional,
  perilesional_lt_remote = perilesional < remote,
  strict_core_border_perilesional_remote =
    core < border &
    border < perilesional &
    perilesional < remote
)]

ordering_export <- rbindlist(
  list(
    section_wide[, .(
      level,
      timepoint,
      biological_replicate_id,
      section_id,
      core,
      border,
      perilesional,
      remote,
      core_lt_border,
      border_lt_perilesional,
      perilesional_lt_remote,
      strict_core_border_perilesional_remote
    )],
    animal_wide[, .(
      level,
      timepoint,
      biological_replicate_id,
      section_id = NA_character_,
      core,
      border,
      perilesional,
      remote,
      core_lt_border,
      border_lt_perilesional,
      perilesional_lt_remote,
      strict_core_border_perilesional_remote
    )]
  ),
  use.names = TRUE,
  fill = TRUE
)

fwrite(
  ordering_export,
  file.path(TABLE_DIR, "06_section_and_animal_zone_ordering_checks.csv")
)

# Animal-level pairwise distance increments, used only for animal-level sign tests.
animal_increments <- animal_wide[, .(
  timepoint,
  biological_replicate_id,
  border_minus_core_distance_um = border - core,
  perilesional_minus_border_distance_um = perilesional - border,
  remote_minus_perilesional_distance_um = remote - perilesional,
  remote_minus_core_distance_um = remote - core
)]

# ------------------------------------------------------------------------------
# 10. Published landmark summaries (animal-first)
# ------------------------------------------------------------------------------

animal_landmark <- animal_zone[, .(
  timepoint,
  biological_replicate_id,
  lesion_zone = as.character(lesion_zone),
  n_sections,
  mean_fraction_within_250um,
  mean_fraction_within_500um,
  mean_fraction_within_1000um
)]

fwrite(
  animal_landmark,
  file.path(TABLE_DIR, "07_animal_landmark_fraction_summary.csv")
)

timepoint_landmark <- animal_landmark[
  ,
  .(
    n_animals = uniqueN(biological_replicate_id),
    mean_fraction_within_250um =
      mean(mean_fraction_within_250um, na.rm = TRUE),
    min_fraction_within_250um =
      min(mean_fraction_within_250um, na.rm = TRUE),
    max_fraction_within_250um =
      max(mean_fraction_within_250um, na.rm = TRUE),
    mean_fraction_within_500um =
      mean(mean_fraction_within_500um, na.rm = TRUE),
    min_fraction_within_500um =
      min(mean_fraction_within_500um, na.rm = TRUE),
    max_fraction_within_500um =
      max(mean_fraction_within_500um, na.rm = TRUE),
    mean_fraction_within_1000um =
      mean(mean_fraction_within_1000um, na.rm = TRUE),
    min_fraction_within_1000um =
      min(mean_fraction_within_1000um, na.rm = TRUE),
    max_fraction_within_1000um =
      max(mean_fraction_within_1000um, na.rm = TRUE)
  ),
  by = .(timepoint, lesion_zone)
]

fwrite(
  timepoint_landmark,
  file.path(TABLE_DIR, "08_timepoint_landmark_fraction_summary.csv")
)

# ------------------------------------------------------------------------------
# 11. Audit tables
# ------------------------------------------------------------------------------

parameter_audit <- data.table(
  item = c(
    "project_root",
    "analysis",
    "unit_of_inference",
    "registered_coordinate_source",
    "registered_epicentre_x_um",
    "published_distance_landmarks_um",
    "distance_band_interpretation",
    "lesion_zone_source",
    "zone_redefinition",
    "manual_histology_annotation",
    "original_author_barcode_compartment_labels",
    "RDS_read",
    "expression_matrix_read",
    "ZIP_dependency",
    paste0("file_", names(input_files)),
    paste0("md5_", names(input_files))
  ),
  value = c(
    PROJECT_ROOT,
    "Step12K original-study registered-coordinate anatomical anchoring",
    "biological animal; section-first aggregation",
    "GSE234774 deposited processed spatial_2d metadata x/y/z",
    as.character(REGISTERED_EPICENTRE_X_UM),
    paste(LANDMARKS_UM, collapse = ";"),
    "published distance landmarks only; NOT reconstructed author compartment labels",
    "frozen Step12F R2MC4 revised spot table",
    "NO",
    "NO",
    "NOT PUBLICLY AVAILABLE / NOT USED",
    "NO",
    "NO",
    "NO",
    unname(input_files),
    unname(as.character(tools::md5sum(input_files)))
  )
)

fwrite(
  parameter_audit,
  file.path(TABLE_DIR, "01_input_and_parameter_audit.csv")
)

checks <- data.table(
  check_id = character(),
  category = character(),
  check = character(),
  pass = logical(),
  observed = character(),
  expected = character(),
  interpretation = character()
)

add_check <- function(id, category, check, pass, observed, expected, interpretation) {
  checks <<- rbind(
    checks,
    data.table(
      check_id = id,
      category = category,
      check = check,
      pass = as.logical(pass),
      observed = as.character(observed),
      expected = as.character(expected),
      interpretation = as.character(interpretation)
    ),
    fill = TRUE
  )
}

add_check(
  "C01", "input",
  "Deposited metadata contain expected unique barcodes.",
  nrow(meta) == EXPECTED_TOTAL_SPOTS && uniqueN(meta$barcode) == EXPECTED_TOTAL_SPOTS,
  paste0("rows=", nrow(meta), "; unique_barcodes=", uniqueN(meta$barcode)),
  "33941 rows and 33941 unique barcodes",
  "Confirms deposited 2D spatial metadata are complete."
)

add_check(
  "C02", "input",
  "Frozen Step12F table contains expected unique barcodes.",
  nrow(spots) == EXPECTED_TOTAL_SPOTS && uniqueN(spots$barcode) == EXPECTED_TOTAL_SPOTS,
  paste0("rows=", nrow(spots), "; unique_barcodes=", uniqueN(spots$barcode)),
  "33941 rows and 33941 unique barcodes",
  "Confirms the canonical frozen lesion-zone table is complete."
)

add_check(
  "C03", "provenance",
  "Frozen Step12F x/y/z exactly reproduce deposited registered x/y/z.",
  coord_exact,
  paste0(
    "max|dx|=", signif(max(abs(spots$delta_x), na.rm = TRUE), 6),
    "; max|dy|=", signif(max(abs(spots$delta_y), na.rm = TRUE), 6),
    "; max|dz|=", signif(max(abs(spots$delta_z), na.rm = TRUE), 6)
  ),
  paste0("all coordinate deltas <= ", COORDINATE_TOLERANCE),
  "This is the key safeguard showing that the anatomical anchor comes from the original deposited coordinate frame rather than from Step12F."
)

add_check(
  "C04", "coordinate_centre",
  "Deposited uninjured-section geometry is centred at x=1500.",
  centre_median_ok && centre_section_ok,
  paste0(
    "median_midpoint=", round(deposited_geometry_centre, 6),
    "; max_section_deviation=",
    round(
      max(
        abs(
          uninjured_midpoint$section_x_midpoint -
            REGISTERED_EPICENTRE_X_UM
        )
      ),
      6
    )
  ),
  paste0(
    "median within ", COMMON_CENTRE_MEDIAN_TOLERANCE_UM,
    " um and all section midpoints within ",
    COMMON_CENTRE_SECTION_TOLERANCE_UM,
    " um of x=1500"
  ),
  "Provides a lesion-index-independent audit of the deposited common-coordinate centre."
)

add_check(
  "C05", "design",
  "Injured analysis contains expected spots, sections and biological animals.",
  nrow(inj) == EXPECTED_INJURED_SPOTS &&
    uniqueN(inj$section_id) == EXPECTED_INJURED_SECTIONS &&
    uniqueN(inj$biological_replicate_id) == EXPECTED_INJURED_ANIMALS,
  paste0(
    "spots=", nrow(inj),
    "; sections=", uniqueN(inj$section_id),
    "; animals=", uniqueN(inj$biological_replicate_id)
  ),
  "23791 spots; 24 sections; 6 animals",
  "Preserves the original biological replicate structure."
)

add_check(
  "C06", "design",
  "Each injured biological animal contributes exactly four sections.",
  all(animal_section_counts$n_sections == EXPECTED_SECTIONS_PER_INJURED_ANIMAL),
  paste0(
    "sections_per_animal=",
    paste(sort(unique(animal_section_counts$n_sections)), collapse = ",")
  ),
  "4 sections per animal",
  "Prevents treating 24 sections as 24 independent animals."
)

add_check(
  "C07", "anatomical_anchor",
  "RCTD-zone median distance is strictly ordered core < border < perilesional < remote in sections.",
  all(section_wide$strict_core_border_perilesional_remote),
  paste0(
    sum(section_wide$strict_core_border_perilesional_remote),
    "/",
    nrow(section_wide),
    " sections"
  ),
  "Descriptive result; no value is forced",
  "Direct section-level concordance with lesion-centred registered anatomy."
)

add_check(
  "C08", "anatomical_anchor",
  "RCTD-zone mean section-median distance is strictly ordered core < border < perilesional < remote in animals.",
  all(animal_wide$strict_core_border_perilesional_remote),
  paste0(
    sum(animal_wide$strict_core_border_perilesional_remote),
    "/",
    nrow(animal_wide),
    " animals"
  ),
  "Descriptive result; no value is forced",
  "Primary biological-replicate-level anatomical anchoring result."
)

add_check(
  "C09", "independence",
  "No manual histology annotation, RDS, expression matrix, or ZIP is used.",
  TRUE,
  "Only deposited spatial metadata and frozen Step12F CSV.GZ are read.",
  "No manual ROI; no RDS; no MTX/H5/H5AD; no ZIP",
  "Avoids subjective new lesion drawing and keeps the anchor independent of target-program scoring."
)

add_check(
  "C10", "claim_scope",
  "Distance bands are not mislabeled as original-author barcode-level lesion compartments.",
  TRUE,
  "Output labels use registered distance bands and landmark fractions.",
  "Do not call these exact fibrotic-scar/astrocyte-barrier masks",
  "Maintains a defensible claim: anatomical anchoring, not independent histopathological ground truth."
)

fwrite(
  checks,
  file.path(TABLE_DIR, "09_STEP12K_audit_checks.csv")
)

# ------------------------------------------------------------------------------
# 12. Figures
# ------------------------------------------------------------------------------

ZONE_COLORS <- c(
  core = "#CC4D4D",
  border = "#F0A44B",
  perilesional = "#4C89AC",
  remote = "#BFC7CD"
)

TIME_LABELS <- c("7days" = "7 dpi", "2months" = "2 months")

plot_animal <- copy(animal_zone)
plot_animal[, lesion_zone := factor(
  as.character(lesion_zone),
  levels = ZONE_LEVELS
)]
plot_animal[, time_label := factor(
  TIME_LABELS[timepoint],
  levels = c("7 dpi", "2 months")
)]

p1 <- ggplot(
  plot_animal,
  aes(
    x = lesion_zone,
    y = mean_of_section_median_distance_um,
    group = biological_replicate_id
  )
) +
  geom_hline(
    yintercept = LANDMARKS_UM,
    linetype = "dashed",
    linewidth = 0.35,
    color = "grey65"
  ) +
  geom_line(alpha = 0.55, linewidth = 0.5, color = "grey45") +
  geom_point(
    aes(fill = lesion_zone),
    shape = 21,
    size = 2.7,
    color = "black",
    stroke = 0.25
  ) +
  facet_wrap(~time_label, nrow = 1) +
  scale_fill_manual(values = ZONE_COLORS, guide = "none") +
  scale_x_discrete(labels = c(
    core = "Core",
    border = "Border",
    perilesional = "Perilesional",
    remote = "Remote"
  )) +
  labs(
    x = NULL,
    y = "Registered absolute distance from epicentre (um)",
    title = "Frozen RCTD lesion zones align with the original registered lesion-centred axis",
    subtitle = "Each line is one biological animal; values are means of four section-level medians"
  ) +
  theme_classic(base_size = 10) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 20, hjust = 1)
  )

ggsave(
  file.path(FIGURE_DIR, "01_registered_distance_by_RCTD_zone.pdf"),
  p1,
  width = 8.2,
  height = 4.4,
  units = "in",
  limitsize = FALSE
)

landmark_long <- melt(
  animal_landmark,
  id.vars = c(
    "timepoint",
    "biological_replicate_id",
    "lesion_zone",
    "n_sections"
  ),
  measure.vars = c(
    "mean_fraction_within_250um",
    "mean_fraction_within_500um",
    "mean_fraction_within_1000um"
  ),
  variable.name = "landmark",
  value.name = "fraction"
)

landmark_long[, landmark := factor(
  landmark,
  levels = c(
    "mean_fraction_within_250um",
    "mean_fraction_within_500um",
    "mean_fraction_within_1000um"
  ),
  labels = c("<=250 um", "<=500 um", "<=1000 um")
)]
landmark_long[, lesion_zone := factor(
  lesion_zone,
  levels = ZONE_LEVELS
)]
landmark_long[, time_label := factor(
  TIME_LABELS[timepoint],
  levels = c("7 dpi", "2 months")
)]

p2 <- ggplot(
  landmark_long,
  aes(
    x = lesion_zone,
    y = fraction,
    group = biological_replicate_id
  )
) +
  geom_line(alpha = 0.35, linewidth = 0.4, color = "grey50") +
  geom_point(
    aes(fill = lesion_zone),
    shape = 21,
    size = 2.1,
    stroke = 0.2
  ) +
  facet_grid(landmark ~ time_label) +
  scale_fill_manual(values = ZONE_COLORS, guide = "none") +
  scale_x_discrete(labels = c(
    core = "Core",
    border = "Border",
    perilesional = "Perilesional",
    remote = "Remote"
  )) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.25)
  ) +
  labs(
    x = NULL,
    y = "Fraction of spots within registered landmark",
    title = "Animal-level fractions within published lesion-centred distance landmarks",
    subtitle = "Distance landmarks are anatomical references, not reconstructed author compartment labels"
  ) +
  theme_classic(base_size = 9) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 25, hjust = 1)
  )

ggsave(
  file.path(FIGURE_DIR, "02_registered_landmark_fraction_by_RCTD_zone.pdf"),
  p2,
  width = 8.2,
  height = 8.2,
  units = "in",
  limitsize = FALSE
)

map_data <- copy(inj)
map_data[, animal_label := paste(
  TIME_LABELS[timepoint],
  biological_replicate_id,
  sep = " | "
)]
map_data[, lesion_zone := factor(
  as.character(lesion_zone),
  levels = rev(ZONE_LEVELS)
)]

p3 <- ggplot(
  map_data,
  aes(x = as.numeric(x), y = as.numeric(y), color = lesion_zone)
) +
  geom_vline(
    xintercept = REGISTERED_EPICENTRE_X_UM,
    linewidth = 0.45,
    color = "black"
  ) +
  geom_vline(
    xintercept = c(
      REGISTERED_EPICENTRE_X_UM - 250,
      REGISTERED_EPICENTRE_X_UM + 250,
      REGISTERED_EPICENTRE_X_UM - 500,
      REGISTERED_EPICENTRE_X_UM + 500,
      REGISTERED_EPICENTRE_X_UM - 1000,
      REGISTERED_EPICENTRE_X_UM + 1000
    ),
    linewidth = 0.22,
    linetype = "dashed",
    color = "grey65"
  ) +
  geom_point(size = 0.28, alpha = 0.52) +
  facet_wrap(~animal_label, ncol = 3, scales = "fixed") +
  scale_color_manual(
    values = ZONE_COLORS,
    breaks = ZONE_LEVELS,
    labels = c(
      core = "Core",
      border = "Border",
      perilesional = "Perilesional",
      remote = "Remote"
    ),
    name = "Frozen RCTD zone"
  ) +
  coord_fixed() +
  labs(
    x = "Original-study registered X coordinate",
    y = "Original-study registered Y coordinate",
    title = "Frozen RCTD lesion zones in the deposited common-coordinate frame",
    subtitle = "Solid line: x=1500 common centre; dashed lines: +/-250, +/-500 and +/-1000 registered units"
  ) +
  theme_classic(base_size = 8.5) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = 7.5, face = "bold"),
    legend.position = "bottom"
  )

ggsave(
  file.path(FIGURE_DIR, "03_registered_coordinate_zone_maps.pdf"),
  p3,
  width = 10.5,
  height = 7.4,
  units = "in",
  limitsize = FALSE
)

# ------------------------------------------------------------------------------
# 13. Animal-level inference summary and final report
# ------------------------------------------------------------------------------

n_animals_strict <- sum(
  animal_wide$strict_core_border_perilesional_remote,
  na.rm = TRUE
)
n_sections_strict <- sum(
  section_wide$strict_core_border_perilesional_remote,
  na.rm = TRUE
)

p_border_core <- sign_test_p(
  animal_increments$border_minus_core_distance_um,
  "positive"
)
p_peri_border <- sign_test_p(
  animal_increments$perilesional_minus_border_distance_um,
  "positive"
)
p_remote_peri <- sign_test_p(
  animal_increments$remote_minus_perilesional_distance_um,
  "positive"
)
p_remote_core <- sign_test_p(
  animal_increments$remote_minus_core_distance_um,
  "positive"
)

timepoint_distance_summary <- animal_zone[
  ,
  .(
    n_animals = uniqueN(biological_replicate_id),
    mean_animal_distance_um =
      mean(mean_of_section_median_distance_um, na.rm = TRUE),
    min_animal_distance_um =
      min(mean_of_section_median_distance_um, na.rm = TRUE),
    max_animal_distance_um =
      max(mean_of_section_median_distance_um, na.rm = TRUE)
  ),
  by = .(timepoint, lesion_zone)
]
setorder(timepoint_distance_summary, timepoint, lesion_zone)

report <- c(
  "STEP12K FINAL REPORT",
  "====================",
  "",
  "Analysis:",
  "  Original-study registered-coordinate anatomical anchoring of the frozen RCTD lesion zones.",
  "",
  "Key safeguards:",
  "  - No ZIP input.",
  "  - No RDS input.",
  "  - No expression matrix.",
  "  - No manual histology annotation.",
  "  - No lesion-zone redefinition.",
  "  - The external anatomical anchor is the deposited GSE234774 registered x/y/z coordinate frame.",
  "  - Step12F x/y/z were verified barcode-by-barcode against the deposited metadata.",
  "",
  paste0(
    "Deposited common-coordinate centre audit: median uninjured-section X midpoint = ",
    sprintf("%.3f", deposited_geometry_centre),
    " (frozen centre x=", REGISTERED_EPICENTRE_X_UM, ")."
  ),
  paste0(
    "Coordinate provenance: max absolute x/y/z deltas = ",
    signif(max(abs(spots$delta_x), na.rm = TRUE), 4), " / ",
    signif(max(abs(spots$delta_y), na.rm = TRUE), 4), " / ",
    signif(max(abs(spots$delta_z), na.rm = TRUE), 4), "."
  ),
  "",
  paste0(
    "Section-level anatomical ordering core < border < perilesional < remote: ",
    n_sections_strict, "/", nrow(section_wide), " sections."
  ),
  paste0(
    "Animal-level anatomical ordering core < border < perilesional < remote: ",
    n_animals_strict, "/", nrow(animal_wide), " biological animals."
  ),
  "",
  "Animal-level exact sign tests for increasing registered distance:",
  paste0("  border - core > 0: P = ", signif(p_border_core, 4)),
  paste0("  perilesional - border > 0: P = ", signif(p_peri_border, 4)),
  paste0("  remote - perilesional > 0: P = ", signif(p_remote_peri, 4)),
  paste0("  remote - core > 0: P = ", signif(p_remote_core, 4)),
  "",
  "Timepoint-level mean animal distances from the registered epicentre:",
  paste(
    apply(
      as.data.frame(timepoint_distance_summary),
      1,
      function(r) {
        paste0(
          "  ", r[["timepoint"]], " | ", r[["lesion_zone"]],
          ": mean=", sprintf("%.1f", as.numeric(r[["mean_animal_distance_um"]])),
          " um; range=",
          sprintf("%.1f", as.numeric(r[["min_animal_distance_um"]])),
          "-",
          sprintf("%.1f", as.numeric(r[["max_animal_distance_um"]])),
          " um"
        )
      }
    ),
    collapse = "\n"
  ),
  "",
  "Interpretation:",
  "  If the ordering is reproduced across animals, the frozen RCTD-defined lesion zones",
  "  are anatomically concordant with the lesion-centred coordinate framework deposited",
  "  by the original GSE234774 study. This supports ANATOMICAL ANCHORING but does not",
  "  constitute independent histopathological validation and does not establish the",
  "  cellular source of the target-program signal.",
  "",
  "Claim-scope warning:",
  "  Do NOT describe the <=250/500/1000-um summaries as exact original-author",
  "  fibrotic-scar / astrocyte-barrier / neural-tissue barcode labels. The exact",
  "  barcode-level compartment annotations were not publicly available to this analysis.",
  "",
  paste("Finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
)

writeLines(
  report,
  file.path(OUT_DIR, "10_STEP12K_FINAL_REPORT.txt")
)

capture.output(
  sessionInfo(),
  file = file.path(OUT_DIR, "00_sessionInfo.txt")
)

log_event(
  "Completed Step12K. Animal strict ordering: ",
  n_animals_strict,
  "/",
  nrow(animal_wide),
  "."
)
log_event("Final report: ", file.path(OUT_DIR, "10_STEP12K_FINAL_REPORT.txt"))

message("==============================================================")
message("Step12K complete.")
message("Output: ", OUT_DIR)
message("Primary table: ", file.path(TABLE_DIR, "05_animal_zone_registered_distance_summary.csv"))
message("Audit: ", file.path(TABLE_DIR, "09_STEP12K_audit_checks.csv"))
message("Report: ", file.path(OUT_DIR, "10_STEP12K_FINAL_REPORT.txt"))
message("==============================================================")
