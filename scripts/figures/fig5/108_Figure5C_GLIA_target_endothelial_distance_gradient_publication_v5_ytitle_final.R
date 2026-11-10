# ============================================================
# 108_Figure5C_GLIA_target_endothelial_distance_gradient_publication_v5_ytitle_final.R
# v5 final visual fix: preserve all v4 content/layout except the C y-axis title,
# which is reduced slightly so the full label fits cleanly at final assembly size.
# Figure 5C | Shared target–endothelial distance-to-core gradient
# ------------------------------------------------------------
# Revised main-figure concept:
#   The previous section-correlation summary looked like an audit/QC plot.
#   Figure 5C is therefore redesigned as a spatial-gradient panel.
#
# Main message:
#   Target-program activity and endothelial-neighborhood signal are both
#   highest near the lesion core and decline with increasing distance,
#   at both 7 dpi and 2 months.
#
# Visual hierarchy:
#   thin lines  = biological-animal trajectories (descriptive replication)
#   thick lines = mean across biological animals
#
# Statistics:
#   - spots are NOT treated as independent replicates
#   - section-level bin means are first calculated
#   - section means are averaged within biological animal
#   - biological animal is the top-level unit
# ============================================================

options(stringsAsFactors = FALSE, warn = 1)

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(scales)
})

# ------------------------------------------------------------
# User settings
# ------------------------------------------------------------
ROOT_CANDIDATES <- unique(c(Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE)), "."))

OUTPUT_DIR_NAME <- "09_final_figures/GLIA/Figure5/single_panels/Figure5C_v5_ytitle_final"

BASE_FAMILY <- "Arial"
DPI_OUT     <- 600

TIMEPOINT_ORDER  <- c("7days", "2months")
TIMEPOINT_LABELS <- c("7days" = "7 dpi", "2months" = "2 months")

# 95% of the uploaded injured spots fall at approximately <= 30 NN units.
# Restricting the display to this range avoids a sparse far-tail dominating
# the horizontal scale while retaining the biologically relevant gradient.
MAX_DISTANCE_NN <- 30
BIN_WIDTH_NN    <- 3

MIN_SPOTS_PER_SECTION_BIN <- 8

TARGET_COLOR <- "#C44E52"
ENDO_COLOR   <- "#4C72B0"

THIN_LINE_ALPHA <- 0.26
THIN_LINE_WIDTH <- 0.55
MEAN_LINE_WIDTH <- 1.45
MEAN_POINT_SIZE <- 2.15

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------
stop2 <- function(...) stop(..., call. = FALSE)

locate_root <- function() {
  for (root in ROOT_CANDIDATES) {
    if (!dir.exists(root)) next

    key_original <- file.path(
      root,
      "04_results/step12F_GSE234774_independent_lesion_zones_target_program_validation/tables/02_spot_independent_lesion_zones_and_scores.csv.gz"
    )
    key_handoff <- file.path(root, "02_spot_independent_lesion_zones_and_scores.csv.gz")

    if (file.exists(key_original) || file.exists(key_handoff)) {
      return(normalizePath(root, winslash = "/", mustWork = TRUE))
    }
  }

  stop2("Could not locate a valid xjz project / Figure5 handoff root.")
}

file_in_root <- function(root, original_rel, handoff_name = basename(original_rel), must = TRUE) {
  candidates <- c(file.path(root, original_rel), file.path(root, handoff_name))
  hit <- candidates[file.exists(candidates)][1]

  if ((is.na(hit) || length(hit) == 0) && must) {
    stop2(
      "Required file not found:\n",
      "  original: ", original_rel, "\n",
      "  handoff:  ", handoff_name, "\n",
      "  root:     ", root
    )
  }

  if (is.na(hit) || length(hit) == 0) return(NA_character_)
  normalizePath(hit, winslash = "/", mustWork = TRUE)
}

extract_endothelial_scores <- function(neighborhood_df) {
  nm <- names(neighborhood_df)

  # Frozen independent recipient-neighborhood abundance.
  # DO NOT use target_neighbor_Endothelial_colocalization because that already
  # contains the target score and would create circularity.
  exact_col <- "neighbor_Endothelial_weight"

  if (exact_col %in% nm) {
    out <- neighborhood_df %>%
      dplyr::transmute(
        barcode = .data$barcode,
        section_id = .data$section_id,
        endothelial_neighborhood_score = .data[[exact_col]]
      )
    attr(out, "source_column") <- exact_col
    return(out)
  }

  if ("recipient_type" %in% nm) {
    candidate_cols <- c(
      "recipient_neighborhood_score",
      "neighborhood_score",
      "recipient_score",
      "score"
    )
    chosen <- candidate_cols[candidate_cols %in% nm][1]

    if (is.na(chosen) || length(chosen) == 0) {
      chosen <- nm[grepl("score", nm, ignore.case = TRUE)][1]
    }

    if (is.na(chosen) || length(chosen) == 0) {
      stop2("Could not identify an endothelial neighborhood score column.")
    }

    out <- neighborhood_df %>%
      dplyr::filter(recipient_type == "Endothelial") %>%
      dplyr::transmute(
        barcode = .data$barcode,
        section_id = .data$section_id,
        endothelial_neighborhood_score = .data[[chosen]]
      )

    attr(out, "source_column") <- chosen
    return(out)
  }

  stop2(
    "Could not find the independent endothelial neighborhood term. ",
    "Expected 'neighbor_Endothelial_weight'."
  )
}

safe_z <- function(x) {
  x <- as.numeric(x)
  s <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(s) || s <= 0) {
    return(rep(NA_real_, length(x)))
  }
  (x - mean(x, na.rm = TRUE)) / s
}

# ------------------------------------------------------------
# Input files
# ------------------------------------------------------------
ROOT <- locate_root()

spots_file <- file_in_root(
  ROOT,
  "04_results/step12F_GSE234774_independent_lesion_zones_target_program_validation/tables/02_spot_independent_lesion_zones_and_scores.csv.gz"
)

section_manifest_file <- file_in_root(
  ROOT,
  "04_results/step12A2_GSE234774_spatial_metadata_asset_linkage_repair/03_condition_specific_section_manifest.csv"
)

neighborhood_spot_file <- file_in_root(
  ROOT,
  "04_results/step12G_GSE234774_recipient_neighborhood_CellChatDB_spatial_closure/tables/04_spot_recipient_neighborhood_scores.csv.gz",
  handoff_name = "04_spot_recipient_neighborhood_scores.csv.gz",
  must = FALSE
)

if (is.na(neighborhood_spot_file)) {
  stop2(
    "Step12G spot-neighborhood output was not found.\n",
    "Expected:\n",
    "04_results/step12G_GSE234774_recipient_neighborhood_CellChatDB_spatial_closure/",
    "tables/04_spot_recipient_neighborhood_scores.csv.gz"
  )
}

spots <- readr::read_csv(spots_file, show_col_types = FALSE)
section_manifest <- readr::read_csv(section_manifest_file, show_col_types = FALSE)
neighborhood_raw <- readr::read_csv(neighborhood_spot_file, show_col_types = FALSE)

required_cols <- c(
  "barcode",
  "section_id",
  "distance_to_core_NN",
  "target_program_adjusted_z_within_section"
)

missing_cols <- setdiff(required_cols, names(spots))
if (length(missing_cols) > 0) {
  stop2("Step12F spot table missing: ", paste(missing_cols, collapse = ", "))
}

endo <- extract_endothelial_scores(neighborhood_raw)
endo_source_col <- attr(endo, "source_column")

# ------------------------------------------------------------
# Merge frozen spatial signals
# ------------------------------------------------------------
dat <- spots %>%
  dplyr::select(
    tidyselect::any_of(c(
      "barcode",
      "section_id",
      "biological_replicate_id",
      "timepoint",
      "injury_condition",
      "distance_to_core_NN",
      "target_program_adjusted_z_within_section"
    ))
  ) %>%
  dplyr::left_join(
    endo,
    by = c("barcode", "section_id")
  )

# Bring in timepoint/animal only when absent from Step12F.
if (!all(c("timepoint", "biological_replicate_id") %in% names(dat))) {
  dat <- dat %>%
    dplyr::left_join(
      section_manifest %>%
        dplyr::select(
          section_id,
          biological_replicate_id_manifest = biological_replicate_id,
          timepoint_manifest = timepoint
        ),
      by = "section_id"
    )

  if (!"timepoint" %in% names(dat)) {
    dat$timepoint <- dat$timepoint_manifest
  }
  if (!"biological_replicate_id" %in% names(dat)) {
    dat$biological_replicate_id <- dat$biological_replicate_id_manifest
  }
}

dat <- dat %>%
  dplyr::filter(
    timepoint %in% TIMEPOINT_ORDER,
    is.finite(distance_to_core_NN),
    is.finite(target_program_adjusted_z_within_section),
    is.finite(endothelial_neighborhood_score)
  )

if ("injury_condition" %in% names(dat) &&
    any(!is.na(dat$injury_condition))) {
  dat <- dat %>% dplyr::filter(injury_condition == "SCI")
}

if (nrow(dat) == 0) stop2("No valid injured 7 dpi / 2 months spots remain.")

# ------------------------------------------------------------
# Section-wise standardization
# ------------------------------------------------------------
# Target is already a within-section adjusted z score.
# Endothelial neighborhood weight is placed on the same visual scale by
# z-standardizing within each section.
dat <- dat %>%
  dplyr::group_by(section_id) %>%
  dplyr::mutate(
    target_z = as.numeric(target_program_adjusted_z_within_section),
    endothelial_z = safe_z(endothelial_neighborhood_score)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::filter(is.finite(target_z), is.finite(endothelial_z))

# ------------------------------------------------------------
# Distance bins
# ------------------------------------------------------------
breaks_use <- seq(0, MAX_DISTANCE_NN, by = BIN_WIDTH_NN)
if (tail(breaks_use, 1) < MAX_DISTANCE_NN) {
  breaks_use <- c(breaks_use, MAX_DISTANCE_NN)
}

dat <- dat %>%
  dplyr::filter(distance_to_core_NN <= MAX_DISTANCE_NN) %>%
  dplyr::mutate(
    distance_bin = cut(
      distance_to_core_NN,
      breaks = breaks_use,
      include.lowest = TRUE,
      right = FALSE
    )
  )

bin_lookup <- tibble::tibble(
  distance_bin = levels(dat$distance_bin),
  distance_mid = head(breaks_use, -1) + diff(breaks_use) / 2
)

# ------------------------------------------------------------
# Hierarchical aggregation
# ------------------------------------------------------------
# Step 1: spot -> section/bin
section_bin <- dat %>%
  dplyr::group_by(
    timepoint,
    biological_replicate_id,
    section_id,
    distance_bin
  ) %>%
  dplyr::summarise(
    n_spots = dplyr::n(),
    target_z = mean(target_z, na.rm = TRUE),
    endothelial_z = mean(endothelial_z, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::filter(n_spots >= MIN_SPOTS_PER_SECTION_BIN) %>%
  dplyr::left_join(bin_lookup, by = "distance_bin")

# Step 2: section/bin -> biological animal/bin
animal_bin <- section_bin %>%
  dplyr::group_by(
    timepoint,
    biological_replicate_id,
    distance_mid
  ) %>%
  dplyr::summarise(
    n_sections = dplyr::n(),
    target_z = mean(target_z, na.rm = TRUE),
    endothelial_z = mean(endothelial_z, na.rm = TRUE),
    .groups = "drop"
  )

# Step 3: long format for plotting
animal_long <- animal_bin %>%
  tidyr::pivot_longer(
    cols = c(target_z, endothelial_z),
    names_to = "signal",
    values_to = "score_z"
  ) %>%
  dplyr::mutate(
    signal = dplyr::recode(
      signal,
      target_z = "Target program",
      endothelial_z = "Endothelial neighborhood"
    ),
    timepoint = factor(timepoint, levels = TIMEPOINT_ORDER)
  )

# Step 4: biological-animal mean curve
group_long <- animal_long %>%
  dplyr::group_by(timepoint, distance_mid, signal) %>%
  dplyr::summarise(
    n_animals = dplyr::n_distinct(biological_replicate_id),
    mean_z = mean(score_z, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::filter(n_animals >= 2)

# ------------------------------------------------------------
# Plot
# ------------------------------------------------------------
signal_cols <- c(
  "Target program" = TARGET_COLOR,
  "Endothelial neighborhood" = ENDO_COLOR
)

p <- ggplot() +

  # Individual biological animals: thin background trajectories
  geom_line(
    data = animal_long,
    aes(
      x = distance_mid,
      y = score_z,
      colour = signal,
      group = interaction(signal, biological_replicate_id)
    ),
    linewidth = THIN_LINE_WIDTH,
    alpha = THIN_LINE_ALPHA
  ) +

  # Group mean across biological animals
  geom_line(
    data = group_long,
    aes(
      x = distance_mid,
      y = mean_z,
      colour = signal,
      group = signal
    ),
    linewidth = MEAN_LINE_WIDTH,
    lineend = "round"
  ) +

  geom_point(
    data = group_long,
    aes(
      x = distance_mid,
      y = mean_z,
      colour = signal
    ),
    size = MEAN_POINT_SIZE,
    stroke = 0
  ) +

  geom_hline(
    yintercept = 0,
    colour = "#B7B7B7",
    linewidth = 0.42,
    linetype = "dashed"
  ) +

  facet_grid(
    . ~ timepoint,
    labeller = as_labeller(TIMEPOINT_LABELS)
  ) +

  scale_colour_manual(
    values = signal_cols,
    breaks = c("Target program", "Endothelial neighborhood"),
    name = NULL
  ) +

  scale_x_continuous(
    breaks = seq(0, MAX_DISTANCE_NN, by = 10),
    limits = c(0, MAX_DISTANCE_NN),
    expand = expansion(mult = c(0.02, 0.04))
  ) +

  labs(
    tag = NULL,
    x = "Distance from lesion core (NN units)",
    y = "Standardized spatial score (z)"
  ) +

  theme_classic(base_family = BASE_FAMILY) +

  theme(
    plot.tag = element_blank(),

    strip.background = element_blank(),
    strip.text = element_text(
      size = 15.0,
      face = "bold",
      family = BASE_FAMILY,
      margin = margin(b = 4)
    ),

    axis.title.x = element_text(
      size = 15.0,
      face = "plain",
      family = BASE_FAMILY
    ),

    axis.title.y = element_text(
      size = 12.6,
      face = "plain",
      family = BASE_FAMILY,
      margin = margin(r = 3)
    ),

    axis.text = element_text(
      size = 13.0,
      family = BASE_FAMILY
    ),

    axis.line = element_line(linewidth = 0.58),

    legend.position = "top",
    legend.justification = "center",
    legend.direction = "horizontal",
    legend.text = element_text(
      size = 13.0,
      family = BASE_FAMILY
    ),
    legend.key.width = unit(6.5, "mm"),
    legend.margin = margin(0, 0, 1, 0),

    panel.spacing.x = unit(4.5, "mm"),

    plot.margin = margin(2, 10, 3, 22)
  )

# ------------------------------------------------------------
# Outputs
# ------------------------------------------------------------
out_dir <- file.path(ROOT, OUTPUT_DIR_NAME)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_base <- file.path(
  out_dir,
  "Figure5C_target_endothelial_distance_to_core_gradient_v5_ytitle_final"
)

pdf_device <- if (capabilities("cairo")) cairo_pdf else grDevices::pdf

ggsave(
  paste0(out_base, ".pdf"),
  p,
  width = 5.8,
  height = 2.78,
  units = "in",
  device = pdf_device,
  bg = "white"
)

ggsave(
  paste0(out_base, ".png"),
  p,
  width = 5.8,
  height = 2.78,
  units = "in",
  dpi = DPI_OUT,
  bg = "white"
)

if (requireNamespace("svglite", quietly = TRUE)) {
  ggsave(
    paste0(out_base, ".svg"),
    p,
    width = 5.8,
    height = 2.78,
    units = "in",
    bg = "white"
  )
}

readr::write_csv(
  section_bin,
  paste0(out_base, "_section_bin_summary.csv")
)

readr::write_csv(
  animal_bin,
  paste0(out_base, "_animal_bin_summary.csv")
)

readr::write_csv(
  group_long,
  paste0(out_base, "_group_mean_curve.csv")
)

audit <- tibble::tibble(
  item = c(
    "n_spots_used",
    "n_sections_used",
    "n_animals_used",
    "max_distance_NN",
    "bin_width_NN",
    "min_spots_per_section_bin",
    "endothelial_source_column",
    "aggregation",
    "interpretation"
  ),
  value = c(
    nrow(dat),
    dplyr::n_distinct(dat$section_id),
    dplyr::n_distinct(paste(dat$timepoint, dat$biological_replicate_id)),
    MAX_DISTANCE_NN,
    BIN_WIDTH_NN,
    MIN_SPOTS_PER_SECTION_BIN,
    endo_source_col,
    "spot -> section/bin -> animal/bin -> group mean",
    "thin curves are biological animals; thick curves are animal-level means"
  )
)

readr::write_csv(
  audit,
  paste0(out_base, "_audit.csv")
)

message("Figure 5C v5 y-title final-fix script completed.")
message("Output directory: ", out_dir)
