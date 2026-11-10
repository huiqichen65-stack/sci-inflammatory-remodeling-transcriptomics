# ============================================================
# 99_Figure5A_GLIA_7dpi_target_endothelial_spatial_coupling_reference_style.R
# Figure 5A | 7 dpi target-program / endothelial spatial coupling (V3 FINAL)
# ------------------------------------------------------------
# Purpose
#   Build Figure 5A as a three-panel row in the style of the
#   provided reference:
#     left   = target-program spatial map
#     middle = endothelial-neighborhood spatial map
#     right  = target vs endothelial spot-level concordance
#
# Scientific boundary
#   - Uses the frozen Step12F target-program spatial validation output
#     and the frozen Step12G endothelial recipient-neighborhood output.
#   - Displays a representative 7 dpi section.
#   - The right panel is descriptive visual concordance at spot level.
#     It is NOT used as inferential evidence.
#   - No lesion-zone polygons are drawn on the tissue.
# ============================================================

options(stringsAsFactors = FALSE, warn = 1)

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(png)
  library(scales)
  library(grid)
  library(MASS)
  library(patchwork)
})

# ------------------------------------------------------------
# User settings
# ------------------------------------------------------------
ROOT_CANDIDATES <- unique(c(Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE)), "."))

TARGET_TIMEPOINT <- "7days"

# Preferred section from the earlier Figure 5A work.
PREFERRED_SECTION_ID <- "106_C_7days_2"
USE_PREFERRED_SECTION_IF_AVAILABLE <- TRUE

# Optional fallback to auto-select if the preferred section is absent.
AUTO_SELECT_IF_PREFERRED_MISSING <- TRUE

OUTPUT_DIR_NAME <- "09_final_figures/GLIA/Figure5/single_panels/Figure5A"

BASE_FAMILY <- "Arial"
DPI_OUT     <- 600

# Reference-like tissue rendering
BACKGROUND_STRENGTH <- 0.18

# Spatial-spot appearance
SPOT_SIZE   <- 1.38
SPOT_STROKE <- 0.035
SPOT_ALPHA  <- 1.00
SPOT_BORDER <- "#202020"

# Tight crop around all tissue spots
CROP_MARGIN_X_PX <- 20
CROP_MARGIN_Y_PX <- 8

# Target-program palette (same visual family as prior Figure 5A)
TARGET_COLORS <- c(
  "#090909",
  "#21102F",
  "#4A176D",
  "#7A287D",
  "#B63D80",
  "#E65B62",
  "#F58B45",
  "#F6D66A"
)

# Endothelial map palette
ENDO_COLORS <- c(
  "#090909",
  "#221034",
  "#4D1871",
  "#822D80",
  "#BF4781",
  "#E9665E",
  "#F59641",
  "#F6D96E"
)

# Correlation scatter palette
SCATTER_COLORS <- c(
  "#2D2C87",
  "#6A2FB3",
  "#B3368A",
  "#F06D4B",
  "#FDE725"
)

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------
stop2 <- function(...) stop(..., call. = FALSE)

locate_root <- function() {
  for (root in ROOT_CANDIDATES) {
    if (!dir.exists(root)) next
    key1 <- file.path(
      root,
      "04_results/step12F_GSE234774_independent_lesion_zones_target_program_validation/tables/02_spot_independent_lesion_zones_and_scores.csv.gz"
    )
    key2 <- file.path(root, "02_spot_independent_lesion_zones_and_scores.csv.gz")
    if (file.exists(key1) || file.exists(key2)) {
      return(normalizePath(root, winslash = "/", mustWork = TRUE))
    }
  }
  stop2("Could not locate project root. Edit ROOT_CANDIDATES.")
}

file_in_root <- function(root, original_rel, handoff_name = basename(original_rel), must = TRUE) {
  candidates <- c(file.path(root, original_rel), file.path(root, handoff_name))
  hit <- candidates[file.exists(candidates)][1]
  if ((is.na(hit) || length(hit) == 0) && must) {
    stop2("Required file not found:\n  ", original_rel, "\nroot = ", root)
  }
  if (is.na(hit) || length(hit) == 0) return(NA_character_)
  normalizePath(hit, winslash = "/", mustWork = TRUE)
}

first_existing <- function(paths) {
  hit <- paths[file.exists(paths)][1]
  if (is.na(hit) || length(hit) == 0) return(NA_character_)
  normalizePath(hit, winslash = "/", mustWork = TRUE)
}

read_any_table <- function(path) {
  if (is.na(path)) return(NULL)
  if (grepl("\\.rds$", path, ignore.case = TRUE)) {
    obj <- readRDS(path)
    if (is.data.frame(obj)) return(as.data.frame(obj))
    stop2("RDS does not contain a data.frame: ", path)
  }
  readr::read_csv(path, show_col_types = FALSE)
}

read_png_gz <- function(path) {
  con <- gzfile(path, open = "rb")
  on.exit(close(con), add = TRUE)
  raw_vec <- readBin(con, what = "raw", n = 1e8)
  png::readPNG(raw_vec)
}

lighten_image <- function(img, strength = 0.18) {
  strength <- max(0, min(1, strength))
  out <- img
  if (length(dim(out)) != 3 || dim(out)[3] < 3) return(out)
  out[, , 1:3] <- 1 - strength * (1 - out[, , 1:3])
  if (dim(out)[3] == 4) out[, , 4] <- 1
  out
}

read_hires_scalef <- function(path) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)
  txt <- paste(readLines(con, warn = FALSE), collapse = "")
  m <- regexec('"tissue_hires_scalef"\\s*:\\s*([0-9eE+\\.-]+)', txt)
  hit <- regmatches(txt, m)[[1]]
  if (length(hit) < 2) stop2("Could not parse tissue_hires_scalef from: ", path)
  as.numeric(hit[2])
}

nice_pretty_breaks <- function(x, n = 4L) {
  b <- pretty(x, n = n)
  b[b >= min(x, na.rm = TRUE) & b <= max(x, na.rm = TRUE)]
}

choose_representative_section <- function(
  target_section_contrast,
  endothelial_section_contrast,
  preferred_section_id = NULL
) {
  if (!is.null(preferred_section_id) &&
      preferred_section_id %in% target_section_contrast$section_id &&
      preferred_section_id %in% endothelial_section_contrast$section_id) {
    return(preferred_section_id)
  }

  targ <- target_section_contrast %>%
    filter(
      timepoint == TARGET_TIMEPOINT,
      metric == "target_adjusted_within_section_z",
      contrast == "core_minus_remote"
    ) %>%
    select(section_id, biological_replicate_id, target_core_effect = zone_minus_remote)

  endo <- endothelial_section_contrast %>%
    filter(timepoint == TARGET_TIMEPOINT)

  # Keep the most plausible endothelial effect column if present.
  if (!"zone_minus_remote" %in% names(endo)) {
    stop2("Could not find zone_minus_remote in Step12G endothelial section-contrast table.")
  }

  # Use explicit Endothelial recipient rows.
  if ("recipient_type" %in% names(endo)) {
    endo <- endo %>% filter(recipient_type == "Endothelial")
  } else {
    stop2("Could not identify recipient_type column in Step12G endothelial table.")
  }

  # Prefer contrast == core_minus_remote when available.
  if ("contrast" %in% names(endo)) {
    endo <- endo %>% filter(contrast == "core_minus_remote")
  }

  endo <- endo %>%
    select(section_id, endothelial_core_effect = zone_minus_remote)

  joint <- inner_join(targ, endo, by = "section_id") %>%
    mutate(
      target_rank = rank(target_core_effect, ties.method = "first"),
      endothelial_rank = rank(endothelial_core_effect, ties.method = "first"),
      rank_sum = target_rank + endothelial_rank
    ) %>%
    arrange(desc(rank_sum))

  if (nrow(joint) == 0) {
    stop2("Failed to auto-select a representative 7 dpi section.")
  }

  # Choose a section near the upper half but not necessarily the single most extreme.
  joint <- joint %>%
    mutate(rank_position = row_number(),
           distance_to_upper_mid = abs(rank_position - ceiling(n() * 0.35))) %>%
    arrange(distance_to_upper_mid, desc(rank_sum), section_id)

  joint$section_id[1]
}

extract_endothelial_scores <- function(neighborhood_df) {
  nm <- names(neighborhood_df)

  # Case 1: long format
  long_score_candidates <- c(
    "recipient_neighborhood_score",
    "neighborhood_score",
    "score",
    "weighted_neighbor_score",
    "recipient_score"
  )

  if ("recipient_type" %in% nm) {
    score_col <- long_score_candidates[long_score_candidates %in% nm][1]
    if (length(score_col) == 0 || is.na(score_col)) {
      # fallback: any column containing "score"
      score_col <- nm[grepl("score", nm, ignore.case = TRUE)][1]
    }
    if (is.na(score_col) || length(score_col) == 0) {
      stop2("Could not identify the endothelial neighborhood score column in long-format Step12G data.")
    }
    out <- neighborhood_df %>%
      filter(recipient_type == "Endothelial") %>%
      transmute(
        barcode   = .data$barcode,
        section_id = .data$section_id,
        endothelial_neighborhood_score = .data[[score_col]]
      )
    attr(out, "source_column") <- score_col
    return(out)
  }

  # Case 2: wide format
  # IMPORTANT: use the independent recipient-neighborhood abundance term.
  # Do NOT use target_neighbor_Endothelial_colocalization here because that
  # already contains the target score and would make the scatter partly circular.
  exact_col <- "neighbor_Endothelial_weight"

  if (exact_col %in% nm) {
    chosen <- exact_col
  } else {
    wide_cols <- nm[
      grepl("^neighbor_.*endothelial.*weight$", nm, ignore.case = TRUE)
    ]
    if (length(wide_cols) == 0) {
      stop2(
        "Could not find the independent Endothelial neighborhood-weight column. ",
        "Expected 'neighbor_Endothelial_weight'."
      )
    }
    chosen <- wide_cols[1]
  }

  out <- neighborhood_df %>%
    transmute(
      barcode = .data$barcode,
      section_id = .data$section_id,
      endothelial_neighborhood_score = .data[[chosen]]
    )
  attr(out, "source_column") <- chosen
  out
}

compute_point_density <- function(x, y, n = 200L) {
  ok <- is.finite(x) & is.finite(y)
  z <- rep(NA_real_, length(x))
  if (sum(ok) < 10L) return(z)

  kd <- MASS::kde2d(x[ok], y[ok], n = n)

  xi <- findInterval(x[ok], kd$x, all.inside = TRUE)
  yi <- findInterval(y[ok], kd$y, all.inside = TRUE)
  z_ok <- kd$z[cbind(xi, yi)]
  z[ok] <- z_ok
  z
}

make_spatial_map <- function(df, fill_var, title_text, img_raster, img_w, img_h,
                             xlim_use, ylim_use, fill_colors) {
  fill_vals <- df[[fill_var]]
  fill_lims <- quantile(fill_vals, probs = c(0.02, 0.98), na.rm = TRUE)
  fill_lims <- unname(fill_lims)

  # Guarantee a non-degenerate range
  if (!all(is.finite(fill_lims)) || diff(fill_lims) <= 0) {
    fill_lims <- range(fill_vals, na.rm = TRUE)
  }
  if (!all(is.finite(fill_lims)) || diff(fill_lims) <= 0) {
    fill_lims <- c(0, 1)
  }

  breaks_use <- nice_pretty_breaks(fill_lims, n = 4)

  ggplot(df, aes(x = x_hires, y = y_hires, fill = .data[[fill_var]])) +
    annotation_raster(
      as.raster(img_raster),
      xmin = 0,
      xmax = img_w,
      ymin = img_h,
      ymax = 0
    ) +
    geom_point(
      shape = 21,
      size = SPOT_SIZE,
      colour = SPOT_BORDER,
      stroke = SPOT_STROKE,
      alpha = SPOT_ALPHA
    ) +
    scale_fill_gradientn(
      colours = fill_colors,
      limits = fill_lims,
      oob = scales::squish,
      breaks = breaks_use,
      labels = label_number(accuracy = 0.01)
    ) +
    coord_equal(
      xlim = xlim_use,
      ylim = ylim_use,
      expand = FALSE,
      clip = "on"
    ) +
    guides(
      fill = guide_colorbar(
        title = NULL,
        barheight = unit(31, "mm"),
        barwidth = unit(3.8, "mm"),
        frame.colour = "#BDBDBD",
        ticks.colour = "#4D4D4D"
      )
    ) +
    labs(title = title_text) +
    theme_void(base_family = BASE_FAMILY) +
    theme(
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.title = element_text(
        size = 9.4,
        face = "plain",
        hjust = 0.5,
        margin = margin(b = 4)
      ),
      legend.position = "right",
      legend.box.spacing = unit(0.2, "mm"),
      legend.margin = margin(0, 0, 0, 0),
      legend.text = element_text(size = 7.6),
      plot.margin = margin(1, 0.5, 1, 0.5)
    )
}

# ------------------------------------------------------------
# Locate inputs
# ------------------------------------------------------------
ROOT <- locate_root()

spots_file <- file_in_root(
  ROOT,
  "04_results/step12F_GSE234774_independent_lesion_zones_target_program_validation/tables/02_spot_independent_lesion_zones_and_scores.csv.gz"
)

target_section_contrast_file <- file_in_root(
  ROOT,
  "04_results/step12F_GSE234774_independent_lesion_zones_target_program_validation/tables/04_section_zone_vs_remote_contrasts.csv"
)

section_manifest_file <- file_in_root(
  ROOT,
  "04_results/step12A2_GSE234774_spatial_metadata_asset_linkage_repair/03_condition_specific_section_manifest.csv"
)

asset_manifest_file <- file_in_root(
  ROOT,
  "04_results/step12A2_GSE234774_spatial_metadata_asset_linkage_repair/02_2d_capture_area_asset_manifest.csv"
)

# Step12G files: prefer original project output; handoff may not contain them.
neighborhood_spot_file <- file_in_root(
  ROOT,
  "04_results/step12G_GSE234774_recipient_neighborhood_CellChatDB_spatial_closure/tables/04_spot_recipient_neighborhood_scores.csv.gz",
  handoff_name = "04_spot_recipient_neighborhood_scores.csv.gz",
  must = FALSE
)

neighborhood_section_file <- file_in_root(
  ROOT,
  "04_results/step12G_GSE234774_recipient_neighborhood_CellChatDB_spatial_closure/tables/05_section_recipient_neighborhood_zone_contrasts.csv",
  handoff_name = "05_section_recipient_neighborhood_zone_contrasts.csv",
  must = FALSE
)

if (is.na(neighborhood_spot_file)) {
  stop2(
    "Step12G spot-neighborhood file was not found.\n",
    "Required original file:\n",
    "  04_results/step12G_GSE234774_recipient_neighborhood_CellChatDB_spatial_closure/tables/04_spot_recipient_neighborhood_scores.csv.gz\n",
    "Please make sure your xjz project contains the frozen Step12G output."
  )
}

# ------------------------------------------------------------
# Read tables
# ------------------------------------------------------------
spots <- read_any_table(spots_file)
target_section_contrast <- read_any_table(target_section_contrast_file)
section_manifest <- read_any_table(section_manifest_file)
asset_manifest <- read_any_table(asset_manifest_file)
neighborhood_spot_raw <- read_any_table(neighborhood_spot_file)
neighborhood_section_contrast <- if (!is.na(neighborhood_section_file)) read_any_table(neighborhood_section_file) else NULL

required_spot_cols <- c(
  "barcode", "section_id", "coordinate_x_fullres", "coordinate_y_fullres",
  "target_program_adjusted_z_within_section"
)
missing_spot_cols <- setdiff(required_spot_cols, names(spots))
if (length(missing_spot_cols) > 0) {
  stop2("Step12F spot table is missing required columns: ",
        paste(missing_spot_cols, collapse = ", "))
}

# ------------------------------------------------------------
# Choose section
# ------------------------------------------------------------
section_ids_7d <- section_manifest %>%
  filter(timepoint == TARGET_TIMEPOINT) %>%
  pull(section_id)

SECTION_ID <- NA_character_

if (USE_PREFERRED_SECTION_IF_AVAILABLE &&
    PREFERRED_SECTION_ID %in% section_ids_7d) {
  SECTION_ID <- PREFERRED_SECTION_ID
} else if (AUTO_SELECT_IF_PREFERRED_MISSING &&
           !is.null(neighborhood_section_contrast)) {
  SECTION_ID <- choose_representative_section(
    target_section_contrast = target_section_contrast,
    endothelial_section_contrast = neighborhood_section_contrast,
    preferred_section_id = NULL
  )
} else {
  if (length(section_ids_7d) == 0) stop2("No 7 dpi sections found.")
  SECTION_ID <- section_ids_7d[1]
}

sec_meta <- section_manifest %>% filter(section_id == !!SECTION_ID)
if (nrow(sec_meta) != 1) stop2("Expected one section_manifest row for ", SECTION_ID)

CAPTURE_AREA_ID <- sec_meta$capture_area_id[1]
BIO_REP_ID      <- sec_meta$biological_replicate_id[1]
TIMEPOINT       <- sec_meta$timepoint[1]

# ------------------------------------------------------------
# Merge Step12F target scores with Step12G endothelial scores
# ------------------------------------------------------------
endo_df <- extract_endothelial_scores(neighborhood_spot_raw)
endo_source_col <- attr(endo_df, "source_column")

section_df <- spots %>%
  filter(section_id == !!SECTION_ID) %>%
  left_join(endo_df, by = c("barcode", "section_id"))

if (!"endothelial_neighborhood_score" %in% names(section_df)) {
  stop2("Failed to merge endothelial scores into the selected section.")
}

n_complete <- sum(
  is.finite(section_df$target_program_adjusted_z_within_section) &
    is.finite(section_df$endothelial_neighborhood_score)
)

if (n_complete < 50) {
  stop2("Too few complete spots after merging target and endothelial scores: n = ", n_complete)
}

# ------------------------------------------------------------
# Locate image + scale factor
# ------------------------------------------------------------
asset_row <- asset_manifest %>% filter(capture_area_id == !!CAPTURE_AREA_ID)
if (nrow(asset_row) != 1) stop2("Expected exactly one asset row for capture_area_id = ", CAPTURE_AREA_ID)

asset_dir <- file.path(
  ROOT,
  "04_results/step12A2_GSE234774_spatial_metadata_asset_linkage_repair/audit_extracted_2d_assets"
)

hires_file <- first_existing(c(
  file.path(asset_dir, asset_row$image_hires_member[1]),
  file.path(ROOT, asset_row$image_hires_member[1])
))

scale_file <- first_existing(c(
  file.path(asset_dir, asset_row$scalefactor_member[1]),
  file.path(ROOT, asset_row$scalefactor_member[1])
))

if (is.na(hires_file) || is.na(scale_file)) {
  stop2("Could not locate hires tissue image or scalefactor file for capture_area_id = ", CAPTURE_AREA_ID)
}

img_original <- read_png_gz(hires_file)
img_light    <- lighten_image(img_original, BACKGROUND_STRENGTH)
hires_scalef <- read_hires_scalef(scale_file)

img_h <- dim(img_original)[1]
img_w <- dim(img_original)[2]

section_df <- section_df %>%
  mutate(
    x_hires = coordinate_x_fullres * hires_scalef,
    y_hires = coordinate_y_fullres * hires_scalef
  )

in_bounds <- with(
  section_df,
  x_hires >= 0 & x_hires <= img_w & y_hires >= 0 & y_hires <= img_h
)

if (mean(in_bounds, na.rm = TRUE) < 0.95) {
  stop2("Coordinate alignment check failed for the selected section.")
}

xmin <- max(0, min(section_df$x_hires, na.rm = TRUE) - CROP_MARGIN_X_PX)
xmax <- min(img_w, max(section_df$x_hires, na.rm = TRUE) + CROP_MARGIN_X_PX)
ymin <- max(0, min(section_df$y_hires, na.rm = TRUE) - CROP_MARGIN_Y_PX)
ymax <- min(img_h, max(section_df$y_hires, na.rm = TRUE) + CROP_MARGIN_Y_PX)

xlim_use <- c(xmin, xmax)
ylim_use <- c(ymax, ymin)

# ------------------------------------------------------------
# Build maps
# ------------------------------------------------------------
p_target <- make_spatial_map(
  df = section_df,
  fill_var = "target_program_adjusted_z_within_section",
  title_text = "Signature score of\nTarget program",
  img_raster = img_light,
  img_w = img_w,
  img_h = img_h,
  xlim_use = xlim_use,
  ylim_use = ylim_use,
  fill_colors = TARGET_COLORS
)

p_endo <- make_spatial_map(
  df = section_df,
  fill_var = "endothelial_neighborhood_score",
  title_text = "Signature score of\nEndothelial neighborhood",
  img_raster = img_light,
  img_w = img_w,
  img_h = img_h,
  xlim_use = xlim_use,
  ylim_use = ylim_use,
  fill_colors = ENDO_COLORS
)

# ------------------------------------------------------------
# Correlation panel
# ------------------------------------------------------------
corr_df <- section_df %>%
  filter(
    is.finite(target_program_adjusted_z_within_section),
    is.finite(endothelial_neighborhood_score)
  ) %>%
  transmute(
    target_score = target_program_adjusted_z_within_section,
    endothelial_score = endothelial_neighborhood_score
  )

rho <- suppressWarnings(
  cor(
    corr_df$target_score,
    corr_df$endothelial_score,
    method = "spearman",
    use = "complete.obs"
  )
)

corr_df <- corr_df %>%
  mutate(local_density = compute_point_density(target_score, endothelial_score))

lab_text <- paste0("\u03C1 = ", formatC(rho, format = "f", digits = 2))

x_rng <- range(corr_df$target_score, na.rm = TRUE)
y_rng <- range(corr_df$endothelial_score, na.rm = TRUE)
x_annot <- x_rng[1] + 0.04 * diff(x_rng)
y_annot <- y_rng[2] - 0.06 * diff(y_rng)

p_corr <- ggplot(corr_df, aes(x = target_score, y = endothelial_score)) +
  geom_point(aes(colour = local_density), size = 1.20, alpha = 0.90) +
  scale_colour_gradientn(colours = SCATTER_COLORS, guide = "none") +
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    colour = "black",
    fill = "grey75",
    linewidth = 0.78
  ) +
  annotate(
    "text",
    x = x_rng[1] + 0.07 * diff(x_rng),
    y = y_rng[2] - 0.10 * diff(y_rng),
    label = lab_text,
    hjust = 0,
    vjust = 1,
    size = 6.2 / .pt,
    family = BASE_FAMILY,
    fontface = "bold"
  ) +
  labs(
    x = "Signature score of\nTarget program",
    y = "Signature score of endothelial\nneighborhood"
  ) +
  theme_classic(base_family = BASE_FAMILY) +
  theme(
    axis.title = element_text(size = 9.4, face = "bold"),
    axis.text = element_text(size = 8.2),
    axis.line = element_line(linewidth = 0.6),
    plot.margin = margin(8, 5, 10, 8)
  )

# ------------------------------------------------------------
# Combine
# ------------------------------------------------------------
# Put the single Figure-panel label "A" and the vertical timepoint label
# in a dedicated narrow strip. This avoids patchwork auto-tagging the
# three internal subplots as B/C/D.
row_lab <- ggplot() +
  annotate(
    "text",
    x = 0.10,
    y = 0.94,
    label = "A",
    size = 17 / .pt,
    family = BASE_FAMILY,
    fontface = "bold",
    hjust = 0
  ) +
  annotate(
    "text",
    x = 0.76,
    y = 0.50,
    label = "7 dpi",
    angle = 90,
    size = 10.0 / .pt,
    family = BASE_FAMILY,
    fontface = "bold"
  ) +
  theme_void() +
  xlim(0, 1) + ylim(0, 1)

final_plot <- row_lab + p_target + p_endo + p_corr +
  plot_layout(widths = c(0.11, 1.08, 1.08, 1.10))

# ------------------------------------------------------------
# Outputs
# ------------------------------------------------------------
out_dir <- file.path(ROOT, OUTPUT_DIR_NAME)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_base <- file.path(
  out_dir,
  "Figure5A_7dpi_target_endothelial_spatial_coupling_reference_style_v3_final"
)

pdf_device <- if (capabilities("cairo")) cairo_pdf else grDevices::pdf

ggsave(
  paste0(out_base, ".pdf"),
  final_plot,
  width = 10.2,
  height = 3.15,
  units = "in",
  device = pdf_device,
  bg = "white"
)

ggsave(
  paste0(out_base, ".png"),
  final_plot,
  width = 10.2,
  height = 3.15,
  units = "in",
  dpi = DPI_OUT,
  bg = "white"
)

if (requireNamespace("svglite", quietly = TRUE)) {
  ggsave(
    paste0(out_base, ".svg"),
    final_plot,
    width = 10.2,
    height = 3.15,
    units = "in",
    bg = "white"
  )
}

audit <- tibble::tibble(
  item = c(
    "section_id",
    "capture_area_id",
    "biological_replicate_id",
    "timepoint",
    "n_spots_total",
    "n_spots_complete_for_correlation",
    "fraction_spots_in_image",
    "hires_scalef",
    "step12G_endothelial_source_column",
    "spearman_r_visual_only",
    "background_strength",
    "spot_size",
    "spot_stroke",
    "display_concept"
  ),
  value = c(
    SECTION_ID,
    CAPTURE_AREA_ID,
    BIO_REP_ID,
    TIMEPOINT,
    nrow(section_df),
    n_complete,
    round(mean(in_bounds, na.rm = TRUE), 5),
    signif(hires_scalef, 7),
    endo_source_col,
    signif(rho, 5),
    BACKGROUND_STRENGTH,
    SPOT_SIZE,
    SPOT_STROKE,
    "reference-style triad: target map + endothelial map + descriptive concordance"
  )
)

readr::write_csv(audit, paste0(out_base, "_audit.csv"))

plot_data_cols <- intersect(
  c(
    "barcode",
    "capture_area_id",
    "section_id",
    "biological_replicate_id",
    "timepoint",
    "coordinate_x_fullres",
    "coordinate_y_fullres",
    "x_hires",
    "y_hires",
    "lesion_zone",
    "distance_to_core_NN",
    "Microglia_RCTD_weight",
    "target_program_adjusted_z_within_section",
    "endothelial_neighborhood_score"
  ),
  names(section_df)
)

plot_data_used <- section_df[, plot_data_cols, drop = FALSE]
readr::write_csv(plot_data_used, paste0(out_base, "_plot_data_used.csv"))

message("Figure 5A script completed.")
message("Selected section: ", SECTION_ID)
message("Output directory: ", out_dir)
