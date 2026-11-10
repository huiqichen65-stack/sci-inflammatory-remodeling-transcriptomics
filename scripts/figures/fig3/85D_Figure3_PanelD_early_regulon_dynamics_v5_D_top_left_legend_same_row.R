# ==============================================================================
# 85D_Figure3_PanelD_early_regulon_dynamics_v5_D_top_left_legend_same_row.R
#
# Figure 3D | Dynamics of the earliest Target-associated regulons
#
# PURPOSE
# -------
# Draw 2x2 branch-progress trajectory panels for the earliest data-derived
# Target-associated regulons identified in Figure 3C.
#
# SCIENTIFIC LOGIC
# ----------------
# Figure 3C first shows, from frozen Step7F timing outputs, that the earliest
# qualified-onset group among Target-associated candidates is:
#   Nfkb1, Nfe2l2, Relb, Nfkb2
#
# Figure 3D now visualizes how these four regulon activities change along branch
# progress:
#   - Shared prebranch (negative progress to 0)
#   - Alternative branch (0 to 1)
#   - Target program (0 to 1)
#
# DATA SOURCES
# ------------
# Uses only frozen Step7F outputs:
#   1) GSE162610_step7F_regulon_fate_timing_results.rds
#      -> regulon_fit_curves
#   2) 04_key_candidate_fate_timing.tsv
#      -> identify earliest group
#   3) 07_key_analyte_sample_balanced_bins.tsv.gz
#      -> sample-balanced means / SD for plotting points + error bars
#
# IMPORTANT INTERPRETATION
# ------------------------
# - The fitted curves are frozen Step7F model predictions.
# - The points/error bars are sample-balanced bin summaries.
# - The vertical dashed line at x = 0 marks the branch origin.
# - This panel shows branch-associated dynamics, not causality.
# - We intentionally use "Branch progress" display language and do NOT show the
#   internal Y_29 label in the figure.
#
# OUTPUT
# ------
# <PROJECT_ROOT>/09_final_figures/GLIA/Figure3D_v5/
# ==============================================================================

rm(list = ls())
gc()
options(stringsAsFactors = FALSE, scipen = 999)

# ------------------------------------------------------------------------------
# 0. USER SETTINGS
# ------------------------------------------------------------------------------
PROJECT_ROOT <- Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
OUT_DIR <- file.path(
  PROJECT_ROOT,
  "09_final_figures", "GLIA", "Figure3D_v5"
)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

FONT_FAMILY <- "Arial"

FIG_WIDTH_MM  <- 170
FIG_HEIGHT_MM <- 126

COL_TEXT   <- "#222222"
COL_GRID   <- "#E8E8E8"
COL_SHARED <- "#7A6AA6"
COL_TARGET <- "#B64B45"
COL_OTHER  <- "#7F8C8D"

STATE_DISPLAY <- c(
  "Shared_prebranch" = "Shared prebranch",
  "Arm_2"            = "Target program",
  "Arm_1"            = "Alternative branch"
)

STATE_COLORS <- c(
  "Shared prebranch" = COL_SHARED,
  "Target program"   = COL_TARGET,
  "Alternative branch"     = COL_OTHER
)

# ------------------------------------------------------------------------------
# 1. PACKAGES
# ------------------------------------------------------------------------------
required_pkgs <- c(
  "data.table",
  "ggplot2",
  "patchwork",
  "svglite",
  "ragg"
)

missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_pkgs) > 0L) {
  stop(
    "Missing required R package(s): ",
    paste(missing_pkgs, collapse = ", "),
    "\nInstall them first, then rerun.",
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

# ------------------------------------------------------------------------------
# 2. INPUTS
# ------------------------------------------------------------------------------
STEP7F_RDS <- file.path(
  PROJECT_ROOT,
  "04_results",
  "step7F_GSE162610_regulon_fate_timing",
  "RDS",
  "GSE162610_step7F_regulon_fate_timing_results.rds"
)

STEP7F_TIMING_TSV <- file.path(
  PROJECT_ROOT,
  "04_results",
  "step7F_GSE162610_regulon_fate_timing",
  "tables",
  "04_key_candidate_fate_timing.tsv"
)

STEP7F_BINS_TSV <- file.path(
  PROJECT_ROOT,
  "04_results",
  "step7F_GSE162610_regulon_fate_timing",
  "tables",
  "07_key_analyte_sample_balanced_bins.tsv.gz"
)

input_audit <- data.frame(
  file = c(
    "Step7F frozen RDS",
    "Step7F key candidate timing table",
    "Step7F sample-balanced bins"
  ),
  path = c(
    STEP7F_RDS,
    STEP7F_TIMING_TSV,
    STEP7F_BINS_TSV
  ),
  exists = c(
    file.exists(STEP7F_RDS),
    file.exists(STEP7F_TIMING_TSV),
    file.exists(STEP7F_BINS_TSV)
  ),
  size_bytes = c(
    if (file.exists(STEP7F_RDS)) file.info(STEP7F_RDS)$size else NA_real_,
    if (file.exists(STEP7F_TIMING_TSV)) file.info(STEP7F_TIMING_TSV)$size else NA_real_,
    if (file.exists(STEP7F_BINS_TSV)) file.info(STEP7F_BINS_TSV)$size else NA_real_
  ),
  stringsAsFactors = FALSE
)

write.csv(
  input_audit,
  file.path(OUT_DIR, "00_input_audit.csv"),
  row.names = FALSE
)

if (!all(input_audit$exists)) {
  stop(
    "One or more Figure 3D inputs are missing. Inspect 00_input_audit.csv.",
    call. = FALSE
  )
}

step7f <- readRDS(STEP7F_RDS)
timing_df <- fread(STEP7F_TIMING_TSV, data.table = FALSE)
bins_df <- fread(STEP7F_BINS_TSV, data.table = FALSE)

# ------------------------------------------------------------------------------
# 3. SCHEMA AUDIT
# ------------------------------------------------------------------------------
rds_required <- c("regulon_fit_curves")
timing_required <- c(
  "regulon",
  "TF",
  "model_success",
  "onset_progress_qualified",
  "priority_group",
  "direction_matches_RSS",
  "seven_dpi_direction_consistent",
  "RSS_preferred_state",
  "preferred_RSS"
)
bins_required <- c(
  "state",
  "bin",
  "axis_value",
  "balanced_mean",
  "balanced_sd",
  "n_samples",
  "total_cells",
  "analyte",
  "analyte_type",
  "axis"
)

schema_audit <- rbind(
  data.frame(
    source = "Step7F RDS",
    item = rds_required,
    present = rds_required %in% names(step7f),
    stringsAsFactors = FALSE
  ),
  data.frame(
    source = "04_key_candidate_fate_timing.tsv",
    item = timing_required,
    present = timing_required %in% colnames(timing_df),
    stringsAsFactors = FALSE
  ),
  data.frame(
    source = "07_key_analyte_sample_balanced_bins.tsv.gz",
    item = bins_required,
    present = bins_required %in% colnames(bins_df),
    stringsAsFactors = FALSE
  )
)

write.csv(
  schema_audit,
  file.path(OUT_DIR, "01_schema_audit.csv"),
  row.names = FALSE
)

if (!all(schema_audit$present)) {
  stop(
    "Figure 3D schema audit failed. Inspect 01_schema_audit.csv.",
    call. = FALSE
  )
}

if (!is.list(step7f$regulon_fit_curves)) {
  stop("step7f$regulon_fit_curves is not a list.", call. = FALSE)
}

# ------------------------------------------------------------------------------
# 4. IDENTIFY THE EARLIEST TARGET-ASSOCIATED GROUP (from frozen data)
# ------------------------------------------------------------------------------
candidate_df <- timing_df[
  timing_df$priority_group == "Arm2_priority" &
    timing_df$model_success %in% TRUE &
    !is.na(timing_df$onset_progress_qualified) &
    timing_df$direction_matches_RSS %in% TRUE &
    timing_df$seven_dpi_direction_consistent %in% TRUE &
    timing_df$RSS_preferred_state == "Arm_2",
  ,
  drop = FALSE
]

if (nrow(candidate_df) == 0L) {
  stop("No frozen Target-associated timing candidates remain.", call. = FALSE)
}

candidate_df$onset_progress_qualified <- as.numeric(candidate_df$onset_progress_qualified)
candidate_df$preferred_RSS <- as.numeric(candidate_df$preferred_RSS)

min_onset <- min(candidate_df$onset_progress_qualified, na.rm = TRUE)
earliest_df <- candidate_df[
  abs(candidate_df$onset_progress_qualified - min_onset) < 1e-12,
  ,
  drop = FALSE
]

# Order the four panels by preferred RSS, highest first.
earliest_df <- earliest_df[
  order(-earliest_df$preferred_RSS, earliest_df$TF),
  ,
  drop = FALSE
]

if (nrow(earliest_df) != 4L) {
  warning(
    "Expected 4 earliest regulons from Figure 3C logic; recovered ",
    nrow(earliest_df),
    ". The script will continue with the data-derived count."
  )
}

panel_regulons <- earliest_df$regulon
panel_tfs <- earliest_df$TF

panel_map <- data.frame(
  regulon = panel_regulons,
  TF = panel_tfs,
  onset_progress_qualified = earliest_df$onset_progress_qualified,
  preferred_RSS = earliest_df$preferred_RSS,
  stringsAsFactors = FALSE
)

write.csv(
  panel_map,
  file.path(OUT_DIR, "02_Figure3D_panel_regulons.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 5. BUILD FROZEN CURVE DATA
# ------------------------------------------------------------------------------
curve_parts <- list()

for (ii in seq_len(nrow(panel_map))) {

  reg <- panel_map$regulon[ii]
  tf  <- panel_map$TF[ii]

  fit_result <- step7f$regulon_fit_curves[[reg]]

  if (is.null(fit_result) || !isTRUE(fit_result$success)) {
    stop(
      "Frozen fit result missing or unsuccessful for regulon: ",
      reg,
      call. = FALSE
    )
  }

  required_fit_fields <- c(
    "pre_grid", "post_grid",
    "pre_prediction", "arm1_prediction", "arm2_prediction"
  )

  if (!all(required_fit_fields %in% names(fit_result))) {
    stop(
      "Frozen fit result for ", reg,
      " is missing one or more required fields: ",
      paste(setdiff(required_fit_fields, names(fit_result)), collapse = ", "),
      call. = FALSE
    )
  }

  curve_parts[[length(curve_parts) + 1L]] <- data.frame(
    TF = tf,
    regulon = reg,
    state = "Shared_prebranch",
    state_display = STATE_DISPLAY[["Shared_prebranch"]],
    branch_progress = as.numeric(fit_result$pre_grid),
    fitted_value = as.numeric(fit_result$pre_prediction),
    stringsAsFactors = FALSE
  )

  curve_parts[[length(curve_parts) + 1L]] <- data.frame(
    TF = tf,
    regulon = reg,
    state = "Arm_1",
    state_display = STATE_DISPLAY[["Arm_1"]],
    branch_progress = as.numeric(fit_result$post_grid),
    fitted_value = as.numeric(fit_result$arm1_prediction),
    stringsAsFactors = FALSE
  )

  curve_parts[[length(curve_parts) + 1L]] <- data.frame(
    TF = tf,
    regulon = reg,
    state = "Arm_2",
    state_display = STATE_DISPLAY[["Arm_2"]],
    branch_progress = as.numeric(fit_result$post_grid),
    fitted_value = as.numeric(fit_result$arm2_prediction),
    stringsAsFactors = FALSE
  )
}

curve_df <- do.call(rbind, curve_parts)

curve_df$state_display <- factor(
  curve_df$state_display,
  levels = c("Shared prebranch", "Target program", "Alternative branch")
)

write.csv(
  curve_df,
  file.path(OUT_DIR, "03_Figure3D_curve_data.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 6. SAMPLE-BALANCED BIN SUMMARIES
# ------------------------------------------------------------------------------
bin_map <- setNames(panel_map$TF, panel_map$regulon)

plot_bins <- bins_df[
  bins_df$axis == "Y29_progress" &
    bins_df$analyte_type == "regulon_AUC" &
    bins_df$analyte %in% names(bin_map),
  ,
  drop = FALSE
]

if (nrow(plot_bins) == 0L) {
  stop("No sample-balanced bin rows found for panel regulons.", call. = FALSE)
}

plot_bins$TF <- unname(bin_map[as.character(plot_bins$analyte)])
plot_bins$state_display <- unname(STATE_DISPLAY[as.character(plot_bins$state)])
plot_bins$axis_value <- as.numeric(plot_bins$axis_value)
plot_bins$balanced_mean <- as.numeric(plot_bins$balanced_mean)
plot_bins$balanced_sd <- as.numeric(plot_bins$balanced_sd)
plot_bins$n_samples <- as.numeric(plot_bins$n_samples)

plot_bins$ymin <- ifelse(
  is.finite(plot_bins$balanced_sd) & plot_bins$n_samples >= 2,
  plot_bins$balanced_mean - plot_bins$balanced_sd,
  NA_real_
)

plot_bins$ymax <- ifelse(
  is.finite(plot_bins$balanced_sd) & plot_bins$n_samples >= 2,
  plot_bins$balanced_mean + plot_bins$balanced_sd,
  NA_real_
)

plot_bins$state_display <- factor(
  plot_bins$state_display,
  levels = c("Shared prebranch", "Target program", "Alternative branch")
)

write.csv(
  plot_bins,
  file.path(OUT_DIR, "04_Figure3D_sample_balanced_bins.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 7. GLOBAL AXIS RANGE
# ------------------------------------------------------------------------------
y_all <- c(curve_df$fitted_value, plot_bins$balanced_mean, plot_bins$ymin, plot_bins$ymax)
y_all <- y_all[is.finite(y_all)]

y_min <- min(y_all, na.rm = TRUE)
y_max <- max(y_all, na.rm = TRUE)
y_pad <- max((y_max - y_min) * 0.08, 0.003)

y_limits <- c(y_min - y_pad, y_max + y_pad)
x_limits <- c(-1, 1)

# ------------------------------------------------------------------------------
# 8. PANEL FUNCTION
# ------------------------------------------------------------------------------
build_panel <- function(tf_name, show_legend = FALSE) {

  reg_name <- panel_map$regulon[panel_map$TF == tf_name][1]
  onset_val <- panel_map$onset_progress_qualified[panel_map$TF == tf_name][1]

  d_curve <- curve_df[curve_df$TF == tf_name, , drop = FALSE]
  d_bins  <- plot_bins[plot_bins$TF == tf_name, , drop = FALSE]

  p <- ggplot() +
    geom_vline(
      xintercept = 0,
      linewidth = 0.40,
      colour = COL_TEXT,
      linetype = 2
    ) +
    geom_errorbar(
      data = d_bins[is.finite(d_bins$ymin) & is.finite(d_bins$ymax), , drop = FALSE],
      aes(
        x = axis_value,
        ymin = ymin,
        ymax = ymax,
        colour = state_display
      ),
      width = 0.025,
      linewidth = 0.32,
      alpha = 0.90
    ) +
    geom_point(
      data = d_bins,
      aes(
        x = axis_value,
        y = balanced_mean,
        fill = state_display
      ),
      shape = 21,
      size = 2.05,
      stroke = 0.38,
      colour = "white",
      alpha = 0.98
    ) +
    geom_line(
      data = d_curve,
      aes(
        x = branch_progress,
        y = fitted_value,
        colour = state_display
      ),
      linewidth = 1.05,
      lineend = "round"
    ) +
    scale_colour_manual(
      values = STATE_COLORS,
      name = NULL
    ) +
    scale_fill_manual(
      values = STATE_COLORS,
      name = NULL
    ) +
    scale_x_continuous(
      limits = x_limits,
      breaks = c(-1, -0.5, 0, 0.5, 1)
    ) +
    coord_cartesian(
      ylim = y_limits,
      clip = "off"
    ) +
    labs(
      title = tf_name,
      x = "Branch progress",
      y = "Regulon activity (AUC)"
    ) +
    theme_classic(
      base_family = FONT_FAMILY,
      base_size = 8
    ) +
    theme(
      plot.title = element_text(
        size = 9,
        face = "bold",
        colour = COL_TEXT,
        hjust = 0.5,
        margin = margin(b = 3)
      ),
      axis.title.x = element_text(
        size = 8,
        colour = COL_TEXT,
        margin = margin(t = 3)
      ),
      axis.title.y = element_text(
        size = 8,
        colour = COL_TEXT,
        margin = margin(r = 4)
      ),
      axis.text.x = element_text(
        size = 7,
        colour = COL_TEXT
      ),
      axis.text.y = element_text(
        size = 7,
        colour = COL_TEXT
      ),
      axis.line = element_line(
        linewidth = 0.45,
        colour = COL_TEXT
      ),
      axis.ticks = element_line(
        linewidth = 0.35,
        colour = COL_TEXT
      ),
      panel.grid.major.y = element_line(
        linewidth = 0.28,
        colour = COL_GRID
      ),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = if (show_legend) "bottom" else "none",
      legend.text = element_text(size = 7, colour = COL_TEXT),
      plot.margin = margin(5, 6, 6, 6)
    )

  p
}

# ------------------------------------------------------------------------------
# 9. BUILD COMPOSITE
# ------------------------------------------------------------------------------
plot_list <- lapply(
  seq_len(nrow(panel_map)),
  function(i) {
    build_panel(
      panel_map$TF[i],
      show_legend = (i == 1L)
    )
  }
)

# ---------------------------------------------------------------------------
# HEADER + 2x2 BODY
# ---------------------------------------------------------------------------
# Put panel label D and the shared legend in the SAME top header row.
# D occupies a dedicated far-left cell and is vertically positioned slightly
# higher than the legend, so it can never appear below the legend.
p_panel_D <- ggplot() +
  annotate(
    "text",
    x = 0.06,
    y = 0.88,
    label = "D",
    hjust = 0,
    vjust = 1,
    family = FONT_FAMILY,
    fontface = "bold",
    size = 4.25,
    colour = COL_TEXT
  ) +
  coord_cartesian(
    xlim = c(0, 1),
    ylim = c(0, 1),
    clip = "off"
  ) +
  theme_void() +
  theme(
    plot.margin = margin(0, 0, 0, 0)
  )

# Top-level patchwork design:
#   row 1: D | collected legend | blank
#   row 2: Nfkb1 | Nfe2l2
#   row 3: Relb  | Nfkb2
#
# Keeping guide_area() on the SAME nesting level as all four plots ensures
# patchwork can place the collected legend exactly in the header.
layout_design <- "
ABBC
DDEE
FFGG
"

Figure3D <- wrap_plots(
  A = p_panel_D,
  B = guide_area(),
  C = plot_spacer(),
  D = plot_list[[1]],
  E = plot_list[[2]],
  F = plot_list[[3]],
  G = plot_list[[4]],
  design = layout_design,
  widths = c(0.13, 1, 1, 0.13),
  heights = c(0.115, 1, 1),
  guides = "collect"
) &
  theme(
    legend.position = "top",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.justification = "center",
    legend.margin = margin(0, 0, 0, 0),
    legend.text = element_text(
      size = 7,
      colour = COL_TEXT
    )
  )

# ------------------------------------------------------------------------------
# 10. EXPORT
# ------------------------------------------------------------------------------
WIDTH_IN  <- FIG_WIDTH_MM / 25.4
HEIGHT_IN <- FIG_HEIGHT_MM / 25.4

PDF_FILE <- file.path(
  OUT_DIR,
  "Figure3D_early_regulon_dynamics_v5.pdf"
)

SVG_FILE <- file.path(
  OUT_DIR,
  "Figure3D_early_regulon_dynamics_v5.svg"
)

TIFF_FILE <- file.path(
  OUT_DIR,
  "Figure3D_early_regulon_dynamics_v5.tiff"
)

PNG_FILE <- file.path(
  OUT_DIR,
  "Figure3D_early_regulon_dynamics_v5.png"
)

ggsave(
  filename = PDF_FILE,
  plot = Figure3D,
  device = function(filename, width, height, ...) {
    grDevices::cairo_pdf(
      filename = filename,
      width = width,
      height = height,
      family = FONT_FAMILY,
      onefile = TRUE
    )
  },
  width = WIDTH_IN,
  height = HEIGHT_IN,
  units = "in",
  bg = "white"
)

ggsave(
  filename = SVG_FILE,
  plot = Figure3D,
  device = svglite::svglite,
  width = WIDTH_IN,
  height = HEIGHT_IN,
  units = "in",
  bg = "white"
)

ggsave(
  filename = TIFF_FILE,
  plot = Figure3D,
  device = ragg::agg_tiff,
  width = WIDTH_IN,
  height = HEIGHT_IN,
  units = "in",
  res = 600,
  compression = "lzw",
  background = "white"
)

ggsave(
  filename = PNG_FILE,
  plot = Figure3D,
  device = ragg::agg_png,
  width = WIDTH_IN,
  height = HEIGHT_IN,
  units = "in",
  res = 300,
  background = "white"
)

# ------------------------------------------------------------------------------
# 11. MD5 / SESSION
# ------------------------------------------------------------------------------
source_md5 <- data.frame(
  file = c(STEP7F_RDS, STEP7F_TIMING_TSV, STEP7F_BINS_TSV),
  md5 = unname(tools::md5sum(c(STEP7F_RDS, STEP7F_TIMING_TSV, STEP7F_BINS_TSV))),
  stringsAsFactors = FALSE
)

write.csv(
  source_md5,
  file.path(OUT_DIR, "source_file_md5.csv"),
  row.names = FALSE
)

writeLines(
  capture.output(sessionInfo()),
  file.path(OUT_DIR, "sessionInfo.txt")
)

# ------------------------------------------------------------------------------
# 12. CONSOLE SUMMARY
# ------------------------------------------------------------------------------
cat("\n============================================================\n")
cat("Figure 3D v5 completed successfully.\n")
cat("Output directory:\n", OUT_DIR, "\n", sep = "")
cat("\nPanel regulons used (data-derived from earliest Figure 3C group):\n")
print(panel_map, row.names = FALSE)
cat("\nDesign:\n")
cat("- 2x2 branch-progress trajectory panels\n- panel label D placed at the true upper-left of the header\n- shared legend centered in the same header row, never above D\n- per-panel qualified-onset subtitles removed\n")
cat("- frozen Step7F fitted curves\n")
cat("- sample-balanced bin means with SD error bars\n")
cat("- dashed line at x = 0 marks the branch origin\n")
cat("- Shared prebranch / Target program / Alternative branch shown together\n")
cat("\nInterpretation boundary:\n")
cat("- Figure 3D shows branch-associated regulon dynamics only.\n")
cat("- It does NOT by itself prove causal regulation.\n")
cat("============================================================\n")
