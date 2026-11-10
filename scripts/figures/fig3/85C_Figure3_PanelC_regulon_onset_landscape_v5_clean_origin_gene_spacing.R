# ==============================================================================
# 85C_Figure3_PanelC_regulon_onset_landscape_v5_clean_origin_gene_spacing.R
#
# Figure 3C | Temporal prioritization landscape of Target-associated regulons
#
# DESIGN GOAL
# -----------
# Replace the visually plain single lollipop plot with a more publication-grade
# composite panel while preserving exactly the same frozen Step7F evidence.
#
# LEFT  : onset of branch-associated divergence (main timeline)
# RIGHT : Target-program regulon specificity (preferred RSS), aligned by TF
#
# The earliest frozen onset group is highlighted as a light red band and with
# red TF labels. This makes the selection of Nfkb1/Nfe2l2/Relb/Nfkb2 visually
# explicit without manually forcing four genes into the panel.
#
# IMPORTANT
# ---------
# - No timing model is rerun.
# - No Y_29 label is shown in the main figure.
# - "Branch origin" is the display language for onset = 0.
# - Earliest group is defined from the minimum frozen onset_progress_qualified.
# - RSS is read directly from frozen Step7F preferred_RSS.
#
# OUTPUT
# ------
# <PROJECT_ROOT>/09_final_figures/GLIA/Figure3C_v5/
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
  "09_final_figures", "GLIA", "Figure3C_v5"
)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

FONT_FAMILY <- "Arial"

FIG_WIDTH_MM  <- 150
FIG_HEIGHT_MM <- 88

COL_TEXT       <- "#222222"
COL_MUTED      <- "#7F8C8D"
COL_LINE       <- "#B9BDC8"
COL_GRID       <- "#E9E9E9"
COL_TARGET     <- "#B64B45"
COL_TARGET_BG  <- "#F7E8E5"
COL_RSS_FILL   <- "#B8A6D8"
COL_RSS_DARK   <- "#6F5B97"
COL_SEPARATOR  <- "#D5D5D5"

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
# 2. INPUT
# ------------------------------------------------------------------------------
STEP7F_TIMING_TSV <- file.path(
  PROJECT_ROOT,
  "04_results",
  "step7F_GSE162610_regulon_fate_timing",
  "tables",
  "04_key_candidate_fate_timing.tsv"
)

input_audit <- data.frame(
  file = "Step7F key candidate timing table",
  path = STEP7F_TIMING_TSV,
  exists = file.exists(STEP7F_TIMING_TSV),
  size_bytes = if (file.exists(STEP7F_TIMING_TSV)) {
    file.info(STEP7F_TIMING_TSV)$size
  } else {
    NA_real_
  },
  stringsAsFactors = FALSE
)

write.csv(
  input_audit,
  file.path(OUT_DIR, "00_input_audit.csv"),
  row.names = FALSE
)

if (!input_audit$exists) {
  stop("Figure 3C input missing. Inspect 00_input_audit.csv.", call. = FALSE)
}

timing_df <- fread(STEP7F_TIMING_TSV, data.table = FALSE)

# ------------------------------------------------------------------------------
# 3. SCHEMA AUDIT
# ------------------------------------------------------------------------------
required_cols <- c(
  "regulon",
  "TF",
  "model_success",
  "onset_progress_qualified",
  "onset_category",
  "priority_group",
  "direction_matches_RSS",
  "seven_dpi_direction_consistent",
  "RSS_preferred_state",
  "preferred_RSS",
  "RSS_margin",
  "CSI_module",
  "terminal_effect_z",
  "timing_classification"
)

schema_audit <- data.frame(
  source = "04_key_candidate_fate_timing.tsv",
  item = required_cols,
  present = required_cols %in% colnames(timing_df),
  stringsAsFactors = FALSE
)

write.csv(
  schema_audit,
  file.path(OUT_DIR, "01_schema_audit.csv"),
  row.names = FALSE
)

if (!all(schema_audit$present)) {
  stop("Figure 3C schema audit failed. Inspect 01_schema_audit.csv.", call. = FALSE)
}

# ------------------------------------------------------------------------------
# 4. FROZEN TARGET-ASSOCIATED CANDIDATE SET
# ------------------------------------------------------------------------------
plot_df <- timing_df[
  timing_df$priority_group == "Arm2_priority" &
    timing_df$model_success %in% TRUE &
    !is.na(timing_df$onset_progress_qualified) &
    timing_df$direction_matches_RSS %in% TRUE &
    timing_df$seven_dpi_direction_consistent %in% TRUE &
    timing_df$RSS_preferred_state == "Arm_2",
  ,
  drop = FALSE
]

if (nrow(plot_df) == 0L) {
  stop("No frozen Target-associated timing candidates remain.", call. = FALSE)
}

plot_df$onset_progress_qualified <- as.numeric(plot_df$onset_progress_qualified)
plot_df$preferred_RSS <- as.numeric(plot_df$preferred_RSS)
plot_df$RSS_margin <- as.numeric(plot_df$RSS_margin)
plot_df$terminal_effect_z <- as.numeric(plot_df$terminal_effect_z)

min_onset <- min(plot_df$onset_progress_qualified, na.rm = TRUE)
plot_df$is_earliest_group <- abs(
  plot_df$onset_progress_qualified - min_onset
) < 1e-12

# Order: onset first, then RSS descending.
plot_df <- plot_df[
  order(
    plot_df$onset_progress_qualified,
    -plot_df$preferred_RSS,
    plot_df$TF
  ),
  ,
  drop = FALSE
]

plot_df$display_label <- as.character(plot_df$TF)
plot_df$display_factor <- factor(
  plot_df$display_label,
  levels = rev(plot_df$display_label)
)

plot_df$onset_text <- sprintf("%.2f", plot_df$onset_progress_qualified)

write.csv(
  plot_df[
    ,
    c(
      "TF",
      "regulon",
      "onset_progress_qualified",
      "onset_category",
      "preferred_RSS",
      "RSS_margin",
      "terminal_effect_z",
      "CSI_module",
      "timing_classification",
      "is_earliest_group"
    )
  ],
  file.path(OUT_DIR, "02_Figure3C_plot_data.csv"),
  row.names = FALSE
)

write.csv(
  plot_df[
    plot_df$is_earliest_group,
    c(
      "TF",
      "regulon",
      "onset_progress_qualified",
      "preferred_RSS",
      "terminal_effect_z",
      "CSI_module",
      "timing_classification"
    ),
    drop = FALSE
  ],
  file.path(OUT_DIR, "03_Figure3C_earliest_group_audit.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 5. COMMON SETTINGS
# ------------------------------------------------------------------------------
n_tf <- nrow(plot_df)

# Earliest rows are at the top after factor reversal.
early_tf <- plot_df$TF[plot_df$is_earliest_group]
early_levels_numeric <- which(
  levels(plot_df$display_factor) %in% early_tf
)

# Continuous y-coordinate band around earliest rows
early_ymin <- min(early_levels_numeric) - 0.48
early_ymax <- max(early_levels_numeric) + 0.48

x_onset_max <- max(plot_df$onset_progress_qualified, na.rm = TRUE)
x_onset_upper <- max(0.60, ceiling((x_onset_max + 0.04) * 20) / 20)

rss_min <- floor((min(plot_df$preferred_RSS) - 0.02) * 20) / 20
rss_max <- ceiling((max(plot_df$preferred_RSS) + 0.015) * 20) / 20

# ------------------------------------------------------------------------------
# 6. MAIN ONSET TIMELINE
# ------------------------------------------------------------------------------
p_onset <- ggplot(
  plot_df,
  aes(
    x = onset_progress_qualified,
    y = display_factor
  )
) +
  # very light highlight behind earliest rows
  annotate(
    "rect",
    xmin = -Inf,
    xmax = Inf,
    ymin = early_ymin,
    ymax = early_ymax,
    fill = COL_TARGET_BG,
    colour = NA
  ) +
  # horizontal guide/lollipop stems
  geom_segment(
    aes(
      x = 0,
      xend = onset_progress_qualified,
      y = display_factor,
      yend = display_factor
    ),
    linewidth = 0.60,
    colour = COL_LINE,
    lineend = "round"
  ) +
  # branch-origin line
  geom_vline(
    xintercept = 0,
    linewidth = 0.55,
    colour = COL_TEXT
  ) +
  # onset dots
  geom_point(
    aes(
      fill = is_earliest_group
    ),
    shape = 21,
    size = 2.9,
    stroke = 0.45,
    colour = "white"
  ) +
  # onset values
  geom_text(
    aes(
      label = onset_text,
      colour = is_earliest_group
    ),
    nudge_x = 0.018,
    hjust = 0,
    family = FONT_FAMILY,
    size = 2.45
  ) +
  # label above origin
  annotate(
    "text",
    x = 0.060,
    y = n_tf + 0.58,
    label = "Branch origin",
    hjust = 0,
    vjust = 0.5,
    family = FONT_FAMILY,
    fontface = "bold",
    size = 2.65,
    colour = COL_TARGET
  ) +
  scale_fill_manual(
    values = c(
      "FALSE" = COL_MUTED,
      "TRUE" = COL_TARGET
    ),
    guide = "none"
  ) +
  scale_colour_manual(
    values = c(
      "FALSE" = COL_MUTED,
      "TRUE" = COL_TARGET
    ),
    guide = "none"
  ) +
  scale_x_continuous(
    limits = c(0, x_onset_upper),
    breaks = seq(0, x_onset_upper, by = 0.1),
    expand = expansion(mult = c(0, 0.06))
  ) +
  scale_y_discrete(
    expand = expansion(add = c(0.50, 1.35))
  ) +
  coord_cartesian(clip = "off") +
  labs(
    tag = "C",
    x = "Branch progress at onset of divergence",
    y = NULL
  ) +
  theme_classic(
    base_family = FONT_FAMILY,
    base_size = 8
  ) +
  theme(
    plot.tag = element_text(
      family = FONT_FAMILY,
      size = 12,
      face = "bold",
      colour = COL_TEXT
    ),
    plot.tag.position = c(0.004, 0.995),
    axis.title.x = element_text(
      size = 8,
      colour = COL_TEXT,
      margin = margin(t = 4)
    ),
    axis.text.x = element_text(
      size = 7,
      colour = COL_TEXT
    ),
    axis.text.y = element_text(
      size = 7.6,
      colour = ifelse(
        levels(plot_df$display_factor) %in% early_tf,
        COL_TARGET,
        COL_TEXT
      ),
      margin = margin(r = 9)
    ),
    axis.line.x = element_line(
      linewidth = 0.45,
      colour = COL_TEXT
    ),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_line(
      linewidth = 0.35,
      colour = COL_TEXT
    ),
    panel.grid.major.x = element_line(
      linewidth = 0.28,
      colour = COL_GRID
    ),
    panel.grid.minor = element_blank(),
    plot.margin = margin(14, 8, 6, 10)
  )

# ------------------------------------------------------------------------------
# 7. ALIGNED RSS MINI-BAR PANEL
# ------------------------------------------------------------------------------
p_rss <- ggplot(
  plot_df,
  aes(
    x = preferred_RSS,
    y = display_factor
  )
) +
  annotate(
    "rect",
    xmin = -Inf,
    xmax = Inf,
    ymin = early_ymin,
    ymax = early_ymax,
    fill = COL_TARGET_BG,
    colour = NA
  ) +
  geom_segment(
    aes(
      x = rss_min,
      xend = preferred_RSS,
      y = display_factor,
      yend = display_factor,
      colour = is_earliest_group
    ),
    linewidth = 2.8,
    lineend = "round"
  ) +
  geom_point(
    aes(
      x = preferred_RSS,
      fill = is_earliest_group
    ),
    shape = 21,
    size = 2.7,
    stroke = 0.45,
    colour = "white"
  ) +
  scale_colour_manual(
    values = c(
      "FALSE" = COL_RSS_FILL,
      "TRUE" = COL_TARGET
    ),
    guide = "none"
  ) +
  scale_fill_manual(
    values = c(
      "FALSE" = COL_RSS_DARK,
      "TRUE" = COL_TARGET
    ),
    guide = "none"
  ) +
  scale_x_continuous(
    limits = c(rss_min, rss_max),
    breaks = pretty(c(rss_min, rss_max), n = 3),
    expand = expansion(mult = c(0, 0.03))
  ) +
  scale_y_discrete(
    expand = expansion(add = c(0.50, 1.35))
  ) +
  labs(
    x = "Target-program RSS",
    y = NULL
  ) +
  theme_classic(
    base_family = FONT_FAMILY,
    base_size = 8
  ) +
  theme(
    axis.title.x = element_text(
      size = 8,
      colour = COL_TEXT,
      margin = margin(t = 4)
    ),
    axis.text.x = element_text(
      size = 7,
      colour = COL_TEXT
    ),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    axis.line.x = element_line(
      linewidth = 0.45,
      colour = COL_TEXT
    ),
    axis.ticks.x = element_line(
      linewidth = 0.35,
      colour = COL_TEXT
    ),
    panel.grid = element_blank(),
    plot.margin = margin(14, 6, 6, 2)
  )

# ------------------------------------------------------------------------------
# 8. COMPOSITE
# ------------------------------------------------------------------------------
Figure3C <- p_onset + p_rss +
  plot_layout(
    widths = c(2.25, 0.95)
  )

# ------------------------------------------------------------------------------
# 9. EXPORT
# ------------------------------------------------------------------------------
WIDTH_IN  <- FIG_WIDTH_MM / 25.4
HEIGHT_IN <- FIG_HEIGHT_MM / 25.4

PDF_FILE <- file.path(
  OUT_DIR,
  "Figure3C_regulon_onset_landscape_v5.pdf"
)

SVG_FILE <- file.path(
  OUT_DIR,
  "Figure3C_regulon_onset_landscape_v5.svg"
)

TIFF_FILE <- file.path(
  OUT_DIR,
  "Figure3C_regulon_onset_landscape_v5.tiff"
)

PNG_FILE <- file.path(
  OUT_DIR,
  "Figure3C_regulon_onset_landscape_v5.png"
)

ggsave(
  filename = PDF_FILE,
  plot = Figure3C,
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
  plot = Figure3C,
  device = svglite::svglite,
  width = WIDTH_IN,
  height = HEIGHT_IN,
  units = "in",
  bg = "white"
)

ggsave(
  filename = TIFF_FILE,
  plot = Figure3C,
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
  plot = Figure3C,
  device = ragg::agg_png,
  width = WIDTH_IN,
  height = HEIGHT_IN,
  units = "in",
  res = 300,
  background = "white"
)

writeLines(
  capture.output(sessionInfo()),
  file.path(OUT_DIR, "sessionInfo.txt")
)

cat("\n============================================================\n")
cat("Figure 3C v5 completed.\n")
cat("Output directory:\n", OUT_DIR, "\n", sep = "")
cat("\nDesign:\n")
cat("- left: branch-divergence onset timeline\n")
cat("- right: aligned Target-program RSS mini-bars\n")
cat("- earliest data-derived group highlighted by a light red band\n")
cat("- no Y_29 terminology in the main figure\n")
cat("============================================================\n")
