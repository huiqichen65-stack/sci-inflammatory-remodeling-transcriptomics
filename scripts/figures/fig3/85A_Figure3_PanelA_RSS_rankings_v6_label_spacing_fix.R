# ==============================================================================
# 85A_Figure3_PanelA_RSS_rankings_v6_label_spacing_fix.R
#
# Figure 3A | Three-state RSS ranking panel in a reference-like style
#
# DESIGN UPDATE AFTER VISUAL REVIEW
# ---------------------------------
# This version keeps the SCI/branching logic of Figure 3A:
#   Shared prebranch | Target program | Alternative branch
# but redraws the panel to look much closer to the reference Fig. 5F:
#
#   - each state is a compact RSS ranking plot
#   - all regulons shown as a grey descending point cloud
#   - only a small set of selected regulons is highlighted
#   - labels are arranged in a clean right-side stack
#   - fine connector segments link labels to their true RSS positions
#   - no hollow circles
#   - Target program is visually central and slightly wider
#
# SCIENTIFIC INTERPRETATION BOUNDARY
# ----------------------------------
# - Figure 3A shows RSS specificity only.
# - It does NOT make causality or direct-regulation claims.
# - Nfkb1 / Nfe2l2 / Relb / Nfkb2 are retained in the Target-program panel
#   without bolding or enlarged font because they are pre-specified downstream branchpoint-coincident
#   regulons used in Figure 3C/D, not because they are necessarily the top four
#   RSS regulons.
# - All RSS values and ranks are read directly from frozen Step7E results.
# - Plot labels use the frozen regulon names (displayed as e.g. Nfkb1 (+)) rather than
#   bare TF names, because Figure 3A ranks regulons, not TF expression.
#
# OUTPUT
# ------
# <PROJECT_ROOT>/09_final_figures/GLIA/Figure3A_v6/
#
#   Figure3A_RSS_rankings_v6.pdf
#   Figure3A_RSS_rankings_v6.svg
#   Figure3A_RSS_rankings_v6.tiff
#   Figure3A_RSS_rankings_v6.png
#   00_input_audit.csv
#   01_schema_audit.csv
#   02_Figure3A_plot_data.csv
#   03_Figure3A_label_audit.csv
#   04_Figure3A_key_regulon_audit.csv
#   source_file_md5.csv
#   sessionInfo.txt
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
  "09_final_figures", "GLIA", "Figure3A_v6"
)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

FONT_FAMILY <- "Arial"

KEY_REGULONS <- c("Nfkb1", "Nfe2l2", "Relb", "Nfkb2")

# Shared/Other are comparison panels; Target is the main panel.
N_LABEL_SHARED <- 4L
N_LABEL_TARGET <- 6L
N_LABEL_OTHER  <- 4L

COL_SHARED <- "#7A6AA6"
COL_TARGET <- "#B64B45"
COL_OTHER  <- "#7F8C8D"
COL_GREY   <- "#CFCFCF"
COL_TEXT   <- "#222222"

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

FIG_WIDTH_MM  <- 182
FIG_HEIGHT_MM <- 70

# ------------------------------------------------------------------------------
# 1. PACKAGES
# ------------------------------------------------------------------------------
required_pkgs <- c(
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
  library(ggplot2)
  library(patchwork)
})

# ------------------------------------------------------------------------------
# 2. INPUT
# ------------------------------------------------------------------------------
STEP7E_RDS <- file.path(
  PROJECT_ROOT,
  "04_results",
  "step7E_GSE162610_RSS_CSI_regulon_specificity_modules",
  "RDS",
  "GSE162610_step7E_RSS_CSI_results.rds"
)

input_audit <- data.frame(
  file = "Step7E frozen RDS",
  path = STEP7E_RDS,
  exists = file.exists(STEP7E_RDS),
  size_bytes = if (file.exists(STEP7E_RDS)) file.info(STEP7E_RDS)$size else NA_real_,
  stringsAsFactors = FALSE
)

write.csv(
  input_audit,
  file.path(OUT_DIR, "00_input_audit.csv"),
  row.names = FALSE
)

if (!input_audit$exists) {
  stop("Step7E frozen RDS not found. Inspect 00_input_audit.csv.", call. = FALSE)
}

step7e <- readRDS(STEP7E_RDS)

# ------------------------------------------------------------------------------
# 3. SCHEMA AUDIT
# ------------------------------------------------------------------------------
required_rds_elements <- c("RSS_primary", "RSS_primary_long", "parameters")
required_long_cols <- c(
  "regulon", "TF", "state", "RSS", "RSS_rank_in_state",
  "preferred_state", "preferred_RSS", "RSS_margin"
)

schema_audit <- rbind(
  data.frame(
    source = "Step7E RDS",
    item = required_rds_elements,
    present = required_rds_elements %in% names(step7e),
    stringsAsFactors = FALSE
  ),
  data.frame(
    source = "Step7E RSS_primary_long",
    item = required_long_cols,
    present = if ("RSS_primary_long" %in% names(step7e)) {
      required_long_cols %in% colnames(step7e$RSS_primary_long)
    } else {
      rep(FALSE, length(required_long_cols))
    },
    stringsAsFactors = FALSE
  )
)

write.csv(
  schema_audit,
  file.path(OUT_DIR, "01_schema_audit.csv"),
  row.names = FALSE
)

if (!all(schema_audit$present)) {
  stop("Schema audit failed. Inspect 01_schema_audit.csv.", call. = FALSE)
}

rss_long <- as.data.frame(step7e$RSS_primary_long, stringsAsFactors = FALSE)
rss_long$RSS <- as.numeric(rss_long$RSS)
rss_long$RSS_rank_in_state <- as.integer(rss_long$RSS_rank_in_state)

expected_states <- c("Shared_prebranch", "Arm_2", "Arm_1")
if (!all(expected_states %in% unique(rss_long$state))) {
  stop("Unexpected state set in Step7E RSS_primary_long.", call. = FALSE)
}

rss_long$state_display <- unname(STATE_DISPLAY[as.character(rss_long$state)])
if (anyNA(rss_long$state_display)) {
  stop("Display-state remapping failed.", call. = FALSE)
}

# Exact frozen target-program audit for the four downstream key regulons.
key_audit <- rss_long[
  rss_long$state == "Arm_2" & rss_long$TF %in% KEY_REGULONS,
  c("regulon", "TF", "RSS", "RSS_rank_in_state", "preferred_state", "preferred_RSS", "RSS_margin"),
  drop = FALSE
]
key_audit <- key_audit[match(KEY_REGULONS, key_audit$TF), , drop = FALSE]

write.csv(
  key_audit,
  file.path(OUT_DIR, "04_Figure3A_key_regulon_audit.csv"),
  row.names = FALSE
)

if (nrow(key_audit) != length(KEY_REGULONS) || anyNA(key_audit$TF)) {
  stop("One or more key regulons are absent from frozen Target-program RSS.", call. = FALSE)
}

# ------------------------------------------------------------------------------
# 4. GLOBAL SCALE SETTINGS
# ------------------------------------------------------------------------------
all_n_regulons <- max(rss_long$RSS_rank_in_state, na.rm = TRUE)
rss_max <- max(rss_long$RSS, na.rm = TRUE)

# Extend x-axis to the right to create in-panel label space.
x_data_max  <- all_n_regulons
x_plot_max  <- all_n_regulons + 70
# State-specific label columns.
# Shared-prebranch labels are shifted left so long names such as Creb3l2 (+)
# remain fully inside the first panel and do not overlap the Target-program y-axis.
X_LABEL_SHARED <- all_n_regulons + 27
X_LABEL_TARGET <- all_n_regulons + 42
X_LABEL_OTHER  <- all_n_regulons + 42

y_min <- 0
y_max <- ceiling((rss_max + 0.02) * 20) / 20

# ------------------------------------------------------------------------------
# 5. THEME
# ------------------------------------------------------------------------------
theme_rank <- function(show_y = TRUE) {
  theme_classic(base_family = FONT_FAMILY, base_size = 8) +
    theme(
      plot.title = element_text(
        size = 9, face = "bold", hjust = 0.5,
        colour = COL_TEXT, margin = margin(b = 2)
      ),
      axis.title.x = element_text(size = 8, colour = COL_TEXT, margin = margin(t = 2)),
      axis.title.y = if (show_y) {
        element_text(size = 8, colour = COL_TEXT, margin = margin(r = 4))
      } else {
        element_blank()
      },
      axis.text.x = element_text(size = 7, colour = COL_TEXT),
      axis.text.y = if (show_y) {
        element_text(size = 7, colour = COL_TEXT)
      } else {
        element_blank()
      },
      axis.ticks.x = element_line(linewidth = 0.35, colour = COL_TEXT),
      axis.ticks.y = if (show_y) {
        element_line(linewidth = 0.35, colour = COL_TEXT)
      } else {
        element_blank()
      },
      axis.line = element_line(linewidth = 0.45, colour = COL_TEXT),
      legend.position = "none",
      plot.margin = margin(5, 5, 5, 5)
    )
}

# ------------------------------------------------------------------------------
# 6. LABEL SELECTION + RIGHT-SIDE STACK LAYOUT
# ------------------------------------------------------------------------------
make_label_table <- function(d, state_name) {

  n_top <- switch(
    state_name,
    "Shared prebranch" = N_LABEL_SHARED,
    "Target program"   = N_LABEL_TARGET,
    "Alternative branch"     = N_LABEL_OTHER,
    stop("Unknown state_name")
  )

  lab <- d[d$RSS_rank_in_state <= n_top, , drop = FALSE]

  if (state_name == "Target program") {
    key_rows <- d[d$TF %in% KEY_REGULONS, , drop = FALSE]
    lab <- rbind(lab, key_rows)
    lab <- lab[!duplicated(lab$TF), , drop = FALSE]
  }

  lab <- lab[order(-lab$RSS, lab$RSS_rank_in_state), , drop = FALSE]

  # Right-side stacked label positions, high to low, with larger vertical spacing.
  n <- nrow(lab)
  if (n == 0L) stop("No labels selected for state: ", state_name, call. = FALSE)

  top_pad <- if (state_name == "Target program") 0.030 else 0.040
  span    <- if (state_name == "Target program") 0.24 else 0.19

  y_top <- min(y_max - top_pad, max(d$RSS, na.rm = TRUE) + 0.008)
  y_bot <- max(y_top - span, 0.10)

  if (n == 1L) {
    y_stack <- y_top
  } else {
    y_stack <- seq(y_top, y_bot, length.out = n)
  }

  # small alternating nudge to avoid visual crowding
  if (n >= 2L) {
    alt <- rep(c(0.000, -0.004), length.out = n)
    y_stack <- y_stack + alt
  }

  lab$label_x <- switch(
    state_name,
    "Shared prebranch" = X_LABEL_SHARED,
    "Target program"   = X_LABEL_TARGET,
    "Alternative branch"     = X_LABEL_OTHER,
    stop("Unknown state_name for label x-position")
  )
  lab$label_y <- y_stack
  lab$point_x <- lab$RSS_rank_in_state
  lab$point_y <- lab$RSS
  lab$state_display <- state_name
  lab$is_key <- lab$TF %in% KEY_REGULONS

  # IMPORTANT: display the frozen regulon name, including "(+)".
  # Selection still uses the TF column, so no ranking/filtering logic changes.
  lab$display_label <- as.character(lab$regulon)

  # Add a small visual gap between the TF name and the positive-regulon suffix:
  # e.g. "Nfkb1 (+)" rather than "Nfkb1(+)".
  lab$display_label <- sub(
    "\\(\\+\\)$",
    " (+)",
    lab$display_label
  )

  lab
}

# ------------------------------------------------------------------------------
# 7. SINGLE PANEL BUILDER
# ------------------------------------------------------------------------------
build_state_plot <- function(state_name, show_y = TRUE) {

  d <- rss_long[rss_long$state_display == state_name, , drop = FALSE]
  d <- d[order(d$RSS_rank_in_state), , drop = FALSE]

  state_col <- STATE_COLORS[[state_name]]
  label_d <- make_label_table(d, state_name)

  # top points highlighted
  hi_tf <- unique(label_d$TF)
  high_d <- d[d$TF %in% hi_tf, , drop = FALSE]

  p <- ggplot(d, aes(x = RSS_rank_in_state, y = RSS)) +
    geom_point(
      colour = COL_GREY,
      size = 1.05,
      alpha = 0.95
    ) +
    geom_point(
      data = high_d,
      colour = state_col,
      size = 1.45,
      alpha = 0.98
    ) +
    geom_segment(
      data = label_d,
      aes(
        x = point_x, y = point_y,
        xend = label_x - 3.6, yend = label_y
      ),
      inherit.aes = FALSE,
      colour = state_col,
      linewidth = 0.24,
      linetype = 3,
      alpha = 0.75
    ) +
    geom_text(
      data = label_d,
      aes(
        x = label_x,
        y = label_y,
        label = display_label
      ),
      inherit.aes = FALSE,
      family = FONT_FAMILY,
      fontface = "plain",
      colour = state_col,
      size = 2.25,
      hjust = 0
    ) +
    coord_cartesian(
      xlim = c(0, x_plot_max),
      ylim = c(y_min, y_max),
      clip = "off"
    ) +
    scale_x_continuous(
      breaks = c(0, 50, 100, 150, 200),
      labels = c("0", "50", "100", "150", "200")
    ) +
    scale_y_continuous(
      breaks = pretty(c(y_min, y_max), n = 5)
    ) +
    labs(
      title = state_name,
      x = "Regulon rank",
      y = if (show_y) "Regulon specificity score (RSS)" else NULL
    ) +
    theme_rank(show_y = show_y)

  list(plot = p, labels = label_d, plot_data = d)
}

# ------------------------------------------------------------------------------
# 8. BUILD PLOTS
# ------------------------------------------------------------------------------
shared_res <- build_state_plot("Shared prebranch", show_y = TRUE)
target_res <- build_state_plot("Target program", show_y = FALSE)
other_res  <- build_state_plot("Alternative branch", show_y = FALSE)

# Target slightly wider because it is the main panel.
Figure3A_core <- (
  shared_res$plot |
  target_res$plot |
  other_res$plot
) +
  plot_layout(widths = c(0.95, 1.18, 0.95))

Figure3A <- wrap_elements(full = Figure3A_core) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(
        family = FONT_FAMILY,
        size = 12,
        face = "bold",
        colour = COL_TEXT
      ),
      plot.tag.position = c(0.004, 0.995)
    )
  )

# ------------------------------------------------------------------------------
# 9. EXPORT AUDITS
# ------------------------------------------------------------------------------
plot_data_export <- rss_long[
  ,
  c(
    "regulon", "TF", "state", "state_display",
    "RSS", "RSS_rank_in_state",
    "preferred_state", "preferred_RSS", "RSS_margin"
  )
]

write.csv(
  plot_data_export,
  file.path(OUT_DIR, "02_Figure3A_plot_data.csv"),
  row.names = FALSE
)

label_audit <- rbind(
  cbind(panel = "Shared prebranch", shared_res$labels),
  cbind(panel = "Target program", target_res$labels),
  cbind(panel = "Alternative branch", other_res$labels)
)

write.csv(
  label_audit,
  file.path(OUT_DIR, "03_Figure3A_label_audit.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 10. EXPORT FIGURE
# ------------------------------------------------------------------------------
WIDTH_IN  <- FIG_WIDTH_MM / 25.4
HEIGHT_IN <- FIG_HEIGHT_MM / 25.4

PDF_FILE  <- file.path(OUT_DIR, "Figure3A_RSS_rankings_v6.pdf")
SVG_FILE  <- file.path(OUT_DIR, "Figure3A_RSS_rankings_v6.svg")
TIFF_FILE <- file.path(OUT_DIR, "Figure3A_RSS_rankings_v6.tiff")
PNG_FILE  <- file.path(OUT_DIR, "Figure3A_RSS_rankings_v6.png")

ggsave(
  filename = PDF_FILE,
  plot = Figure3A,
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
  plot = Figure3A,
  device = svglite::svglite,
  width = WIDTH_IN,
  height = HEIGHT_IN,
  units = "in",
  bg = "white"
)

ggsave(
  filename = TIFF_FILE,
  plot = Figure3A,
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
  plot = Figure3A,
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
  file = STEP7E_RDS,
  md5 = unname(tools::md5sum(STEP7E_RDS)),
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
# 12. SUMMARY
# ------------------------------------------------------------------------------
cat("\n============================================================\n")
cat("Figure 3A v6 completed successfully.\n")
cat("Output directory:\n", OUT_DIR, "\n", sep = "")
cat("\nPanel design:\n")
cat("- three-state comparative RSS ranking panel\n")
cat("- Shared prebranch and Alternative branch are comparison panels\n")
cat("- Target program is the central/main panel\n")
cat("- labels are stacked on the right to mimic the reference style\n")
cat("- no hollow circles; no large emphasis symbols\n")
cat("- Nfkb1/Nfe2l2/Relb/Nfkb2 labels use the same plain font and size as other labels\n")
cat("- regulon labels are displayed with a space before (+), e.g. Nfkb1 (+)\n")
cat("- Shared-prebranch label column is shifted left to prevent Creb3l2 (+) overlap\n")

cat("\nTarget-program key-regulon frozen audit:\n")
print(
  key_audit[, c("TF", "RSS", "RSS_rank_in_state")],
  row.names = FALSE
)

cat("\nInterpretation boundary:\n")
cat("- Figure 3A = RSS specificity landscape only.\n")
cat("- Nfkb1/Nfe2l2/Relb/Nfkb2 are shown at their true frozen ranks.\n")
cat("- No direct-regulation claim is made in this panel.\n")
cat("============================================================\n")
