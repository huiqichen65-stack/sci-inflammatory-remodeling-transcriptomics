# ==============================================================================
# 83B_Figure2_PanelB_shared_transition_gate_v3_mastertag_safe.R
# v3 visual-only fix: remove the internal B tag so the master can overlay
# a uniform panel letter without masking/clipping the y-axis title.
# Figure 2B | Shared activated / DAM-like transition gate
#
# Biological purpose
# ------------------
# Show that Activated and Early-DAM transcriptional programs are not independent;
# rather, a subset of resident microglia concurrently shows high Activated AUC
# and high Early-DAM AUC, defining a shared activated/DAM-like transitional state.
#
# IMPORTANT
# ---------
# - Read frozen Step3C outputs only.
# - Do NOT rerun AUCell, clustering, UMAP, thresholding, or reclassification.
# - Highlighted cells are exactly the frozen 3 dpi high-confidence transition
#   cells from Step3C:
#       shared_transition_3dpi_highconf == TRUE AND timepoint == "3 dpi"
# - Threshold lines are exactly the frozen Step3C thresholds:
#       Activated_high_q90_cluster3
#       EarlyDAM_high_q90_cluster3
#
# Output
# ------
# <PROJECT_ROOT>/09_final_figures/GLIA/Figure2_panelB_v2/
#   Figure2B_shared_transition_gate_v2.pdf
#   Figure2B_shared_transition_gate_v2.svg
#   Figure2B_shared_transition_gate_v2.tiff
#   Figure2B_shared_transition_gate_v2.png
#   00_input_audit.csv
#   01_threshold_audit.csv
#   02_transition_count_audit.csv
#   sessionInfo.txt
# ==============================================================================

rm(list = ls())
gc()
options(stringsAsFactors = FALSE)
set.seed(20260813)

# ------------------------------------------------------------------------------
# 0. USER SETTINGS
# ------------------------------------------------------------------------------
project_root <- Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
out_dir <- file.path(
  project_root,
  "09_final_figures", "GLIA", "Figure2_panelB_v3_mastertag_safe"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Standalone Panel B size
panel_width_mm  <- 76
panel_height_mm <- 60

# ------------------------------------------------------------------------------
# 1. PACKAGES
# ------------------------------------------------------------------------------
required_pkgs <- c(
  "Seurat",
  "ggplot2",
  "dplyr",
  "scales",
  "svglite",
  "ragg"
)

missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_pkgs) > 0) {
  stop(
    "Missing required package(s): ",
    paste(missing_pkgs, collapse = ", "),
    "\nInstall them first, then rerun."
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(scales)
})

# ------------------------------------------------------------------------------
# 2. FROZEN INPUTS
# ------------------------------------------------------------------------------
input_files <- c(
  step3C_rds = file.path(
    project_root,
    "04_results/step3C_GSE162610_activated_DAM_transition_mapping",
    "GSE162610_step3C_transition_mapped.rds"
  ),
  step3C_thresholds = file.path(
    project_root,
    "04_results/step3C_GSE162610_activated_DAM_transition_mapping",
    "tables", "02_main_transition_thresholds.csv"
  )
)

input_audit <- data.frame(
  key = names(input_files),
  path = unname(input_files),
  exists = file.exists(input_files),
  size_MB = ifelse(
    file.exists(input_files),
    round(file.info(input_files)$size / 1024^2, 3),
    NA_real_
  ),
  stringsAsFactors = FALSE
)

write.csv(
  input_audit,
  file.path(out_dir, "00_input_audit.csv"),
  row.names = FALSE
)

if (!all(input_audit$exists)) {
  print(input_audit)
  stop("Required frozen input missing. See 00_input_audit.csv.")
}

message("Reading frozen Step3C inputs...")
obj3C <- readRDS(input_files[["step3C_rds"]])
thr   <- read.csv(input_files[["step3C_thresholds"]], check.names = FALSE)

if (!inherits(obj3C, "Seurat")) {
  stop("Step3C frozen input must be a Seurat object.")
}

# ------------------------------------------------------------------------------
# 3. SCHEMA CHECKS
# ------------------------------------------------------------------------------
required_meta <- c(
  "TransAUC_Activated_injury_core",
  "TransAUC_Early_DAM_transition",
  "shared_transition_3dpi_highconf",
  "timepoint"
)

missing_meta <- setdiff(required_meta, colnames(obj3C[[]]))
if (length(missing_meta) > 0) {
  stop("Step3C metadata missing: ", paste(missing_meta, collapse = ", "))
}

required_thr_cols <- c("threshold", "value")
missing_thr <- setdiff(required_thr_cols, colnames(thr))
if (length(missing_thr) > 0) {
  stop("Threshold table missing required columns: ", paste(missing_thr, collapse = ", "))
}

act_thr <- thr$value[thr$threshold == "Activated_high_q90_cluster3"]
early_thr <- thr$value[thr$threshold == "EarlyDAM_high_q90_cluster3"]

if (length(act_thr) != 1 || length(early_thr) != 1) {
  stop(
    "Could not uniquely recover frozen Step3C thresholds: ",
    "Activated_high_q90_cluster3 and/or EarlyDAM_high_q90_cluster3."
  )
}

threshold_audit <- data.frame(
  threshold = c("Activated_high_q90_cluster3", "EarlyDAM_high_q90_cluster3"),
  value = c(act_thr, early_thr),
  stringsAsFactors = FALSE
)

write.csv(
  threshold_audit,
  file.path(out_dir, "01_threshold_audit.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 4. MASTER PLOT DATA
# ------------------------------------------------------------------------------
meta <- obj3C[[]]

plot_df <- data.frame(
  cell_id = rownames(meta),
  Activated_AUCell = as.numeric(meta$TransAUC_Activated_injury_core),
  EarlyDAM_AUCell = as.numeric(meta$TransAUC_Early_DAM_transition),
  timepoint = as.character(meta$timepoint),
  highlight = as.logical(meta$shared_transition_3dpi_highconf) &
    as.character(meta$timepoint) == "3 dpi",
  stringsAsFactors = FALSE
)

if (
  any(!is.finite(plot_df$Activated_AUCell)) ||
  any(!is.finite(plot_df$EarlyDAM_AUCell))
) {
  stop("Non-finite AUCell values detected.")
}

transition_audit <- data.frame(
  group = c(
    "All resident microglia",
    "3 dpi high-confidence transition"
  ),
  n_cells = c(
    nrow(plot_df),
    sum(plot_df$highlight, na.rm = TRUE)
  ),
  stringsAsFactors = FALSE
)

write.csv(
  transition_audit,
  file.path(out_dir, "02_transition_count_audit.csv"),
  row.names = FALSE
)

# Deterministic plotting order so highlighted cells are not buried.
set.seed(20260813)
bg_df <- plot_df %>%
  filter(!highlight) %>%
  slice_sample(prop = 1)

fg_df <- plot_df %>%
  filter(highlight) %>%
  slice_sample(prop = 1)

# ------------------------------------------------------------------------------
# 5. STYLE
# ------------------------------------------------------------------------------
base_family <- "Arial"
# V2 restrained publication palette.
# Neutral background cells + muted brick-red transition signal.
transition_col <- "#B85C57"
bg_col <- "#B9BDC3"
quad_fill <- "#F3E7E5"

theme_pub <- theme_classic(base_family = base_family, base_size = 7.2) +
  theme(
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),

    axis.title = element_text(
      family = base_family,
      size = 7.0,
      colour = "black"
    ),
    axis.text = element_text(
      family = base_family,
      size = 6.2,
      colour = "black"
    ),

    plot.tag = element_text(
      family = base_family,
      face = "bold",
      size = 11.5,
      colour = "black",
      hjust = 0,
      vjust = 1
    ),
    plot.tag.position = c(0.012, 0.988),

    plot.margin = margin(4.5, 3.2, 2.8, 4.2, unit = "mm")
  )

geom_cells <- function(...) {
  if (requireNamespace("ggrastr", quietly = TRUE)) {
    ggrastr::geom_point_rast(
      ...,
      raster.dpi = 500,
      stroke = 0
    )
  } else {
    geom_point(...)
  }
}

# ------------------------------------------------------------------------------
# 6. AXIS / ANNOTATION GEOMETRY
# ------------------------------------------------------------------------------
xr <- range(plot_df$Activated_AUCell, finite = TRUE)
yr <- range(plot_df$EarlyDAM_AUCell, finite = TRUE)
dx <- diff(xr)
dy <- diff(yr)

x_limits <- c(
  xr[1] - 0.04 * dx,
  xr[2] + 0.05 * dx
)
y_limits <- c(
  yr[1] - 0.04 * dy,
  yr[2] + 0.05 * dy
)

# Upper-right region shading to quietly indicate the dual-high gate.
rect_xmin <- act_thr
rect_xmax <- x_limits[2]
rect_ymin <- early_thr
rect_ymax <- y_limits[2]

# Label location within the shaded dual-high region.
label_x <- rect_xmax - 0.045 * dx
label_y <- rect_ymax - 0.060 * dy

# ------------------------------------------------------------------------------
# 7. PLOT
# ------------------------------------------------------------------------------
pB <- ggplot() +
  annotate(
    "rect",
    xmin = rect_xmin,
    xmax = rect_xmax,
    ymin = rect_ymin,
    ymax = rect_ymax,
    fill = quad_fill,
    alpha = 0.50
  ) +

  geom_cells(
    data = bg_df,
    aes(x = Activated_AUCell, y = EarlyDAM_AUCell),
    colour = bg_col,
    size = 0.25,
    alpha = 0.20
  ) +

  geom_cells(
    data = fg_df,
    aes(x = Activated_AUCell, y = EarlyDAM_AUCell),
    colour = transition_col,
    size = 0.40,
    alpha = 0.78
  ) +

  geom_vline(
    xintercept = act_thr,
    linetype = "22",
    linewidth = 0.31,
    colour = "#707070"
  ) +

  geom_hline(
    yintercept = early_thr,
    linetype = "22",
    linewidth = 0.31,
    colour = "#707070"
  ) +

  annotate(
    "text",
    x = label_x,
    y = label_y,
    label = "Shared activated /\nDAM-like transition",
    family = base_family,
    fontface = "bold",
    size = 1.85,
    hjust = 1,
    vjust = 1,
    lineheight = 0.92,
    colour = transition_col
  ) +

  scale_x_continuous(
    limits = x_limits,
    expand = expansion(mult = 0),
    labels = label_number(accuracy = 0.1)
  ) +

  scale_y_continuous(
    limits = y_limits,
    expand = expansion(mult = 0),
    labels = label_number(accuracy = 0.1)
  ) +

  labs(
    x = "Activated program (AUCell AUC)",
    y = "Early-DAM program (AUCell AUC)",
    tag = NULL
  ) +
  theme_pub +
  theme(
    axis.line = element_line(linewidth = 0.40, colour = "#333333"),
    plot.tag = element_blank()
  )

# ------------------------------------------------------------------------------
# 8. SAVE
# ------------------------------------------------------------------------------
stem <- "Figure2B_shared_transition_gate_v3_mastertag_safe"

ggsave(
  file.path(out_dir, paste0(stem, ".pdf")),
  pB,
  device = cairo_pdf,
  width = panel_width_mm,
  height = panel_height_mm,
  units = "mm",
  bg = "white"
)

ggsave(
  file.path(out_dir, paste0(stem, ".svg")),
  pB,
  device = svglite::svglite,
  width = panel_width_mm,
  height = panel_height_mm,
  units = "mm",
  bg = "white"
)

ggsave(
  file.path(out_dir, paste0(stem, ".tiff")),
  pB,
  device = ragg::agg_tiff,
  width = panel_width_mm,
  height = panel_height_mm,
  units = "mm",
  dpi = 600,
  compression = "lzw",
  background = "white"
)

ggsave(
  file.path(out_dir, paste0(stem, ".png")),
  pB,
  device = ragg::agg_png,
  width = panel_width_mm,
  height = panel_height_mm,
  units = "mm",
  dpi = 300,
  background = "white"
)

capture.output(
  sessionInfo(),
  file = file.path(out_dir, "sessionInfo.txt")
)

message("")
message("============================================================")
message("Figure 2B v3 master-tag-safe COMPLETE")
message("Biological message:")
message("Concurrent elevation of Activated and Early-DAM programs defines")
message("a shared activated/DAM-like transitional state.")
message("")
message("Primary review file:")
message(file.path(out_dir, paste0(stem, '.pdf')))
message("============================================================")
