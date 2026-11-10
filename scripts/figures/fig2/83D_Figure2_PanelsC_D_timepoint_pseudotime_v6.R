# ==============================================================================
# 83D_Figure2_PanelsC_D_timepoint_pseudotime_v6.R
# Figure 2C-D | Split injury-time maps (C) and corrected pseudotime (D)
#
# Biological purpose
# ------------------
# Revised Figure 2 layout:
#   (1) Panel C: show 1 dpi, 3 dpi, and 7 dpi separately on the SAME frozen trajectory;
#   (2) Panel D: show corrected pseudotime as a standalone panel.
# Other-timepoint cells are retained as light-gray background in Panel C so that
# all three subplots share exactly the same UMAP frame and principal graph.
#
# IMPORTANT
# ---------
# - Read the frozen corrected Step5BR2 object only.
# - Do NOT rerun preprocess_cds(), reduce_dimension(), cluster_cells(),
#   learn_graph(), or order_cells().
# - Use the existing corrected UMAP / principal graph from Step5BR2.
# - Do NOT display Root, Y_29, leaf labels, branch labels, graph-node IDs, or
#   any audit-style node annotations.
#
# Output
# ------
# <PROJECT_ROOT>/09_final_figures/GLIA/Figure2_panelsC_D_v6/
#   Figure2C_timepoint_split_v6.{pdf,svg,tiff,png}
#   Figure2D_pseudotime_v6.{pdf,svg,tiff,png}
#   Figure2CD_preview_v6.{pdf,png}
#   00_input_audit.csv
#   01_schema_audit.csv
#   02_timepoint_counts.csv
#   03_display_palette.csv
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
  "09_final_figures", "GLIA", "Figure2_panelsC_D_v6"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

panelC_width_mm  <- 138
panelD_width_mm  <- 62
panel_height_mm  <- 54
preview_width_mm <- 180

# ------------------------------------------------------------------------------
# 1. PACKAGES
# ------------------------------------------------------------------------------
required_pkgs <- c(
  "monocle3",
  "ggplot2",
  "patchwork",
  "SingleCellExperiment",
  "SummarizedExperiment",
  "svglite",
  "ragg",
  "scales"
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
  library(monocle3)
  library(ggplot2)
  library(patchwork)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(scales)
})

# ------------------------------------------------------------------------------
# 2. INPUT
# ------------------------------------------------------------------------------
input_files <- c(
  br2_cds = file.path(
    project_root,
    "04_results/step5BR2_GSE162610_single_root_branch_audit",
    "GSE162610_step5BR2_single_root_cds.rds"
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

message("Reading frozen Step5BR2 CDS...")
cds <- readRDS(input_files[["br2_cds"]])

# ------------------------------------------------------------------------------
# 3. BASIC SCHEMA CHECKS
# ------------------------------------------------------------------------------
if (!inherits(cds, "cell_data_set")) {
  stop("Input is not a monocle3 cell_data_set.")
}

required_meta <- c("timepoint", "step5BR2_pseudotime")
missing_meta <- setdiff(required_meta, colnames(SummarizedExperiment::colData(cds)))
if (length(missing_meta) > 0) {
  stop("Frozen Step5BR2 CDS metadata missing: ", paste(missing_meta, collapse = ", "))
}

if (!("UMAP" %in% names(reducedDims(cds)))) {
  stop("Frozen Step5BR2 CDS does not contain reducedDims(cds)$UMAP.")
}

if (is.null(principal_graph(cds)$UMAP)) {
  stop("Frozen Step5BR2 CDS does not contain a UMAP principal graph.")
}

schema_audit <- data.frame(
  item = c(
    "n_cells",
    "has_timepoint",
    "has_step5BR2_pseudotime",
    "has_UMAP_reduction",
    "has_UMAP_principal_graph"
  ),
  value = c(
    ncol(cds),
    "TRUE",
    "TRUE",
    "TRUE",
    "TRUE"
  ),
  stringsAsFactors = FALSE
)

write.csv(
  schema_audit,
  file.path(out_dir, "01_schema_audit.csv"),
  row.names = FALSE
)

timepoint_audit <- as.data.frame(
  table(
    factor(
      as.character(SummarizedExperiment::colData(cds)$timepoint),
      levels = c("1 dpi", "3 dpi", "7 dpi")
    )
  ),
  stringsAsFactors = FALSE
)
colnames(timepoint_audit) <- c("timepoint", "n_cells")

if (sum(timepoint_audit$n_cells) != ncol(cds)) {
  stop(
    "Figure2D timepoint audit failed: not all cells map to 1/3/7 dpi."
  )
}

write.csv(
  timepoint_audit,
  file.path(out_dir, "02_timepoint_counts.csv"),
  row.names = FALSE
)

palette_audit <- data.frame(
  role = c(
    "1 dpi",
    "3 dpi",
    "7 dpi",
    "Pseudotime low",
    "Pseudotime high"
  ),
  color = c(
    "#8C1D9C",
    "#D84A73",
    "#F6B73C",
    "#0D0887",
    "#F0F921"
  ),
  stringsAsFactors = FALSE
)

write.csv(
  palette_audit,
  file.path(out_dir, "03_display_palette.csv"),
  row.names = FALSE
)

# Standardized factor for timepoint plotting.
SummarizedExperiment::colData(cds)$Figure2_time <- factor(
  as.character(SummarizedExperiment::colData(cds)$timepoint),
  levels = c("1 dpi", "3 dpi", "7 dpi")
)

# ------------------------------------------------------------------------------
# 4. VISUAL STYLE
# ------------------------------------------------------------------------------
base_family <- "Arial"
graph_col <- "#554A63"   # harmonize with Figure 2C final lavender scheme

# Discrete injury-time colors are used for the three Panel C subplots.
# Panel D uses the same frozen continuous pseudotime palette as the previous right half.
# No biological quantities are recomputed.
TIME_COLORS <- c(
  "1 dpi" = "#8C1D9C",
  "3 dpi" = "#D84A73",
  "7 dpi" = "#F6B73C"
)

PSEUDOTIME_COLS <- c(
  "#0D0887",
  "#5B02A3",
  "#9A179B",
  "#CC4778",
  "#ED7953",
  "#FCA636",
  "#F0F921"
)

theme_pub <- theme_classic(base_family = base_family, base_size = 7.2) +
  theme(
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),

    axis.title = element_text(
      family = base_family,
      size = 6.8,
      colour = "black"
    ),
    axis.text = element_text(
      family = base_family,
      size = 5.9,
      colour = "black"
    ),

    legend.title = element_text(
      family = base_family,
      size = 5.4,
      colour = "black"
    ),
    legend.text = element_text(
      family = base_family,
      size = 5.2,
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

    plot.margin = margin(2.5, 2.5, 2.5, 2.5, unit = "mm")
  )

# ------------------------------------------------------------------------------
# 5. MONOCLE3-COMPATIBLE PLOT WRAPPER
# ------------------------------------------------------------------------------
make_traj_plot <- function(
  cds_obj,
  color_by,
  color_scale = NULL,
  show_legend = TRUE,
  cell_size = 0.31
) {
  plot_args <- list(
    cds = cds_obj,
    reduction_method = "UMAP",
    color_cells_by = color_by,
    show_trajectory_graph = TRUE,
    label_cell_groups = FALSE,
    label_leaves = FALSE,
    label_branch_points = FALSE,
    label_roots = FALSE,
    graph_label_size = 2.4,
    cell_size = cell_size,
    trajectory_graph_color = graph_col,
    trajectory_graph_segment_size = 0.74
  )

  supported_args <- names(formals(monocle3::plot_cells))
  keep_names <- names(plot_args) %in% supported_args

  dropped_args <- names(plot_args)[!keep_names]
  if (length(dropped_args) > 0) {
    message(
      "plot_cells(): skipping unsupported plotting argument(s): ",
      paste(dropped_args, collapse = ", ")
    )
  }

  p <- do.call(monocle3::plot_cells, plot_args[keep_names])

  if (!is.null(color_scale)) {
    p <- suppressMessages(p + color_scale)
  }

  p +
    labs(
      x = "UMAP 1",
      y = "UMAP 2",
      colour = NULL
    ) +
    theme_pub +
    theme(
      legend.position = if (show_legend) "right" else "none",
      axis.line = element_line(linewidth = 0.40, colour = graph_col),
      plot.title = element_blank()
    )
}

# ------------------------------------------------------------------------------
# 6. BUILD REVISED PANEL C | 1 dpi / 3 dpi / 7 dpi separately
# ------------------------------------------------------------------------------
# IMPORTANT:
# - All three subplots use the FULL frozen CDS, the SAME UMAP coordinates,
#   and the SAME principal graph.
# - Only the focal timepoint is highlighted; all other cells are light gray.
# - This is display-only and does not subset/relearn/reorder the trajectory.

time_vec <- as.character(SummarizedExperiment::colData(cds)$timepoint)
time_levels <- c("1 dpi", "3 dpi", "7 dpi")
background_col <- "#D9D9D9"

make_timepoint_focus_plot <- function(tp, add_panel_tag = FALSE) {

  meta_name <- paste0(
    "Figure2_focus_",
    gsub("[^A-Za-z0-9]+", "_", tp)
  )

  SummarizedExperiment::colData(cds)[[meta_name]] <- factor(
    ifelse(time_vec == tp, tp, "Other timepoints"),
    levels = c("Other timepoints", tp)
  )

  focus_values <- c(
    "Other timepoints" = background_col,
    TIME_COLORS[[tp]]
  )
  names(focus_values)[2] <- tp

  p <- make_traj_plot(
    cds,
    color_by = meta_name,
    color_scale = scale_color_manual(
      values = focus_values,
      drop = FALSE,
      guide = "none"
    ),
    show_legend = FALSE,
    cell_size = 0.30
  ) +
    labs(
      title = tp,
      tag = if (add_panel_tag) "C" else NULL
    ) +
    theme(
      plot.title = element_text(
        family = base_family,
        size = 7.2,
        face = "plain",
        hjust = 0.5,
        colour = "black",
        margin = margin(0, 0, 1.0, 0, unit = "mm")
      ),
      plot.margin = margin(2.0, 1.0, 2.0, 1.0, unit = "mm")
    )

  p
}

pC_1dpi <- make_timepoint_focus_plot("1 dpi", add_panel_tag = TRUE)
pC_3dpi <- make_timepoint_focus_plot("3 dpi", add_panel_tag = FALSE)
pC_7dpi <- make_timepoint_focus_plot("7 dpi", add_panel_tag = FALSE)

pC <- patchwork::wrap_plots(
  pC_1dpi,
  pC_3dpi,
  pC_7dpi,
  nrow = 1,
  widths = c(1, 1, 1)
)

# ------------------------------------------------------------------------------
# 7. BUILD REVISED PANEL D | corrected pseudotime only
# ------------------------------------------------------------------------------
pD_core <- make_traj_plot(
  cds,
  color_by = "step5BR2_pseudotime",
  color_scale = scale_color_gradientn(
    colours = PSEUDOTIME_COLS,
    name = "Pseudotime",
    guide = guide_colorbar(
      direction = "vertical",
      title.position = "top",
      title.hjust = 0,
      barwidth = grid::unit(2.5, "mm"),
      barheight = grid::unit(23, "mm"),
      ticks = TRUE,
      frame.colour = "#666666",
      frame.linewidth = 0.25
    )
  ),
  show_legend = TRUE,
  cell_size = 0.31
) +
  labs(tag = "D") +
  theme(
    plot.title = element_blank(),
    legend.position = "right",
    legend.title = element_text(
      family = base_family,
      size = 5.5,
      face = "plain",
      hjust = 0
    ),
    legend.text = element_text(
      family = base_family,
      size = 5.2
    ),
    legend.margin = margin(0, 0, 0, 0),
    plot.margin = margin(2.0, 0.5, 2.0, 1.0, unit = "mm")
  )

# Give Panel D a dedicated guide column so the map body is not squeezed.
pD <- patchwork::wrap_plots(
  pD_core,
  patchwork::guide_area(),
  nrow = 1,
  widths = c(1.00, 0.24),
  guides = "collect"
) &
  theme(legend.position = "right")

# Optional review preview: the new C occupies the former left side,
# while the former pseudotime right side is now standalone D.
pCD_preview <- patchwork::wrap_plots(
  pC,
  pD,
  nrow = 1,
  widths = c(3.0, 1.15)
)

# ------------------------------------------------------------------------------
# 8. SAVE HELPERS
# ------------------------------------------------------------------------------
save_panel_multiformat <- function(plot_obj, stem, width_mm, height_mm) {

  ggsave(
    file.path(out_dir, paste0(stem, ".pdf")),
    plot_obj,
    device = cairo_pdf,
    width = width_mm,
    height = height_mm,
    units = "mm",
    bg = "white"
  )

  ggsave(
    file.path(out_dir, paste0(stem, ".svg")),
    plot_obj,
    device = svglite::svglite,
    width = width_mm,
    height = height_mm,
    units = "mm",
    bg = "white"
  )

  ggsave(
    file.path(out_dir, paste0(stem, ".tiff")),
    plot_obj,
    device = ragg::agg_tiff,
    width = width_mm,
    height = height_mm,
    units = "mm",
    dpi = 600,
    compression = "lzw",
    background = "white"
  )

  ggsave(
    file.path(out_dir, paste0(stem, ".png")),
    plot_obj,
    device = ragg::agg_png,
    width = width_mm,
    height = height_mm,
    units = "mm",
    dpi = 300,
    background = "white"
  )
}

message("Saving revised Figure 2C...")
save_panel_multiformat(
  pC,
  "Figure2C_timepoint_split_v6",
  panelC_width_mm,
  panel_height_mm
)

message("Saving revised Figure 2D...")
save_panel_multiformat(
  pD,
  "Figure2D_pseudotime_v6",
  panelD_width_mm,
  panel_height_mm
)

# Review-only combined preview.
ggsave(
  file.path(out_dir, "Figure2CD_preview_v6.pdf"),
  pCD_preview,
  device = cairo_pdf,
  width = preview_width_mm,
  height = panel_height_mm,
  units = "mm",
  bg = "white"
)

ggsave(
  file.path(out_dir, "Figure2CD_preview_v6.png"),
  pCD_preview,
  device = ragg::agg_png,
  width = preview_width_mm,
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
message("Figure 2C-D v6 COMPLETE")
message("Panel C: 1 dpi / 3 dpi / 7 dpi shown separately on the SAME frozen trajectory.")
message("Panel D: corrected Step5BR2 pseudotime shown as a standalone panel.")
message("No trajectory learning, rooting, ordering, or topology was rerun.")
message("Panel C:")
message(file.path(out_dir, "Figure2C_timepoint_split_v6.pdf"))
message("Panel D:")
message(file.path(out_dir, "Figure2D_pseudotime_v6.pdf"))
message("Preview:")
message(file.path(out_dir, "Figure2CD_preview_v6.pdf"))
message("============================================================")
