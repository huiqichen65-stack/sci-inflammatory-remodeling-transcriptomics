# ==============================================================================
# 83A_Figure2_PanelA_resident_AUCell_bridge_v6_GLIA_countfix.R
# v6 keeps the v5 visual revision and explicitly uses dplyr::count()
# to avoid count() namespace masking/conflicts during master assembly.
# Figure 2A | Resident-state landscape with Activated and Early-DAM programs
#
# Biological purpose
# ------------------
# Bridge Figure 1 -> Figure 2:
#   A1. Figure-1-defined resident microglial states
#   A2. Activated transcriptional program (frozen AUCell AUC)
#   A3. Early-DAM transcriptional program (frozen AUCell AUC)
#
# IMPORTANT
# ---------
# - Read frozen RDS only.
# - Never rerun QC, UMAP, clustering, state annotation, or AUCell.
# - A1 uses EXACT Figure 1 state mapping:
#       cluster 0   -> Canonical resident
#       clusters 3/4 -> Activated resident
#       cluster 9   -> IFN-responsive resident
#       clusters 5/7 -> Cycling resident
# - All three maps use the SAME frozen "umap.microglia.step3A" coordinates.
# - Step3B is the Figure 1 anchor; Step3C supplies only frozen AUCell scores.
# - If Step3B and Step3C cell sets or UMAP coordinates disagree, STOP.
#
# Output
# ------
# <PROJECT_ROOT>/09_final_figures/GLIA/Figure2_panelA_v4/
#   Figure2A_resident_AUCell_bridge_v4.pdf
#   Figure2A_resident_AUCell_bridge_v4.svg
#   Figure2A_resident_AUCell_bridge_v4.tiff
#   Figure2A_resident_AUCell_bridge_v4.png
#   00_input_audit.csv
#   01_alignment_audit.csv
#   02_state_counts.csv
#   03_AUCell_ranges.csv
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
  "09_final_figures", "GLIA", "Figure2_panelA_v6_GLIA_countfix"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Standalone Panel A size.
# Designed from the start as one publication panel containing 3 matched UMAPs.
panel_width_mm  <- 164
panel_height_mm <- 48

# ------------------------------------------------------------------------------
# 1. PACKAGES
# ------------------------------------------------------------------------------
required_pkgs <- c(
  "Seurat",
  "ggplot2",
  "dplyr",
  "patchwork",
  "scales",
  "ggrastr",
  "ggrepel",
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
  library(patchwork)
  library(scales)
  library(ggrepel)
})

# ------------------------------------------------------------------------------
# 2. FROZEN INPUTS
# ------------------------------------------------------------------------------
input_files <- c(
  # Exact Figure 1G source object.
  step3B_microglia = file.path(
    project_root,
    "04_results/step3B_GSE162610_state_scoring",
    "GSE162610_step3B_clean_state_scored.rds"
  ),

  # Frozen AUCell scores used by Figure 2.
  step3C_transition = file.path(
    project_root,
    "04_results/step3C_GSE162610_activated_DAM_transition_mapping",
    "GSE162610_step3C_transition_mapped.rds"
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

message("Reading frozen Figure 1 / Figure 2 bridge inputs...")
obj3B <- readRDS(input_files[["step3B_microglia"]])
obj3C <- readRDS(input_files[["step3C_transition"]])

if (!inherits(obj3B, "Seurat") || !inherits(obj3C, "Seurat")) {
  stop("Both frozen inputs must be Seurat objects.")
}

# ------------------------------------------------------------------------------
# 3. SCHEMA + PROVENANCE AUDIT
# ------------------------------------------------------------------------------
reduction_name <- "umap.microglia.step3A"

if (!reduction_name %in% Reductions(obj3B)) {
  stop("Step3B object lacks frozen reduction: ", reduction_name)
}
if (!reduction_name %in% Reductions(obj3C)) {
  stop("Step3C object lacks frozen reduction: ", reduction_name)
}

required_3B <- c(
  "step2_final_cluster",
  "sample_id",
  "timepoint"
)

required_3C <- c(
  "TransAUC_Activated_injury_core",
  "TransAUC_Early_DAM_transition"
)

miss3B <- setdiff(required_3B, colnames(obj3B[[]]))
miss3C <- setdiff(required_3C, colnames(obj3C[[]]))

if (length(miss3B) > 0) {
  stop("Step3B metadata missing: ", paste(miss3B, collapse = ", "))
}
if (length(miss3C) > 0) {
  stop("Step3C metadata missing: ", paste(miss3C, collapse = ", "))
}

cells3B <- colnames(obj3B)
cells3C <- colnames(obj3C)

same_cell_set <- identical(sort(cells3B), sort(cells3C))

if (!same_cell_set) {
  only3B <- setdiff(cells3B, cells3C)
  only3C <- setdiff(cells3C, cells3B)

  alignment_audit <- data.frame(
    item = c(
      "n_Step3B",
      "n_Step3C",
      "Step3B_only",
      "Step3C_only"
    ),
    value = c(
      length(cells3B),
      length(cells3C),
      length(only3B),
      length(only3C)
    )
  )

  write.csv(
    alignment_audit,
    file.path(out_dir, "01_alignment_audit.csv"),
    row.names = FALSE
  )

  stop(
    "Step3B and Step3C cell sets are not identical. ",
    "Figure 2A must use the exact Figure 1 resident landscape. ",
    "See 01_alignment_audit.csv."
  )
}

# Reorder Step3C to Figure1/Step3B order.
obj3C <- obj3C[, cells3B]

emb3B <- Embeddings(obj3B, reduction = reduction_name)
emb3C <- Embeddings(obj3C, reduction = reduction_name)

emb3B <- emb3B[cells3B, 1:2, drop = FALSE]
emb3C <- emb3C[cells3B, 1:2, drop = FALSE]

coord_abs_diff <- abs(emb3B - emb3C)
max_coord_diff <- max(coord_abs_diff, na.rm = TRUE)

alignment_audit <- data.frame(
  item = c(
    "n_cells",
    "same_cell_set",
    "max_absolute_UMAP_coordinate_difference"
  ),
  value = c(
    as.character(length(cells3B)),
    as.character(same_cell_set),
    format(max_coord_diff, scientific = TRUE)
  ),
  stringsAsFactors = FALSE
)

write.csv(
  alignment_audit,
  file.path(out_dir, "01_alignment_audit.csv"),
  row.names = FALSE
)

# Downstream object should preserve the exact frozen embedding.
# A tiny floating-point tolerance is allowed.
if (!is.finite(max_coord_diff) || max_coord_diff > 1e-8) {
  stop(
    "Step3B and Step3C frozen UMAP coordinates differ (max |delta| = ",
    signif(max_coord_diff, 4),
    "). Stop rather than silently mixing embeddings."
  )
}

# ------------------------------------------------------------------------------
# 4. EXACT FIGURE 1 RESIDENT-STATE MAPPING
# ------------------------------------------------------------------------------
resident_state_cols <- c(
  "Canonical resident"      = "#4C78A8",
  "Activated resident"      = "#E45756",
  "IFN-responsive resident" = "#1B9E77",
  "Cycling resident"        = "#B279A2"
)

collapse_resident_state <- function(step2_cluster) {
  x <- as.character(step2_cluster)

  dplyr::case_when(
    x == "0"           ~ "Canonical resident",
    x %in% c("3", "4") ~ "Activated resident",
    x == "9"           ~ "IFN-responsive resident",
    x %in% c("5", "7") ~ "Cycling resident",
    TRUE               ~ NA_character_
  )
}

meta3B <- obj3B[[]]
meta3C <- obj3C[[]]

resident_state <- collapse_resident_state(meta3B$step2_final_cluster)

if (anyNA(resident_state)) {
  bad <- sort(
    unique(
      as.character(
        meta3B$step2_final_cluster[is.na(resident_state)]
      )
    )
  )

  stop(
    "Unexpected Step2 final cluster(s): ",
    paste(bad, collapse = ", "),
    ". No Figure 1 state label was guessed."
  )
}

resident_state <- factor(
  resident_state,
  levels = names(resident_state_cols)
)

# ------------------------------------------------------------------------------
# 5. MASTER PLOT DATA
# ------------------------------------------------------------------------------
plot_df <- data.frame(
  cell_id = cells3B,
  UMAP_1 = emb3B[, 1],
  UMAP_2 = emb3B[, 2],

  resident_state = resident_state,

  Activated_AUCell = as.numeric(
    meta3C[cells3B, "TransAUC_Activated_injury_core"]
  ),
  EarlyDAM_AUCell = as.numeric(
    meta3C[cells3B, "TransAUC_Early_DAM_transition"]
  ),

  stringsAsFactors = FALSE
)

if (
  any(!is.finite(plot_df$Activated_AUCell)) ||
  any(!is.finite(plot_df$EarlyDAM_AUCell))
) {
  stop("Non-finite AUCell values detected in frozen Step3C metadata.")
}

# Same deterministic point order for ALL three maps.
set.seed(20260813)
point_order <- sample(seq_len(nrow(plot_df)))
plot_df <- plot_df[point_order, , drop = FALSE]

state_counts <- plot_df %>%
  dplyr::count(resident_state, name = "n_cells")

write.csv(
  state_counts,
  file.path(out_dir, "02_state_counts.csv"),
  row.names = FALSE
)

score_ranges <- data.frame(
  program = c("Activated", "Early-DAM"),
  minimum = c(
    min(plot_df$Activated_AUCell),
    min(plot_df$EarlyDAM_AUCell)
  ),
  median = c(
    median(plot_df$Activated_AUCell),
    median(plot_df$EarlyDAM_AUCell)
  ),
  maximum = c(
    max(plot_df$Activated_AUCell),
    max(plot_df$EarlyDAM_AUCell)
  ),
  stringsAsFactors = FALSE
)

write.csv(
  score_ranges,
  file.path(out_dir, "03_AUCell_ranges.csv"),
  row.names = FALSE
)

layout_audit <- data.frame(
  subpanel = c("A1", "A2", "A3"),
  UMAP_core_width_unit = c(1, 1, 1),
  UMAP_core_height_unit = c(1, 1, 1),
  right_legend_width_unit = c(0, 0.30, 0.30),
  stringsAsFactors = FALSE
)

write.csv(
  layout_audit,
  file.path(out_dir, "04_equal_UMAP_layout_audit.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 6. FIGURE 1-MATCHED GLOBAL STYLE
# ------------------------------------------------------------------------------
base_family <- "Arial"

theme_pub <- theme_classic(
  base_family = base_family,
  base_size = 8.0
) +
  theme(
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),

    axis.title = element_text(
      family = base_family,
      size = 7.8,
      colour = "black"
    ),
    axis.text = element_text(
      family = base_family,
      size = 7.2,
      colour = "black"
    ),

    legend.title = element_text(
      family = base_family,
      size = 6.6,
      colour = "black"
    ),
    legend.text = element_text(
      family = base_family,
      size = 6.2,
      colour = "black"
    ),

    plot.tag = element_text(
      family = base_family,
      face = "bold",
      size = 16,
      colour = "black"
    ),
    plot.tag.position = c(-0.035, 1.035),

    plot.margin = margin(2.5, 2.0, 2.5, 2.0, unit = "mm")
  )

theme_umap <- theme_pub +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank()
  )

# ------------------------------------------------------------------------------
# 7. UMAP HELPERS
# ------------------------------------------------------------------------------

# Same Figure 1 lower-left vector axes, slightly tuned for the smaller Figure 2A maps.
add_umap_axes <- function(
  p,
  df,
  x_fraction = 0.285,
  y_fraction = 0.285
) {
  xr <- range(df$UMAP_1, na.rm = TRUE)
  yr <- range(df$UMAP_2, na.rm = TRUE)
  dx <- diff(xr)
  dy <- diff(yr)

  # V2:
  # Move the coordinate pocket clearly away from the cell cloud.
  # This prevents the Figure1-style axes from touching the lower-left UMAP body.
  x0 <- xr[1] - 0.235 * dx
  y0 <- yr[1] - 0.255 * dy

  x1 <- x0 + x_fraction * dx
  y1 <- y0 + y_fraction * dy

  # More white space at lower-left so the enlarged arrows/text remain isolated.
  x_limits <- c(
    xr[1] - 0.405 * dx,
    xr[2] + 0.025 * dx
  )
  y_limits <- c(
    yr[1] - 0.435 * dy,
    yr[2] + 0.045 * dy
  )

  p +
    scale_x_continuous(
      limits = x_limits,
      expand = expansion(mult = 0)
    ) +
    scale_y_continuous(
      limits = y_limits,
      expand = expansion(mult = 0)
    ) +
    coord_equal(clip = "off") +

    annotate(
      "segment",
      x = x0, xend = x1,
      y = y0, yend = y0,
      linewidth = 0.82,
      colour = "black",
      lineend = "round",
      arrow = grid::arrow(
        type = "closed",
        length = grid::unit(2.05, "mm")
      )
    ) +

    annotate(
      "segment",
      x = x0, xend = x0,
      y = y0, yend = y1,
      linewidth = 0.82,
      colour = "black",
      lineend = "round",
      arrow = grid::arrow(
        type = "closed",
        length = grid::unit(2.05, "mm")
      )
    ) +

    annotate(
      "text",
      x = (x0 + x1) / 2,
      y = y0 - 0.088 * dy,
      label = "UMAP_1",
      family = base_family,
      size = 2.15,
      hjust = 0.5,
      vjust = 1
    ) +

    annotate(
      "text",
      x = x0 - 0.094 * dx,
      y = (y0 + y1) / 2,
      label = "UMAP_2",
      family = base_family,
      size = 2.15,
      angle = 90,
      hjust = 0.5,
      vjust = 0.5
    )
}

label_centers <- function(df) {
  # Use the largest contributing frozen Step2 cluster for each display state,
  # matching the Figure 1 logic and avoiding labels in empty geometric centers.
  tmp <- df %>%
    mutate(
      step2_final_cluster = as.character(
        meta3B[cell_id, "step2_final_cluster"]
      )
    )

  largest_subgroup <- tmp %>%
    dplyr::count(
      resident_state,
      step2_final_cluster,
      name = "n_cells"
    ) %>%
    group_by(resident_state) %>%
    slice_max(
      order_by = n_cells,
      n = 1,
      with_ties = FALSE
    ) %>%
    ungroup() %>%
    select(resident_state, step2_final_cluster)

  tmp %>%
    inner_join(
      largest_subgroup,
      by = c("resident_state", "step2_final_cluster")
    ) %>%
    group_by(resident_state) %>%
    summarise(
      UMAP_1 = median(UMAP_1),
      UMAP_2 = median(UMAP_2),
      .groups = "drop"
    ) %>%
    mutate(
      label = dplyr::recode(
        as.character(resident_state),
        "Canonical resident" = "Canonical\nresident",
        "Activated resident" = "Activated\nresident",
        "IFN-responsive resident" = "IFN-responsive\nresident",
        "Cycling resident" = "Cycling\nresident"
      )
    )
}

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
# 8. A1 | FIGURE-1-DEFINED RESIDENT STATES
# ------------------------------------------------------------------------------
centers <- label_centers(plot_df)

pA1 <- ggplot(
  plot_df,
  aes(
    x = UMAP_1,
    y = UMAP_2,
    colour = resident_state
  )
) +
  geom_cells(
    size = 0.22,
    alpha = 0.88
  ) +
  scale_colour_manual(
    values = resident_state_cols,
    drop = FALSE
  ) +
  ggrepel::geom_text_repel(
    data = centers,
    aes(
      x = UMAP_1,
      y = UMAP_2,
      label = label
    ),
    inherit.aes = FALSE,
    family = base_family,
    size = 2.18,
    colour = "black",
    fontface = "plain",
    seed = 20260813,
    max.overlaps = Inf,
    box.padding = 0.18,
    point.padding = 0.04,
    min.segment.length = Inf,
    segment.colour = NA,
    force = 0.65
  ) +
  labs(
    tag = "A"
  ) +
  theme_umap +
  theme(
    legend.position = "none"
  )

pA1 <- add_umap_axes(pA1, plot_df)

# ------------------------------------------------------------------------------
# 9. A2/A3 | FROZEN AUCELL PROGRAM MAPS
# ------------------------------------------------------------------------------

# Same restrained palette for both programs.
# Each program keeps its OWN numeric range, avoiding an implied direct
# cross-gene-set magnitude comparison.
aucell_cols <- c(
  "#ECECEC",
  "#D8D0EA",
  "#B39AD4",
  "#8064B6",
  "#4D347F"
)

make_aucell_umap <- function(
  df,
  score_col,
  program_name
) {
  score <- df[[score_col]]
  lims <- range(score, finite = TRUE)

  if (diff(lims) == 0) {
    lims <- lims + c(-1e-6, 1e-6)
  }

  ggplot(
    df,
    aes(
      x = UMAP_1,
      y = UMAP_2,
      colour = .data[[score_col]]
    )
  ) +
    geom_cells(
      size = 0.22,
      alpha = 0.92
    ) +
    scale_colour_gradientn(
      colours = aucell_cols,
      limits = lims,
      oob = scales::squish,
      name = paste0(program_name, "\nAUCell AUC"),
      guide = guide_colourbar(
        direction = "vertical",
        title.position = "top",
        title.hjust = 0,
        barwidth = grid::unit(2.0, "mm"),
        barheight = grid::unit(19.0, "mm"),
        ticks = TRUE,
        frame.colour = "#6A6A6A",
        frame.linewidth = 0.25
      )
    ) +
    theme_umap +
    theme(
      # No top plot title in V2.
      plot.title = element_blank(),

      legend.position = "right",
      legend.justification = c(0.5, 0.5),

      legend.title = element_text(
        family = base_family,
        size = 6.2,
        face = "plain",
        colour = "black",
        hjust = 0,
        lineheight = 0.95,
        margin = margin(b = 1.0, unit = "mm")
      ),

      legend.text = element_text(
        family = base_family,
        size = 5.7,
        colour = "black"
      ),

      legend.margin = margin(0.5, 1.5, 0.5, 0.6, unit = "mm"),
      legend.box.margin = margin(0, 0, 0, 0, unit = "mm")
    )
}

pA2_core <- make_aucell_umap(
  plot_df,
  score_col = "Activated_AUCell",
  program_name = "Activated program"
)
pA2_core <- add_umap_axes(pA2_core, plot_df)

pA3_core <- make_aucell_umap(
  plot_df,
  score_col = "EarlyDAM_AUCell",
  program_name = "Early-DAM program"
)
pA3_core <- add_umap_axes(pA3_core, plot_df)

# ------------------------------------------------------------------------------
# 10. COMPOSITE PANEL A
# ------------------------------------------------------------------------------
# V3 key layout rule:
#   - A1, A2_core and A3_core each occupy EXACTLY one equal UMAP-width unit.
#   - A2/A3 legends live in their own narrow guide columns.
#   - Therefore all three UMAP bodies have the same width AND the same height.

pA2 <- patchwork::wrap_plots(
  pA2_core,
  patchwork::guide_area(),
  nrow = 1,
  widths = c(1.00, 0.30),
  guides = "collect"
) &
  theme(legend.position = "right")

pA3 <- patchwork::wrap_plots(
  pA3_core,
  patchwork::guide_area(),
  nrow = 1,
  widths = c(1.00, 0.30),
  guides = "collect"
) &
  theme(legend.position = "right")

pA <- patchwork::wrap_plots(
  pA1,
  pA2,
  pA3,
  nrow = 1,
  # Total widths account for the narrow legend columns.
  # The UMAP body within each block remains exactly 1.00 unit wide.
  widths = c(1.00, 1.30, 1.30)
)

# ------------------------------------------------------------------------------
# 11. SAVE
# ------------------------------------------------------------------------------
stem <- "Figure2A_resident_AUCell_bridge_v6_GLIA_countfix"

ggsave(
  file.path(out_dir, paste0(stem, ".pdf")),
  pA,
  device = cairo_pdf,
  width = panel_width_mm,
  height = panel_height_mm,
  units = "mm",
  bg = "white"
)

ggsave(
  file.path(out_dir, paste0(stem, ".svg")),
  pA,
  device = svglite::svglite,
  width = panel_width_mm,
  height = panel_height_mm,
  units = "mm",
  bg = "white"
)

ggsave(
  file.path(out_dir, paste0(stem, ".tiff")),
  pA,
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
  pA,
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
message("Figure 2A v6 GLIA count-fix COMPLETE")
message("Biological message:")
message("Figure-1-defined resident states -> Activated / Early-DAM AUCell programs")
message("")
message("Primary review file:")
message(file.path(out_dir, paste0(stem, ".pdf")))
message("============================================================")
