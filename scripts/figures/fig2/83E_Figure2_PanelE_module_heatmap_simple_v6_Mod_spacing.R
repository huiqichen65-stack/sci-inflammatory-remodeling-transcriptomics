# ==============================================================================
# 83E_Figure2_PanelE_module_heatmap_simple_v6_Mod_spacing.R
# v5 visual-only revision: display module labels as Mod1-Mod5 (instead of M1-M5).
# Figure 2E | Simple module heatmap in the style of the provided reference
#
# Goal
# ----
# Draw a clean publication-style module heatmap similar to the reference:
# - left module labels (Mod1, Mod2, ...)
# - top 2-color stage bar with a proper legend
# - blue-white-red continuous heatmap
# - representative genes on the right
# - NO pathway/term text
#
# Frozen-data contract
# --------------------
# - Read frozen Step6B dynamic-module results only.
# - Do NOT rerun pseudotime, module detection, graph_test, or enrichment.
# - Only perform display interpolation across the x-axis to make the heatmap
#   visually smooth; this does not alter the frozen biological results.
# ==============================================================================

rm(list = ls())
gc()
options(stringsAsFactors = FALSE)

# ------------------------------------------------------------------------------
# 0. USER SETTINGS
# ------------------------------------------------------------------------------
project_root <- Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
out_dir <- file.path(
  project_root,
  "09_final_figures", "GLIA", "Figure2_panelE_simple_v6_Mod_spacing"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

panel_width_mm  <- 155
panel_height_mm <- 84

genes_per_module <- 4
display_columns <- 120

# ------------------------------------------------------------------------------
# 1. PACKAGES
# ------------------------------------------------------------------------------
required_pkgs <- c(
  "ComplexHeatmap",
  "circlize",
  "dplyr",
  "readr",
  "stringr",
  "grid",
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
  library(ComplexHeatmap)
  library(circlize)
  library(dplyr)
  library(readr)
  library(stringr)
  library(grid)
  library(scales)
})

# ------------------------------------------------------------------------------
# 2. INPUT
# ------------------------------------------------------------------------------
input_files <- c(
  step6B_rds = file.path(
    project_root,
    "04_results/step6B_GSE162610_dynamic_branch_modules",
    "GSE162610_step6B_dynamic_modules_results.rds"
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

message("Reading frozen Step6B input...")
step6B <- readRDS(input_files[["step6B_rds"]])

# ------------------------------------------------------------------------------
# 3. HELPERS
# ------------------------------------------------------------------------------
first_existing_col <- function(df, candidates, required = TRUE, context = "table") {
  hit <- intersect(candidates, colnames(df))
  if (length(hit) == 0) {
    if (required) {
      stop(
        context, " is missing expected column(s): ",
        paste(candidates, collapse = ", "),
        ". Observed columns: ",
        paste(colnames(df), collapse = ", ")
      )
    } else {
      return(NULL)
    }
  }
  hit[1]
}

safe_numeric <- function(x) {
  out <- suppressWarnings(as.numeric(x))
  if (all(is.na(out))) return(rep(NA_real_, length(x)))
  out
}

natural_order <- function(x) {
  if (length(x) <= 1) return(x)
  num <- suppressWarnings(as.numeric(stringr::str_extract(x, "(?<!\\d)\\d+(?!\\d)")))
  ord <- order(is.na(num), num, x)
  x[ord]
}

# ------------------------------------------------------------------------------
# 4. RECOVER FROZEN STEP6B OBJECTS
# ------------------------------------------------------------------------------
profile_z <- step6B$sample_balanced_profile_z
gene_modules <- as.data.frame(step6B$gene_modules)

if (is.null(profile_z) || is.null(gene_modules)) {
  stop("Step6B result lacks sample_balanced_profile_z and/or gene_modules.")
}

profile_z <- as.matrix(profile_z)

if (is.null(rownames(profile_z)) || nrow(profile_z) == 0 || ncol(profile_z) == 0) {
  stop("Frozen Step6B profile matrix is empty or lacks row names.")
}

gene_col <- first_existing_col(
  gene_modules,
  c("gene_id", "gene", "id"),
  context = "step6B$gene_modules"
)
module_col <- first_existing_col(
  gene_modules,
  c("raw_module", "module"),
  context = "step6B$gene_modules"
)
membership_col <- first_existing_col(
  gene_modules,
  c("module_membership_cor", "membership_cor", "membership"),
  required = FALSE,
  context = "step6B$gene_modules"
)
highconf_col <- first_existing_col(
  gene_modules,
  c("source_Step6A_highconf", "strict_pseudobulk_supported"),
  required = FALSE,
  context = "step6B$gene_modules"
)

gene_modules2 <- gene_modules %>%
  mutate(
    gene_id = as.character(.data[[gene_col]]),
    raw_module = as.character(.data[[module_col]]),
    membership = if (!is.null(membership_col))
      safe_numeric(.data[[membership_col]]) else NA_real_,
    highconf = if (!is.null(highconf_col))
      as.logical(.data[[highconf_col]]) else FALSE
  ) %>%
  filter(gene_id %in% rownames(profile_z)) %>%
  distinct(gene_id, .keep_all = TRUE)

if (nrow(gene_modules2) == 0) {
  stop("No frozen Step6B module genes were found in the frozen profile matrix.")
}

pre_cols <- grep("^Pre_", colnames(profile_z), value = TRUE)
target_cols <- grep("^(Arm2_|Target_)", colnames(profile_z), value = TRUE)

pre_cols <- natural_order(pre_cols)
target_cols <- natural_order(target_cols)
keep_cols <- c(pre_cols, target_cols)

if (length(pre_cols) == 0 || length(target_cols) == 0) {
  stop(
    "Frozen Step6B profile matrix must contain both Pre_ and Arm2_/Target_ bins. ",
    "Observed columns: ",
    paste(colnames(profile_z), collapse = ", ")
  )
}

# ------------------------------------------------------------------------------
# 5. PRESERVE FROZEN MODULE IDENTITIES
# ------------------------------------------------------------------------------
module_order <- natural_order(unique(as.character(gene_modules2$raw_module)))

# Preserve frozen module identities; change DISPLAY LABELS only.
# Examples: M1 -> Mod1, M2 -> Mod2; numeric raw labels 1 -> Mod1.
module_short_raw <- ifelse(
  grepl("^M", module_order, ignore.case = TRUE),
  module_order,
  paste0("M", module_order)
)

module_display <- sub(
  "^M([0-9]+)$",
  "Mod\\1",
  module_short_raw,
  perl = TRUE
)

# Fallback for any unexpected non-M<number> label:
# keep a stable Mod1, Mod2, ... display label in frozen module order.
bad_display <- !grepl("^Mod[0-9]+$", module_display)
if (any(bad_display)) {
  module_display[bad_display] <- paste0(
    "Mod",
    match(module_order[bad_display], module_order)
  )
}

module_map <- data.frame(
  raw_module = module_order,
  module_short_raw = module_short_raw,
  module_short = module_display,
  stringsAsFactors = FALSE
)

module_short_map <- module_map$module_short
names(module_short_map) <- module_map$raw_module

write.csv(
  module_map,
  file.path(out_dir, "01_module_map.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 6. SELECT REPRESENTATIVE GENES
# ------------------------------------------------------------------------------
gene_modules_ranked <- gene_modules2 %>%
  mutate(
    highconf_sort = ifelse(is.na(highconf), FALSE, highconf),
    membership_sort = ifelse(is.na(membership), -Inf, membership)
  ) %>%
  group_by(raw_module) %>%
  arrange(desc(highconf_sort), desc(membership_sort), gene_id, .by_group = TRUE) %>%
  ungroup()

display_genes <- gene_modules_ranked %>%
  group_by(raw_module) %>%
  slice_head(n = genes_per_module) %>%
  ungroup() %>%
  mutate(
    raw_module = factor(raw_module, levels = module_order)
  ) %>%
  arrange(raw_module, desc(membership_sort), gene_id) %>%
  mutate(
    module_short = module_short_map[as.character(raw_module)]
  )

if (nrow(display_genes) == 0) {
  stop("No display genes were selected for Figure 2E.")
}

write.csv(
  display_genes,
  file.path(out_dir, "02_gene_selection_audit.csv"),
  row.names = FALSE
)

heat_order <- display_genes$gene_id

# ------------------------------------------------------------------------------
# 7. BUILD HEATMAP MATRIX (DISPLAY-ONLY INTERPOLATION)
# ------------------------------------------------------------------------------
E_mat_raw <- as.matrix(
  profile_z[
    heat_order,
    keep_cols,
    drop = FALSE
  ]
)

E_mat_raw[!is.finite(E_mat_raw)] <- 0

if (
  is.null(dim(E_mat_raw)) ||
  length(dim(E_mat_raw)) != 2 ||
  nrow(E_mat_raw) == 0 ||
  ncol(E_mat_raw) == 0
) {
  stop("Figure 2E frozen heatmap matrix is invalid or empty.")
}

x_old <- seq_len(ncol(E_mat_raw))
x_new <- seq(
  from = 1,
  to = ncol(E_mat_raw),
  length.out = display_columns
)

E_mat <- t(
  apply(
    E_mat_raw,
    1,
    function(y) {
      stats::approx(
        x = x_old,
        y = as.numeric(y),
        xout = x_new,
        method = "linear",
        ties = "ordered"
      )$y
    }
  )
)

rownames(E_mat) <- rownames(E_mat_raw)
colnames(E_mat) <- paste0("display_", seq_len(ncol(E_mat)))

E_mat[E_mat >  2.5] <-  2.5
E_mat[E_mat < -2.5] <- -2.5

E_row_split <- factor(
  display_genes$module_short,
  levels = module_short_map[module_order]
)

if (length(E_row_split) != nrow(E_mat)) {
  stop("Figure 2E row split length does not match heatmap rows.")
}

write.csv(
  data.frame(
    item = c("frozen_columns", "display_columns"),
    value = c(ncol(E_mat_raw), ncol(E_mat))
  ),
  file.path(out_dir, "03_display_interpolation_audit.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 8. DISPLAY STYLE
# ------------------------------------------------------------------------------
BASE_FAMILY <- "Arial"

# Top 2-color annotation and legend
STAGE_COLORS <- c(
  "Shared prebranch" = "#8B86B8",
  "Target-program branch" = "#C9675C"
)

module_fill_palette <- c(
  "#8794B8",  # Mod1 slate blue
  "#8794B8",  # Mod2 same tone family
  "#A88BAF",  # Mod3 mauve
  "#A88BAF",  # Mod4 mauve
  "#8794B8",  # Mod5 slate blue
  "#8794B8",
  "#A88BAF"
)

module_fill_map <- setNames(
  module_fill_palette[seq_along(levels(E_row_split))],
  levels(E_row_split)
)

# smooth blue-white-red
z_col_fun <- circlize::colorRamp2(
  c(-2.5, 0, 2.5),
  c("#4F6DAF", "#F7F7F7", "#C85B55")
)

stage_factor <- factor(
  c(rep("Shared prebranch", length(pre_cols)),
    rep("Target-program branch", length(target_cols))),
  levels = c("Shared prebranch", "Target-program branch")
)

# expand stage annotation to display_columns
stage_display <- c(
  rep("Shared prebranch", round(display_columns * length(pre_cols) / length(keep_cols))),
  rep("Target-program branch", display_columns - round(display_columns * length(pre_cols) / length(keep_cols)))
)
stage_display <- factor(
  stage_display,
  levels = c("Shared prebranch", "Target-program branch")
)

top_anno <- HeatmapAnnotation(
  Stage = stage_display,
  col = list(Stage = STAGE_COLORS),
  show_annotation_name = FALSE,
  show_legend = FALSE,
  simple_anno_size = unit(4.0, "mm"),
  gp = gpar(col = NA)
)

# Explicit legends are drawn manually in a dedicated right-side area.
# This guarantees that the Stage text is complete and that Stage sits beside Row z.
rowz_legend <- ComplexHeatmap::Legend(
  title = "Row z",
  col_fun = z_col_fun,
  at = c(-2.5, 0, 2.5),
  labels = c("-2.5", "0", "2.5"),
  direction = "vertical",
  legend_height = unit(24, "mm"),
  grid_width = unit(3.2, "mm"),
  labels_gp = gpar(
    fontsize = 7.2,
    fontfamily = BASE_FAMILY
  ),
  title_gp = gpar(
    fontsize = 7.8,
    fontfamily = BASE_FAMILY,
    fontface = "bold"
  ),
  border = NA
)

stage_legend <- ComplexHeatmap::Legend(
  title = "Stage",
  labels = c("Shared prebranch", "Target-program branch"),
  legend_gp = gpar(
    fill = unname(STAGE_COLORS[
      c("Shared prebranch", "Target-program branch")
    ]),
    col = NA
  ),
  ncol = 1,
  grid_width = unit(4.0, "mm"),
  grid_height = unit(2.8, "mm"),
  labels_gp = gpar(
    fontsize = 7.2,
    fontfamily = BASE_FAMILY
  ),
  title_gp = gpar(
    fontsize = 7.8,
    fontfamily = BASE_FAMILY,
    fontface = "bold"
  ),
  gap = unit(1.6, "mm")
)

combined_legend <- ComplexHeatmap::packLegend(
  rowz_legend,
  stage_legend,
  direction = "vertical",
  gap = unit(5.0, "mm")
)

left_anno <- rowAnnotation(
  Module = anno_block(
    # No colored blocks: keep only small module labels.
    gp = gpar(
      fill = NA,
      col = NA
    ),
    labels = levels(E_row_split),
    labels_gp = gpar(
      fontsize = 8.0,
      fontfamily = BASE_FAMILY,
      fontface = "bold",
      col = "black"
    ),
    labels_rot = 0
  ),
  show_annotation_name = FALSE,
  width = unit(9.0, "mm")
)

right_anno <- rowAnnotation(
  Gene = anno_text(
    heat_order,
    gp = gpar(
      fontsize = 7.6,
      fontfamily = BASE_FAMILY,
      fontface = "italic",
      col = "black"
    ),
    just = "left",
    location = unit(0, "npc"),
    width = unit(30, "mm")
  ),
  show_annotation_name = FALSE
)

E_ht <- Heatmap(
  E_mat,
  name = "Row z",
  col = z_col_fun,
  show_heatmap_legend = FALSE,

  cluster_rows = FALSE,
  cluster_columns = FALSE,
  cluster_row_slices = FALSE,

  row_split = E_row_split,
  row_gap = unit(1.2, "mm"),

  top_annotation = top_anno,
  left_annotation = left_anno,
  right_annotation = right_anno,

  show_row_names = FALSE,
  show_column_names = FALSE,

  row_title_gp = gpar(fontsize = 0.1, col = "transparent"),

  border = FALSE,
  rect_gp = gpar(col = NA),

  use_raster = TRUE,
  raster_quality = 5,
  raster_resize_mat = FALSE,

  width = unit(70, "mm"),

  heatmap_legend_param = list(
    at = c(-2.5, 0, 2.5),
    labels = c("-2.5", "0", "2.5"),
    title = "Row z",
    title_position = "topcenter",
    title_gp = gpar(
      fontsize = 8.0,
      fontfamily = BASE_FAMILY,
      fontface = "bold"
    ),
    labels_gp = gpar(
      fontsize = 7.2,
      fontfamily = BASE_FAMILY
    ),
    grid_width = unit(3.2, "mm"),
    legend_height = unit(22, "mm"),
    border = NA
  )
)

# ------------------------------------------------------------------------------
# 9. DRAW PANEL
# ------------------------------------------------------------------------------
draw_panel_E <- function() {

  # --------------------------------------------------------------------------
  # Fixed heatmap viewport
  # --------------------------------------------------------------------------
  # Do not rely on ComplexHeatmap padding for the outer page geometry.
  # The whole heatmap is deliberately shifted right and down, leaving:
  #   - a true blank margin on the left;
  #   - a true blank upper-left corner for panel tag E;
  #   - a dedicated right zone for the two legends.
  pushViewport(
    viewport(
      x = unit(11, "mm"),
      y = unit(4.0, "mm"),
      width = unit(108, "mm"),
      height = unit(panel_height_mm - 10, "mm"),
      just = c("left", "bottom"),
      clip = "off",
      name = "figure2E_heatmap_region"
    )
  )

  ht_drawn <- draw(
    E_ht,
    heatmap_legend_side = "right",
    annotation_legend_side = "right",
    merge_legend = FALSE,
    padding = unit(c(1.0, 1.0, 1.0, 1.0), "mm"),
    newpage = FALSE
  )

  popViewport()

  # --------------------------------------------------------------------------
  # Manual legend zone
  # --------------------------------------------------------------------------
  # Row z is placed above Stage. The fixed viewport ensures
  # "Shared prebranch" and "Target-program branch" are never cropped.
  pushViewport(
    viewport(
      x = unit(116, "mm"),
      y = unit(42, "mm"),
      width = unit(37, "mm"),
      height = unit(58, "mm"),
      just = c("left", "center"),
      clip = "off",
      name = "figure2E_legend_region"
    )
  )

  ComplexHeatmap::draw(
    combined_legend,
    x = unit(0, "npc"),
    y = unit(0.5, "npc"),
    just = c("left", "center")
  )

  popViewport()

  # --------------------------------------------------------------------------
  # Panel tag
  # --------------------------------------------------------------------------
  # E is now in page-level whitespace, not inside the heatmap viewport.
  grid.text(
    "E",
    x = unit(2.0, "mm"),
    y = unit(1, "npc") - unit(2.0, "mm"),
    just = c("left", "top"),
    gp = gpar(
      fontsize = 16,
      fontface = "bold",
      fontfamily = BASE_FAMILY,
      col = "black"
    )
  )

  invisible(ht_drawn)
}

# ------------------------------------------------------------------------------
# 10. SAVE
# ------------------------------------------------------------------------------
stem <- "Figure2E_module_heatmap_simple_v6_Mod_spacing"

message(
  "Rendering with font family: ", BASE_FAMILY,
  "; Cairo support: ", capabilities("cairo")
)

# PDF
if (!capabilities("cairo")) {
  stop(
    "This R installation does not report Cairo support. ",
    "Figure 2E requires cairo_pdf() for reliable Arial rendering."
  )
}

grDevices::cairo_pdf(
  filename = file.path(out_dir, paste0(stem, ".pdf")),
  width = panel_width_mm / 25.4,
  height = panel_height_mm / 25.4,
  family = BASE_FAMILY,
  onefile = FALSE,
  bg = "white"
)
grid.newpage()
draw_panel_E()
dev.off()

# SVG
svglite::svglite(
  file = file.path(out_dir, paste0(stem, ".svg")),
  width = panel_width_mm / 25.4,
  height = panel_height_mm / 25.4,
  bg = "white"
)
grid.newpage()
draw_panel_E()
dev.off()

# TIFF
ragg::agg_tiff(
  filename = file.path(out_dir, paste0(stem, ".tiff")),
  width = panel_width_mm,
  height = panel_height_mm,
  units = "mm",
  res = 600,
  compression = "lzw",
  background = "white"
)
grid.newpage()
draw_panel_E()
dev.off()

# PNG
ragg::agg_png(
  filename = file.path(out_dir, paste0(stem, ".png")),
  width = panel_width_mm,
  height = panel_height_mm,
  units = "mm",
  res = 300,
  background = "white"
)
grid.newpage()
draw_panel_E()
dev.off()

capture.output(
  sessionInfo(),
  file = file.path(out_dir, "sessionInfo.txt")
)

message("")
message("============================================================")
message("Figure 2E simple module heatmap v5 Mod-label revision COMPLETE")
message("Primary review file:")
message(file.path(out_dir, paste0(stem, ".pdf")))
message("============================================================")
