# ============================================================
# v7 GLIA readability revision: enlarged axis/gene/legend text; violin data and grouping unchanged.
# Figure 4D v6 | Side-by-side FULL violin plots (not split halves)
#
# USER-REQUESTED FIXES
#   1) Left and right should each be a COMPLETE violin, not half violins
#   2) Legend should be on the RIGHT side of the whole figure
#   3) Left y-axis numbers should not overlap
#
# DATA
#   - RNA normalized expression from Step1F Seurat atlas
#   - standardized metadata from frozen Step8A manifest
#
# OUTPUT
#   - one-page PDF
#   - one PNG preview
# ============================================================

options(stringsAsFactors = FALSE)

# ------------------------------------------------------------
# 0. Packages
# ------------------------------------------------------------
pkgs <- c("Seurat", "ggplot2", "dplyr", "tibble", "grid", "scales")
miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]

if (length(miss) > 0) {
  stop(
    "Missing package(s): ",
    paste(miss, collapse = ", "),
    "\nInstall with:\ninstall.packages(c(",
    paste(sprintf('"%s"', miss), collapse = ", "),
    "))"
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(grid)
  library(scales)
})

# ------------------------------------------------------------
# 1. Paths / settings
# ------------------------------------------------------------
ROOT <- Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
ATLAS_RDS_FILE <- file.path(
  ROOT,
  "04_results",
  "step1F_GSE162610_broad_annotation",
  "GSE162610_broad_annotation_audit.rds"
)

STEP8A_MANIFEST_FILE <- file.path(
  ROOT,
  "04_results",
  "step8A_GSE162610_CellChat_input_preflight",
  "RDS",
  "GSE162610_step8A_CellChat_input_manifest.rds"
)

OUTDIR <- file.path(
  ROOT,
  "09_final_figures",
  "GLIA",
  "Figure4",
  "single_panels",
  "Figure4D"
)

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

PDF_FILE <- file.path(
  OUTDIR,
  "Figure4D_reference_style_full_violin_v7_GLIA_readability.pdf"
)

PNG_FILE <- file.path(
  OUTDIR,
  "Figure4D_reference_style_full_violin_v7_GLIA_readability.png"
)

if (file.exists(PDF_FILE)) unlink(PDF_FILE)
if (file.exists(PNG_FILE)) unlink(PNG_FILE)

GENES_TO_PLOT <- c("Tnf", "Icam1", "Ccl2")
TIMEPOINTS_TO_PLOT <- c("1dpi", "7dpi")

TIMEPOINT_COLORS <- c(
  "1dpi" = "#EA766A",
  "7dpi" = "#11A7B0"
)

GROUP_ORDER_RAW <- c(
  "Microglia_Arm_2",
  "Other_myeloid",
  "Microglia_other",
  "Neutrophil",
  "T_cell",
  "Astrocyte",
  "Endothelial",
  "Fibroblast",
  "Pericyte_VSMC",
  "OPC"
)

GROUP_DISPLAY <- c(
  "Microglia_Arm_2" = "Target program",
  "Other_myeloid"   = "Other myeloid",
  "Microglia_other" = "Microglia other",
  "Neutrophil"      = "Neutrophil",
  "T_cell"          = "T cell",
  "Astrocyte"       = "Astrocyte",
  "Endothelial"     = "Endothelial",
  "Fibroblast"      = "Fibroblast",
  "Pericyte_VSMC"   = "Pericyte/VSMC",
  "OPC"             = "OPC"
)

FONT <- "Arial"

if (.Platform$OS.type == "windows") {
  try(
    grDevices::windowsFonts(
      Arial = grDevices::windowsFont("Arial")
    ),
    silent = TRUE
  )
}

# ------------------------------------------------------------
# 2. Load atlas + standardized metadata
# ------------------------------------------------------------
if (!file.exists(ATLAS_RDS_FILE)) {
  stop("Atlas not found:\n", ATLAS_RDS_FILE)
}

if (!file.exists(STEP8A_MANIFEST_FILE)) {
  stop("Step8A manifest not found:\n", STEP8A_MANIFEST_FILE)
}

atlas <- readRDS(ATLAS_RDS_FILE)
if (!inherits(atlas, "Seurat")) {
  stop("Atlas RDS is not a Seurat object.")
}

manifest <- readRDS(STEP8A_MANIFEST_FILE)
if (!is.list(manifest) || is.null(manifest$cell_metadata)) {
  stop("Step8A manifest$cell_metadata is missing.")
}

metadata <- as.data.frame(manifest$cell_metadata, stringsAsFactors = FALSE)

required_meta <- c(
  "cell",
  "sample_id",
  "timepoint_standard",
  "cellchat_group",
  "include_step8_timecourse"
)

missing_meta <- setdiff(required_meta, colnames(metadata))
if (length(missing_meta) > 0) {
  stop(
    "Step8A cell_metadata missing: ",
    paste(missing_meta, collapse = ", ")
  )
}

as_logical_safe <- function(x) {
  if (is.logical(x)) return(x)
  if (is.numeric(x) || is.integer(x)) return(!is.na(x) & x != 0)

  y <- trimws(tolower(as.character(x)))
  out <- rep(NA, length(y))
  out[y %in% c("true", "t", "1", "yes", "y")] <- TRUE
  out[y %in% c("false", "f", "0", "no", "n")] <- FALSE
  out
}

metadata <- metadata %>%
  mutate(
    cell = as.character(cell),
    sample_id = as.character(sample_id),
    timepoint_standard = as.character(timepoint_standard),
    cellchat_group = as.character(cellchat_group),
    include_step8_timecourse = as_logical_safe(include_step8_timecourse),
    group_display = unname(GROUP_DISPLAY[cellchat_group])
  ) %>%
  filter(
    include_step8_timecourse %in% TRUE,
    timepoint_standard %in% TIMEPOINTS_TO_PLOT,
    cellchat_group %in% GROUP_ORDER_RAW
  )

if (nrow(metadata) == 0) {
  stop("No cells remained after metadata filtering.")
}

# ------------------------------------------------------------
# 3. RNA normalized expression
# ------------------------------------------------------------
Seurat::DefaultAssay(atlas) <- "RNA"

expr_mat <- try(
  Seurat::GetAssayData(
    atlas,
    assay = "RNA",
    layer = "data"
  ),
  silent = TRUE
)

if (inherits(expr_mat, "try-error")) {
  expr_mat <- Seurat::GetAssayData(
    atlas,
    assay = "RNA",
    slot = "data"
  )
}

missing_genes <- setdiff(GENES_TO_PLOT, rownames(expr_mat))
if (length(missing_genes) > 0) {
  stop(
    "Missing requested gene(s) from RNA assay: ",
    paste(missing_genes, collapse = ", ")
  )
}

cells_use <- intersect(metadata$cell, colnames(expr_mat))
if (length(cells_use) == 0) {
  stop("No Step8A cell IDs matched the atlas RNA matrix.")
}

metadata <- metadata[match(cells_use, metadata$cell), , drop = FALSE]
expr_sub <- expr_mat[GENES_TO_PLOT, cells_use, drop = FALSE]

plot_df <- lapply(
  GENES_TO_PLOT,
  function(g) {
    data.frame(
      cell = cells_use,
      gene = g,
      expression = as.numeric(expr_sub[g, ]),
      stringsAsFactors = FALSE
    )
  }
) %>%
  bind_rows() %>%
  left_join(
    metadata %>%
      select(cell, timepoint_standard, group_display),
    by = "cell"
  ) %>%
  filter(!is.na(group_display), !is.na(timepoint_standard)) %>%
  mutate(
    gene = factor(gene, levels = GENES_TO_PLOT),
    group_display = factor(
      group_display,
      levels = unname(GROUP_DISPLAY[GROUP_ORDER_RAW])
    ),
    timepoint_standard = factor(
      timepoint_standard,
      levels = TIMEPOINTS_TO_PLOT
    )
  )

write.csv(
  plot_df,
  file.path(OUTDIR, "Figure4D_v7_plot_data_used.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 4. Plot
# ------------------------------------------------------------
# Use FULL violins for both groups, placed side-by-side.
dodge_width <- 0.72
violin_width <- 0.60

p <- ggplot(
  plot_df,
  aes(
    x = group_display,
    y = expression,
    fill = timepoint_standard
  )
) +
  geom_violin(
    position = position_dodge(width = dodge_width),
    width = violin_width,
    scale = "width",
    trim = TRUE,
    colour = "#333333",
    linewidth = 0.32,
    alpha = 0.98
  ) +
  facet_grid(
    rows = vars(gene),
    scales = "free_y",
    switch = "y"
  ) +
  scale_fill_manual(
    values = TIMEPOINT_COLORS,
    breaks = TIMEPOINTS_TO_PLOT,
    labels = c("1 dpi", "7 dpi"),
    name = NULL,
    drop = FALSE
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.02, 0.06)),
    n.breaks = 2,
    labels = scales::label_number(accuracy = 0.1)
  ) +
  geom_hline(
    yintercept = 0,
    colour = "#333333",
    linewidth = 0.34
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_classic(base_family = FONT) +
  theme(
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text.y.left = element_text(
      angle = 0,
      size = 11.5,
      margin = margin(r = 8),
      colour = "black"
    ),

    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 9.3,
      colour = "black"
    ),
    axis.text.y = element_text(
      size = 8.0,
      colour = "black"
    ),

    axis.ticks = element_line(
      colour = "#333333",
      linewidth = 0.32
    ),
    axis.line.x = element_line(
      colour = "#333333",
      linewidth = 0.42
    ),
    axis.line.y = element_line(
      colour = "#333333",
      linewidth = 0.42
    ),

    panel.spacing.y = unit(0.06, "lines"),

    # Legend on the right side of the whole figure
    legend.position = "right",
    legend.justification = c("center", "center"),
    legend.text = element_text(size = 9.4),
    legend.key.height = unit(0.46, "cm"),
    legend.key.width = unit(0.30, "cm"),
    legend.background = element_rect(fill = "white", colour = NA),

    plot.margin = margin(
      t = 5,
      r = 12,
      b = 5,
      l = 10
    )
  )

# ------------------------------------------------------------
# 5. Save helper
# ------------------------------------------------------------
draw_with_panel_letter <- function(plot_obj) {
  grid.newpage()
  grid.draw(ggplotGrob(plot_obj))

  grid.text(
    "D",
    x = unit(0.010, "npc"),
    y = unit(0.985, "npc"),
    just = c("left", "top"),
    gp = gpar(
      fontsize = 24,
      fontface = "bold",
      fontfamily = FONT
    )
  )
}

# ------------------------------------------------------------
# 6. Save outputs
# ------------------------------------------------------------
if (capabilities("cairo")) {
  grDevices::cairo_pdf(
    filename = PDF_FILE,
    width = 8.3,
    height = 4.25,
    family = FONT,
    onefile = FALSE,
    bg = "white"
  )
} else {
  grDevices::pdf(
    file = PDF_FILE,
    width = 8.3,
    height = 4.25,
    family = "Helvetica",
    onefile = FALSE,
    useDingbats = FALSE,
    bg = "white"
  )
}

draw_with_panel_letter(p)
dev.off()

grDevices::png(
  filename = PNG_FILE,
  width = 2656,
  height = 1360,
  res = 320,
  bg = "white"
)

draw_with_panel_letter(p)
dev.off()

# ------------------------------------------------------------
# 7. Audit
# ------------------------------------------------------------
write.csv(
  data.frame(
    atlas = ATLAS_RDS_FILE,
    manifest = STEP8A_MANIFEST_FILE,
    genes = paste(GENES_TO_PLOT, collapse = " | "),
    timepoints = paste(TIMEPOINTS_TO_PLOT, collapse = " vs "),
    n_cells = length(cells_use),
    pdf_exists = file.exists(PDF_FILE),
    png_exists = file.exists(PNG_FILE),
    stringsAsFactors = FALSE
  ),
  file.path(OUTDIR, "Figure4D_v7_output_audit.csv"),
  row.names = FALSE
)

cat("\n============================================================\n")
cat("Figure 4D v7 GLIA readability revision completed.\n")
cat("Full side-by-side violin plots (NOT split halves).\n")
cat("Legend placed on the right side of the whole figure.\n")
cat("Left y-axis labels reduced to avoid overlap.\n\n")
cat("PDF:\n  ", PDF_FILE, "\n", sep = "")
cat("PNG:\n  ", PNG_FILE, "\n", sep = "")
cat("============================================================\n")
