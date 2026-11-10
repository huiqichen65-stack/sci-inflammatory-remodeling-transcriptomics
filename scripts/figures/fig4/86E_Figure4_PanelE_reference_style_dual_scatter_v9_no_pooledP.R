# ==============================================================================
# v9 GLIA revision: remove pooled cell-level P values from the displayed annotations;
# retain Spearman rho as a descriptive association measure; analysis/statistics unchanged.
# Figure 4E | Reference-style dual scatter plots
#
# PURPOSE
# -------
# Generate Figure 4E as two publication-style correlation panels that connect
# the prioritized TNF axis from upstream regulatory state to downstream
# endothelial response.
#
#   LEFT  : recipient-side response coupling
#           TNF_NFKB_RESPONSE score  vs  ENDOTHELIAL_ACTIVATION score
#           in Tnfrsf1a-positive Endothelial cells at 7 dpi
#
#   RIGHT : sender-side regulatory-to-output coupling
#           frozen Step7D Nfe2l2(+) AUCell regulon activity  vs  Tnf expression
#           in Target program (Microglia_Arm_2) cells at 7 dpi
#
# IMPORTANT DISPLAY RULE FOR THE RIGHT PANEL
# ------------------------------------------
# - Correlation statistics are taken from the frozen Step8D closure table and
#   therefore reflect the COMPLETE target-program cell set.
# - The black regression line and grey confidence band are ALSO fitted using the
#   complete target-program cell set.
# - For visual clarity only, cells with undetected Tnf expression (Tnf == 0)
#   are omitted from the displayed point layer.
# - The right-panel viewport is zoomed to the Tnf-positive range, but this does
#   not change the all-cell fit or reported statistics.
#
# FROZEN INPUTS
# -------------
# 1) Step1F atlas Seurat object: normalized RNA expression
# 2) Step8A frozen manifest: standardized cell metadata / cell groups
# 3) Step8C recipient-program gene-coverage table: signature genes
# 4) Step7D AUCell RDS: real Nfe2l2(+) regulon activity by cell
# 5) Step8D closure table: frozen Nfe2l2-Tnf association statistics
#
# No CellChat analysis is rerun.
# ==============================================================================

rm(list = ls())
gc()
options(stringsAsFactors = FALSE, scipen = 999)

# ------------------------------------------------------------------------------
# 0. SETTINGS
# ------------------------------------------------------------------------------
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

STEP8C_GENESET_FILE <- file.path(
  ROOT,
  "04_results",
  "step8C_GSE162610_recipient_program_integration",
  "tables",
  "02_recipient_program_gene_coverage.tsv"
)

STEP7D_AUC_FILE <- file.path(
  ROOT,
  "04_results",
  "step7D_GSE162610_AUCell_regulon_activity_scoring",
  "RDS",
  "GSE162610_step7D_AUC_regulons_by_cells.rds"
)

STEP8D_CLOSURE_FILE <- file.path(
  ROOT,
  "04_results",
  "step8D_GSE162610_mechanism_closure_QC",
  "tables",
  "06_focus_TF_ligand_closure_summary.tsv"
)

OUTDIR <- file.path(
  ROOT,
  "09_final_figures",
  "GLIA",
  "Figure4",
  "single_panels",
  "Figure4E"
)

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

PDF_FILE <- file.path(
  OUTDIR,
  "Figure4E_reference_style_dual_scatter_v9_no_pooledP.pdf"
)
PNG_FILE <- file.path(
  OUTDIR,
  "Figure4E_reference_style_dual_scatter_v9_no_pooledP.png"
)
LEFT_DATA_FILE <- file.path(
  OUTDIR,
  "Figure4E_v6_leftpanel_TNFNFKB_vs_EndothelialActivation.csv"
)
RIGHT_DATA_FILE_FULL <- file.path(
  OUTDIR,
  "Figure4E_v6_rightpanel_Nfe2l2AUCell_vs_Tnf_allcells.csv"
)
RIGHT_DATA_FILE_DISPLAY <- file.path(
  OUTDIR,
  "Figure4E_v6_rightpanel_Nfe2l2AUCell_vs_Tnf_display_TnfPositiveOnly.csv"
)
AUDIT_FILE <- file.path(
  OUTDIR,
  "Figure4E_reference_style_dual_scatter_v9_no_pooledP_audit.csv"
)

FONT_FAMILY <- "Arial"
TIMEPOINT_USE <- "7dpi"
RECIPIENT_GROUP <- "Endothelial"
SENDER_GROUP <- "Microglia_Arm_2"
RECIPIENT_RECEPTOR_GENE <- "Tnfrsf1a"
LEFT_PROGRAM_X <- "TNF_NFKB_RESPONSE"
LEFT_PROGRAM_Y <- "ENDOTHELIAL_ACTIVATION"
RIGHT_TF <- "Nfe2l2"
RIGHT_REGULON <- "Nfe2l2(+)"
RIGHT_LIGAND_GENE <- "Tnf"
MAX_POINTS_PER_PANEL <- 1600
POINT_SIZE <- 1.60
POINT_ALPHA <- 0.70

if (.Platform$OS.type == "windows") {
  try(grDevices::windowsFonts(Arial = grDevices::windowsFont("Arial")), silent = TRUE)
}

# ------------------------------------------------------------------------------
# 1. PACKAGES
# ------------------------------------------------------------------------------
required_pkgs <- c("Seurat", "ggplot2", "dplyr", "data.table", "grid", "gridExtra")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(
    "Missing package(s): ", paste(missing_pkgs, collapse = ", "),
    "\nInstall with:\ninstall.packages(c(",
    paste(sprintf('"%s"', missing_pkgs), collapse = ", "),
    "))",
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(data.table)
  library(grid)
  library(gridExtra)
})

# ------------------------------------------------------------------------------
# 2. HELPERS
# ------------------------------------------------------------------------------
find_exact_file <- function(root, filename) {
  if (!dir.exists(root)) return(NA_character_)
  hits <- list.files(
    root,
    pattern = paste0("^", gsub("\\.", "\\\\.", filename), "$"),
    recursive = TRUE,
    full.names = TRUE
  )
  if (length(hits) == 0L) return(NA_character_)
  hits[1]
}

resolve_file <- function(primary_path, fallback_name) {
  if (file.exists(primary_path)) return(primary_path)
  hit <- find_exact_file(file.path(ROOT, "04_results"), fallback_name)
  if (!is.na(hit) && file.exists(hit)) return(hit)
  stop("Could not find required file:\n", primary_path, call. = FALSE)
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

format_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) return(format(p, scientific = TRUE, digits = 2))
  sprintf("%.4f", p)
}

cor_label <- function(x, y, method = "spearman") {
  ct <- suppressWarnings(stats::cor.test(x, y, method = method, exact = FALSE))
  paste0("Spearman rho = ", sprintf("%.3f", unname(ct$estimate)))
}

signature_score_from_expr <- function(expr_sub) {
  zmat <- t(scale(t(as.matrix(expr_sub))))
  zmat[!is.finite(zmat)] <- 0
  colMeans(zmat, na.rm = TRUE)
}

sample_for_plot <- function(df, max_n) {
  if (nrow(df) <= max_n) return(df)
  per_sample <- ceiling(max_n / length(unique(df$sample_id)))
  out <- df %>%
    group_by(sample_id) %>%
    slice_sample(n = min(n(), per_sample)) %>%
    ungroup()
  if (nrow(out) > max_n) out <- out %>% slice_sample(n = max_n)
  out
}

normalize_name <- function(x) {
  x <- as.character(x)
  x <- gsub("\\s+", "", x)
  x <- gsub("_", "", x)
  x <- gsub("-", "", x)
  tolower(x)
}

extract_auc_matrix <- function(obj, cell_ids) {
  # Accept matrix/data.frame directly
  candidates <- list(obj)

  # Expand one level if list-like
  if (is.list(obj)) {
    candidates <- c(candidates, obj)
    for (nm in names(obj)) {
      if (is.list(obj[[nm]])) {
        candidates <- c(candidates, obj[[nm]])
      }
    }
  }

  best <- NULL
  best_overlap <- -1

  for (cand in candidates) {
    ok <- is.matrix(cand) || is.data.frame(cand) || inherits(cand, "dgCMatrix") || inherits(cand, "Matrix")
    if (!ok) next

    rn <- rownames(cand)
    cn <- colnames(cand)
    if (is.null(rn) || is.null(cn)) next

    row_overlap <- sum(rn %in% cell_ids)
    col_overlap <- sum(cn %in% cell_ids)
    overlap <- max(row_overlap, col_overlap)
    if (overlap > best_overlap) {
      best <- cand
      best_overlap <- overlap
    }
  }

  if (is.null(best)) {
    stop("Could not identify a usable AUCell matrix in the Step7D RDS.", call. = FALSE)
  }

  rn <- rownames(best)
  cn <- colnames(best)
  row_overlap <- sum(rn %in% cell_ids)
  col_overlap <- sum(cn %in% cell_ids)

  if (col_overlap >= row_overlap) {
    # rows = regulons, cols = cells
    out <- best
  } else {
    # rows = cells, cols = regulons
    out <- t(as.matrix(best))
  }

  out
}

# ------------------------------------------------------------------------------
# 3. LOAD INPUTS
# ------------------------------------------------------------------------------
ATLAS_RDS_FILE <- resolve_file(ATLAS_RDS_FILE, "GSE162610_broad_annotation_audit.rds")
STEP8A_MANIFEST_FILE <- resolve_file(STEP8A_MANIFEST_FILE, "GSE162610_step8A_CellChat_input_manifest.rds")
STEP8C_GENESET_FILE <- resolve_file(STEP8C_GENESET_FILE, "02_recipient_program_gene_coverage.tsv")
STEP7D_AUC_FILE <- resolve_file(STEP7D_AUC_FILE, "GSE162610_step7D_AUC_regulons_by_cells.rds")
STEP8D_CLOSURE_FILE <- resolve_file(STEP8D_CLOSURE_FILE, "06_focus_TF_ligand_closure_summary.tsv")

atlas <- readRDS(ATLAS_RDS_FILE)
if (!inherits(atlas, "Seurat")) stop("Atlas RDS is not a Seurat object.", call. = FALSE)

manifest <- readRDS(STEP8A_MANIFEST_FILE)
if (!is.list(manifest) || is.null(manifest$cell_metadata)) {
  stop("Step8A manifest$cell_metadata is missing.", call. = FALSE)
}

metadata <- as.data.frame(manifest$cell_metadata, stringsAsFactors = FALSE)
required_meta <- c("cell", "sample_id", "timepoint_standard", "cellchat_group", "include_step8_timecourse")
missing_meta <- setdiff(required_meta, colnames(metadata))
if (length(missing_meta) > 0) {
  stop("Step8A cell_metadata missing: ", paste(missing_meta, collapse = ", "), call. = FALSE)
}

metadata <- metadata %>%
  mutate(
    cell = as.character(cell),
    sample_id = as.character(sample_id),
    timepoint_standard = as.character(timepoint_standard),
    cellchat_group = as.character(cellchat_group),
    include_step8_timecourse = as_logical_safe(include_step8_timecourse)
  ) %>%
  filter(include_step8_timecourse %in% TRUE, timepoint_standard == TIMEPOINT_USE)

if (nrow(metadata) == 0) stop("No cells remained after metadata filtering.", call. = FALSE)

# ------------------------------------------------------------------------------
# 4. RNA EXPRESSION
# ------------------------------------------------------------------------------
DefaultAssay(atlas) <- "RNA"
expr_mat <- try(GetAssayData(atlas, assay = "RNA", layer = "data"), silent = TRUE)
if (inherits(expr_mat, "try-error")) {
  expr_mat <- GetAssayData(atlas, assay = "RNA", slot = "data")
}

required_genes <- c(RECIPIENT_RECEPTOR_GENE, RIGHT_LIGAND_GENE)
missing_genes <- setdiff(required_genes, rownames(expr_mat))
if (length(missing_genes) > 0) {
  stop("Missing required gene(s) in RNA assay: ", paste(missing_genes, collapse = ", "), call. = FALSE)
}

# ------------------------------------------------------------------------------
# 5. LEFT PANEL: program-response coupling in endothelial cells
# ------------------------------------------------------------------------------
geneset_tbl <- data.table::fread(STEP8C_GENESET_FILE, data.table = FALSE, showProgress = FALSE)
need_gene_cols <- c("response_program", "present_genes")
miss_gene_cols <- setdiff(need_gene_cols, colnames(geneset_tbl))
if (length(miss_gene_cols) > 0) {
  stop("Gene-set table missing columns: ", paste(miss_gene_cols, collapse = ", "), call. = FALSE)
}

get_program_genes <- function(program_name) {
  row <- geneset_tbl %>% filter(response_program == program_name)
  if (nrow(row) != 1) {
    stop("Could not uniquely identify response program: ", program_name, call. = FALSE)
  }
  genes <- unlist(strsplit(as.character(row$present_genes[1]), "\\s*\\|\\s*"))
  genes <- unique(trimws(genes))
  genes <- genes[nchar(genes) > 0]
  genes <- intersect(genes, rownames(expr_mat))
  if (length(genes) == 0) {
    stop("No genes for program ", program_name, " were found in the RNA matrix.", call. = FALSE)
  }
  genes
}

left_x_genes <- get_program_genes(LEFT_PROGRAM_X)
left_y_genes <- get_program_genes(LEFT_PROGRAM_Y)

left_meta <- metadata %>%
  filter(cellchat_group == RECIPIENT_GROUP, cell %in% colnames(expr_mat))
if (nrow(left_meta) == 0) stop("No endothelial cells found for left panel.", call. = FALSE)

left_cells_all <- left_meta$cell
left_expr_all <- expr_mat[, left_cells_all, drop = FALSE]
left_tnfrsf1a <- as.numeric(left_expr_all[RECIPIENT_RECEPTOR_GENE, left_cells_all])
left_score_x <- signature_score_from_expr(left_expr_all[left_x_genes, , drop = FALSE])
left_score_y <- signature_score_from_expr(left_expr_all[left_y_genes, , drop = FALSE])

left_df <- data.frame(
  cell = left_cells_all,
  sample_id = left_meta$sample_id,
  timepoint_standard = left_meta$timepoint_standard,
  Tnfrsf1a = left_tnfrsf1a,
  TNF_NFKB_RESPONSE = as.numeric(left_score_x),
  ENDOTHELIAL_ACTIVATION = as.numeric(left_score_y),
  stringsAsFactors = FALSE
) %>%
  filter(Tnfrsf1a > 0)

if (nrow(left_df) < 10) {
  stop("Too few Tnfrsf1a-positive endothelial cells for the left panel.", call. = FALSE)
}

left_label <- cor_label(left_df$TNF_NFKB_RESPONSE, left_df$ENDOTHELIAL_ACTIVATION, method = "spearman")
left_label_sub <- "Descriptive cell-level association"
write.csv(left_df, LEFT_DATA_FILE, row.names = FALSE)
left_plot_df <- sample_for_plot(left_df, MAX_POINTS_PER_PANEL)

# ------------------------------------------------------------------------------
# 6. RIGHT PANEL: real Step7D AUCell activity vs Tnf
# ------------------------------------------------------------------------------
right_meta <- metadata %>%
  filter(cellchat_group == SENDER_GROUP, cell %in% colnames(expr_mat))
if (nrow(right_meta) == 0) stop("No target-program cells found for right panel.", call. = FALSE)

right_cells <- right_meta$cell
right_expr <- expr_mat[, right_cells, drop = FALSE]

auc_obj <- readRDS(STEP7D_AUC_FILE)
auc_mat <- extract_auc_matrix(auc_obj, cell_ids = metadata$cell)

# find the regulon row robustly
row_candidates <- rownames(auc_mat)
row_norm <- normalize_name(row_candidates)
target_norm <- normalize_name(c(RIGHT_REGULON, RIGHT_TF, paste0(RIGHT_TF, "(+)") ))
idx <- which(row_norm %in% target_norm)
if (length(idx) == 0) {
  # fallback: partial match on nfe2l2 and plus sign or exact TF
  idx <- grep("nfe2l2", row_norm)
}
if (length(idx) == 0) {
  stop("Could not locate an Nfe2l2 regulon in the Step7D AUCell matrix.", call. = FALSE)
}
if (length(idx) > 1) {
  # prefer an explicit positive regulon if available
  plus_idx <- idx[grepl("\\+", row_candidates[idx], fixed = FALSE) | grepl("positive", row_norm[idx])]
  if (length(plus_idx) >= 1) idx <- plus_idx[1] else idx <- idx[1]
}
regulon_name_used <- row_candidates[idx[1]]

right_auc_values <- as.numeric(auc_mat[idx[1], right_cells])
right_tnf_values <- as.numeric(right_expr[RIGHT_LIGAND_GENE, right_cells])

right_df_full <- data.frame(
  cell = right_cells,
  sample_id = right_meta$sample_id,
  timepoint_standard = right_meta$timepoint_standard,
  Nfe2l2_AUCell = right_auc_values,
  Tnf = right_tnf_values,
  stringsAsFactors = FALSE
)

if (nrow(right_df_full) < 10) {
  stop("Too few target-program cells for the right panel.", call. = FALSE)
}

# Display only Tnf-positive cells; statistics remain frozen/all-cells
right_df_display <- right_df_full %>% filter(Tnf > 0)
if (nrow(right_df_display) < 10) {
  stop("Too few Tnf-positive target-program cells for displayed right-panel scatter.", call. = FALSE)
}

write.csv(right_df_full, RIGHT_DATA_FILE_FULL, row.names = FALSE)
write.csv(right_df_display, RIGHT_DATA_FILE_DISPLAY, row.names = FALSE)
right_plot_df <- sample_for_plot(right_df_display, MAX_POINTS_PER_PANEL)

# Frozen Step8D stats for annotation
closure_tbl <- data.table::fread(STEP8D_CLOSURE_FILE, data.table = FALSE, showProgress = FALSE)
need_closure_cols <- c("TF", "regulon", "ligand", "pooled_spearman_rho", "pooled_spearman_p_value", "sample_rhos")
miss_closure_cols <- setdiff(need_closure_cols, colnames(closure_tbl))
if (length(miss_closure_cols) > 0) {
  stop("Closure table missing columns: ", paste(miss_closure_cols, collapse = ", "), call. = FALSE)
}

closure_row <- closure_tbl %>%
  filter(TF == RIGHT_TF, ligand == RIGHT_LIGAND_GENE)
if (nrow(closure_row) == 0) {
  stop("Could not find the frozen Nfe2l2-Tnf closure row in Step8D.", call. = FALSE)
}
closure_row <- closure_row[1, , drop = FALSE]

right_label_main <- paste0(
  "Spearman rho = ", sprintf("%.3f", as.numeric(closure_row$pooled_spearman_rho))
)
right_label_sub <- paste0(
  "Descriptive cell-level association",
  "\nTnf = 0 cells omitted from display; fit/statistic use all cells"
)

# ------------------------------------------------------------------------------
# 7. COLORS
# ------------------------------------------------------------------------------
all_samples <- unique(c(left_df$sample_id, right_df_full$sample_id))
palette_base <- c(
  "#6A3D9A", "#1F78B4", "#E66101", "#FDB863",
  "#B2ABD2", "#F46D43", "#3288BD", "#9E0142",
  "#FEE08B", "#66C2A5", "#5E3C99", "#F1A340"
)
if (length(all_samples) > length(palette_base)) {
  palette_base <- grDevices::colorRampPalette(palette_base)(length(all_samples))
}
sample_colors <- setNames(palette_base[seq_along(all_samples)], all_samples)

# ------------------------------------------------------------------------------
# 8. PLOT FUNCTIONS
# ------------------------------------------------------------------------------
base_theme <- theme_classic(base_family = FONT_FAMILY) +
  theme(
    plot.title = element_text(size = 12.2, face = "plain", hjust = 0.5, margin = margin(b = 4)),
    axis.title = element_text(size = 10.6, colour = "black"),
    axis.text = element_text(size = 9.1, colour = "black"),
    axis.line = element_line(linewidth = 0.45, colour = "#333333"),
    axis.ticks = element_line(linewidth = 0.35, colour = "#333333"),
    plot.margin = margin(14, 8, 6, 8),
    legend.position = "none"
  )

make_scatter <- function(df, xvar, yvar, xlab, ylab, title_text,
                         ann_text_main, ann_text_sub = NULL,
                         smooth_df = NULL, zoom_y_to_display = FALSE) {
  if (is.null(smooth_df)) smooth_df <- df

  xvals <- df[[xvar]]
  yvals <- df[[yvar]]
  x_rng <- range(xvals, na.rm = TRUE)
  y_rng <- range(yvals, na.rm = TRUE)
  ann_x <- x_rng[1] + 0.06 * diff(x_rng)
  ann_y <- y_rng[2] - 0.06 * diff(y_rng)

  p <- ggplot(df, aes_string(x = xvar, y = yvar)) +
    geom_smooth(
      data = smooth_df,
      mapping = aes_string(x = xvar, y = yvar),
      method = "lm",
      se = TRUE,
      inherit.aes = FALSE,
      linewidth = 0.70,
      colour = "black",
      fill = "#BFBFBF",
      alpha = 0.35
    ) +
    geom_point(
      aes(colour = sample_id),
      size = POINT_SIZE,
      alpha = POINT_ALPHA,
      show.legend = FALSE
    ) +
    scale_colour_manual(values = sample_colors) +
    scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0.02, 0.14))
    ) +
    annotate(
      "text",
      x = ann_x,
      y = ann_y,
      label = ann_text_main,
      hjust = 0,
      vjust = 1,
      size = 3.9,
      family = FONT_FAMILY,
      fontface = "italic"
    ) +
    labs(title = title_text, x = xlab, y = ylab) +
    base_theme

  if (!is.null(ann_text_sub)) {
    p <- p + annotate(
      "text",
      x = ann_x,
      y = ann_y - 0.08 * diff(y_rng),
      label = ann_text_sub,
      hjust = 0,
      vjust = 1,
      size = 3.3,
      family = FONT_FAMILY,
      colour = "#4A4A4A"
    )
  }

  # For the right panel, the regression is fit on ALL cells, but the viewport
  # is zoomed to the displayed Tnf-positive range so hidden zero-expression
  # cells do not force the y-axis back to zero. This changes display only, not
  # the fitted line or reported statistics.
  if (isTRUE(zoom_y_to_display)) {
    yr <- range(yvals, na.rm = TRUE)
    ypad <- max(0.05, 0.04 * diff(yr))
    p <- p + coord_cartesian(ylim = c(yr[1] - ypad, yr[2] + ypad))
  }

  p
}

p_left <- make_scatter(
  df = left_plot_df,
  xvar = "TNF_NFKB_RESPONSE",
  yvar = "ENDOTHELIAL_ACTIVATION",
  xlab = "Signature score of TNF/NFκB response",
  ylab = "Signature score of endothelial activation",
  title_text = "GSE162610",
  ann_text_main = left_label,
  ann_text_sub = left_label_sub
)

p_right <- make_scatter(
  df = right_plot_df,
  xvar = "Nfe2l2_AUCell",
  yvar = "Tnf",
  xlab = "Activity of Nfe2l2(+) regulon",
  ylab = "Expression of Tnf",
  title_text = "GSE162610",
  ann_text_main = right_label_main,
  ann_text_sub = right_label_sub,
  smooth_df = right_df_full,
  zoom_y_to_display = TRUE
)

# ------------------------------------------------------------------------------
# 9. ARRANGE AND SAVE
# ------------------------------------------------------------------------------
page_grob <- function() {
  gridExtra::arrangeGrob(p_left, p_right, ncol = 2, widths = c(1, 1))
}

draw_page <- function() {
  g <- page_grob()
  grid.newpage()
  grid.draw(g)
  grid.text(
    "E",
    x = unit(0.012, "npc"),
    y = unit(0.985, "npc"),
    just = c("left", "top"),
    gp = gpar(fontsize = 30, fontface = "bold", fontfamily = FONT_FAMILY)
  )
}

if (file.exists(PDF_FILE)) unlink(PDF_FILE)
if (file.exists(PNG_FILE)) unlink(PNG_FILE)

if (capabilities("cairo")) {
  grDevices::cairo_pdf(filename = PDF_FILE, width = 7.4, height = 3.9, family = FONT_FAMILY, onefile = FALSE, bg = "white")
} else {
  grDevices::pdf(file = PDF_FILE, width = 7.4, height = 3.9, family = "Helvetica", onefile = FALSE, useDingbats = FALSE, bg = "white")
}
draw_page()
dev.off()

grDevices::png(filename = PNG_FILE, width = 2368, height = 1248, res = 320, bg = "white")
draw_page()
dev.off()

# ------------------------------------------------------------------------------
# 10. AUDIT
# ------------------------------------------------------------------------------
audit <- data.frame(
  item = c(
    "atlas_rds",
    "step8a_manifest",
    "step8c_geneset_table",
    "step7d_auc_file",
    "step8d_closure_file",
    "timepoint_use",
    "left_panel_group",
    "left_panel_filter",
    "left_panel_x_program",
    "left_panel_y_program",
    "right_panel_group",
    "right_panel_auc_regulon_requested",
    "right_panel_auc_regulon_used",
    "right_panel_stats_source",
    "right_panel_display_filter",
    "displayed_correlation_annotation",
    "left_n_cells_full",
    "right_n_cells_all",
    "right_n_cells_display",
    "left_n_points_plotted",
    "right_n_points_plotted",
    "pdf",
    "png"
  ),
  value = c(
    ATLAS_RDS_FILE,
    STEP8A_MANIFEST_FILE,
    STEP8C_GENESET_FILE,
    STEP7D_AUC_FILE,
    STEP8D_CLOSURE_FILE,
    TIMEPOINT_USE,
    RECIPIENT_GROUP,
    paste0(RECIPIENT_RECEPTOR_GENE, " > 0"),
    LEFT_PROGRAM_X,
    LEFT_PROGRAM_Y,
    SENDER_GROUP,
    RIGHT_REGULON,
    regulon_name_used,
    "Frozen Step8D pooled all-cell statistics",
    paste0(RIGHT_LIGAND_GENE, " > 0 (display only)"),
    "Spearman rho only; pooled cell-level P values intentionally omitted from figure",
    as.character(nrow(left_df)),
    as.character(nrow(right_df_full)),
    as.character(nrow(right_df_display)),
    as.character(nrow(left_plot_df)),
    as.character(nrow(right_plot_df)),
    PDF_FILE,
    PNG_FILE
  ),
  stringsAsFactors = FALSE
)
write.csv(audit, AUDIT_FILE, row.names = FALSE)

cat("\n============================================================\n")
cat("Figure 4E v9 GLIA revision completed successfully.\n")
cat("Left : TNF_NFKB_RESPONSE vs ENDOTHELIAL_ACTIVATION\n")
cat("       (Tnfrsf1a-positive endothelial cells, 7 dpi)\n")
cat("Right: real Step7D Nfe2l2(+) AUCell activity vs Tnf\n")
cat("       (all-cell regression/statistics; displayed points: Tnf > 0)\n\n")
cat("PDF:\n  ", PDF_FILE, "\n", sep = "")
cat("PNG:\n  ", PNG_FILE, "\n", sep = "")
cat("============================================================\n")
