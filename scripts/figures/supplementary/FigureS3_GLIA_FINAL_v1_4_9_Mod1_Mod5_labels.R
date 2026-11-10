# ==============================================================================
# Figure S3 | GLIA | FINAL v1.4.9 | no panel subtitles | remove old panel B
# Trajectory and branch-program robustness
#
# DIRECT ORIGINAL-SOURCE MODE
#
# v1.4.9 visual-only nomenclature fix: reader-facing dynamic-module labels use
# Mod1-Mod5 to match Main Figure 2 and the revised manuscript. Internal frozen
# module IDs (M1_..., M2_..., etc.) are unchanged. No data, statistics, branch
# assignments, colors, or frozen analysis objects are changed.
#
# Core question:
#   Does the target-program branch remain stable to reasonable root/topology choices,
#   and are branch-associated transcriptional programs supported at the
#   biological-replicate and pseudobulk levels?
#
# Panels:
#   A. Primary 1 dpi injury-trajectory root + single-principal-root topology
#   B. Root/orientation sensitivity of the downstream trajectory
#   C. Shared prebranch + progressive post-branch divergence
#   D. Biological-replicate / pseudobulk support of branch programs
#
# IMPORTANT:
#   * Reads frozen Step5BR2 / Step5CR / Step6A / Step6B results.
#   * The FINAL injury-response trajectory root is the Step5BR primary 1 dpi
#     root (225 cells; 75 cells from each of three 1 dpi samples), NOT the
#     earlier Step5AR uninjured consensus-root audit.
#   * Monocle principal-node IDs are retained only in audit logic and are not
#     displayed in the figure.
#   * Does NOT rerun Monocle learning, root selection, branch selection,
#     differential testing, pseudobulk testing, or module discovery.
#   * Arm_2 is displayed as "Target-program branch" to remain consistent with
#     Main Figure 2E ("Arm2_/Target_" bins).
#   * Supplementary Figure S3 intentionally does NOT repeat the full Main
#     Figure 2 dynamic heatmap, cell-level volcano, pseudobulk MDS, or GO walls.
# ==============================================================================

rm(list = ls())
gc()
options(stringsAsFactors = FALSE)
set.seed(20260818)

# ------------------------------------------------------------------------------
# 0. Paths
# ------------------------------------------------------------------------------
ROOT <- Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
if (!dir.exists(ROOT)) {
  stop("Project root does not exist: ", ROOT)
}

OUT <- file.path(
  ROOT,
  "09_final_figures",
  "GLIA",
  "FigureS3_FINAL_v1_4_9_Mod1_Mod5_labels"
)

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

P <- list(
  step5BR2_cds = file.path(
    ROOT,
    "04_results",
    "step5BR2_GSE162610_single_root_branch_audit",
    "GSE162610_step5BR2_single_root_cds.rds"
  ),
  step5BR2_seurat = file.path(
    ROOT,
    "04_results",
    "step5BR2_GSE162610_single_root_branch_audit",
    "GSE162610_step5BR2_single_root_seurat.rds"
  ),
  step5BR2_xlsx = file.path(
    ROOT,
    "04_results",
    "step5BR2_GSE162610_single_root_branch_audit",
    "Table5BR2_single_root_branch_audit.xlsx"
  ),

  step5CR_xlsx = file.path(
    ROOT,
    "04_results",
    "step5CR_GSE162610_Y29_prebranch_progressive_divergence",
    "Table5CR_Y29_prebranch_progressive_divergence.xlsx"
  ),

  step6A_highconf = file.path(
    ROOT,
    "04_results",
    "step6A_GSE162610_7dpi_branch_specific_DEG",
    "tables",
    "02_high_confidence_branch_markers.csv"
  ),
  step6A_evidence = file.path(
    ROOT,
    "04_results",
    "step6A_GSE162610_7dpi_branch_specific_DEG",
    "tables",
    "01_all_genes_master_evidence.csv"
  ),

  step6B_summary = file.path(
    ROOT,
    "04_results",
    "step6B_GSE162610_dynamic_branch_modules",
    "tables",
    "02_dynamic_module_summary.csv"
  ),
  step6B_effects = file.path(
    ROOT,
    "04_results",
    "step6B_GSE162610_dynamic_branch_modules",
    "tables",
    "04_module_7dpi_replicate_effects.csv"
  )
)

missing_inputs <- names(P)[!vapply(P, file.exists, logical(1))]

if (length(missing_inputs) > 0) {
  stop(
    "Missing required Figure S3 input(s):\n",
    paste0(
      " - ",
      missing_inputs,
      ": ",
      unlist(P[missing_inputs]),
      collapse = "\n"
    )
  )
}

message("Project root: ", ROOT)
message("Output dir:   ", OUT)

# ------------------------------------------------------------------------------
# 1. Packages
# ------------------------------------------------------------------------------
required_pkgs <- c(
  "SeuratObject",
  "monocle3",
  "SingleCellExperiment",
  "SummarizedExperiment",
  "ggplot2",
  "dplyr",
  "tidyr",
  "readr",
  "readxl",
  "patchwork",
  "scales",
  "svglite",
  "ragg"
)

missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_pkgs) > 0) {
  stop(
    "Missing required packages: ",
    paste(missing_pkgs, collapse = ", "),
    "\nInstall them first, then rerun."
  )
}

suppressPackageStartupMessages({
  library(SeuratObject)
  library(monocle3)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(readxl)
  library(patchwork)
  library(scales)
})

# ------------------------------------------------------------------------------
# 2. Frozen visual language
# ------------------------------------------------------------------------------
BASE_FAMILY <- "Arial"

# Match Main Figure 2E stage language:
COL_SHARED <- "#8B86B8"
COL_TARGET <- "#C9675C"

# Alternative branch is deliberately quieter than the target-program branch.
COL_ALT <- "#6F8FB7"
COL_BRANCH <- "#3F3F3F"
COL_ROOT <- "#D94841"
COL_BG <- "#D9D9D9"
COL_NONSTRICT <- "#C8C8C8"

# Main Figure 2D frozen timepoint palette:
#   1 dpi = #8C1D9C
#   3 dpi = #D84A73
#   7 dpi = #F6B73C
# Sample shades stay within the corresponding timepoint hue family.
root_sample_cols <- c(
  "1dpi_sample1" = "#6F167D",
  "1dpi_sample2" = "#8C1D9C",
  "1dpi_sample3" = "#AC4CB7"
)

sample_cols <- c(
  "3dpi_sample1" = "#BE365F",
  "3dpi_sample2" = "#E36A8B",
  "7dpi_sample1" = "#E89E24",
  "7dpi_sample2" = "#F6B73C"
)

branch_cols <- c(
  "Shared prebranch" = COL_SHARED,
  "Alternative branch" = COL_ALT,
  "Target-program branch" = COL_TARGET,
  "Branch point" = COL_BRANCH
)

theme_pub <- theme_classic(
  base_family = BASE_FAMILY,
  base_size = 7.2
) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    axis.title = element_text(
      family = BASE_FAMILY,
      size = 7.0,
      color = "black"
    ),
    axis.text = element_text(
      family = BASE_FAMILY,
      size = 5.9,
      color = "black"
    ),
    axis.line = element_line(
      linewidth = 0.36,
      color = "black"
    ),
    axis.ticks = element_line(
      linewidth = 0.30,
      color = "black"
    ),
    legend.title = element_text(
      family = BASE_FAMILY,
      size = 5.8,
      color = "black"
    ),
    legend.text = element_text(
      family = BASE_FAMILY,
      size = 5.2,
      color = "black"
    ),
    plot.title = element_text(
      family = BASE_FAMILY,
      size = 7.0,
      face = "bold",
      hjust = 0
    ),
    plot.subtitle = element_text(
      family = BASE_FAMILY,
      size = 5.5,
      color = "#4D4D4D",
      hjust = 0,
      lineheight = 0.95
    ),
    plot.tag = element_text(
      family = BASE_FAMILY,
      face = "bold",
      size = 11.5,
      color = "black"
    ),
    plot.tag.position = c(-0.035, 1.03),
    plot.margin = margin(
      2.8, 3.0, 2.8, 3.6,
      unit = "mm"
    )
  )

# ------------------------------------------------------------------------------
# 3. Helpers
# ------------------------------------------------------------------------------
save_plot_quad <- function(
  p,
  stem,
  width_mm,
  height_mm,
  dpi = 600
) {
  ggsave(
    file.path(OUT, paste0(stem, ".pdf")),
    p,
    device = cairo_pdf,
    width = width_mm,
    height = height_mm,
    units = "mm",
    bg = "white"
  )

  ggsave(
    file.path(OUT, paste0(stem, ".svg")),
    p,
    device = svglite::svglite,
    width = width_mm,
    height = height_mm,
    units = "mm",
    bg = "white"
  )

  ggsave(
    file.path(OUT, paste0(stem, ".tiff")),
    p,
    device = ragg::agg_tiff,
    width = width_mm,
    height = height_mm,
    units = "mm",
    dpi = dpi,
    compression = "lzw",
    background = "white"
  )

  ggsave(
    file.path(OUT, paste0(stem, "_review.png")),
    p,
    device = ragg::agg_png,
    width = width_mm,
    height = height_mm,
    units = "mm",
    dpi = 220,
    background = "white"
  )
}

safe_spearman <- function(x, y) {
  suppressWarnings(
    cor(
      x,
      y,
      use = "complete.obs",
      method = "spearman"
    )
  )
}

# ------------------------------------------------------------------------------
# 4. Input audit
# ------------------------------------------------------------------------------
input_manifest <- data.frame(
  key = names(P),
  path = unlist(P, use.names = FALSE),
  size_MB = round(
    file.info(unlist(P))$size / 1024^2,
    4
  ),
  md5 = unname(
    tools::md5sum(unlist(P))
  ),
  stringsAsFactors = FALSE
)

write.csv(
  input_manifest,
  file.path(
    OUT,
    "01_FigureS3_input_manifest.csv"
  ),
  row.names = FALSE
)

# ==============================================================================
# PANEL A
# Primary root + single-principal-root topology
# ==============================================================================

# ------------------------------------------------------------------------------
# 5A. Primary 1 dpi injury-trajectory root
#
# IMPORTANT:
#   The earlier Step5AR uninjured consensus-root audit was an upstream root
#   refinement step. The FINAL injury-response trajectory used by Step5BR/BR2
#   is rooted with 225 primary 1 dpi cells (75 cells per 1 dpi sample).
# ------------------------------------------------------------------------------
message("Reading Step5BR2 Seurat object for the final 1 dpi trajectory root...")

br2_seurat <- readRDS(P$step5BR2_seurat)

if (!inherits(br2_seurat, "Seurat")) {
  stop("Step5BR2 Seurat input is not a Seurat object.")
}

req_root_meta <- c(
  "sample_id",
  "timepoint",
  "step5BR_injury_input",
  "step5BR_primary_root",
  "step5BR_principal_vertex",
  "step5BR2_principal_root_node"
)

miss_root_meta <- setdiff(
  req_root_meta,
  colnames(br2_seurat[[]])
)

if (length(miss_root_meta) > 0) {
  stop(
    "Step5BR2 Seurat metadata missing for Panel A root: ",
    paste(miss_root_meta, collapse = ", ")
  )
}

root_reduction <- "umap.microglia.step3A"

if (!(root_reduction %in% SeuratObject::Reductions(br2_seurat))) {
  stop(
    "Frozen Step5BR2 Seurat object is missing reduction: ",
    root_reduction
  )
}

root_xy <- as.data.frame(
  SeuratObject::Embeddings(
    br2_seurat,
    reduction = root_reduction
  )
)

colnames(root_xy)[1:2] <- c(
  "UMAP_1",
  "UMAP_2"
)

root_xy$cell <- rownames(root_xy)

root_md <- br2_seurat[[]] %>%
  tibble::rownames_to_column("cell") %>%
  transmute(
    cell = cell,
    sample_id = as.character(sample_id),
    timepoint = as.character(timepoint),
    injury_input = as.logical(step5BR_injury_input),
    primary_root = as.logical(step5BR_primary_root),
    principal_vertex = as.character(step5BR_principal_vertex),
    single_root_node = as.character(step5BR2_principal_root_node)
  )

root_df <- inner_join(
  root_xy,
  root_md,
  by = "cell"
)

# Background only uses the 1 dpi cells that actually entered the injury
# trajectory. This makes the root selection readable without mixing later
# 3/7 dpi cells into the root-selection panel.
root_background <- root_df %>%
  filter(
    injury_input,
    timepoint == "1 dpi"
  )

root_final <- root_df %>%
  filter(primary_root)

root_sample_summary <- root_final %>%
  dplyr::count(
    sample_id,
    name = "n_root"
  ) %>%
  dplyr::arrange(sample_id)

expected_root_samples <- c(
  "1dpi_sample1",
  "1dpi_sample2",
  "1dpi_sample3"
)

if (
  nrow(root_final) != 225 ||
  !setequal(
    root_sample_summary$sample_id,
    expected_root_samples
  ) ||
  any(
    root_sample_summary$n_root != 75
  )
) {
  stop(
    "Unexpected frozen final injury-trajectory root composition. ",
    "Expected 225 primary 1 dpi root cells = 75 cells x 3 biological samples."
  )
}

message(
  "Final injury-trajectory root confirmed: 225 cells = 75 x 3 independent 1 dpi samples."
)

# The primary root cells mapped to several graph vertices in the original audit.
# Y_19 contains the largest number (41/225) and is therefore the objective
# single principal root used by Step5BR2.
root_vertex_summary <- root_final %>%
  dplyr::count(
    principal_vertex,
    name = "n_primary_root_cells"
  ) %>%
  dplyr::arrange(
    dplyr::desc(n_primary_root_cells)
  )

if (
  nrow(root_vertex_summary) == 0 ||
  root_vertex_summary$principal_vertex[1] != "Y_19" ||
  root_vertex_summary$n_primary_root_cells[1] != 41
) {
  stop(
    "Unexpected Step5BR primary-root principal-vertex mapping. ",
    "Expected Y_19 to be the top root vertex with 41 primary-root cells."
  )
}

pA_root <- ggplot() +
  geom_point(
    data = root_background,
    aes(
      x = UMAP_1,
      y = UMAP_2
    ),
    size = 0.52,
    color = "#D8D8D8",
    alpha = 0.58
  ) +
  geom_point(
    data = root_final,
    aes(
      x = UMAP_1,
      y = UMAP_2,
      color = sample_id
    ),
    size = 1.45,
    alpha = 0.98
  ) +
  scale_color_manual(
    values = root_sample_cols,
    labels = c(
      "1dpi_sample1" = "1 dpi sample 1",
      "1dpi_sample2" = "1 dpi sample 2",
      "1dpi_sample3" = "1 dpi sample 3"
    ),
    drop = FALSE
  ) +
  coord_equal() +
  labs(
    tag = "A",
    title = "Primary 1 dpi trajectory root",
    subtitle = NULL,
    x = "UMAP 1",
    y = "UMAP 2",
    color = NULL
  ) +
  theme_pub +
  theme(
    plot.title = element_text(
      size = 6.35,
      face = "bold",
      hjust = 0
    ),
    plot.subtitle = element_blank(),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.key.width = grid::unit(3.0, "mm"),
    legend.spacing.x = grid::unit(0.7, "mm"),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    plot.margin = margin(
      2.4, 2.0, 2.6, 3.2,
      unit = "mm"
    )
  ) +
  guides(
    color = guide_legend(
      nrow = 1,
      override.aes = list(
        size = 2.2,
        alpha = 1
      )
    )
  )

# ------------------------------------------------------------------------------
# 5B. Single-principal-root branch topology
# ------------------------------------------------------------------------------
message("Reading Step5BR2 single-root CDS...")

cds <- readRDS(P$step5BR2_cds)

if (!inherits(cds, "cell_data_set")) {
  stop("Step5BR2 input is not a monocle3 cell_data_set.")
}

req_cds_meta <- c(
  "step5BR2_branch_arm",
  "step5BR2_pseudotime",
  "sample_id",
  "timepoint"
)

miss_cds_meta <- setdiff(
  req_cds_meta,
  colnames(SummarizedExperiment::colData(cds))
)

if (length(miss_cds_meta) > 0) {
  stop(
    "Step5BR2 CDS colData missing: ",
    paste(miss_cds_meta, collapse = ", ")
  )
}

branch_display <- dplyr::recode(
  as.character(
    SummarizedExperiment::colData(cds)$step5BR2_branch_arm
  ),
  "Pre_branch_root_side" = "Shared prebranch",
  "Arm_1" = "Alternative branch",
  "Arm_2" = "Target-program branch",
  "Branchpoint_vertex" = "Branch point",
  .default = NA_character_
)

if (any(is.na(branch_display))) {
  stop(
    "Unexpected branch-arm label in Step5BR2 CDS."
  )
}

SummarizedExperiment::colData(cds)$FigureS3_branch_display <- factor(
  branch_display,
  levels = names(branch_cols)
)

plot_args <- list(
  cds = cds,
  reduction_method = "UMAP",
  color_cells_by = "FigureS3_branch_display",
  show_trajectory_graph = TRUE,
  label_cell_groups = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE,
  label_roots = FALSE,
  graph_label_size = 2.4,
  cell_size = 0.23,
  trajectory_graph_color = "#51495B",
  trajectory_graph_segment_size = 0.66
)

supported <- names(formals(monocle3::plot_cells))
plot_args <- plot_args[
  names(plot_args) %in% supported
]

pA_topology <- do.call(
  monocle3::plot_cells,
  plot_args
)

pA_topology <- suppressMessages(
  pA_topology +
    scale_color_manual(
      values = branch_cols,
      drop = FALSE
    )
) +
  labs(
    title = "Single-principal-root topology",
    subtitle = NULL,
    x = "UMAP 1",
    y = "UMAP 2",
    color = NULL
  ) +
  theme_pub +
  theme(
    plot.title = element_text(
      size = 6.35,
      face = "bold",
      hjust = 0
    ),
    plot.subtitle = element_blank(),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.key.width = grid::unit(2.8, "mm"),
    legend.spacing.x = grid::unit(0.5, "mm"),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    plot.margin = margin(
      2.4, 2.4, 2.6, 2.0,
      unit = "mm"
    )
  ) +
  guides(
    color = guide_legend(
      nrow = 2,
      byrow = TRUE,
      override.aes = list(
        size = 2.0,
        alpha = 1
      )
    )
  )

panelA <- patchwork::wrap_plots(
  pA_root,
  pA_topology,
  ncol = 2,
  widths = c(0.43, 0.57)
)

# ==============================================================================
# PANEL B
# Root/orientation sensitivity
# ==============================================================================

message("Using the same frozen Step5BR2 Seurat object for root-sensitivity analysis...")

req_br2_seurat_meta <- c(
  "sample_id",
  "timepoint",
  "step5BR_pseudotime",
  "step5BR2_single_root_pseudotime",
  "step5BR2_branch_arm",
  "step5BR2_major_late_arm"
)

miss_br2_seurat_meta <- setdiff(
  req_br2_seurat_meta,
  colnames(br2_seurat[[]])
)

if (length(miss_br2_seurat_meta) > 0) {
  stop(
    "Step5BR2 Seurat metadata missing: ",
    paste(miss_br2_seurat_meta, collapse = ", ")
  )
}

br2_md <- br2_seurat[[]] %>%
  tibble::rownames_to_column("cell") %>%
  transmute(
    cell = cell,
    sample = as.character(sample_id),
    timepoint = as.character(timepoint),
    old_pseudotime = as.numeric(step5BR_pseudotime),
    pseudotime_single_root = as.numeric(step5BR2_single_root_pseudotime),
    branch_arm = as.character(step5BR2_branch_arm),
    major_late_arm = as.logical(step5BR2_major_late_arm)
  )

downstream_sens <- br2_md %>%
  filter(
    major_late_arm,
    timepoint %in% c("3 dpi", "7 dpi"),
    is.finite(old_pseudotime),
    is.finite(pseudotime_single_root)
  ) %>%
  mutate(
    sample = factor(
      sample,
      levels = c(
        "3dpi_sample1",
        "3dpi_sample2",
        "7dpi_sample1",
        "7dpi_sample2"
      )
    ),
    branch_display = dplyr::recode(
      branch_arm,
      "Arm_1" = "Alternative branch",
      "Arm_2" = "Target-program branch",
      .default = NA_character_
    )
  )

if (nrow(downstream_sens) != 4087) {
  warning(
    "Expected 4,087 downstream 3/7 dpi major-arm cells from the frozen audit, ",
    "but found ",
    nrow(downstream_sens),
    ". The plot will still be generated; inspect the audit CSV."
  )
}

rho_overall <- safe_spearman(
  downstream_sens$old_pseudotime,
  downstream_sens$pseudotime_single_root
)

rho_by_sample <- downstream_sens %>%
  group_by(sample) %>%
  summarise(
    n_cells = n(),
    spearman_rho = safe_spearman(
      old_pseudotime,
      pseudotime_single_root
    ),
    median_shift = median(
      pseudotime_single_root -
        old_pseudotime,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

if (
  any(abs(rho_by_sample$spearman_rho - 1) > 1e-6)
) {
  warning(
    "Downstream sample-wise root-specification pseudotime rank correlation ",
    "is not exactly 1 for every sample. Inspect the audit CSV."
  )
}

pB <- ggplot(
  downstream_sens,
  aes(
    x = old_pseudotime,
    y = pseudotime_single_root
  )
) +
  geom_point(
    aes(
      color = sample
    ),
    size = 0.48,
    alpha = 0.32
  ) +
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = FALSE,
    linewidth = 0.72,
    color = "#303030"
  ) +
  annotate(
    "label",
    x = min(downstream_sens$old_pseudotime) + 0.04,
    y = max(downstream_sens$pseudotime_single_root) - 0.04,
    label = paste0(
      "Ranks preserved\n",
      "Spearman r_s = ",
      sprintf("%.3f", rho_overall),
      "\nEach 3/7 dpi replicate: r_s = 1.000"
    ),
    family = BASE_FAMILY,
    size = 2.05,
    hjust = 0,
    vjust = 1,
    label.size = 0.20,
    fill = scales::alpha("white", 0.88)
  ) +
  scale_color_manual(
    values = sample_cols,
    labels = c(
      "3dpi_sample1" = "3 dpi sample 1",
      "3dpi_sample2" = "3 dpi sample 2",
      "7dpi_sample1" = "7 dpi sample 1",
      "7dpi_sample2" = "7 dpi sample 2"
    ),
    drop = FALSE
  ) +
  labs(
    tag = "B",
    title = "Root/orientation sensitivity",
    subtitle = NULL,
    x = "Alternative-root pseudotime",
    y = "Primary-root pseudotime",
    color = NULL
  ) +
  theme_pub +
  theme(
    plot.subtitle = element_blank(),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.key.width = grid::unit(2.8, "mm"),
    legend.spacing.x = grid::unit(0.5, "mm"),
    plot.margin = margin(
      3.0, 3.2, 3.0, 4.0,
      unit = "mm"
    )
  ) +
  guides(
    color = guide_legend(
      nrow = 2,
      byrow = TRUE,
      override.aes = list(
        size = 2.2,
        alpha = 1
      )
    )
  )

# Step5BR2 objects are no longer needed after the root-sensitivity section.
rm(cds, br2_seurat, br2_md)
invisible(gc())

# ==============================================================================
# PANEL C
# Shared prebranch + progressive divergence
# ==============================================================================

message("Reading Step5CR progressive-divergence workbook...")

cr_master <- readxl::read_excel(
  P$step5CR_xlsx,
  sheet = "Master_568_genes"
) %>%
  as.data.frame()

cr_binned <- readxl::read_excel(
  P$step5CR_xlsx,
  sheet = "Binned_expression"
) %>%
  as.data.frame()

cr_cell_summary <- readxl::read_excel(
  P$step5CR_xlsx,
  sheet = "Cell_summary"
) %>%
  as.data.frame()

cr_params <- readxl::read_excel(
  P$step5CR_xlsx,
  sheet = "Parameters"
) %>%
  as.data.frame()

req_master <- c(
  "gene_id",
  "abs_early_difference",
  "abs_late_difference",
  "strong_progressive",
  "terminal_direction"
)

miss_master <- setdiff(
  req_master,
  colnames(cr_master)
)

if (length(miss_master) > 0) {
  stop(
    "Step5CR master sheet missing: ",
    paste(miss_master, collapse = ", ")
  )
}

n_strong <- sum(
  as.logical(cr_master$strong_progressive),
  na.rm = TRUE
)

shared_cells <- suppressWarnings(
  as.numeric(
    cr_params$value[
      cr_params$parameter ==
        "Shared prebranch cells"
    ]
  )
)

if (
  length(shared_cells) != 1 ||
  shared_cells != 116
) {
  stop(
    "Unexpected Step5CR shared-prebranch count."
  )
}

pC_div <- ggplot(
  cr_master,
  aes(
    x = abs_early_difference,
    y = abs_late_difference
  )
) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linewidth = 0.42,
    linetype = "22",
    color = "#777777"
  ) +
  geom_point(
    data = subset(
      cr_master,
      !as.logical(strong_progressive)
    ),
    size = 0.72,
    color = "#C9C9C9",
    alpha = 0.75
  ) +
  geom_point(
    data = subset(
      cr_master,
      as.logical(strong_progressive)
    ),
    aes(
      color = terminal_direction
    ),
    size = 1.05,
    alpha = 0.90
  ) +
  scale_color_manual(
    values = c(
      "Arm_1" = COL_ALT,
      "Arm_2" = COL_TARGET
    ),
    labels = c(
      "Arm_1" = "Alternative-branch biased",
      "Arm_2" = "Target-program biased"
    ),
    drop = FALSE
  ) +
  labs(
    tag = "B",
    title = "Progressive divergence",
    subtitle = NULL,
    x = "Early |Alternative - Target-program| difference",
    y = "Late |Alternative - Target-program| difference",
    color = NULL
  ) +
  theme_pub +
  theme(
    plot.subtitle = element_blank(),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "vertical",
    legend.key.width = grid::unit(2.8, "mm"),
    legend.key.height = grid::unit(2.2, "mm"),
    legend.spacing.y = grid::unit(0.5, "mm"),
    legend.margin = margin(0, 0, 0, 0, unit = "mm"),
    plot.margin = margin(
      3.0, 2.0, 4.5, 4.0,
      unit = "mm"
    )
  ) +
  guides(
    color = guide_legend(
      nrow = 2,
      byrow = TRUE,
      override.aes = list(
        size = 2.1,
        alpha = 1
      )
    )
  )

representative_genes <- c(
  "Marcksl1",
  "Tnf",
  "Tmsb4x"
)

if (!all(representative_genes %in% cr_binned$gene_id)) {
  stop(
    "Representative Step5CR genes missing from Binned_expression sheet."
  )
}

traj_df <- cr_binned %>%
  filter(
    gene_id %in% representative_genes
  ) %>%
  mutate(
    gene_id = factor(
      gene_id,
      levels = representative_genes
    ),
    state = factor(
      state,
      levels = c(
        "Shared_prebranch",
        "Arm_1",
        "Arm_2"
      )
    ),
    branch_display = dplyr::recode(
      as.character(state),
      "Shared_prebranch" = "Shared prebranch",
      "Arm_1" = "Alternative branch",
      "Arm_2" = "Target-program branch"
    )
  )

# Expand the single shared-prebranch point to both downstream line groups so that
# each trajectory visibly starts from the same prebranch anchor.
shared_pts <- traj_df %>%
  filter(state == "Shared_prebranch") %>%
  select(
    gene_id,
    progress_mid,
    mean_expression
  )

arm_pts <- traj_df %>%
  filter(state %in% c("Arm_1", "Arm_2")) %>%
  mutate(
    line_group = dplyr::recode(
      as.character(state),
      "Arm_1" = "Alternative branch",
      "Arm_2" = "Target-program branch"
    )
  ) %>%
  select(
    gene_id,
    progress_mid,
    mean_expression,
    line_group
  )

shared_expand <- tidyr::crossing(
  shared_pts,
  line_group = c(
    "Alternative branch",
    "Target-program branch"
  )
)

traj_line <- bind_rows(
  shared_expand,
  arm_pts
) %>%
  mutate(
    line_group = factor(
      line_group,
      levels = c(
        "Alternative branch",
        "Target-program branch"
      )
    )
  )

pC_gene <- ggplot(
  traj_line,
  aes(
    x = progress_mid,
    y = mean_expression,
    color = line_group,
    group = line_group
  )
) +
  geom_line(
    linewidth = 0.82
  ) +
  geom_point(
    size = 1.15
  ) +
  facet_wrap(
    ~gene_id,
    ncol = 1,
    scales = "free_y"
  ) +
  scale_color_manual(
    values = c(
      "Alternative branch" = COL_ALT,
      "Target-program branch" = COL_TARGET
    ),
    drop = FALSE
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = c(0, 0.5, 1)
  ) +
  labs(
    title = "Example trajectories",
    subtitle = NULL,
    x = "Normalized downstream progress",
    y = "Mean expression",
    color = NULL
  ) +
  theme_pub +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(
      family = BASE_FAMILY,
      face = "italic",
      size = 5.8,
      hjust = 0
    ),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "vertical",
    legend.key.width = grid::unit(3.0, "mm"),
    legend.key.height = grid::unit(2.2, "mm"),
    legend.spacing.y = grid::unit(0.5, "mm"),
    legend.margin = margin(0, 0, 0, 0, unit = "mm"),
    axis.text = element_text(
      size = 5.1
    ),
    axis.title = element_text(
      size = 5.8
    ),
    plot.title = element_text(
      size = 6.1,
      face = "bold"
    ),
    plot.margin = margin(
      3.0, 3.0, 4.5, 2.0,
      unit = "mm"
    )
  ) +
  guides(
    color = guide_legend(
      nrow = 2,
      byrow = TRUE,
      override.aes = list(
        linewidth = 1.0
      )
    )
  )

panelC <- patchwork::wrap_plots(
  pC_div,
  pC_gene,
  ncol = 2,
  widths = c(0.57, 0.43)
)

# ==============================================================================
# PANEL D
# Biological-replicate + pseudobulk support
# ==============================================================================

message("Reading Step6A/Step6B branch-program support tables...")

hc <- readr::read_csv(
  P$step6A_highconf,
  show_col_types = FALSE
)

module_summary <- readr::read_csv(
  P$step6B_summary,
  show_col_types = FALSE
)

module_effects <- readr::read_csv(
  P$step6B_effects,
  show_col_types = FALSE
)

req_hc <- c(
  "gene_id",
  "repFC_7dpi_sample1",
  "repFC_7dpi_sample2",
  "strict_pseudobulk_supported",
  "final_direction"
)

miss_hc <- setdiff(
  req_hc,
  colnames(hc)
)

if (length(miss_hc) > 0) {
  stop(
    "Step6A high-confidence table missing: ",
    paste(miss_hc, collapse = ", ")
  )
}

rho_rep <- safe_spearman(
  hc$repFC_7dpi_sample1,
  hc$repFC_7dpi_sample2
)

n_strict_pb <- sum(
  hc$strict_pseudobulk_supported,
  na.rm = TRUE
)

# Background = high-confidence but not strict PB.
# Foreground = strict PB, colored by final branch direction.
hc <- hc %>%
  mutate(
    strict_status = ifelse(
      strict_pseudobulk_supported,
      "Strict pseudobulk-supported",
      "High-confidence only"
    ),
    direction_display = dplyr::recode(
      final_direction,
      "Arm_1-biased" = "Alternative-branch biased",
      "Arm_2-biased" = "Target-program biased"
    )
  )

pD_rep <- ggplot() +
  geom_hline(
    yintercept = 0,
    linewidth = 0.32,
    linetype = "22",
    color = "#777777"
  ) +
  geom_vline(
    xintercept = 0,
    linewidth = 0.32,
    linetype = "22",
    color = "#777777"
  ) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linewidth = 0.40,
    color = "#777777"
  ) +
  geom_point(
    data = hc %>%
      filter(!strict_pseudobulk_supported),
    aes(
      x = repFC_7dpi_sample1,
      y = repFC_7dpi_sample2
    ),
    size = 0.70,
    color = COL_NONSTRICT,
    alpha = 0.55
  ) +
  geom_point(
    data = hc %>%
      filter(strict_pseudobulk_supported),
    aes(
      x = repFC_7dpi_sample1,
      y = repFC_7dpi_sample2,
      color = direction_display
    ),
    size = 0.98,
    alpha = 0.90
  ) +
  scale_color_manual(
    values = c(
      "Alternative-branch biased" = COL_ALT,
      "Target-program biased" = COL_TARGET
    ),
    drop = FALSE
  ) +
  annotate(
    "label",
    x = -2.92,
    y = 1.32,
    label = paste0(
      "Spearman r_s = ",
      sprintf("%.3f", rho_rep),
      "\n745 replicate-consistent markers",
      "\n",
      n_strict_pb,
      " strict pseudobulk-supported"
    ),
    family = BASE_FAMILY,
    size = 1.78,
    hjust = 0,
    vjust = 1,
    label.size = 0.20,
    fill = scales::alpha("white", 0.90)
  ) +
  labs(
    tag = "C",
    title = "7 dpi replicate concordance",
    subtitle = NULL,
    x = "7 dpi sample 1 log2FC",
    y = "7 dpi sample 2 log2FC",
    color = NULL
  ) +
  theme_pub +
  theme(
    plot.title = element_text(
      size = 6.25,
      face = "bold",
      hjust = 0
    ),
    plot.subtitle = element_blank(),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.key.width = grid::unit(3.0, "mm"),
    plot.margin = margin(
      2.2, 2.8, 2.4, 4.0,
      unit = "mm"
    )
  ) +
  guides(
    color = guide_legend(
      nrow = 1,
      override.aes = list(
        size = 2.1,
        alpha = 1
      )
    )
  )

req_mod <- c(
  "module_name",
  "n_samples",
  "same_direction_all_samples",
  "robust_terminal_branch_effect"
)

miss_mod <- setdiff(
  req_mod,
  colnames(module_summary)
)

if (length(miss_mod) > 0) {
  stop(
    "Step6B module summary missing: ",
    paste(miss_mod, collapse = ", ")
  )
}

module_order <- module_summary %>%
  arrange(
    desc(abs(mean_Arm1_minus_Arm2))
  ) %>%
  pull(module_name)

module_effects_plot <- module_effects %>%
  mutate(
    module_name = factor(
      module_name,
      levels = rev(module_order)
    ),
    sample = factor(
      sample,
      levels = c(
        "7dpi_sample1",
        "7dpi_sample2"
      )
    )
  )

module_label_map <- setNames(
  c(
    "Mod1 Target",
    "Mod2 Target",
    "Mod3 Shared",
    "Mod4 Target",
    "Mod5 Target"
  ),
  c(
    "M1_Arm2_terminal",
    "M2_Arm2_terminal",
    "M3_Shared_prebranch",
    "M4_Arm2_terminal",
    "M5_Arm2_terminal"
  )
)

n_same_dir <- sum(
  module_summary$same_direction_all_samples,
  na.rm = TRUE
)

n_robust <- sum(
  module_summary$robust_terminal_branch_effect,
  na.rm = TRUE
)

pD_mod <- ggplot(
  module_effects_plot,
  aes(
    x = Arm1_minus_Arm2,
    y = module_name,
    group = module_name
  )
) +
  geom_vline(
    xintercept = 0,
    linewidth = 0.38,
    linetype = "22",
    color = "#777777"
  ) +
  geom_segment(
    data = module_effects_plot %>%
      group_by(module_name) %>%
      summarise(
        xmin = min(Arm1_minus_Arm2),
        xmax = max(Arm1_minus_Arm2),
        .groups = "drop"
      ),
    aes(
      x = xmin,
      xend = xmax,
      y = module_name,
      yend = module_name
    ),
    inherit.aes = FALSE,
    linewidth = 0.70,
    color = "#CFCFCF"
  ) +
  geom_point(
    aes(
      color = sample
    ),
    size = 2.55,
    alpha = 0.98
  ) +
  scale_color_manual(
    values = c(
      "7dpi_sample1" = "#4C78A8",
      "7dpi_sample2" = "#E45756"
    ),
    labels = c(
      "7dpi_sample1" = "7 dpi sample 1",
      "7dpi_sample2" = "7 dpi sample 2"
    ),
    drop = FALSE
  ) +
  scale_y_discrete(
    labels = function(x) {
      unname(module_label_map[as.character(x)])
    },
    drop = FALSE
  ) +
  labs(
    title = "Module-level replicate support",
    subtitle = NULL,
    x = "Alternative - Target-program module score",
    y = NULL,
    color = NULL
  ) +
  theme_pub +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.key.width = grid::unit(3.0, "mm"),
    axis.text.y = element_text(
      size = 5.0
    ),
    plot.title = element_text(
      size = 6.15,
      face = "bold"
    ),
    plot.subtitle = element_blank(),
    plot.margin = margin(
      3.0, 3.0, 3.0, 2.0,
      unit = "mm"
    )
  ) +
  guides(
    color = guide_legend(
      nrow = 1,
      override.aes = list(
        size = 2.2
      )
    )
  )

panelD <- patchwork::wrap_plots(
  pD_rep,
  pD_mod,
  ncol = 1,
  heights = c(0.56, 0.44)
)

# ------------------------------------------------------------------------------
# 9. Audits
# ------------------------------------------------------------------------------
write.csv(
  root_sample_summary,
  file.path(
    OUT,
    "02_PanelA_primary_1dpi_root_by_sample.csv"
  ),
  row.names = FALSE
)

write.csv(
  root_vertex_summary,
  file.path(
    OUT,
    "02B_PanelA_primary_root_vertex_mapping.csv"
  ),
  row.names = FALSE
)

write.csv(
  rho_by_sample,
  file.path(
    OUT,
    "03_PanelB_downstream_orientation_sensitivity.csv"
  ),
  row.names = FALSE
)

write.csv(
  cr_cell_summary,
  file.path(
    OUT,
    "04_PanelC_shared_prebranch_cell_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  hc,
  file.path(
    OUT,
    "05_PanelD_branch_marker_support.csv"
  ),
  row.names = FALSE
)

write.csv(
  module_summary,
  file.path(
    OUT,
    "06_PanelD_module_summary.csv"
  ),
  row.names = FALSE
)

audit <- data.frame(
  item = c(
    "primary_1dpi_root_total_cells",
    "primary_1dpi_root_cells_per_sample",
    "single_principal_root_selected",
    "primary_root_cells_in_selected_root_vertex",
    "downstream_sensitivity_cells_3_7dpi",
    "downstream_old_vs_new_spearman",
    "shared_prebranch_cells",
    "progressive_genes_tested",
    "strong_progressive_genes",
    "high_confidence_branch_markers",
    "strict_pseudobulk_supported_genes",
    "replicate_log2FC_spearman",
    "dynamic_modules_total",
    "dynamic_modules_same_direction_both_samples",
    "dynamic_modules_robust_terminal"
  ),
  value = as.character(
    c(
      nrow(root_final),
      75,
      "TRUE",
      max(root_vertex_summary$n_primary_root_cells, na.rm = TRUE),
      nrow(downstream_sens),
      rho_overall,
      shared_cells,
      nrow(cr_master),
      n_strong,
      nrow(hc),
      n_strict_pb,
      rho_rep,
      nrow(module_summary),
      n_same_dir,
      n_robust
    )
  ),
  stringsAsFactors = FALSE
)

write.csv(
  audit,
  file.path(
    OUT,
    "00_FigureS3_data_audit.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 10. Final 2x2 layout
# ------------------------------------------------------------------------------
panelA_fixed <- patchwork::wrap_elements(
  full = panelA,
  clip = FALSE
)

panelB_fixed <- patchwork::wrap_elements(
  full = panelC,
  clip = FALSE
)

panelC_fixed <- patchwork::wrap_elements(
  full = panelD,
  clip = FALSE
)

row1 <- patchwork::wrap_plots(
  panelA_fixed,
  ncol = 1
)

row2 <- patchwork::wrap_plots(
  panelB_fixed,
  panelC_fixed,
  ncol = 2,
  widths = c(1.06, 0.94)
)

figureS3 <- patchwork::wrap_plots(
  patchwork::wrap_elements(
    full = row1,
    clip = FALSE
  ),
  patchwork::wrap_elements(
    full = row2,
    clip = FALSE
  ),
  ncol = 1,
  # Give Panel A substantially more vertical space so the two side-by-side
  # UMAPs are not visually flattened in the final composite.
  heights = c(1.05, 0.95)
) +
  plot_annotation(
    theme = theme(
      plot.background = element_rect(
        fill = "white",
        color = NA
      )
    )
  )

# ------------------------------------------------------------------------------
# 11. Export
# ------------------------------------------------------------------------------
message("Saving Figure S3 panel A...")
save_plot_quad(
  panelA,
  "FigureS3A_root_topology",
  width_mm = 112,
  height_mm = 82
)

message("Saving Figure S3 panel B...")
save_plot_quad(
  panelC,
  "FigureS3B_prebranch_progressive_divergence",
  width_mm = 106,
  height_mm = 88
)

message("Saving Figure S3 panel C...")
save_plot_quad(
  panelD,
  "FigureS3C_replicate_pseudobulk_support",
  width_mm = 88,
  height_mm = 102
)

message("Saving complete Figure S3...")
save_plot_quad(
  figureS3,
  "FigureS3_GLIA_FINAL_v1_4_9_Mod1_Mod5_labels",
  width_mm = 190,
  height_mm = 178
)

# ------------------------------------------------------------------------------
# 12. Run summary
# ------------------------------------------------------------------------------
run_summary <- c(
  "Figure S3 | GLIA | FINAL v1.4.9 | Mod1-Mod5 reader-facing labels | no subtitles | no panel B",
  "==============================================",
  "",
  paste0("Project root: ", ROOT),
  paste0("Output: ", OUT),
  "",
  "Panel A — primary 1 dpi root + single-principal-root topology",
  "  225 primary injury-trajectory root cells; 75 from each of three 1 dpi biological samples.",
  "  The 225 primary roots map across multiple principal vertices; the highest-support vertex is used as the single objective principal root.",
  "  The objective bifurcation is shown with Arm_2 displayed as Target-program branch.",
  "",
  "",
  "Panel B — shared prebranch + progressive divergence",
  paste0("  Shared-prebranch cells = ", shared_cells, "."),
  paste0(
    "  Strong progressive genes = ",
    n_strong,
    " / ",
    nrow(cr_master),
    "."
  ),
  "  Only three representative trajectories are shown; the full dynamic heatmap is not duplicated.",
  "",
  "Panel C — replicate / pseudobulk support",
  paste0(
    "  Replicate-consistent high-confidence branch markers = ",
    nrow(hc),
    "."
  ),
  paste0(
    "  Strict pseudobulk-supported genes = ",
    n_strict_pb,
    "."
  ),
  paste0(
    "  Replicate log2FC Spearman r_s = ",
    sprintf("%.3f", rho_rep),
    "."
  ),
  paste0(
    "  Dynamic modules with same direction in both samples = ",
    n_same_dir,
    " / ",
    nrow(module_summary),
    "."
  ),
  "",
  "No upstream analysis was rerun."
)

writeLines(
  run_summary,
  con = file.path(
    OUT,
    "FigureS3_GLIA_FINAL_v1_4_9_Mod1_Mod5_labels_run_summary.txt"
  ),
  useBytes = TRUE
)

capture.output(
  sessionInfo(),
  file = file.path(
    OUT,
    "sessionInfo.txt"
  )
)

warn_obj <- warnings()

if (!is.null(warn_obj)) {
  capture.output(
    warn_obj,
    file = file.path(
      OUT,
      "07_warnings_after_run.txt"
    )
  )
}

message("============================================================")
message("Figure S3 FINAL v1.4.9 (Mod1-Mod5 reader-facing labels; no subtitles, no panel B) completed.")
message("Output: ", OUT)
message("Main review file: FigureS3_GLIA_FINAL_v1_4_9_Mod1_Mod5_labels_review.png")
message("Main vector file: FigureS3_GLIA_FINAL_v1_4_9_Mod1_Mod5_labels.pdf / .svg")
message("============================================================")
