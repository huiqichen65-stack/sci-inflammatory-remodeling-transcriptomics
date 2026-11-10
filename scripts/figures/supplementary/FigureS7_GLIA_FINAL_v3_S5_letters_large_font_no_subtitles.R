# =============================================================================
# FigureS7_GLIA_FINAL_v3_S5_letters_large_font_no_subtitles.R
# GLIA | Supplementary Figure S7 | Reviewer 2, Major Comment 2
#
# Project convention (identical in structure to FigureS7_GLIA_FINAL_v1_FINAL.R):
#   Script:  <PROJECT_ROOT>/script/FigureS7_GLIA_FINAL_v3_S5_letters_large_font_no_subtitles.R
#   Inputs:  xjz/04_results/step6F_GSE162610_dissociation_activation_sensitivity/
#            xjz/04_results/step6FRH_GSE162610_Harmony_matched_trajectory_sensitivity/
#                 31N_source_matched_Marsh_paired_trajectory/
#   Outputs: xjz/09_final_figures/GLIA/FigureS7_FINAL_v3_S5_letters_large_font_no_subtitles/
#
# Design: use the FINAL main-manuscript Figure 2/3 palette:
#   Target: #B64B45 | general Figure 2 blue: #4C78A8 (Alternative)
#   Figure 2 orange: #F28E2B (new, additional downstream branch)
#   Figure 2 shared-prebranch purple: #7A6AA6 (NOT used for a new branch)
#   Root-side and unassigned cells: separate neutral grays.
# The role-based colors in panels A-C are identical.
#
# Panel A: Published Marsh activation-signature scores, sample x branch.
# Panel B: Source-matched full-gene and Marsh-masked principal graph, with the
#          THIRD downstream arm colored explicitly (rather than gray).
# Panel C: Per-sample correspondence of original branch cells after masking.
#
# v2 FIX: Panel B uses fixed facet scales with coord_fixed(ratio=1);
#          avoid ggplot2 incompatible free facets / fixed aspect ratio.
# PLOTTING ONLY. Input matrices, pseudotime, original labels, and CDS objects are
# READ-ONLY. No Harmony, Monocle relearning, feature re-scoring, or new tests.
# Every publication file and audit goes ONLY into 09_final_figures/GLIA.
# =============================================================================

options(stringsAsFactors = FALSE)
root_dir <- Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
if (!dir.exists(root_dir)) stop("Project directory not found: ", root_dir)

step6f <- file.path(root_dir, "04_results",
  "step6F_GSE162610_dissociation_activation_sensitivity")
step31n <- file.path(root_dir, "04_results",
  "step6FRH_GSE162610_Harmony_matched_trajectory_sensitivity",
  "31N_source_matched_Marsh_paired_trajectory")
out_dir <- file.path(root_dir, "09_final_figures", "GLIA", "FigureS7_FINAL_v3_S5_letters_large_font_no_subtitles")

base_name <- "FigureS7_GLIA_FINAL_v3_S5_letters_large_font_no_subtitles"
input_paths <- c(
  score_cells = file.path(step6f, "tables", "03_cell_level_Marsh_activation_scores.csv"),
  score_summary = file.path(step6f, "tables", "04_activation_score_sample_by_branch_summary.csv"),
  control_cells = file.path(step31n, "tables", "11_all_genes_control_cell_level_correspondence.csv"),
  masked_cells = file.path(step31n, "tables", "11_Marsh_activation_masked_cell_level_correspondence.csv"),
  recovery = file.path(step31n, "tables", "14_paired_sample_level_recovery_comparison.csv"),
  topology = file.path(step31n, "tables", "16_paired_final_topology_summary.csv"),
  control_cds = file.path(step31n, "RDS", "31N_all_genes_control_relearned_Monocle_CDS.rds"),
  masked_cds = file.path(step31n, "RDS", "31N_Marsh_activation_masked_relearned_Monocle_CDS.rds")
)
missing_files <- input_paths[!file.exists(input_paths)]
if (length(missing_files)) stop(
  "Missing frozen Step6F / Step31N input(s):\n",
  paste(missing_files, collapse = "\n"),
  "\nDo not rerun the analyses; check the original files and paths.")

packages <- c("ggplot2", "patchwork", "monocle3", "igraph", "SingleCellExperiment")
not_installed <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(not_installed)) stop("Required packages are missing: ",
                                paste(not_installed, collapse = ", "))
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(grid)
})

read_step <- function(p) utils::read.csv(p, check.names = FALSE, stringsAsFactors = FALSE)
assert_cols <- function(d, fields, label) {
  absent <- setdiff(fields, names(d))
  if (length(absent)) stop(label, " missing columns: ", paste(absent, collapse = ", "))
}
near <- function(x, y, tolerance = 1e-6) all(is.finite(x), is.finite(y),
                                              abs(x-y) <= tolerance)

scores <- read_step(input_paths[["score_cells"]])
scores_summary <- read_step(input_paths[["score_summary"]])
control <- read_step(input_paths[["control_cells"]])
masked <- read_step(input_paths[["masked_cells"]])
recovery <- read_step(input_paths[["recovery"]])
topology <- read_step(input_paths[["topology"]])

assert_cols(scores, c("cell_id", "sample_id", "branch",
  "score_all_CNS_activation", "score_micro_myeloid_activation",
  "score_activation_union"), "Cell-level Marsh score table")
assert_cols(scores_summary, c("sample_id", "branch", "signature",
  "n_cells", "mapped_genes", "median", "q25", "q75"), "Marsh score summary")
required_cell_columns <- c("barcode", "sample_id", "timepoint", "original_arm",
                           "new_arm", "aligned_biological_arm", "UMAP_1", "UMAP_2")
assert_cols(control, required_cell_columns, "31N full-gene cell correspondence")
assert_cols(masked, required_cell_columns, "31N Marsh-masked cell correspondence")
assert_cols(recovery, c("condition", "sample_id", "original_arm",
                        "n_original", "n_recovered", "recall"), "31N recovery table")
assert_cols(topology, c("condition", "branchpoint", "graph_nodes", "graph_edges"),
            "31N topology table")

if (anyDuplicated(control$barcode) || anyDuplicated(masked$barcode) ||
    nrow(control) != 14005L || !setequal(control$barcode, masked$barcode)) {
  stop("31N input barcodes do not match the frozen 14,005 trajectory cells")
}
masked <- masked[match(control$barcode, masked$barcode), , drop = FALSE]
if (!identical(control$barcode, masked$barcode) ||
    !identical(control$original_arm, masked$original_arm) ||
    !identical(control$sample_id, masked$sample_id) ||
    !identical(control$timepoint, masked$timepoint)) {
  stop("Full-gene and Marsh-masked results do not have identical frozen labels")
}
if (any(!is.finite(control$UMAP_1)) || any(!is.finite(control$UMAP_2)) ||
    any(!is.finite(masked$UMAP_1)) || any(!is.finite(masked$UMAP_2))) {
  stop("Missing or nonfinite UMAP coordinates in 31N data")
}

# Main-Figure-2/3 color system, fixed in every Figure S7 panel.
COL <- list(
  target = "#B64B45",           # exact Fig 2/3 Target program red
  alternative = "#4C78A8",      # exact Fig 2 blue
  third = "#F28E2B",            # exact Fig 2 orange; new branch ONLY
  shared = "#7A6AA6",           # exact Fig 2/3 shared-prebranch purple (reserved)
  root = "#B8BFC5",             # neutral, not equivalent to shared-prebranch
  unassigned = "#E6E6E6",       # neutral lighter gray
  ink = "#202020",
  grid = "#E4E4E4"
)
PAL <- c("Alternative" = COL$alternative,
         "Target" = COL$target,
         "Additional downstream" = COL$third,
         "Prebranch/root side" = COL$root,
         "Unassigned" = COL$unassigned)
ARM_LABS <- c(Arm_1 = "Alternative", Arm_2 = "Target")
SAMPLE_IDS <- c("7dpi_sample1", "7dpi_sample2")
SAMPLE_LABS <- c("7 dpi: sample 1", "7 dpi: sample 2")
SIGNATURE_LABELS <- c(
  all_CNS_activation = "All CNS activation",
  micro_myeloid_activation = "Micro/myeloid activation",
  activation_union = "Activation union"
)

# --------------------------------------- Panel A: true cell-level scores.
sa <- scores[scores$sample_id %in% SAMPLE_IDS &
               scores$branch %in% names(ARM_LABS), , drop = FALSE]
if (nrow(sa) != 2541L || anyDuplicated(sa$cell_id)) {
  stop("Expected 2,541 original 7-dpi terminal-branch cells in Step6F score data")
}
score_long <- do.call(rbind, lapply(names(SIGNATURE_LABELS), function(sig) {
  data.frame(cell_id = sa$cell_id,
    sample_id = sa$sample_id,
    branch = unname(ARM_LABS[sa$branch]),
    signature = unname(SIGNATURE_LABELS[[sig]]),
    value = as.numeric(sa[[paste0("score_", sig)]]))
}))
if (any(!is.finite(score_long$value))) stop("Some Marsh score values are nonfinite")
score_long$sample <- factor(score_long$sample_id,
                            levels = SAMPLE_IDS, labels = SAMPLE_LABS)
score_long$signature <- factor(score_long$signature,
                               levels = unname(SIGNATURE_LABELS))
score_long$branch <- factor(score_long$branch, levels = c("Alternative", "Target"))

# Cross-check raw scores against frozen summary; check medians, IQR and counts.
for (i in seq_len(nrow(scores_summary))) {
  s <- scores_summary[i, , drop = FALSE]
  if (!s$signature %in% names(SIGNATURE_LABELS) ||
      !s$sample_id %in% SAMPLE_IDS || !s$branch %in% names(ARM_LABS)) next
  idx <- sa$sample_id == s$sample_id & sa$branch == s$branch
  v <- as.numeric(sa[[paste0("score_", s$signature)]][idx])
  if (length(v) != s$n_cells ||
      !near(stats::median(v), as.numeric(s$median)) ||
      !near(unname(stats::quantile(v, .25)), as.numeric(s$q25)) ||
      !near(unname(stats::quantile(v, .75)), as.numeric(s$q75))) {
    stop("Marsh score summary disagrees with cell-level data at row ", i)
  }
}

# ------------------------ Panel B: actual frozen Monocle3 principal graphs.
# Never draw an invented tree or use the other condition's vertex positions.
extract_graph <- function(path, df, condition_id, condition_label) {
  cds <- readRDS(path)
  graph <- monocle3::principal_graph(cds)[["UMAP"]]
  if (is.null(graph)) stop("UMAP principal graph unavailable: ", path)
  vertices <- igraph::V(graph)$name
  xy <- cds@principal_graph_aux[["UMAP"]]$dp_mst
  if (is.null(xy)) stop("Principal-graph node positions missing in: ", path)
  xy <- as.matrix(xy)
  if (nrow(xy) == 2L && all(vertices %in% colnames(xy))) {
    pos <- t(xy[, vertices, drop = FALSE])
  } else if (ncol(xy) == 2L && all(vertices %in% rownames(xy))) {
    pos <- xy[vertices, , drop = FALSE]
  } else {
    stop("Cannot match Monocle3 graph vertices to their positions: ", path)
  }
  colnames(pos) <- c("x", "y")
  if (any(!is.finite(pos))) stop("Nonfinite principal graph vertex positions")

  edges <- igraph::as_data_frame(graph, what = "edges")
  ix <- match(edges$from, vertices)
  iy <- match(edges$to, vertices)
  if (anyNA(ix) || anyNA(iy)) stop("Graph edges have unknown vertices")
  segments <- data.frame(condition = condition_label,
                         x = pos[ix, 1], y = pos[ix, 2],
                         xend = pos[iy, 1], yend = pos[iy, 2])
  topo <- topology[topology$condition == condition_id, , drop = FALSE]
  if (nrow(topo) != 1L) stop("Missing/duplicate 31N topology: ", condition_id)
  bp <- topo$branchpoint[1]
  if (!bp %in% vertices) stop("Objective branchpoint not found in CDS: ", bp)
  marks <- data.frame(condition = condition_label,
                      x = pos[bp, 1], y = pos[bp, 2], label = "BP")
  if (igraph::vcount(graph) != topo$graph_nodes[1] ||
      igraph::ecount(graph) != topo$graph_edges[1]) {
    stop("Graph nodes/edges disagree with 31N topology table: ", condition_id)
  }

  u <- SingleCellExperiment::reducedDim(cds, "UMAP")
  if (is.null(rownames(u))) rownames(u) <- colnames(cds)
  if (!setequal(df$barcode, rownames(u))) {
    stop("The selected Monocle CDS has different cell barcodes: ", path)
  }
  idx <- match(df$barcode, rownames(u))
  difference <- max(abs(u[idx, 1] - df$UMAP_1),
                    abs(u[idx, 2] - df$UMAP_2))
  if (!is.finite(difference) || difference > 1e-5) {
    stop("Cell UMAP coordinates differ between CSV and Monocle CDS: ", path,
         " (max absolute difference=", signif(difference, 5), ")")
  }
  list(edges = segments, marks = marks,
       audit = data.frame(condition = condition_id,
                          branchpoint = bp,
                          nodes = igraph::vcount(graph),
                          edges = igraph::ecount(graph)))
}
g0 <- extract_graph(input_paths[["control_cds"]], control,
                    "all_genes_control", "Full-gene control")
g1 <- extract_graph(input_paths[["masked_cds"]], masked,
                    "Marsh_activation_masked", "Marsh genes excluded")
segments <- rbind(g0$edges, g1$edges)
branch_marks <- rbind(g0$marks, g1$marks)
graph_audit <- rbind(g0$audit, g1$audit)

# Cell labels assigned by 31N, not inferred anew at plotting time.
cell_state <- function(x, condition) {
  x <- x[x$timepoint == "7 dpi", , drop = FALSE]
  x$condition <- condition
  x$fate <- ifelse(x$aligned_biological_arm == "Aligned_Alternative", "Alternative",
    ifelse(x$aligned_biological_arm == "Aligned_Target", "Target",
    ifelse(x$new_arm == "Other_downstream", "Additional downstream",
    ifelse(x$new_arm == "Root_side", "Prebranch/root side", "Unassigned"))))
  if (any(!x$fate %in% names(PAL))) stop("Unknown cell category in 31N output")
  x
}
plot_cells <- rbind(cell_state(control, "Full-gene control"),
                    cell_state(masked, "Marsh genes excluded"))
plot_cells$condition <- factor(plot_cells$condition,
                               levels = c("Full-gene control", "Marsh genes excluded"))
plot_cells$fate <- factor(plot_cells$fate, levels = names(PAL))
# Gray first, then new downstream arm, then alternative/target on top.
z_order <- c("Unassigned", "Prebranch/root side", "Additional downstream",
             "Alternative", "Target")
plot_cells <- plot_cells[order(match(as.character(plot_cells$fate), z_order)), ]
segments$condition <- factor(segments$condition,
                             levels = levels(plot_cells$condition))
branch_marks$condition <- factor(branch_marks$condition,
                                 levels = levels(plot_cells$condition))

# ------------------------- Panel C: per-sample 7-dpi reassignment.
mask7 <- cell_state(masked, "Marsh genes excluded")
mask7 <- mask7[mask7$original_arm %in% names(ARM_LABS), , drop = FALSE]
if (nrow(mask7) != 2541L) stop("Unexpected 7-dpi original-arm cell count")
mask7$original_branch <- unname(ARM_LABS[mask7$original_arm])
mask7$sample <- factor(mask7$sample_id,
                       levels = SAMPLE_IDS, labels = c("Sample 1", "Sample 2"))
row_names <- c("Sample 1 / Alternative", "Sample 1 / Target",
               "Sample 2 / Alternative", "Sample 2 / Target")
mask7$group <- factor(paste(mask7$sample, mask7$original_branch, sep = " / "),
                      levels = row_names)
if (anyNA(mask7$group)) stop("Unexpected original sample or branch labels")
count_table <- as.data.frame(table(group = mask7$group,
                                   fate = factor(mask7$fate, levels = names(PAL))),
                             stringsAsFactors = FALSE)
colnames(count_table)[3] <- "n"
count_table$n <- as.integer(count_table$n)
expected_totals <- c("Sample 1 / Alternative" = 848L,
                     "Sample 1 / Target" = 644L,
                     "Sample 2 / Alternative" = 551L,
                     "Sample 2 / Target" = 498L)
observed_totals <- tapply(count_table$n, count_table$group, sum)
if (!identical(as.integer(observed_totals[names(expected_totals)]),
               unname(expected_totals))) {
  stop("31N recovery totals do not match the original sample/arm populations")
}
count_table$denominator <- unname(expected_totals[as.character(count_table$group)])
count_table$proportion <- count_table$n / count_table$denominator
count_table$label <- ifelse(count_table$proportion >= .10,
  paste0(count_table$n, "\n", sprintf("%.1f%%", 100 * count_table$proportion)), "")
count_table$fate <- factor(count_table$fate, levels = names(PAL))
count_table$group <- factor(count_table$group, levels = rev(row_names))

# Validate retention counts against the independent 31N per-sample audit.
for (i in seq_len(nrow(recovery))) {
  rr <- recovery[i, , drop = FALSE]
  if (rr$condition != "Marsh_activation_masked") next
  x <- mask7[mask7$original_arm == rr$original_arm, , drop = FALSE]
  if (rr$sample_id != "ALL_7dpi") x <- x[x$sample_id == rr$sample_id, , drop = FALSE]
  retained_category <- unname(ARM_LABS[rr$original_arm])
  if (nrow(x) != rr$n_original ||
      sum(x$fate == retained_category) != rr$n_recovered) {
    stop("Panel C count disagreeing with 31N recovery row ", i)
  }
}
all_control <- recovery[recovery$condition == "all_genes_control" &
                        recovery$sample_id != "ALL_7dpi", , drop = FALSE]
if (nrow(all_control) != 4L || any(all_control$recall < 1 - 1e-8)) {
  stop("Full-gene control fails original source-matched recovery audit")
}

# ----------------------------------------------- Unified GLIA styling.
# Supplementary Figure 7 display revision:
#   * Panel letters use the SAME geometry as Supplementary Figure 5.
#   * All panel subtitles are removed.
#   * All reader-facing fonts are enlarged.
#   * Analysis values, graph topology, cell labels, and statistics are unchanged.

font_family <- "Arial"
BASE_FAMILY <- "Arial"
COL_TEXT <- COL$ink

base_theme <- ggplot2::theme_bw(
  base_size = 12.2,
  base_family = font_family
) +
  ggplot2::theme(
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    panel.border = ggplot2::element_rect(
      color = "#D4D4D4",
      fill = NA,
      linewidth = .46
    ),
    strip.background = ggplot2::element_rect(
      fill = "#F0F0F0",
      color = "#D4D4D4",
      linewidth = .40
    ),
    strip.text = ggplot2::element_text(
      face = "bold",
      size = 11.5,
      color = COL$ink
    ),
    plot.title = ggplot2::element_blank(),
    plot.subtitle = ggplot2::element_blank(),
    axis.title = ggplot2::element_text(
      size = 11.8,
      color = COL$ink
    ),
    axis.text = ggplot2::element_text(
      size = 10.8,
      color = COL$ink
    ),
    legend.title = ggplot2::element_text(
      size = 11.2
    ),
    legend.text = ggplot2::element_text(
      size = 10.6
    ),
    plot.tag = ggplot2::element_blank(),
    plot.margin = ggplot2::margin(7, 10, 8, 10)
  )

# EXACT Supplementary Figure 5 panel-letter geometry.
make_panel_header <- function(letter, title) {
  patchwork::wrap_elements(
    full = grid::grobTree(
      grid::textGrob(
        label = letter,
        x = grid::unit(0.022, "npc"),
        y = grid::unit(0.78, "npc"),
        just = c("left", "center"),
        gp = grid::gpar(
          fontfamily = BASE_FAMILY,
          fontface = "bold",
          fontsize = 14.5,
          col = "black"
        )
      ),
      grid::textGrob(
        label = title,
        x = grid::unit(0.50, "npc"),
        y = grid::unit(0.42, "npc"),
        just = c("center", "center"),
        gp = grid::gpar(
          fontfamily = BASE_FAMILY,
          fontface = "bold",
          fontsize = 11.8,
          col = COL_TEXT
        )
      )
    )
  )
}

assemble_panel <- function(body, letter, title) {
  patchwork::wrap_plots(
    make_panel_header(letter, title),
    body,
    ncol = 1,
    heights = c(0.16, 1.00)
  )
}

# ---------------- Panel A
pA_body <- ggplot2::ggplot(
  score_long,
  ggplot2::aes(
    x = branch,
    y = value,
    fill = branch
  )
) +
  ggplot2::geom_violin(
    trim = FALSE,
    width = .80,
    color = COL$ink,
    linewidth = .34,
    alpha = .82
  ) +
  ggplot2::geom_boxplot(
    width = .17,
    outlier.shape = NA,
    color = COL$ink,
    fill = "white",
    linewidth = .44,
    alpha = .95
  ) +
  ggplot2::facet_grid(
    sample ~ signature
  ) +
  ggplot2::scale_fill_manual(
    values = PAL[c("Alternative", "Target")]
  ) +
  ggplot2::labs(
    x = NULL,
    y = "Mean log1p(CP10k) expression"
  ) +
  base_theme +
  ggplot2::theme(
    legend.position = "none",
    axis.text.x = ggplot2::element_text(
      size = 11.0
    ),
    axis.text.y = ggplot2::element_text(
      size = 10.8
    )
  )

panelA <- assemble_panel(
  body = pA_body,
  letter = "A",
  title = "Dissociation-associated activation signatures"
)

# ---------------- Panel B
pB_body <- ggplot2::ggplot() +
  ggplot2::geom_point(
    data = plot_cells,
    ggplot2::aes(
      x = UMAP_1,
      y = UMAP_2,
      colour = fate
    ),
    size = .48,
    alpha = .85,
    stroke = 0
  ) +
  ggplot2::geom_segment(
    data = segments,
    ggplot2::aes(
      x = x,
      y = y,
      xend = xend,
      yend = yend
    ),
    inherit.aes = FALSE,
    linewidth = .62,
    color = COL$ink,
    alpha = .9
  ) +
  ggplot2::geom_point(
    data = branch_marks,
    ggplot2::aes(
      x = x,
      y = y
    ),
    inherit.aes = FALSE,
    shape = 21,
    size = 3.8,
    stroke = .95,
    fill = "white",
    color = COL$ink
  ) +
  ggplot2::geom_text(
    data = branch_marks,
    ggplot2::aes(
      x = x,
      y = y,
      label = label
    ),
    inherit.aes = FALSE,
    vjust = -1.0,
    size = 3.9,
    color = COL$ink,
    family = font_family,
    fontface = "bold"
  ) +
  ggplot2::facet_wrap(
    ~ condition,
    nrow = 1,
    scales = "fixed"
  ) +
  ggplot2::scale_colour_manual(
    values = PAL,
    drop = FALSE
  ) +
  ggplot2::coord_fixed(
    ratio = 1
  ) +
  ggplot2::labs(
    x = "UMAP 1",
    y = "UMAP 2",
    colour = "Trajectory region"
  ) +
  base_theme +
  ggplot2::theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.key.width = grid::unit(8, "mm"),
    legend.spacing.x = grid::unit(2.5, "mm")
  ) +
  ggplot2::guides(
    colour = ggplot2::guide_legend(
      nrow = 1,
      byrow = TRUE,
      override.aes = list(
        size = 3.4,
        alpha = 1
      )
    )
  )

panelB <- assemble_panel(
  body = pB_body,
  letter = "B",
  title = "Trajectory reconstruction with and without Marsh activation genes"
)

# ---------------- Panel C
pC_body <- ggplot2::ggplot(
  count_table,
  ggplot2::aes(
    x = group,
    y = proportion,
    fill = fate
  )
) +
  ggplot2::geom_col(
    width = .63,
    colour = "white",
    linewidth = .40
  ) +
  ggplot2::geom_text(
    ggplot2::aes(
      label = label
    ),
    position = ggplot2::position_stack(
      vjust = .5
    ),
    colour = COL$ink,
    family = font_family,
    lineheight = .95,
    size = 3.9
  ) +
  ggplot2::coord_flip() +
  ggplot2::scale_fill_manual(
    values = PAL,
    drop = FALSE
  ) +
  ggplot2::scale_y_continuous(
    limits = c(0, 1),
    breaks = c(
      0,
      .25,
      .5,
      .75,
      1
    ),
    labels = function(x) {
      paste0(
        round(x * 100),
        "%"
      )
    },
    expand = ggplot2::expansion(
      mult = c(0, 0)
    )
  ) +
  ggplot2::labs(
    x = NULL,
    y = "Cells within original sample/branch"
  ) +
  base_theme +
  ggplot2::theme(
    legend.position = "none",
    axis.text.y = ggplot2::element_text(
      size = 11.4
    ),
    axis.text.x = ggplot2::element_text(
      size = 10.8
    )
  )

panelC <- assemble_panel(
  body = pC_body,
  letter = "C",
  title = "Reassignment of original branch cells after Marsh gene exclusion"
)

combined <- patchwork::wrap_plots(
  patchwork::wrap_elements(
    full = panelA,
    clip = FALSE
  ),
  patchwork::wrap_elements(
    full = panelB,
    clip = FALSE
  ),
  patchwork::wrap_elements(
    full = panelC,
    clip = FALSE
  ),
  ncol = 1,
  heights = c(
    1.28,
    1.22,
    0.78
  )
) +
  patchwork::plot_annotation(
    theme = ggplot2::theme(
      plot.background = ggplot2::element_rect(
        fill = "white",
        color = NA
      )
    )
  )

# ----------------------------------------------- Export in S7 naming scheme.
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
pdf_out <- file.path(out_dir, paste0(base_name, ".pdf"))
svg_out <- file.path(out_dir, paste0(base_name, ".svg"))
tiff_out <- file.path(out_dir, paste0(base_name, ".tiff"))
png_out <- file.path(out_dir, paste0(base_name, "_review.png"))

# Capture (without suppressing) device/font warnings to make a large warning
# batch actionable rather than treating every warning as a data failure.
export_warnings <- character()
observe_export <- function(expr, label) {
  withCallingHandlers(expr, warning = function(w) {
    export_warnings <<- c(export_warnings, paste0(label, ": ", conditionMessage(w)))
  })
}
pdf_device <- function(filename, width, height, ...) {
  grDevices::cairo_pdf(filename, width = width, height = height, ...)
}
# Cairo embeds Arial where supported; use plain PDF fallback if unavailable.
observe_export(tryCatch(
  ggplot2::ggsave(pdf_out, combined, width = 17, height = 16.8,
                  units = "in", device = pdf_device, bg = "white"),
  error = function(e) {
    warning("Cairo PDF failed; retrying with default PDF device: ", conditionMessage(e))
    ggplot2::ggsave(pdf_out, combined, width = 17, height = 16.8,
                    units = "in", device = "pdf", bg = "white")
  }), "PDF")
observe_export(tryCatch(
  ggplot2::ggsave(svg_out, combined, width = 17, height = 16.8,
                  units = "in", device = "svg", bg = "white"),
  error = function(e) warning("SVG export unavailable: ", conditionMessage(e))), "SVG")
observe_export(tryCatch(
  ggplot2::ggsave(tiff_out, combined, width = 17, height = 16.8,
                  units = "in", dpi = 300, compression = "lzw",
                  device = "tiff", bg = "white"),
  error = function(e) warning("TIFF export unavailable: ", conditionMessage(e))), "TIFF")
observe_export(ggplot2::ggsave(png_out, combined, width = 17, height = 16.8,
                units = "in", dpi = 170, device = "png", bg = "white"), "PNG")
if (length(export_warnings)) {
  wt <- sort(table(export_warnings), decreasing = TRUE)
  writeLines(c(
    paste("Figure S7 v2: export warnings, total =", length(export_warnings)),
    "Warnings are not suppressed in the R console; grouped here for diagnosis.",
    paste(as.integer(wt), names(wt), sep = " x ")
  ), file.path(out_dir, "02_FigureS7_v3_export_warnings.txt"), useBytes = TRUE)
} else {
  writeLines("No warnings occurred during figure export.",
    file.path(out_dir, "02_FigureS7_v3_export_warnings.txt"), useBytes = TRUE)
}

# Publication-oriented audit (not supplementary tables S14-S17).
mask_extra_n <- sum(mask7$fate == "Additional downstream")
mask_alt_to_third <- sum(mask7$original_arm == "Arm_1" &
                         mask7$fate == "Additional downstream")
mask_target_to_alt <- sum(mask7$original_arm == "Arm_2" &
                          mask7$fate == "Alternative")
mask_alt_keep <- sum(mask7$original_arm == "Arm_1" & mask7$fate == "Alternative")
mask_target_keep <- sum(mask7$original_arm == "Arm_2" & mask7$fate == "Target")
if (mask_alt_keep != 1001L || mask_target_keep != 927L ||
    mask_alt_to_third != 341L || mask_target_to_alt != 174L) {
  stop("Correspondence numbers diverged from the approved 31N final report")
}

audit <- data.frame(
  check = c("figure_A_original_branch_cells", "figure_A_two_biological_samples",
    "figure_A_published_activation_sets_scored", "figure_B_trajectory_cell_count",
    "figure_B_fullgene_vertices", "figure_B_fullgene_edges",
    "figure_B_masked_vertices", "figure_B_masked_edges",
    "figure_C_original_alt_cells", "figure_C_original_target_cells",
    "figure_C_original_alt_retained", "figure_C_original_target_retained",
    "figure_C_alt_to_new_third_arm", "figure_C_target_to_new_alt_arm",
    "figure_C_all_original_arm_cells", "figure_B_masked_third_arm_separately_colored",
    "only_two_biological_samples", "no_new_statistical_tests"),
  value = c(nrow(sa), length(unique(sa$sample_id)), length(SIGNATURE_LABELS),
    nrow(control), g0$audit$nodes, g0$audit$edges, g1$audit$nodes,
    g1$audit$edges, 1399, 1142, mask_alt_keep, mask_target_keep,
    mask_alt_to_third, mask_target_to_alt, nrow(mask7), TRUE, TRUE, TRUE))
utils::write.csv(audit, file.path(out_dir, "00_FigureS7_v3_data_audit.csv"),
                 row.names = FALSE)
manifest <- data.frame(source = names(input_paths),
  project_relative_path = substring(
    normalizePath(input_paths, winslash = "/", mustWork = TRUE),
    nchar(gsub("\\\\", "/", root_dir)) + 2),
  md5 = unname(tools::md5sum(input_paths)), stringsAsFactors = FALSE)
utils::write.csv(manifest, file.path(out_dir, "01_FigureS7_v3_input_manifest.csv"),
                 row.names = FALSE)

figure_legend <- paste(
  "Figure S7. Dissociation-associated activation and trajectory sensitivity to exclusion of Marsh activation genes.",
  "(A) Cell-level distributions of the published All CNS Cells Activation and Micro/Myeloid Shared Activation signatures of Marsh et al. (2022), and their nonredundant union, in the original target and alternative branches at 7 dpi. Two biological samples are shown separately. Scores are mean log1p-normalized counts per 10,000 UMIs across the detected genes in each signature. Violins represent cell distributions and inset boxplots indicate median and interquartile range; cells are not independent biological replicates.",
  "(B) Matched full-gene and Marsh-gene-excluded trajectories reconstructed using the original source-matched workflow. Both conditions contained the same 14,005 trajectory cells and used the same 225 root cells. A total of 26 detectable genes from the published 34-gene activation union were excluded before normalization and re-inference. Only 7-dpi cells are colored for clarity, with saved Monocle3 principal graphs overlaid. Objective graph branchpoints are marked BP. The additional downstream branch in the excluded-gene condition is colored separately from prebranch/root-side and unassigned cells. Internal vertex identities are not assumed to correspond between graphs.",
  "(C) Proportions of original 7-dpi target- and alternative-branch cells reassigned to each graph region after gene exclusion, calculated independently within each original branch and biological sample. Values are descriptive. Across the two samples, 1,001/1,399 alternative and 927/1,142 target cells retained corresponding assignments; 341 alternative cells entered an additional downstream branch and 174 target cells were reassigned to the newly defined alternative branch.",
  "Colors follow the main Figures 2/3 visual system: red for target-program, blue for alternative, orange for the new downstream branch, and neutral gray for prebranch/root-side or unassigned cells. These analyses do not establish that either branch is exclusively a dissociation artifact or that the original binary topology is preserved.",
  sep = " "
)
writeLines(figure_legend,
           file.path(out_dir, "FigureS7_GLIA_FINAL_v3_S5_letters_large_font_no_subtitles_legend.txt"), useBytes = TRUE)
writeLines(c(
  "Figure S7 | GLIA | v3 | S5-format panel letters | enlarged fonts | no subtitles | Reviewer 2 Major Comment 2",
  paste("Run date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste("Project root:", root_dir),
  paste("Input Step6F:", step6f),
  paste("Input Step31N:", step31n),
  paste("Output:", out_dir),
  "Final palette source: Figure 2/3 main scripts (Target #B64B45, blue #4C78A8, orange #F28E2B).",
  paste("Full / masked graph BP:", g0$audit$branchpoint, "/", g1$audit$branchpoint),
  paste("Original alternative / target retained:", mask_alt_keep, "/", mask_target_keep),
  paste("Original alternative -> third downstream:", mask_alt_to_third),
  paste("Original target -> new alternative:", mask_target_to_alt),
  "Panel B: fixed facet scales + coord_fixed(ratio=1), avoids free-scale/fixed-aspect error.",
  paste("Export warnings:", length(export_warnings),
        "(see 02_FigureS7_v3_export_warnings.txt)"),
  "No dataset, graph, root, trajectory, pseudotime, or per-cell score recomputed.",
  paste("Publication files:", paste(basename(c(pdf_out, svg_out, tiff_out, png_out)), collapse = "; "))
), file.path(out_dir, "FigureS7_GLIA_FINAL_v3_S5_letters_large_font_no_subtitles_run_summary.txt"), useBytes = TRUE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))

message("Figure S7 output directory: ", out_dir)
message("Publication PDF: ", pdf_out)
message("Publication TIFF: ", tiff_out)
message("Review PNG: ", png_out)
