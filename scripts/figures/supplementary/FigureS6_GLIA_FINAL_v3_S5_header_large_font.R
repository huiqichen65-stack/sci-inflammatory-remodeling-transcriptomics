# =============================================================================
# FigureS6_GLIA_FINAL_v3_S5_header_large_font.R
# GLIA | Supplementary Figure S6 | Reviewer 2, Major Comment 1
#
# Project convention:
#   Script:  <PROJECT_ROOT>/script/FigureS6_GLIA_FINAL_v3_S5_header_large_font.R
#   Inputs:  xjz/04_results/step2D_GSE162610_target_program_lineage_specificity_audit/tables/
#   Outputs: xjz/09_final_figures/GLIA/FigureS6_FINAL_v3_S5_header_large_font/
#
# Figure purpose: cross-lineage evidence for the identity of the 7-dpi target
# program. No cell-level pseudoreplicate inference and no new statistical tests.
# This script uses only frozen Step 2D results; it does NOT read/write any ZIP.
# =============================================================================

options(stringsAsFactors = FALSE)
root_dir <- Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
step_dir <- file.path(root_dir, "04_results",
  "step2D_GSE162610_target_program_lineage_specificity_audit")
input_dir <- file.path(step_dir, "tables")
out_dir <- file.path(root_dir, "09_final_figures", "GLIA", "FigureS6_FINAL_v3_S5_header_large_font")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

packages <- c("ggplot2", "patchwork")
not_installed <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(not_installed)) {
  stop("Install packages first: install.packages(c(",
       paste(shQuote(not_installed), collapse = ", "), "))")
}
library(ggplot2)
library(patchwork)
library(grid)

input_paths <- c(
  score = file.path(input_dir, "03_target_signature_by_sample_lineage_group.csv"),
  marker = file.path(input_dir, "05_marker_expression_by_sample_lineage_group.csv"),
  auc = file.path(input_dir, "07_lineage_AUCell_by_sample_group.csv"),
  counts = file.path(input_dir, "02_sample_by_lineage_group_counts.csv")
)
missing <- input_paths[!file.exists(input_paths)]
if (length(missing)) {
  stop("Required Step 2D input files are missing:\n", paste(missing, collapse = "\n"))
}
read_step <- function(p) read.csv(p, check.names = FALSE, stringsAsFactors = FALSE)
score <- read_step(input_paths[["score"]])
marker <- read_step(input_paths[["marker"]])
auc <- read_step(input_paths[["auc"]])
counts <- read_step(input_paths[["counts"]])

assert_cols <- function(x, fields, label) {
  absent <- setdiff(fields, names(x))
  if (length(absent)) stop(label, " is missing required columns: ",
                           paste(absent, collapse = ", "))
}
assert_cols(score, c("sample_id","group","n_cells","median_score","q25","q75"), "Score table")
assert_cols(marker, c("sample_id","group","gene","n_cells","pct_detected","mean_logCP10k"), "Marker table")
assert_cols(auc, c("sample_id","group","lineage_signature","n_cells","n_with_score","median_AUCell"), "AUCell table")
assert_cols(counts, c("sample_id","group","n_cells"), "Count table")

sample_ids <- c("7dpi_sample1", "7dpi_sample2")
sample_names <- c("7 dpi: sample 1", "7 dpi: sample 2")
group_ids <- c(
  "Arm_2_target", "Arm_1_alternative", "Excluded_infiltrating_myeloid",
  "Excluded_BAM_like", "Removed_residual_monocyte_like",
  "Neutrophil_author_annotation"
)
group_names <- c(
  "Target", "Alternative", "Excluded infiltrating\nmyeloid",
  "Excluded BAM-like", "Removed residual\nmonocyte-like",
  "Author-annotated\nneutrophils"
)
names(group_names) <- group_ids

group_colors <- c(
  "Target" = "#B12B34", "Alternative" = "#356C9B",
  "Excluded infiltrating\nmyeloid" = "#7665A4",
  "Excluded BAM-like" = "#438C78",
  "Removed residual\nmonocyte-like" = "#AF7540",
  "Author-annotated\nneutrophils" = "#767676"
)

assert_groups <- function(x, what) {
  if (!setequal(unique(x$sample_id), sample_ids) ||
      !setequal(unique(x$group), group_ids)) {
    stop(what, ": unexpected sample IDs or comparison groups. Verify input versions.")
  }
}
assert_groups(score, "Score table")
assert_groups(marker, "Marker table")
assert_groups(auc, "AUCell table")
assert_groups(counts, "Count table")

# Freeze display order across all panels and data audits.
prep_group <- function(df) {
  df$group_display <- factor(unname(group_names[df$group]), levels = unname(group_names))
  df$sample_display <- factor(df$sample_id, levels = sample_ids, labels = sample_names)
  df
}
score <- prep_group(score)
marker <- prep_group(marker)
auc <- prep_group(auc)
counts <- prep_group(counts)

if (nrow(score) != 12L || anyDuplicated(score[c("sample_id", "group")]))
  stop("Score table must contain exactly 12 unique sample-group rows.")
if (any(!is.finite(score$median_score) | !is.finite(score$q25) |
        !is.finite(score$q75) | score$q25 > score$median_score |
        score$median_score > score$q75))
  stop("Score table has invalid or missing median/IQR values.")
if (any(!is.finite(marker$pct_detected) | marker$pct_detected < 0 |
        marker$pct_detected > 1))
  stop("pct_detected must be a fraction between 0 and 1.")
if (anyDuplicated(marker[c("sample_id", "group", "gene")]))
  stop("Marker expression table has duplicate sample/group/gene rows.")
if (anyDuplicated(auc[c("sample_id", "group", "lineage_signature")]))
  stop("AUCell table has duplicate sample/group/signature rows.")

n_check <- merge(score[c("sample_id", "group", "n_cells")],
                 counts[c("sample_id", "group", "n_cells")],
                 by = c("sample_id", "group"), suffixes = c("_score", "_count"))
if (nrow(n_check) != 12L || any(n_check$n_cells_score != n_check$n_cells_count))
  stop("Counts do not match between signature and group-count files.")

# All 18 displayed genes are from the Step 2D 31-gene marker analysis.
# They span the precise marker types requested by Reviewer 2 (Major Comment 1).
genes <- c(
  "Tnf", "C5ar1", "C3ar1", "Myd88", "Icam1",                  # inflammatory
  "Hexb", "Sall1", "P2ry12", "Tmem119", "Fcrls",             # resident-like
  "Ms4a7", "Ccr2", "Ly6c2", "S100a8", "S100a9",             # infiltrating
  "Lyve1", "Cd163", "Mrc1"                                   # BAM-like
)
missing_genes <- setdiff(genes, unique(marker$gene))
if (length(missing_genes))
  stop("Marker table is missing required genes: ", paste(missing_genes, collapse = ", "))
marker <- marker[marker$gene %in% genes, , drop = FALSE]
marker$gene <- factor(marker$gene, levels = genes)
if (nrow(marker) != length(genes) * 12L)
  stop("Marker table does not contain all 18 genes in both samples and all 6 groups.")

expected_auc <- c("AUCell_Resident_microglia", "AUCell_Infiltrating_mono_macro", "AUCell_BAM_like")
if (!all(expected_auc %in% unique(auc$lineage_signature)))
  stop("Lineage AUCell table is missing one of the three required signatures.")

# IMPORTANT: the author-annotated neutrophil group overlaps with other groups,
# rather than constituting a disjoint lineage. AUCell coverage is incomplete in
# that group (15/184 and 7/16 cells in the supplied Step2D outputs). Exclude this
# group from Panel B ONLY, not from Panels A and C. Never show its partial AUCell
# subset as if it characterized all author-annotated neutrophils.
auc_neutro <- auc[auc$group == "Neutrophil_author_annotation", , drop = FALSE]
if (nrow(auc_neutro) && any(auc_neutro$n_with_score < auc_neutro$n_cells)) {
  message("Panel B: author-neutrophil AUCell omitted because of incomplete coverage.")
}
auc_plot <- auc[auc$group != "Neutrophil_author_annotation", , drop = FALSE]
auc_plot$group_display <- factor(as.character(auc_plot$group_display),
                                levels = unname(group_names[1:5]))
auc_plot$lineage_display <- factor(auc_plot$lineage_signature,
  levels = expected_auc,
  labels = c("Resident microglia", "Infiltrating mono/macrophage", "BAM-like"))
if (any(auc_plot$n_with_score < auc_plot$n_cells)) {
  stop("AUCell coverage is also incomplete in one of the five Panel B groups; review data.")
}

# =============================================================================
# S5-matched panel-letter geometry + globally enlarged typography
# =============================================================================

BASE_FAMILY <- "Arial"
COL_TEXT <- "#202020"
DPI_OUT <- 600

# Keep the original biological-group colours.
# Typography is intentionally larger than the original S7 script.
theme_s6_large <- function(base_size = 10.8) {
  theme_classic(
    base_family = BASE_FAMILY,
    base_size = base_size
  ) +
    theme(
      text = element_text(
        family = BASE_FAMILY,
        colour = COL_TEXT
      ),
      axis.title = element_text(
        family = BASE_FAMILY,
        size = 11.2,
        colour = COL_TEXT
      ),
      axis.text = element_text(
        family = BASE_FAMILY,
        size = 9.6,
        colour = COL_TEXT
      ),
      axis.line = element_line(
        linewidth = 0.55,
        colour = "#333333"
      ),
      axis.ticks = element_line(
        linewidth = 0.50,
        colour = "#333333"
      ),
      strip.background = element_blank(),
      strip.text = element_text(
        family = BASE_FAMILY,
        face = "bold",
        size = 10.8,
        colour = COL_TEXT,
        margin = margin(b = 3)
      ),
      legend.title = element_text(
        family = BASE_FAMILY,
        size = 10.2,
        colour = COL_TEXT
      ),
      legend.text = element_text(
        family = BASE_FAMILY,
        size = 9.5,
        colour = COL_TEXT
      ),
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      plot.tag = element_blank(),
      plot.margin = margin(
        3.2, 4.0, 3.6, 4.0,
        unit = "mm"
      )
    )
}

# -------------------------------------------------------------------------
# EXACT Supplementary Figure 5 panel-header geometry
# The panel letter position / alignment / size matches Figure S5.
# -------------------------------------------------------------------------
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
          fontsize = 10.2,
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
    heights = c(0.18, 1.00)
  )
}

save_plot_multi <- function(
  plot,
  stem,
  width,
  height,
  dpi = DPI_OUT
) {
  pdf_device <- if (capabilities("cairo")) cairo_pdf else grDevices::pdf

  ggsave(
    paste0(stem, ".pdf"),
    plot,
    width = width,
    height = height,
    units = "in",
    device = pdf_device,
    bg = "white"
  )

  ggsave(
    paste0(stem, "_review.png"),
    plot,
    width = width,
    height = height,
    units = "in",
    dpi = 220,
    bg = "white"
  )

  if (requireNamespace("svglite", quietly = TRUE)) {
    ggsave(
      paste0(stem, ".svg"),
      plot,
      width = width,
      height = height,
      units = "in",
      bg = "white"
    )
  }

  if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_tiff(
      filename = paste0(stem, ".tiff"),
      width = width,
      height = height,
      units = "in",
      res = dpi,
      compression = "lzw",
      background = "white"
    )
    print(plot)
    grDevices::dev.off()
  } else {
    grDevices::tiff(
      filename = paste0(stem, ".tiff"),
      width = width,
      height = height,
      units = "in",
      res = dpi,
      compression = "lzw",
      bg = "white"
    )
    print(plot)
    grDevices::dev.off()
  }
}

# --- Panel A: frozen target-program signature, descriptive median + IQR ---
score$n_label <- paste0("n=", score$n_cells)

pA_body <- ggplot(
  score,
  aes(
    x = group_display,
    y = median_score,
    color = group_display
  )
) +
  geom_linerange(
    aes(ymin = q25, ymax = q75),
    linewidth = 1.25,
    show.legend = FALSE
  ) +
  geom_point(
    size = 3.0,
    show.legend = FALSE
  ) +
  geom_text(
    aes(label = n_label),
    color = "#303030",
    family = BASE_FAMILY,
    size = 3.0,
    vjust = -0.90,
    show.legend = FALSE
  ) +
  facet_wrap(
    ~ sample_display,
    nrow = 1
  ) +
  scale_color_manual(
    values = group_colors,
    drop = FALSE
  ) +
  scale_y_continuous(
    expand = expansion(
      mult = c(0.10, 0.23)
    )
  ) +
  labs(
    y = "Median mean log1p(CP10k) (IQR)",
    x = NULL
  ) +
  theme_s6_large() +
  theme(
    axis.text.x = element_text(
      angle = 34,
      hjust = 1,
      vjust = 1,
      size = 9.1
    ),
    axis.text.y = element_text(
      size = 9.8
    ),
    panel.spacing.x = unit(4.5, "mm")
  )

panelA <- assemble_panel(
  body = pA_body,
  letter = "A",
  title = "Frozen 496-gene target-program score"
)

# --- Panel B: lineage AUCell scores ---
pB_body <- ggplot(
  auc_plot,
  aes(
    x = group_display,
    y = lineage_display,
    fill = median_AUCell
  )
) +
  geom_tile(
    color = "white",
    linewidth = 0.85
  ) +
  geom_text(
    aes(label = sprintf("%.2f", median_AUCell)),
    family = BASE_FAMILY,
    size = 3.25,
    colour = COL_TEXT
  ) +
  facet_wrap(
    ~ sample_display,
    nrow = 1
  ) +
  scale_fill_gradient(
    low = "#F3F6FA",
    high = "#3E78A4",
    name = "Median\nAUCell"
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_s6_large() +
  theme(
    axis.text.x = element_text(
      angle = 36,
      hjust = 1,
      vjust = 1,
      size = 8.8
    ),
    axis.text.y = element_text(
      size = 9.7
    ),
    legend.position = "top",
    legend.justification = "center",
    legend.direction = "horizontal",
    legend.key.width = unit(8, "mm"),
    legend.key.height = unit(4, "mm"),
    panel.spacing.x = unit(5.0, "mm")
  )

panelB <- assemble_panel(
  body = pB_body,
  letter = "B",
  title = "Lineage-program scores"
)

# --- Panel C: marker expression ---
marker$group_rev <- factor(
  as.character(marker$group_display),
  levels = rev(unname(group_names))
)

pC_body <- ggplot(
  marker,
  aes(
    x = gene,
    y = group_rev
  )
) +
  geom_point(
    aes(
      size = pct_detected,
      fill = mean_logCP10k
    ),
    shape = 21,
    stroke = 0.32,
    color = "#303030"
  ) +
  facet_wrap(
    ~ sample_display,
    nrow = 1
  ) +
  scale_size_continuous(
    range = c(0.35, 6.6),
    breaks = c(
      0.0,
      0.25,
      0.5,
      0.75,
      1.0
    ),
    labels = function(x) {
      paste0(round(x * 100), "%")
    },
    name = "Cells detected"
  ) +
  scale_fill_gradient(
    low = "#F2F5F9",
    high = "#A22A35",
    name = "Mean\nlogCP10k"
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_s6_large() +
  theme(
    axis.text.x = element_text(
      angle = 50,
      hjust = 1,
      vjust = 1,
      size = 9.4
    ),
    axis.text.y = element_text(
      size = 10.0
    ),
    legend.position = "right",
    legend.box = "vertical",
    legend.spacing.y = unit(2.5, "mm"),
    panel.spacing.x = unit(5.0, "mm")
  )

panelC <- assemble_panel(
  body = pC_body,
  letter = "C",
  title = "Inflammatory and lineage-associated markers"
)

# Keep the original S7 layout:
# A + B on the top row; C spans the full lower row.
row1 <- patchwork::wrap_plots(
  patchwork::wrap_elements(
    full = panelA,
    clip = FALSE
  ),
  patchwork::wrap_elements(
    full = panelB,
    clip = FALSE
  ),
  ncol = 2,
  widths = c(1.02, 0.98)
)

row2 <- patchwork::wrap_plots(
  patchwork::wrap_elements(
    full = panelC,
    clip = FALSE
  ),
  ncol = 1
)

combined <- patchwork::wrap_plots(
  patchwork::wrap_elements(
    full = row1,
    clip = FALSE
  ),
  patchwork::wrap_elements(
    full = row2,
    clip = FALSE
  ),
  ncol = 1,
  heights = c(
    1.02,
    1.08
  )
) +
  patchwork::plot_annotation(
    theme = theme(
      plot.background = element_rect(
        fill = "white",
        color = NA
      )
    )
  )

base_name <- "FigureS6_GLIA_FINAL_v3_S5_header_large_font"

pdf_out  <- file.path(out_dir, paste0(base_name, ".pdf"))
svg_out  <- file.path(out_dir, paste0(base_name, ".svg"))
tiff_out <- file.path(out_dir, paste0(base_name, ".tiff"))
png_out  <- file.path(out_dir, paste0(base_name, "_review.png"))

# Larger canvas prevents clipping after globally enlarging typography.
save_plot_multi(
  panelA,
  file.path(
    out_dir,
    "FigureS6A_target_program_score_v3_large_font"
  ),
  width = 8.4,
  height = 6.3
)

save_plot_multi(
  panelB,
  file.path(
    out_dir,
    "FigureS6B_lineage_program_scores_v3_large_font"
  ),
  width = 8.4,
  height = 6.3
)

save_plot_multi(
  panelC,
  file.path(
    out_dir,
    "FigureS6C_lineage_marker_expression_v3_large_font"
  ),
  width = 10.0,
  height = 6.8
)

save_plot_multi(
  combined,
  file.path(
    out_dir,
    base_name
  ),
  width = 18.0,
  height = 13.0
)

# Match the existing xjz FigureS1-S6 figure audit naming convention.
audit <- data.frame(
  check = c("source_files_present", "sample_groups_in_all_inputs", "signature_gene_count",
            "n_unique_groups_per_sample", "marker_genes_displayed",
            "n_target_7dpi_sample1", "n_target_7dpi_sample2",
            "neutrophil_AUCell_excluded_due_to_incomplete_coverage",
            "two_biological_samples_only", "no_significance_testing"),
  result = c(TRUE, TRUE, 496, 6, length(genes),
             score$n_cells[score$sample_id == sample_ids[1] & score$group == group_ids[1]],
             score$n_cells[score$sample_id == sample_ids[2] & score$group == group_ids[1]],
             TRUE, TRUE, TRUE)
)
write.csv(audit, file.path(out_dir, "00_FigureS6_data_audit.csv"), row.names = FALSE)
manifest <- data.frame(source = names(input_paths),
                       project_relative_path = substring(normalizePath(input_paths,
                         winslash = "/", mustWork = TRUE), nchar(gsub("\\\\", "/", root_dir)) + 2),
                       md5 = unname(tools::md5sum(input_paths)),
                       stringsAsFactors = FALSE)
write.csv(manifest, file.path(out_dir, "01_FigureS6_input_manifest.csv"), row.names = FALSE)

figure_legend <- paste(
  "Figure S6. Cross-lineage evaluation of the 7-dpi target-program signature and myeloid identity.",
  "(A) Median and interquartile range of the frozen 496-gene Arm_2-specific program score",
  "(mean log1p counts per 10,000 UMIs), shown separately for the two 7-dpi biological samples.",
  "Numbers indicate cells per sample/group. (B) Median AUCell scores for resident microglial,",
  "infiltrating monocyte/macrophage, and BAM-like lineage signatures. The author-annotated",
  "neutrophil comparator is omitted from this panel because original AUCell scores were",
  "available for only 15/184 and 7/16 cells in the two samples, respectively; it remains",
  "included in panels A and C. (C) Expression of representative inflammatory, resident",
  "microglial, infiltrating myeloid, and BAM-associated markers. Dot size indicates the",
  "fraction of cells with nonzero raw counts; color denotes mean log1p-normalized expression",
  "per 10,000 UMIs. Author-annotated neutrophils constitute an overlapping comparator",
  "rather than an exclusive lineage group. All comparisons are descriptive; only two",
  "biological samples were available at 7 dpi. Target-program score differences are not",
  "independent confirmation of lineage identity because the signature was selected from",
  "the discovery branch; transcriptomic profiles do not establish developmental origin.",
  sep = " "
)
writeLines(figure_legend, file.path(out_dir, "FigureS6_GLIA_FINAL_v3_S5_header_large_font_legend.txt"))
writeLines(c(
  "Figure S6 | GLIA | v3 | S5-format panel letters + enlarged typography | Reviewer 2 Major Comment 1",
  paste0("Run date: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste0("Project root: ", root_dir),
  paste0("Input folder: ", input_dir),
  paste0("Output folder: ", out_dir),
  paste0("Target branch: n=", audit$result[6], " / n=", audit$result[7], " at the two 7-dpi samples"),
  "Panels A and C have six groups, panel B has five groups owing to partial neutrophil AUCell coverage.",
  "No statistical tests or new expression estimates were calculated; this is a figure-only assembly.",
  "Panel letters use the same header geometry as Supplementary Figure 5; all reader-facing typography was enlarged.",
  paste0("Output files: ", paste(basename(c(pdf_out, svg_out, tiff_out, png_out)), collapse = "; "))
), file.path(out_dir, "FigureS6_GLIA_FINAL_v3_S5_header_large_font_run_summary.txt"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))

message("Figure S6 output directory: ", out_dir)
message("Main PDF: ", pdf_out)
message("Main TIFF: ", tiff_out)
message("Review PNG: ", png_out)
