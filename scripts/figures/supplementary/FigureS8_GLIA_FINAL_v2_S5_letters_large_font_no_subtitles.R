# ==============================================================================
# FigureS8_GLIA_FINAL_v2_S5_letters_large_font_no_subtitles.R
# GLIA | Supplementary Figure S8 | Reviewer 2, Major Comment 3
#
# Project convention (same as Figure S7 and Figure S8):
#   Script:  <PROJECT_ROOT>/script/FigureS8_GLIA_FINAL_v2_S5_letters_large_font_no_subtitles.R
#   Inputs:  xjz/04_results/step11B_.../RDS/ and tables/
#            xjz/09_final_figures/GLIA/Figure4/single_panels/Figure4F/
#   Outputs: xjz/09_final_figures/GLIA/FigureS8_FINAL_v2_S5_letters_large_font_no_subtitles/
#
# Panel A: Joint nucleus-level distributions of the two frozen signatures, by time.
# Panel B: Marginal nucleus-level scores; library medians overlaid to show replication.
# Panel C: Exploratory, injured-library-only endothelial associations.
#
# Scientific safeguards:
#   * Reads already-scored Step11B Seurat metadata; does NOT recompute AUCell,
#     pseudobulk ssGSEA, Harmony, UMAP or any cell state/trajectory.
#   * Target = Arm_2_specific; Alternative = Arm_1_specific.
#   * Both gene sets have different sizes; their absolute AUCell values MUST NOT
#     be directly compared or used to define an arbitrary high/low population.
#   * Density plots are descriptive, not evidence of discrete branches.
#   * Correlation unit = biological library (n = 11 injured libraries), NOT nuclei.
#   * Injury-only restriction does NOT eliminate timepoint confounding.
#   * All existing analysis results and main-text figures are READ ONLY.
# ==============================================================================

options(stringsAsFactors = FALSE)
ROOT <- Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
STEP11B <- file.path(ROOT, "04_results",
  "step11B_GSE234774_refined_main_Arm_specific_shared_AUCell_ssGSEA")
SCORED_RDS <- file.path(STEP11B, "RDS",
  "GSE234774_step11B_refined_main_Arm_scored_BPCells.rds")
LIBRARY_SCORES <- file.path(STEP11B, "tables", "06_refined_library_scores_wide.csv")
FIGURE4F <- file.path(ROOT, "09_final_figures", "GLIA", "Figure4",
  "single_panels", "Figure4F", "Figure4F_v8_library_level_summary.csv")
OUTDIR <- file.path(ROOT, "09_final_figures", "GLIA", "FigureS8_FINAL_v2_S5_letters_large_font_no_subtitles")

inputs <- c(Step11B_scored_object = SCORED_RDS,
            Step11B_library_scores = LIBRARY_SCORES,
            Figure4F_library_summary = FIGURE4F)
missing <- inputs[!file.exists(inputs)]
if (length(missing)) {
  stop("Figure S8 requires these existing files:\n",
       paste(names(missing), missing, sep = ": ", collapse = "\n"),
       "\nCheck the frozen xjz output folders; do not regenerate the analysis just to plot.")
}

required_packages <- c("ggplot2", "patchwork", "SeuratObject")
not_installed <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(not_installed)) {
  stop("Install the required package(s): ", paste(not_installed, collapse = ", "))
}
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(grid)
})
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

check_cols <- function(df, cols, label) {
  absent <- setdiff(cols, names(df))
  if (length(absent)) stop(label, " missing columns: ", paste(absent, collapse = ", "))
}
read_csv_frozen <- function(path) read.csv(path, check.names = FALSE,
                                           stringsAsFactors = FALSE)

# ------------------------------------------------------------------------------
# 1. Read frozen scores and perform exact sample- and library-level crosschecks.
# ------------------------------------------------------------------------------
message("Reading frozen Step11B scored Seurat metadata (no expression-matrix access)...")
scored <- readRDS(SCORED_RDS)
if (!inherits(scored, "Seurat")) stop("Step11B saved object is not a Seurat object.")
needed_metadata <- c("biological_library", "timepoint",
                     "Arm_1_specific_AUCell", "Arm_2_specific_AUCell")
check_cols(scored@meta.data, needed_metadata, "Step11B object metadata")
meta <- as.data.frame(scored@meta.data[, needed_metadata, drop = FALSE])
meta$nucleus_barcode <- rownames(scored@meta.data)
rm(scored)
invisible(gc())

meta$biological_library <- as.character(meta$biological_library)
meta$timepoint <- as.character(meta$timepoint)
meta$alternative <- as.numeric(meta$Arm_1_specific_AUCell)
meta$target <- as.numeric(meta$Arm_2_specific_AUCell)
if (anyNA(meta[, c("biological_library", "timepoint", "alternative", "target")]) ||
    any(!is.finite(meta$alternative) | !is.finite(meta$target))) {
  stop("Missing or non-finite library, timepoint, or AUCell scores in frozen metadata.")
}
if (anyDuplicated(meta$nucleus_barcode)) stop("Duplicate nucleus barcodes in Step11B.")

TP <- c("Uninjured", "7 dpi", "14 dpi", "1 month", "2 months")
if (!setequal(unique(meta$timepoint), TP)) {
  stop("Unexpected Step11B timepoint labels: ",
       paste(sort(unique(meta$timepoint)), collapse = ", "))
}
meta$timepoint <- factor(meta$timepoint, levels = TP)
if (nrow(meta) != 19617L || length(unique(meta$biological_library)) != 13L) {
  stop("Unexpected Step11B population version: expected 19,617 nuclei / 13 libraries, got ",
       nrow(meta), " nuclei / ", length(unique(meta$biological_library)), " libraries.")
}

wide <- read_csv_frozen(LIBRARY_SCORES)
check_cols(wide, c("biological_library", "timepoint", "n_cells",
  "Arm_1_specific_AUCell_median", "Arm_2_specific_AUCell_median"), "Step11B library scores")
endo <- read_csv_frozen(FIGURE4F)
check_cols(endo, c("biological_library", "timepoint", "target_program",
  "endothelial_activation", "endothelial_tnf_nfkB"), "Figure4F library summary")
if (nrow(wide) != 13L || nrow(endo) != 13L ||
    anyDuplicated(wide$biological_library) || anyDuplicated(endo$biological_library)) {
  stop("Each library summary input must contain 13 unique biological libraries.")
}
if (!setequal(unique(meta$biological_library), wide$biological_library) ||
    !setequal(wide$biological_library, endo$biological_library)) {
  stop("Library identifiers differ between Step11B nuclei, Step11B wide CSV, and Figure4F.")
}

library_audit <- do.call(rbind, lapply(sort(unique(meta$biological_library)), function(lib) {
  d <- meta[meta$biological_library == lib, , drop = FALSE]
  if (length(unique(d$timepoint)) != 1L) stop("Multiple timepoints for library: ", lib)
  data.frame(biological_library = lib, timepoint = as.character(d$timepoint[1]),
             n_cells = nrow(d),
             target_AUCell_median_from_nuclei = median(d$target),
             alternative_AUCell_median_from_nuclei = median(d$alternative))
}))
wi <- match(library_audit$biological_library, wide$biological_library)
ei <- match(library_audit$biological_library, endo$biological_library)
if (anyNA(wi) || anyNA(ei)) stop("Internal library matching failed.")
library_audit$target_AUCell_median_Step11B <- wide$Arm_2_specific_AUCell_median[wi]
library_audit$alternative_AUCell_median_Step11B <- wide$Arm_1_specific_AUCell_median[wi]
library_audit$target_AUCell_median_Figure4F <- endo$target_program[ei]
library_audit$Figure4F_timepoint <- endo$timepoint[ei]
library_audit$Step11B_timepoint <- wide$timepoint[wi]
library_audit$Step11B_n_cells <- wide$n_cells[wi]
library_audit$target_abs_difference <- abs(
  library_audit$target_AUCell_median_from_nuclei -
    library_audit$target_AUCell_median_Step11B)
library_audit$alternative_abs_difference <- abs(
  library_audit$alternative_AUCell_median_from_nuclei -
    library_audit$alternative_AUCell_median_Step11B)
library_audit$Figure4F_target_abs_difference <- abs(
  library_audit$target_AUCell_median_from_nuclei -
    library_audit$target_AUCell_median_Figure4F)
if (any(library_audit$target_abs_difference > 1e-6) ||
    any(library_audit$alternative_abs_difference > 1e-6) ||
    any(library_audit$Figure4F_target_abs_difference > 1e-6) ||
    any(library_audit$n_cells != library_audit$Step11B_n_cells) ||
    any(library_audit$timepoint != library_audit$Step11B_timepoint) ||
    any(library_audit$timepoint != library_audit$Figure4F_timepoint)) {
  stop("Nucleus-level scores do not match the frozen Step11B/Figure4F library outputs.")
}

# ------------------------------------------------------------------------------
# 2. Injured-only associations: two endpoints, biological-library unit.
# ------------------------------------------------------------------------------
endo$timepoint <- as.character(endo$timepoint)
if (!setequal(unique(endo$timepoint), TP)) stop("Unexpected Figure4F timepoints.")
injured <- endo[endo$timepoint != "Uninjured", , drop = FALSE]
if (nrow(injured) != 11L) stop("Expected 11 injured biological libraries.")
if (any(!is.finite(injured$target_program) |
        !is.finite(injured$endothelial_activation) |
        !is.finite(injured$endothelial_tnf_nfkB))) {
  stop("Non-finite library-level program or endothelial scores.")
}

test_one <- function(y, endpoint_name) {
  v <- suppressWarnings(cor.test(injured$target_program, y,
                                 method = "spearman", exact = FALSE))
  data.frame(endpoint = endpoint_name, n_libraries = length(y),
             spearman_rho = unname(v$estimate), nominal_p = v$p.value,
             adjustment = "None (two-sided nominal P)")
}
correlation_results <- rbind(
  test_one(injured$endothelial_activation, "Overall endothelial activation"),
  test_one(injured$endothelial_tnf_nfkB, "Endothelial TNF/NF-kB response")
)
if (abs(correlation_results$spearman_rho[1] - 0.45454545) > 0.00001 ||
    abs(correlation_results$spearman_rho[2] - 0.8) > 0.00001) {
  stop("Injured-only correlations differ from the previously audited source data.")
}

# ------------------------------------------------------------------------------
# 3. Plotting data: fixed identical axes for joint score maps.
#    Percent-per-bin is normalized within each timepoint; no cell-level tests.
# ------------------------------------------------------------------------------
count_nuclei <- as.data.frame(table(meta$timepoint), stringsAsFactors = FALSE)
names(count_nuclei) <- c("timepoint", "n_nuclei")
count_libraries <- aggregate(biological_library ~ timepoint, data = meta,
  FUN = function(x) length(unique(x)))
names(count_libraries)[2] <- "n_libraries"
counts <- merge(count_nuclei, count_libraries, by = "timepoint", sort = FALSE)
counts <- counts[match(TP, counts$timepoint), , drop = FALSE]
if (sum(counts$n_nuclei) != 19617L || sum(counts$n_libraries) != 13L) {
  stop("Timepoint nucleus/library counts do not reconcile with Step11B.")
}
facet_labels <- setNames(paste0(counts$timepoint, "\n",
                                format(counts$n_nuclei, big.mark = ",", trim = TRUE),
                                " nuclei / ", counts$n_libraries, " libraries"), TP)
meta$timepoint_label <- factor(unname(facet_labels[as.character(meta$timepoint)]),
                               levels = unname(facet_labels[TP]))

nx <- 35L
ny <- 36L
xrange <- range(meta$alternative)
yrange <- range(meta$target)
xx <- seq(xrange[1] - 1e-7, xrange[2] + 1e-7, length.out = nx + 1L)
yy <- seq(yrange[1] - 1e-7, yrange[2] + 1e-7, length.out = ny + 1L)
xb <- cut(meta$alternative, xx, include.lowest = TRUE, labels = FALSE)
yb <- cut(meta$target, yy, include.lowest = TRUE, labels = FALSE)
if (anyNA(xb) || anyNA(yb)) stop("Joint score-bin assignment failed.")
binned <- aggregate(rep(1L, nrow(meta)),
  by = list(timepoint_label = meta$timepoint_label,
            x_bin = xb, y_bin = yb), FUN = sum)
names(binned)[4] <- "n_nuclei_in_bin"
binned$timepoint_label <- factor(binned$timepoint_label,
                                levels = levels(meta$timepoint_label))
bin_total <- setNames(counts$n_nuclei, unname(facet_labels[TP]))
binned$pct_within_timepoint <- 100 * binned$n_nuclei_in_bin /
  unname(bin_total[as.character(binned$timepoint_label)])
binned$x <- (xx[binned$x_bin] + xx[binned$x_bin + 1L]) / 2
binned$y <- (yy[binned$y_bin] + yy[binned$y_bin + 1L]) / 2

# Separate vertical scales are essential: Arm_1 (100 genes) and Arm_2 (496
# genes) produce AUCell scores that are not directly comparable in magnitude.
plot_nucleus <- rbind(
  data.frame(timepoint = meta$timepoint, program = "Target program",
             AUCell = meta$target),
  data.frame(timepoint = meta$timepoint, program = "Alternative program",
             AUCell = meta$alternative)
)
plot_nucleus$program <- factor(plot_nucleus$program,
                               levels = c("Target program", "Alternative program"))
library_medians <- rbind(
  data.frame(timepoint = factor(library_audit$timepoint, levels = TP),
             program = "Target program",
             AUCell = library_audit$target_AUCell_median_from_nuclei),
  data.frame(timepoint = factor(library_audit$timepoint, levels = TP),
             program = "Alternative program",
             AUCell = library_audit$alternative_AUCell_median_from_nuclei)
)
library_medians$program <- factor(library_medians$program,
                                 levels = levels(plot_nucleus$program))

scatter <- rbind(
  data.frame(biological_library = injured$biological_library,
             timepoint = injured$timepoint, target_program = injured$target_program,
             endothelial_score = injured$endothelial_activation,
             endpoint = "Overall endothelial activation"),
  data.frame(biological_library = injured$biological_library,
             timepoint = injured$timepoint, target_program = injured$target_program,
             endothelial_score = injured$endothelial_tnf_nfkB,
             endpoint = "Endothelial TNF/NF-kB response")
)
format_cor_label <- function(rho, p) {
  paste0("Spearman rho = ", sprintf("%.3f", rho),
         "; nominal P = ", formatC(p, digits = 3, format = "f"))
}
correlation_results$facet <- paste0(
  correlation_results$endpoint, "\n",
  mapply(format_cor_label, correlation_results$spearman_rho,
         correlation_results$nominal_p, USE.NAMES = FALSE))
label_map <- setNames(correlation_results$facet, correlation_results$endpoint)
scatter$endpoint_label <- factor(unname(label_map[scatter$endpoint]),
                                 levels = correlation_results$facet)
scatter$timepoint <- factor(scatter$timepoint, levels = TP[-1])

# ------------------------------------------------------------------------------
# 4. GLIA style: display revision for Supplementary Figure 8.
#    - panel letters match Supplementary Figure 5 header geometry
#    - all panel subtitles removed
#    - reader-facing fonts enlarged
# ------------------------------------------------------------------------------

COL_TARGET <- "#B64B45"
COL_ALTERNATIVE <- "#4C78A8"
COL_INK <- "#202020"
PROGRAM_COLORS <- c("Target program" = COL_TARGET,
                    "Alternative program" = COL_ALTERNATIVE)
INJURY_TIME_COLORS <- c("7 dpi" = "#D62728", "14 dpi" = "#F28E2B",
                        "1 month" = "#7A6AA6", "2 months" = "#438C78")

BASE_FAMILY <- "Arial"
COL_TEXT <- COL_INK

base_theme <- theme_bw(base_size = 12.2, base_family = BASE_FAMILY) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "#F0F0F0", colour = "#CFCFCF"),
    strip.text = element_text(face = "bold", size = 11.2, colour = COL_INK),
    axis.text = element_text(colour = COL_INK, size = 10.6),
    axis.title = element_text(size = 11.3, colour = COL_INK),
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    legend.text = element_text(size = 10.5),
    legend.title = element_text(size = 10.8),
    plot.tag = element_blank(),
    plot.margin = margin(7, 12, 7, 12)
  )

# Exact Supplementary Figure 5 panel-letter geometry.
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
          fontsize = 11.5,
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

# Panel A.
pA_body <- ggplot(binned, aes(x, y, fill = pct_within_timepoint)) +
  geom_tile(width = diff(xx)[1] * 0.99, height = diff(yy)[1] * 0.99,
            colour = NA) +
  facet_wrap(~ timepoint_label, nrow = 1, scales = "fixed") +
  scale_fill_gradient(low = "#E7EBEF", high = "#375F7F",
                      name = "Nuclei per bin\n(% within time)") +
  scale_x_continuous(limits = range(xx), expand = expansion(mult = 0.01)) +
  scale_y_continuous(limits = range(yy), expand = expansion(mult = 0.01)) +
  labs(x = "Alternative-program AUCell", y = "Target-program AUCell") +
  base_theme +
  theme(
    panel.spacing.x = grid::unit(3.5, "mm"),
    strip.text = element_text(size = 10.4),
    legend.position = "right"
  )

panelA <- assemble_panel(
  body = pA_body,
  letter = "A",
  title = "Joint target- and alternative-program scores"
)

# Panel B.
pB_body <- ggplot(plot_nucleus, aes(timepoint, AUCell, fill = program)) +
  geom_violin(trim = TRUE, width = 0.76, alpha = 0.75,
              colour = "#575757", linewidth = 0.28, show.legend = FALSE) +
  geom_boxplot(width = 0.13, outlier.shape = NA, fill = "white",
               colour = "#3B3B3B", linewidth = 0.33, show.legend = FALSE) +
  geom_point(data = library_medians,
             mapping = aes(x = timepoint, y = AUCell), inherit.aes = FALSE,
             position = position_jitter(width = 0.065, height = 0, seed = 16),
             size = 2.2, stroke = 0.62, shape = 21, fill = "white",
             colour = COL_INK) +
  facet_grid(program ~ ., scales = "free_y", switch = "y") +
  scale_fill_manual(values = PROGRAM_COLORS) +
  labs(x = NULL, y = "AUCell (program-specific scale)") +
  base_theme +
  theme(
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 90, size = 11.0),
    axis.text.x = element_text(size = 10.8),
    panel.spacing.y = grid::unit(3.8, "mm")
  )

panelB <- assemble_panel(
  body = pB_body,
  letter = "B",
  title = "Time-resolved nucleus-level score distributions"
)

# Panel C.
pC_body <- ggplot(scatter,
             aes(x = target_program, y = endothelial_score, colour = timepoint)) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
              linewidth = 0.62, linetype = "dashed",
              colour = "#8A8A8A", show.legend = FALSE) +
  geom_point(size = 3.2, alpha = 0.92) +
  facet_wrap(~ endpoint_label, nrow = 1, scales = "free_y") +
  scale_colour_manual(values = INJURY_TIME_COLORS, name = "Injury timepoint") +
  labs(x = "Resident-microglial target-program median AUCell",
       y = "Endothelial response score") +
  base_theme +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    strip.text = element_text(size = 10.8),
    panel.spacing.x = grid::unit(7.5, "mm")
  )

panelC <- assemble_panel(
  body = pC_body,
  letter = "C",
  title = "Microglial target-program and endothelial responses in injured libraries"
)

combined <- patchwork::wrap_plots(
  patchwork::wrap_elements(full = panelA, clip = FALSE),
  patchwork::wrap_elements(full = panelB, clip = FALSE),
  patchwork::wrap_elements(full = panelC, clip = FALSE),
  ncol = 1,
  heights = c(1.45, 1.58, 1.24)
) +
  patchwork::plot_annotation(
    theme = theme(
      plot.background = element_rect(fill = "white", color = NA)
    )
  )

# ------------------------------------------------------------------------------
# 5. Export publication-quality outputs and provenance/audit files.
# ------------------------------------------------------------------------------
BASE <- "FigureS8_GLIA_FINAL_v2_S5_letters_large_font_no_subtitles"
pdf_file <- file.path(OUTDIR, paste0(BASE, ".pdf"))
svg_file <- file.path(OUTDIR, paste0(BASE, ".svg"))
tiff_file <- file.path(OUTDIR, paste0(BASE, ".tiff"))
png_file <- file.path(OUTDIR, paste0(BASE, "_review.png"))
WIDTH_IN <- 16.0
HEIGHT_IN <- 14.6
export_warnings <- character(0)
export_plot <- function(path, device, ..., optional = FALSE) {
  result <- tryCatch(withCallingHandlers(
    ggsave(filename = path, plot = combined, width = WIDTH_IN,
           height = HEIGHT_IN, units = "in", device = device, bg = "white", ...),
    warning = function(w) {
      export_warnings <<- c(export_warnings,
        paste(basename(path), conditionMessage(w), sep = ": "))
      invokeRestart("muffleWarning")
    }), error = function(e) e)
  if (inherits(result, "error")) {
    if (optional) {
      export_warnings <<- c(export_warnings,
        paste("Optional export failed:", basename(path), conditionMessage(result)))
      warning("Optional export unavailable: ", basename(path), " - ",
              conditionMessage(result))
      return(invisible(FALSE))
    }
    stop("Required export failed: ", basename(path), ": ", conditionMessage(result))
  }
  invisible(TRUE)
}

message("Exporting Figure S8 PDF, SVG, TIFF, and preview PNG...")
export_plot(pdf_file, "pdf")
export_plot(svg_file, "svg", optional = TRUE)
export_plot(tiff_file, "tiff", dpi = 300, compression = "lzw", optional = TRUE)
export_plot(png_file, "png", dpi = 180)

# 0. Audit and input hashes. Do not expose or serialize expression matrices.
checks <- data.frame(check = c(
  "Step11B scored nuclei", "Step11B unique libraries", "Timepoints",
  "Two signatures present", "No missing per-nucleus AUCell",
  "Step11B/Figure4F library medians agree", "Injured libraries",
  "Correlation unit is biological library", "No cell-level significance tests",
  "No new trajectory or scoring"),
  value = c(nrow(meta), length(unique(meta$biological_library)),
            paste(TP, collapse = "; "), TRUE, TRUE, TRUE,
            nrow(injured), TRUE, TRUE, TRUE), stringsAsFactors = FALSE)
write.csv(checks, file.path(OUTDIR, "00_FigureS8_data_audit.csv"), row.names = FALSE)

project_relative <- function(path) {
  a <- normalizePath(path, winslash = "/", mustWork = TRUE)
  root <- normalizePath(ROOT, winslash = "/", mustWork = TRUE)
  prefix <- paste0(root, "/")
  if (!startsWith(a, prefix)) stop("Input is outside project root: ", a)
  substring(a, nchar(prefix) + 1L)
}
manifest <- data.frame(source = names(inputs),
                       project_relative_path = vapply(inputs, project_relative, character(1)),
                       md5 = unname(tools::md5sum(inputs)),
                       stringsAsFactors = FALSE)
write.csv(manifest, file.path(OUTDIR, "01_FigureS8_input_manifest.csv"),
          row.names = FALSE)
write.csv(correlation_results[, c("endpoint", "n_libraries", "spearman_rho",
                                "nominal_p", "adjustment")],
          file.path(OUTDIR, "02_FigureS8_injured_library_correlations.csv"),
          row.names = FALSE)
write.csv(counts, file.path(OUTDIR, "03_FigureS8_timepoint_sample_sizes.csv"),
          row.names = FALSE)
write.csv(library_audit, file.path(OUTDIR, "04_FigureS8_library_score_audit.csv"),
          row.names = FALSE)
writeLines(if (length(export_warnings)) export_warnings else "No plotting/export warnings.",
           file.path(OUTDIR, "05_FigureS8_export_warnings.txt"))

legend_text <- paste(
  "Figure S8. Time-resolved target- and alternative-program signatures and",
  "injured-library endothelial associations in independent GSE234774 snRNA-seq.",
  "(A) Joint nucleus-level distributions of the frozen Arm_2-specific (target)",
  "and Arm_1-specific (alternative) AUCell scores across five timepoints.",
  "Heatmap intensity represents the percentage of nuclei within each timepoint",
  "in a common score bin; axes and bins are held constant across timepoints.",
  "(B) Marginal nucleus-level AUCell distributions for each signature, shown",
  "on separate vertical scales because the two gene sets contain different",
  "numbers of genes and their absolute AUCell values are not directly comparable.",
  "Violin distributions summarize nuclei, boxes indicate median and interquartile",
  "range, and white symbols denote individual biological-library medians.",
  "(C) Spearman correlations between the library-level median microglial",
  "target-program AUCell score and overall endothelial activation or endothelial",
  "TNF/NF-kappaB response, restricted to 11 injured biological libraries.",
  "Points correspond to independent biological libraries and are colored by",
  "injury timepoint; dashed lines are descriptive linear guides only.",
  "Displayed P values are nominal two-sided values without multiplicity correction.",
  "Nucleus-level plots are descriptive; statistical inference is not based on",
  "individual nuclei. Continued signature detection does not establish a",
  "discrete target-high/alternative-low population or reproduce the original",
  "Monocle3 branching topology. Restricting correlations to injured libraries",
  "does not eliminate residual confounding by injury timepoint.",
  sep = " ")
writeLines(legend_text, file.path(OUTDIR, "FigureS8_GLIA_FINAL_v2_S5_letters_large_font_no_subtitles_legend.txt"),
           useBytes = TRUE)
writeLines(c(
  "Figure S8 | GLIA | v2 | S5-format panel letters | enlarged fonts | no subtitles | Reviewer 2 Major Comment 3",
  paste0("Date: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste0("Project root: ", ROOT),
  paste0("Output folder: ", OUTDIR),
  paste0("Nuclei: ", nrow(meta), "; biological libraries: ", nrow(wide)),
  paste0("Injured-only libraries: ", nrow(injured)),
  paste0("Overall endothelial activation rho/P: ",
         sprintf("%.3f", correlation_results$spearman_rho[1]), " / ",
         signif(correlation_results$nominal_p[1], 4)),
  paste0("Endothelial TNF/NF-kB rho/P: ",
         sprintf("%.3f", correlation_results$spearman_rho[2]), " / ",
         signif(correlation_results$nominal_p[2], 4)),
  "Inputs are read-only, no expression or AUCell recalculation performed.",
  "Panel C P values are nominal; no correction for injury timepoint was performed.",
  "Panel letters follow the Supplementary Figure 5 header format; panel subtitles were removed and all reader-facing fonts were enlarged.",
  paste0("Exports: ", paste(basename(c(pdf_file, svg_file,
                                       tiff_file, png_file)), collapse = "; "))
), file.path(OUTDIR, "FigureS8_GLIA_FINAL_v2_S5_letters_large_font_no_subtitles_run_summary.txt"),
useBytes = TRUE)
writeLines(capture.output(sessionInfo()), file.path(OUTDIR, "sessionInfo.txt"))
message("Figure S8 output directory: ", OUTDIR)
message("Main PDF: ", pdf_file)
message("TIFF: ", tiff_file)
message("Review PNG: ", png_file)
message("Analysis finished; original Step11B/Figure4F results were not modified.")
