# ==============================================================================
# Figure S9 | GLIA | Major Comment 5 | FINAL v4 aligned large font
# Endothelial receptor-detection sensitivity and the cellular Tnf source at 7 dpi
#
# One supplementary figure with FOUR panels:
#   A  Tnfrsf1a detection and depth imbalance in 7-dpi endothelial cells
#   B  Original endothelial/TNF-NF-kB program effects before/after depth correction
#      (point estimates only; no de-overlapped group contrasts)
#   C  Original versus nonoverlapping endothelial-program correlation
#   D  Tnf detection and average expression across ALL 7-dpi populations,
#      including refined myeloid subgroups
#
# DIRECT ORIGINAL-RESULT MODE (same principle as Figure S1 v8):
#   - Reads ONLY the previously completed xjz/04_results/step8G/... outputs.
#   - Does NOT read any evidence ZIP, upload bundle, or temporary payload.
#   - Does NOT rerun normalization, scoring, gene set selection or statistics.
#   - Reuses the exact numerical results produced by Step 8G.
#
# Inputs (tables in Step 8G):
#   03_original_step8C_reproduction_audit.csv
#   04_receptor_depth_adjusted_within_sample_models.csv
#   05_receptor_positive_negative_depth_audit.csv
#   07_program_pair_correlations_before_after_deoverlap.csv
#   10_Tnf_refined_myeloid_and_other_groups_by_sample.csv
#   12_STEP8G_audit_checks.csv
#
# Output:
#   xjz/09_final_figures/GLIA/FigureS9_FINAL_v4_aligned_large_font/
#     FigureS9A_receptor_detection_and_depth.{pdf,svg,tiff}
#     FigureS9B_endothelial_effect_sensitivity.{pdf,svg,tiff}
#     FigureS9C_gene_overlap_correlation.{pdf,svg,tiff}
#     FigureS9D_Tnf_cellular_sources.{pdf,svg,tiff}
#     FigureS9_GLIA_FINAL_v4_aligned_large_font.{pdf,svg,tiff}
#     FigureS9_GLIA_FINAL_v4_aligned_large_font_review.png
#     00_FigureS9_data_audit.csv
#     01_FigureS9_input_manifest.csv
#     02_FigureS9_*_plot_data.csv (frozen plotting inputs)
#     FigureS9_GLIA_FINAL_v4_aligned_large_font_run_summary.txt
#     Table_S21_revised.xlsx and Table_S21_revised.csv
#     sessionInfo.txt
#
# Figure style follows the supplied FigureS1_GLIA_FINAL_v8_FINAL.R:
#   Arial; compact academic styling; consistent restrained red/grey/blue;
#   patchwork with wrapped nested layouts; PDF/SVG/600-dpi TIFF/PNG preview.
# ==============================================================================

options(stringsAsFactors = FALSE)
set.seed(20261026)

# ------------------------------------------------------------------------------
# 0. Direct original project paths
# ------------------------------------------------------------------------------
ROOT <- Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
if (!dir.exists(ROOT)) stop("Missing project root: ", ROOT)

STEP8G <- file.path(
  ROOT, "04_results",
  "step8G_GSE162610_receptor_depth_overlap_Tnf_source_sensitivity"
)
OUT <- file.path(ROOT, "09_final_figures", "GLIA", "FigureS9_FINAL_v4_aligned_large_font")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

message("Figure S9 direct-source mode")
message("Input directory:  ", STEP8G)
message("Output directory: ", OUT)

# ------------------------------------------------------------------------------
# 1. Packages (same visualization pipeline as reference Figure S1)
# ------------------------------------------------------------------------------
required_pkgs <- c("ggplot2", "dplyr", "readr", "patchwork", "scales",
                   "svglite", "ragg", "openxlsx")
missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_pkgs)) {
  stop("Missing R packages: ", paste(missing_pkgs, collapse = ", "),
       "\nInstall using: install.packages(c(",
       paste(sprintf('"%s"', missing_pkgs), collapse = ", "), "))")
}
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(patchwork)
  library(scales)
})

# ------------------------------------------------------------------------------
# 2. Exact upstream Step 8G outputs, no re-analysis
# ------------------------------------------------------------------------------
P <- list(
  original_reproduction = file.path(STEP8G, "tables",
    "03_original_step8C_reproduction_audit.csv"),
  models = file.path(STEP8G, "tables",
    "04_receptor_depth_adjusted_within_sample_models.csv"),
  depth = file.path(STEP8G, "tables",
    "05_receptor_positive_negative_depth_audit.csv"),
  correlations = file.path(STEP8G, "tables",
    "07_program_pair_correlations_before_after_deoverlap.csv"),
  tnf_groups = file.path(STEP8G, "tables",
    "10_Tnf_refined_myeloid_and_other_groups_by_sample.csv"),
  step8G_checks = file.path(STEP8G, "tables",
    "12_STEP8G_audit_checks.csv")
)
missing_inputs <- names(P)[!vapply(P, file.exists, logical(1))]
if (length(missing_inputs)) {
  stop("Step 8G original output(s) missing:\n",
       paste(sprintf(" - %s: %s", missing_inputs,
                     unlist(P[missing_inputs])), collapse = "\n"),
       "\nRun the existing Step 8G analysis first; do not use a ZIP instead.")
}

manifest <- data.frame(
  input_key = names(P),
  absolute_path = unname(unlist(P)),
  bytes = as.numeric(file.info(unlist(P))$size),
  md5 = unname(tools::md5sum(unlist(P))),
  stringsAsFactors = FALSE
)
write.csv(manifest, file.path(OUT, "01_FigureS9_input_manifest.csv"),
          row.names = FALSE)

read_one <- function(path) {
  as.data.frame(readr::read_csv(path, show_col_types = FALSE,
                                 progress = FALSE),
                stringsAsFactors = FALSE)
}
reproduction <- read_one(P$original_reproduction)
mod0 <- read_one(P$models)
depth0 <- read_one(P$depth)
rho0 <- read_one(P$correlations)
tnf0 <- read_one(P$tnf_groups)
checks <- read_one(P$step8G_checks)

# ------------------------------------------------------------------------------
# 3. Hard audit before plotting; stop rather than silently drop data
# ------------------------------------------------------------------------------
need_columns <- function(df, columns, key) {
  miss <- setdiff(columns, names(df))
  if (length(miss)) {
    stop(key, " lacks required column(s): ", paste(miss, collapse = ", "))
  }
}
need_columns(reproduction, c("sample_id", "program", "PASS"), "reproduction")
need_columns(mod0, c("sample_id", "program", "gene_set", "covariate_model",
                     "n_cells", "n_pos", "n_neg", "beta_receptor_positive",
                     "se_cell_level_descriptive",
                     "ci95_low_cell_level_descriptive",
                     "ci95_high_cell_level_descriptive",
                     "model_r_squared"), "models")
need_columns(depth0, c("sample_id", "receptor_positive", "n_cells",
                       "fraction_of_sample", "median_nUMI", "median_nFeature"),
             "depth")
need_columns(rho0, c("sample_id", "original_program_rho",
                     "nonoverlap_program_rho", "n_cells"), "correlations")
need_columns(tnf0, c("sample_id", "cell_group", "n_cells", "n_Tnf_detected",
                     "pct_Tnf_detected", "mean_logRNA_all_cells"), "Tnf")
need_columns(checks, c("check", "PASS"), "Step 8G checks")

SAMPLES <- c("7dpi_sample1", "7dpi_sample2")
sample_short <- c("7dpi_sample1" = "7 dpi 1",
                  "7dpi_sample2" = "7 dpi 2")
expected_programs <- c("ENDOTHELIAL_ACTIVATION",
                       "TNF_NFKB_RESPONSE",
                       "ENDOTHELIAL_ACTIVATION_NONOVERLAP",
                       "TNF_NFKB_RESPONSE_NONOVERLAP")
expected_models <- c("unadjusted", "adjust_log_nUMI",
                     "adjust_log_nFeature", "adjust_both")

pair_complete <- function(d, expected_n, cols) {
  nrow(d) == expected_n && !anyDuplicated(d[, cols, drop = FALSE])
}

all_samples <- function(d) setequal(unique(as.character(d$sample_id)), SAMPLES)
valid_bools <- function(v) all(!is.na(v) & as.logical(v))

qc <- data.frame(
  check = c(
    "Upstream Step8G audits all PASS",
    "Original Step8C reproduction all PASS",
    "Exactly 2 intended 7dpi biological samples in every input",
    "32 model rows: 2 samples x 4 programs x 4 covariate models",
    "All four original/nonoverlap gene-program identifiers present",
    "Both baseline and dual-depth-corrected models available",
    "Depth summary has four unique receptor-detection rows",
    "Program correlation summary has two sample rows",
    "Tnf summary has 36 unique sample/group rows",
    "Refined Tnf summary covers all 14505 7dpi cells",
    "Tnf proportions, per-group cell counts, and model CIs valid"
  ),
  PASS = FALSE, details = "", stringsAsFactors = FALSE
)
qc$PASS[1] <- valid_bools(checks$PASS)
qc$PASS[2] <- valid_bools(reproduction$PASS)
qc$PASS[3] <- all(vapply(list(reproduction, mod0, depth0, rho0, tnf0),
                         all_samples, logical(1)))
qc$PASS[4] <- pair_complete(mod0, 32,
  c("sample_id", "program", "covariate_model")) &&
  setequal(unique(mod0$covariate_model), expected_models) &&
  all(mod0$n_cells == mod0$n_pos + mod0$n_neg)
qc$PASS[5] <- setequal(unique(mod0$program), expected_programs) &&
  setequal(unique(mod0$gene_set), c("ORIGINAL", "NONOVERLAP"))
qc$PASS[6] <- all(c("unadjusted", "adjust_both") %in% mod0$covariate_model)
qc$PASS[7] <- pair_complete(depth0, 4,
  c("sample_id", "receptor_positive")) &&
  setequal(as.character(depth0$receptor_positive), c("TRUE", "FALSE"))
qc$PASS[8] <- pair_complete(rho0, 2, "sample_id") &&
  all(abs(rho0$original_program_rho) <= 1) &&
  all(abs(rho0$nonoverlap_program_rho) <= 1)
qc$PASS[9] <- pair_complete(tnf0, 36, c("sample_id", "cell_group")) &&
  all(table(tnf0$sample_id) == 18)
qc$PASS[10] <- sum(tnf0$n_cells) == 14505 &&
  all(tnf0$n_cells > 0) &&
  all(tnf0$n_Tnf_detected <= tnf0$n_cells)
qc$PASS[11] <- all(is.finite(mod0$beta_receptor_positive)) &&
  all(mod0$ci95_low_cell_level_descriptive <=
      mod0$beta_receptor_positive) &&
  all(mod0$beta_receptor_positive <=
      mod0$ci95_high_cell_level_descriptive) &&
  all(is.finite(tnf0$pct_Tnf_detected)) &&
  all(tnf0$pct_Tnf_detected >= 0 & tnf0$pct_Tnf_detected <= 100) &&
  all(is.finite(tnf0$mean_logRNA_all_cells))
qc$details[10] <- paste("n=", sum(tnf0$n_cells))
write.csv(qc, file.path(OUT, "00_FigureS9_data_audit.csv"),
          row.names = FALSE)
if (!all(qc$PASS)) {
  stop("Figure S9 input audit FAILED:\n - ",
       paste(qc$check[!qc$PASS], collapse = "\n - "),
       "\nReview: ", file.path(OUT, "00_FigureS9_data_audit.csv"))
}
message("All 11 Figure S9 audit checks passed.")

# ------------------------------------------------------------------------------
# 4. Frozen style — consistent with the supplied Figure S1
# ------------------------------------------------------------------------------
base_family <- "Arial"
COLOR_S1 <- "#B33F3F"
COLOR_S2 <- "#E45756"
COLOR_ADJUSTED <- "#4C78A8"
COLOR_GREY <- "#A5AAB2"
COLOR_DARK <- "#30343A"
sample_colors <- c("7dpi_sample1" = COLOR_S1,
                   "7dpi_sample2" = COLOR_S2)

pub_theme <- ggplot2::theme_classic(base_family = base_family,
                                    base_size = 9.2) +
  ggplot2::theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    axis.title = element_text(size = 9.0, color = "black"),
    axis.text = element_text(size = 7.8, color = "black"),
    axis.line = element_line(linewidth = 0.42, color = COLOR_DARK),
    axis.ticks = element_line(linewidth = 0.36, color = COLOR_DARK),
    legend.title = element_text(size = 8.0, color = "black"),
    legend.text = element_text(size = 7.4, color = "black"),
    legend.key.height = grid::unit(3.2, "mm"),
    legend.key.width = grid::unit(3.4, "mm"),
    plot.title = element_text(size = 9.2, face = "bold", hjust = 0,
                              margin = margin(b = 2.2, unit = "mm")),
    plot.subtitle = element_blank(),
    plot.tag = element_text(size = 13.0, face = "bold", color = "black"),
    plot.tag.position = c(-0.035, 1.05),
    plot.margin = margin(3.4, 3.6, 3.4, 4.0, unit = "mm")
  )

# Shared top-row header: Panel A and Panel B use the SAME grob geometry,
# guaranteeing that A and B are exactly level in the assembled top row.
make_top_panel_header <- function(letter, title) {
  patchwork::wrap_elements(
    full = grid::grobTree(
      grid::textGrob(
        label = letter,
        x = grid::unit(0.018, "npc"),
        y = grid::unit(0.68, "npc"),
        just = c("left", "center"),
        gp = grid::gpar(
          fontfamily = base_family,
          fontface = "bold",
          fontsize = 13.0,
          col = "black"
        )
      ),
      grid::textGrob(
        label = title,
        x = grid::unit(0.50, "npc"),
        y = grid::unit(0.48, "npc"),
        just = c("center", "center"),
        gp = grid::gpar(
          fontfamily = base_family,
          fontface = "bold",
          fontsize = 9.2,
          col = COLOR_DARK
        )
      )
    )
  )
}

assemble_top_panel <- function(body, letter, title) {
  patchwork::wrap_plots(
    make_top_panel_header(letter, title),
    body,
    ncol = 1,
    heights = c(0.14, 1.00)
  )
}

# ------------------------------------------------------------------------------
# 5. Panel A — Tnfrsf1a detection and sequencing-depth imbalance
# ------------------------------------------------------------------------------
depth <- depth0 %>%
  mutate(
    sample_id = factor(sample_id, levels = SAMPLES),
    sample_label = factor(unname(sample_short[as.character(sample_id)]),
                          levels = unname(sample_short[SAMPLES])),
    detection = ifelse(as.logical(receptor_positive),
                       "Tnfrsf1a detected", "Not detected"),
    detection = factor(detection,
                       levels = c("Not detected", "Tnfrsf1a detected"))
  )
plot_a_detection_data <- depth %>%
  filter(detection == "Tnfrsf1a detected") %>%
  mutate(percent_detected = 100 * fraction_of_sample)

pA1 <- ggplot(plot_a_detection_data,
              aes(x = sample_label, y = percent_detected,
                  fill = as.character(sample_id))) +
  geom_col(width = 0.53, color = NA) +
  geom_text(aes(label = sprintf("%.1f%%", percent_detected)),
            vjust = -0.52, size = 2.75, family = base_family,
            color = COLOR_DARK) +
  scale_fill_manual(values = sample_colors) +
  scale_y_continuous(limits = c(0, 59), breaks = c(0, 20, 40),
                     labels = function(x) paste0(x, "%"),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = NULL, y = "Endothelial cells (%)") +
  pub_theme +
  theme(legend.position = "none", axis.ticks.x = element_blank(),
        plot.margin = margin(4.5, 3.2, 2.0, 4.2, unit = "mm"))

pA2 <- ggplot(depth,
              aes(x = sample_label, y = median_nUMI, fill = detection)) +
  geom_col(position = position_dodge(width = 0.77),
           width = 0.68, color = NA) +
  geom_text(aes(label = scales::comma(round(median_nUMI, 0))),
            position = position_dodge(width = 0.77),
            vjust = -0.35, size = 2.35, family = base_family,
            color = COLOR_DARK) +
  scale_fill_manual(values = c("Not detected" = COLOR_GREY,
                               "Tnfrsf1a detected" = COLOR_S1)) +
  scale_y_continuous(limits = c(0, 8750),
                     breaks = c(0, 2000, 4000, 6000, 8000),
                     labels = scales::comma,
                     expand = expansion(mult = c(0, 0))) +
  labs(title = "Median UMI counts by receptor detection",
       x = NULL, y = "Median nUMI", fill = NULL) +
  pub_theme +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        legend.key.width = grid::unit(2.6, "mm"),
        legend.margin = margin(t = 0.5, unit = "mm"),
        plot.title = element_text(size = 8.8, face = "bold", hjust = 0),
        plot.margin = margin(1.2, 3.2, 2.0, 4.2, unit = "mm")) +
  guides(fill = guide_legend(nrow = 1))

panelA_body <- patchwork::wrap_plots(
  pA1, pA2,
  ncol = 1,
  heights = c(0.86, 1.12)
)

panelA <- assemble_top_panel(
  body = panelA_body,
  letter = "A",
  title = "Tnfrsf1a detection"
)

# ------------------------------------------------------------------------------
# 6. Panel B — ORIGINAL gene sets only; no model is refitted
# ------------------------------------------------------------------------------
# Deliberately exclude NONOVERLAP receptor-positive vs receptor-undetected
# contrasts from Figure S9B. Panel C below retains the separate score-correlation
# sensitivity analysis before/after removing shared genes.
original_program_ids <- c("ENDOTHELIAL_ACTIVATION", "TNF_NFKB_RESPONSE")
original_program_labels <- c(
  "Endothelial activation",
  "TNF/NF-kappaB response"
)

plot_b_data <- mod0 %>%
  filter(
    gene_set == "ORIGINAL",
    program %in% original_program_ids,
    covariate_model %in% c("unadjusted", "adjust_both")
  ) %>%
  mutate(
    sample_id = factor(sample_id, levels = SAMPLES),
    sample_label = factor(unname(sample_short[as.character(sample_id)]),
                          levels = unname(sample_short[SAMPLES])),
    adjustment = ifelse(covariate_model == "unadjusted",
                        "Unadjusted", "Adjusted: nUMI + nFeature"),
    adjustment = factor(adjustment,
                        levels = c("Unadjusted", "Adjusted: nUMI + nFeature")),
    program_display = case_when(
      program == "ENDOTHELIAL_ACTIVATION" ~ original_program_labels[1],
      program == "TNF_NFKB_RESPONSE" ~ original_program_labels[2],
      TRUE ~ NA_character_
    )
  )

if (anyNA(plot_b_data$program_display) || nrow(plot_b_data) != 8) {
  stop("Panel B must contain only 2 original programs x 2 samples x 2 models (8 rows).")
}
plot_b_data$program_display <- factor(
  plot_b_data$program_display,
  levels = rev(original_program_labels)
)
plot_b_data <- plot_b_data %>%
  mutate(
    y_index = as.numeric(program_display),
    y_plot = y_index + ifelse(adjustment == "Unadjusted", 0.135, -0.135)
  )

pB <- ggplot(plot_b_data,
             aes(x = beta_receptor_positive, y = y_plot,
                 color = adjustment)) +
  geom_vline(xintercept = 0, linetype = "22", color = "#777777",
             linewidth = 0.42) +
  geom_point(size = 2.65, stroke = 0) +
  facet_grid(sample_label ~ ., switch = "y") +
  scale_color_manual(values = c("Unadjusted" = COLOR_S1,
                                "Adjusted: nUMI + nFeature" = COLOR_ADJUSTED)) +
  scale_y_continuous(
    breaks = seq_along(original_program_labels),
    labels = rev(original_program_labels),
    limits = c(0.55, 2.45),
    expand = expansion(mult = 0)
  ) +
  scale_x_continuous(breaks = seq(-0.05, 0.20, 0.05),
                     labels = scales::label_number(accuracy = 0.01)) +
  coord_cartesian(xlim = c(-0.075, 0.242), clip = "off") +
  labs(
    x = "Tnfrsf1a transcript detected vs undetected (score difference)",
    y = NULL, color = NULL
  ) +
  pub_theme +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 8.2, face = "bold"),
    panel.spacing.y = grid::unit(3.0, "mm"),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(size = 7.6, lineheight = 0.98),
    axis.text.x = element_text(size = 7.4),
    axis.title.x = element_text(size = 8.4),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.text = element_text(size = 7.2),
    legend.key.width = grid::unit(3.0, "mm"),
    plot.margin = margin(4.5, 3.5, 2.4, 4.4, unit = "mm")
  ) +
  guides(color = guide_legend(nrow = 1))

panelB <- assemble_top_panel(
  body = pB,
  letter = "B",
  title = "Endothelial score sensitivity"
)

# Export a revised Table S21 from the exact Step8G model output.
# Keep the original 14-gene and 20-gene programs only; exclude all NONOVERLAP
# positive-versus-undetected effect rows. All four model specifications remain.
table_s21_raw <- mod0 %>%
  filter(
    gene_set == "ORIGINAL",
    program %in% original_program_ids,
    covariate_model %in% expected_models
  ) %>%
  mutate(
    sample_order = match(as.character(sample_id), SAMPLES),
    program_order = match(program, original_program_ids),
    model_order = match(covariate_model, expected_models)
  ) %>%
  arrange(sample_order, program_order, model_order)

if (nrow(table_s21_raw) != 16 ||
    any(table_s21_raw$gene_set != "ORIGINAL") ||
    any(grepl("NONOVERLAP", table_s21_raw$program))) {
  stop("Revised Table S21 must contain 16 original-program model rows only.")
}

model_labels <- c(
  "unadjusted" = "Unadjusted",
  "adjust_log_nUMI" = "Adjusted for log(nUMI + 1)",
  "adjust_log_nFeature" = "Adjusted for log(nFeature + 1)",
  "adjust_both" = "Adjusted for log(nUMI + 1) and log(nFeature + 1)"
)
program_labels_table <- c(
  "ENDOTHELIAL_ACTIVATION" = "Endothelial activation (14 genes)",
  "TNF_NFKB_RESPONSE" = "TNF/NF-kB response (20 genes)"
)

table_s21_export <- table_s21_raw %>%
  mutate(
    sample = unname(sample_short[as.character(sample_id)]),
    program = unname(program_labels_table[as.character(program)]),
    covariate_model = unname(model_labels[as.character(covariate_model)])
  ) %>%
  select(
    sample, program, covariate_model,
    n_cells, n_pos, n_neg,
    beta_receptor_positive, se_cell_level_descriptive,
    ci95_low_cell_level_descriptive, ci95_high_cell_level_descriptive,
    positive_direction, model_r_squared
  )

write.csv(
  table_s21_export,
  file.path(OUT, "Table_S21_revised.csv"),
  row.names = FALSE,
  na = ""
)

wb_s21 <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb_s21, "Table S21")
openxlsx::writeData(
  wb_s21, "Table S21", table_s21_export,
  withFilter = TRUE, keepNA = FALSE
)
openxlsx::freezePane(wb_s21, "Table S21", firstRow = TRUE)
header_style_s21 <- openxlsx::createStyle(
  textDecoration = "bold", fgFill = "#D9E2F3",
  border = "bottom", wrapText = TRUE, valign = "center"
)
openxlsx::addStyle(
  wb_s21, "Table S21", header_style_s21,
  rows = 1, cols = seq_len(ncol(table_s21_export)),
  gridExpand = TRUE
)
openxlsx::setRowHeights(wb_s21, "Table S21", rows = 1, heights = 32)
openxlsx::setColWidths(
  wb_s21, "Table S21", cols = seq_len(ncol(table_s21_export)),
  widths = "auto"
)
numeric_style_s21 <- openxlsx::createStyle(numFmt = "0.0000")
openxlsx::addStyle(
  wb_s21, "Table S21", numeric_style_s21,
  rows = 2:(nrow(table_s21_export) + 1),
  cols = c(7:10, 12), gridExpand = TRUE, stack = TRUE
)

openxlsx::addWorksheet(wb_s21, "Notes")
notes_s21 <- data.frame(
  Item = c(
    "Scope",
    "Tnfrsf1a status",
    "Effect estimate",
    "Covariate adjustment",
    "Confidence intervals",
    "Biological replication",
    "Excluded results"
  ),
  Description = c(
    "Only the original endothelial activation (14-gene) and TNF/NF-kB response (20-gene) scores are included.",
    "Detected/undetected refers to Tnfrsf1a transcript detection, not receptor protein measurement.",
    "beta_receptor_positive is the within-sample score difference for transcript-detected versus transcript-undetected endothelial cells.",
    "Models are fit separately within each sample: unadjusted; adjusted for log(nUMI + 1); adjusted for log(nFeature + 1); and adjusted for both covariates.",
    "Intervals are descriptive cell-level 95% intervals and should not be interpreted as biological-replicate confidence intervals.",
    "Two 7-dpi samples are shown separately; cells are not independent biological replicates.",
    "All positive-versus-undetected effect rows calculated from the nonoverlapping/de-overlapped gene sets were removed from this revised table."
  ),
  stringsAsFactors = FALSE
)
openxlsx::writeData(wb_s21, "Notes", notes_s21, withFilter = FALSE)
openxlsx::addStyle(
  wb_s21, "Notes", header_style_s21,
  rows = 1, cols = seq_len(ncol(notes_s21)), gridExpand = TRUE
)
openxlsx::setColWidths(wb_s21, "Notes", cols = 1, widths = 24)
openxlsx::setColWidths(wb_s21, "Notes", cols = 2, widths = 105)
openxlsx::saveWorkbook(
  wb_s21, file.path(OUT, "Table_S21_revised.xlsx"), overwrite = TRUE
)

# ------------------------------------------------------------------------------
# 7. Panel C — shared-gene dependence of correlation (precomputed Spearman rho)
# ------------------------------------------------------------------------------
plot_c_data <- rbind(
  data.frame(sample_id = rho0$sample_id,
             version = "Original programs",
             rho = rho0$original_program_rho),
  data.frame(sample_id = rho0$sample_id,
             version = "8 shared genes removed",
             rho = rho0$nonoverlap_program_rho)
)
plot_c_data$sample_id <- factor(plot_c_data$sample_id, levels = SAMPLES)
plot_c_data$version <- factor(plot_c_data$version,
                              levels = c("Original programs",
                                         "8 shared genes removed"))

pC <- ggplot(plot_c_data,
             aes(x = version, y = rho,
                 group = sample_id, color = as.character(sample_id))) +
  geom_line(linewidth = 0.86) +
  geom_point(size = 3.6) +
  geom_text(aes(label = sprintf("%.3f", rho)),
            size = 2.85, color = COLOR_DARK,
            vjust = -1.05, family = base_family) +
  scale_color_manual(values = sample_colors,
                     labels = unname(sample_short[SAMPLES])) +
  scale_y_continuous(breaks = seq(0, 0.8, 0.2),
                     limits = c(0, 0.91),
                     expand = expansion(mult = c(0, 0))) +
  labs(tag = "C", title = "Program-score correlation",
       x = NULL, y = "Spearman rho", color = NULL) +
  pub_theme +
  theme(
    axis.text.x = element_text(size = 7.7),
    axis.text.y = element_text(size = 7.7),
    axis.ticks.x = element_blank(),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.text = element_text(size = 7.6),
    plot.margin = margin(4.0, 4.0, 2.2, 4.0, unit = "mm")
  ) +
  guides(color = guide_legend(nrow = 1))

# ------------------------------------------------------------------------------
# 8. Panel D — the ONE nonduplicated DotPlot, all refined Step8G populations
# ------------------------------------------------------------------------------
priority_groups <- c(
  "Target-program microglia (Arm_2)",
  "Other microglia",
  "Alternative-branch microglia (Arm_1)",
  "Infiltrating macrophage-like",
  "Infiltrating monocyte",
  "IFN-responsive infiltrating mono/mac",
  "BAM-like macrophage",
  "Neutrophil (Step8A exclusive)",
  "Other myeloid (unresolved)",
  "T_cell", "Endothelial", "Astrocyte", "Ependymal",
  "Fibroblast", "Pericyte_VSMC", "OPC",
  "Oligodendrocyte", "Neuron"
)
if (!setequal(unique(tnf0$cell_group), priority_groups)) {
  stop("Unexpected refined myeloid cell groups; update Figure S9 labels only ",
       "after checking the original Step8G data.")
}

labels_D <- c(
  "Target-program microglia (Arm_2)" = "Target microglia",
  "Other microglia" = "Other microglia",
  "Alternative-branch microglia (Arm_1)" = "Alternative microglia",
  "Infiltrating macrophage-like" = "Infiltrating macrophage-like",
  "Infiltrating monocyte" = "Infiltrating monocyte",
  "IFN-responsive infiltrating mono/mac" = "IFN-responsive mono/mac",
  "BAM-like macrophage" = "BAM-like macrophage",
  "Neutrophil (Step8A exclusive)" = "Neutrophils",
  "Other myeloid (unresolved)" = "Other myeloid (unresolved)",
  "T_cell" = "T cells", "Endothelial" = "Endothelial",
  "Astrocyte" = "Astrocyte", "Ependymal" = "Ependymal",
  "Fibroblast" = "Fibroblast", "Pericyte_VSMC" = "Pericyte/VSMC",
  "OPC" = "OPC", "Oligodendrocyte" = "Oligodendrocyte",
  "Neuron" = "Neuron"
)

plot_d_data <- tnf0 %>%
  mutate(
    sample_label = factor(unname(sample_short[sample_id]),
                          levels = unname(sample_short[SAMPLES])),
    cell_group = factor(cell_group, levels = rev(priority_groups))
  )

pD <- ggplot(plot_d_data,
             aes(x = sample_label, y = cell_group,
                 size = pct_Tnf_detected,
                 fill = mean_logRNA_all_cells)) +
  geom_point(shape = 21, color = "#676767", stroke = 0.24) +
  scale_y_discrete(labels = labels_D,
                   expand = expansion(add = 0.65)) +
  scale_x_discrete(expand = expansion(add = 0.5)) +
  scale_size_continuous(name = "Tnf detected (%)", range = c(0.8, 6.6),
                        limits = c(0, 100),
                        breaks = c(0, 25, 50, 75, 100)) +
  scale_fill_gradient(name = "Mean normalized Tnf",
                      low = "#E8E8E8", high = COLOR_S1,
                      limits = c(0, 1.85),
                      breaks = c(0, 0.6, 1.2, 1.8),
                      labels = scales::label_number(accuracy = 0.1)) +
  labs(tag = "D", title = "Cellular sources of Tnf at 7 dpi",
       x = NULL, y = NULL) +
  pub_theme +
  theme(
    axis.line = element_blank(), axis.ticks = element_blank(),
    axis.text.y = element_text(size = 7.4, color = COLOR_DARK),
    axis.text.x = element_text(size = 8.0, face = "bold"),
    panel.grid.major.y = element_line(color = "#F0F0F0", linewidth = 0.22),
    panel.grid.major.x = element_line(color = "#F2F2F2", linewidth = 0.22),
    legend.position = "bottom", legend.box = "vertical",
    legend.margin = margin(t = 0.2, unit = "mm"),
    legend.text = element_text(size = 7.2),
    legend.title = element_text(size = 7.5),
    plot.margin = margin(4.0, 3.8, 2.0, 4.1, unit = "mm")
  ) +
  guides(
    size = guide_legend(order = 1, nrow = 1,
                        override.aes = list(fill = "#B3BBC4")),
    fill = guide_colorbar(order = 2, barwidth = grid::unit(28, "mm"),
                          barheight = grid::unit(3.2, "mm"))
  )

# Key caveat documented in legend/run summary rather than drawing an ambiguous
# warning badge over the data: 7dpi_sample2 has only 13 neutrophils.
neut_row <- tnf0[tnf0$sample_id == "7dpi_sample2" &
                   tnf0$cell_group == "Neutrophil (Step8A exclusive)", ]
if (nrow(neut_row) != 1 || neut_row$n_cells != 13) {
  warning("Neutrophil cell count differs from the original Step8G audit; ",
          "check the updated figure legend.")
}

# ------------------------------------------------------------------------------
# 9. Frozen plotting data and complete assembled figure
# ------------------------------------------------------------------------------
write.csv(plot_a_detection_data,
          file.path(OUT, "02_FigureS9A_detection_plot_data.csv"), row.names = FALSE)
write.csv(depth,
          file.path(OUT, "02_FigureS9A_depth_plot_data.csv"), row.names = FALSE)
write.csv(plot_b_data,
          file.path(OUT, "02_FigureS9B_effect_plot_data.csv"), row.names = FALSE)
write.csv(plot_c_data,
          file.path(OUT, "02_FigureS9C_correlation_plot_data.csv"),
          row.names = FALSE)
write.csv(plot_d_data,
          file.path(OUT, "02_FigureS9D_Tnf_source_plot_data.csv"),
          row.names = FALSE)

# Protect against patchwork flattening nested subplots (same technique as S1).
panelA_fixed <- patchwork::wrap_elements(full = panelA, clip = FALSE)
panelB_fixed <- patchwork::wrap_elements(full = panelB, clip = FALSE)
panelC_fixed <- patchwork::wrap_elements(full = pC, clip = FALSE)
panelD_fixed <- patchwork::wrap_elements(full = pD, clip = FALSE)

row_top <- patchwork::wrap_plots(panelA_fixed, panelB_fixed, ncol = 2,
                                widths = c(0.405, 0.595))
row_bottom <- patchwork::wrap_plots(panelC_fixed, panelD_fixed, ncol = 2,
                                   widths = c(0.405, 0.595))
figureS9 <- patchwork::wrap_plots(
  patchwork::wrap_elements(full = row_top, clip = FALSE),
  patchwork::wrap_elements(full = row_bottom, clip = FALSE),
  ncol = 1, heights = c(0.91, 1.09)
) + patchwork::plot_annotation(
  theme = theme(plot.background = element_rect(fill = "white", color = NA))
)

# ------------------------------------------------------------------------------
# 10. High-resolution publication exports: each panel + full figure
# ------------------------------------------------------------------------------
save_plot_quad <- function(plot, stem, width_mm, height_mm, dpi = 600,
                           preview = FALSE) {
  filename <- file.path(OUT, stem)
  ggplot2::ggsave(
    filename = paste0(filename, ".pdf"), plot = plot,
    device = grDevices::cairo_pdf,
    width = width_mm, height = height_mm, units = "mm", bg = "white"
  )
  ggplot2::ggsave(
    filename = paste0(filename, ".svg"), plot = plot,
    device = svglite::svglite,
    width = width_mm, height = height_mm, units = "mm", bg = "white"
  )
  ggplot2::ggsave(
    filename = paste0(filename, ".tiff"), plot = plot,
    device = ragg::agg_tiff,
    width = width_mm, height = height_mm, units = "mm",
    dpi = dpi, compression = "lzw", background = "white"
  )
  if (preview) {
    ggplot2::ggsave(
      filename = paste0(filename, "_review.png"), plot = plot,
      device = ragg::agg_png,
      width = width_mm, height = height_mm, units = "mm",
      dpi = 220, background = "white"
    )
  }
}

message("Exporting Figure S9 panel A...")
save_plot_quad(panelA, "FigureS9A_receptor_detection_and_depth",
               width_mm = 83, height_mm = 104)
message("Exporting Figure S9 panel B...")
save_plot_quad(panelB, "FigureS9B_endothelial_effect_sensitivity",
               width_mm = 118, height_mm = 104)
message("Exporting Figure S9 panel C...")
save_plot_quad(pC, "FigureS9C_gene_overlap_correlation",
               width_mm = 83, height_mm = 122)
message("Exporting Figure S9 panel D...")
save_plot_quad(pD, "FigureS9D_Tnf_cellular_sources",
               width_mm = 118, height_mm = 122)
message("Exporting complete Figure S9...")
save_plot_quad(figureS9, "FigureS9_GLIA_FINAL_v4_aligned_large_font",
               width_mm = 200, height_mm = 232, preview = TRUE)

# ------------------------------------------------------------------------------
# 11. Reproducibility notes and caveats
# ------------------------------------------------------------------------------
summary_lines <- c(
  "Figure S9 | GLIA | Major Comment 5 | FINAL v2",
  "================================================",
  paste0("Original project root: ", ROOT),
  paste0("Original Step8G result: ", STEP8G),
  paste0("Final figure output:  ", OUT), "",
  "Panel A: Tnfrsf1a detected fractions and nUMI medians by receptor detection.",
  "Panel B: Original 14-gene endothelial activation and 20-gene TNF/NF-kB",
  "         receptor-detected vs undetected contrasts only; no NONOVERLAP contrasts.",
  "         Shows unadjusted and dual-depth-adjusted point estimates; descriptive cell-level",
  "         intervals remain available in Table S21 but are not drawn in the panel.",
  "Table S21 revised: original programs only; all four covariate models retained.",
  "         NONOVERLAP positive-vs-undetected effect rows were excluded.",
  "Panel C: Original vs. overlap-removed cross-program Spearman correlations.",
  "         The eight shared genes are removed from BOTH programs.",
  "Panel D: All 14505 original Step8G 7-dpi cells, grouped into 18",
  "         nonoverlapping categories for each of two biological samples.",
  "         Dot size = fraction with detected Tnf mRNA; fill = average",
  "         log-normalized Tnf RNA across all cells in that group.",
  "         Sample 2 contains only 13 neutrophils; interpret cautiously.", "",
  "All statistical values were imported directly from Step8G CSV outputs.",
  "No inference of ligand secretion, TNFRSF1A signaling flux, or causality is made.",
  "DO NOT describe receptor-undetected cells as receptor-negative by protein assay.",
  "Figure S9 is a four-panel supplemental figure; individual panels are also exported.",
  "Panel subtitles were removed; Panel B horizontal CI bars were removed; Panel D display labels omit Arm_1/Arm_2.",
  "Panels A and B use one shared external-header geometry so their panel letters are exactly level; all reader-facing fonts were enlarged.",
  paste0("Step8G input validation passed: ", sum(qc$PASS), "/", nrow(qc)),
  ""
)
writeLines(summary_lines,
           file.path(OUT, "FigureS9_GLIA_FINAL_v4_aligned_large_font_run_summary.txt"),
           useBytes = TRUE)
capture.output(sessionInfo(), file = file.path(OUT, "sessionInfo.txt"))
warning_obj <- warnings()
if (!is.null(warning_obj)) {
  capture.output(warning_obj,
                 file = file.path(OUT, "03_warnings_after_run.txt"))
}

message("===============================================================")
message("Figure S9 FINAL v2 complete. All 11 data audits passed.")
message("Publication vector file: FigureS9_GLIA_FINAL_v4_aligned_large_font.pdf / .svg")
message("Publication 600-dpi TIFF: FigureS9_GLIA_FINAL_v4_aligned_large_font.tiff")
message("Preview file: FigureS9_GLIA_FINAL_v4_aligned_large_font_review.png")
message("Output: ", OUT)
message("===============================================================")
