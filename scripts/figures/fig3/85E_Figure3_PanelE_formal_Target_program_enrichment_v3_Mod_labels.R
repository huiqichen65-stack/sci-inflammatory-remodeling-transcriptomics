# ==============================================================================
# 85E_Figure3_PanelE_formal_Target_program_enrichment_v3_Mod_labels.R
#
# Figure 3E | Candidate-regulon target enrichment in the formal Target program
#
# FINAL DESIGN
# ------------
# v3: display-only terminology update: M1/M2 -> Mod1/Mod2;
#     frozen data, formal set IDs, statistics, and visual encodings are unchanged.
# Figure 3C/D identified four earliest branch-associated regulons:
#   Nfkb1(+), Nfe2l2(+), Relb(+), Nfkb2(+)
#
# Figure 3E asks one focused question:
#   Are the motif-supported targets of these four regulons enriched in the
#   formal Target-program gene set defined by Step6B Mod1 + Mod2?
#
# IMPORTANT:
# - This panel does NOT rerun target inference, module discovery, enrichment,
#   pseudotime, or any statistical test.
# - It reads the already-frozen Step7G target-program closure table.
# - The displayed enrichment statistics are specifically those for
#   Step6B_formal_Arm2_M1_M2.
# - Mod1 and Mod2 are shown as one formal Target-program gene set.
# - Progressive / strong-progressive / high-confidence marker domains are NOT
#   displayed in the main panel.
# - No causal or "driver" claim is made.
#
# VISUAL ENCODING
# ---------------
# x position   = frozen formal-module fold enrichment
# point size   = number of overlapping target genes
# point fill   = -log10(BH q), capped at 20 for display
# dashed x = 1 = no enrichment
#
# OUTPUT
# ------
# <PROJECT_ROOT>/09_final_figures/GLIA/Figure3E_v2/
#
#   Figure3E_formal_Target_program_enrichment_v2.pdf
#   Figure3E_formal_Target_program_enrichment_v2.svg
#   Figure3E_formal_Target_program_enrichment_v2.tiff
#   Figure3E_formal_Target_program_enrichment_v2.png
#   00_input_audit.csv
#   01_schema_audit.csv
#   02_Figure3E_plot_data.csv
#   03_Figure3E_key_result_audit.csv
#   source_file_md5.csv
#   sessionInfo.txt
# ==============================================================================

rm(list = ls())
gc()
options(stringsAsFactors = FALSE, scipen = 999)

# ------------------------------------------------------------------------------
# 0. USER SETTINGS
# ------------------------------------------------------------------------------
PROJECT_ROOT <- Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
OUT_DIR <- file.path(
  PROJECT_ROOT,
  "09_final_figures", "GLIA", "Figure3E_v3_Mod_labels"
)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

FONT_FAMILY <- "Arial"

# The same four regulons already prioritized visually in Figure 3C/D.
KEY_TFS <- c("Nfkb1", "Nfe2l2", "Relb", "Nfkb2")

# Frozen formal Target-program set used by Step7G.
EXPECTED_FORMAL_SET <- "Step6B_formal_Arm2_M1_M2"

# Figure palette consistent with accepted Figure 3 panels.
COL_TARGET      <- "#B64B45"
COL_TARGET_MID  <- "#D98C84"
COL_LIGHT       <- "#F3F3F3"
COL_GREY        <- "#B8B8B8"
COL_DARK        <- "#222222"
COL_NULL        <- "#8A8A8A"

# Display cap only; underlying q values remain untouched in the exported audit.
NEGLOG10Q_CAP <- 20

FIG_WIDTH_MM  <- 92
FIG_HEIGHT_MM <- 68

# ------------------------------------------------------------------------------
# 1. PACKAGES
# ------------------------------------------------------------------------------
required_pkgs <- c(
  "data.table",
  "ggplot2",
  "svglite",
  "ragg",
  "scales"
)

missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_pkgs) > 0L) {
  stop(
    "Missing required R package(s): ",
    paste(missing_pkgs, collapse = ", "),
    "\nInstall them first, then rerun.",
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# ------------------------------------------------------------------------------
# 2. FROZEN INPUT
# ------------------------------------------------------------------------------
STEP7G_CLOSURE <- file.path(
  PROJECT_ROOT,
  "04_results",
  "step7G_GSE162610_regulon_target_program_closure",
  "tables",
  "03_timing_candidate_target_closure.tsv"
)

input_audit <- data.frame(
  file = "Step7G frozen candidate target-closure table",
  path = STEP7G_CLOSURE,
  exists = file.exists(STEP7G_CLOSURE),
  size_bytes = if (file.exists(STEP7G_CLOSURE)) {
    file.info(STEP7G_CLOSURE)$size
  } else {
    NA_real_
  },
  stringsAsFactors = FALSE
)

write.csv(
  input_audit,
  file.path(OUT_DIR, "00_input_audit.csv"),
  row.names = FALSE
)

if (!input_audit$exists) {
  stop(
    "Frozen Step7G closure table not found. Inspect 00_input_audit.csv.",
    call. = FALSE
  )
}

closure <- data.table::fread(STEP7G_CLOSURE, data.table = FALSE)

# ------------------------------------------------------------------------------
# 3. SCHEMA AUDIT
# ------------------------------------------------------------------------------
required_cols <- c(
  "regulon",
  "TF",
  "n_targets",
  "onset_progress_qualified",
  "seven_dpi_direction_consistent",
  "formal_module_set",
  "formal_module_overlap_n",
  "formal_module_fold_enrichment",
  "formal_module_q",
  "formal_module_support",
  "target_closure_strength",
  "mechanistic_priority"
)

schema_audit <- data.frame(
  source = "03_timing_candidate_target_closure.tsv",
  item = required_cols,
  present = required_cols %in% colnames(closure),
  stringsAsFactors = FALSE
)

write.csv(
  schema_audit,
  file.path(OUT_DIR, "01_schema_audit.csv"),
  row.names = FALSE
)

if (!all(schema_audit$present)) {
  stop(
    "Schema audit failed. Inspect 01_schema_audit.csv.",
    call. = FALSE
  )
}

# ------------------------------------------------------------------------------
# 4. SELECT THE SAME FOUR REGULONS USED IN FIGURE 3C/D
# ------------------------------------------------------------------------------
d <- closure[
  closure$TF %in% KEY_TFS,
  required_cols,
  drop = FALSE
]

# Require exactly one frozen row per requested TF.
tf_counts <- table(d$TF)

if (
  length(tf_counts) != length(KEY_TFS) ||
  any(!KEY_TFS %in% names(tf_counts)) ||
  any(tf_counts[KEY_TFS] != 1L)
) {
  stop(
    "Expected exactly one frozen Step7G row for each of: ",
    paste(KEY_TFS, collapse = ", "),
    call. = FALSE
  )
}

# Put rows back into narrative order.
d <- d[match(KEY_TFS, d$TF), , drop = FALSE]

# Confirm that Figure 3E is using the intended M1+M2 formal gene set.
if (any(is.na(d$formal_module_set)) ||
    any(d$formal_module_set != EXPECTED_FORMAL_SET)) {
  stop(
    "Unexpected formal_module_set detected. Expected: ",
    EXPECTED_FORMAL_SET,
    call. = FALSE
  )
}

# Numerical safety checks.
numeric_cols <- c(
  "n_targets",
  "formal_module_overlap_n",
  "formal_module_fold_enrichment",
  "formal_module_q"
)

for (nm in numeric_cols) {
  d[[nm]] <- as.numeric(d[[nm]])
}

if (
  anyNA(d$formal_module_overlap_n) ||
  anyNA(d$formal_module_fold_enrichment) ||
  anyNA(d$formal_module_q)
) {
  stop(
    "Missing frozen formal-module enrichment statistic(s) for Figure 3E.",
    call. = FALSE
  )
}

if (any(d$formal_module_q < 0 | d$formal_module_q > 1)) {
  stop("formal_module_q contains values outside [0, 1].", call. = FALSE)
}

# ------------------------------------------------------------------------------
# 5. DISPLAY VARIABLES
# ------------------------------------------------------------------------------
# Use the actual frozen regulon name and add a visual gap before "(+)".
d$display_label <- sub(
  "\\(\\+\\)$",
  " (+)",
  as.character(d$regulon)
)

# Keep the C/D narrative order from top to bottom.
d$display_label <- factor(
  d$display_label,
  levels = rev(sub("\\(\\+\\)$", " (+)", paste0(KEY_TFS, "(+)")))
)

# Display significance; q itself remains unchanged in output CSV.
d$neglog10_q_raw <- -log10(pmax(d$formal_module_q, 1e-300))
d$neglog10_q_display <- pmin(d$neglog10_q_raw, NEGLOG10Q_CAP)

# Useful audit-only categorical label.
d$BH_q_lt_0_05 <- d$formal_module_q < 0.05

# ------------------------------------------------------------------------------
# 6. EXPORT FROZEN PLOT DATA / AUDIT
# ------------------------------------------------------------------------------
plot_export <- d[, c(
  "regulon",
  "TF",
  "display_label",
  "n_targets",
  "onset_progress_qualified",
  "seven_dpi_direction_consistent",
  "formal_module_set",
  "formal_module_overlap_n",
  "formal_module_fold_enrichment",
  "formal_module_q",
  "formal_module_support",
  "neglog10_q_raw",
  "neglog10_q_display",
  "BH_q_lt_0_05",
  "target_closure_strength",
  "mechanistic_priority"
)]

write.csv(
  plot_export,
  file.path(OUT_DIR, "02_Figure3E_plot_data.csv"),
  row.names = FALSE
)

key_result_audit <- d[, c(
  "regulon",
  "TF",
  "n_targets",
  "formal_module_overlap_n",
  "formal_module_fold_enrichment",
  "formal_module_q",
  "formal_module_support",
  "target_closure_strength",
  "mechanistic_priority"
)]

write.csv(
  key_result_audit,
  file.path(OUT_DIR, "03_Figure3E_key_result_audit.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 7. PLOT
# ------------------------------------------------------------------------------
# Dynamic x-axis upper bound so labels do not get cramped.
x_max_data <- max(d$formal_module_fold_enrichment, na.rm = TRUE)
x_upper <- ceiling((x_max_data + 0.35) * 2) / 2
x_upper <- max(x_upper, 6.5)

# Stems extend from the null expectation (fold enrichment = 1).
p <- ggplot(
  d,
  aes(
    x = formal_module_fold_enrichment,
    y = display_label
  )
) +
  # Null expectation.
  geom_vline(
    xintercept = 1,
    linetype = 2,
    linewidth = 0.45,
    colour = COL_NULL
  ) +
  # Lollipop stems.
  geom_segment(
    aes(
      x = 1,
      xend = formal_module_fold_enrichment,
      yend = display_label
    ),
    linewidth = 0.75,
    colour = COL_GREY,
    lineend = "round"
  ) +
  # Main points.
  geom_point(
    aes(
      size = formal_module_overlap_n,
      fill = neglog10_q_display
    ),
    shape = 21,
    colour = "white",
    stroke = 0.55
  ) +
  scale_size_continuous(
    range = c(3.1, 7.8),
    breaks = scales::pretty_breaks(n = 3),
    name = "Overlap genes"
  ) +
  scale_fill_gradientn(
    colours = c(
      "#F2F2F2",
      "#E8C8C3",
      COL_TARGET_MID,
      COL_TARGET
    ),
    limits = c(0, NEGLOG10Q_CAP),
    breaks = c(0, 5, 10, 15, 20),
    labels = c("0", "5", "10", "15", "\u226520"),
    oob = scales::squish,
    name = expression(-log[10]("BH q"))
  ) +
  scale_x_continuous(
    limits = c(0.75, x_upper),
    breaks = seq(1, floor(x_upper), by = 1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    x = "Fold enrichment in target-terminal modules (Mod1/Mod2)",
    y = NULL,
    tag = "E"
  ) +
  theme_classic(
    base_family = FONT_FAMILY,
    base_size = 8
  ) +
  theme(
    plot.tag = element_text(
      family = FONT_FAMILY,
      size = 12,
      face = "bold",
      colour = COL_DARK
    ),
    plot.tag.position = c(0.006, 0.995),

    axis.title.x = element_text(
      family = FONT_FAMILY,
      size = 8,
      colour = COL_DARK,
      margin = margin(t = 5)
    ),
    axis.text.x = element_text(
      family = FONT_FAMILY,
      size = 7,
      colour = COL_DARK
    ),
    axis.text.y = element_text(
      family = FONT_FAMILY,
      size = 8,
      colour = COL_DARK,
      face = "plain",
      margin = margin(r = 4)
    ),
    axis.ticks = element_line(
      linewidth = 0.35,
      colour = COL_DARK
    ),
    axis.line = element_line(
      linewidth = 0.45,
      colour = COL_DARK
    ),

    legend.position = "right",
    legend.box = "vertical",
    legend.box.spacing = grid::unit(1.2, "mm"),
    legend.margin = margin(0, 0, 0, 0),
    legend.spacing.y = grid::unit(1.0, "mm"),
    legend.title = element_text(
      family = FONT_FAMILY,
      size = 7,
      colour = COL_DARK
    ),
    legend.text = element_text(
      family = FONT_FAMILY,
      size = 6.5,
      colour = COL_DARK
    ),
    legend.key.height = grid::unit(3.0, "mm"),

    plot.margin = margin(7, 2, 5, 7)
  ) +
  guides(
    size = guide_legend(
      order = 1,
      override.aes = list(
        fill = COL_TARGET_MID,
        colour = "white"
      )
    ),
    fill = guide_colourbar(
      order = 2,
      barheight = grid::unit(17, "mm"),
      barwidth = grid::unit(3.0, "mm"),
      ticks.colour = COL_DARK,
      frame.colour = COL_DARK
    )
  )

Figure3E <- p

# ------------------------------------------------------------------------------
# 8. EXPORT FIGURE
# ------------------------------------------------------------------------------
WIDTH_IN  <- FIG_WIDTH_MM / 25.4
HEIGHT_IN <- FIG_HEIGHT_MM / 25.4

PDF_FILE <- file.path(
  OUT_DIR,
  "Figure3E_formal_Target_program_enrichment_v3_Mod_labels.pdf"
)

SVG_FILE <- file.path(
  OUT_DIR,
  "Figure3E_formal_Target_program_enrichment_v3_Mod_labels.svg"
)

TIFF_FILE <- file.path(
  OUT_DIR,
  "Figure3E_formal_Target_program_enrichment_v3_Mod_labels.tiff"
)

PNG_FILE <- file.path(
  OUT_DIR,
  "Figure3E_formal_Target_program_enrichment_v3_Mod_labels.png"
)

ggsave(
  filename = PDF_FILE,
  plot = Figure3E,
  device = function(filename, width, height, ...) {
    grDevices::cairo_pdf(
      filename = filename,
      width = width,
      height = height,
      family = FONT_FAMILY,
      onefile = TRUE
    )
  },
  width = WIDTH_IN,
  height = HEIGHT_IN,
  units = "in",
  bg = "white"
)

ggsave(
  filename = SVG_FILE,
  plot = Figure3E,
  device = svglite::svglite,
  width = WIDTH_IN,
  height = HEIGHT_IN,
  units = "in",
  bg = "white"
)

ggsave(
  filename = TIFF_FILE,
  plot = Figure3E,
  device = ragg::agg_tiff,
  width = WIDTH_IN,
  height = HEIGHT_IN,
  units = "in",
  res = 600,
  compression = "lzw",
  background = "white"
)

ggsave(
  filename = PNG_FILE,
  plot = Figure3E,
  device = ragg::agg_png,
  width = WIDTH_IN,
  height = HEIGHT_IN,
  units = "in",
  res = 300,
  background = "white"
)

# ------------------------------------------------------------------------------
# 9. SOURCE MD5 / SESSION INFO
# ------------------------------------------------------------------------------
source_md5 <- data.frame(
  file = STEP7G_CLOSURE,
  md5 = unname(tools::md5sum(STEP7G_CLOSURE)),
  stringsAsFactors = FALSE
)

write.csv(
  source_md5,
  file.path(OUT_DIR, "source_file_md5.csv"),
  row.names = FALSE
)

writeLines(
  capture.output(sessionInfo()),
  file.path(OUT_DIR, "sessionInfo.txt")
)

# ------------------------------------------------------------------------------
# 10. CONSOLE SUMMARY
# ------------------------------------------------------------------------------
cat("\n============================================================\n")
cat("Figure 3E v3 completed successfully.\n")
cat("Output directory:\n", OUT_DIR, "\n", sep = "")

cat("\nFigure logic:\n")
cat("- same four branch-associated regulons as Figure 3C/D\n")
cat("- formal Target-program gene set = frozen Step6B Mod1 + Mod2 (display terminology only)\n")
cat("- x = fold enrichment\n")
cat("- point size = overlapping targets\n")
cat("- point fill = -log10(BH q), display-capped at 20\n")
cat("- dashed x = 1 marks no enrichment\n")

cat("\nFrozen formal Target-program enrichment:\n")
print(
  d[, c(
    "regulon",
    "n_targets",
    "formal_module_overlap_n",
    "formal_module_fold_enrichment",
    "formal_module_q",
    "formal_module_support"
  )],
  row.names = FALSE
)

cat("\nInterpretation boundary:\n")
cat("- Figure 3E shows target-program overlap/enrichment only.\n")
cat("- It does not establish direct causal control of the Target program.\n")
cat("- Regulon targets are frozen motif-supported predicted targets from upstream analysis.\n")
cat("============================================================\n")
