# ==============================================================================
# 85F_Figure3_PanelF_v9_black_borders.R
#
# Figure 3F | Top-ranked Target-program targets of Nfkb1 and Nfe2l2
#
# FINAL DESIGN
# ------------
# - Step7G only; NO CellChat input.
# - Keep the same objective ranking rule for both regulons:
#       evidence_count (descending) -> pySCENIC weight (descending)
# - Display ONLY the Top 10 frozen priority targets per regulon.
# - x = within-regulon pySCENIC target-weight percentile
# - point size = evidence_count
# - point FILL = formal Target-program module (M1 / M2)
# - point BORDER = shared priority target vs TF-prioritized only
# - gene labels are placed in a compact far-left in-panel gutter.
# - horizontal stems begin immediately to the right of the longest labels,
#   leaving only a small visual gap and never crossing the text.
# - bubble size is reduced so adjacent ranks do not overlap.
# - extra right-side plotting room prevents 100% points from being clipped.
# - NO gene-specific highlighting.
# - NO special treatment of Tnf / Ccl3 / Icam1.
#
# INTERPRETATION
# --------------
# This panel visualizes the highest-ranked motif-supported predicted targets
# that connect Nfkb1(+) and Nfe2l2(+) to the formal Target program (M1/M2).
# "Evidence count" is the number of supported Step7G target-closure evidence
# layers (maximum 5); it is NOT regulatory strength.
#
# OUTPUT
# ------
# <PROJECT_ROOT>/09_final_figures/GLIA/Figure3F_v7/
# ==============================================================================

rm(list = ls())
gc()
options(stringsAsFactors = FALSE, scipen = 999)

# ------------------------------------------------------------------------------
# 0. SETTINGS
# ------------------------------------------------------------------------------
PROJECT_ROOT <- Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
OUT_DIR <- file.path(
  PROJECT_ROOT,
  "09_final_figures", "GLIA", "Figure3F_v9_black_borders"
)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

FONT_FAMILY <- "Arial"
KEY_TFS <- c("Nfkb1", "Nfe2l2")
TOP_N <- 10L
EXPECTED_PRIORITY_EDGES_PER_TF <- 30L

# Stronger fills than v2 so M1/M2 remain clearly visible in PDF/SVG/PPT.
COL_M1 <- "#8F7CC3"
COL_M2 <- "#C96F63"

# Border encodes shared status only.
COL_SHARED_BORDER <- "#222222"
COL_UNIQUE_BORDER <- "#B9B9B9"

COL_TEXT <- "#222222"
COL_MUTED <- "#666666"
COL_STEM <- "#D8D8D8"
COL_GRID <- "#ECECEC"

FIG_WIDTH_MM <- 118
FIG_HEIGHT_MM <- 80

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
    "Missing required package(s): ",
    paste(missing_pkgs, collapse = ", "),
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
STEP7G_EDGES <- file.path(
  PROJECT_ROOT,
  "04_results",
  "step7G_GSE162610_regulon_target_program_closure",
  "tables",
  "05_priority_TF_target_network_edges.tsv"
)

input_audit <- data.frame(
  file = "05_priority_TF_target_network_edges.tsv",
  path = STEP7G_EDGES,
  exists = file.exists(STEP7G_EDGES),
  size_bytes = if (file.exists(STEP7G_EDGES)) {
    file.info(STEP7G_EDGES)$size
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

if (!file.exists(STEP7G_EDGES)) {
  stop(
    "Frozen Step7G priority-edge table not found. ",
    "Inspect 00_input_audit.csv.",
    call. = FALSE
  )
}

d0 <- data.table::fread(
  STEP7G_EDGES,
  data.table = FALSE,
  showProgress = FALSE
)

# ------------------------------------------------------------------------------
# 3. SCHEMA AUDIT
# ------------------------------------------------------------------------------
required_cols <- c(
  "TF",
  "regulon",
  "target",
  "weight",
  "weight_percentile_within_regulon",
  "evidence_count",
  "Step6B_module",
  "formal_Arm2_module",
  "mechanistic_priority"
)

schema_audit <- data.frame(
  source = "05_priority_TF_target_network_edges.tsv",
  column = required_cols,
  present = required_cols %in% colnames(d0),
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
# 4. SELECT Nfkb1 / Nfe2l2
# ------------------------------------------------------------------------------
d <- d0[
  d0$TF %in% KEY_TFS,
  required_cols,
  drop = FALSE
]

d$TF <- as.character(d$TF)
d$regulon <- as.character(d$regulon)
d$target <- as.character(d$target)
d$Step6B_module <- as.character(d$Step6B_module)
d$weight <- as.numeric(d$weight)
d$weight_percentile_within_regulon <- as.numeric(
  d$weight_percentile_within_regulon
)
d$evidence_count <- as.integer(d$evidence_count)

if (anyNA(d$weight) || any(d$weight <= 0)) {
  stop("Invalid frozen pySCENIC target weight detected.", call. = FALSE)
}

if (anyNA(d$weight_percentile_within_regulon)) {
  stop("Missing within-regulon target-weight percentile.", call. = FALSE)
}

if (anyNA(d$evidence_count)) {
  stop("Missing evidence_count.", call. = FALSE)
}

if (anyDuplicated(d[, c("TF", "target")])) {
  stop("Duplicate TF-target pair detected.", call. = FALSE)
}

# Frozen network-size check.
edge_counts <- table(d$TF)

if (
  any(!KEY_TFS %in% names(edge_counts)) ||
  any(
    as.integer(edge_counts[KEY_TFS]) !=
      rep(EXPECTED_PRIORITY_EDGES_PER_TF, length(KEY_TFS))
  )
) {
  stop(
    "Expected exactly 30 frozen Step7G priority targets per regulon.",
    call. = FALSE
  )
}

# ------------------------------------------------------------------------------
# 5. FORMAL TARGET-PROGRAM MODULE + SHARED STATUS
# ------------------------------------------------------------------------------
formal_flag <- as.logical(d$formal_Arm2_module)

if (anyNA(formal_flag) || any(!formal_flag)) {
  stop(
    "At least one selected frozen priority target is not in formal M1/M2.",
    call. = FALSE
  )
}

d$module_display <- ifelse(
  grepl("^M1_", d$Step6B_module),
  "Mod1",
  ifelse(
    grepl("^M2_", d$Step6B_module),
    "Mod2",
    NA_character_
  )
)

if (anyNA(d$module_display)) {
  stop(
    "Unexpected Step6B module among Nfkb1/Nfe2l2 priority targets.",
    call. = FALSE
  )
}

nfkb1_targets <- unique(d$target[d$TF == "Nfkb1"])
nfe2l2_targets <- unique(d$target[d$TF == "Nfe2l2"])
shared_targets <- intersect(nfkb1_targets, nfe2l2_targets)

d$shared_priority_target <- d$target %in% shared_targets

d$target_status <- ifelse(
  d$shared_priority_target,
  "Shared priority target",
  "TF-prioritized only"
)

# ------------------------------------------------------------------------------
# 6. OBJECTIVE RANKING + TOP 10 ONLY
# ------------------------------------------------------------------------------
# Same ranking rule for both regulons:
#   1) evidence_count descending
#   2) frozen pySCENIC weight descending
#   3) gene symbol only as deterministic tie-break
rank_parts <- lapply(KEY_TFS, function(tf_name) {

  x <- d[d$TF == tf_name, , drop = FALSE]

  x <- x[
    order(
      -x$evidence_count,
      -x$weight,
      x$target
    ),
    ,
    drop = FALSE
  ]

  x$full_rank <- seq_len(nrow(x))

  x <- x[
    seq_len(min(TOP_N, nrow(x))),
    ,
    drop = FALSE
  ]

  x$display_rank <- seq_len(nrow(x))

  x
})

plot_d <- do.call(rbind, rank_parts)
rownames(plot_d) <- NULL

if (
  any(
    table(plot_d$TF)[KEY_TFS] != TOP_N
  )
) {
  stop(
    "Top-10 selection failed for one or both regulons.",
    call. = FALSE
  )
}

# Convert percentile robustly to 0-100 for display.
pct_max <- max(
  plot_d$weight_percentile_within_regulon,
  na.rm = TRUE
)

plot_d$weight_percentile_display <- if (pct_max <= 1.000001) {
  plot_d$weight_percentile_within_regulon * 100
} else {
  plot_d$weight_percentile_within_regulon
}

if (
  any(
    plot_d$weight_percentile_display < -1e-6 |
      plot_d$weight_percentile_display > 100.0001
  )
) {
  stop(
    "Unexpected within-regulon target-weight percentile range.",
    call. = FALSE
  )
}

# Build a dedicated left-side gene-name gutter INSIDE each facet.
# This avoids ggplot dropping text outside scale limits and guarantees that
# horizontal stems never pass through gene names.
LABEL_X <- 1.5
STEM_START_X <- 27
PERCENTILE_MAX_X <- 100
PLOT_X_MAX <- 108

plot_d$label_x <- LABEL_X
plot_d$x_point <- STEM_START_X +
  (plot_d$weight_percentile_display / 100) *
  (PERCENTILE_MAX_X - STEM_START_X)

plot_d$TF_display <- factor(
  plot_d$TF,
  levels = c("Nfkb1", "Nfe2l2"),
  labels = c("Nfkb1 (+)", "Nfe2l2 (+)")
)

# ------------------------------------------------------------------------------
# 7. AUDIT EXPORTS
# ------------------------------------------------------------------------------
write.csv(
  plot_d,
  file.path(OUT_DIR, "02_Figure3F_top10_plot_data.csv"),
  row.names = FALSE
)

full_rank_audit <- do.call(rbind, lapply(KEY_TFS, function(tf_name) {

  x <- d[d$TF == tf_name, , drop = FALSE]

  x <- x[
    order(
      -x$evidence_count,
      -x$weight,
      x$target
    ),
    ,
    drop = FALSE
  ]

  x$full_rank <- seq_len(nrow(x))

  x
}))

write.csv(
  full_rank_audit,
  file.path(OUT_DIR, "03_Figure3F_full30_rank_audit.csv"),
  row.names = FALSE
)

summary_table <- data.frame(
  metric = c(
    "Nfkb1_frozen_priority_targets",
    "Nfe2l2_frozen_priority_targets",
    "shared_priority_targets_in_full30",
    "displayed_targets_per_regulon",
    "displayed_total_points"
  ),
  value = c(
    length(nfkb1_targets),
    length(nfe2l2_targets),
    length(shared_targets),
    TOP_N,
    nrow(plot_d)
  ),
  stringsAsFactors = FALSE
)

write.csv(
  summary_table,
  file.path(OUT_DIR, "04_Figure3F_summary.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 8. PLOT
# ------------------------------------------------------------------------------
# One single shape-21 point layer maps:
#   fill   -> Mod1/Mod2
#   border -> uniform black
# This keeps the panel visually simple while preserving module identity.

p <- ggplot(
  plot_d,
  aes(
    x = x_point,
    y = display_rank
  )
) +

  # Quiet lollipop stems.
  geom_segment(
    aes(
      x = STEM_START_X,
      xend = x_point,
      yend = display_rank
    ),
    linewidth = 0.38,
    colour = COL_STEM,
    lineend = "round"
  ) +

  # Single point layer:
  # fill = Mod1/Mod2; border = uniform black.
  geom_point(
    aes(
      size = evidence_count,
      fill = module_display
    ),
    shape = 21,
    colour = "black",
    stroke = 0.95
  ) +

  # Gene labels use the same rule for both regulons: all displayed Top 10.
# Labels are placed in a compact far-left gutter; stems start at STEM_START_X.
  geom_text(
    aes(
      x = label_x,
      label = target
    ),
    hjust = 0,
    family = FONT_FAMILY,
    size = 2.65,
    colour = COL_TEXT
  ) +

  facet_wrap(
    ~ TF_display,
    nrow = 1
  ) +

  scale_y_reverse(
    breaks = 1:10,
    expand = expansion(mult = c(0.01, 0.02))
  ) +

  scale_x_continuous(
    limits = c(0, PLOT_X_MAX),
    breaks = STEM_START_X +
      (c(0, 25, 50, 75, 100) / 100) *
      (PERCENTILE_MAX_X - STEM_START_X),
    labels = c("0%", "25%", "50%", "75%", "100%"),
    expand = expansion(mult = c(0, 0))
  ) +

  scale_size_continuous(
    range = c(2.7, 4.5),
    breaks = sort(unique(plot_d$evidence_count)),
    name = "Evidence count",
    guide = "none"
  ) +

  scale_fill_manual(
    values = c(
      Mod1 = COL_M1,
      Mod2 = COL_M2
    ),
    breaks = c("Mod1", "Mod2"),
    name = "Target-program module"
  ) +

  labs(
    tag = "F",
    x = "Within-regulon target-weight percentile",
    y = "Priority rank"
  ) +

  coord_cartesian(
    xlim = c(0, PLOT_X_MAX),
    ylim = c(10.75, 0.25),
    clip = "on"
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
      colour = COL_TEXT
    ),
    plot.tag.position = c(0.005, 0.995),

    strip.background = element_blank(),
    strip.text = element_text(
      family = FONT_FAMILY,
      size = 9.4,
      face = "bold",
      colour = COL_TEXT,
      margin = margin(b = 4)
    ),

    axis.title = element_text(
      family = FONT_FAMILY,
      size = 7.7,
      colour = COL_TEXT
    ),
    axis.text = element_text(
      family = FONT_FAMILY,
      size = 6.6,
      colour = COL_TEXT
    ),

    axis.line = element_line(
      colour = "#444444",
      linewidth = 0.45
    ),
    axis.ticks = element_line(
      colour = "#444444",
      linewidth = 0.35
    ),

    panel.grid.major.y = element_line(
      colour = COL_GRID,
      linewidth = 0.28
    ),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),

    panel.spacing.x = grid::unit(8, "mm"),

    legend.position = "bottom",
    legend.box = "horizontal",
    legend.title = element_text(
      family = FONT_FAMILY,
      size = 6.9,
      colour = COL_TEXT
    ),
    legend.text = element_text(
      family = FONT_FAMILY,
      size = 6.5,
      colour = COL_TEXT
    ),
    legend.key.width = grid::unit(4.4, "mm"),
    legend.key.height = grid::unit(3.8, "mm"),
    legend.spacing.x = grid::unit(2.2, "mm"),
    legend.box.spacing = grid::unit(0.8, "mm"),
    legend.margin = margin(1, 0, 0, 0),

    plot.margin = margin(7, 10, 4, 7)
  ) +

  guides(
    fill = guide_legend(
      order = 1,
      title.position = "top",
      override.aes = list(
        shape = 21,
        size = 4.1,
        colour = "#666666",
        stroke = 0.7
      )
    )
  )

Figure3F <- p

# ------------------------------------------------------------------------------
# 9. EXPORT
# ------------------------------------------------------------------------------
WIDTH_IN <- FIG_WIDTH_MM / 25.4
HEIGHT_IN <- FIG_HEIGHT_MM / 25.4

PDF_FILE <- file.path(
  OUT_DIR,
  "Figure3F_top10_ranked_target_landscape_v9_black_borders.pdf"
)

SVG_FILE <- file.path(
  OUT_DIR,
  "Figure3F_top10_ranked_target_landscape_v9_black_borders.svg"
)

TIFF_FILE <- file.path(
  OUT_DIR,
  "Figure3F_top10_ranked_target_landscape_v9_black_borders.tiff"
)

PNG_FILE <- file.path(
  OUT_DIR,
  "Figure3F_top10_ranked_target_landscape_v9_black_borders.png"
)

ggsave(
  PDF_FILE,
  Figure3F,
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
  SVG_FILE,
  Figure3F,
  device = svglite::svglite,
  width = WIDTH_IN,
  height = HEIGHT_IN,
  units = "in",
  bg = "white"
)

ggsave(
  TIFF_FILE,
  Figure3F,
  device = ragg::agg_tiff,
  width = WIDTH_IN,
  height = HEIGHT_IN,
  units = "in",
  res = 600,
  compression = "lzw",
  background = "white"
)

ggsave(
  PNG_FILE,
  Figure3F,
  device = ragg::agg_png,
  width = WIDTH_IN,
  height = HEIGHT_IN,
  units = "in",
  res = 300,
  background = "white"
)

# ------------------------------------------------------------------------------
# 10. REPRODUCIBILITY
# ------------------------------------------------------------------------------
write.csv(
  data.frame(
    file = STEP7G_EDGES,
    md5 = unname(tools::md5sum(STEP7G_EDGES)),
    stringsAsFactors = FALSE
  ),
  file.path(OUT_DIR, "source_file_md5.csv"),
  row.names = FALSE
)

writeLines(
  capture.output(sessionInfo()),
  file.path(OUT_DIR, "sessionInfo.txt")
)

# ------------------------------------------------------------------------------
# 11. CONSOLE SUMMARY
# ------------------------------------------------------------------------------
cat("\n============================================================\n")
cat("Figure 3F v9 completed successfully.\n")
cat("Output directory:\n", OUT_DIR, "\n", sep = "")

cat("\nDisplayed Top 10 Nfkb1 targets:\n")
print(
  plot_d[
    plot_d$TF == "Nfkb1",
    c(
      "display_rank",
      "target",
      "evidence_count",
      "weight",
      "weight_percentile_display",
      "module_display",
      "shared_priority_target"
    )
  ],
  row.names = FALSE
)

cat("\nDisplayed Top 10 Nfe2l2 targets:\n")
print(
  plot_d[
    plot_d$TF == "Nfe2l2",
    c(
      "display_rank",
      "target",
      "evidence_count",
      "weight",
      "weight_percentile_display",
      "module_display",
      "shared_priority_target"
    )
  ],
  row.names = FALSE
)

cat("\nFigure rules:\n")
cat("- Step7G only; no CellChat input.\n")
cat("- Only objective Top 10 targets per regulon are displayed.\n")
cat("- Ranking = evidence_count descending, then pySCENIC weight descending.\n")
cat("- Point fill = formal M1/M2 module.\n")
cat("- Point border = shared priority target vs TF-prioritized only.\n")
cat("- Gene labels are placed in a compact far-left in-panel gutter.\n")
cat("- Stems begin immediately to the right of labels with a small fixed gap.\n")
cat("- Bubble size is reduced and figure height increased to prevent overlap.\n")
cat("- Extra right-side plotting room prevents 100% points from clipping.\n")
cat("- Point size = evidence count; this is evidence-layer count, not regulatory strength.\n")
cat("- No gene-specific highlighting is used.\n")
cat("============================================================\n")
