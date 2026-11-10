# ==============================================================================
# Figure S4 | GLIA | FINAL v6.5 | no panel subtitles
# Regulatory mechanism robustness and Nfe2l2-Tnf support
#
# DIRECT ORIGINAL-SOURCE MODE
#
# v6.5 visual-only change: panel-title subtitles removed; data/statistics/layout unchanged.
#
# Core question:
#   Are the key Figure 3 regulons robust to the RSS sensitivity definition,
#   do their temporal dynamics persist after biological-sample balancing,
#   do their regulon targets close onto the target program across
#   multiple independent evidence layers, and is the Nfe2l2-associated -> Tnf
#   relationship supported by multi-layer but NON-CAUSAL evidence?
#
# Panels:
#   A. RSS sensitivity concordance
#   B. Sample-balanced regulon timing on retained pseudotime
#   C. Target-program regulatory closure across five evidence domains
#   D. Compact Nfe2l2-Tnf multi-layer evidence
#
# IMPORTANT:
#   * Reads frozen Step7E / Step7F / Step7G / Step8D tables directly.
#   * Does NOT rerun pySCENIC, AUCell, RSS, CSI, trajectory modeling,
#     regulon-target enrichment, motif pruning, or CellChat.
#   * Reader-facing terminology follows Main Figures 3-4:
#       Shared prebranch | Alternative branch | Target program
#     Internal Arm_1 / Arm_2 labels are used only to map frozen source data.
#   * Panel B deliberately uses SAMPLE-BALANCED RETAINED-PSEUDOTIME summaries
#     rather than copying Main Figure 3D's fitted branch-progress curves.
#   * Panel C deliberately summarizes FIVE closure domains rather than
#     copying Main Figure 3E's single formal-module enrichment panel.
#   * Panel D deliberately avoids the Main Figure 4E Nfe2l2-AUCell vs Tnf
#     scatter and instead shows the underlying multi-layer evidence audit.
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
  "FigureS4_FINAL_v6_5_no_subtitles"
)

dir.create(
  OUT,
  recursive = TRUE,
  showWarnings = FALSE
)

P <- list(
  rss_primary = file.path(
    ROOT,
    "04_results",
    "step7E_GSE162610_RSS_CSI_regulon_specificity_modules",
    "tables",
    "03_RSS_primary_long.tsv"
  ),
  rss_sensitivity = file.path(
    ROOT,
    "04_results",
    "step7E_GSE162610_RSS_CSI_regulon_specificity_modules",
    "tables",
    "03c_RSS_sensitivity_informative_long.tsv"
  ),
  rss_audit = file.path(
    ROOT,
    "04_results",
    "step7E_GSE162610_RSS_CSI_regulon_specificity_modules",
    "tables",
    "03d_RSS_primary_vs_sensitivity_audit.tsv"
  ),
  rss_key_tf = file.path(
    ROOT,
    "04_results",
    "step7E_GSE162610_RSS_CSI_regulon_specificity_modules",
    "tables",
    "09_key_TF_RSS_CSI_audit.tsv"
  ),

  timing_key = file.path(
    ROOT,
    "04_results",
    "step7F_GSE162610_regulon_fate_timing",
    "tables",
    "04_key_candidate_fate_timing.tsv"
  ),
  timing_balanced_bins = file.path(
    ROOT,
    "04_results",
    "step7F_GSE162610_regulon_fate_timing",
    "tables",
    "07_key_analyte_sample_balanced_bins.tsv.gz"
  ),

  closure_timing = file.path(
    ROOT,
    "04_results",
    "step7G_GSE162610_regulon_target_program_closure",
    "tables",
    "03_timing_candidate_target_closure.tsv"
  ),

  motif_audit = file.path(
    ROOT,
    "04_results",
    "step8D_GSE162610_mechanism_closure_QC",
    "tables",
    "04_focus_TF_ligand_motif_target_audit.tsv"
  ),
  correlation_long = file.path(
    ROOT,
    "04_results",
    "step8D_GSE162610_mechanism_closure_QC",
    "tables",
    "05_focus_TF_ligand_two_sample_correlation_long.tsv"
  ),
  closure_summary = file.path(
    ROOT,
    "04_results",
    "step8D_GSE162610_mechanism_closure_QC",
    "tables",
    "06_focus_TF_ligand_closure_summary.tsv"
  )
)

missing_inputs <- names(P)[
  !vapply(
    P,
    file.exists,
    logical(1)
  )
]

if (length(missing_inputs) > 0) {
  stop(
    "Missing required Figure S4 input(s):\n",
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
  "ggplot2",
  "dplyr",
  "tidyr",
  "readr",
  "patchwork",
  "scales",
  "svglite",
  "ragg"
)

missing_pkgs <- required_pkgs[
  !vapply(
    required_pkgs,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_pkgs) > 0) {
  stop(
    "Missing required packages: ",
    paste(
      missing_pkgs,
      collapse = ", "
    ),
    "\nInstall them first, then rerun."
  )
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(patchwork)
  library(scales)
})

# ------------------------------------------------------------------------------
# 2. Frozen visual language from Main Figure 3
# ------------------------------------------------------------------------------
BASE_FAMILY <- "Arial"

COL_SHARED <- "#7A6AA6"
COL_TARGET <- "#B64B45"
COL_OTHER  <- "#7F8C8D"

COL_TEXT <- "#222222"
COL_GRID <- "#E8E8E8"
COL_GREY <- "#CFCFCF"
COL_LIGHT_GREY <- "#F2F2F2"

STATE_COLORS <- c(
  "Shared prebranch" = COL_SHARED,
  "Alternative branch" = COL_OTHER,
  "Target program" = COL_TARGET
)

CORE_TFS <- c(
  "Nfkb1",
  "Nfe2l2",
  "Relb",
  "Nfkb2"
)

CORE_REGULONS <- paste0(
  CORE_TFS,
  "(+)"
)

theme_pub <- theme_classic(
  base_family = BASE_FAMILY,
  base_size = 7.2
) +
  theme(
    plot.background = element_rect(
      fill = "white",
      color = NA
    ),
    panel.background = element_rect(
      fill = "white",
      color = NA
    ),
    axis.title = element_text(
      family = BASE_FAMILY,
      size = 7.0,
      color = COL_TEXT
    ),
    axis.text = element_text(
      family = BASE_FAMILY,
      size = 5.9,
      color = COL_TEXT
    ),
    axis.line = element_line(
      linewidth = 0.36,
      color = COL_TEXT
    ),
    axis.ticks = element_line(
      linewidth = 0.30,
      color = COL_TEXT
    ),
    legend.title = element_text(
      family = BASE_FAMILY,
      size = 5.8,
      color = COL_TEXT
    ),
    legend.text = element_text(
      family = BASE_FAMILY,
      size = 5.2,
      color = COL_TEXT
    ),
    plot.title = element_text(
      family = BASE_FAMILY,
      size = 7.1,
      face = "bold",
      hjust = 0,
      color = COL_TEXT
    ),
    plot.subtitle = element_blank(),
    plot.tag = element_text(
      family = BASE_FAMILY,
      face = "bold",
      size = 11.5,
      color = "black"
    ),
    plot.tag.position = c(
      -0.035,
      1.03
    ),
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
    file.path(
      OUT,
      paste0(
        stem,
        ".pdf"
      )
    ),
    p,
    device = cairo_pdf,
    width = width_mm,
    height = height_mm,
    units = "mm",
    bg = "white"
  )

  ggsave(
    file.path(
      OUT,
      paste0(
        stem,
        ".svg"
      )
    ),
    p,
    device = svglite::svglite,
    width = width_mm,
    height = height_mm,
    units = "mm",
    bg = "white"
  )

  ggsave(
    file.path(
      OUT,
      paste0(
        stem,
        ".tiff"
      )
    ),
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
    file.path(
      OUT,
      paste0(
        stem,
        "_review.png"
      )
    ),
    p,
    device = ragg::agg_png,
    width = width_mm,
    height = height_mm,
    units = "mm",
    dpi = 220,
    background = "white"
  )
}

state_display <- function(x) {
  dplyr::recode(
    as.character(x),
    "Shared_prebranch" = "Shared prebranch",
    "Arm_1" = "Alternative branch",
    "Arm_2" = "Target program",
    .default = as.character(x)
  )
}

# ------------------------------------------------------------------------------
# 4. Input audit
# ------------------------------------------------------------------------------
input_manifest <- data.frame(
  key = names(P),
  path = unlist(
    P,
    use.names = FALSE
  ),
  size_MB = round(
    file.info(
      unlist(P)
    )$size / 1024^2,
    4
  ),
  md5 = unname(
    tools::md5sum(
      unlist(P)
    )
  ),
  stringsAsFactors = FALSE
)

write.csv(
  input_manifest,
  file.path(
    OUT,
    "01_FigureS4_input_manifest.csv"
  ),
  row.names = FALSE
)

# ==============================================================================
# PANEL A
# RSS sensitivity concordance
# ==============================================================================

message("Reading Step7E RSS primary/sensitivity tables...")

rss_primary <- readr::read_tsv(
  P$rss_primary,
  show_col_types = FALSE
)

rss_sens <- readr::read_tsv(
  P$rss_sensitivity,
  show_col_types = FALSE
)

rss_audit <- readr::read_tsv(
  P$rss_audit,
  show_col_types = FALSE
)

rss_key <- readr::read_tsv(
  P$rss_key_tf,
  show_col_types = FALSE
)

req_audit <- c(
  "state",
  "primary_regulons",
  "sensitivity_regulons",
  "Spearman_rank_correlation_on_primary",
  "top20_overlap_n",
  "top20_overlap_fraction"
)

if (
  length(
    setdiff(
      req_audit,
      colnames(rss_audit)
    )
  ) > 0
) {
  stop(
    "Step7E RSS audit schema is incomplete."
  )
}

# Frozen audit expectations:
# 215 primary regulons vs 284 informative-sensitivity regulons;
# common-set rank correlation = 1 and top20 overlap = 20/20 in each state.
if (
  !all(
    rss_audit$primary_regulons == 215
  ) ||
  !all(
    rss_audit$sensitivity_regulons == 284
  ) ||
  !all(
    abs(
      rss_audit$Spearman_rank_correlation_on_primary -
        1
    ) < 1e-8
  ) ||
  !all(
    rss_audit$top20_overlap_n == 20
  )
) {
  warning(
    "Frozen Step7E RSS sensitivity summary differs from the expected ",
    "215-vs-284 / rho=1 / top20=20/20 audit. Inspect exported CSV."
  )
}

rss_audit_plot <- rss_audit %>%
  mutate(
    state_display = factor(
      state_display(state),
      levels = c(
        "Shared prebranch",
        "Alternative branch",
        "Target program"
      )
    )
  )

rss_metric_long <- bind_rows(
  rss_audit_plot %>%
    transmute(
      state_display,
      metric = "Rank\nr_s",
      label = sprintf(
        "%.2f",
        Spearman_rank_correlation_on_primary
      )
    ),
  rss_audit_plot %>%
    transmute(
      state_display,
      metric = "Top-20\noverlap",
      label = paste0(
        top20_overlap_n,
        "/20"
      )
    )
) %>%
  mutate(
    metric = factor(
      metric,
      levels = c(
        "Rank\nr_s",
        "Top-20\noverlap"
      )
    )
  )

pA_summary <- ggplot(
  rss_metric_long,
  aes(
    x = metric,
    y = state_display
  )
) +
  geom_tile(
    aes(
      fill = state_display
    ),
    color = "white",
    linewidth = 1.0,
    alpha = 0.20
  ) +
  geom_text(
    aes(
      label = label,
      color = state_display
    ),
    family = BASE_FAMILY,
    fontface = "bold",
    size = 2.55
  ) +
  scale_fill_manual(
    values = STATE_COLORS,
    guide = "none"
  ) +
  scale_color_manual(
    values = STATE_COLORS,
    guide = "none"
  ) +
  labs(
    tag = "A",
    title = "Global sensitivity",
    subtitle = NULL,
    x = NULL,
    y = NULL
  ) +
  theme_minimal(
    base_family = BASE_FAMILY,
    base_size = 7
  ) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(
      size = 4.7,
      color = COL_TEXT,
      lineheight = 0.90
    ),
    axis.text.y = element_text(
      size = 5.6,
      color = COL_TEXT
    ),
    plot.title = element_text(
      size = 6.4,
      face = "bold",
      hjust = 0
    ),
    plot.subtitle = element_blank(),
    plot.tag = element_text(
      family = BASE_FAMILY,
      face = "bold",
      size = 11.5,
      color = "black"
    ),
    plot.tag.position = c(
      -0.10,
      1.06
    ),
    plot.margin = margin(
      3.0, 2.0, 2.5, 3.8,
      unit = "mm"
    )
  )

rss_key_rank <- rss_primary %>%
  filter(
    state == "Arm_2",
    TF %in% CORE_TFS
  ) %>%
  select(
    regulon,
    TF,
    primary_rank = RSS_rank_in_state
  ) %>%
  inner_join(
    rss_sens %>%
      filter(
        state == "Arm_2",
        TF %in% CORE_TFS
      ) %>%
      select(
        regulon,
        TF,
        sensitivity_rank = RSS_rank_in_state
      ),
    by = c(
      "regulon",
      "TF"
    )
  ) %>%
  mutate(
    TF = factor(
      TF,
      levels = rev(
        CORE_TFS
      )
    ),
    rank_label = paste0(
      primary_rank,
      " \u2192 ",
      sensitivity_rank
    )
  )

if (
  nrow(rss_key_rank) != length(CORE_TFS)
) {
  stop(
    "Could not recover all four core Target-program regulons from Step7E."
  )
}

rss_key_long <- rss_key_rank %>%
  select(
    TF,
    primary_rank,
    sensitivity_rank
  ) %>%
  pivot_longer(
    cols = c(
      primary_rank,
      sensitivity_rank
    ),
    names_to = "analysis",
    values_to = "rank"
  ) %>%
  mutate(
    analysis = factor(
      analysis,
      levels = c(
        "primary_rank",
        "sensitivity_rank"
      ),
      labels = c(
        "Primary",
        "Sensitivity"
      )
    )
  )

pA_rank <- ggplot() +
  geom_segment(
    data = rss_key_rank,
    aes(
      x = primary_rank,
      xend = sensitivity_rank,
      y = TF,
      yend = TF
    ),
    linewidth = 0.70,
    color = "#AFAFAF"
  ) +
  geom_point(
    data = rss_key_long,
    aes(
      x = rank,
      y = TF,
      shape = analysis,
      fill = analysis
    ),
    size = 2.15,
    stroke = 0.55,
    color = COL_TARGET
  ) +

  scale_x_log10(
    breaks = c(
      1, 5, 25, 100
    ),
    limits = c(
      0.85,
      165
    )
  ) +
  scale_shape_manual(
    values = c(
      "Primary" = 21,
      "Sensitivity" = 16
    )
  ) +
  scale_fill_manual(
    values = c(
      "Primary" = "white",
      "Sensitivity" = COL_TARGET
    )
  ) +
  scale_y_discrete(
    labels = c(
      "Nfkb1" = "Nfkb1   1→1",
      "Nfe2l2" = "Nfe2l2   5→5",
      "Relb" = "Relb   26→26",
      "Nfkb2" = "Nfkb2   128→129"
    )
  ) +
  labs(
    title = "Core rank robustness",
    subtitle = NULL,
    x = "Regulon rank (log scale)",
    y = NULL,
    shape = NULL,
    fill = NULL
  ) +
  theme_pub +
  theme(
    plot.title = element_text(
      size = 5.9,
      face = "bold"
    ),
    plot.subtitle = element_blank(),
    axis.text.y = element_text(
      size = 5.15
    ),
    axis.text.x = element_text(
      size = 5.0
    ),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.key.width = grid::unit(
      2.8,
      "mm"
    ),
    plot.margin = margin(
      2.5, 3.5, 2.5, 4.5,
      unit = "mm"
    )
  ) +
  guides(
    fill = "none",
    shape = guide_legend(
      nrow = 1,
      override.aes = list(
        size = 2.0
      )
    )
  )

panelA <- patchwork::wrap_plots(
  pA_summary,
  pA_rank,
  ncol = 2,
  widths = c(
    0.37,
    0.63
  )
) &
  theme(
    plot.margin = margin(
      3.0, 3.2, 2.8, 3.2,
      unit = "mm"
    )
  )

# ==============================================================================
# PANEL B
# Sample-balanced regulon timing
# ==============================================================================

message("Reading Step7F sample-balanced bins...")

timing_key <- readr::read_tsv(
  P$timing_key,
  show_col_types = FALSE
)

timing_bins <- readr::read_tsv(
  P$timing_balanced_bins,
  show_col_types = FALSE
)

core_timing_key <- timing_key %>%
  filter(
    TF %in% CORE_TFS
  )

if (
  nrow(core_timing_key) !=
    length(CORE_TFS)
) {
  stop(
    "Could not recover the four core Step7F regulon timing rows."
  )
}

# All four final core regulons must have qualified onset at branch origin and
# reproducible 7 dpi branch direction in the frozen analysis.
if (
  !all(
    core_timing_key$onset_progress_qualified == 0
  ) ||
  !all(
    core_timing_key$seven_dpi_direction_consistent
  )
) {
  warning(
    "Core Step7F timing summary differs from the expected ",
    "branch-origin onset + 7 dpi direction-consistent result."
  )
}

timing_plot <- timing_bins %>%
  filter(
    analyte %in% CORE_REGULONS,
    analyte_type == "regulon_AUC",
    axis == "retained_pseudotime"
  ) %>%
  mutate(
    TF = sub(
      "\\(\\+\\)$",
      "",
      analyte
    ),
    TF = factor(
      TF,
      levels = CORE_TFS
    ),
    state_display = factor(
      state_display(state),
      levels = c(
        "Shared prebranch",
        "Alternative branch",
        "Target program"
      )
    )
  )

# Frozen balanced-pseudotime audit:
# Shared prebranch = 2 biological samples; postbranch arms = 4 samples.
sample_audit_B <- timing_plot %>%
  group_by(
    TF,
    state_display
  ) %>%
  summarise(
    min_n_samples = min(
      n_samples,
      na.rm = TRUE
    ),
    max_n_samples = max(
      n_samples,
      na.rm = TRUE
    ),
    n_bins = n(),
    .groups = "drop"
  )

pB <- ggplot(
  timing_plot,
  aes(
    x = axis_value,
    y = balanced_mean,
    color = state_display,
    group = state_display
  )
) +
  geom_line(
    linewidth = 0.72,
    alpha = 0.95
  ) +
  geom_point(
    size = 1.10,
    alpha = 0.96
  ) +
  facet_wrap(
    ~TF,
    ncol = 2,
    scales = "free_y"
  ) +
  scale_color_manual(
    values = STATE_COLORS,
    drop = FALSE
  ) +
  labs(
    tag = "B",
    title = "Sample-balanced regulon timing",
    subtitle = NULL,
    x = "Retained pseudotime",
    y = "Sample-balanced regulon AUC",
    color = NULL
  ) +
  theme_pub +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(
      family = BASE_FAMILY,
      face = "bold",
      size = 5.9,
      hjust = 0
    ),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.key.width = grid::unit(
      3.0,
      "mm"
    ),
    axis.text = element_text(
      size = 4.9
    ),
    axis.title = element_text(
      size = 5.8
    ),
    panel.spacing = grid::unit(
      2.6,
      "mm"
    ),
    plot.margin = margin(
      3.0, 2.5, 2.8, 3.8,
      unit = "mm"
    )
  ) +
  guides(
    color = guide_legend(
      nrow = 1,
      override.aes = list(
        linewidth = 1.1,
        size = 1.8
      )
    )
  )

# ==============================================================================
# PANEL C
# Target-program regulatory closure
# ==============================================================================

message("Reading Step7G Target-program closure table...")

closure_timing <- readr::read_tsv(
  P$closure_timing,
  show_col_types = FALSE
)

core_closure <- closure_timing %>%
  filter(
    TF %in% CORE_TFS
  )

if (
  nrow(core_closure) !=
    length(CORE_TFS)
) {
  stop(
    "Could not recover all four core Step7G closure rows."
  )
}

make_domain <- function(
  df,
  domain,
  overlap_col,
  fold_col,
  q_col,
  support_col
) {
  df %>%
    transmute(
      TF,
      domain = domain,
      overlap_n = .data[[overlap_col]],
      fold_enrichment = .data[[fold_col]],
      q_value = .data[[q_col]],
      support = as.logical(
        .data[[support_col]]
      ),
      target_closure_strength
    )
}

closure_long <- bind_rows(
  make_domain(
    core_closure,
    "Progressive\nTarget",
    "progressive_overlap_n",
    "progressive_fold_enrichment",
    "progressive_q",
    "progressive_support"
  ),
  make_domain(
    core_closure,
    "Strong-progressive\nTarget",
    "strong_progressive_overlap_n",
    "strong_progressive_fold_enrichment",
    "strong_progressive_q",
    "strong_progressive_support"
  ),
  make_domain(
    core_closure,
    "High-conf.\nmarkers",
    "highconf_marker_overlap_n",
    "highconf_marker_fold_enrichment",
    "highconf_marker_q",
    "highconf_marker_support"
  ),
  make_domain(
    core_closure,
    "Strict PB\nmarkers",
    "strict_marker_overlap_n",
    "strict_marker_fold_enrichment",
    "strict_marker_q",
    "strict_marker_support"
  ),
  make_domain(
    core_closure,
    "Target-terminal\nmodules",
    "formal_module_overlap_n",
    "formal_module_fold_enrichment",
    "formal_module_q",
    "formal_module_support"
  )
) %>%
  mutate(
    TF = factor(
      TF,
      levels = rev(
        CORE_TFS
      )
    ),
    domain = factor(
      domain,
      levels = c(
        "Progressive\nTarget",
        "Strong-progressive\nTarget",
        "High-conf.\nmarkers",
        "Strict PB\nmarkers",
        "Target-terminal\nmodules"
      )
    ),
    log2_fold = log2(
      pmax(
        fold_enrichment,
        1e-6
      )
    ),
    fill_value = ifelse(
      support,
      log2_fold,
      NA_real_
    ),
    support_display = ifelse(
      support,
      "Supported",
      "Not significant"
    )
  )

closure_strength <- core_closure %>%
  transmute(
    TF = factor(
      TF,
      levels = rev(
        CORE_TFS
      )
    ),
    target_closure_strength
  )

# Preserve the exact final hierarchy:
# Nfkb1/Nfe2l2 = strong closure; Relb/Nfkb2 = weak closure.
expected_strength <- c(
  "Nfkb1" = "strong",
  "Nfe2l2" = "strong",
  "Relb" = "weak",
  "Nfkb2" = "weak"
)

observed_strength <- setNames(
  as.character(
    core_closure$target_closure_strength
  ),
  core_closure$TF
)

if (
  !all(
    observed_strength[
      names(
        expected_strength
      )
    ] ==
      expected_strength
  )
) {
  warning(
    "Step7G closure-strength hierarchy differs from the frozen ",
    "Nfkb1/Nfe2l2 strong; Relb/Nfkb2 weak expectation."
  )
}

tf_label_map <- setNames(
  paste0(
    CORE_TFS,
    "  (",
    expected_strength[
      CORE_TFS
    ],
    ")"
  ),
  CORE_TFS
)

pC <- ggplot(
  closure_long,
  aes(
    x = domain,
    y = TF
  )
) +
  geom_point(
    aes(
      size = overlap_n,
      fill = fill_value,
      color = support_display
    ),
    shape = 21,
    stroke = 0.65
  ) +
  scale_fill_gradient(
    low = "#F3E4E2",
    high = COL_TARGET,
    na.value = "white",
    name = "log2 fold\nenrichment"
  ) +
  scale_color_manual(
    values = c(
      "Supported" = COL_TARGET,
      "Not significant" = "#B8B8B8"
    ),
    guide = "none"
  ) +
  scale_size_continuous(
    range = c(
      2.0,
      7.2
    ),
    breaks = c(
      5, 20, 50, 100
    ),
    name = "Overlap genes"
  ) +
  scale_y_discrete(
    labels = function(x) {
      unname(
        tf_label_map[
          as.character(x)
        ]
      )
    }
  ) +
  labs(
    tag = "C",
    title = "Target-program regulatory closure",
    subtitle = NULL,
    x = NULL,
    y = NULL
  ) +
  theme_pub +
  theme(
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(
      size = 4.65,
      lineheight = 0.90,
      angle = 0,
      hjust = 0.5
    ),
    axis.text.y = element_text(
      size = 5.6
    ),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.direction = "horizontal",
    legend.spacing.x = grid::unit(1.4, "mm"),
    panel.grid.major.x = element_line(
      color = "#EFEFEF",
      linewidth = 0.35
    ),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.text = element_text(
      size = 4.7
    ),
    legend.title = element_text(
      size = 4.9
    ),
    plot.subtitle = element_blank(),
    plot.margin = margin(
      3.0, 2.8, 3.2, 3.8,
      unit = "mm"
    )
  ) +
  guides(
    size = guide_legend(
      order = 1,
      nrow = 1
    ),
    fill = guide_colorbar(
      order = 2,
      barwidth = grid::unit(
        16,
        "mm"
      ),
      barheight = grid::unit(
        2.2,
        "mm"
      )
    )
  )

# ==============================================================================
# PANEL D
# Compact Nfe2l2-Tnf evidence
# ==============================================================================

message("Reading Step8D Nfe2l2-Tnf evidence tables...")

motif_audit <- readr::read_tsv(
  P$motif_audit,
  show_col_types = FALSE
)

correlation_long <- readr::read_tsv(
  P$correlation_long,
  show_col_types = FALSE
)

closure_summary <- readr::read_tsv(
  P$closure_summary,
  show_col_types = FALSE
)

focus_pairs <- c(
  "Nfe2l2",
  "Nfkb1"
)

motif_tnf <- motif_audit %>%
  filter(
    ligand == "Tnf",
    TF %in% focus_pairs
  )

closure_tnf <- closure_summary %>%
  filter(
    ligand == "Tnf",
    TF %in% focus_pairs
  )

corr_tnf <- correlation_long %>%
  filter(
    ligand == "Tnf",
    TF %in% focus_pairs,
    analysis_level == "WITHIN_7DPI_SAMPLE"
  )

if (
  nrow(motif_tnf) != 2 ||
  nrow(closure_tnf) != 2 ||
  nrow(corr_tnf) != 4
) {
  stop(
    "Step8D Nfe2l2/Nfkb1-Tnf evidence rows are incomplete."
  )
}

# Exact interpretation boundary:
# Nfe2l2-Tnf = direct motif-supported candidate with replicated association.
# Nfkb1-Tnf = indirect replicated association; no direct motif edge.
nfe_row <- closure_tnf %>%
  filter(
    TF == "Nfe2l2"
  )

nfkb1_row <- closure_tnf %>%
  filter(
    TF == "Nfkb1"
  )

if (
  !isTRUE(
    nfe_row$motif_supported_direct_target[[1]]
  ) ||
  isTRUE(
    nfkb1_row$motif_supported_direct_target[[1]]
  ) ||
  !isTRUE(
    nfe_row$two_sample_state_concordance[[1]]
  ) ||
  !isTRUE(
    nfe_row$two_sample_positive_correlation[[1]]
  )
) {
  stop(
    "Frozen Step8D Nfe2l2-Tnf interpretation does not match ",
    "the expected motif-supported / replicated-association audit."
  )
}

target_counts <- rss_key %>%
  filter(
    TF %in% focus_pairs
  ) %>%
  select(
    TF,
    n_targets
  )

rho_wide <- corr_tnf %>%
  select(
    TF,
    sample_id,
    spearman_rho
  ) %>%
  mutate(
    sample_id = dplyr::recode(
      sample_id,
      "7dpi_sample1" = "rho_sample1",
      "7dpi_sample2" = "rho_sample2"
    )
  ) %>%
  pivot_wider(
    names_from = sample_id,
    values_from = spearman_rho
  )

d_summary <- closure_tnf %>%
  select(
    TF,
    motif_supported_direct_target,
    target_weight,
    target_weight_rank_within_regulon,
    target_weight_percentile_within_regulon,
    n_samples_state_concordant,
    n_samples_positive_correlation,
    regulatory_closure_status,
    manuscript_interpretation
  ) %>%
  left_join(
    target_counts,
    by = "TF"
  ) %>%
  left_join(
    rho_wide,
    by = "TF"
  ) %>%
  mutate(
    pair = ifelse(
      TF == "Nfe2l2",
      "Nfe2l2 \u2192 Tnf",
      "Nfkb1 \u2192 Tnf\n(context)"
    ),
    pair = factor(
      pair,
      levels = c(
        "Nfkb1 \u2192 Tnf\n(context)",
        "Nfe2l2 \u2192 Tnf"
      )
    ),
    direct_motif_label = ifelse(
      motif_supported_direct_target,
      "Yes",
      "No"
    ),
    weight_label = ifelse(
      motif_supported_direct_target,
      paste0(
        round(
          target_weight_rank_within_regulon
        ),
        "/",
        n_targets,
        "\n",
        sprintf(
          "%.1f",
          100 *
            target_weight_percentile_within_regulon
        ),
        "th pct"
      ),
      "\u2014"
    ),
    rho1_label = paste0(
      "r_s = ",
      sprintf(
        "%.3f",
        rho_sample1
      )
    ),
    rho2_label = paste0(
      "r_s = ",
      sprintf(
        "%.3f",
        rho_sample2
      )
    ),
    concordance_label = paste0(
      n_samples_state_concordant,
      "/2 state\n",
      n_samples_positive_correlation,
      "/2 +r"
    ),
    tier_label = ifelse(
      motif_supported_direct_target,
      "Motif-supported\n+ replicated",
      "Replicated\nassociation"
    ),
    row_type = ifelse(
      TF == "Nfe2l2",
      "Nfe2l2-Tnf",
      "Nfkb1-Tnf"
    )
  )

evidence_cells <- bind_rows(
  d_summary %>%
    transmute(
      pair,
      row_type,
      metric = "Motif\ntarget",
      label = direct_motif_label,
      positive = motif_supported_direct_target
    ),
  d_summary %>%
    transmute(
      pair,
      row_type,
      metric = "Target\nrank",
      label = weight_label,
      positive = motif_supported_direct_target
    ),
  d_summary %>%
    transmute(
      pair,
      row_type,
      metric = "7 dpi S1\nassociation",
      label = rho1_label,
      positive = rho_sample1 > 0
    ),
  d_summary %>%
    transmute(
      pair,
      row_type,
      metric = "7 dpi S2\nassociation",
      label = rho2_label,
      positive = rho_sample2 > 0
    ),
  d_summary %>%
    transmute(
      pair,
      row_type,
      metric = "Replicate\nsupport",
      label = concordance_label,
      positive =
        n_samples_state_concordant == 2 &
        n_samples_positive_correlation == 2
    ),
  d_summary %>%
    transmute(
      pair,
      row_type,
      metric = "Evidence tier",
      label = tier_label,
      positive = TF == "Nfe2l2"
    )
) %>%
  mutate(
    metric = factor(
      metric,
      levels = c(
        "Motif\ntarget",
        "Target\nrank",
        "7 dpi S1\nassociation",
        "7 dpi S2\nassociation",
        "Replicate\nsupport",
        "Evidence tier"
      )
    )
  )

pD <- ggplot(
  evidence_cells,
  aes(
    x = metric,
    y = pair
  )
) +
  geom_tile(
    aes(
      fill = row_type
    ),
    color = "white",
    linewidth = 1.2
  ) +
  geom_text(
    aes(
      label = label,
      color = ifelse(
        positive,
        "positive",
        "neutral"
      )
    ),
    family = BASE_FAMILY,
    size = 1.92,
    lineheight = 0.90
  ) +
  scale_fill_manual(
    values = c(
      "Nfe2l2-Tnf" = "#F5E5E3",
      "Nfkb1-Tnf" = "#F1F1F1"
    ),
    guide = "none"
  ) +
  scale_color_manual(
    values = c(
      "positive" = COL_TARGET,
      "neutral" = "#777777"
    ),
    guide = "none"
  ) +
  labs(
    tag = "D",
    title = "Nfe2l2\u2013Tnf multi-layer support",
    subtitle = NULL,
    x = NULL,
    y = NULL
  ) +
  theme_minimal(
    base_family = BASE_FAMILY,
    base_size = 7
  ) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(
      size = 4.45,
      color = COL_TEXT,
      lineheight = 0.90
    ),
    axis.text.y = element_text(
      size = 5.45,
      color = COL_TEXT
    ),
    plot.title = element_text(
      family = BASE_FAMILY,
      size = 7.1,
      face = "bold",
      hjust = 0,
      color = COL_TEXT
    ),
    plot.subtitle = element_blank(),
    plot.tag = element_text(
      family = BASE_FAMILY,
      face = "bold",
      size = 11.5
    ),
    plot.tag.position = c(
      -0.035,
      1.03
    ),
    plot.margin = margin(
      3.0, 3.0, 3.0, 3.8,
      unit = "mm"
    )
  )

# ==============================================================================
# AUDITS
# ==============================================================================

write.csv(
  rss_audit_plot,
  file.path(
    OUT,
    "02_PanelA_RSS_sensitivity_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  rss_key_rank,
  file.path(
    OUT,
    "02B_PanelA_core_Target_regulon_ranks.csv"
  ),
  row.names = FALSE
)

write.csv(
  sample_audit_B,
  file.path(
    OUT,
    "03_PanelB_sample_balancing_audit.csv"
  ),
  row.names = FALSE
)

readr::write_tsv(
  timing_plot,
  file.path(
    OUT,
    "03B_PanelB_plot_data.tsv"
  )
)

write.csv(
  closure_long,
  file.path(
    OUT,
    "04_PanelC_closure_matrix.csv"
  ),
  row.names = FALSE
)

write.csv(
  d_summary,
  file.path(
    OUT,
    "05_PanelD_Nfe2l2_Tnf_evidence_summary.csv"
  ),
  row.names = FALSE
)

# Key final data audit
nfe_tnf <- d_summary %>%
  filter(
    TF == "Nfe2l2"
  )

audit <- data.frame(
  item = c(
    "RSS_primary_regulons",
    "RSS_sensitivity_regulons",
    "RSS_states_with_rank_r_s_1",
    "RSS_states_with_top20_overlap_20",
    "core_regulons",
    "core_regulons_branch_origin_onset",
    "Nfkb1_target_closure_strength",
    "Nfe2l2_target_closure_strength",
    "Relb_target_closure_strength",
    "Nfkb2_target_closure_strength",
    "Nfe2l2_Tnf_direct_motif_target",
    "Nfe2l2_Tnf_target_weight_rank",
    "Nfe2l2_Tnf_target_weight_percentile",
    "Nfe2l2_Tnf_7dpi_sample1_r_s",
    "Nfe2l2_Tnf_7dpi_sample2_r_s",
    "Nfe2l2_Tnf_two_sample_state_concordance",
    "Nfe2l2_Tnf_two_sample_positive_correlation"
  ),
  value = as.character(
    c(
      unique(
        rss_audit$primary_regulons
      )[1],
      unique(
        rss_audit$sensitivity_regulons
      )[1],
      sum(
        abs(
          rss_audit$Spearman_rank_correlation_on_primary -
            1
        ) < 1e-8
      ),
      sum(
        rss_audit$top20_overlap_n == 20
      ),
      length(
        CORE_TFS
      ),
      sum(
        core_timing_key$onset_progress_qualified == 0
      ),
      observed_strength[["Nfkb1"]],
      observed_strength[["Nfe2l2"]],
      observed_strength[["Relb"]],
      observed_strength[["Nfkb2"]],
      nfe_tnf$motif_supported_direct_target[[1]],
      nfe_tnf$target_weight_rank_within_regulon[[1]],
      nfe_tnf$target_weight_percentile_within_regulon[[1]],
      nfe_tnf$rho_sample1[[1]],
      nfe_tnf$rho_sample2[[1]],
      nfe_tnf$n_samples_state_concordant[[1]],
      nfe_tnf$n_samples_positive_correlation[[1]]
    )
  ),
  stringsAsFactors = FALSE
)

write.csv(
  audit,
  file.path(
    OUT,
    "00_FigureS4_data_audit.csv"
  ),
  row.names = FALSE
)

# ==============================================================================
# FINAL 2x2 LAYOUT
# ==============================================================================

panelA_fixed <- patchwork::wrap_elements(
  full = panelA,
  clip = FALSE
)

panelB_fixed <- patchwork::wrap_elements(
  full = pB,
  clip = FALSE
)

panelC_fixed <- patchwork::wrap_elements(
  full = pC,
  clip = FALSE
)

panelD_fixed <- patchwork::wrap_elements(
  full = pD,
  clip = FALSE
)

row1 <- patchwork::wrap_plots(
  panelA_fixed,
  panelB_fixed,
  ncol = 2,
  widths = c(
    0.95,
    1.05
  )
)

row2 <- patchwork::wrap_plots(
  panelC_fixed,
  panelD_fixed,
  ncol = 2,
  widths = c(
    0.96,
    1.04
  )
)

figureS4 <- patchwork::wrap_plots(
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
    0.94,
    1.06
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

# ------------------------------------------------------------------------------
# 6. Export
# ------------------------------------------------------------------------------
message("Saving Figure S4 panel A...")
save_plot_quad(
  panelA,
  "FigureS4A_RSS_sensitivity_concordance",
  width_mm = 94,
  height_mm = 78
)

message("Saving Figure S4 panel B...")
save_plot_quad(
  pB,
  "FigureS4B_sample_balanced_regulon_timing",
  width_mm = 100,
  height_mm = 78
)

message("Saving Figure S4 panel C...")
save_plot_quad(
  pC,
  "FigureS4C_Target_program_regulatory_closure",
  width_mm = 94,
  height_mm = 86
)

message("Saving Figure S4 panel D...")
save_plot_quad(
  pD,
  "FigureS4D_Nfe2l2_Tnf_multilayer_support",
  width_mm = 104,
  height_mm = 86
)

message("Saving complete Figure S4...")
save_plot_quad(
  figureS4,
  "FigureS4_GLIA_FINAL_v6_5_no_panel_subtitles",
  width_mm = 190,
  height_mm = 188
)

# ------------------------------------------------------------------------------
# 7. Run summary
# ------------------------------------------------------------------------------
nfe_s1 <- nfe_tnf$rho_sample1[[1]]
nfe_s2 <- nfe_tnf$rho_sample2[[1]]

run_summary <- c(
  "Figure S4 | GLIA | FINAL v6.5 | no panel subtitles",
  "==============================================",
  "",
  paste0(
    "Project root: ",
    ROOT
  ),
  paste0(
    "Output: ",
    OUT
  ),
  "",
  "Panel A — RSS sensitivity concordance",
  paste0(
    "  Primary regulons = ",
    unique(
      rss_audit$primary_regulons
    )[1],
    "; informative sensitivity set = ",
    unique(
      rss_audit$sensitivity_regulons
    )[1],
    "."
  ),
  "  Common-set rank concordance is 1.00 and top-20 overlap is 20/20 in all three states.",
  "  Core Target-program rank labels are shown as primary -> full sensitivity-set rank.",
  "",
  "Panel B — sample-balanced regulon timing",
  "  Uses retained pseudotime sample-balanced means only; no Step7F fitted branch-progress curve is copied.",
  "  Nfkb1, Nfe2l2, Relb and Nfkb2 all have qualified onset at the branch origin and consistent 7 dpi Target direction.",
  "",
  "Panel C — Target-program regulatory closure",
  "  Five evidence domains are summarized: progressive, strong-progressive, high-confidence markers, strict pseudobulk markers, and target-terminal modules.",
  paste0(
    "  Closure strength: Nfkb1=",
    observed_strength[["Nfkb1"]],
    "; Nfe2l2=",
    observed_strength[["Nfe2l2"]],
    "; Relb=",
    observed_strength[["Relb"]],
    "; Nfkb2=",
    observed_strength[["Nfkb2"]],
    "."
  ),
  "",
  "Panel D — compact Nfe2l2-Tnf evidence",
  paste0(
    "  Nfe2l2(+) -> Tnf is motif-supported; Tnf target rank = ",
    round(
      nfe_tnf$target_weight_rank_within_regulon[[1]]
    ),
    "/",
    nfe_tnf$n_targets[[1]],
    " (",
    sprintf(
      "%.1f",
      100 *
        nfe_tnf$target_weight_percentile_within_regulon[[1]]
    ),
    "th percentile)."
  ),
  paste0(
    "  Within-sample Spearman r_s: 7 dpi sample 1 = ",
    sprintf(
      "%.3f",
      nfe_s1
    ),
    "; 7 dpi sample 2 = ",
    sprintf(
      "%.3f",
      nfe_s2
    ),
    "."
  ),
  "  Both 7 dpi samples show the same Target-state direction and positive within-sample association.",
  "  Interpretation is associative / motif-supported candidate evidence, not perturbational proof of direct causality.",
  "",
  "No upstream analysis was rerun."
)

writeLines(
  run_summary,
  con = file.path(
    OUT,
    "FigureS4_GLIA_FINAL_v6_5_no_panel_subtitles_run_summary.txt"
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
      "06_warnings_after_run.txt"
    )
  )
}

message("============================================================")
message("Figure S4 FINAL v6.5 (no panel subtitles) completed.")
message("Output: ", OUT)
message("Main review file: FigureS4_GLIA_FINAL_v6_5_no_panel_subtitles_review.png")
message("Main vector file: FigureS4_GLIA_FINAL_v6_5_no_panel_subtitles.pdf / .svg")
message("============================================================")
