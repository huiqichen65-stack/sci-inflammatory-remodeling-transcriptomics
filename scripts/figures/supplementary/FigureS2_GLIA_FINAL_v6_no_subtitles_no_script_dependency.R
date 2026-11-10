# ==============================================================================
# Figure S2 | GLIA | FINAL v4 FINAL
# Resident microglial state and Activated/Early-DAM transition robustness
#
# DIRECT ORIGINAL-SOURCE MODE
#
# Core question:
#   Are the resident microglial state definitions program-supported, represented
#   across biological samples/timepoints, and is the Activated/Early-DAM
#   transition robust to cell-cycle exclusion and pseudobulk analysis?
#
# Panels:
#   A. State-definition evidence:
#      key marker DotPlot + compact AUCell state-program heatmap
#   B. State representation across biological samples/timepoints
#   C. Activated/Early-DAM transition robustness:
#      sample-level temporal transition fraction + 3 dpi replicate validation
#   D. Non-cycling and pseudobulk control:
#      refined-state retention + cell-vs-pseudobulk concordance summary
#
# IMPORTANT:
#   * Reads frozen Step3B/Step3C RDS and Step3E tables only.
#   * Source .R scripts are NOT required for Figure S2 rendering.
#   * Does NOT rerun QC, Harmony, clustering, marker discovery, AUCell,
#     transition mapping, DEG, pseudobulk, or trajectory analysis.
#   * Resident-state display classes use the EXACT Figure 1 collapse:
#       Step2 final cluster 0   -> Canonical resident
#       Step2 final clusters 3/4 -> Activated resident
#       Step2 final cluster 9   -> IFN-responsive resident
#       Step2 final clusters 5/7 -> Cycling resident
#   * Biological sample is the preferred display/replicate unit.
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
  "FigureS2_FINAL_v6_no_subtitles"
)

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

message("Project root: ", ROOT)
message("Output dir:   ", OUT)

P <- list(
  step3B_rds = file.path(
    ROOT,
    "04_results",
    "step3B_GSE162610_state_scoring",
    "GSE162610_step3B_clean_state_scored.rds"
  ),
  step3C_rds = file.path(
    ROOT,
    "04_results",
    "step3C_GSE162610_activated_DAM_transition_mapping",
    "GSE162610_step3C_transition_mapped.rds"
  ),

  step3B_program_availability = file.path(
    ROOT,
    "04_results",
    "step3B_GSE162610_state_scoring",
    "tables",
    "03_state_program_gene_set_availability.csv"
  ),
  step3B_auc_cluster = file.path(
    ROOT,
    "04_results",
    "step3B_GSE162610_state_scoring",
    "tables",
    "04_AUCell_state_program_cluster_summary.csv"
  ),
  step3B_time = file.path(
    ROOT,
    "04_results",
    "step3B_GSE162610_state_scoring",
    "tables",
    "08_cluster_by_timepoint_after_cleanup.csv"
  ),
  step3B_sample = file.path(
    ROOT,
    "04_results",
    "step3B_GSE162610_state_scoring",
    "tables",
    "09_cluster_by_sample_after_cleanup.csv"
  ),

  step3C_thresholds = file.path(
    ROOT,
    "04_results",
    "step3C_GSE162610_activated_DAM_transition_mapping",
    "tables",
    "02_main_transition_thresholds.csv"
  ),
  step3C_time = file.path(
    ROOT,
    "04_results",
    "step3C_GSE162610_activated_DAM_transition_mapping",
    "tables",
    "04_transition_status_by_timepoint.csv"
  ),
  step3C_3dpi = file.path(
    ROOT,
    "04_results",
    "step3C_GSE162610_activated_DAM_transition_mapping",
    "tables",
    "10_3dpi_transition_by_biological_sample.csv"
  ),
  step3C_sensitivity = file.path(
    ROOT,
    "04_results",
    "step3C_GSE162610_activated_DAM_transition_mapping",
    "tables",
    "12_transition_threshold_sensitivity.csv"
  ),

  step3E_state_summary = file.path(
    ROOT,
    "04_results",
    "step3E_GSE162610_noncycling_transition_refinement",
    "tables",
    "01_step3E_state_summary.csv"
  ),
  step3E_state_audit = file.path(
    ROOT,
    "04_results",
    "step3E_GSE162610_noncycling_transition_refinement",
    "tables",
    "02_step3E_state_cell_audit.csv"
  ),
  step3E_concordance = file.path(
    ROOT,
    "04_results",
    "step3E_GSE162610_noncycling_transition_refinement",
    "tables",
    "14_cell_vs_pseudobulk_concordance.csv"
  ),
  step3E_pb_meta = file.path(
    ROOT,
    "04_results",
    "step3E_GSE162610_noncycling_transition_refinement",
    "pseudobulk",
    "00_pseudobulk_sample_state_metadata.csv"
  ),

  sample_manifest = file.path(
    ROOT,
    "00_admin",
    "step0_GSE162610_final",
    "02_GSE162610_sample_manifest_FINAL.csv"
  )
)

missing_inputs <- names(P)[!vapply(P, file.exists, logical(1))]
if (length(missing_inputs) > 0) {
  stop(
    "Missing required Figure S2 input(s):\n",
    paste0(
      " - ",
      missing_inputs,
      ": ",
      unlist(P[missing_inputs]),
      collapse = "\n"
    )
  )
}

# ------------------------------------------------------------------------------
# 1. Packages
# ------------------------------------------------------------------------------
required_pkgs <- c(
  "Seurat",
  "SeuratObject",
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
  library(Seurat)
  library(SeuratObject)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(patchwork)
  library(scales)
})

# ------------------------------------------------------------------------------
# 2. Frozen colors / display labels
# ------------------------------------------------------------------------------
base_family <- "Arial"

time_levels <- c(
  "Uninjured",
  "1 dpi",
  "3 dpi",
  "7 dpi"
)

time_cols <- c(
  "Uninjured" = "#4C78A8",
  "1 dpi"     = "#EDC948",
  "3 dpi"     = "#F28E2B",
  "7 dpi"     = "#E45756"
)

resident_levels <- c(
  "Canonical resident",
  "Activated resident",
  "IFN-responsive resident",
  "Cycling resident"
)

resident_cols <- c(
  "Canonical resident"      = "#4C78A8",
  "Activated resident"      = "#E45756",
  "IFN-responsive resident" = "#1B9E77",
  "Cycling resident"        = "#B279A2"
)

resident_short <- c(
  "Canonical resident"      = "Canonical",
  "Activated resident"      = "Activated",
  "IFN-responsive resident" = "IFN-responsive",
  "Cycling resident"        = "Cycling"
)

state3E_levels <- c(
  "Homeostatic_root",
  "Early_activated_1dpi",
  "Shared_transition_3dpi",
  "Established_DAM_7dpi"
)

state3E_cols <- c(
  "Homeostatic_root"       = time_cols[["Uninjured"]],
  "Early_activated_1dpi"   = time_cols[["1 dpi"]],
  "Shared_transition_3dpi" = time_cols[["3 dpi"]],
  "Established_DAM_7dpi"   = time_cols[["7 dpi"]]
)

state3E_labels <- c(
  "Homeostatic_root"       = "Homeostatic\nroot",
  "Early_activated_1dpi"   = "Early activated\n1 dpi",
  "Shared_transition_3dpi" = "Shared transition\n3 dpi",
  "Established_DAM_7dpi"   = "Established DAM\n7 dpi"
)

state3E_labels_D <- c(
  "Homeostatic_root"       = "Homeostatic",
  "Early_activated_1dpi"   = "Early activated",
  "Shared_transition_3dpi" = "Shared transition",
  "Established_DAM_7dpi"   = "Established DAM"
)

comparison_levels <- c(
  "C1_Early_vs_Homeostatic",
  "C2_Transition_vs_Early",
  "C3_DAM_vs_Transition"
)

comparison_labels <- c(
  "C1_Early_vs_Homeostatic" = "Early vs\nhomeostatic",
  "C2_Transition_vs_Early"  = "Transition vs\nearly",
  "C3_DAM_vs_Transition"    = "DAM vs\ntransition"
)

comparison_cols <- c(
  "C1_Early_vs_Homeostatic" = "#D95F02",
  "C2_Transition_vs_Early"  = "#7570B3",
  "C3_DAM_vs_Transition"    = "#1B9E77"
)

# ------------------------------------------------------------------------------
# 3. Publication theme
# ------------------------------------------------------------------------------
theme_pub <- theme_classic(
  base_family = base_family,
  base_size = 7.2
) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    axis.title = element_text(
      family = base_family,
      size = 7.0,
      color = "black"
    ),
    axis.text = element_text(
      family = base_family,
      size = 6.0,
      color = "black"
    ),
    axis.line = element_line(
      linewidth = 0.35,
      color = "black"
    ),
    axis.ticks = element_line(
      linewidth = 0.30,
      color = "black"
    ),
    legend.title = element_text(
      family = base_family,
      size = 6.0,
      color = "black"
    ),
    legend.text = element_text(
      family = base_family,
      size = 5.5,
      color = "black"
    ),
    plot.title = element_text(
      family = base_family,
      size = 7.1,
      face = "bold",
      hjust = 0
    ),
    plot.subtitle = element_text(
      family = base_family,
      size = 5.7,
      color = "#4D4D4D",
      hjust = 0
    ),
    plot.tag = element_text(
      family = base_family,
      face = "bold",
      size = 11.5,
      color = "black"
    ),
    plot.tag.position = c(-0.035, 1.035),
    plot.margin = margin(
      3.2, 3.0, 3.2, 3.8,
      unit = "mm"
    )
  )

# ------------------------------------------------------------------------------
# 4. Helpers
# ------------------------------------------------------------------------------
collapse_resident_state <- function(step2_cluster) {
  x <- as.character(step2_cluster)

  dplyr::case_when(
    x == "0"           ~ "Canonical resident",
    x %in% c("3", "4") ~ "Activated resident",
    x == "9"           ~ "IFN-responsive resident",
    x %in% c("5", "7") ~ "Cycling resident",
    TRUE               ~ NA_character_
  )
}

sample_display <- function(x) {
  x <- as.character(x)

  dplyr::case_when(
    x == "uninj_sample1" ~ "U1",
    x == "uninj_sample2" ~ "U2",
    x == "uninj_sample3" ~ "U3",
    x == "1dpi_sample1"  ~ "1d1",
    x == "1dpi_sample2"  ~ "1d2",
    x == "1dpi_sample3"  ~ "1d3",
    x == "3dpi_sample1"  ~ "3d1",
    x == "3dpi_sample2"  ~ "3d2",
    x == "7dpi_sample1"  ~ "7d1",
    x == "7dpi_sample2"  ~ "7d2",
    TRUE ~ x
  )
}

sample_order <- c(
  "uninj_sample1",
  "uninj_sample2",
  "uninj_sample3",
  "1dpi_sample1",
  "1dpi_sample2",
  "1dpi_sample3",
  "3dpi_sample1",
  "3dpi_sample2",
  "7dpi_sample1",
  "7dpi_sample2"
)

row_zscore <- function(x) {
  x <- as.numeric(x)

  if (all(is.na(x))) {
    return(rep(NA_real_, length(x)))
  }

  s <- stats::sd(x, na.rm = TRUE)

  if (!is.finite(s) || s == 0) {
    return(rep(0, length(x)))
  }

  (x - mean(x, na.rm = TRUE)) / s
}

save_plot_quad <- function(
  p,
  stem,
  width_mm,
  height_mm,
  dpi = 600
) {
  ggsave(
    filename = file.path(OUT, paste0(stem, ".pdf")),
    plot = p,
    device = cairo_pdf,
    width = width_mm,
    height = height_mm,
    units = "mm",
    bg = "white"
  )

  ggsave(
    filename = file.path(OUT, paste0(stem, ".svg")),
    plot = p,
    device = svglite::svglite,
    width = width_mm,
    height = height_mm,
    units = "mm",
    bg = "white"
  )

  ggsave(
    filename = file.path(OUT, paste0(stem, ".tiff")),
    plot = p,
    device = ragg::agg_tiff,
    width = width_mm,
    height = height_mm,
    units = "mm",
    dpi = dpi,
    compression = "lzw",
    background = "white"
  )

  ggsave(
    filename = file.path(OUT, paste0(stem, "_review.png")),
    plot = p,
    device = ragg::agg_png,
    width = width_mm,
    height = height_mm,
    units = "mm",
    dpi = 220,
    background = "white"
  )
}

# ------------------------------------------------------------------------------
# 5. Input manifest
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
    "01_FigureS2_input_manifest.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 6. PANEL A + PANEL B | Step3B frozen resident-microglia object
# ------------------------------------------------------------------------------
message("Reading frozen Step3B state-scored object...")

micro <- readRDS(P$step3B_rds)

if (!inherits(micro, "Seurat")) {
  stop("Step3B input is not a Seurat object.")
}

DefaultAssay(micro) <- "RNA"

required_meta_B <- c(
  "step2_final_cluster",
  "sample_id",
  "timepoint",
  "AUCell_Homeostatic",
  "AUCell_Activated_injury",
  "AUCell_DAM_like_candidate",
  "AUCell_IFN_response",
  "AUCell_Cycling_S",
  "AUCell_Cycling_G2M"
)

missing_meta_B <- setdiff(
  required_meta_B,
  colnames(micro[[]])
)

if (length(missing_meta_B) > 0) {
  stop(
    "Step3B object is missing required metadata: ",
    paste(missing_meta_B, collapse = ", ")
  )
}

micro$resident_state_display <- collapse_resident_state(
  micro$step2_final_cluster
)

if (any(is.na(micro$resident_state_display))) {
  bad <- sort(
    unique(
      as.character(
        micro$step2_final_cluster[
          is.na(micro$resident_state_display)
        ]
      )
    )
  )

  stop(
    "Unexpected Step2 final cluster(s): ",
    paste(bad, collapse = ", "),
    ". No resident-state label was guessed."
  )
}

micro$resident_state_display <- factor(
  micro$resident_state_display,
  levels = resident_levels
)

micro$timepoint <- factor(
  as.character(micro$timepoint),
  levels = time_levels
)

# ------------------------------------------------------------------------------
# 6A. State-definition marker DotPlot
# ------------------------------------------------------------------------------
marker_df <- tibble::tribble(
  ~state,                     ~gene,
  "Canonical resident",       "P2ry12",
  "Canonical resident",       "Tmem119",
  "Canonical resident",       "Sall1",
  "Canonical resident",       "Fcrls",

  "Activated resident",       "Apoe",
  "Activated resident",       "Cst7",
  "Activated resident",       "Cd68",
  "Activated resident",       "Trem2",

  "IFN-responsive resident",  "Ifit1",
  "IFN-responsive resident",  "Isg15",
  "IFN-responsive resident",  "Mx1",
  "IFN-responsive resident",  "Usp18",

  "Cycling resident",         "Mki67",
  "Cycling resident",         "Top2a",
  "Cycling resident",         "Birc5",
  "Cycling resident",         "Mcm5"
)

marker_df <- marker_df %>%
  filter(gene %in% rownames(micro))

marker_avail <- marker_df %>%
  count(state, name = "n_available")

if (
  nrow(marker_avail) != 4 ||
  any(marker_avail$n_available < 3)
) {
  stop(
    "Panel A marker availability is insufficient. ",
    "Need at least 3 available marker genes per resident state."
  )
}

marker_order <- marker_df$gene

pA_dot <- Seurat::DotPlot(
  object = micro,
  features = marker_order,
  group.by = "resident_state_display",
  assay = "RNA",
  dot.scale = 5.0,
  scale = TRUE
) +
  scale_color_gradient2(
    low = "#2166AC",
    mid = "#F7F7F7",
    high = "#B2182B",
    midpoint = 0,
    breaks = c(-1, 0, 1),
    labels = c("-1", "0", "1"),
    name = "Average\nexpression"
  ) +
  scale_size(
    range = c(0.2, 3.2),
    name = "Percent\nexpressed"
  ) +
  scale_y_discrete(
    labels = resident_short,
    drop = FALSE
  ) +
  labs(
    tag = "A",
    title = "State-marker evidence",
    x = NULL,
    y = NULL
  ) +
  theme_pub +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    panel.grid.major = element_line(
      linewidth = 0.20,
      color = "#E8E8E8"
    ),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(
      size = 5.3,
      angle = 52,
      hjust = 1,
      vjust = 1
    ),
    axis.text.y = element_text(
      size = 5.8
    ),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.margin = margin(
      t = 0.5,
      unit = "mm"
    ),
    legend.key.height = grid::unit(
      2.4,
      "mm"
    ),
    legend.key.width = grid::unit(
      3.8,
      "mm"
    ),
    plot.margin = margin(
      2.4, 2.0, 2.5, 3.5,
      unit = "mm"
    )
  ) +
  guides(
    size = guide_legend(
      order = 1,
      title.position = "top"
    ),
    color = guide_colorbar(
      order = 2,
      title.position = "top",
      barwidth = grid::unit(13, "mm"),
      barheight = grid::unit(2.0, "mm")
    )
  )

# ------------------------------------------------------------------------------
# 6B. Compact AUCell state-program heatmap
# ------------------------------------------------------------------------------
auc_meta <- micro[[]] %>%
  transmute(
    resident_state_display = factor(
      as.character(resident_state_display),
      levels = resident_levels
    ),
    Homeostatic = AUCell_Homeostatic,
    Activated_injury = AUCell_Activated_injury,
    DAM_like = AUCell_DAM_like_candidate,
    IFN_response = AUCell_IFN_response,
    Cycling_S = AUCell_Cycling_S,
    Cycling_G2M = AUCell_Cycling_G2M
  )

auc_state <- auc_meta %>%
  group_by(resident_state_display) %>%
  summarise(
    across(
      Homeostatic:Cycling_G2M,
      ~ median(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = Homeostatic:Cycling_G2M,
    names_to = "program",
    values_to = "median_auc"
  ) %>%
  group_by(program) %>%
  mutate(
    z = row_zscore(median_auc)
  ) %>%
  ungroup()

program_levels <- c(
  "Homeostatic",
  "Activated_injury",
  "DAM_like",
  "IFN_response",
  "Cycling_S",
  "Cycling_G2M"
)

program_labels <- c(
  "Homeostatic"     = "Homeostatic",
  "Activated_injury"= "Activated injury",
  "DAM_like"        = "DAM-like",
  "IFN_response"    = "IFN response",
  "Cycling_S"       = "Cycling S",
  "Cycling_G2M"     = "Cycling G2/M"
)

auc_state <- auc_state %>%
  mutate(
    program = factor(
      program,
      levels = rev(program_levels)
    ),
    resident_state_display = factor(
      resident_state_display,
      levels = resident_levels
    )
  )

pA_auc <- ggplot(
  auc_state,
  aes(
    x = resident_state_display,
    y = program,
    fill = z
  )
) +
  geom_tile(
    linewidth = 0.32,
    color = "white"
  ) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "#F7F7F7",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-2.2, 2.2),
    oob = scales::squish,
    name = "Row Z"
  ) +
  scale_x_discrete(
    labels = resident_short,
    drop = FALSE
  ) +
  scale_y_discrete(
    labels = program_labels,
    drop = FALSE
  ) +
  labs(
    title = "AUCell programs",
    x = NULL,
    y = NULL
  ) +
  coord_fixed(ratio = 0.74) +
  theme_pub +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(
      size = 5.2,
      angle = 35,
      hjust = 1
    ),
    axis.text.y = element_text(
      size = 5.5
    ),
    legend.position = "bottom",
    legend.title = element_text(
      size = 5.4
    ),
    legend.text = element_text(
      size = 5.0
    ),
    legend.key.height = grid::unit(
      2.0,
      "mm"
    ),
    legend.key.width = grid::unit(
      12,
      "mm"
    ),
    plot.margin = margin(
      2.4, 2.4, 2.5, 1.0,
      unit = "mm"
    )
  ) +
  guides(
    fill = guide_colorbar(
      title.position = "top",
      barwidth = grid::unit(15, "mm"),
      barheight = grid::unit(2.0, "mm")
    )
  )

panelA <- patchwork::wrap_plots(
  pA_dot,
  pA_auc,
  ncol = 2,
  widths = c(0.62, 0.38)
)

# ------------------------------------------------------------------------------
# 7. PANEL B | resident-state representation across biological samples
# ------------------------------------------------------------------------------
comp_B <- micro[[]] %>%
  transmute(
    sample_id = as.character(sample_id),
    timepoint = factor(
      as.character(timepoint),
      levels = time_levels
    ),
    resident_state_display = factor(
      as.character(resident_state_display),
      levels = resident_levels
    )
  ) %>%
  filter(
    !is.na(sample_id),
    !is.na(timepoint),
    !is.na(resident_state_display)
  ) %>%
  count(
    sample_id,
    timepoint,
    resident_state_display,
    name = "n_cells"
  ) %>%
  group_by(
    sample_id,
    timepoint
  ) %>%
  tidyr::complete(
    resident_state_display = factor(
      resident_levels,
      levels = resident_levels
    ),
    fill = list(
      n_cells = 0
    )
  ) %>%
  mutate(
    sample_total = sum(n_cells),
    proportion = n_cells / sample_total
  ) %>%
  ungroup() %>%
  mutate(
    sample_id = factor(
      sample_id,
      levels = sample_order
    ),
    sample_label = factor(
      sample_display(sample_id),
      levels = sample_display(sample_order)
    )
  )

# Dot-size matrix: directly displays all 10 biological samples.
pB <- ggplot(
  comp_B,
  aes(
    x = sample_label,
    y = resident_state_display,
    size = proportion,
    color = resident_state_display
  )
) +
  geom_point(
    alpha = 0.95
  ) +
  geom_vline(
    xintercept = c(3.5, 6.5, 8.5),
    linewidth = 0.32,
    color = "#A5A5A5",
    linetype = "22"
  ) +
  annotate(
    "text",
    x = c(2, 5, 7.5, 9.5),
    y = 4.52,
    label = time_levels,
    family = base_family,
    size = 1.90,
    fontface = "bold",
    color = unname(time_cols)
  ) +
  scale_color_manual(
    values = resident_cols,
    labels = resident_short,
    drop = FALSE
  ) +
  scale_size_area(
    max_size = 7.0,
    limits = c(0, 1),
    breaks = c(0.1, 0.3, 0.5, 0.7),
    labels = scales::label_percent(
      accuracy = 1
    ),
    name = "Within-sample\nproportion"
  ) +
  scale_y_discrete(
    labels = resident_short,
    drop = FALSE
  ) +
  scale_x_discrete(
    drop = FALSE
  ) +
  coord_cartesian(
    ylim = c(0.55, 4.72),
    clip = "off"
  ) +
  labs(
    tag = "B",
    title = "State representation across biological samples",
    x = "Biological sample",
    y = NULL
  ) +
  theme_pub +
  theme(
    panel.grid.major = element_line(
      linewidth = 0.20,
      color = "#EAEAEA"
    ),
    panel.grid.minor = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(
      size = 5.4,
      angle = 45,
      hjust = 1
    ),
    axis.text.y = element_text(
      size = 5.5
    ),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.margin = margin(
      t = 0.4,
      unit = "mm"
    ),
    plot.margin = margin(
      3.2, 4.0, 3.0, 4.0,
      unit = "mm"
    )
  ) +
  guides(
    color = "none",
    size = guide_legend(
      title.position = "top",
      nrow = 1
    )
  )

# No longer need the full Step3B object in memory.
rm(micro)
invisible(gc())

# ------------------------------------------------------------------------------
# 8. PANEL C | Step3C transition robustness
# ------------------------------------------------------------------------------
message("Reading frozen Step3C transition-mapped object...")

trans <- readRDS(P$step3C_rds)

if (!inherits(trans, "Seurat")) {
  stop("Step3C input is not a Seurat object.")
}

required_meta_C <- c(
  "sample_id",
  "timepoint",
  "transition_status",
  "shared_transition_3dpi_highconf"
)

missing_meta_C <- setdiff(
  required_meta_C,
  colnames(trans[[]])
)

if (length(missing_meta_C) > 0) {
  stop(
    "Step3C object is missing required metadata: ",
    paste(missing_meta_C, collapse = ", ")
  )
}

trans_md <- trans[[]] %>%
  transmute(
    sample_id = as.character(sample_id),
    timepoint = factor(
      as.character(timepoint),
      levels = time_levels
    ),
    transition_status = as.character(
      transition_status
    ),
    shared_transition_3dpi_highconf = as.logical(
      shared_transition_3dpi_highconf
    )
  )

# Sample-level broad Activated/Early-DAM transition fraction.
trans_sample <- trans_md %>%
  group_by(
    sample_id,
    timepoint
  ) %>%
  summarise(
    n_cells = n(),
    transition_fraction = mean(
      transition_status ==
        "Activated_EarlyDAM_transition",
      na.rm = TRUE
    ),
    .groups = "drop"
  )

trans_mean <- trans_sample %>%
  group_by(timepoint) %>%
  summarise(
    mean_fraction = mean(
      transition_fraction,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

pC_time <- ggplot(
  trans_sample,
  aes(
    x = timepoint,
    y = transition_fraction
  )
) +
  geom_line(
    data = trans_mean,
    aes(
      y = mean_fraction,
      group = 1
    ),
    linewidth = 0.72,
    color = "#3F3F3F"
  ) +
  geom_point(
    data = trans_mean,
    aes(
      y = mean_fraction
    ),
    shape = 21,
    size = 2.7,
    stroke = 0.42,
    fill = "white",
    color = "black"
  ) +
  geom_point(
    aes(
      fill = timepoint
    ),
    shape = 21,
    size = 2.0,
    stroke = 0.32,
    color = "black",
    position = position_jitter(
      width = 0.08,
      height = 0,
      seed = 20260818
    )
  ) +
  scale_fill_manual(
    values = time_cols,
    drop = FALSE
  ) +
  scale_y_continuous(
    limits = c(0, 0.50),
    breaks = seq(
      0,
      0.5,
      by = 0.1
    ),
    labels = scales::label_percent(
      accuracy = 1
    ),
    expand = expansion(
      mult = c(0, 0.02)
    )
  ) +
  labs(
    tag = "C",
    title = "Transition over time",
    x = NULL,
    y = "Transition fraction"
  ) +
  theme_pub +
  theme(
    legend.position = "none",
    plot.title = element_text(
      size = 6.2,
      face = "bold",
      hjust = 0
    ),
    plot.subtitle = element_blank(),
    axis.text.x = element_text(
      size = 5.7
    ),
    axis.text.y = element_text(
      size = 5.5
    ),
    plot.margin = margin(
      3.0, 2.0, 3.0, 3.8,
      unit = "mm"
    )
  )

# Exact main thresholds from Step3C.
thr <- readr::read_csv(
  P$step3C_thresholds,
  show_col_types = FALSE
)

get_thr <- function(key) {
  v <- thr$value[
    thr$threshold == key
  ]

  if (length(v) != 1) {
    stop(
      "Could not uniquely resolve Step3C threshold: ",
      key
    )
  }

  as.numeric(v)
}

act_thr <- get_thr(
  "Activated_high_q90_cluster3"
)

early_thr <- get_thr(
  "EarlyDAM_high_q90_cluster3"
)

mature_thr <- get_thr(
  "MatureDAM_established_q25_cluster1"
)

three <- trans_md %>%
  filter(
    timepoint == "3 dpi"
  ) %>%
  group_by(sample_id) %>%
  summarise(
    n_cells = n(),
    n_highconf = sum(
      shared_transition_3dpi_highconf,
      na.rm = TRUE
    ),
    fraction_highconf = mean(
      shared_transition_3dpi_highconf,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  mutate(
    sample_id = factor(
      sample_id,
      levels = c(
        "3dpi_sample1",
        "3dpi_sample2"
      )
    ),
    sample_label = c(
      "Sample 1",
      "Sample 2"
    )
  )

pC_rep <- ggplot(
  three,
  aes(
    x = sample_label,
    y = fraction_highconf
  )
) +
  geom_col(
    width = 0.62,
    fill = time_cols[["3 dpi"]],
    alpha = 0.94
  ) +
  geom_text(
    aes(
      label = paste0(
        scales::percent(
          fraction_highconf,
          accuracy = 0.1
        ),
        "\n(n=",
        scales::comma(n_highconf),
        ")"
      )
    ),
    family = base_family,
    size = 2.0,
    vjust = -0.25,
    lineheight = 0.95,
    color = "#303030"
  ) +
  scale_y_continuous(
    limits = c(
      0,
      max(
        0.46,
        max(three$fraction_highconf) * 1.18
      )
    ),
    labels = scales::label_percent(
      accuracy = 10
    ),
    expand = expansion(
      mult = c(0, 0)
    )
  ) +
  labs(
    title = "3 dpi replicate support",
    x = NULL,
    y = "Fraction of 3 dpi cells"
  ) +
  theme_pub +
  theme(
    axis.text.x = element_text(
      size = 5.2,
      angle = 0,
      hjust = 0.5
    ),
    axis.text.y = element_text(
      size = 5.3
    ),
    plot.title = element_text(
      size = 5.95,
      face = "bold",
      hjust = 0
    ),
    plot.subtitle = element_blank(),
    plot.margin = margin(
      3.0, 2.5, 3.0, 2.0,
      unit = "mm"
    )
  )

panelC <- patchwork::wrap_plots(
  pC_time,
  patchwork::wrap_elements(full = pC_rep, clip = FALSE),
  ncol = 2,
  widths = c(0.55, 0.45)
)

rm(trans)
invisible(gc())

# ------------------------------------------------------------------------------
# 9. PANEL D | non-cycling refinement + pseudobulk support
# ------------------------------------------------------------------------------
state3E <- readr::read_csv(
  P$step3E_state_summary,
  show_col_types = FALSE
) %>%
  mutate(
    step3E_state = factor(
      as.character(step3E_state),
      levels = state3E_levels
    )
  )

if (
  any(is.na(state3E$step3E_state)) ||
  nrow(state3E) != 4
) {
  stop(
    "Unexpected Step3E state summary. ",
    "Expected exactly four frozen non-cycling states."
  )
}

pD_state <- ggplot(
  state3E,
  aes(
    x = n_cells,
    y = step3E_state,
    fill = step3E_state
  )
) +
  geom_col(
    width = 0.62,
    alpha = 0.94
  ) +
  geom_text(
    aes(
      x = n_cells + 90,
      label = paste0(
        scales::comma(n_cells),
        " | ",
        n_samples,
        " samples"
      )
    ),
    family = base_family,
    size = 1.72,
    hjust = 0,
    color = "#333333"
  ) +
  scale_fill_manual(
    values = state3E_cols,
    drop = FALSE
  ) +
  scale_y_discrete(
    labels = state3E_labels_D,
    drop = FALSE
  ) +
  scale_x_continuous(
    limits = c(
      0,
      max(state3E$n_cells) * 1.53
    ),
    breaks = c(0, 2000, 4000),
    labels = scales::label_comma(),
    expand = expansion(
      mult = c(0, 0)
    )
  ) +
  labs(
    tag = "D",
    title = "Non-cycling state support",
    x = "Cells retained",
    y = NULL
  ) +
  theme_pub +
  theme(
    legend.position = "none",
    axis.text.y = element_text(
      size = 5.2
    ),
    axis.text.x = element_text(
      size = 5.0
    ),
    plot.title = element_text(
      size = 6.3,
      face = "bold",
      hjust = 0
    ),
    plot.subtitle = element_blank(),
    plot.margin = margin(
      3.0, 2.0, 3.0, 3.8,
      unit = "mm"
    )
  )

conc <- readr::read_csv(
  P$step3E_concordance,
  show_col_types = FALSE
)

required_conc <- c(
  "gene",
  "cell_log2FC",
  "pb_log2FC",
  "Sig",
  "comparison",
  "direction_concordant"
)

missing_conc <- setdiff(
  required_conc,
  names(conc)
)

if (length(missing_conc) > 0) {
  stop(
    "Step3E concordance table is missing: ",
    paste(missing_conc, collapse = ", ")
  )
}

conc_summary <- conc %>%
  filter(
    comparison %in% comparison_levels
  ) %>%
  group_by(comparison) %>%
  summarise(
    n_genes = n(),
    pearson_r = suppressWarnings(
      cor(
        cell_log2FC,
        pb_log2FC,
        use = "complete.obs",
        method = "pearson"
      )
    ),
    n_cell_sig = sum(
      Sig == "Sig",
      na.rm = TRUE
    ),
    direction_concordance_sig = mean(
      direction_concordant[
        Sig == "Sig"
      ],
      na.rm = TRUE
    ),
    .groups = "drop"
  )

if (
  nrow(conc_summary) != 3 ||
  any(!is.finite(conc_summary$pearson_r)) ||
  any(!is.finite(conc_summary$direction_concordance_sig))
) {
  stop(
    "Could not calculate complete Step3E cell-vs-pseudobulk concordance metrics."
  )
}

conc_long <- conc_summary %>%
  select(
    comparison,
    pearson_r,
    direction_concordance_sig
  ) %>%
  pivot_longer(
    cols = c(
      pearson_r,
      direction_concordance_sig
    ),
    names_to = "metric_raw",
    values_to = "value"
  ) %>%
  mutate(
    comparison = factor(
      comparison,
      levels = rev(comparison_levels)
    ),
    metric = factor(
      metric_raw,
      levels = c(
        "pearson_r",
        "direction_concordance_sig"
      ),
      labels = c(
        "Pearson r",
        "Directional concordance"
      )
    ),
    comparison_num = as.numeric(comparison),
    y_pos = comparison_num + ifelse(
      metric_raw == "pearson_r",
      -0.17,
      0.17
    ),
    value_label = ifelse(
      metric_raw == "pearson_r",
      sprintf("r=%.3f", value),
      sprintf("%.1f%%", 100 * value)
    ),
    label_x = ifelse(
      value >= 0.985,
      value - 0.005,
      value + 0.005
    ),
    label_hjust = ifelse(
      value >= 0.985,
      1,
      0
    )
  )

pD_conc <- ggplot(
  conc_long,
  aes(
    x = value,
    y = y_pos,
    shape = metric
  )
) +
  geom_vline(
    xintercept = 0.90,
    linewidth = 0.36,
    color = "#8A8A8A",
    linetype = "22"
  ) +
  geom_segment(
    aes(
      x = 0.80,
      xend = value,
      y = y_pos,
      yend = y_pos
    ),
    linewidth = 0.50,
    color = "#D4D4D4"
  ) +
  geom_point(
    aes(
      color = comparison
    ),
    size = 2.55,
    stroke = 0.5
  ) +
  geom_text(
    aes(
      x = label_x,
      y = y_pos,
      hjust = label_hjust,
      label = value_label
    ),
    family = base_family,
    size = 1.68,
    vjust = 0.5,
    color = "#303030",
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = comparison_cols,
    drop = FALSE
  ) +
  scale_shape_manual(
    values = c(16, 17)
  ) +
  scale_y_continuous(
    breaks = 1:3,
    labels = comparison_labels[rev(comparison_levels)],
    limits = c(0.55, 3.45),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_x_continuous(
    limits = c(0.80, 1.005),
    breaks = c(
      0.80,
      0.85,
      0.90,
      0.95,
      1.00
    ),
    labels = scales::label_percent(
      accuracy = 1
    ),
    expand = expansion(
      mult = c(0, 0)
    )
  ) +
  labs(
    title = "Pseudobulk concordance",
    x = "Correlation / concordance",
    y = NULL,
    shape = NULL
  ) +
  theme_pub +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.text = element_text(
      size = 4.65
    ),
    legend.key.width = grid::unit(3.0, "mm"),
    legend.spacing.x = grid::unit(0.8, "mm"),
    axis.text.y = element_text(
      size = 4.9
    ),
    axis.text.x = element_text(
      size = 5.0
    ),
    plot.title = element_text(
      size = 6.3,
      face = "bold",
      hjust = 0
    ),
    plot.subtitle = element_blank(),
    plot.margin = margin(
      2.2, 3.2, 2.4, 3.8,
      unit = "mm"
    )
  ) +
  guides(
    color = "none",
    shape = guide_legend(
      nrow = 1,
      byrow = TRUE,
      override.aes = list(size = 2.3)
    )
  )

panelD <- patchwork::wrap_plots(
  pD_state,
  pD_conc,
  ncol = 1,
  heights = c(0.41, 0.59)
)

# ------------------------------------------------------------------------------
# 10. Audits
# ------------------------------------------------------------------------------
audit <- data.frame(
  item = c(
    "Step3B_cells",
    "Canonical_resident_cells",
    "Activated_resident_cells",
    "IFN_responsive_resident_cells",
    "Cycling_resident_cells",
    "Step3C_cells",
    "3dpi_highconf_sample1_fraction",
    "3dpi_highconf_sample2_fraction",
    "Step3E_Homeostatic_root_cells",
    "Step3E_Early_activated_cells",
    "Step3E_Shared_transition_cells",
    "Step3E_Established_DAM_cells",
    "C1_Pearson_r",
    "C1_direction_concordance_sig",
    "C2_Pearson_r",
    "C2_direction_concordance_sig",
    "C3_Pearson_r",
    "C3_direction_concordance_sig"
  ),
  value = as.character(
    c(
      sum(comp_B$n_cells),
      sum(
        comp_B$n_cells[
          comp_B$resident_state_display ==
            "Canonical resident"
        ]
      ),
      sum(
        comp_B$n_cells[
          comp_B$resident_state_display ==
            "Activated resident"
        ]
      ),
      sum(
        comp_B$n_cells[
          comp_B$resident_state_display ==
            "IFN-responsive resident"
        ]
      ),
      sum(
        comp_B$n_cells[
          comp_B$resident_state_display ==
            "Cycling resident"
        ]
      ),
      sum(trans_sample$n_cells),
      three$fraction_highconf[
        three$sample_id ==
          "3dpi_sample1"
      ],
      three$fraction_highconf[
        three$sample_id ==
          "3dpi_sample2"
      ],
      state3E$n_cells[
        state3E$step3E_state ==
          "Homeostatic_root"
      ],
      state3E$n_cells[
        state3E$step3E_state ==
          "Early_activated_1dpi"
      ],
      state3E$n_cells[
        state3E$step3E_state ==
          "Shared_transition_3dpi"
      ],
      state3E$n_cells[
        state3E$step3E_state ==
          "Established_DAM_7dpi"
      ],
      conc_summary$pearson_r[
        conc_summary$comparison ==
          "C1_Early_vs_Homeostatic"
      ],
      conc_summary$direction_concordance_sig[
        conc_summary$comparison ==
          "C1_Early_vs_Homeostatic"
      ],
      conc_summary$pearson_r[
        conc_summary$comparison ==
          "C2_Transition_vs_Early"
      ],
      conc_summary$direction_concordance_sig[
        conc_summary$comparison ==
          "C2_Transition_vs_Early"
      ],
      conc_summary$pearson_r[
        conc_summary$comparison ==
          "C3_DAM_vs_Transition"
      ],
      conc_summary$direction_concordance_sig[
        conc_summary$comparison ==
          "C3_DAM_vs_Transition"
      ]
    )
  ),
  stringsAsFactors = FALSE
)

write.csv(
  audit,
  file.path(
    OUT,
    "00_FigureS2_data_audit.csv"
  ),
  row.names = FALSE
)

write.csv(
  comp_B,
  file.path(
    OUT,
    "02_PanelB_sample_state_proportions.csv"
  ),
  row.names = FALSE
)

write.csv(
  trans_sample,
  file.path(
    OUT,
    "03_PanelC_sample_transition_fractions.csv"
  ),
  row.names = FALSE
)

write.csv(
  conc_summary,
  file.path(
    OUT,
    "04_PanelD_pseudobulk_concordance_summary.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 11. Final 2 x 2 layout
#
# A and C are intentionally wider, as specified in the supplementary plan.
# ------------------------------------------------------------------------------
panelA_fixed <- patchwork::wrap_elements(
  full = panelA,
  clip = FALSE
)

panelB_fixed <- patchwork::wrap_elements(
  full = pB,
  clip = FALSE
)

panelC_fixed <- patchwork::wrap_elements(
  full = panelC,
  clip = FALSE
)

panelD_fixed <- patchwork::wrap_elements(
  full = panelD,
  clip = FALSE
)

row1 <- patchwork::wrap_plots(
  panelA_fixed,
  panelB_fixed,
  ncol = 2,
  widths = c(1.14, 0.86)
)

row2 <- patchwork::wrap_plots(
  panelC_fixed,
  panelD_fixed,
  ncol = 2,
  widths = c(1.06, 0.94)
)

figureS2 <- patchwork::wrap_plots(
  patchwork::wrap_elements(
    full = row1,
    clip = FALSE
  ),
  patchwork::wrap_elements(
    full = row2,
    clip = FALSE
  ),
  ncol = 1,
  heights = c(1, 1)
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
# 12. Export
# ------------------------------------------------------------------------------
message("Saving Figure S2 panel A...")
save_plot_quad(
  panelA,
  "FigureS2A_state_definition_evidence",
  width_mm = 112,
  height_mm = 82
)

message("Saving Figure S2 panel B...")
save_plot_quad(
  pB,
  "FigureS2B_state_representation_samples",
  width_mm = 78,
  height_mm = 82
)

message("Saving Figure S2 panel C...")
save_plot_quad(
  panelC,
  "FigureS2C_transition_robustness",
  width_mm = 112,
  height_mm = 82
)

message("Saving Figure S2 panel D...")
save_plot_quad(
  panelD,
  "FigureS2D_noncycling_pseudobulk_control",
  width_mm = 88,
  height_mm = 92
)

message("Saving complete Figure S2...")
save_plot_quad(
  figureS2,
  "FigureS2_GLIA_FINAL_v6_no_subtitles_no_script_dependency",
  width_mm = 190,
  height_mm = 186
)

# ------------------------------------------------------------------------------
# 13. Run summary
# ------------------------------------------------------------------------------
run_summary <- c(
  "Figure S2 | GLIA | v6 | no panel subtitles | frozen-input only",
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
  "Panel A:",
  "  Key fixed marker genes across Figure-1 resident-state display classes.",
  "  Compact AUCell heatmap uses six frozen programs and row Z scores across the four state classes.",
  "",
  "Panel B:",
  "  All 10 biological samples are shown directly.",
  "  Dot size is within-sample state proportion; no cell-level P values are used.",
  "",
  "Panel C:",
  "  Activated/Early-DAM transition fraction is calculated separately for each biological sample.",
  "  Open circles are timepoint means.",
  paste0(
    "  High-confidence 3 dpi gate thresholds: Activated >",
    sprintf("%.3f", act_thr),
    "; Early-DAM >",
    sprintf("%.3f", early_thr),
    "; Mature-DAM <",
    sprintf("%.3f", mature_thr),
    "."
  ),
  "",
  "Panel D:",
  "  Non-cycling refined trajectory states remain supported by >=2 biological samples.",
  "  Pseudobulk support is summarized by all-gene Pearson r and directional concordance among cell-level significant genes.",
  "",
  "No upstream analysis was rerun."
)

writeLines(
  run_summary,
  con = file.path(
    OUT,
    "FigureS2_GLIA_FINAL_v6_no_subtitles_no_script_dependency_run_summary.txt"
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
      "05_warnings_after_run.txt"
    )
  )
}

message("============================================================")
message("Figure S2 FINAL v4 FINAL completed.")
message("Output: ", OUT)
message("Main review file: FigureS2_GLIA_FINAL_v6_no_subtitles_no_script_dependency_review.png")
message("Main vector file: FigureS2_GLIA_FINAL_v6_no_subtitles_no_script_dependency.pdf / .svg")
message("============================================================")
