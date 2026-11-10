# ==============================================================================
# Figure S1 | GLIA | FINAL v8 FINAL
# GSE162610 data quality, doublet control and Harmony integration
#
# Input:
#   FigureS1_GLIA_plotting_payload/
#     plotting_payload/
#     source_files/
#
# Output:
#   FigureS1_GLIA_FINAL_v8_FINAL/
#     FigureS1A_QC_filtering_retention.{pdf,svg,tiff,png}
#     FigureS1B_doublet_rate_by_sample.{pdf,svg,tiff,png}
#     FigureS1C_pre_post_Harmony_integration.{pdf,svg,tiff,png}
#     FigureS1_GLIA_FINAL_v8_FINAL.{pdf,svg,tiff,png}
#     00_FigureS1_data_audit.csv
#     01_FigureS1_input_manifest.csv
#     sessionInfo.txt
#
# Design principles:
#   * concise but complete: 3 panels only
#   * no re-normalization, re-Harmony, re-clustering or re-annotation
#   * exact exported cell metadata and exact UMAP coordinates are used
#   * Arial throughout
#   * timepoint colors are frozen to Main Figure 1
#   * compact UMAP axes; avoid excessive white margins
# ==============================================================================

options(stringsAsFactors = FALSE)
set.seed(20260818)

# ------------------------------------------------------------------------------
# 0. Project root — DIRECT ORIGINAL-SOURCE MODE
# ------------------------------------------------------------------------------
# This script does NOT depend on:
#   * FigureS1_GLIA_plotting_payload.zip
#   * an extracted FigureS1_GLIA_plotting_payload folder
#   * any temporary upload-package files
#
# It reads the original Step1A–Step1E results directly from xjz.
# ------------------------------------------------------------------------------
ROOT <- Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
if (!dir.exists(ROOT)) {
  stop("Project root does not exist: ", ROOT)
}

OUT <- file.path(ROOT, "09_final_figures", "GLIA", "FigureS1_FINAL_v8")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

message("Project root: ", ROOT)
message("Output dir:   ", OUT)

# ------------------------------------------------------------------------------
# 1. Packages
# ------------------------------------------------------------------------------
required_pkgs <- c(
  "ggplot2",
  "dplyr",
  "readr",
  "patchwork",
  "scales",
  "ggrastr",
  "svglite",
  "ragg"
)

missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_pkgs) > 0) {
  stop(
    "Missing required R packages: ",
    paste(missing_pkgs, collapse = ", "),
    "\nInstall them first, e.g.:\n",
    "install.packages(c(",
    paste(sprintf('"%s"', missing_pkgs), collapse = ", "),
    "))"
  )
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(patchwork)
  library(scales)
})

# ------------------------------------------------------------------------------
# 2. Frozen style
# ------------------------------------------------------------------------------
base_family <- "Arial"

time_cols <- c(
  "Uninjured" = "#4C78A8",
  "1 dpi"     = "#EDC948",
  "3 dpi"     = "#F28E2B",
  "7 dpi"     = "#E45756"
)

time_levels <- names(time_cols)

# Sample colors are deliberately nested within the four Figure-1 timepoint
# color families. This keeps sample identities distinguishable without
# introducing an unrelated rainbow palette.
sample_cols <- c(
  "uninj_sample1" = "#2F5D8C",
  "uninj_sample2" = "#4C78A8",
  "uninj_sample3" = "#7EA6C9",
  "1dpi_sample1"  = "#A87700",
  "1dpi_sample2"  = "#D6A900",
  "1dpi_sample3"  = "#EDC948",
  "3dpi_sample1"  = "#C96817",
  "3dpi_sample2"  = "#F28E2B",
  "7dpi_sample1"  = "#B33F3F",
  "7dpi_sample2"  = "#E45756"
)

sample_order <- names(sample_cols)

sample_labels <- c(
  "uninj_sample1" = "Uninjured 1",
  "uninj_sample2" = "Uninjured 2",
  "uninj_sample3" = "Uninjured 3",
  "1dpi_sample1"  = "1 dpi 1",
  "1dpi_sample2"  = "1 dpi 2",
  "1dpi_sample3"  = "1 dpi 3",
  "3dpi_sample1"  = "3 dpi 1",
  "3dpi_sample2"  = "3 dpi 2",
  "7dpi_sample1"  = "7 dpi 1",
  "7dpi_sample2"  = "7 dpi 2"
)

theme_pub <- theme_classic(
  base_family = base_family,
  base_size = 7.2
) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    axis.title = element_text(size = 7.0, color = "black"),
    axis.text = element_text(size = 6.1, color = "black"),
    axis.line = element_line(linewidth = 0.35, color = "black"),
    axis.ticks = element_line(linewidth = 0.30, color = "black"),
    legend.title = element_text(size = 6.2, color = "black"),
    legend.text = element_text(size = 5.45, color = "black"),
    legend.key.height = grid::unit(2.6, "mm"),
    legend.key.width = grid::unit(2.6, "mm"),
    plot.title = element_text(
      size = 7.4, face = "bold", hjust = 0.5,
      margin = margin(b = 1.6, unit = "mm")
    ),
    plot.subtitle = element_text(
      size = 6.2, color = "#4D4D4D", hjust = 0.5,
      margin = margin(b = 1.2, unit = "mm")
    ),
    plot.tag = element_text(
      family = base_family,
      face = "bold",
      size = 11.5,
      color = "black"
    ),
    plot.tag.position = c(-0.035, 1.03),
    plot.margin = margin(3.2, 3.2, 3.2, 3.8, unit = "mm")
  )

theme_umap <- theme_pub +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    plot.margin = margin(2.2, 2.2, 2.2, 2.2, unit = "mm")
  )

# ------------------------------------------------------------------------------
# 3. Exact original-source inputs
# ------------------------------------------------------------------------------
P <- list(
  step1A_rds = file.path(
    ROOT, "04_results", "step1A_GSE162610_QC_audit",
    "GSE162610_raw_with_QCmetrics.rds"
  ),
  step1B_rds = file.path(
    ROOT, "04_results", "step1B_GSE162610_QC_filtered",
    "GSE162610_QC_filtered.rds"
  ),
  step1C_rds = file.path(
    ROOT, "04_results", "step1C_GSE162610_scDblFinder",
    "GSE162610_QC_with_scDblFinder_labels.rds"
  ),
  step1D_rds = file.path(
    ROOT, "04_results", "step1D_GSE162610_LogNorm_PCA_preHarmony",
    "GSE162610_LogNorm_PCA_preHarmony.rds"
  ),
  step1E_rds = file.path(
    ROOT, "04_results", "step1E_GSE162610_Harmony_cluster",
    "GSE162610_Harmony_multires_clustered.rds"
  ),
  qc_thresholds = file.path(
    ROOT, "04_results", "step1B_GSE162610_QC_filtered", "tables",
    "01_QC_thresholds.csv"
  ),
  retention_sample = file.path(
    ROOT, "04_results", "step1B_GSE162610_QC_filtered", "tables",
    "02_QC_retention_by_sample.csv"
  ),
  retention_time = file.path(
    ROOT, "04_results", "step1B_GSE162610_QC_filtered", "tables",
    "03_QC_retention_by_timepoint.csv"
  ),
  doublet_sample = file.path(
    ROOT, "04_results", "step1C_GSE162610_scDblFinder", "tables",
    "02_scDblFinder_summary_by_sample.csv"
  ),
  doublet_time = file.path(
    ROOT, "04_results", "step1C_GSE162610_scDblFinder", "tables",
    "03_scDblFinder_summary_by_timepoint.csv"
  ),
  harmony_params = file.path(
    ROOT, "04_results", "step1E_GSE162610_Harmony_cluster", "tables",
    "01_step1E_parameters.csv"
  ),
  sample_manifest = file.path(
    ROOT, "00_admin", "step0_GSE162610_final",
    "02_GSE162610_sample_manifest_FINAL.csv"
  )
)

missing_inputs <- names(P)[!vapply(P, file.exists, logical(1))]
if (length(missing_inputs) > 0) {
  stop(
    "Missing required original Figure S1 inputs:\n",
    paste0(" - ", missing_inputs, ": ", unlist(P[missing_inputs]), collapse = "\n")
  )
}

input_manifest <- data.frame(
  key = names(P),
  path = unlist(P, use.names = FALSE),
  size_MB = round(file.info(unlist(P))$size / 1024^2, 4),
  md5 = unname(tools::md5sum(unlist(P))),
  stringsAsFactors = FALSE
)
write.csv(
  input_manifest,
  file.path(OUT, "01_FigureS1_input_manifest.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 4. Read original RDS objects and extract only what Figure S1 needs
# ------------------------------------------------------------------------------
if (!requireNamespace("SeuratObject", quietly = TRUE)) {
  stop(
    "Package 'SeuratObject' is required.\n",
    "Install it first, then rerun:\n",
    "install.packages('SeuratObject')"
  )
}

safe_metadata <- function(obj) {
  md <- tryCatch(obj[[]], error = function(e) NULL)
  if (is.null(md)) {
    md <- tryCatch(obj@meta.data, error = function(e) NULL)
  }
  if (is.null(md)) stop("Could not extract Seurat metadata.")
  md <- as.data.frame(md)
  md$cell_barcode <- rownames(md)
  md <- md[, c("cell_barcode", setdiff(colnames(md), "cell_barcode")), drop = FALSE]
  md
}

get_embedding <- function(obj, preferred, fallback_pattern = "umap") {
  reds <- tryCatch(
    SeuratObject::Reductions(obj),
    error = function(e) names(obj@reductions)
  )

  if (preferred %in% reds) {
    rn <- preferred
  } else {
    hit <- reds[grepl(fallback_pattern, reds, ignore.case = TRUE)]
    if (length(hit) == 0) {
      stop(
        "Could not find reduction '", preferred,
        "' or any reduction matching '", fallback_pattern, "'.\n",
        "Available reductions: ", paste(reds, collapse = ", ")
      )
    }
    rn <- hit[1]
    warning(
      "Preferred reduction '", preferred,
      "' not found; using '", rn, "' instead."
    )
  }

  emb <- as.data.frame(SeuratObject::Embeddings(obj, reduction = rn))
  emb$cell_barcode <- rownames(emb)
  emb <- emb[, c("cell_barcode", setdiff(colnames(emb), "cell_barcode")), drop = FALSE]

  list(name = rn, data = emb)
}

message("Reading original Step1A RDS...")
obj1A <- readRDS(P$step1A_rds)
before <- safe_metadata(obj1A)
rm(obj1A)
invisible(gc())

message("Reading original Step1B RDS...")
obj1B <- readRDS(P$step1B_rds)
after <- safe_metadata(obj1B)
rm(obj1B)
invisible(gc())

message("Reading original Step1C RDS...")
obj1C <- readRDS(P$step1C_rds)
dblmeta <- safe_metadata(obj1C)
rm(obj1C)
invisible(gc())

message("Reading original Step1D pre-Harmony RDS...")
obj1D <- readRDS(P$step1D_rds)
premeta <- safe_metadata(obj1D)
pre_emb <- get_embedding(
  obj1D,
  preferred = "umap.preharmony",
  fallback_pattern = "umap"
)
rm(obj1D)
invisible(gc())

message("Reading original Step1E post-Harmony RDS...")
obj1E <- readRDS(P$step1E_rds)
postmeta <- safe_metadata(obj1E)
post_emb <- get_embedding(
  obj1E,
  preferred = "umap.harmony",
  fallback_pattern = "umap"
)
rm(obj1E)
invisible(gc())

message("Pre-Harmony reduction:  ", pre_emb$name)
message("Post-Harmony reduction: ", post_emb$name)

pre_xy <- pre_emb$data
post_xy <- post_emb$data
rm(pre_emb, post_emb)
invisible(gc())

# Read frozen audit/statistics tables.
qc_thresholds <- readr::read_csv(P$qc_thresholds, show_col_types = FALSE)
ret_s <- readr::read_csv(P$retention_sample, show_col_types = FALSE)
ret_t <- readr::read_csv(P$retention_time, show_col_types = FALSE)
dbl_s <- readr::read_csv(P$doublet_sample, show_col_types = FALSE)
dbl_t <- readr::read_csv(P$doublet_time, show_col_types = FALSE)
harmony_params <- readr::read_csv(P$harmony_params, show_col_types = FALSE)
sample_manifest <- readr::read_csv(P$sample_manifest, show_col_types = FALSE)

# Validate core metadata columns.
required_qc <- c(
  "cell_barcode", "sample_id", "timepoint",
  "nFeature_RNA", "nCount_RNA", "percent.mt"
)

for (nm in c("before", "after")) {
  dat <- get(nm)
  miss <- setdiff(required_qc, colnames(dat))
  if (length(miss) > 0) {
    stop(nm, " metadata is missing columns: ", paste(miss, collapse = ", "))
  }
}

if (!all(c("cell_barcode", "sample_id", "timepoint") %in% names(premeta))) {
  stop("pre-Harmony metadata lacks barcode/sample/timepoint fields.")
}
if (!all(c("cell_barcode", "sample_id", "timepoint") %in% names(postmeta))) {
  stop("post-Harmony metadata lacks barcode/sample/timepoint fields.")
}

# Standardize UMAP coordinate column names without assuming the reduction's
# internal key. Exact rownames/barcodes are preserved.
standardize_xy <- function(x, label) {
  coord_cols <- setdiff(colnames(x), "cell_barcode")
  if (length(coord_cols) < 2) {
    stop(label, " reduction has fewer than two coordinate columns.")
  }
  x <- x[, c("cell_barcode", coord_cols[1:2]), drop = FALSE]
  colnames(x) <- c("cell_barcode", "UMAP_1", "UMAP_2")
  x
}

pre_xy <- standardize_xy(pre_xy, "Pre-Harmony")
post_xy <- standardize_xy(post_xy, "Post-Harmony")

# Standardize factors.
standardize_meta <- function(x) {
  x %>%
    mutate(
      timepoint = factor(as.character(timepoint), levels = time_levels),
      sample_id = factor(as.character(sample_id), levels = sample_order)
    )
}

before <- standardize_meta(before)
after <- standardize_meta(after)
dblmeta <- standardize_meta(dblmeta)
premeta <- standardize_meta(premeta)
postmeta <- standardize_meta(postmeta)

if (any(is.na(before$timepoint)) || any(is.na(after$timepoint))) {
  stop("Unexpected timepoint labels detected in QC metadata.")
}
if (any(is.na(premeta$sample_id)) || any(is.na(postmeta$sample_id))) {
  stop("Unexpected sample IDs detected in pre/post-Harmony metadata.")
}

# UMAP joins use exact cell barcodes.
pre_df <- premeta %>%
  select(cell_barcode, sample_id, timepoint) %>%
  inner_join(pre_xy, by = "cell_barcode")

post_df <- postmeta %>%
  select(cell_barcode, sample_id, timepoint) %>%
  inner_join(post_xy, by = "cell_barcode")

if (nrow(pre_df) != nrow(premeta)) {
  stop("Pre-Harmony metadata/UMAP barcode join lost cells.")
}
if (nrow(post_df) != nrow(postmeta)) {
  stop("Post-Harmony metadata/UMAP barcode join lost cells.")
}

same_cells <- setequal(pre_df$cell_barcode, post_df$cell_barcode)
if (!same_cells) {
  stop("Pre- and post-Harmony UMAPs do not contain the same cell set.")
}

# Retention and doublet summaries.
ret_s <- ret_s %>%
  mutate(
    timepoint = factor(as.character(timepoint), levels = time_levels),
    sample_id = factor(as.character(sample_id), levels = sample_order),
    sample_display = unname(sample_labels[as.character(sample_id)])
  )

dbl_s <- dbl_s %>%
  mutate(
    timepoint = factor(as.character(timepoint), levels = time_levels),
    sample_id = factor(as.character(sample_id), levels = sample_order),
    sample_display = unname(sample_labels[as.character(sample_id)])
  )

if (nrow(ret_s) != 10 || nrow(dbl_s) != 10) {
  stop("Expected exactly 10 biological samples in retention/doublet tables.")
}

# ------------------------------------------------------------------------------
# 5. Audit
# ------------------------------------------------------------------------------
overall_doublet_pct <- 100 * sum(dbl_s$n_doublet) / sum(dbl_s$n_total)

audit <- data.frame(
  item = c(
    "cells_before_QC",
    "cells_after_QC",
    "cells_removed_by_QC",
    "overall_QC_retention_pct",
    "cells_after_doublet_removal",
    "doublets_detected",
    "overall_doublet_pct",
    "n_biological_samples",
    "preHarmony_cells",
    "postHarmony_cells",
    "pre_post_same_cell_set",
    "min_sample_QC_retention_pct",
    "max_sample_QC_retention_pct",
    "Harmony_group",
    "Harmony_PCs",
    "Harmony_theta",
    "provisional_resolution"
  ),
  value = as.character(c(
    nrow(before),
    nrow(after),
    nrow(before) - nrow(after),
    round(100 * nrow(after) / nrow(before), 3),
    nrow(premeta),
    sum(dbl_s$n_doublet),
    round(overall_doublet_pct, 3),
    length(unique(as.character(premeta$sample_id))),
    nrow(pre_df),
    nrow(post_df),
    same_cells,
    min(ret_s$retention_pct),
    max(ret_s$retention_pct),
    harmony_params$value[harmony_params$parameter == "Harmony group"],
    harmony_params$value[harmony_params$parameter == "PCs"],
    harmony_params$value[harmony_params$parameter == "Harmony theta"],
    harmony_params$value[harmony_params$parameter == "Provisional resolution"]
  )),
  stringsAsFactors = FALSE
)

write.csv(
  audit,
  file.path(OUT, "00_FigureS1_data_audit.csv"),
  row.names = FALSE
)

message("Validated:")
message("  QC: ", nrow(before), " -> ", nrow(after), " cells")
message("  Doublets: ", sum(dbl_s$n_doublet),
        " (", sprintf("%.2f", overall_doublet_pct), "%)")
message("  Harmony UMAP: ", nrow(pre_df), " identical singlet cells pre/post")

# ------------------------------------------------------------------------------
# 6. Helper: compact vector UMAP axes
#
# Unlike the Main Figure 1 helper, this S1-specific helper uses only a small
# external pocket (~7-8% of data span), preventing the four UMAPs from
# acquiring large white borders.
# ------------------------------------------------------------------------------
add_umap_axes_compact <- function(
  p,
  df,
  arrow_fraction = 0.110
) {
  xr <- range(df$UMAP_1, na.rm = TRUE)
  yr <- range(df$UMAP_2, na.rm = TRUE)
  dx <- diff(xr)
  dy <- diff(yr)

  # Dedicated lower-left pocket: arrows stay close to the embedding, while
  # axis labels sit clearly outside the arrow shafts.
  x0 <- xr[1] - 0.020 * dx
  y0 <- yr[1] - 0.030 * dy
  x1 <- x0 + arrow_fraction * dx
  y1 <- y0 + arrow_fraction * dy

  x_limits <- c(xr[1] - 0.085 * dx, xr[2] + 0.008 * dx)
  y_limits <- c(yr[1] - 0.090 * dy, yr[2] + 0.010 * dy)

  p +
    scale_x_continuous(
      limits = x_limits,
      expand = expansion(mult = 0)
    ) +
    scale_y_continuous(
      limits = y_limits,
      expand = expansion(mult = 0)
    ) +
    coord_equal(clip = "off") +
    annotate(
      "segment",
      x = x0, xend = x1,
      y = y0, yend = y0,
      linewidth = 0.44,
      color = "black",
      lineend = "round",
      arrow = grid::arrow(
        type = "closed",
        length = grid::unit(1.08, "mm")
      )
    ) +
    annotate(
      "segment",
      x = x0, xend = x0,
      y = y0, yend = y1,
      linewidth = 0.44,
      color = "black",
      lineend = "round",
      arrow = grid::arrow(
        type = "closed",
        length = grid::unit(1.08, "mm")
      )
    ) +
    annotate(
      "text",
      x = x0 + 0.50 * arrow_fraction * dx,
      y = y0 - 0.045 * dy,
      label = "UMAP_1",
      family = base_family,
      size = 1.23,
      hjust = 0.5,
      vjust = 1
    ) +
    annotate(
      "text",
      x = x0 - 0.045 * dx,
      y = y0 + 0.50 * arrow_fraction * dy,
      label = "UMAP_2",
      family = base_family,
      size = 1.23,
      angle = 90,
      hjust = 0.5,
      vjust = 0.5
    )
}

# ------------------------------------------------------------------------------
# 7. Panel A: QC filtering + sample-level retention
# ------------------------------------------------------------------------------

# QC thresholds from the frozen Step1B table
thr_feature <- suppressWarnings(
  as.numeric(gsub("[^0-9.]+", "", qc_thresholds$rule[
    qc_thresholds$metric == "nFeature_RNA"
  ]))
)
thr_count <- suppressWarnings(
  as.numeric(gsub("[^0-9.]+", "", qc_thresholds$rule[
    qc_thresholds$metric == "nCount_RNA"
  ]))
)
thr_mt <- suppressWarnings(
  as.numeric(gsub("[^0-9.]+", "", qc_thresholds$rule[
    qc_thresholds$metric == "percent.mt"
  ]))
)

qc_both <- bind_rows(
  before %>% mutate(stage = "Before QC"),
  after  %>% mutate(stage = "After QC")
) %>%
  mutate(stage = factor(stage, levels = c("Before QC", "After QC")))

make_qc_violin <- function(metric, title, threshold, upper_q = 0.997) {
  vals <- qc_both[[metric]]
  ymax <- as.numeric(stats::quantile(vals, probs = upper_q, na.rm = TRUE))
  ymax <- max(ymax, threshold * 1.08, na.rm = TRUE)
  # Add a small top margin so extreme observations/violin tips never touch
  # the upper panel boundary. This is particularly important for 1 dpi UMI counts.
  ymax <- ymax * 1.035

  ggplot(qc_both, aes(
    x = timepoint,
    y = .data[[metric]],
    fill = timepoint
  )) +
    geom_violin(
      scale = "width",
      trim = TRUE,
      linewidth = 0.26,
      color = "#424242",
      alpha = 0.88
    ) +
    stat_summary(
      fun = median,
      geom = "point",
      size = 0.62,
      color = "black"
    ) +
    geom_hline(
      yintercept = threshold,
      linetype = "22",
      linewidth = 0.38,
      color = "#6B6B6B"
    ) +
    facet_wrap(~stage, ncol = 1) +
    scale_fill_manual(values = time_cols, drop = FALSE) +
    coord_cartesian(ylim = c(0, ymax), clip = "off") +
    labs(
      title = title,
      x = NULL,
      y = NULL
    ) +
    theme_pub +
    theme(
      legend.position = "none",
      strip.background = element_blank(),
      strip.text = element_text(
        family = base_family,
        size = 5.9,
        color = "#404040",
        face = "plain",
        margin = margin(0.5, 0, 0.8, 0, unit = "mm")
      ),
      axis.text.x = element_text(
        size = 5.8,
        angle = 0,
        hjust = 0.5,
        margin = margin(t = 0.6, unit = "mm")
      ),
      axis.text.y = element_text(size = 5.5),
      axis.ticks.x = element_blank(),
      plot.margin = margin(2.0, 2.0, 1.5, 2.4, unit = "mm")
    )
}

pA_feature <- make_qc_violin(
  "nFeature_RNA",
  "Detected genes",
  thr_feature
) +
  labs(tag = "A")

pA_count <- make_qc_violin(
  "nCount_RNA",
  "UMI counts",
  thr_count,
  upper_q = 1.000
) +
  scale_y_continuous(
    labels = scales::label_number(big.mark = ","),
    expand = expansion(mult = c(0, 0.02))
  )

pA_mt <- make_qc_violin(
  "percent.mt",
  "Mitochondrial reads (%)",
  thr_mt,
  upper_q = 0.999
)

# Retention: keep this directly beneath the QC distributions with minimal gap.
ret_plot <- ret_s %>%
  mutate(
    sample_display = factor(
      sample_display,
      levels = rev(unname(sample_labels[sample_order]))
    )
  )

pA_retention <- ggplot(
  ret_plot,
  aes(y = sample_display, x = retention_pct)
) +
  geom_segment(
    aes(x = 97.4, xend = retention_pct, yend = sample_display),
    linewidth = 0.45,
    color = "#D5D5D5",
    lineend = "round"
  ) +
  geom_point(
    aes(color = timepoint),
    size = 1.65,
    alpha = 0.98
  ) +
  geom_text(
    aes(
      x = pmin(retention_pct + 0.08, 100.14),
      label = sprintf("%.2f%%", retention_pct)
    ),
    family = base_family,
    size = 1.75,
    hjust = 0,
    color = "#333333"
  ) +
  scale_color_manual(values = time_cols, drop = FALSE) +
  scale_x_continuous(
    limits = c(97.35, 100.42),
    breaks = c(97.5, 98.5, 99.5, 100),
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = 0)
  ) +
  labs(
    title = "Sample-level retention",
    x = "Cells retained after QC",
    y = NULL
  ) +
  theme_pub +
  theme(
    legend.position = "none",
    plot.title = element_text(
      size = 6.8,
      face = "bold",
      hjust = 0,
      margin = margin(b = 0.8, unit = "mm")
    ),
    axis.text.y = element_text(size = 5.5),
    axis.text.x = element_text(size = 5.4),
    axis.title.x = element_text(size = 5.8),
    axis.ticks.y = element_blank(),
    plot.margin = margin(0.8, 3.0, 2.0, 3.8, unit = "mm")
  )

panelA_top <- pA_feature | pA_count | pA_mt
panelA <- (panelA_top / pA_retention) +
  plot_layout(heights = c(1.00, 0.58))

# ------------------------------------------------------------------------------
# 8. Panel B: scDblFinder rate by biological sample
# ------------------------------------------------------------------------------
dbl_plot <- dbl_s %>%
  mutate(
    sample_display = factor(
      sample_display,
      levels = rev(unname(sample_labels[sample_order]))
    )
  )

xmax_dbl <- ceiling(max(dbl_plot$detected_doublet_pct) + 1.0)

pB <- ggplot(
  dbl_plot,
  aes(
    y = sample_display,
    x = detected_doublet_pct,
    fill = timepoint
  )
) +
  geom_col(
    width = 0.68,
    linewidth = 0,
    alpha = 0.93
  ) +
  geom_vline(
    xintercept = overall_doublet_pct,
    linetype = "22",
    linewidth = 0.40,
    color = "#606060"
  ) +
  geom_text(
    aes(
      x = detected_doublet_pct + 0.08,
      label = sprintf("%.2f%%", detected_doublet_pct)
    ),
    family = base_family,
    size = 1.82,
    hjust = 0,
    color = "#222222"
  ) +
  annotate(
    "text",
    x = overall_doublet_pct + 0.08,
    y = 10.45,
    label = paste0("Overall ", sprintf("%.2f%%", overall_doublet_pct)),
    family = base_family,
    size = 1.78,
    hjust = 0,
    vjust = 0,
    color = "#555555"
  ) +
  scale_fill_manual(values = time_cols, drop = FALSE) +
  scale_x_continuous(
    limits = c(0, xmax_dbl),
    breaks = seq(0, xmax_dbl, by = 2),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    tag = "B",
    x = "Detected doublets (%)",
    y = NULL,
    fill = "Time after SCI"
  ) +
  theme_pub +
  theme(
    axis.text.y = element_text(size = 5.45),
    axis.text.x = element_text(size = 5.7),
    axis.title.x = element_text(size = 6.3),
    axis.ticks.y = element_blank(),
    legend.position = "none",
    plot.margin = margin(3.0, 3.0, 3.0, 4.0, unit = "mm")
  )

# ------------------------------------------------------------------------------
# 9. Panel C: pre- vs post-Harmony integration
# ------------------------------------------------------------------------------
make_umap <- function(
  df,
  color_col,
  color_values,
  color_labels = waiver(),
  title,
  subtitle,
  point_size = 0.21,
  alpha = 0.80,
  tag = NULL,
  guide_nrow = 1
) {
  p <- ggplot(
    df,
    aes(
      x = UMAP_1,
      y = UMAP_2,
      color = .data[[color_col]]
    )
  ) +
    ggrastr::geom_point_rast(
      size = point_size,
      alpha = alpha,
      raster.dpi = 500,
      stroke = 0
    ) +
    scale_color_manual(
      values = color_values,
      labels = color_labels,
      drop = FALSE
    ) +
    labs(
      title = title,
      subtitle = subtitle,
      tag = tag,
      color = NULL
    ) +
    theme_umap +
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.margin = margin(t = 0.5, unit = "mm"),
      legend.key.height = grid::unit(2.0, "mm"),
      legend.key.width = grid::unit(2.0, "mm"),
      legend.text = element_text(size = 4.55),
      plot.title = element_text(size = 6.8, face = "bold", hjust = 0.5, margin = margin(b = 0.4, unit = "mm")),
      plot.subtitle = element_blank(),
      plot.margin = margin(0.8, 0.9, 0.8, 0.9, unit = "mm")
    ) +
    guides(
      color = guide_legend(
        nrow = guide_nrow,
        byrow = TRUE,
        override.aes = list(size = 1.95, alpha = 1)
      )
    )

  add_umap_axes_compact(p, df)
}

pC_pre_sample <- make_umap(
  pre_df,
  color_col = "sample_id",
  color_values = sample_cols,
  color_labels = sample_labels,
  title = "Pre-Harmony",
  subtitle = NULL,
  point_size = 0.20,
  alpha = 0.72,
  tag = "C",
  guide_nrow = 2
)

pC_post_sample <- make_umap(
  post_df,
  color_col = "sample_id",
  color_values = sample_cols,
  color_labels = sample_labels,
  title = "Post-Harmony",
  subtitle = NULL,
  point_size = 0.20,
  alpha = 0.72,
  guide_nrow = 2
)

pC_pre_time <- make_umap(
  pre_df,
  color_col = "timepoint",
  color_values = time_cols,
  title = "Pre-Harmony",
  subtitle = NULL,
  point_size = 0.20,
  alpha = 0.78,
  guide_nrow = 1
)

pC_post_time <- make_umap(
  post_df,
  color_col = "timepoint",
  color_values = time_cols,
  title = "Post-Harmony",
  subtitle = NULL,
  point_size = 0.20,
  alpha = 0.78,
  guide_nrow = 1
)

# Each row collects its own guide. A single row-level descriptor replaces the
# duplicated subtitle under both UMAPs, giving the actual UMAP bodies more room.
panelC_sample <- ((pC_pre_sample | pC_post_sample) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Biological sample",
    theme = theme(
      plot.title = element_text(
        family = base_family,
        size = 6.0,
        face = "plain",
        color = "#4D4D4D",
        hjust = 0.5,
        margin = margin(b = 0.2, unit = "mm")
      )
    )
  )) &
  theme(legend.position = "bottom")

panelC_time <- ((pC_pre_time | pC_post_time) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Time after SCI",
    theme = theme(
      plot.title = element_text(
        family = base_family,
        size = 6.0,
        face = "plain",
        color = "#4D4D4D",
        hjust = 0.5,
        margin = margin(b = 0.2, unit = "mm")
      )
    )
  )) &
  theme(legend.position = "bottom")

panelC <- (panelC_sample / panelC_time) +
  plot_layout(heights = c(1, 1))

# ------------------------------------------------------------------------------
# 10. Final layout — explicit nesting fix
# ------------------------------------------------------------------------------
# IMPORTANT:
# patchwork may flatten nested patchworks when a new layout specification is
# applied at the parent level. In the previous version, row2 could therefore be
# interpreted as three children (pB + the two rows inside panelC), while
# widths = c(0.27, 0.73) only defined two columns. This caused:
#   "Need 3 panels, but together nrow and ncol only provide 2."
#
# Freeze each completed sub-composition as ONE patchwork-compliant element
# before assembling the final page.
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

row2 <- patchwork::wrap_plots(
  panelB_fixed,
  panelC_fixed,
  ncol = 2,
  widths = c(0.27, 0.73)
)

figureS1 <- patchwork::wrap_plots(
  panelA_fixed,
  patchwork::wrap_elements(full = row2, clip = FALSE),
  ncol = 1,
  heights = c(0.98, 1.04)
) +
  plot_annotation(
    theme = theme(
      plot.background = element_rect(fill = "white", color = NA)
    )
  )

# ------------------------------------------------------------------------------
# 11. Export helpers
# ------------------------------------------------------------------------------
save_plot_quad <- function(
  p,
  stem,
  width_mm,
  height_mm,
  dpi = 600
) {
  pdf_file <- file.path(OUT, paste0(stem, ".pdf"))
  svg_file <- file.path(OUT, paste0(stem, ".svg"))
  tif_file <- file.path(OUT, paste0(stem, ".tiff"))
  png_file <- file.path(OUT, paste0(stem, "_review.png"))

  ggsave(
    filename = pdf_file,
    plot = p,
    device = cairo_pdf,
    width = width_mm,
    height = height_mm,
    units = "mm",
    bg = "white"
  )

  ggsave(
    filename = svg_file,
    plot = p,
    device = svglite::svglite,
    width = width_mm,
    height = height_mm,
    units = "mm",
    bg = "white"
  )

  ggsave(
    filename = tif_file,
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
    filename = png_file,
    plot = p,
    device = ragg::agg_png,
    width = width_mm,
    height = height_mm,
    units = "mm",
    dpi = 220,
    background = "white"
  )
}

message("Saving Panel A...")
save_plot_quad(
  panelA,
  "FigureS1A_QC_filtering_retention",
  width_mm = 190,
  height_mm = 103
)

message("Saving Panel B...")
save_plot_quad(
  pB,
  "FigureS1B_doublet_rate_by_sample",
  width_mm = 66,
  height_mm = 104
)

message("Saving Panel C...")
save_plot_quad(
  panelC,
  "FigureS1C_pre_post_Harmony_integration",
  width_mm = 124,
  height_mm = 104
)

message("Saving complete Figure S1...")
save_plot_quad(
  figureS1,
  "FigureS1_GLIA_FINAL_v8_FINAL",
  width_mm = 190,
  height_mm = 204
)

# ------------------------------------------------------------------------------
# 12. Run summary
# ------------------------------------------------------------------------------
run_summary <- c(
  "Figure S1 GLIA FINAL v8 FINAL",
  "========================",
  paste0("Project root: ", ROOT),
  paste0("Output: ", OUT),
  "",
  paste0("QC cells: ", format(nrow(before), big.mark = ","),
         " -> ", format(nrow(after), big.mark = ",")),
  paste0(
    "Overall QC retention: ",
    sprintf("%.2f%%", 100 * nrow(after) / nrow(before))
  ),
  paste0(
    "Sample retention range: ",
    sprintf("%.2f%%", min(ret_s$retention_pct)),
    " to ",
    sprintf("%.2f%%", max(ret_s$retention_pct))
  ),
  paste0(
    "Detected doublets: ",
    format(sum(dbl_s$n_doublet), big.mark = ","),
    " / ",
    format(sum(dbl_s$n_total), big.mark = ","),
    " (",
    sprintf("%.2f%%", overall_doublet_pct),
    ")"
  ),
  paste0(
    "Pre/post-Harmony cells: ",
    format(nrow(pre_df), big.mark = ","),
    "; identical cell set = ",
    same_cells
  ),
  "",
  "Panel A: nFeature_RNA, nCount_RNA and percent.mt distributions before/after QC; sample-level retention directly beneath.",
  "Panel B: sample-level scDblFinder rate; redundant timepoint legend removed because sample labels already encode timepoint.",
  "Panel C: 2x2 pre/post-Harmony UMAP; sample and timepoint views; longer vector UMAP axes.",
  "",
  "Original Step1A-Step1E RDS objects were read directly; only metadata and existing UMAP embeddings were extracted. No upstream analysis was rerun."
)

writeLines(
  run_summary,
  con = file.path(OUT, "FigureS1_GLIA_FINAL_v8_FINAL_run_summary.txt"),
  useBytes = TRUE
)

capture.output(
  sessionInfo(),
  file = file.path(OUT, "sessionInfo.txt")
)

# Record non-fatal warnings for review without interrupting export.
warn_obj <- warnings()
if (!is.null(warn_obj)) {
  capture.output(
    warn_obj,
    file = file.path(OUT, "02_warnings_after_run.txt")
  )
}

message("============================================================")
message("Figure S1 FINAL v8 FINAL completed.")
message("Output: ", OUT)
message("Main review file: FigureS1_GLIA_FINAL_v8_FINAL_review.png")
message("Main vector file: FigureS1_GLIA_FINAL_v8_FINAL.pdf / .svg")
message("============================================================")
