# ============================================================
# 118_Figure5E_GLIA_target_endothelial_response_spatial_coupling_v4_final.R
# Figure 5E | Spatial coupling of target program with endothelial responses
# ------------------------------------------------------------
# Main message
#   Spatial niches with higher target-program activity also show
#   higher endothelial response states predicted in Figure 4:
#     1) ENDOTHELIAL_ACTIVATION
#     2) TNF_NFKB_RESPONSE
#
# Layout
#                    Endothelial activation     TNF/NFkB response
#       7 dpi           scatter/density            scatter/density
#       2 months        scatter/density            scatter/density
#
# X-axis
#   Target-program spatial score
#
# Y-axis
#   Endothelial-weighted response-program score
#
# Important
#   - Uses the SAME representative sections as Figure 5A/B/D:
#       106_C_7days_2
#       106_C_2months_2
#   - Response scores are computed exactly as in Figure 5D:
#       log-normalized spatial counts ->
#       gene-wise z ->
#       program mean ->
#       multiplied by normalized Endothelial RCTD weight ->
#       z-standardized within section
#   - Target-program score comes from frozen Step12F:
#       target_program_adjusted_z_within_section
# ============================================================

options(stringsAsFactors = FALSE, warn = 1)

suppressPackageStartupMessages({
  library(BPCells)
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(patchwork)
})

# ------------------------------------------------------------
# 0. User settings
# ------------------------------------------------------------
ROOT_CANDIDATES <- unique(c(Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE)), "."))

SECTION_7D <- "106_C_7days_2"
SECTION_2M <- "106_C_2months_2"

OUTPUT_DIR_NAME <- "09_final_figures/GLIA/Figure5/single_panels/Figure5E"

BASE_FAMILY <- "Arial"
DPI_OUT <- 600
SCALE_FACTOR <- 10000

# Point-density palette matching the Figure 5A/B scatter style
DENSITY_COLORS <- c(
  "#24327A",
  "#5D2E9C",
  "#A63A99",
  "#F08A65",
  "#F6E945"
)

POINT_SIZE <- 1.15
POINT_ALPHA <- 0.92
LM_LINEWIDTH <- 0.78

# ------------------------------------------------------------
# 1. Helpers
# ------------------------------------------------------------
stop2 <- function(...) stop(..., call. = FALSE)

locate_root <- function() {
  for (root in ROOT_CANDIDATES) {
    if (!dir.exists(root)) next
    key <- file.path(
      root,
      "04_results",
      "step12B_fast_BPCells_spatial_section_QC",
      "BPCells",
      "GSE234774_spatial_2d_full_counts"
    )
    if (dir.exists(key)) {
      return(normalizePath(root, winslash = "/", mustWork = TRUE))
    }
  }

  stop2(
    "Could not locate the xjz project with Step12B BPCells spatial counts.\n",
    "Please edit ROOT_CANDIDATES."
  )
}

match_genes_case_insensitive <- function(requested, available) {
  requested <- unique(trimws(as.character(requested)))
  requested <- requested[!is.na(requested) & nzchar(requested)]

  mp <- setNames(available, tolower(available))
  hit <- unname(mp[tolower(requested)])

  unique(hit[!is.na(hit)])
}

extract_program_genes <- function(gs_tbl, program_name) {
  row <- gs_tbl %>%
    dplyr::filter(response_program == program_name)

  if (nrow(row) != 1) {
    stop2(
      "Could not uniquely identify response program in Step8C gene-set table: ",
      program_name
    )
  }

  genes <- unlist(
    strsplit(
      as.character(row$present_genes[1]),
      "\\s*\\|\\s*"
    )
  )

  genes <- unique(trimws(genes))
  genes[nchar(genes) > 0]
}

safe_z <- function(x) {
  x <- as.numeric(x)
  s <- sd(x, na.rm = TRUE)

  if (!is.finite(s) || s <= 0) {
    return(rep(0, length(x)))
  }

  (x - mean(x, na.rm = TRUE)) / s
}

program_score <- function(log_expr_gene_by_spot) {
  z <- t(scale(t(as.matrix(log_expr_gene_by_spot))))
  z[!is.finite(z)] <- 0
  colMeans(z, na.rm = TRUE)
}

make_density_colors <- function(x, y, palette = DENSITY_COLORS) {
  ok <- is.finite(x) & is.finite(y)
  cols <- rep("#7F7F7F", length(x))

  if (sum(ok) < 10) {
    cols[ok] <- palette[1]
    return(cols)
  }

  z <- grDevices::densCols(
    x[ok],
    y[ok],
    colramp = grDevices::colorRampPalette(palette)
  )

  cols[ok] <- z
  cols
}

format_p_value <- function(p) {
  if (!is.finite(p)) return("NA")
  if (p < 2.2e-16) return("< 2.2e-16")
  if (p < 0.001) return(sprintf("%.2e", p))
  sprintf("%.3f", p)
}

# ------------------------------------------------------------
# 2. Frozen input paths
# ------------------------------------------------------------
ROOT <- locate_root()
ROOT_RESULTS <- file.path(ROOT, "04_results")

STEP12B_DIR <- file.path(
  ROOT_RESULTS,
  "step12B_fast_BPCells_spatial_section_QC"
)

SPATIAL_BPCELLS_DIR <- file.path(
  STEP12B_DIR,
  "BPCells",
  "GSE234774_spatial_2d_full_counts"
)

STEP12D_WEIGHTS_RDS <- file.path(
  ROOT_RESULTS,
  "step12D_GSE234774_full_mode_RCTD_per_capture_area",
  "RDS",
  "GSE234774_step12D_all_spots_normalized_RCTD_weights.rds"
)

STEP12F_SPOT_RDS <- file.path(
  ROOT_RESULTS,
  "step12F_GSE234774_independent_lesion_zones_target_program_validation",
  "RDS",
  "GSE234774_step12F_spot_independent_lesion_zones_and_scores.rds"
)

STEP8C_GENESET_TSV <- file.path(
  ROOT_RESULTS,
  "step8C_GSE162610_recipient_program_integration",
  "tables",
  "02_recipient_program_gene_coverage.tsv"
)

required_files <- c(
  STEP12D_WEIGHTS_RDS,
  STEP12F_SPOT_RDS,
  STEP8C_GENESET_TSV
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop2(
    "Required frozen input(s) are missing:\n",
    paste0("  - ", missing_files, collapse = "\n")
  )
}

if (!dir.exists(SPATIAL_BPCELLS_DIR)) {
  stop2("Spatial BPCells directory is missing:\n  ", SPATIAL_BPCELLS_DIR)
}

# ------------------------------------------------------------
# 3. Load frozen inputs
# ------------------------------------------------------------
message("Loading frozen Figure 5E inputs ...")

spatial_counts <- BPCells::open_matrix_dir(SPATIAL_BPCELLS_DIR)
weights <- as.matrix(readRDS(STEP12D_WEIGHTS_RDS))
spot_tbl <- as.data.frame(readRDS(STEP12F_SPOT_RDS))
gs_tbl <- read.delim(
  STEP8C_GENESET_TSV,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# ------------------------------------------------------------
# 4. Endothelial RCTD column
# ------------------------------------------------------------
endo_candidates <- colnames(weights)[
  grepl("^Endothelial$|endothelial", colnames(weights), ignore.case = TRUE)
]

if (length(endo_candidates) == 0) {
  stop2(
    "No Endothelial column was found in normalized Step12D RCTD weights.\n",
    "Observed columns:\n  ",
    paste(colnames(weights), collapse = ", ")
  )
}

ENDO_COL <- if ("Endothelial" %in% endo_candidates) {
  "Endothelial"
} else {
  endo_candidates[1]
}

message("Using RCTD endothelial column: ", ENDO_COL)

# ------------------------------------------------------------
# 5. Program gene sets
# ------------------------------------------------------------
genes_activation_requested <- extract_program_genes(
  gs_tbl,
  "ENDOTHELIAL_ACTIVATION"
)

genes_tnf_requested <- extract_program_genes(
  gs_tbl,
  "TNF_NFKB_RESPONSE"
)

available_genes <- rownames(spatial_counts)

genes_activation <- match_genes_case_insensitive(
  genes_activation_requested,
  available_genes
)

genes_tnf <- match_genes_case_insensitive(
  genes_tnf_requested,
  available_genes
)

if (length(genes_activation) < 3) {
  stop2(
    "Too few ENDOTHELIAL_ACTIVATION genes were found in spatial counts: ",
    length(genes_activation)
  )
}

if (length(genes_tnf) < 3) {
  stop2(
    "Too few TNF_NFKB_RESPONSE genes were found in spatial counts: ",
    length(genes_tnf)
  )
}

# ------------------------------------------------------------
# 6. Score one selected section
# ------------------------------------------------------------
score_section <- function(section_id) {
  sec <- spot_tbl %>%
    dplyr::filter(.data$section_id == .env$section_id)

  if (nrow(sec) == 0) {
    stop2("No Step12F spots found for section: ", section_id)
  }

  observed_sections <- unique(as.character(sec$section_id))
  if (length(observed_sections) != 1L || observed_sections[1] != section_id) {
    stop2(
      "Section filtering failed inside score_section(). Requested: ",
      section_id,
      "; observed: ",
      paste(observed_sections, collapse = " | ")
    )
  }

  message(
    "Scoring section ",
    section_id,
    " | spots = ",
    nrow(sec)
  )

  if (!"target_program_adjusted_z_within_section" %in% names(sec)) {
    stop2("Step12F spot table is missing target_program_adjusted_z_within_section.")
  }

  barcodes <- as.character(sec$barcode)

  missing_counts <- setdiff(barcodes, colnames(spatial_counts))
  missing_weights <- setdiff(barcodes, rownames(weights))

  if (length(missing_counts) > 0) {
    stop2(
      "Spatial BPCells counts are missing ",
      length(missing_counts),
      " barcodes for section ",
      section_id
    )
  }

  if (length(missing_weights) > 0) {
    stop2(
      "RCTD weights are missing ",
      length(missing_weights),
      " barcodes for section ",
      section_id
    )
  }

  all_program_genes <- unique(c(genes_activation, genes_tnf))

  section_counts_all <- spatial_counts[, barcodes, drop = FALSE]
  library_size <- as.numeric(colSums(section_counts_all))
  names(library_size) <- barcodes

  if (any(!is.finite(library_size)) || any(library_size <= 0)) {
    stop2("Invalid spatial library size encountered in section: ", section_id)
  }

  prog_counts <- as.matrix(
    spatial_counts[all_program_genes, barcodes, drop = FALSE]
  )

  log_expr <- log1p(
    sweep(prog_counts, 2, library_size, "/") * SCALE_FACTOR
  )

  activation_score <- program_score(
    log_expr[genes_activation, , drop = FALSE]
  )

  tnf_score <- program_score(
    log_expr[genes_tnf, , drop = FALSE]
  )

  endo_weight <- as.numeric(weights[barcodes, ENDO_COL])
  names(endo_weight) <- barcodes

  sec$target_program_score <- as.numeric(
    sec$target_program_adjusted_z_within_section
  )

  sec$endothelial_weight <- endo_weight[sec$barcode]

  sec$endothelial_activation_weighted <- safe_z(
    activation_score[sec$barcode] * endo_weight[sec$barcode]
  )

  sec$endothelial_tnf_nfkb_weighted <- safe_z(
    tnf_score[sec$barcode] * endo_weight[sec$barcode]
  )

  sec$scored_section_id <- section_id

  sec
}

sec_7d <- score_section(SECTION_7D) %>%
  dplyr::mutate(display_timepoint = "7 dpi")

sec_2m <- score_section(SECTION_2M) %>%
  dplyr::mutate(display_timepoint = "2 months")

plot_data <- dplyr::bind_rows(sec_7d, sec_2m)

# Explicit sanity check: the two timepoints must come from different sections.
used_sections <- unique(plot_data[, c("display_timepoint", "section_id")])

if (nrow(used_sections) != 2 ||
    SECTION_7D == SECTION_2M ||
    identical(
      sort(unique(sec_7d$barcode)),
      sort(unique(sec_2m$barcode))
    )) {
  stop2(
    "Figure 5E sanity check failed: 7 dpi and 2 months do not appear to use distinct sections."
  )
}

message("Figure 5E section audit:")
print(used_sections)

# ------------------------------------------------------------
# 7. Shared display limits for fair cross-timepoint comparison
# ------------------------------------------------------------
common_x_limits <- range(
  c(
    sec_7d$target_program_score,
    sec_2m$target_program_score
  ),
  na.rm = TRUE
)

activation_y_limits <- range(
  c(
    sec_7d$endothelial_activation_weighted,
    sec_2m$endothelial_activation_weighted
  ),
  na.rm = TRUE
)

tnf_y_limits <- range(
  c(
    sec_7d$endothelial_tnf_nfkb_weighted,
    sec_2m$endothelial_tnf_nfkb_weighted
  ),
  na.rm = TRUE
)

# Add a small visual margin while preserving the same limits across panels.
expand_limits <- function(lim, frac = 0.035) {
  span <- diff(lim)
  if (!is.finite(span) || span <= 0) span <- 1
  c(lim[1] - frac * span, lim[2] + frac * span)
}

common_x_limits <- expand_limits(common_x_limits)
activation_y_limits <- expand_limits(activation_y_limits)
tnf_y_limits <- expand_limits(tnf_y_limits)

message(
  "Common x limits: ",
  paste(round(common_x_limits, 3), collapse = " to ")
)
message(
  "Activation y limits: ",
  paste(round(activation_y_limits, 3), collapse = " to ")
)
message(
  "TNF/NFkB y limits: ",
  paste(round(tnf_y_limits, 3), collapse = " to ")
)

# ------------------------------------------------------------
# 8. Scatter-panel builder
# ------------------------------------------------------------
make_scatter <- function(
  df,
  x_col,
  y_col,
  y_lab,
  x_limits,
  y_limits,
  show_y = TRUE,
  show_x = TRUE
) {
  dat <- df %>%
    dplyr::select(
      barcode,
      x = all_of(x_col),
      y = all_of(y_col)
    ) %>%
    dplyr::filter(is.finite(x), is.finite(y))

  if (nrow(dat) < 10) {
    stop2("Too few finite points for plotting.")
  }

  rho <- suppressWarnings(
    cor(
      dat$x,
      dat$y,
      method = "spearman",
      use = "complete.obs"
    )
  )

  dat$density_col <- make_density_colors(dat$x, dat$y)

  label_txt <- paste0(
    "\u03C1 = ",
    format(round(rho, 2), nsmall = 2)
  )

  p <- ggplot(dat, aes(x = x, y = y)) +
    geom_point(
      colour = dat$density_col,
      size = POINT_SIZE,
      alpha = POINT_ALPHA
    ) +
    geom_smooth(
      method = "lm",
      formula = y ~ x,
      se = FALSE,
      linewidth = LM_LINEWIDTH,
      colour = "black"
    ) +
    annotate(
      "text",
      x = x_limits[1] + 0.06 * diff(x_limits),
      y = y_limits[2] - 0.07 * diff(y_limits),
      label = label_txt,
      hjust = 0,
      vjust = 1,
      family = BASE_FAMILY,
      size = 6.0 / .pt,
      fontface = "italic"
    ) +
    scale_x_continuous(
      limits = x_limits,
      expand = c(0, 0),
      breaks = scales::pretty_breaks(n = 4)
    ) +
    scale_y_continuous(
      limits = y_limits,
      expand = c(0, 0),
      breaks = scales::pretty_breaks(n = 4)
    ) +
    labs(
      x = if (show_x) "Target-program spatial score" else NULL,
      y = if (show_y) y_lab else NULL
    ) +
    theme_classic(base_family = BASE_FAMILY) +
    theme(
      axis.title.x = element_text(size = 9.4, face = "plain"),
      axis.title.y = element_text(size = 9.4, face = "plain"),
      axis.text = element_text(size = 8.0),
      axis.line = element_line(linewidth = 0.56),
      plot.margin = margin(2, 3, 2, 2)
    )

  if (!show_y) {
    p <- p + theme(axis.title.y = element_blank())
  }

  if (!show_x) {
    p <- p + theme(axis.title.x = element_blank())
  }

  list(
    plot = p,
    rho = rho,
    n = nrow(dat)
  )
}

# ------------------------------------------------------------
# 9. Build four panels
# ------------------------------------------------------------
sc_7d_activation <- make_scatter(
  sec_7d,
  "target_program_score",
  "endothelial_activation_weighted",
  "Endothelial activation score",
  common_x_limits,
  activation_y_limits,
  TRUE,
  TRUE
)

sc_7d_tnf <- make_scatter(
  sec_7d,
  "target_program_score",
  "endothelial_tnf_nfkb_weighted",
  "TNF/NF\u03baB response score",
  common_x_limits,
  tnf_y_limits,
  TRUE,
  TRUE
)

sc_2m_activation <- make_scatter(
  sec_2m,
  "target_program_score",
  "endothelial_activation_weighted",
  "Endothelial activation score",
  common_x_limits,
  activation_y_limits,
  TRUE,
  TRUE
)

sc_2m_tnf <- make_scatter(
  sec_2m,
  "target_program_score",
  "endothelial_tnf_nfkb_weighted",
  "TNF/NF\u03baB response score",
  common_x_limits,
  tnf_y_limits,
  TRUE,
  TRUE
)

# ------------------------------------------------------------
# 9. Labels + layout
# ------------------------------------------------------------
make_header <- function(label) {
  ggplot() +
    annotate(
      "text",
      x = 0.5,
      y = 0.5,
      label = label,
      size = 9.8 / .pt,
      family = BASE_FAMILY,
      fontface = "plain"
    ) +
    xlim(0, 1) + ylim(0, 1) +
    theme_void()
}

make_blank_header <- function(panel_letter = NULL) {
  ggplot() +
    {
      if (!is.null(panel_letter)) {
        annotate(
          "text",
          x = 0.16,
          y = 0.72,
          label = panel_letter,
          size = 16 / .pt,
          family = BASE_FAMILY,
          fontface = "bold",
          hjust = 0
        )
      }
    } +
    xlim(0, 1) + ylim(0, 1) +
    theme_void()
}

make_row_label <- function(label) {
  ggplot() +
    annotate(
      "text",
      x = 0.80,
      y = 0.50,
      label = label,
      angle = 90,
      size = 9.3 / .pt,
      family = BASE_FAMILY,
      fontface = "bold"
    ) +
    xlim(0, 1) + ylim(0, 1) +
    theme_void()
}

header_row <- make_blank_header("E") +
  make_header("Endothelial activation") +
  make_header("TNF/NF\u03baB response") +
  patchwork::plot_layout(widths = c(0.08, 1, 1))

row_7d <- make_row_label("7 dpi") +
  sc_7d_activation$plot +
  sc_7d_tnf$plot +
  patchwork::plot_layout(widths = c(0.08, 1, 1))

row_2m <- make_row_label("2 months") +
  sc_2m_activation$plot +
  sc_2m_tnf$plot +
  patchwork::plot_layout(widths = c(0.08, 1, 1))

final_plot <- header_row / row_7d / row_2m +
  patchwork::plot_layout(heights = c(0.11, 1, 1))

# ------------------------------------------------------------
# 10. Outputs
# ------------------------------------------------------------
out_dir <- file.path(ROOT, OUTPUT_DIR_NAME)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_base <- file.path(
  out_dir,
  "Figure5E_target_endothelial_response_spatial_coupling_v4_final"
)

pdf_device <- if (capabilities("cairo")) cairo_pdf else grDevices::pdf

ggsave(
  paste0(out_base, ".pdf"),
  final_plot,
  width = 8.2,
  height = 4.35,
  units = "in",
  device = pdf_device,
  bg = "white"
)

ggsave(
  paste0(out_base, ".png"),
  final_plot,
  width = 8.2,
  height = 4.35,
  units = "in",
  dpi = DPI_OUT,
  bg = "white"
)

if (requireNamespace("svglite", quietly = TRUE)) {
  ggsave(
    paste0(out_base, ".svg"),
    final_plot,
    width = 8.2,
    height = 4.35,
    units = "in",
    bg = "white"
  )
}

# ------------------------------------------------------------
# 11. Plot-data and audit
# ------------------------------------------------------------
plot_data_out <- plot_data %>%
  dplyr::select(
    tidyselect::any_of(
      c(
        "barcode",
        "section_id",
        "capture_area_id",
        "biological_replicate_id",
        "timepoint",
        "display_timepoint",
        "coordinate_x_fullres",
        "coordinate_y_fullres",
        "target_program_score",
        "target_program_adjusted_z_within_section",
        "endothelial_weight",
        "endothelial_activation_weighted",
        "endothelial_tnf_nfkb_weighted",
        "distance_to_core_NN"
      )
    )
  )

readr::write_csv(
  plot_data_out,
  paste0(out_base, "_plot_data.csv")
)

summary_tbl <- tibble::tibble(
  panel = c(
    "7dpi_activation",
    "7dpi_TNF_NFKB",
    "2months_activation",
    "2months_TNF_NFKB"
  ),
  section_id = c(
    unique(sec_7d$section_id)[1],
    unique(sec_7d$section_id)[1],
    unique(sec_2m$section_id)[1],
    unique(sec_2m$section_id)[1]
  ),
  n_spots = c(
    sc_7d_activation$n,
    sc_7d_tnf$n,
    sc_2m_activation$n,
    sc_2m_tnf$n
  ),
  spearman_r = c(
    sc_7d_activation$rho,
    sc_7d_tnf$rho,
    sc_2m_activation$rho,
    sc_2m_tnf$rho
  )
)

readr::write_csv(
  summary_tbl,
  paste0(out_base, "_correlation_summary.csv")
)

gene_audit <- tibble::tibble(
  program = c(
    "ENDOTHELIAL_ACTIVATION",
    "TNF_NFKB_RESPONSE"
  ),
  requested_n = c(
    length(genes_activation_requested),
    length(genes_tnf_requested)
  ),
  spatial_present_n = c(
    length(genes_activation),
    length(genes_tnf)
  ),
  spatial_present_genes = c(
    paste(genes_activation, collapse = " | "),
    paste(genes_tnf, collapse = " | ")
  )
)

readr::write_csv(
  gene_audit,
  paste0(out_base, "_gene_audit.csv")
)

audit <- tibble::tibble(
  item = c(
    "section_7dpi",
    "section_2months",
    "spatial_counts",
    "RCTD_weights",
    "Step8C_genesets",
    "endothelial_RCTD_column",
    "scoring",
    "endothelial_weighting",
    "display_standardization",
    "common_x_limits",
    "activation_y_limits",
    "TNF_NFkB_y_limits",
    "interpretation_boundary"
  ),
  value = c(
    SECTION_7D,
    SECTION_2M,
    SPATIAL_BPCELLS_DIR,
    STEP12D_WEIGHTS_RDS,
    STEP8C_GENESET_TSV,
    ENDO_COL,
    "log-normalized spot expression -> gene-wise z -> program mean",
    "program score multiplied by frozen normalized Endothelial RCTD weight",
    "weighted response signal z-standardized within section",
    paste(signif(common_x_limits, 5), collapse = " | "),
    paste(signif(activation_y_limits, 5), collapse = " | "),
    paste(signif(tnf_y_limits, 5), collapse = " | "),
    "descriptive spot-level Spearman R only; no spot-level inferential P value; not proof of direct ligand-receptor binding"
  )
)

readr::write_csv(
  audit,
  paste0(out_base, "_audit.csv")
)

message("============================================================")
message("Figure 5E completed.")
message("7 dpi section : ", SECTION_7D)
message("2 months      : ", SECTION_2M)
message("Output PDF    : ", paste0(out_base, ".pdf"))
message("============================================================")
