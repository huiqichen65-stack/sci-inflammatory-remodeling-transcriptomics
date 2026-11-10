# ============================================================
# 117_Figure5D_GLIA_endothelial_response_spatial_maps_v5_rowlabel_fix.R
# v5 crop-cleanup only: preserve v4 map/legend typography and widen the
# left row-label gutter so '7 dpi' and '2 months' render fully.
# Figure 5D | Spatial endothelial response-state validation
# ------------------------------------------------------------
# Main message
#   The endothelial-rich lesion niche identified in Figure 5A-C
#   also shows the recipient response programs predicted by Figure 4:
#     1) ENDOTHELIAL_ACTIVATION
#     2) TNF_NFKB_RESPONSE
#
# Layout
#                 Endothelial activation     TNF/NF-kB response
#       7 dpi             map                      map
#       2 months          map                      map
#
# Important
#   - Uses the SAME representative sections as Figure 5A/B:
#       106_C_7days_2
#       106_C_2months_2
#   - Uses frozen Step8C response-program gene sets.
#   - Uses frozen Step12D normalized RCTD Endothelial weights.
#   - Uses frozen Step12B BPCells spatial counts.
#   - The displayed signal is an ENDOTHELIAL-WEIGHTED spatial
#     response-program score, not a claim of single-cell identity.
#   - No new clustering, no LR rerun, no new gene-set discovery.
# ============================================================

options(stringsAsFactors = FALSE, warn = 1)

suppressPackageStartupMessages({
  library(BPCells)
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(png)
  library(scales)
  library(grid)
  library(patchwork)
})

# ------------------------------------------------------------
# 0. User settings
# ------------------------------------------------------------
ROOT_CANDIDATES <- unique(c(Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE)), "."))

SECTION_7D <- "106_C_7days_2"
SECTION_2M <- "106_C_2months_2"

OUTPUT_DIR_NAME <- "09_final_figures/GLIA/Figure5/single_panels/Figure5D_v5_rowlabel_fix"

BASE_FAMILY <- "Arial"
DPI_OUT <- 600

# H&E visibility
BACKGROUND_STRENGTH <- 0.15

# Spatial spots
SPOT_SIZE   <- 1.36
SPOT_STROKE <- 0.035
SPOT_BORDER <- "#202020"

# Crop
CROP_MARGIN_X_PX <- 20
CROP_MARGIN_Y_PX <- 8

# Deterministic log-normalization
SCALE_FACTOR <- 10000

# Continuous spatial palette, matching Figure 5A/B visual family
PROGRAM_COLORS <- c(
  "#090909",
  "#21102F",
  "#4A176D",
  "#7A287D",
  "#B63D80",
  "#E65B62",
  "#F58B45",
  "#F6D66A"
)

# Robust display range shared across 7 dpi and 2 months
DISPLAY_QUANTILES <- c(0.02, 0.98)

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

first_existing <- function(paths) {
  hit <- paths[file.exists(paths)][1]
  if (length(hit) == 0 || is.na(hit)) return(NA_character_)
  normalizePath(hit, winslash = "/", mustWork = TRUE)
}

read_png_gz <- function(path) {
  con <- gzfile(path, open = "rb")
  on.exit(close(con), add = TRUE)
  raw_vec <- readBin(con, what = "raw", n = 1e8)
  png::readPNG(raw_vec)
}

lighten_image <- function(img, strength = 0.16) {
  strength <- max(0, min(1, strength))
  out <- img

  if (length(dim(out)) != 3 || dim(out)[3] < 3) return(out)

  out[, , 1:3] <- 1 - strength * (1 - out[, , 1:3])

  if (dim(out)[3] == 4) {
    out[, , 4] <- 1
  }

  out
}

read_hires_scalef <- function(path) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)

  txt <- paste(readLines(con, warn = FALSE), collapse = "")

  m <- regexec(
    '"tissue_hires_scalef"\\s*:\\s*([0-9eE+\\.-]+)',
    txt
  )

  hit <- regmatches(txt, m)[[1]]

  if (length(hit) < 2) {
    stop2("Could not parse tissue_hires_scalef from: ", path)
  }

  as.numeric(hit[2])
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

# Same deterministic gene-wise z-score concept used for recipient
# response programs in Figure 4, but applied to spatial log-normalized
# spot expression.
program_score <- function(log_expr_gene_by_spot) {
  z <- t(scale(t(as.matrix(log_expr_gene_by_spot))))
  z[!is.finite(z)] <- 0
  colMeans(z, na.rm = TRUE)
}

robust_limits <- function(x) {
  lim <- quantile(
    x,
    probs = DISPLAY_QUANTILES,
    na.rm = TRUE
  )

  lim <- unname(lim)

  if (!all(is.finite(lim)) || diff(lim) <= 0) {
    lim <- range(x, na.rm = TRUE)
  }

  if (!all(is.finite(lim)) || diff(lim) <= 0) {
    lim <- c(-1, 1)
  }

  lim
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

STEP12A2_DIR <- file.path(
  ROOT_RESULTS,
  "step12A2_GSE234774_spatial_metadata_asset_linkage_repair"
)

ASSET_MANIFEST_FILE <- file.path(
  STEP12A2_DIR,
  "02_2d_capture_area_asset_manifest.csv"
)

ASSET_DIR <- file.path(
  STEP12A2_DIR,
  "audit_extracted_2d_assets"
)

required_files <- c(
  STEP12D_WEIGHTS_RDS,
  STEP12F_SPOT_RDS,
  STEP8C_GENESET_TSV,
  ASSET_MANIFEST_FILE
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop2(
    "Required frozen input(s) are missing:\n",
    paste0("  - ", missing_files, collapse = "\n")
  )
}

if (!dir.exists(SPATIAL_BPCELLS_DIR)) {
  stop2(
    "Spatial BPCells directory is missing:\n  ",
    SPATIAL_BPCELLS_DIR
  )
}

# ------------------------------------------------------------
# 3. Load frozen inputs
# ------------------------------------------------------------
message("Loading frozen Figure 5D inputs ...")

spatial_counts <- BPCells::open_matrix_dir(
  SPATIAL_BPCELLS_DIR
)

weights <- as.matrix(
  readRDS(STEP12D_WEIGHTS_RDS)
)

spot_tbl <- as.data.frame(
  readRDS(STEP12F_SPOT_RDS)
)

gs_tbl <- read.delim(
  STEP8C_GENESET_TSV,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

asset_manifest <- readr::read_csv(
  ASSET_MANIFEST_FILE,
  show_col_types = FALSE
)

# ------------------------------------------------------------
# 4. Validate / recover endothelial RCTD column
# ------------------------------------------------------------
endo_candidates <- colnames(weights)[
  grepl(
    "^Endothelial$|endothelial",
    colnames(weights),
    ignore.case = TRUE
  )
]

if (length(endo_candidates) == 0) {
  stop2(
    "No Endothelial column was found in normalized Step12D RCTD weights.\n",
    "Observed columns:\n  ",
    paste(colnames(weights), collapse = ", ")
  )
}

# Prefer exact match.
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

message(
  "ENDOTHELIAL_ACTIVATION genes in spatial data: ",
  length(genes_activation),
  "/",
  length(genes_activation_requested)
)

message(
  "TNF_NFKB_RESPONSE genes in spatial data: ",
  length(genes_tnf),
  "/",
  length(genes_tnf_requested)
)

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

  barcodes <- as.character(sec$barcode)

  # Strict barcode overlap with spatial counts and RCTD weights.
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

  all_program_genes <- unique(
    c(
      genes_activation,
      genes_tnf
    )
  )

  # Library sizes from the full spatial matrix.
  # BPCells evaluates this lazily without loading the complete matrix into RAM.
  section_counts_all <- spatial_counts[
    ,
    barcodes,
    drop = FALSE
  ]

  library_size <- as.numeric(
    colSums(section_counts_all)
  )

  names(library_size) <- barcodes

  if (any(!is.finite(library_size)) || any(library_size <= 0)) {
    stop2(
      "Invalid spatial library size encountered in section: ",
      section_id
    )
  }

  # Only materialize the small program-gene subset.
  prog_counts <- as.matrix(
    spatial_counts[
      all_program_genes,
      barcodes,
      drop = FALSE
    ]
  )

  # Deterministic log-normalization.
  log_expr <- log1p(
    sweep(
      prog_counts,
      2,
      library_size,
      "/"
    ) * SCALE_FACTOR
  )

  activation_score <- program_score(
    log_expr[
      genes_activation,
      ,
      drop = FALSE
    ]
  )

  tnf_score <- program_score(
    log_expr[
      genes_tnf,
      ,
      drop = FALSE
    ]
  )

  endo_weight <- as.numeric(
    weights[
      barcodes,
      ENDO_COL
    ]
  )

  names(endo_weight) <- barcodes

  # Endothelial-weighted response signal.
  #
  # We first calculate the deterministic response-program score in each spot,
  # then attenuate it by the frozen normalized endothelial RCTD fraction.
  # Finally, z-standardize within section for display so that the maps show
  # the spatial distribution of response activity rather than differences
  # in arbitrary score scale between sections.
  activation_weighted_raw <- activation_score * endo_weight
  tnf_weighted_raw <- tnf_score * endo_weight

  sec$endothelial_weight <- endo_weight[sec$barcode]

  sec$endothelial_activation_weighted <- safe_z(
    activation_weighted_raw[sec$barcode]
  )

  sec$endothelial_tnf_nfkb_weighted <- safe_z(
    tnf_weighted_raw[sec$barcode]
  )

  sec
}

sec_7d <- score_section(SECTION_7D)
sec_2m <- score_section(SECTION_2M)

plot_data <- bind_rows(
  sec_7d %>% mutate(display_timepoint = "7 dpi"),
  sec_2m %>% mutate(display_timepoint = "2 months")
)

# ------------------------------------------------------------
# 7. Shared display limits across timepoints, separately by program
# ------------------------------------------------------------
activation_lims <- robust_limits(
  plot_data$endothelial_activation_weighted
)

tnf_lims <- robust_limits(
  plot_data$endothelial_tnf_nfkb_weighted
)

message(
  "Activation display limits: ",
  paste(round(activation_lims, 3), collapse = " to ")
)

message(
  "TNF/NF-kB display limits: ",
  paste(round(tnf_lims, 3), collapse = " to ")
)

# ------------------------------------------------------------
# 8. Image / coordinate helpers
# ------------------------------------------------------------
prepare_section_image <- function(sec_df) {

  capture_area <- unique(sec_df$capture_area_id)

  if (length(capture_area) != 1) {
    stop2(
      "Expected one capture_area_id for section ",
      unique(sec_df$section_id)
    )
  }

  asset_row <- asset_manifest %>%
    dplyr::filter(
      .data$capture_area_id == capture_area
    )

  if (nrow(asset_row) != 1) {
    stop2(
      "Expected exactly one asset row for capture area: ",
      capture_area
    )
  }

  hires_file <- first_existing(
    c(
      file.path(
        ASSET_DIR,
        asset_row$image_hires_member[1]
      ),
      file.path(
        ROOT,
        asset_row$image_hires_member[1]
      )
    )
  )

  scale_file <- first_existing(
    c(
      file.path(
        ASSET_DIR,
        asset_row$scalefactor_member[1]
      ),
      file.path(
        ROOT,
        asset_row$scalefactor_member[1]
      )
    )
  )

  if (is.na(hires_file)) {
    stop2(
      "Could not locate H&E hires image for capture area: ",
      capture_area
    )
  }

  if (is.na(scale_file)) {
    stop2(
      "Could not locate scalefactor JSON.gz for capture area: ",
      capture_area
    )
  }

  img_original <- read_png_gz(
    hires_file
  )

  img_light <- lighten_image(
    img_original,
    BACKGROUND_STRENGTH
  )

  hires_scalef <- read_hires_scalef(
    scale_file
  )

  img_h <- dim(img_original)[1]
  img_w <- dim(img_original)[2]

  sec_df <- sec_df %>%
    dplyr::mutate(
      x_hires = coordinate_x_fullres * hires_scalef,
      y_hires = coordinate_y_fullres * hires_scalef
    )

  xmin <- max(
    0,
    min(sec_df$x_hires, na.rm = TRUE) - CROP_MARGIN_X_PX
  )

  xmax <- min(
    img_w,
    max(sec_df$x_hires, na.rm = TRUE) + CROP_MARGIN_X_PX
  )

  ymin <- max(
    0,
    min(sec_df$y_hires, na.rm = TRUE) - CROP_MARGIN_Y_PX
  )

  ymax <- min(
    img_h,
    max(sec_df$y_hires, na.rm = TRUE) + CROP_MARGIN_Y_PX
  )

  list(
    df = sec_df,
    image = img_light,
    img_w = img_w,
    img_h = img_h,
    xlim = c(xmin, xmax),
    ylim = c(ymax, ymin),
    hires_file = hires_file,
    scale_file = scale_file
  )
}

img_7d <- prepare_section_image(
  plot_data %>%
    dplyr::filter(section_id == SECTION_7D)
)

img_2m <- prepare_section_image(
  plot_data %>%
    dplyr::filter(section_id == SECTION_2M)
)

# ------------------------------------------------------------
# 9. Spatial map function
# ------------------------------------------------------------
make_map <- function(
  prep,
  score_col,
  score_limits,
  title_text = NULL,
  show_legend = TRUE
) {

  p <- ggplot(
    prep$df,
    aes(
      x = x_hires,
      y = y_hires,
      fill = .data[[score_col]]
    )
  ) +

    annotation_raster(
      as.raster(prep$image),
      xmin = 0,
      xmax = prep$img_w,
      ymin = prep$img_h,
      ymax = 0
    ) +

    geom_point(
      shape = 21,
      size = SPOT_SIZE,
      colour = SPOT_BORDER,
      stroke = SPOT_STROKE,
      alpha = 1
    ) +

    scale_fill_gradientn(
      colours = PROGRAM_COLORS,
      limits = score_limits,
      oob = scales::squish,
      breaks = scales::pretty_breaks(n = 4),
      name = "Weighted score (z)"
    ) +

    coord_equal(
      xlim = prep$xlim,
      ylim = prep$ylim,
      expand = FALSE,
      clip = "on"
    ) +

    labs(
      title = title_text
    ) +

    theme_void(
      base_family = BASE_FAMILY
    ) +

    theme(
      plot.background = element_rect(
        fill = "white",
        colour = NA
      ),

      panel.background = element_rect(
        fill = "white",
        colour = NA
      ),

      plot.title = element_text(
        family = BASE_FAMILY,
        size = 11.5,
        face = "plain",
        hjust = 0.5,
        margin = margin(b = 2)
      ),

      legend.position = if (show_legend) "right" else "none",

      legend.title = element_text(
        family = BASE_FAMILY,
        size = 10.5,
        face = "plain"
      ),

      legend.text = element_text(
        family = BASE_FAMILY,
        size = 9.8
      ),

      legend.key.height = unit(4.0, "mm"),

      plot.margin = margin(0.15, 0.15, 0.15, 0.15)
    )

  if (show_legend) {
    p <- p +
      guides(
        fill = guide_colorbar(
          title.position = "top",
          title.hjust = 0,
          barheight = unit(19, "mm"),
          barwidth = unit(3.5, "mm"),
          frame.colour = "#BDBDBD",
          ticks.colour = "#555555"
        )
      )
  }

  p
}

# ------------------------------------------------------------
# 10. Four maps
# ------------------------------------------------------------
p_7d_activation <- make_map(
  img_7d,
  "endothelial_activation_weighted",
  activation_lims,
  NULL,
  TRUE
)

p_7d_tnf <- make_map(
  img_7d,
  "endothelial_tnf_nfkb_weighted",
  tnf_lims,
  NULL,
  TRUE
)

p_2m_activation <- make_map(
  img_2m,
  "endothelial_activation_weighted",
  activation_lims,
  NULL,
  TRUE
)

p_2m_tnf <- make_map(
  img_2m,
  "endothelial_tnf_nfkb_weighted",
  tnf_lims,
  NULL,
  TRUE
)

# ------------------------------------------------------------
# 11. Titles + row labels + final layout
# ------------------------------------------------------------
make_header <- function(label) {
  ggplot() +
    annotate(
      "text",
      x = 0.5,
      y = 0.5,
      label = label,
      size = 13.0 / .pt,
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
          x = 0.10,
          y = 0.70,
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
      x = 0.62,
      y = 0.50,
      label = label,
      angle = 90,
      size = 12.0 / .pt,
      family = BASE_FAMILY,
      fontface = "bold"
    ) +
    xlim(0, 1) + ylim(0, 1) +
    theme_void()
}

header_row <- make_blank_header("D") +
  make_header("Endothelial activation") +
  make_header("TNF/NF\u03baB response") +
  patchwork::plot_layout(widths = c(0.14, 1, 1))

label_7d <- make_row_label("7 dpi")
label_2m <- make_row_label("2 months")

row_7d <- label_7d + p_7d_activation + p_7d_tnf +
  patchwork::plot_layout(widths = c(0.14, 1, 1))

row_2m <- label_2m + p_2m_activation + p_2m_tnf +
  patchwork::plot_layout(widths = c(0.14, 1, 1))

final_plot <- header_row / row_7d / row_2m +
  patchwork::plot_layout(
    heights = c(0.10, 1, 1)
  )

# ------------------------------------------------------------
# 12. Outputs
# ------------------------------------------------------------
out_dir <- file.path(
  ROOT,
  OUTPUT_DIR_NAME
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

out_base <- file.path(
  out_dir,
  "Figure5D_endothelial_response_spatial_maps_v5_rowlabel_fix"
)

pdf_device <- if (capabilities("cairo")) {
  cairo_pdf
} else {
  grDevices::pdf
}

ggsave(
  paste0(out_base, ".pdf"),
  final_plot,
  width = 7.8,
  height = 3.2,
  units = "in",
  device = pdf_device,
  bg = "white"
)

ggsave(
  paste0(out_base, ".png"),
  final_plot,
  width = 7.8,
  height = 3.2,
  units = "in",
  dpi = DPI_OUT,
  bg = "white"
)

if (requireNamespace("svglite", quietly = TRUE)) {
  ggsave(
    paste0(out_base, ".svg"),
    final_plot,
    width = 7.8,
    height = 3.2,
    units = "in",
    bg = "white"
  )
}

# ------------------------------------------------------------
# 13. Plot data + audit
# ------------------------------------------------------------
plot_data_out <- bind_rows(
  img_7d$df,
  img_2m$df
) %>%
  dplyr::select(
    tidyselect::any_of(
      c(
        "barcode",
        "section_id",
        "capture_area_id",
        "biological_replicate_id",
        "timepoint",
        "coordinate_x_fullres",
        "coordinate_y_fullres",
        "x_hires",
        "y_hires",
        "endothelial_weight",
        "endothelial_activation_weighted",
        "endothelial_tnf_nfkb_weighted",
        "target_program_adjusted_z_within_section",
        "distance_to_core_NN"
      )
    )
  )

readr::write_csv(
  plot_data_out,
  paste0(
    out_base,
    "_plot_data.csv"
  )
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
  paste0(
    out_base,
    "_gene_audit.csv"
  )
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
    "activation_display_low",
    "activation_display_high",
    "TNF_NFkB_display_low",
    "TNF_NFkB_display_high",
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
    "weighted signal z-standardized within section",
    signif(activation_lims[1], 5),
    signif(activation_lims[2], 5),
    signif(tnf_lims[1], 5),
    signif(tnf_lims[2], 5),
    "endothelial-weighted spatial response signal; not single-cell proof and not direct TNF-TNFRSF1A binding"
  )
)

readr::write_csv(
  audit,
  paste0(
    out_base,
    "_audit.csv"
  )
)

message("============================================================")
message("Figure 5D v5 row-label fix completed.")
message("7 dpi section : ", SECTION_7D)
message("2 months      : ", SECTION_2M)
message("Output PDF    : ", paste0(out_base, ".pdf"))
message("============================================================")
