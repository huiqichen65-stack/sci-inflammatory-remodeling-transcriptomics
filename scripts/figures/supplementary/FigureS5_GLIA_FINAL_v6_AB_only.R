# ==============================================================================
# FigureS5_GLIA_FINAL_v6_AB_only.R
# GLIA Supplementary Figure S5
# Spatial robustness of lesion-associated Target-endothelial coupling
#
# Revised 1 x 2 layout
#   A  All-animal lesion-zone validation with Shared control
#   B  Section/animal-level endothelial response coupling
#
# IMPORTANT statistical hierarchy
#   - spots are used only to construct section-level quantities
#   - sections are summarized within biological animals
#   - biological animal is the top-level replicate
#
# Panel B is NOT the two representative Figure 5E sections.
# It recomputes ENDOTHELIAL_ACTIVATION and TNF_NFKB_RESPONSE across ALL injured
# sections (7 dpi and 2 months), then summarizes section correlations within
# biological animals.
# ==============================================================================

options(stringsAsFactors = FALSE, warn = 1)

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(patchwork)
  library(scales)
  library(grid)
  library(BPCells)
})

# ------------------------------------------------------------------------------
# 0. User settings
# ------------------------------------------------------------------------------

ROOT_CANDIDATES <- unique(c(Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE)), "."))

OUTPUT_DIR_NAME <- "09_final_figures/GLIA/FigureS5_FINAL_v6_AB_only"

BASE_FAMILY <- "Arial"
DPI_OUT <- 600
SCALE_FACTOR <- 10000

TIMEPOINT_ORDER <- c("7days", "2months")
TIMEPOINT_LABELS <- c(
  "7days" = "7 dpi",
  "2months" = "2 months"
)

ZONE_ORDER <- c("core", "border", "perilesional")
ZONE_LABELS <- c(
  "core" = "Core",
  "border" = "Border",
  "perilesional" = "Peri."
)

# Main-figure-compatible semantic colours
COL_TARGET <- "#C44E52"
COL_SHARED <- "#8172B2"
COL_ENDO   <- "#4C72B0"
COL_TNF    <- "#D95F02"
COL_TEXT   <- "#202020"
COL_ZERO   <- "#A8A8A8"
COL_GRID   <- "#E7E7E7"
COL_LIGHT  <- "#F4F5F7"
COL_QC     <- "#4C72B0"

# ------------------------------------------------------------------------------
# 1. Helpers
# ------------------------------------------------------------------------------

stop2 <- function(...) stop(..., call. = FALSE)

locate_root <- function() {
  for (root in ROOT_CANDIDATES) {
    if (!dir.exists(root)) next

    key1 <- file.path(
      root,
      "04_results",
      "step12F_GSE234774_independent_lesion_zones_target_program_validation",
      "tables",
      "05_animal_zone_vs_remote_contrasts.csv"
    )

    if (file.exists(key1)) {
      return(normalizePath(root, winslash = "/", mustWork = TRUE))
    }
  }

  stop2(
    "Could not locate the xjz project.\n",
    "Edit ROOT_CANDIDATES near the top of this script."
  )
}

must_file <- function(path) {
  if (!file.exists(path)) {
    stop2("Required file is missing:\n  ", path)
  }
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

must_dir <- function(path) {
  if (!dir.exists(path)) {
    stop2("Required directory is missing:\n  ", path)
  }
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

safe_num <- function(x) suppressWarnings(as.numeric(x))

safe_z <- function(x) {
  x <- as.numeric(x)
  s <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(s) || s <= 0) return(rep(0, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

program_score <- function(log_expr_gene_by_spot) {
  z <- t(scale(t(as.matrix(log_expr_gene_by_spot))))
  z[!is.finite(z)] <- 0
  colMeans(z, na.rm = TRUE)
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
    dplyr::filter(.data$response_program == .env$program_name)

  if (nrow(row) != 1) {
    stop2(
      "Could not uniquely identify response program: ",
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

pretty_pct <- function(x) paste0(round(100 * x, 1), "%")

theme_s6 <- function(base_size = 8.2) {
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
        size = 8.5,
        colour = COL_TEXT
      ),
      axis.text = element_text(
        family = BASE_FAMILY,
        size = 7.5,
        colour = COL_TEXT
      ),
      axis.line = element_line(
        linewidth = 0.45,
        colour = "#333333"
      ),
      axis.ticks = element_line(
        linewidth = 0.40,
        colour = "#333333"
      ),
      strip.background = element_blank(),
      strip.text = element_text(
        family = BASE_FAMILY,
        face = "bold",
        size = 8.3,
        colour = COL_TEXT,
        margin = margin(b = 2)
      ),
      legend.text = element_text(
        family = BASE_FAMILY,
        size = 7.4
      ),
      legend.title = element_blank(),
      plot.margin = margin(2.5, 3.5, 3.0, 3.5, unit = "mm")
    )
}

# Figure-S4-like header geometry:
# - panel letter is complete and at the far upper-left of the WHOLE panel
# - the letter sits slightly ABOVE the centered title
# - both are inside a full-width header strip, so child-axis widths cannot move them
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
          fontsize = 7.8,
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

# ------------------------------------------------------------------------------
# 2. Frozen input paths
# ------------------------------------------------------------------------------

ROOT <- locate_root()
OUT <- file.path(ROOT, OUTPUT_DIR_NAME)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# ---- Shared section manifest (needed by new Panel B metadata fallback)
A_SECTION_MANIFEST <- must_file(file.path(
  ROOT,
  "04_results/step12A2_GSE234774_spatial_metadata_asset_linkage_repair",
  "03_condition_specific_section_manifest.csv"
))


# ---- New Panel A
B_ANIMAL_CONTRAST <- must_file(file.path(
  ROOT,
  "04_results/step12F_GSE234774_independent_lesion_zones_target_program_validation",
  "tables",
  "05_animal_zone_vs_remote_contrasts.csv"
))

B_ANIMAL_INFERENCE <- must_file(file.path(
  ROOT,
  "04_results/step12F_GSE234774_independent_lesion_zones_target_program_validation",
  "tables",
  "06_animal_level_inference_summary.csv"
))


# ---- New Panel B: exact all-section response rescoring
D_SPATIAL_COUNTS_DIR <- must_dir(file.path(
  ROOT,
  "04_results/step12B_fast_BPCells_spatial_section_QC",
  "BPCells",
  "GSE234774_spatial_2d_full_counts"
))

D_RCTD_WEIGHTS_RDS <- must_file(file.path(
  ROOT,
  "04_results/step12D_GSE234774_full_mode_RCTD_per_capture_area",
  "RDS",
  "GSE234774_step12D_all_spots_normalized_RCTD_weights.rds"
))

D_STEP12F_SPOT_RDS <- must_file(file.path(
  ROOT,
  "04_results/step12F_GSE234774_independent_lesion_zones_target_program_validation",
  "RDS",
  "GSE234774_step12F_spot_independent_lesion_zones_and_scores.rds"
))

D_STEP8C_GENESET <- must_file(file.path(
  ROOT,
  "04_results/step8C_GSE162610_recipient_program_integration",
  "tables",
  "02_recipient_program_gene_coverage.tsv"
))

input_manifest <- tibble::tibble(
  panel = c(
    "Shared metadata",
    "A", "A",
    "B", "B", "B", "B"
  ),
  input = c(
    A_SECTION_MANIFEST,
    B_ANIMAL_CONTRAST,
    B_ANIMAL_INFERENCE,
    D_SPATIAL_COUNTS_DIR,
    D_RCTD_WEIGHTS_RDS,
    D_STEP12F_SPOT_RDS,
    D_STEP8C_GENESET
  )
)

readr::write_csv(
  input_manifest,
  file.path(OUT, "01_FigureS5_input_manifest.csv")
)

# ==============================================================================
# 3. SHARED METADATA | section manifest retained for Panel B metadata fallback
# ==============================================================================

section_manifest <- readr::read_csv(
  A_SECTION_MANIFEST,
  show_col_types = FALSE
)

required_section_manifest <- c(
  "section_id",
  "timepoint",
  "replicate_index",
  "biological_replicate_id"
)

miss_manifest <- setdiff(required_section_manifest, names(section_manifest))
if (length(miss_manifest) > 0) {
  stop2(
    "Section manifest missing columns: ",
    paste(miss_manifest, collapse = ", ")
  )
}

# ==============================================================================
# 4. PANEL A | All-animal lesion-zone validation with Shared control
# ==============================================================================

animal_contrast <- readr::read_csv(
  B_ANIMAL_CONTRAST,
  show_col_types = FALSE
)

animal_inference <- readr::read_csv(
  B_ANIMAL_INFERENCE,
  show_col_types = FALSE
)

b_dat <- animal_contrast %>%
  dplyr::filter(
    .data$timepoint %in% TIMEPOINT_ORDER,
    .data$metric %in% c(
      "target_adjusted_within_section_z",
      "shared_adjusted_z"
    ),
    .data$zone %in% ZONE_ORDER
  ) %>%
  dplyr::mutate(
    program = dplyr::recode(
      .data$metric,
      target_adjusted_within_section_z = "Target program",
      shared_adjusted_z = "Shared prebranch"
    ),
    zone = factor(
      .data$zone,
      levels = ZONE_ORDER,
      labels = unname(ZONE_LABELS[ZONE_ORDER])
    ),
    timepoint = factor(
      .data$timepoint,
      levels = TIMEPOINT_ORDER
    )
  )

if (nrow(b_dat) == 0) {
  stop2("Panel A filtering returned no rows.")
}

b_mean <- b_dat %>%
  dplyr::group_by(
    .data$timepoint,
    .data$zone,
    .data$program
  ) %>%
  dplyr::summarise(
    n_animals = dplyr::n_distinct(.data$biological_replicate_id),
    mean_contrast = mean(.data$animal_mean_contrast, na.rm = TRUE),
    .groups = "drop"
  )

program_cols <- c(
  "Target program" = COL_TARGET,
  "Shared prebranch" = COL_SHARED
)

pB <- ggplot() +
  geom_hline(
    yintercept = 0,
    linewidth = 0.45,
    linetype = "dashed",
    colour = COL_ZERO
  ) +
  geom_line(
    data = b_dat,
    aes(
      x = zone,
      y = animal_mean_contrast,
      colour = program,
      group = interaction(program, biological_replicate_id)
    ),
    linewidth = 0.55,
    alpha = 0.28
  ) +
  geom_point(
    data = b_dat,
    aes(
      x = zone,
      y = animal_mean_contrast,
      colour = program
    ),
    size = 1.65,
    alpha = 0.55
  ) +
  geom_line(
    data = b_mean,
    aes(
      x = zone,
      y = mean_contrast,
      colour = program,
      group = program
    ),
    linewidth = 1.35,
    lineend = "round"
  ) +
  geom_point(
    data = b_mean,
    aes(
      x = zone,
      y = mean_contrast,
      colour = program
    ),
    size = 2.55
  ) +
  facet_grid(
    . ~ timepoint,
    labeller = as_labeller(TIMEPOINT_LABELS)
  ) +
  scale_colour_manual(
    values = program_cols,
    breaks = c("Target program", "Shared prebranch")
  ) +
  labs(
    x = NULL,
    y = "Animal mean zone − remote contrast (z)"
  ) +
  theme_s6() +
  theme(
    legend.position = "top",
    legend.justification = "center",
    legend.direction = "horizontal",
    legend.margin = margin(0, 0, 1, 0),
    panel.spacing.x = unit(5, "mm"),
    axis.text.x = element_text(size = 6.8)
  )

panelA <- assemble_panel(
  body = pB,
  letter = "A",
  title = "All-animal lesion-zone validation"
)

# ==============================================================================
# 5. PANEL B | Section/animal-level endothelial response coupling
# ==============================================================================

message("Panel B: loading frozen spatial inputs ...")

spatial_counts <- BPCells::open_matrix_dir(
  D_SPATIAL_COUNTS_DIR
)

weights <- as.matrix(
  readRDS(D_RCTD_WEIGHTS_RDS)
)

spot_tbl <- as.data.frame(
  readRDS(D_STEP12F_SPOT_RDS)
)

gs_tbl <- read.delim(
  D_STEP8C_GENESET,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# ----- Required Step12F fields
required_D_spot <- c(
  "barcode",
  "section_id",
  "target_program_adjusted_z_within_section"
)

miss_D <- setdiff(required_D_spot, names(spot_tbl))
if (length(miss_D) > 0) {
  stop2(
    "Panel B Step12F spot table missing: ",
    paste(miss_D, collapse = ", ")
  )
}

# Bring timepoint/animal from manifest if they are not already present.
if (!all(c("timepoint", "biological_replicate_id") %in% names(spot_tbl))) {
  spot_tbl <- spot_tbl %>%
    dplyr::left_join(
      section_manifest %>%
        dplyr::select(
          section_id,
          timepoint_manifest = timepoint,
          biological_replicate_id_manifest = biological_replicate_id
        ),
      by = "section_id"
    )

  if (!"timepoint" %in% names(spot_tbl)) {
    spot_tbl$timepoint <- spot_tbl$timepoint_manifest
  }

  if (!"biological_replicate_id" %in% names(spot_tbl)) {
    spot_tbl$biological_replicate_id <-
      spot_tbl$biological_replicate_id_manifest
  }
}

# ----- Endothelial RCTD column
endo_candidates <- colnames(weights)[
  grepl(
    "^Endothelial$|endothelial",
    colnames(weights),
    ignore.case = TRUE
  )
]

if (length(endo_candidates) == 0) {
  stop2(
    "No Endothelial column found in Step12D RCTD weights."
  )
}

ENDO_COL <- if ("Endothelial" %in% endo_candidates) {
  "Endothelial"
} else {
  endo_candidates[1]
}

# ----- Frozen response programs from Step8C
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
    "Too few ENDOTHELIAL_ACTIVATION genes in spatial counts: ",
    length(genes_activation)
  )
}

if (length(genes_tnf) < 3) {
  stop2(
    "Too few TNF_NFKB_RESPONSE genes in spatial counts: ",
    length(genes_tnf)
  )
}

d_spots <- spot_tbl %>%
  dplyr::filter(
    .data$timepoint %in% TIMEPOINT_ORDER
  ) %>%
  dplyr::mutate(
    barcode = as.character(.data$barcode),
    section_id = as.character(.data$section_id),
    biological_replicate_id =
      as.character(.data$biological_replicate_id)
  )

if (nrow(d_spots) == 0) {
  stop2("Panel B found no 7 dpi / 2 months Step12F spots.")
}

all_barcodes <- unique(d_spots$barcode)

missing_counts <- setdiff(
  all_barcodes,
  colnames(spatial_counts)
)

missing_weights <- setdiff(
  all_barcodes,
  rownames(weights)
)

if (length(missing_counts) > 0) {
  stop2(
    "Panel B spatial counts missing ",
    length(missing_counts),
    " Step12F barcodes."
  )
}

if (length(missing_weights) > 0) {
  stop2(
    "Panel B RCTD weights missing ",
    length(missing_weights),
    " Step12F barcodes."
  )
}

all_program_genes <- unique(
  c(genes_activation, genes_tnf)
)

message(
  "Panel B | injured spots = ",
  length(all_barcodes),
  " | program genes = ",
  length(all_program_genes)
)

# ----- Efficient one-pass normalization
#
# Library size is computed once for all injured spots.
# Only the small frozen program-gene matrix is materialized.
message("Panel B: computing spatial library sizes ...")

library_size <- as.numeric(
  colSums(
    spatial_counts[, all_barcodes, drop = FALSE]
  )
)

names(library_size) <- all_barcodes

if (any(!is.finite(library_size)) ||
    any(library_size <= 0)) {
  stop2("Panel B encountered invalid spatial library sizes.")
}

message("Panel B: loading frozen response-program genes ...")

prog_counts <- as.matrix(
  spatial_counts[
    all_program_genes,
    all_barcodes,
    drop = FALSE
  ]
)

log_expr <- log1p(
  sweep(
    prog_counts,
    2,
    library_size[colnames(prog_counts)],
    "/"
  ) * SCALE_FACTOR
)

# ----- Score each section using the SAME logic as Figure 5E
section_ids <- sort(unique(d_spots$section_id))

section_results <- vector(
  "list",
  length(section_ids)
)

for (ii in seq_along(section_ids)) {

  sid <- section_ids[ii]

  sec <- d_spots %>%
    dplyr::filter(.data$section_id == .env$sid)

  barcodes <- sec$barcode

  if (length(barcodes) < 20) {
    warning(
      "Skipping very small section in Panel D: ",
      sid,
      " | n = ",
      length(barcodes)
    )
    next
  }

  sec_log <- log_expr[, barcodes, drop = FALSE]

  activation_raw <- program_score(
    sec_log[
      genes_activation,
      ,
      drop = FALSE
    ]
  )

  tnf_raw <- program_score(
    sec_log[
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

  target_score <- as.numeric(
    sec$target_program_adjusted_z_within_section
  )

  activation_weighted <- safe_z(
    activation_raw[barcodes] * endo_weight[barcodes]
  )

  tnf_weighted <- safe_z(
    tnf_raw[barcodes] * endo_weight[barcodes]
  )

  rho_activation <- suppressWarnings(
    cor(
      target_score,
      activation_weighted,
      method = "spearman",
      use = "complete.obs"
    )
  )

  rho_tnf <- suppressWarnings(
    cor(
      target_score,
      tnf_weighted,
      method = "spearman",
      use = "complete.obs"
    )
  )

  section_results[[ii]] <- tibble::tibble(
    section_id = sid,
    timepoint = unique(sec$timepoint)[1],
    biological_replicate_id =
      unique(sec$biological_replicate_id)[1],
    n_spots = nrow(sec),
    response_program = c(
      "Endothelial activation",
      "TNF/NFκB response"
    ),
    spearman_rho = c(
      rho_activation,
      rho_tnf
    )
  )
}

d_section <- dplyr::bind_rows(section_results) %>%
  dplyr::filter(
    is.finite(.data$spearman_rho)
  ) %>%
  dplyr::mutate(
    timepoint = factor(
      .data$timepoint,
      levels = TIMEPOINT_ORDER
    ),
    response_program = factor(
      .data$response_program,
      levels = c(
        "Endothelial activation",
        "TNF/NFκB response"
      )
    )
  )

if (nrow(d_section) == 0) {
  stop2("Panel B produced no valid section-level correlations.")
}

# Fisher-z average of section correlations within each biological animal.
d_animal <- d_section %>%
  dplyr::mutate(
    fisher_z = atanh(
      pmin(
        pmax(.data$spearman_rho, -0.999999),
        0.999999
      )
    )
  ) %>%
  dplyr::group_by(
    .data$timepoint,
    .data$biological_replicate_id,
    .data$response_program
  ) %>%
  dplyr::summarise(
    n_sections = dplyr::n(),
    animal_rho = tanh(
      mean(.data$fisher_z, na.rm = TRUE)
    ),
    min_section_rho = min(
      .data$spearman_rho,
      na.rm = TRUE
    ),
    max_section_rho = max(
      .data$spearman_rho,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

d_group <- d_animal %>%
  dplyr::group_by(
    .data$timepoint,
    .data$response_program
  ) %>%
  dplyr::summarise(
    n_animals = dplyr::n_distinct(
      .data$biological_replicate_id
    ),
    mean_animal_rho = mean(
      .data$animal_rho,
      na.rm = TRUE
    ),
    positive_animals = sum(
      .data$animal_rho > 0,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

response_cols <- c(
  "Endothelial activation" = COL_ENDO,
  "TNF/NFκB response" = COL_TNF
)

# section-level dots are descriptive; larger points are biological animals.
pD <- ggplot() +
  geom_hline(
    yintercept = 0,
    colour = COL_ZERO,
    linewidth = 0.45,
    linetype = "dashed"
  ) +
  geom_point(
    data = d_section,
    aes(
      x = response_program,
      y = spearman_rho,
      colour = response_program
    ),
    size = 1.15,
    alpha = 0.22,
    position = position_jitter(
      width = 0.07,
      height = 0
    )
  ) +
  geom_line(
    data = d_animal,
    aes(
      x = response_program,
      y = animal_rho,
      group = biological_replicate_id
    ),
    colour = "#A0A0A0",
    linewidth = 0.50,
    alpha = 0.55
  ) +
  geom_point(
    data = d_animal,
    aes(
      x = response_program,
      y = animal_rho,
      colour = response_program
    ),
    size = 2.45
  ) +
  geom_point(
    data = d_group,
    aes(
      x = response_program,
      y = mean_animal_rho,
      colour = response_program
    ),
    shape = 18,
    size = 3.9
  ) +
  facet_grid(
    . ~ timepoint,
    labeller = as_labeller(TIMEPOINT_LABELS)
  ) +
  scale_colour_manual(
    values = response_cols,
    breaks = c(
      "Endothelial activation",
      "TNF/NFκB response"
    ),
    labels = c(
      "Endothelial activation",
      "TNF/NFκB response"
    )
  ) +
  scale_x_discrete(
    labels = NULL
  ) +
  labs(
    x = NULL,
    y = "Target–response Spearman ρ"
  ) +
  theme_s6() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "top",
    legend.justification = "center",
    legend.direction = "horizontal",
    legend.key.width = unit(5, "mm"),
    legend.spacing.x = unit(2, "mm"),
    legend.margin = margin(0, 0, 1, 0),
    panel.spacing.x = unit(6, "mm")
  )

panelB <- assemble_panel(
  body = pD,
  letter = "B",
  title = "Animal-level Target–endothelial response coupling"
)

# ==============================================================================
# 6. Final 1 x 2 assembly
# ==============================================================================

final_plot <- patchwork::wrap_plots(
  panelA,
  panelB,
  ncol = 2,
  widths = c(1, 1)
)

# ==============================================================================
# 8. Frozen audit / plot data outputs
# ==============================================================================


readr::write_csv(
  b_dat,
  file.path(
    OUT,
    "02_PanelA_animal_zone_target_shared_plot_data.csv"
  )
)

readr::write_csv(
  b_mean,
  file.path(
    OUT,
    "03_PanelA_group_mean_summary.csv"
  )
)


readr::write_csv(
  d_section,
  file.path(
    OUT,
    "04_PanelB_section_level_response_correlations.csv"
  )
)

readr::write_csv(
  d_animal,
  file.path(
    OUT,
    "05_PanelB_animal_level_response_correlations.csv"
  )
)

readr::write_csv(
  d_group,
  file.path(
    OUT,
    "06_PanelB_group_response_summary.csv"
  )
)

audit <- tibble::tibble(
  item = c(
    "n_spatial_sections_manifest",
    "n_spatial_animals",
    "PanelA_animals_7dpi",
    "PanelA_animals_2months",
    "PanelB_sections_7dpi",
    "PanelB_sections_2months",
    "PanelB_animals_7dpi",
    "PanelB_animals_2months",
    "PanelB_ENDOTHELIAL_ACTIVATION_genes_present",
    "PanelB_TNF_NFKB_RESPONSE_genes_present",
    "PanelB_endothelial_RCTD_column"
  ),
  value = c(
    nrow(section_manifest),
    dplyr::n_distinct(section_manifest$replicate_index),
    dplyr::n_distinct(
      b_dat$biological_replicate_id[
        as.character(b_dat$timepoint) == "7days"
      ]
    ),
    dplyr::n_distinct(
      b_dat$biological_replicate_id[
        as.character(b_dat$timepoint) == "2months"
      ]
    ),
    dplyr::n_distinct(
      d_section$section_id[
        as.character(d_section$timepoint) == "7days"
      ]
    ),
    dplyr::n_distinct(
      d_section$section_id[
        as.character(d_section$timepoint) == "2months"
      ]
    ),
    dplyr::n_distinct(
      d_animal$biological_replicate_id[
        as.character(d_animal$timepoint) == "7days"
      ]
    ),
    dplyr::n_distinct(
      d_animal$biological_replicate_id[
        as.character(d_animal$timepoint) == "2months"
      ]
    ),
    length(genes_activation),
    length(genes_tnf),
    ENDO_COL
  )
)

readr::write_csv(
  audit,
  file.path(
    OUT,
    "00_FigureS6_data_audit.csv"
  )
)

# ==============================================================================
# 8. Export individual panels and combined figure
# ==============================================================================

save_plot_multi(
  panelA,
  file.path(
    OUT,
    "FigureS5A_all_animal_lesion_zone_validation"
  ),
  width = 6.6,
  height = 3.25
)

save_plot_multi(
  panelB,
  file.path(
    OUT,
    "FigureS5B_section_animal_endothelial_response_coupling"
  ),
  width = 6.6,
  height = 3.25
)

save_plot_multi(
  final_plot,
  file.path(
    OUT,
    "FigureS5_GLIA_FINAL_v6_AB_only"
  ),
  width = 7.35,
  height = 3.55
)

# ==============================================================================
# 10. Run summary
# ==============================================================================

summary_lines <- c(
  "GLIA Supplementary Figure S5 | FINAL v6 | A/B only",
  "===================================================",
  "",
  "Core question:",
  paste(
    "Does lesion-associated Target-endothelial coupling remain robust",
    "at the biological-animal level?"
  ),
  "",
  "Panel A | All-animal lesion-zone validation with Shared control",
  "  Thin lines/points = individual biological animals; bold line = mean across animals.",
  "  Biological animal is the displayed replicate.",
  "  Target and Shared-prebranch adjusted scores are shown against remote tissue.",
  "",
  "Panel B | Section/animal-level endothelial response coupling",
  paste0(
    "  ENDOTHELIAL_ACTIVATION spatial genes: ",
    length(genes_activation)
  ),
  paste0(
    "  TNF_NFKB_RESPONSE spatial genes: ",
    length(genes_tnf)
  ),
  paste0(
    "  7 dpi sections: ",
    dplyr::n_distinct(
      d_section$section_id[
        as.character(d_section$timepoint) == "7days"
      ]
    ),
    " | animals: ",
    dplyr::n_distinct(
      d_animal$biological_replicate_id[
        as.character(d_animal$timepoint) == "7days"
      ]
    )
  ),
  paste0(
    "  2 months sections: ",
    dplyr::n_distinct(
      d_section$section_id[
        as.character(d_section$timepoint) == "2months"
      ]
    ),
    " | animals: ",
    dplyr::n_distinct(
      d_animal$biological_replicate_id[
        as.character(d_animal$timepoint) == "2months"
      ]
    )
  ),
  "  Small translucent points = section-level Spearman rho.",
  "  Large circles = Fisher-z-averaged biological-animal rho.",
  "  Diamonds = mean across biological animals.",
  "  No spot-level p-value is used as independent inference.",
  "",
  paste0("Output directory: ", OUT),
  paste0("Created: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
)

writeLines(
  summary_lines,
  con = file.path(
    OUT,
    "FigureS5_GLIA_FINAL_v6_AB_only_run_summary.txt"
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

message(
  "\nDONE\nFigure S5 v6 A/B-only outputs:\n  ",
  OUT
)
