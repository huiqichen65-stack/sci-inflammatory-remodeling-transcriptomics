# ==============================================================================
# 90_assemble_Figure3_GLIA_A_to_F_v13_Ev3_Fv9.R
#
# Figure 3 final assembly | updated E v3 + F v9 | Figure-4-matched letter spacing + GLIA typography
#
# v6 fixes the issues seen in the v5 review PDF:
#
# 1) Keep A/B right edges safe.
# 2) Keep B top-left clean.
# 3) Further enlarge the D-letter mask so the extra original D is fully removed.
# 4) Keep C/E top-left gene labels safe.
# 5) Shrink E further.
# 6) Keep standardized A-F panel letters on the same rows/columns.
#
# STRATEGY
# --------
# - The original embedded letters are covered using SMALLER panel-specific masks,
#   instead of one large universal mask that could cover real content in C/E.
# - Panel bodies are horizontally centered within each column cell.
# - Panel pairs are vertically center-aligned within each row.
# - A/B use a modest upscale factor; E/F use a modest downscale factor.
# - Final canvas keeps a bottom safety margin to avoid clipping.
#
# OUTPUT
# ------
# <PROJECT_ROOT>/09_final_figures/GLIA/
#   Figure3_final_v8_Dmask_Esmaller/
# ==============================================================================

rm(list = ls())
gc()
options(stringsAsFactors = FALSE, scipen = 999)

# ------------------------------------------------------------------------------
# 0. SETTINGS
# ------------------------------------------------------------------------------
PROJECT_ROOT <- Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
SCRIPT_DIR   <- PROJECT_ROOT

OUT_DIR <- file.path(
  PROJECT_ROOT,
  "09_final_figures",
  "GLIA",
  "Figure3_final_v13_Ev3_Fv9"
)

COMBINED_DIR <- file.path(OUT_DIR, "combined")
SINGLE_DIR   <- file.path(OUT_DIR, "single_panels")
LOG_DIR      <- file.path(OUT_DIR, "logs")

dir.create(COMBINED_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(SINGLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)

RUN_PANEL_SCRIPTS_IF_OUTPUT_MISSING <- TRUE

COMBINED_WIDTH_MM <- 180
DPI <- 600L

COMBINED_WIDTH_PX <- as.integer(round(COMBINED_WIDTH_MM / 25.4 * DPI))
GAP_X_PX <- 68L
GAP_Y_PX <- 38L
BOTTOM_SAFETY_PX <- 85L

CELL_WIDTH_PX <- as.integer((COMBINED_WIDTH_PX - GAP_X_PX) / 2)

# Reserve a small internal horizontal safety margin in each column cell so
# enlarged panels never clip at the right edge.
CELL_INNER_X_PAD <- 18L

# A small panel border protects content without creating large blank areas.
PANEL_BORDER_X <- 10L
PANEL_BORDER_Y <- 10L

# Standardized panel letters.
LETTER_FONT   <- "Arial"
LETTER_SIZE   <- 108
LETTER_WEIGHT <- 700
LETTER_COLOR  <- "black"

LETTER_X_LEFT  <- 70L
LETTER_X_RIGHT <- CELL_WIDTH_PX + GAP_X_PX + 12L
LETTER_Y       <- 6L
LETTER_HEADROOM_PX <- 128L
LEFT_PANEL_NUDGE_PX <- 24L

# Panel-specific masks to remove the original embedded panel letters without
# covering nearby labels/genes. C and E intentionally use narrower masks.
LETTER_MASK <- list(
  A = list(w = 118L, h = 92L),
  B = list(w = 154L, h = 112L),
  C = list(w = 78L,  h = 72L),
  D = list(w = 210L, h = 126L),
  E = list(w = 84L,  h = 72L),
  F = list(w = 118L, h = 92L)
)

# Row-/panel-level rescaling:
# - A/B a bit larger
# - C/D baseline
# - E/F a bit smaller
PANEL_SCALE <- c(
  A = 1.02,
  B = 1.02,
  C = 1.00,
  D = 1.00,
  E = 0.90,
  F = 0.95
)

# ------------------------------------------------------------------------------
# 1. REQUIRED PACKAGE
# ------------------------------------------------------------------------------
required_pkgs <- c("magick")

missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_pkgs) > 0L) {
  stop(
    "Missing required package(s): ",
    paste(missing_pkgs, collapse = ", "),
    "\nInstall with install.packages('magick').",
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(magick)
})

# ------------------------------------------------------------------------------
# 2. ACCEPTED PANEL OUTPUT DEFINITIONS
# ------------------------------------------------------------------------------
panel_files <- data.frame(
  panel = c("A", "B", "C", "D", "E", "F"),

  output_dir = file.path(
    PROJECT_ROOT,
    "09_final_figures",
    "GLIA",
    c(
      "Figure3A_v6",
      "Figure3B_v10_GLIA_readability",
      "Figure3C_v5",
      "Figure3D_v5",
      "Figure3E_v3_Mod_labels",
      "Figure3F_v9_black_borders"
    )
  ),

  pdf_name = c(
    "Figure3A_RSS_rankings_v6.pdf",
    "Figure3B_CSI_module_heatmap_v10_GLIA_readability.pdf",
    "Figure3C_regulon_onset_landscape_v5.pdf",
    "Figure3D_early_regulon_dynamics_v5.pdf",
    "Figure3E_formal_Target_program_enrichment_v3_Mod_labels.pdf",
    "Figure3F_top10_ranked_target_landscape_v9_black_borders.pdf"
  ),

  svg_name = c(
    "Figure3A_RSS_rankings_v6.svg",
    "Figure3B_CSI_module_heatmap_v10_GLIA_readability.svg",
    "Figure3C_regulon_onset_landscape_v5.svg",
    "Figure3D_early_regulon_dynamics_v5.svg",
    "Figure3E_formal_Target_program_enrichment_v3_Mod_labels.svg",
    "Figure3F_top10_ranked_target_landscape_v9_black_borders.svg"
  ),

  tiff_name = c(
    "Figure3A_RSS_rankings_v6.tiff",
    "Figure3B_CSI_module_heatmap_v10_GLIA_readability.tiff",
    "Figure3C_regulon_onset_landscape_v5.tiff",
    "Figure3D_early_regulon_dynamics_v5.tiff",
    "Figure3E_formal_Target_program_enrichment_v3_Mod_labels.tiff",
    "Figure3F_top10_ranked_target_landscape_v9_black_borders.tiff"
  ),

  png_name = c(
    "Figure3A_RSS_rankings_v6.png",
    "Figure3B_CSI_module_heatmap_v10_GLIA_readability.png",
    "Figure3C_regulon_onset_landscape_v5.png",
    "Figure3D_early_regulon_dynamics_v5.png",
    "Figure3E_formal_Target_program_enrichment_v3_Mod_labels.png",
    "Figure3F_top10_ranked_target_landscape_v9_black_borders.png"
  ),

  script_pattern = c(
    "^85A_Figure3_PanelA_RSS_rankings_v6_label_spacing_fix.*\\.R$",
    "^85B_Figure3_PanelB_CSI_module_heatmap_v10_GLIA_readability.*\\.R$",
    "^85C_Figure3_PanelC_regulon_onset_landscape_v5_clean_origin_gene_spacing.*\\.R$",
    "^85D_Figure3_PanelD_early_regulon_dynamics_v5_D_top_left_legend_same_row.*\\.R$",
    "^85E_Figure3_PanelE_formal_Target_program_enrichment_v3_Mod_labels.*\\.R$",
    "^85F_Figure3_PanelF_v9_black_borders.*\\.R$"
  ),

  stringsAsFactors = FALSE
)

panel_files$pdf_path  <- file.path(panel_files$output_dir, panel_files$pdf_name)
panel_files$svg_path  <- file.path(panel_files$output_dir, panel_files$svg_name)
panel_files$tiff_path <- file.path(panel_files$output_dir, panel_files$tiff_name)
panel_files$png_path  <- file.path(panel_files$output_dir, panel_files$png_name)

# ------------------------------------------------------------------------------
# 3. VERIFY PANEL BUNDLES
# ------------------------------------------------------------------------------
refresh_status <- function(tbl) {
  tbl$pdf_exists  <- file.exists(tbl$pdf_path)
  tbl$svg_exists  <- file.exists(tbl$svg_path)
  tbl$tiff_exists <- file.exists(tbl$tiff_path)
  tbl$png_exists  <- file.exists(tbl$png_path)
  tbl$all_required_exist <-
    tbl$pdf_exists &
    tbl$svg_exists &
    tbl$tiff_exists &
    tbl$png_exists
  tbl
}

find_panel_script <- function(pattern, panel_id) {
  hits <- sort(list.files(
    SCRIPT_DIR,
    pattern = pattern,
    full.names = TRUE,
    recursive = TRUE
  ))
  if (length(hits) == 0L) {
    stop(
      "Panel ", panel_id,
      " output bundle is incomplete and its source R script was not found in:\n",
      SCRIPT_DIR,
      call. = FALSE
    )
  }
  hits[1]
}

run_panel_script_safely <- function(script_path, panel_id) {
  message("Running missing panel ", panel_id, ": ", basename(script_path))
  panel_env <- new.env(parent = globalenv())
  tryCatch(
    sys.source(script_path, envir = panel_env, chdir = FALSE),
    error = function(e) {
      stop(
        "Panel ", panel_id, " failed while sourcing:\n",
        script_path,
        "\n\nOriginal error:\n",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )
  invisible(TRUE)
}

panel_files <- refresh_status(panel_files)

write.csv(
  panel_files,
  file.path(LOG_DIR, "00_panel_output_status_before.csv"),
  row.names = FALSE
)

if (!all(panel_files$all_required_exist)) {

  if (!RUN_PANEL_SCRIPTS_IF_OUTPUT_MISSING) {
    stop("One or more final panel bundles are incomplete.", call. = FALSE)
  }

  for (i in which(!panel_files$all_required_exist)) {
    script_path <- find_panel_script(
      panel_files$script_pattern[i],
      panel_files$panel[i]
    )
    run_panel_script_safely(script_path, panel_files$panel[i])
  }

  panel_files <- refresh_status(panel_files)
}

write.csv(
  panel_files,
  file.path(LOG_DIR, "01_panel_output_status_after.csv"),
  row.names = FALSE
)

if (!all(panel_files$all_required_exist)) {
  stop("At least one final A-F panel bundle remains incomplete.", call. = FALSE)
}

# ------------------------------------------------------------------------------
# 4. COPY STANDARDIZED SINGLE-PANEL OUTPUTS
# ------------------------------------------------------------------------------
single_manifest <- vector("list", nrow(panel_files))

for (i in seq_len(nrow(panel_files))) {

  panel_id <- paste0("Figure3", panel_files$panel[i])
  panel_dir <- file.path(SINGLE_DIR, panel_id)
  dir.create(panel_dir, recursive = TRUE, showWarnings = FALSE)

  dest_pdf  <- file.path(panel_dir, paste0(panel_id, ".pdf"))
  dest_svg  <- file.path(panel_dir, paste0(panel_id, ".svg"))
  dest_tiff <- file.path(panel_dir, paste0(panel_id, ".tiff"))
  dest_png  <- file.path(panel_dir, paste0(panel_id, ".png"))

  ok <- c(
    file.copy(panel_files$pdf_path[i],  dest_pdf,  overwrite = TRUE),
    file.copy(panel_files$svg_path[i],  dest_svg,  overwrite = TRUE),
    file.copy(panel_files$tiff_path[i], dest_tiff, overwrite = TRUE),
    file.copy(panel_files$png_path[i],  dest_png,  overwrite = TRUE)
  )

  if (!all(ok)) {
    stop("Failed to copy standardized files for ", panel_id, call. = FALSE)
  }

  single_manifest[[i]] <- data.frame(
    panel = panel_files$panel[i],
    pdf = dest_pdf,
    svg = dest_svg,
    tiff = dest_tiff,
    png = dest_png,
    stringsAsFactors = FALSE
  )
}

single_manifest <- do.call(rbind, single_manifest)

write.csv(
  single_manifest,
  file.path(LOG_DIR, "02_single_panel_manifest.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 5. IMAGE HELPERS
# ------------------------------------------------------------------------------
tiff_for <- function(panel_id) {
  panel_files$tiff_path[match(panel_id, panel_files$panel)]
}

erase_original_letter <- function(img, panel_id) {
  mask <- image_blank(
    width = LETTER_MASK[[panel_id]]$w,
    height = LETTER_MASK[[panel_id]]$h,
    color = "white"
  )
  image_composite(
    img,
    mask,
    gravity = "northwest",
    offset = "+0+0"
  )
}

prepare_panel_body <- function(panel_id) {

  img <- image_read(tiff_for(panel_id))
  img <- image_trim(img)

  # Remove only the old letter region with a small panel-specific mask.
  img <- erase_original_letter(img, panel_id)

  # Restore a tiny safety border.
  img <- image_border(
    img,
    color = "white",
    geometry = paste0(PANEL_BORDER_X, "x", PANEL_BORDER_Y)
  )

  # Row-/panel-specific rescaling, but always keep the scaled panel inside
  # the column cell with a small left/right safety margin.
  target_width <- max(
    1L,
    min(
      as.integer(round(CELL_WIDTH_PX * PANEL_SCALE[[panel_id]])),
      CELL_WIDTH_PX - 2L * CELL_INNER_X_PAD
    )
  )

  img <- image_resize(img, geometry = paste0(target_width, "x"))
  img
}

place_on_canvas <- function(canvas, img, x, y) {
  image_composite(
    canvas,
    img,
    gravity = "northwest",
    offset = paste0("+", x, "+", y)
  )
}

annotate_row_letters <- function(row_canvas, left_letter, right_letter) {

  row_canvas <- image_annotate(
    row_canvas,
    text = left_letter,
    gravity = "northwest",
    location = paste0("+", LETTER_X_LEFT, "+", LETTER_Y),
    size = LETTER_SIZE,
    font = LETTER_FONT,
    weight = LETTER_WEIGHT,
    color = LETTER_COLOR
  )

  row_canvas <- image_annotate(
    row_canvas,
    text = right_letter,
    gravity = "northwest",
    location = paste0("+", LETTER_X_RIGHT, "+", LETTER_Y),
    size = LETTER_SIZE,
    font = LETTER_FONT,
    weight = LETTER_WEIGHT,
    color = LETTER_COLOR
  )

  row_canvas
}

build_centered_row <- function(left_panel_id, right_panel_id) {

  left_img  <- prepare_panel_body(left_panel_id)
  right_img <- prepare_panel_body(right_panel_id)

  left_info  <- image_info(left_img)
  right_info <- image_info(right_img)

  # dynamic body height based on the larger of the two scaled panels
  body_height <- max(left_info$height, right_info$height)

  # Reserve a true top strip for the standardized panel letters.
  # This prevents the larger manuscript-level letters from covering real content.
  row_height <- body_height + LETTER_HEADROOM_PX

  # exact vertical center alignment within the body region, below the letter strip
  left_y  <- LETTER_HEADROOM_PX +
    as.integer(round((body_height - left_info$height) / 2))
  right_y <- LETTER_HEADROOM_PX +
    as.integer(round((body_height - right_info$height) / 2))

  # exact horizontal center alignment within each column cell
  left_x_base <- as.integer(round((CELL_WIDTH_PX - left_info$width) / 2))
  left_x <- min(
    left_x_base + LEFT_PANEL_NUDGE_PX,
    CELL_WIDTH_PX - left_info$width
  )

  right_x <- CELL_WIDTH_PX + GAP_X_PX +
    as.integer(round((CELL_WIDTH_PX - right_info$width) / 2))

  row_canvas <- image_blank(
    width = COMBINED_WIDTH_PX,
    height = row_height,
    color = "white"
  )

  row_canvas <- place_on_canvas(row_canvas, left_img,  left_x,  left_y)
  row_canvas <- place_on_canvas(row_canvas, right_img, right_x, right_y)
  row_canvas <- annotate_row_letters(row_canvas, left_panel_id, right_panel_id)

  list(
    image = row_canvas,
    audit = data.frame(
      left_panel = left_panel_id,
      right_panel = right_panel_id,
      row_height_px = row_height,
      left_width_px = left_info$width,
      right_width_px = right_info$width,
      left_height_px = left_info$height,
      right_height_px = right_info$height,
      left_x_offset_px = left_x,
      right_x_offset_px = right_x,
      left_y_offset_px = left_y,
      right_y_offset_px = right_y,
      left_center_y_px = left_y + left_info$height / 2,
      right_center_y_px = right_y + right_info$height / 2,
      stringsAsFactors = FALSE
    )
  )
}

# ------------------------------------------------------------------------------
# 6. BUILD ROWS
# ------------------------------------------------------------------------------
row1_obj <- build_centered_row("A", "B")
row2_obj <- build_centered_row("C", "D")
row3_obj <- build_centered_row("E", "F")

row1 <- row1_obj$image
row2 <- row2_obj$image
row3 <- row3_obj$image

row_audit <- do.call(rbind, list(row1_obj$audit, row2_obj$audit, row3_obj$audit))
row_audit$center_difference_px <- abs(row_audit$left_center_y_px - row_audit$right_center_y_px)

write.csv(
  row_audit,
  file.path(LOG_DIR, "03_center_alignment_audit.csv"),
  row.names = FALSE
)

if (any(row_audit$center_difference_px > 1)) {
  stop("Internal center-alignment audit failed.", call. = FALSE)
}

# ------------------------------------------------------------------------------
# 7. STACK ROWS
# ------------------------------------------------------------------------------
make_gap <- function(height) {
  image_blank(width = COMBINED_WIDTH_PX, height = height, color = "white")
}

Figure3_core <- image_append(
  c(row1, make_gap(GAP_Y_PX), row2, make_gap(GAP_Y_PX), row3),
  stack = TRUE
)

core_info <- image_info(Figure3_core)

Figure3_combined <- image_extent(
  Figure3_core,
  geometry = paste0(COMBINED_WIDTH_PX, "x", core_info$height + BOTTOM_SAFETY_PX),
  gravity = "northwest",
  color = "white"
)

final_info <- image_info(Figure3_combined)

# ------------------------------------------------------------------------------
# 8. EXPORT
# ------------------------------------------------------------------------------
combined_tiff <- file.path(
  COMBINED_DIR,
  "Figure3_GLIA_A_to_F_v13_Ev3_Fv9_600dpi.tiff"
)

combined_png <- file.path(
  COMBINED_DIR,
  "Figure3_GLIA_A_to_F_v13_Ev3_Fv9_600dpi.png"
)

combined_pdf <- file.path(
  COMBINED_DIR,
  "Figure3_GLIA_A_to_F_v13_Ev3_Fv9_review.pdf"
)

image_write(
  Figure3_combined,
  path = combined_tiff,
  format = "tiff",
  density = paste0(DPI),
  compression = "LZW"
)

image_write(
  Figure3_combined,
  path = combined_png,
  format = "png",
  density = paste0(DPI)
)

PDF_WIDTH_IN  <- final_info$width / DPI
PDF_HEIGHT_IN <- final_info$height / DPI

grDevices::cairo_pdf(
  filename = combined_pdf,
  width = PDF_WIDTH_IN,
  height = PDF_HEIGHT_IN,
  family = "Arial",
  onefile = TRUE
)

par(mar = c(0, 0, 0, 0), xaxs = "i", yaxs = "i")
plot.new()
rasterImage(
  as.raster(Figure3_combined),
  xleft = 0,
  ybottom = 0,
  xright = 1,
  ytop = 1,
  interpolate = TRUE
)
dev.off()

# ------------------------------------------------------------------------------
# 9. AUDIT
# ------------------------------------------------------------------------------
combined_manifest <- data.frame(
  format = c("TIFF_600dpi", "PNG_600dpi", "PDF_review"),
  path = c(combined_tiff, combined_png, combined_pdf),
  exists = file.exists(c(combined_tiff, combined_png, combined_pdf)),
  stringsAsFactors = FALSE
)

write.csv(
  combined_manifest,
  file.path(LOG_DIR, "04_combined_manifest.csv"),
  row.names = FALSE
)

letter_audit <- data.frame(
  property = c("font", "size", "weight", "left_column_x", "right_column_x", "letter_y"),
  value = c(LETTER_FONT, LETTER_SIZE, LETTER_WEIGHT, LETTER_X_LEFT, LETTER_X_RIGHT, LETTER_Y),
  stringsAsFactors = FALSE
)

write.csv(
  letter_audit,
  file.path(LOG_DIR, "05_panel_letter_alignment.csv"),
  row.names = FALSE
)

scale_audit <- data.frame(
  panel = names(PANEL_SCALE),
  scale = as.numeric(PANEL_SCALE),
  letter_mask_w = vapply(LETTER_MASK, function(x) x$w, integer(1)),
  letter_mask_h = vapply(LETTER_MASK, function(x) x$h, integer(1)),
  stringsAsFactors = FALSE
)

write.csv(
  scale_audit,
  file.path(LOG_DIR, "06_scale_and_mask_audit.csv"),
  row.names = FALSE
)

geometry_audit <- data.frame(
  metric = c(
    "combined_width_mm",
    "combined_width_px",
    "combined_height_px",
    "cell_width_px",
    "horizontal_gap_px",
    "vertical_gap_px",
    "letter_headroom_px",
    "left_panel_nudge_px",
    "bottom_safety_px",
    "dpi"
  ),
  value = c(
    COMBINED_WIDTH_MM,
    final_info$width,
    final_info$height,
    CELL_WIDTH_PX,
    GAP_X_PX,
    GAP_Y_PX,
    LETTER_HEADROOM_PX,
    LEFT_PANEL_NUDGE_PX,
    BOTTOM_SAFETY_PX,
    DPI
  ),
  stringsAsFactors = FALSE
)

write.csv(
  geometry_audit,
  file.path(LOG_DIR, "07_geometry_audit.csv"),
  row.names = FALSE
)

writeLines(
  capture.output(sessionInfo()),
  file.path(LOG_DIR, "sessionInfo.txt")
)

# ------------------------------------------------------------------------------
# 10. SUMMARY
# ------------------------------------------------------------------------------
cat("\n============================================================\n")
cat("Figure 3 v13 (E v3 + F v9) / Figure-4 visual-match assembly completed successfully.\n\n")

cat("Key fixes:\n")
cat("- Panel E uses Mod1/Mod2 display terminology (v3).\n")
cat("- Panel F uses Mod1/Mod2 and uniform black point borders (v9).\n")
cat("- D uses a much larger mask so the extra original D is fully removed.\n")
cat("- B keeps a larger mask; C and E keep smaller masks to protect nearby labels.\n")
cat("- A/B are slightly enlarged but clipped safely inside the cell width.\n")
cat("- E/F are enlarged relative to v8 for better final-size readability.\n")
cat("- Panel bodies are centered horizontally within each cell.\n")
cat("- Paired panels remain vertically center-aligned.\n\n")

cat("Combined outputs:\n")
print(combined_manifest, row.names = FALSE)

cat("\nReview PDF:\n")
cat(combined_pdf, "\n")

cat("\nPrimary 600-dpi TIFF:\n")
cat(combined_tiff, "\n")
cat("============================================================\n")
