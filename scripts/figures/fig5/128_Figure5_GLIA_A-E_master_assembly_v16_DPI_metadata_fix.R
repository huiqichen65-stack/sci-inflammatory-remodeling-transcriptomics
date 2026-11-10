# ============================================================
# 128_Figure5_GLIA_A-E_master_assembly_v16_DPI_metadata_fix.R
# v16: preserves all v15 analysis/layout logic and explicitly writes 600-dpi metadata
# into exported PNG/TIFF files; analysis/statistics/layout are unchanged.
# Figure 5 | TRUE immutable-PDF assembly
# ------------------------------------------------------------
# FIXES THE v9 RUNTIME ERRORS:
#   - removes magick::image_colorspace() entirely
#     (not exported in the user's magick version)
#   - fixes pdftools::pdf_convert() filename-template warnings
#
# MORE IMPORTANTLY, THIS VERSION IMPLEMENTS THE REQUEST EXACTLY:
#
#   A/B:
#     - remove the ORIGINAL A/B letters in the source plotting scripts
#       by changing only label = "A"/"B" to label = ""
#     - keep the original timepoint-label strip and all layout geometry
#     - add ONE new standardized external letter after final placement
#
#   C/D/E:
#     - remove the original panel letters in temporary script copies
#       while preserving the original layout geometry
#     - add ONE new standardized external letter after final placement
#
#   C/D row:
#     - C is smaller
#     - D is larger
#     - C artwork is vertically CENTERED relative to D
#     - C/D letters themselves remain top-aligned at the row-box top
#
#   D COLOR:
#     - each letterless PDF is rendered DIRECTLY at its final effective
#       pixel size; no 1200->600 downsampling and no colorspace conversion
#     - this avoids the extra resampling step that made D look lighter
#
# INTERNAL PANEL CONTENT:
#   - no ggplot/patchwork re-layout during final assembly
#   - final panel PDFs are treated as finished artworks
#   - only the original panel-letter text is blanked in temporary scripts
#   - whole-panel scaling only
#
# FINAL LAYOUT:
#   A  full width
#   B  full width
#   C | D   (C smaller, D larger; vertically center-aligned artwork)
#   E  full width
# ============================================================

options(stringsAsFactors = FALSE, warn = 1)

# ------------------------------------------------------------
# 0. USER SETTINGS
# ------------------------------------------------------------
ROOT <- Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
SCRIPT_DIR <- ROOT

OUT_DIR <- file.path(
  ROOT,
  "09_final_figures",
  "GLIA",
  "Figure5",
  "FINAL_ASSEMBLY_v16_DPI_METADATA_FIX"
)

dir.create(
  OUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

# Final output resolution
DPI_OUT <- 600

# Final page / content widths
PAGE_WIDTH_MM <- 180
CONTENT_WIDTH_MM <- 176

# External spacing only
TOP_MARGIN_MM <- 3
BOTTOM_MARGIN_MM <- 3
ROW_GAP_MM <- 1.5
CD_GAP_MM <- 2.5

# D larger, C smaller
# 68 + 3 gap + 105 = 176 mm
CD_WIDTH_MM <- c(
  C = 74.5,
  D = 99.0
)

# Standardized external letters
LETTER_FAMILY <- "Arial"
LETTER_SIZE_PT <- 16
LETTER_WEIGHT <- 700

# Letter inset from the upper-left corner of each PANEL BOX
LETTER_INSET_X_MM <- 1.2
LETTER_INSET_Y_MM <- 1.0

# ------------------------------------------------------------
# 1. REQUIRED PACKAGES
# ------------------------------------------------------------
required_pkgs <- c(
  "pdftools",
  "magick",
  "png"
)

missing_pkgs <- required_pkgs[
  !vapply(
    required_pkgs,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_pkgs) > 0) {
  stop(
    "Missing required R package(s): ",
    paste(missing_pkgs, collapse = ", "),
    "\nInstall them first:\n",
    'install.packages(c("pdftools","magick","png"))',
    call. = FALSE
  )
}

# ------------------------------------------------------------
# 2. HELPERS
# ------------------------------------------------------------
stop2 <- function(...) stop(..., call. = FALSE)

first_existing <- function(paths) {
  hit <- paths[file.exists(paths)][1]
  if (length(hit) == 0 || is.na(hit)) return(NA_character_)
  normalizePath(
    hit,
    winslash = "/",
    mustWork = TRUE
  )
}

resolve_script <- function(candidates) {

  # Prefer exact filenames in xjz root, then search xjz recursively.
  direct <- first_existing(file.path(SCRIPT_DIR, candidates))
  if (!is.na(direct)) return(direct)

  all_r <- list.files(
    SCRIPT_DIR,
    pattern = "\\.R$",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )

  for (cand in candidates) {
    hit <- all_r[basename(all_r) == cand]
    if (length(hit) > 0) {
      return(normalizePath(hit[1], winslash = "/", mustWork = TRUE))
    }
  }

  stop2(
    "Could not locate required Figure 5 script anywhere under xjz.\n",
    "Required filename(s):\n  ",
    paste(candidates, collapse = "\n  ")
  )
}

mm_to_px <- function(mm) {
  round(
    mm / 25.4 * DPI_OUT
  )
}

# ------------------------------------------------------------
# 3. LOCATE FROZEN PANEL SCRIPTS
# ------------------------------------------------------------
SCRIPT_A <- resolve_script(
  c(
    "101_Figure5A_GLIA_7dpi_target_endothelial_spatial_coupling_reference_style_v4_rho_fix.R",
    "101_Figure5A_GLIA_7dpi_target_endothelial_spatial_coupling_reference_style_v3_final.R",
    "101_Figure5A_GLIA_7dpi_target_endothelial_spatial_coupling_reference_style_v3_final(1).R"
  )
)

SCRIPT_B <- resolve_script(
  c(
    "102_Figure5B_GLIA_2months_target_endothelial_spatial_coupling_reference_style_v2_rho_fix.R",
    "102_Figure5B_GLIA_2months_target_endothelial_spatial_coupling_reference_style.R",
    "102_Figure5B_GLIA_2months_target_endothelial_spatial_coupling_reference_style(1).R"
  )
)

SCRIPT_C <- resolve_script(
  c(
    "108_Figure5C_GLIA_target_endothelial_distance_gradient_publication_v5_ytitle_final.R"
  )
)

SCRIPT_D <- resolve_script(
  c(
    "117_Figure5D_GLIA_endothelial_response_spatial_maps_v5_rowlabel_fix.R"
  )
)

SCRIPT_E <- resolve_script(
  c(
    "118_Figure5E_GLIA_target_endothelial_response_spatial_coupling_v5_rho_fix.R",
    "118_Figure5E_GLIA_target_endothelial_response_spatial_coupling_v4_final.R",
    "118_Figure5E_GLIA_target_endothelial_response_spatial_coupling_v4_final(1).R"
  )
)

# ------------------------------------------------------------
# 4. CREATE LETTERLESS COPIES OF A-E SCRIPTS
# ------------------------------------------------------------
# The replacement patterns below remove ONLY the original figure-panel
# letter text. They do NOT alter the layout object that contains the letter.
# Therefore all original blank space / alignment geometry is preserved.

LETTERLESS_SUBDIR <- "09_final_figures/GLIA/Figure5/FINAL_ASSEMBLY_v16_DPI_METADATA_FIX/letterless_source"

panel_specs <- list(
  A = list(
    src = SCRIPT_A,
    remove_pattern = 'label\\s*=\\s*"A"',
    remove_replacement = 'label = ""',
    old_basename = "Figure5A_7dpi_target_endothelial_spatial_coupling_reference_style_v3_final",
    new_basename = "Figure5A_LETTERLESS"
  ),

  B = list(
    src = SCRIPT_B,
    remove_pattern = 'label\\s*=\\s*"B"',
    remove_replacement = 'label = ""',
    old_basename = "Figure5B_2months_target_endothelial_spatial_coupling_reference_style",
    new_basename = "Figure5B_LETTERLESS"
  ),

  C = list(
    src = SCRIPT_C,
    remove_pattern = 'tag\\s*=\\s*NULL',
    remove_replacement = 'tag = NULL',
    old_basename = "Figure5C_target_endothelial_distance_to_core_gradient_v5_ytitle_final",
    new_basename = "Figure5C_LETTERLESS"
  ),

  D = list(
    src = SCRIPT_D,
    remove_pattern = 'make_blank_header\\("D"\\)',
    remove_replacement = 'make_blank_header("")',
    old_basename = "Figure5D_endothelial_response_spatial_maps_v5_rowlabel_fix",
    new_basename = "Figure5D_LETTERLESS"
  ),

  E = list(
    src = SCRIPT_E,
    remove_pattern = 'make_blank_header\\("E"\\)',
    remove_replacement = 'make_blank_header("")',
    old_basename = "Figure5E_target_endothelial_response_spatial_coupling_v4_final",
    new_basename = "Figure5E_LETTERLESS"
  )
)

make_letterless_pdf <- function(
  panel_letter,
  spec
) {

  x <- readLines(
    spec$src,
    warn = FALSE,
    encoding = "UTF-8"
  )

  n_remove <- sum(
    grepl(
      spec$remove_pattern,
      x,
      perl = TRUE
    )
  )

  if (n_remove != 1L) {
    stop2(
      "Panel ",
      panel_letter,
      ": expected exactly one original-letter target, found ",
      n_remove,
      ". No temporary script was run."
    )
  }

  # Remove ONLY original panel-letter text.
  x <- sub(
    spec$remove_pattern,
    spec$remove_replacement,
    x,
    perl = TRUE
  )

  # Redirect output folder.
  x <- sub(
    'OUTPUT_DIR_NAME\\s*<-\\s*"[^"]+"',
    paste0(
      'OUTPUT_DIR_NAME <- "',
      LETTERLESS_SUBDIR,
      '"'
    ),
    x,
    perl = TRUE
  )

  # Redirect output basename.
  old_quoted <- paste0(
    '"',
    spec$old_basename,
    '"'
  )

  new_quoted <- paste0(
    '"',
    spec$new_basename,
    '"'
  )

  n_base <- sum(
    grepl(
      old_quoted,
      x,
      fixed = TRUE
    )
  )

  if (n_base < 1L) {
    stop2(
      "Panel ",
      panel_letter,
      ": could not find frozen output basename in source script."
    )
  }

  x <- gsub(
    old_quoted,
    new_quoted,
    x,
    fixed = TRUE
  )

  temp_script <- file.path(
    OUT_DIR,
    paste0(
      "TEMP_Figure5",
      panel_letter,
      "_LETTERLESS.R"
    )
  )

  writeLines(
    x,
    temp_script,
    useBytes = TRUE
  )

  e <- new.env(
    parent = globalenv()
  )

  cat(
    "Re-running Figure 5",
    panel_letter,
    " with ONLY the original panel letter blanked ...\n",
    sep = ""
  )

  sys.source(
    temp_script,
    envir = e,
    keep.source = TRUE
  )

  expected_pdf <- file.path(
    ROOT,
    LETTERLESS_SUBDIR,
    paste0(
      spec$new_basename,
      ".pdf"
    )
  )

  if (!file.exists(expected_pdf)) {
    stop2(
      "Panel ",
      panel_letter,
      ": expected letterless PDF was not created:\n  ",
      expected_pdf
    )
  }

  normalizePath(
    expected_pdf,
    winslash = "/",
    mustWork = TRUE
  )
}

LETTERLESS_PDF <- setNames(
  character(5),
  c(
    "A",
    "B",
    "C",
    "D",
    "E"
  )
)

for (lab in names(panel_specs)) {
  LETTERLESS_PDF[lab] <- make_letterless_pdf(
    lab,
    panel_specs[[lab]]
  )
}

cat(
  "\nLetterless frozen PDFs:\n",
  paste0(
    names(LETTERLESS_PDF),
    ": ",
    LETTERLESS_PDF,
    collapse = "\n"
  ),
  "\n\n",
  sep = ""
)

# ------------------------------------------------------------
# 5. READ EXACT PDF PAGE RATIOS
# ------------------------------------------------------------
page_info <- lapply(
  LETTERLESS_PDF,
  pdftools::pdf_pagesize
)

aspect <- vapply(
  page_info,
  function(z) {
    z$width[1] / z$height[1]
  },
  FUN.VALUE = numeric(1)
)

source_width_in <- vapply(
  page_info,
  function(z) {
    z$width[1] / 72
  },
  FUN.VALUE = numeric(1)
)

# ------------------------------------------------------------
# 6. FINAL PHYSICAL PANEL SIZES
# ------------------------------------------------------------
# A/B/E full width.
# C/D fixed external widths, with D larger.

panel_width_mm <- c(
  A = CONTENT_WIDTH_MM,
  B = CONTENT_WIDTH_MM,
  C = unname(CD_WIDTH_MM["C"]),
  D = unname(CD_WIDTH_MM["D"]),
  E = CONTENT_WIDTH_MM
)

panel_height_mm <- panel_width_mm / aspect

# C/D artwork vertically centered inside a common row box.
CD_ROW_HEIGHT_MM <- max(
  panel_height_mm["C"],
  panel_height_mm["D"]
)

PAGE_HEIGHT_MM <-
  TOP_MARGIN_MM +
  panel_height_mm["A"] +
  ROW_GAP_MM +
  panel_height_mm["B"] +
  ROW_GAP_MM +
  CD_ROW_HEIGHT_MM +
  ROW_GAP_MM +
  panel_height_mm["E"] +
  BOTTOM_MARGIN_MM

cat(
  "Final panel sizes:\n",
  paste0(
    names(panel_width_mm),
    ": ",
    sprintf(
      "%.2f x %.2f mm",
      panel_width_mm,
      panel_height_mm
    ),
    collapse = "\n"
  ),
  "\n",
  sprintf(
    "Complete page: %.2f x %.2f mm\n\n",
    PAGE_WIDTH_MM,
    PAGE_HEIGHT_MM
  ),
  sep = ""
)

# ------------------------------------------------------------
# 7. RENDER EACH PDF DIRECTLY AT FINAL EFFECTIVE SIZE
# ------------------------------------------------------------
# Instead of:
#     render huge -> resize down
# we calculate the PDF render DPI needed to produce approximately the
# final 600-dpi pixel width directly.
#
# This removes an extra interpolation step and is the safest raster
# route for keeping D visually close to its frozen PDF.

RENDER_DIR <- file.path(
  OUT_DIR,
  "_direct_final_render"
)

dir.create(
  RENDER_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

rendered_files <- setNames(
  character(5),
  names(LETTERLESS_PDF)
)

for (lab in names(LETTERLESS_PDF)) {

  desired_width_px <- mm_to_px(
    panel_width_mm[lab]
  )

  render_dpi <- round(
    desired_width_px /
      source_width_in[lab]
  )

  render_dpi <- max(
    72,
    render_dpi
  )

  # pdftools::pdf_convert internally calls sprintf(filenames, pages, format).
  # Therefore the filename MUST contain placeholders.
  filename_template <- file.path(
    RENDER_DIR,
    paste0(
      "panel_",
      lab,
      "_page_%d.%s"
    )
  )

  converted <- pdftools::pdf_convert(
    pdf = LETTERLESS_PDF[lab],
    format = "png",
    pages = 1,
    filenames = filename_template,
    dpi = render_dpi,
    verbose = FALSE
  )

  if (length(converted) < 1L || !file.exists(converted[1])) {
    stop2(
      "Direct rendering failed for panel ",
      lab
    )
  }

  rendered_files[lab] <- normalizePath(
    converted[1],
    winslash = "/",
    mustWork = TRUE
  )

  cat(
    "Panel ",
    lab,
    ": direct render DPI = ",
    render_dpi,
    "\n",
    sep = ""
  )
}

# ------------------------------------------------------------
# 8. READ DIRECTLY RENDERED PANEL ARTWORKS
# ------------------------------------------------------------
img <- lapply(
  rendered_files,
  magick::image_read
)

names(img) <- names(rendered_files)

# Tiny rounding differences are permitted because render DPI is integer.
# No image_resize() is used: this is deliberate.

img_info <- lapply(
  img,
  magick::image_info
)

# ------------------------------------------------------------
# 9. CREATE FINAL PAGE
# ------------------------------------------------------------
canvas_w_px <- mm_to_px(
  PAGE_WIDTH_MM
)

canvas_h_px <- mm_to_px(
  PAGE_HEIGHT_MM
)

content_w_px <- mm_to_px(
  CONTENT_WIDTH_MM
)

left_full_px <- round(
  (
    canvas_w_px -
      content_w_px
  ) / 2
)

row_gap_px <- mm_to_px(
  ROW_GAP_MM
)

cd_gap_px <- mm_to_px(
  CD_GAP_MM
)

canvas <- magick::image_blank(
  width = canvas_w_px,
  height = canvas_h_px,
  color = "white"
)

# ------------------------------------------------------------
# 10. PLACE A
# ------------------------------------------------------------
x_A <- left_full_px
y_A <- mm_to_px(
  TOP_MARGIN_MM
)

canvas <- magick::image_composite(
  canvas,
  img$A,
  operator = "over",
  offset = paste0(
    "+",
    x_A,
    "+",
    y_A
  )
)

# ------------------------------------------------------------
# 11. PLACE B
# ------------------------------------------------------------
x_B <- left_full_px
y_B <- y_A +
  img_info$A$height[1] +
  row_gap_px

canvas <- magick::image_composite(
  canvas,
  img$B,
  operator = "over",
  offset = paste0(
    "+",
    x_B,
    "+",
    y_B
  )
)

# ------------------------------------------------------------
# 12. PLACE C / D
# ------------------------------------------------------------
y_CD_top <- y_B +
  img_info$B$height[1] +
  row_gap_px

# Use actual rendered widths.
cd_total_width_px <-
  img_info$C$width[1] +
  cd_gap_px +
  img_info$D$width[1]

x_C <- round(
  (
    canvas_w_px -
      cd_total_width_px
  ) / 2
)

x_D <- x_C +
  img_info$C$width[1] +
  cd_gap_px

row_cd_height_px <- max(
  img_info$C$height[1],
  img_info$D$height[1]
)

# USER REQUEST:
# C artwork goes DOWN and shares the same vertical CENTER as D.
y_C <- y_CD_top +
  round(
    (
      row_cd_height_px -
        img_info$C$height[1]
    ) / 2
  )

y_D <- y_CD_top +
  round(
    (
      row_cd_height_px -
        img_info$D$height[1]
    ) / 2
  )

canvas <- magick::image_composite(
  canvas,
  img$C,
  operator = "over",
  offset = paste0(
    "+",
    x_C,
    "+",
    y_C
  )
)

canvas <- magick::image_composite(
  canvas,
  img$D,
  operator = "over",
  offset = paste0(
    "+",
    x_D,
    "+",
    y_D
  )
)

# ------------------------------------------------------------
# 13. PLACE E
# ------------------------------------------------------------
x_E <- left_full_px

y_E <- y_CD_top +
  row_cd_height_px +
  row_gap_px

canvas <- magick::image_composite(
  canvas,
  img$E,
  operator = "over",
  offset = paste0(
    "+",
    x_E,
    "+",
    y_E
  )
)

# ------------------------------------------------------------
# 14. ADD ONE STANDARDIZED LETTER PER PANEL
# ------------------------------------------------------------
# All letters are drawn AFTER final panel scaling and placement.
# Therefore they are physically the same 16 pt in the final figure.
#
# C/D letters use y_CD_top, so their letters are TOP-ALIGNED,
# even though the smaller C artwork itself is vertically centered lower.

pointsize_px <- round(
  LETTER_SIZE_PT *
    DPI_OUT /
    72
)

inset_x_px <- mm_to_px(
  LETTER_INSET_X_MM
)

inset_y_px <- mm_to_px(
  LETTER_INSET_Y_MM
)

letter_xy <- list(
  A = c(
    x_A + inset_x_px,
    y_A + inset_y_px
  ),

  B = c(
    x_B + inset_x_px,
    y_B + inset_y_px
  ),

  C = c(
    x_C + inset_x_px,
    y_CD_top + inset_y_px
  ),

  D = c(
    x_D + inset_x_px,
    y_CD_top + inset_y_px
  ),

  E = c(
    x_E + inset_x_px,
    y_E + inset_y_px
  )
)

# image_annotate supports gravity/location relative to the WHOLE canvas.
# To guarantee exact coordinates and full visibility, annotate one temporary
# transparent full-size layer for each letter and composite it.

add_letter_to_canvas <- function(
  base,
  letter,
  xy
) {

  overlay <- magick::image_blank(
    width = canvas_w_px,
    height = canvas_h_px,
    color = "transparent"
  )

  overlay <- magick::image_annotate(
    overlay,
    text = letter,
    size = pointsize_px,
    gravity = "northwest",
    location = paste0(
      "+",
      round(xy[1]),
      "+",
      round(xy[2])
    ),
    font = LETTER_FAMILY,
    weight = LETTER_WEIGHT,
    color = "black"
  )

  magick::image_composite(
    base,
    overlay,
    operator = "over"
  )
}

for (lab in names(letter_xy)) {
  canvas <- add_letter_to_canvas(
    canvas,
    lab,
    letter_xy[[lab]]
  )
}

# ------------------------------------------------------------
# 15. EXPORT COMPLETE FIGURE
# ------------------------------------------------------------
FULL_PNG <- file.path(
  OUT_DIR,
  "Figure5_GLIA_complete_A-E_v16_DPI_METADATA_FIX_600dpi.png"
)

FULL_TIF <- file.path(
  OUT_DIR,
  "Figure5_GLIA_complete_A-E_v16_DPI_METADATA_FIX_600dpi.tif"
)

FULL_PDF <- file.path(
  OUT_DIR,
  "Figure5_GLIA_complete_A-E_v16_DPI_METADATA_FIX.pdf"
)

magick::image_write(
  canvas,
  path = FULL_PNG,
  format = "png",
  density = paste0(DPI_OUT, "x", DPI_OUT)
)

magick::image_write(
  canvas,
  path = FULL_TIF,
  format = "tiff",
  density = paste0(DPI_OUT, "x", DPI_OUT),
  compression = "lzw"
)

arr <- png::readPNG(
  FULL_PNG
)

grDevices::pdf(
  FULL_PDF,
  width = PAGE_WIDTH_MM / 25.4,
  height = PAGE_HEIGHT_MM / 25.4,
  useDingbats = FALSE,
  bg = "white"
)

grid::grid.newpage()

grid::grid.raster(
  arr,
  width = unit(1, "npc"),
  height = unit(1, "npc"),
  interpolate = FALSE
)

grDevices::dev.off()

# ------------------------------------------------------------
# 16. EXPORT STANDARDIZED SINGLE PANELS
# ------------------------------------------------------------
# These are constructed from the same letterless final artwork.
# No second letter can occur.

single_dir <- file.path(
  OUT_DIR,
  "single_panels"
)

dir.create(
  single_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

make_single_panel <- function(
  lab,
  panel_img
) {

  info <- magick::image_info(
    panel_img
  )

  single_canvas <- magick::image_blank(
    width = info$width[1],
    height = info$height[1],
    color = "white"
  )

  single_canvas <- magick::image_composite(
    single_canvas,
    panel_img,
    operator = "over",
    offset = "+0+0"
  )

  overlay <- magick::image_blank(
    width = info$width[1],
    height = info$height[1],
    color = "transparent"
  )

  overlay <- magick::image_annotate(
    overlay,
    text = lab,
    size = pointsize_px,
    gravity = "northwest",
    location = paste0(
      "+",
      inset_x_px,
      "+",
      inset_y_px
    ),
    font = LETTER_FAMILY,
    weight = LETTER_WEIGHT,
    color = "black"
  )

  single_canvas <- magick::image_composite(
    single_canvas,
    overlay,
    operator = "over"
  )

  png_file <- file.path(
    single_dir,
    paste0(
      "Figure5",
      lab,
      "_standardized_v16.png"
    )
  )

  tif_file <- file.path(
    single_dir,
    paste0(
      "Figure5",
      lab,
      "_standardized_v16.tif"
    )
  )

  magick::image_write(
    single_canvas,
    path = png_file,
    format = "png",
    density = paste0(DPI_OUT, "x", DPI_OUT)
  )

  magick::image_write(
    single_canvas,
    path = tif_file,
    format = "tiff",
    density = paste0(DPI_OUT, "x", DPI_OUT),
    compression = "lzw"
  )
}

for (lab in names(img)) {
  make_single_panel(
    lab,
    img[[lab]]
  )
}

# ------------------------------------------------------------
# 17. AUDIT
# ------------------------------------------------------------
audit <- data.frame(
  panel = c(
    "A",
    "B",
    "C",
    "D",
    "E"
  ),

  original_letter_removed_before_assembly = TRUE,

  replacement_letter_count = 1,

  replacement_letter_family = LETTER_FAMILY,

  replacement_letter_size_pt = LETTER_SIZE_PT,

  replacement_letter_bold = TRUE,

  whole_panel_uniform_scaling = TRUE,

  internal_panel_reflow_during_assembly = FALSE,

  final_width_mm = unname(
    panel_width_mm
  ),

  final_nominal_height_mm = unname(
    panel_height_mm
  ),

  stringsAsFactors = FALSE
)

write.csv(
  audit,
  file.path(
    OUT_DIR,
    "Figure5_GLIA_complete_A-E_v16_DPI_METADATA_FIX_audit.csv"
  ),
  row.names = FALSE
)

cat(
  "\n============================================================\n",
  "Figure 5 v16 DPI-metadata-fix assembly completed successfully.\n",
  "============================================================\n",
  "Runtime fixes:\n",
  "  - no magick::image_colorspace()\n",
  "  - pdf_convert uses valid %d / %s filename templates\n\n",
  "Requested visual fixes:\n",
  "  - A/B original PDF letters removed before rendering\n",
  "  - only ONE new standardized A/B letter remains\n",
  "  - C original letter removed and ONE complete C is added\n",
  "  - C artwork is moved down to center-align with D\n",
  "  - C/D letters remain top-aligned\n",
  "  - D is larger and C is smaller\n",
  "  - D is directly rendered at final scale (no extra downsampling)\n",
  "  - PNG/TIFF outputs explicitly store ", DPI_OUT, " x ", DPI_OUT, " dpi metadata\n\n",
  "Output directory:\n  ",
  OUT_DIR,
  "\n\n",
  "Main PDF:\n  ",
  FULL_PDF,
  "\n",
  "============================================================\n",
  sep = ""
)
