# ==============================================================================
# Figure 4 GLIA final assembly v16 - revised Figure 4E (no pooled cell-level P) + BH-corrected Figure 4F
#
# Based on the supplied v14_ROOT_EF_fix assembly script.
# A-D: run the previously reviewed panel scripts.
# E: preferentially run the revised v9 panel script with pooled cell-level P values removed from display.
# F: reuse the already-reviewed v8 TNF/NF-kB/BH PDF without recalculating it.
# No silent fallback to an old Figure 4F file or source script.
# Layout, rasterization and panel placement are unchanged; the E source-letter mask is tightened to protect the longer annotation.
# All assembly outputs and cleaned panels use v16 names, preserving earlier versions.
# ==============================================================================

# Do not clear an interactive RStudio workspace.
gc()
options(stringsAsFactors = FALSE, scipen = 999)

# ------------------------------------------------------------------------------
# 0. SETTINGS
# ------------------------------------------------------------------------------
ROOT <- Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
OUTDIR <- file.path(
  ROOT,
  "09_final_figures",
  "GLIA",
  "Figure4",
  "final_assembly"
)
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# The standalone Figure4F v8 script must be run first; the assembler does not rerun it.
F_V8_PDF <- file.path(
  ROOT, "09_final_figures", "GLIA", "Figure4", "single_panels", "Figure4F",
  "Figure4F_GSE234774_library_level_target_TNF_NFKB_correlation_v8_BH.pdf"
)
if (!file.exists(F_V8_PDF)) {
  stop(
    "The revised Figure4F v8 PDF was not found:\n", F_V8_PDF,
    "\nRun 86F_Figure4_PanelF_ONE_FILE_v8_TNF_NFKB_BH.R first.",
    call. = FALSE
  )
}


# Final two-column width
FIG_WIDTH_MM  <- 180
FIG_HEIGHT_MM <- 230

FIG_WIDTH_IN  <- FIG_WIDTH_MM / 25.4
FIG_HEIGHT_IN <- FIG_HEIGHT_MM / 25.4

FONT_FAMILY <- "Arial"
LETTER_FONTFACE <- "bold"
LETTER_FONTSIZE_FINAL <- 16
LETTER_FONTSIZE_SINGLE <- 22

# Standardized new letter location
LETTER_X <- 0.012
LETTER_Y <- 0.992

# Raster erase rectangles are specified in SOURCE-IMAGE normalized coordinates:
# x0/x1 from LEFT, y0/y1 from TOP.
#
# These are deliberately panel-specific.
# B is erased much deeper because its old B sits below the top-left corner.
# F is also erased deeper to remove the small remaining mark.
ERASE_RECTS <- list(
  A = list(c(x0 = 0.000, x1 = 0.135, y0 = 0.000, y1 = 0.205)),
  B = list(
    # v8: tight erase box around the source panel letter only.
    # The previous x1=0.155/y1=0.360 box cut into the leading text of B row labels.
    c(x0 = 0.000, x1 = 0.060, y0 = 0.000, y1 = 0.120)
  ),
  C = list(
    # v12: tighten the C erase box; the old wide/deep mask clipped the leading text of the top ligand-receptor pairs.
    c(x0 = 0.000, x1 = 0.055, y0 = 0.000, y1 = 0.115)
  ),
  D = list(c(x0 = 0.000, x1 = 0.045, y0 = 0.000, y1 = 0.085)),
  E = list(
    # v16: tighter mask around the source E only.
    # Figure4E v9 uses a longer "Spearman rho = ..." annotation, so keep the mask
    # shallow/narrow to avoid clipping the annotation while still removing the source letter.
    c(x0 = 0.000, x1 = 0.060, y0 = 0.000, y1 = 0.115)
  ),
  F = list(
    # v9: tight mask around the source F only; the old deep mask clipped the correlation label and y-axis title.
    c(x0 = 0.000, x1 = 0.090, y0 = 0.000, y1 = 0.180)
  )
)

# After removing old letters, trim white margins.
TRIM_THRESHOLD <- 0.985
TRIM_PAD_FRAC <- 0.018

# Panel-specific extra padding AFTER trimming.
# D gets extra TOP whitespace to protect the top violin row.
PAD_AFTER_TRIM <- list(
  A = c(top = 0.015, right = 0.010, bottom = 0.015, left = 0.010),
  B = c(top = 0.020, right = 0.008, bottom = 0.015, left = 0.008),
  C = c(top = 0.015, right = 0.010, bottom = 0.015, left = 0.010),
  D = c(top = 0.105, right = 0.016, bottom = 0.030, left = 0.040),
  E = c(top = 0.055, right = 0.010, bottom = 0.020, left = 0.010),
  F = c(top = 0.085, right = 0.010, bottom = 0.020, left = 0.010)
)

FINAL_PDF <- file.path(OUTDIR, "Figure4_GLIA_final_v16_E_noP_TNF_NFKB_BH.pdf")
FINAL_PNG <- file.path(OUTDIR, "Figure4_GLIA_final_v16_E_noP_TNF_NFKB_BH.png")
FINAL_TIF <- file.path(OUTDIR, "Figure4_GLIA_final_v16_E_noP_TNF_NFKB_BH.tif")
AUDIT_CSV <- file.path(OUTDIR, "Figure4_GLIA_final_v16_E_noP_TNF_NFKB_BH_audit.csv")

if (.Platform$OS.type == "windows") {
  try(
    grDevices::windowsFonts(
      Arial = grDevices::windowsFont("Arial")
    ),
    silent = TRUE
  )
}

# ------------------------------------------------------------------------------
# 1. PACKAGES
# ------------------------------------------------------------------------------
required_pkgs <- c("pdftools", "png", "grid", "abind")
missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_pkgs) > 0) {
  stop(
    "Missing package(s): ",
    paste(missing_pkgs, collapse = ", "),
    "\nInstall with:\ninstall.packages(c(",
    paste(sprintf('"%s"', missing_pkgs), collapse = ", "),
    "))",
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(pdftools)
  library(png)
  library(grid)
  library(abind)
})

# ------------------------------------------------------------------------------
# 2. GENERIC HELPERS
# ------------------------------------------------------------------------------
get_script_dir <- function() {
  cmd <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd, value = TRUE)

  if (length(file_arg) > 0) {
    return(
      dirname(
        normalizePath(
          sub("^--file=", "", file_arg[1])
        )
      )
    )
  }

  fr <- sys.frames()

  for (i in rev(seq_along(fr))) {
    if (!is.null(fr[[i]]$ofile)) {
      return(
        dirname(
          normalizePath(fr[[i]]$ofile)
        )
      )
    }
  }

  getwd()
}

find_first_existing <- function(paths) {
  hit <- paths[file.exists(paths)][1]

  if (length(hit) == 0 || is.na(hit)) {
    return(NA_character_)
  }

  normalizePath(
    hit,
    winslash = "/",
    mustWork = TRUE
  )
}

source_panel_script <- function(script_path) {

  env <- new.env(
    parent = globalenv()
  )

  sys.source(
    script_path,
    envir = env,
    keep.source = FALSE
  )

  pdf_path <- NA_character_
  outdir <- NA_character_

  if (
    exists(
      "PDF_FILE",
      envir = env,
      inherits = FALSE
    )
  ) {
    pdf_path <- get(
      "PDF_FILE",
      envir = env
    )
  }

  if (
    exists(
      "OUTDIR",
      envir = env,
      inherits = FALSE
    )
  ) {
    outdir <- get(
      "OUTDIR",
      envir = env
    )
  }

  if (
    !is.na(pdf_path) &&
    file.exists(pdf_path)
  ) {
    return(
      normalizePath(
        pdf_path,
        winslash = "/",
        mustWork = TRUE
      )
    )
  }

  if (
    !is.na(outdir) &&
    dir.exists(outdir)
  ) {

    pdfs <- list.files(
      outdir,
      pattern = "\\.pdf$",
      full.names = TRUE
    )

    pdfs <- pdfs[
      file.info(pdfs)$isdir %in% FALSE
    ]

    if (length(pdfs) > 0) {

      pdfs <- pdfs[
        order(
          file.info(pdfs)$mtime,
          decreasing = TRUE
        )
      ]

      return(
        normalizePath(
          pdfs[1],
          winslash = "/",
          mustWork = TRUE
        )
      )
    }
  }

  stop(
    "Could not recover PDF output after sourcing:\n",
    script_path,
    call. = FALSE
  )
}

pdf_first_page_to_image <- function(
  pdf_path,
  dpi = 450
) {

  tmp_png <- tempfile(
    fileext = ".png"
  )

  on.exit(
    unlink(tmp_png),
    add = TRUE
  )

  pdftools::pdf_convert(
    pdf_path,
    format = "png",
    pages = 1,
    dpi = dpi,
    filenames = tmp_png,
    verbose = FALSE
  )

  png::readPNG(tmp_png)
}

# ------------------------------------------------------------------------------
# 3. IMAGE CLEANUP
# ------------------------------------------------------------------------------
erase_rect_from_top <- function(
  img,
  rect
) {

  h <- dim(img)[1]
  w <- dim(img)[2]

  x0 <- max(
    1,
    floor(rect["x0"] * w) + 1
  )

  x1 <- min(
    w,
    ceiling(rect["x1"] * w)
  )

  y0 <- max(
    1,
    floor(rect["y0"] * h) + 1
  )

  y1 <- min(
    h,
    ceiling(rect["y1"] * h)
  )

  if (
    x1 >= x0 &&
    y1 >= y0
  ) {

    if (length(dim(img)) == 3) {

      # RGB channels
      img[
        y0:y1,
        x0:x1,
        1:min(3, dim(img)[3])
      ] <- 1

      # Keep alpha opaque if alpha exists
      if (dim(img)[3] >= 4) {
        img[
          y0:y1,
          x0:x1,
          4
        ] <- 1
      }

    } else {

      img[
        y0:y1,
        x0:x1
      ] <- 1
    }
  }

  img
}

erase_old_letter <- function(
  img,
  letter
) {

  rects <- ERASE_RECTS[[letter]]

  if (is.null(rects)) {
    return(img)
  }

  for (r in rects) {
    img <- erase_rect_from_top(
      img,
      r
    )
  }

  img
}

trim_white_image <- function(
  img,
  threshold = TRIM_THRESHOLD,
  pad_frac = TRIM_PAD_FRAC
) {

  h <- dim(img)[1]
  w <- dim(img)[2]

  if (length(dim(img)) == 3) {

    rgb_n <- min(
      3,
      dim(img)[3]
    )

    # pixel is non-white if at least one RGB channel is darker than threshold
    nonwhite <- apply(
      img[, , seq_len(rgb_n), drop = FALSE],
      c(1, 2),
      function(z) {
        any(z < threshold)
      }
    )

  } else {

    nonwhite <- img < threshold
  }

  rows <- which(
    rowSums(nonwhite) > 0
  )

  cols <- which(
    colSums(nonwhite) > 0
  )

  if (
    length(rows) == 0 ||
    length(cols) == 0
  ) {
    return(img)
  }

  pad_y <- ceiling(
    h * pad_frac
  )

  pad_x <- ceiling(
    w * pad_frac
  )

  r0 <- max(
    1,
    min(rows) - pad_y
  )

  r1 <- min(
    h,
    max(rows) + pad_y
  )

  c0 <- max(
    1,
    min(cols) - pad_x
  )

  c1 <- min(
    w,
    max(cols) + pad_x
  )

  img[
    r0:r1,
    c0:c1,
    ,
    drop = FALSE
  ]
}

add_white_padding <- function(
  img,
  top = 0,
  right = 0,
  bottom = 0,
  left = 0
) {

  h <- dim(img)[1]
  w <- dim(img)[2]
  nc <- dim(img)[3]

  pt <- round(h * top)
  pb <- round(h * bottom)
  pl <- round(w * left)
  pr <- round(w * right)

  out <- array(
    1,
    dim = c(
      h + pt + pb,
      w + pl + pr,
      nc
    )
  )

  out[
    (pt + 1):(pt + h),
    (pl + 1):(pl + w),
    seq_len(nc)
  ] <- img

  out
}


resize_width_array <- function(
  img,
  new_w
) {

  h <- dim(img)[1]
  w <- dim(img)[2]
  nc <- dim(img)[3]

  if (new_w <= 1 || w <= 1) {
    return(img)
  }

  x_old <- seq(0, 1, length.out = w)
  x_new <- seq(0, 1, length.out = new_w)

  out <- array(
    1,
    dim = c(h, new_w, nc)
  )

  for (ch in seq_len(nc)) {
    for (r in seq_len(h)) {
      out[r, , ch] <- approx(
        x = x_old,
        y = img[r, , ch],
        xout = x_new,
        method = "linear",
        ties = "ordered"
      )$y
    }
  }

  out
}

compress_B_legend_expand_body <- function(
  img,
  split_frac = 0.905,
  new_legend_frac = 0.050,
  gap_frac = 0.012
) {

  h <- dim(img)[1]
  w <- dim(img)[2]
  nc <- dim(img)[3]

  main_end <- max(
    2,
    min(
      w - 2,
      round(w * split_frac)
    )
  )

  legend_start <- main_end + 1

  main_img <- img[
    ,
    seq_len(main_end),
    ,
    drop = FALSE
  ]

  legend_img <- img[
    ,
    legend_start:w,
    ,
    drop = FALSE
  ]

  gap_w <- max(
    2,
    round(w * gap_frac)
  )

  legend_w <- max(
    8,
    round(w * new_legend_frac)
  )

  main_w <- w - gap_w - legend_w

  main_resized <- resize_width_array(
    main_img,
    main_w
  )

  legend_resized <- resize_width_array(
    legend_img,
    legend_w
  )

  gap_img <- array(
    1,
    dim = c(h, gap_w, nc)
  )

  out <- abind::abind(
    main_resized,
    gap_img,
    legend_resized,
    along = 2
  )

  out
}

prepare_panel_image <- function(
  img,
  letter
) {

  # 1. erase source letter before any scaling
  img <- erase_old_letter(
    img,
    letter
  )

  # 2. remove unnecessary white margins
  img <- trim_white_image(
    img
  )

  # 3. panel-specific internal composition tweaks
  if (letter == "B") {
    # v7: preserve the source panel's native aspect ratio.
    # The B source script now handles its own fixed geometry and readable legend.
    # Do NOT horizontally rescale the main body/legend after rasterization.
    img <- img
  }

  # 4. protective / aesthetic panel-specific padding
  p <- PAD_AFTER_TRIM[[letter]]

  img <- add_white_padding(
    img,
    top = p["top"],
    right = p["right"],
    bottom = p["bottom"],
    left = p["left"]
  )

  img
}

# ------------------------------------------------------------------------------
# 4. DRAW / EXPORT HELPERS
# ------------------------------------------------------------------------------
draw_raster_with_new_letter <- function(
  img,
  letter,
  fontsize
) {

  grid::grid.newpage()

  grid::grid.raster(
    image = img,
    x = unit(0.5, "npc"),
    y = unit(0.5, "npc"),
    width = unit(0.97, "npc"),
    height = unit(0.97, "npc"),
    just = c(
      "center",
      "center"
    ),
    interpolate = TRUE
  )

  grid::grid.text(
    label = letter,
    x = unit(0.015, "npc"),
    y = unit(0.985, "npc"),
    just = c(
      "left",
      "top"
    ),
    gp = gpar(
      fontfamily = FONT_FAMILY,
      fontface = LETTER_FONTFACE,
      fontsize = fontsize,
      col = "black"
    )
  )
}

save_clean_panel <- function(
  img,
  letter
) {

  base_name <- paste0(
    "Figure4",
    letter,
    "_clean_v16_E_noP_TNF_NFKB_BH"
  )

  pdf_out <- file.path(
    OUTDIR,
    paste0(base_name, ".pdf")
  )

  png_out <- file.path(
    OUTDIR,
    paste0(base_name, ".png")
  )

  tif_out <- file.path(
    OUTDIR,
    paste0(base_name, ".tif")
  )

  h <- dim(img)[1]
  w <- dim(img)[2]

  width_in <- 7.2
  height_in <- width_in * h / w

  if (capabilities("cairo")) {

    grDevices::cairo_pdf(
      filename = pdf_out,
      width = width_in,
      height = height_in,
      family = FONT_FAMILY,
      onefile = FALSE,
      bg = "white"
    )

  } else {

    grDevices::pdf(
      file = pdf_out,
      width = width_in,
      height = height_in,
      onefile = FALSE,
      family = "Helvetica",
      useDingbats = FALSE,
      bg = "white"
    )
  }

  draw_raster_with_new_letter(
    img,
    letter,
    LETTER_FONTSIZE_SINGLE
  )

  grDevices::dev.off()

  grDevices::png(
    filename = png_out,
    width = round(
      width_in * 400
    ),
    height = round(
      height_in * 400
    ),
    res = 400,
    bg = "white"
  )

  draw_raster_with_new_letter(
    img,
    letter,
    LETTER_FONTSIZE_SINGLE
  )

  grDevices::dev.off()

  grDevices::tiff(
    filename = tif_out,
    width = round(
      width_in * 600
    ),
    height = round(
      height_in * 600
    ),
    res = 600,
    compression = "lzw",
    bg = "white"
  )

  draw_raster_with_new_letter(
    img,
    letter,
    LETTER_FONTSIZE_SINGLE
  )

  grDevices::dev.off()

  c(
    pdf = pdf_out,
    png = png_out,
    tif = tif_out
  )
}

draw_panel_into_viewport <- function(
  img,
  letter,
  vp,
  x_shift = 0,
  y_shift = 0,
  scale_mult = 1
) {

  grid::pushViewport(vp)

  img_h <- dim(img)[1]
  img_w <- dim(img)[2]

  img_aspect <- img_w / img_h

  vp_w <- convertWidth(
    unit(1, "npc"),
    "in",
    valueOnly = TRUE
  )

  vp_h <- convertHeight(
    unit(1, "npc"),
    "in",
    valueOnly = TRUE
  )

  vp_aspect <- vp_w / vp_h

  if (
    img_aspect >= vp_aspect
  ) {

    draw_w <- 0.985
    draw_h <- (
      0.985 *
      vp_aspect /
      img_aspect
    )

  } else {

    draw_h <- 0.985
    draw_w <- (
      0.985 *
      img_aspect /
      vp_aspect
    )
  }

  draw_w <- min(
    1.02,
    draw_w * scale_mult
  )

  draw_h <- min(
    1.02,
    draw_h * scale_mult
  )

  grid::grid.raster(
    image = img,
    x = unit(
      0.5 + x_shift,
      "npc"
    ),
    y = unit(
      0.5 + y_shift,
      "npc"
    ),
    width = unit(
      draw_w,
      "npc"
    ),
    height = unit(
      draw_h,
      "npc"
    ),
    just = c(
      "center",
      "center"
    ),
    interpolate = TRUE
  )

  grid::grid.text(
    label = letter,
    x = unit(
      LETTER_X,
      "npc"
    ),
    y = unit(
      LETTER_Y,
      "npc"
    ),
    just = c(
      "left",
      "top"
    ),
    gp = gpar(
      fontfamily = FONT_FAMILY,
      fontface = LETTER_FONTFACE,
      fontsize = LETTER_FONTSIZE_FINAL,
      col = "black"
    )
  )

  grid::popViewport()
}

# ------------------------------------------------------------------------------
# 5. FINAL MANUAL LAYOUT
# ------------------------------------------------------------------------------
draw_final_figure4 <- function(
  panel_imgs,
  letters
) {

  grid::grid.newpage()

  # page edges
  x0 <- 0.018
  x1 <- 0.985
  y0 <- 0.015
  y1 <- 0.985

  W <- x1 - x0
  H <- y1 - y0

  # ---------------------------------------------------------------------------
  # TOP ROW: A + B
  #
  # B receives the larger half because its two heatmaps need more width.
  # ---------------------------------------------------------------------------
  y_top_1 <- y1
  y_top_0 <- y1 - H * 0.305

  xA0 <- x0
  xA1 <- x0 + W * 0.385

  xB0 <- xA1 + W * 0.010
  xB1 <- x1

  # ---------------------------------------------------------------------------
  # MIDDLE BLOCK
  #
  # C large on left.
  # D/E stacked on right.
  # Right side is widened to enlarge D.
  # ---------------------------------------------------------------------------
  y_mid_1 <- y_top_0 - H * 0.010
  y_mid_0 <- y_mid_1 - H * 0.340

  xC0 <- x0
  xC1 <- x0 + W * 0.495

  xR0 <- xC1 + W * 0.012
  xR1 <- x1

  # Give D substantially MORE vertical room.
  yD1 <- y_mid_1
  yD0 <- y_mid_0 + H * 0.150

  yE1 <- yD0 - H * 0.008
  yE0 <- y_mid_0

  # ---------------------------------------------------------------------------
  # F full width
  # ---------------------------------------------------------------------------
  yF1 <- y_mid_0 - H * 0.010
  yF0 <- y0

  # A
  draw_panel_into_viewport(
    panel_imgs[[1]],
    letters[1],
    viewport(
      x = unit(
        (xA0 + xA1) / 2,
        "npc"
      ),
      y = unit(
        (y_top_0 + y_top_1) / 2,
        "npc"
      ),
      width = unit(
        xA1 - xA0,
        "npc"
      ),
      height = unit(
        y_top_1 - y_top_0,
        "npc"
      )
    ),
    scale_mult = 0.93
  )

  # B — deliberately enlarged
  draw_panel_into_viewport(
    panel_imgs[[2]],
    letters[2],
    viewport(
      x = unit(
        (xB0 + xB1) / 2,
        "npc"
      ),
      y = unit(
        (y_top_0 + y_top_1) / 2,
        "npc"
      ),
      width = unit(
        xB1 - xB0,
        "npc"
      ),
      height = unit(
        y_top_1 - y_top_0,
        "npc"
      )
    ),
    scale_mult = 1.16
  )

  # C
  draw_panel_into_viewport(
    panel_imgs[[3]],
    letters[3],
    viewport(
      x = unit(
        (xC0 + xC1) / 2,
        "npc"
      ),
      y = unit(
        (y_mid_0 + y_mid_1) / 2,
        "npc"
      ),
      width = unit(
        xC1 - xC0,
        "npc"
      ),
      height = unit(
        y_mid_1 - y_mid_0,
        "npc"
      )
    ),
    scale_mult = 1.00
  )

  # D — deliberately enlarged; slightly downward to protect top row
  draw_panel_into_viewport(
    panel_imgs[[4]],
    letters[4],
    viewport(
      x = unit(
        (xR0 + xR1) / 2,
        "npc"
      ),
      y = unit(
        (yD0 + yD1) / 2,
        "npc"
      ),
      width = unit(
        xR1 - xR0,
        "npc"
      ),
      height = unit(
        yD1 - yD0,
        "npc"
      )
    ),
    x_shift = 0.018,
    y_shift = -0.022,
    scale_mult = 1.02
  )

  # E
  draw_panel_into_viewport(
    panel_imgs[[5]],
    letters[5],
    viewport(
      x = unit(
        (xR0 + xR1) / 2,
        "npc"
      ),
      y = unit(
        (yE0 + yE1) / 2,
        "npc"
      ),
      width = unit(
        xR1 - xR0,
        "npc"
      ),
      height = unit(
        yE1 - yE0,
        "npc"
      )
    ),
    scale_mult = 1.04
  )

  # F
  draw_panel_into_viewport(
    panel_imgs[[6]],
    letters[6],
    viewport(
      x = unit(
        (x0 + x1) / 2,
        "npc"
      ),
      y = unit(
        (yF0 + yF1) / 2,
        "npc"
      ),
      width = unit(
        x1 - x0,
        "npc"
      ),
      height = unit(
        yF1 - yF0,
        "npc"
      )
    ),
    scale_mult = 1.00
  )
}

# ------------------------------------------------------------------------------
# 6. FIND SIX PANEL SCRIPTS
# ------------------------------------------------------------------------------
SCRIPT_DIR <- get_script_dir()

script_candidates <- list(

  A = c(
    file.path(
      SCRIPT_DIR,
      "86A_Figure4_PanelA_global_CellChat_network_v8_stronger.R"
    ),
    file.path(
      "<DOWNLOADS>",
      "86A_Figure4_PanelA_global_CellChat_network_v8_stronger.R"
    )
  ),

  B = c(
    file.path(
      SCRIPT_DIR,
      "86B_Figure4_PanelB_7dpi_reference_style_nondifferential_onepage_v8_GLIA_readability.R"
    ),
    file.path(
      SCRIPT_DIR,
      "86B_Figure4_PanelB_7dpi_reference_style_nondifferential_onepage_v7_manual_grid(1).R"
    ),
    file.path(
      SCRIPT_DIR,
      "86B_Figure4_PanelB_7dpi_reference_style_nondifferential_onepage_v7_manual_grid.R"
    ),
    file.path(
      "<DOWNLOADS>",
      "86B_Figure4_PanelB_7dpi_reference_style_nondifferential_onepage_v7_manual_grid(1).R"
    ),
    file.path(
      "<DOWNLOADS>",
      "86B_Figure4_PanelB_7dpi_reference_style_nondifferential_onepage_v7_manual_grid.R"
    )
  ),

  C = c(
    file.path(
      SCRIPT_DIR,
      "86C_Figure4_PanelC_reference_style_bubbleplot_v4_GLIA_readability.R"
    ),
    file.path(
      SCRIPT_DIR,
      "86C_Figure4_PanelC_reference_style_bubbleplot_v3_rewrite.R"
    ),
    file.path(
      "<DOWNLOADS>",
      "86C_Figure4_PanelC_reference_style_bubbleplot_v3_rewrite.R"
    )
  ),

  D = c(
    file.path(
      SCRIPT_DIR,
      "86D_Figure4_PanelD_full_violin_v7_GLIA_readability.R"
    ),
    file.path(
      SCRIPT_DIR,
      "86D_Figure4_PanelD_full_violin_v6_side_by_side_fix.R"
    ),
    file.path(
      "<DOWNLOADS>",
      "86D_Figure4_PanelD_full_violin_v6_side_by_side_fix.R"
    )
  ),

  E = c(
    file.path(
      ROOT,
      "86E_Figure4_PanelE_reference_style_dual_scatter_v9_no_pooledP.R"
    ),
    file.path(
      SCRIPT_DIR,
      "86E_Figure4_PanelE_reference_style_dual_scatter_v9_no_pooledP.R"
    ),
    file.path(
      "<DOWNLOADS>",
      "86E_Figure4_PanelE_reference_style_dual_scatter_v9_no_pooledP.R"
    ),
    file.path(
      ROOT,
      "86E_Figure4_PanelE_reference_style_dual_scatter_v8_topTick_fix.R"
    ),
    file.path(
      SCRIPT_DIR,
      "86E_Figure4_PanelE_reference_style_dual_scatter_v8_topTick_fix.R"
    ),
    file.path(
      SCRIPT_DIR,
      "86E_Figure4_PanelE_reference_style_dual_scatter_v7_GLIA_readability_FIXED.R"
    ),
    file.path(
      SCRIPT_DIR,
      "86E_Figure4_PanelE_reference_style_dual_scatter_v7_GLIA_readability(1).R"
    ),
    file.path(
      SCRIPT_DIR,
      "86E_Figure4_PanelE_reference_style_dual_scatter_v7_GLIA_readability.R"
    ),
    file.path(
      SCRIPT_DIR,
      "86E_Figure4_PanelE_reference_style_dual_scatter_v6_AUCell_allcellFit(1).R"
    ),
    file.path(
      SCRIPT_DIR,
      "86E_Figure4_PanelE_reference_style_dual_scatter_v6_AUCell_allcellFit.R"
    ),
    file.path(
      "<DOWNLOADS>",
      "86E_Figure4_PanelE_reference_style_dual_scatter_v6_AUCell_allcellFit(1).R"
    ),
    file.path(
      "<DOWNLOADS>",
      "86E_Figure4_PanelE_reference_style_dual_scatter_v6_AUCell_allcellFit.R"
    )
  ),

  F = c(
    file.path(
      SCRIPT_DIR,
      "86F_Figure4_PanelF_ONE_FILE_v8_TNF_NFKB_BH.R"
    ),
    file.path(
      ROOT,
      "86F_Figure4_PanelF_ONE_FILE_v8_TNF_NFKB_BH.R"
    ),
    file.path(
      "<DOWNLOADS>",
      "86F_Figure4_PanelF_ONE_FILE_v8_TNF_NFKB_BH.R"
    )
  )
)

script_paths <- vapply(
  script_candidates,
  find_first_existing,
  character(1)
)

# Safety check: the final GLIA assembly must use the revised Figure4E v9 script.
# This prevents silent reversion to an older panel that displays pooled cell-level P values.
if (!is.na(script_paths["E"]) &&
    !grepl("v9_no_pooledP", basename(script_paths["E"]), fixed = TRUE)) {
  stop(
    "Figure4E v9 no-pooled-P script was not found.\n",
    "Resolved E script instead to:\n  ", script_paths["E"], "\n",
    "Place 86E_Figure4_PanelE_reference_style_dual_scatter_v9_no_pooledP.R ",
    "beside this merge script, under ROOT, or in Downloads, then rerun.",
    call. = FALSE
  )
}

if (
  any(
    is.na(script_paths)
  )
) {

  miss <- names(
    script_paths
  )[
    is.na(script_paths)
  ]

  stop(
    "Could not find panel script(s): ",
    paste(
      miss,
      collapse = ", "
    ),
    "\nPlace all six scripts in Downloads or beside this merge script.",
    call. = FALSE
  )
}

# ------------------------------------------------------------------------------
# 7. RUN + PREPARE PANELS
# ------------------------------------------------------------------------------
letters <- c(
  "A",
  "B",
  "C",
  "D",
  "E",
  "F"
)

panel_names <- names(
  script_paths
)

source_pdf_paths <- setNames(
  rep(
    NA_character_,
    6
  ),
  panel_names
)

panel_imgs <- vector(
  "list",
  6
)

for (
  i in seq_along(
    panel_names
  )
) {

  nm <- panel_names[i]

  cat(
    "\n============================================================\n"
  )

  cat(
    "Running panel ",
    nm,
    ":\n",
    sep = ""
  )

  cat(
    "  ",
    script_paths[nm],
    "\n",
    sep = ""
  )

  if (nm == "F") {
    # The vetted Figure4F v8 standalone plot is reused as-is.
    # Do not source F here: this avoids needless recomputation or an old v7 fallback.
    source_pdf_paths[nm] <- normalizePath(
      F_V8_PDF, winslash = "/", mustWork = TRUE
    )
    cat("Reusing revised Figure4F v8 PDF (TNF/NF-kB; BH-adjusted):\n  ",
        source_pdf_paths[nm], "\n", sep = "")
  } else {
    source_pdf_paths[nm] <- source_panel_script(
      script_paths[nm]
    )
  }

  raw_img <- pdf_first_page_to_image(
    source_pdf_paths[nm],
    dpi = 450
  )

  panel_imgs[[i]] <- prepare_panel_image(
    raw_img,
    letters[i]
  )

  cat(
    "Prepared panel ",
    letters[i],
    " after source-letter erasure + white-margin trim.\n",
    sep = ""
  )
}

names(
  panel_imgs
) <- panel_names

# ------------------------------------------------------------------------------
# 8. SAVE CLEAN SINGLE PANELS
# ------------------------------------------------------------------------------
single_records <- list()

for (
  i in seq_along(
    letters
  )
) {

  outs <- save_clean_panel(
    panel_imgs[[i]],
    letters[i]
  )

  single_records[[i]] <- data.frame(
    panel = letters[i],
    source_script = script_paths[i],
    source_pdf = source_pdf_paths[i],
    clean_pdf = outs["pdf"],
    clean_png = outs["png"],
    clean_tif = outs["tif"],
    stringsAsFactors = FALSE
  )
}

# ------------------------------------------------------------------------------
# 9. SAVE FINAL FIGURE
# ------------------------------------------------------------------------------
for (
  f in c(
    FINAL_PDF,
    FINAL_PNG,
    FINAL_TIF
  )
) {
  if (
    file.exists(f)
  ) {
    unlink(f)
  }
}

if (
  capabilities("cairo")
) {

  grDevices::cairo_pdf(
    filename = FINAL_PDF,
    width = FIG_WIDTH_IN,
    height = FIG_HEIGHT_IN,
    family = FONT_FAMILY,
    onefile = FALSE,
    bg = "white"
  )

} else {

  grDevices::pdf(
    file = FINAL_PDF,
    width = FIG_WIDTH_IN,
    height = FIG_HEIGHT_IN,
    family = "Helvetica",
    onefile = FALSE,
    useDingbats = FALSE,
    bg = "white"
  )
}

draw_final_figure4(
  panel_imgs,
  letters
)

grDevices::dev.off()

grDevices::png(
  filename = FINAL_PNG,
  width = round(
    FIG_WIDTH_IN * 450
  ),
  height = round(
    FIG_HEIGHT_IN * 450
  ),
  res = 450,
  bg = "white"
)

draw_final_figure4(
  panel_imgs,
  letters
)

grDevices::dev.off()

grDevices::tiff(
  filename = FINAL_TIF,
  width = round(
    FIG_WIDTH_IN * 600
  ),
  height = round(
    FIG_HEIGHT_IN * 600
  ),
  res = 600,
  compression = "lzw",
  bg = "white"
)

draw_final_figure4(
  panel_imgs,
  letters
)

grDevices::dev.off()

# ------------------------------------------------------------------------------
# 10. AUDIT
# ------------------------------------------------------------------------------
audit_df <- do.call(
  rbind,
  single_records
)

audit_extra <- data.frame(
  panel = "Figure4",
  source_script = NA,
  source_pdf = NA,
  clean_pdf = FINAL_PDF,
  clean_png = FINAL_PNG,
  clean_tif = FINAL_TIF,
  stringsAsFactors = FALSE
)

audit_all <- rbind(
  audit_df,
  audit_extra
)

write.csv(
  audit_all,
  AUDIT_CSV,
  row.names = FALSE
)

cat(
  "\n============================================================\n"
)

cat(
  "Figure 4 v16 GLIA assembly completed: revised E (no pooled cell-level P) + TNF/NF-kB/BH F.\n"
)

cat(
  "- source letters erased BEFORE scaling\n"
)

cat(
  "- white margins trimmed\n"
)

cat(
  "- A reduced; B enlarged\n"
)

cat(
  "- B native aspect ratio preserved; revised source controls legend/text\n"
)

cat(
  "- D old-letter mask minimized further; extra top/left padding added\n"
)

cat(
  "- E uses revised v9 display: Spearman rho only; pooled cell-level P values omitted\n"
)

cat(
  "- F reuses the revised v8 TNF/NF-kB PDF; source letter erased\n\n"
)

cat(
  "PDF : ",
  FINAL_PDF,
  "\n",
  sep = ""
)

cat(
  "PNG : ",
  FINAL_PNG,
  "\n",
  sep = ""
)

cat(
  "TIFF: ",
  FINAL_TIF,
  "\n",
  sep = ""
)

cat(
  "============================================================\n"
)
