# ==============================================================================
# 84_Figure2_GLIA_complete_AtoF_master_v10_CDsplit_Ev6.R
# v10 adapts the master to the revised Figure 2 structure:
# - old Panel C is removed;
# - revised Panel C contains three side-by-side timepoint maps (1/3/7 dpi);
# - revised Panel D contains pseudotime only;
# - revised Panel E uses Mod1-Mod5 labels with increased left spacing;
# - Panel B remains tag-free and is not masked.
#
# Figure 2 master assembly | APPROVED panels A-F | submission layout
#
# DESIGN PRINCIPLE
# ----------------
# - Five approved source scripts generate six final panels A-F.
# - The revised C/D source script is run once and exports both Panel C and Panel D.
# - Each standalone panel is exported first by its own approved script.
# - The COMPLETE Figure 2 is assembled from those standalone PNG outputs.
# - White outer margins are trimmed for assembly only.
# - Every panel is drawn inside its own clipped viewport -> NO panel overlap.
#
# JOURNAL-ORIENTED OUTPUT
# -----------------------
# GLIA currently requests high-quality figures and accepts individual figures
# at 300 dpi for submission. The journal does not specify one mandatory figure
# width in mm in its current public Author Guidelines.
#
# Therefore this master uses a conservative full-width production canvas:
#   180 mm x 160 mm
#
# Complete outputs:
#   PDF  : 180 x 160 mm
#   SVG  : 180 x 160 mm
#   TIFF : 600 dpi
#   PNG  : 300 dpi
#
# DEFAULT LAYOUT
# --------------
# Row 1: A | B
# Row 2: C | D
# Row 3: E | F
#
# OUTPUT
# ------
# <PROJECT_ROOT>/09_final_figures/GLIA/
#   Figure2_complete_A-F_v10_CDsplit_Ev6/
#
#       Figure2_complete_A-F_v10_CDsplit_Ev6.pdf
#       Figure2_complete_A-F_v10_CDsplit_Ev6.svg
#       Figure2_complete_A-F_v10_CDsplit_Ev6.tiff
#       Figure2_complete_A-F_v10_CDsplit_Ev6.png
#
#       panels/
#           Figure2A_final.pdf/svg/tiff/png
#           ...
#           Figure2F_final.pdf/svg/tiff/png
#
#       00_panel_script_audit.csv
#       01_panel_output_manifest.csv
#       02_master_layout_manifest.csv
#       03_trim_manifest.csv
#       04_submission_preflight.csv
#       sessionInfo.txt
# ==============================================================================

rm(list = ls())
gc()
options(stringsAsFactors = FALSE)

# ------------------------------------------------------------------------------
# 0. USER SETTINGS
# ------------------------------------------------------------------------------
project_root <- Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
panel_script_dir <- project_root

out_dir <- file.path(
  project_root,
  "09_final_figures", "GLIA",
  "Figure2_complete_A-F_v10_CDsplit_Ev6"
)

panel_copy_dir <- file.path(out_dir, "panels")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(panel_copy_dir, recursive = TRUE, showWarnings = FALSE)

figure_width_mm  <- 180
figure_height_mm <- 160

margin_left_mm   <- 3.0
margin_right_mm  <- 3.0
margin_top_mm    <- 3.5
margin_bottom_mm <- 3.5

gap_x_mm <- 3.0
gap_y_mm <- 4

row1_h_mm <- 42
row2_h_mm <- 49
row3_h_mm <- 54

row1_A_w_mm <- 116.5
row1_B_w_mm <- 174 - gap_x_mm - row1_A_w_mm  # 54.5

row2_C_w_mm <- 118.0
row2_D_w_mm <- 174 - gap_x_mm - row2_C_w_mm  # 53.0

row3_E_w_mm <- 87.0
row3_F_w_mm <- 174 - gap_x_mm - row3_E_w_mm  # 84

trim_white_threshold <- 0.992
trim_padding_px <- 12
b_extra_top_padding_px <- 24

# ------------------------------------------------------------------------------
# 1. PACKAGES
# ------------------------------------------------------------------------------
required_pkgs <- c("grid", "png", "svglite", "ragg")

missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_pkgs) > 0) {
  stop(
    "Missing required master-assembly package(s): ",
    paste(missing_pkgs, collapse = ", "),
    "\nInstall them first, then rerun."
  )
}

suppressPackageStartupMessages({
  library(grid)
})

# ------------------------------------------------------------------------------
# 2. RESOLVE THE FIVE APPROVED SOURCE SCRIPTS
# ------------------------------------------------------------------------------
resolve_script <- function(script_dir, candidates) {

  # 1) Prefer files directly under the xjz project root.
  candidate_paths <- file.path(script_dir, candidates)
  hit <- candidate_paths[file.exists(candidate_paths)]

  if (length(hit) > 0) {
    return(normalizePath(hit[1], winslash = "/", mustWork = TRUE))
  }

  # 2) Search recursively, but still require an exact approved filename.
  all_r <- list.files(
    script_dir,
    pattern = "\\.R$",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )

  for (cand in candidates) {
    hit2 <- all_r[basename(all_r) == cand]
    if (length(hit2) > 0) {
      return(normalizePath(hit2[1], winslash = "/", mustWork = TRUE))
    }
  }

  stop(
    "Could not find approved Figure 2 source script under xjz.\n",
    "Project root: ", script_dir, "\n",
    "Required filename(s):\n  ",
    paste(candidates, collapse = "\n  "),
    "\n\nPut the required scripts anywhere inside <PROJECT_ROOT> and rerun."
  )
}

source_scripts <- c(
  A = resolve_script(
    panel_script_dir,
    c(
      "83A_Figure2_PanelA_resident_AUCell_bridge_v6_GLIA_countfix.R"
    )
  ),

  B = resolve_script(
    panel_script_dir,
    c(
      "83B_Figure2_PanelB_shared_transition_gate_v3_mastertag_safe.R"
    )
  ),

  CD = resolve_script(
    panel_script_dir,
    c(
      "83D_Figure2_PanelsC_D_timepoint_pseudotime_v6.R"
    )
  ),

  E = resolve_script(
    panel_script_dir,
    c(
      "83E_Figure2_PanelE_module_heatmap_simple_v6_Mod_spacing.R"
    )
  ),

  F = resolve_script(
    panel_script_dir,
    c(
      "83F_Figure2_PanelF_GO_KEGG_dotplot_v2_GLIA_readability.R"
    )
  )
)

panel_script_audit <- data.frame(
  source_key = names(source_scripts),
  produces_panel = c("A", "B", "C + D", "E", "F"),
  script = unname(source_scripts),
  exists = file.exists(source_scripts),
  stringsAsFactors = FALSE
)

write.csv(
  panel_script_audit,
  file.path(out_dir, "00_panel_script_audit.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 3. RUN EACH APPROVED SOURCE SCRIPT IN ITS OWN ENVIRONMENT
# ------------------------------------------------------------------------------
run_source_script <- function(source_key, script_path) {

  message("")
  message("============================================================")
  message("Running approved Figure 2 source: ", source_key)
  message(script_path)
  message("============================================================")

  ee <- new.env(parent = globalenv())

  sys.source(
    script_path,
    envir = ee,
    keep.source = FALSE
  )

  ee
}

envA  <- run_source_script("A",  source_scripts[["A"]])
envB  <- run_source_script("B",  source_scripts[["B"]])
envCD <- run_source_script("C/D", source_scripts[["CD"]])
envE  <- run_source_script("E",  source_scripts[["E"]])
envF  <- run_source_script("F",  source_scripts[["F"]])

# Standard single-panel scripts must expose out_dir + stem.
for (nm in c("A", "B", "E", "F")) {

  ee <- get(paste0("env", nm))

  miss <- c("out_dir", "stem")[
    !vapply(
      c("out_dir", "stem"),
      exists,
      logical(1),
      envir = ee,
      inherits = FALSE
    )
  ]

  if (length(miss) > 0) {
    stop(
      "Figure 2 source ", nm,
      " did not leave expected object(s): ",
      paste(miss, collapse = ", ")
    )
  }
}

# Revised C/D source intentionally exports two independent final panel stems.
if (!exists("out_dir", envir = envCD, inherits = FALSE)) {
  stop("Revised C/D source did not leave expected object: out_dir")
}

panel_defs <- list(
  A = list(
    out_dir = envA$out_dir,
    stem = envA$stem
  ),
  B = list(
    out_dir = envB$out_dir,
    stem = envB$stem
  ),
  C = list(
    out_dir = envCD$out_dir,
    stem = "Figure2C_timepoint_split_v6"
  ),
  D = list(
    out_dir = envCD$out_dir,
    stem = "Figure2D_pseudotime_v6"
  ),
  E = list(
    out_dir = envE$out_dir,
    stem = envE$stem
  ),
  F = list(
    out_dir = envF$out_dir,
    stem = envF$stem
  )
)

# ------------------------------------------------------------------------------
# 4. VALIDATE STANDALONE OUTPUTS AND COPY THEM
# ------------------------------------------------------------------------------
copy_panel_outputs <- function(panel_letter, def) {

  exts <- c(".pdf", ".svg", ".tiff", ".png")

  srcs <- file.path(
    def$out_dir,
    paste0(def$stem, exts)
  )

  dsts <- file.path(
    panel_copy_dir,
    paste0("Figure2", panel_letter, "_final", exts)
  )

  exists_vec <- file.exists(srcs)
  copied <- rep(FALSE, length(srcs))

  if (any(exists_vec)) {
    copied[exists_vec] <- file.copy(
      srcs[exists_vec],
      dsts[exists_vec],
      overwrite = TRUE
    )
  }

  data.frame(
    panel = panel_letter,
    extension = exts,
    source = srcs,
    source_exists = exists_vec,
    central_copy = dsts,
    copied = copied,
    stringsAsFactors = FALSE
  )
}

copy_manifest <- do.call(
  rbind,
  lapply(
    names(panel_defs),
    function(pp) copy_panel_outputs(pp, panel_defs[[pp]])
  )
)

write.csv(
  copy_manifest,
  file.path(out_dir, "01_panel_output_manifest.csv"),
  row.names = FALSE
)

if (any(!copy_manifest$source_exists)) {
  print(copy_manifest[!copy_manifest$source_exists, ])
  stop("At least one standalone final panel output is missing.")
}

if (any(copy_manifest$source_exists & !copy_manifest$copied)) {
  print(copy_manifest[copy_manifest$source_exists & !copy_manifest$copied, ])
  stop("At least one standalone final panel output could not be copied.")
}

# ------------------------------------------------------------------------------
# 5. READ + TRIM STANDALONE PNG PANELS FOR MASTER ASSEMBLY
# ------------------------------------------------------------------------------
read_panel_png <- function(panel_letter, def) {

  png_path <- file.path(
    def$out_dir,
    paste0(def$stem, ".png")
  )

  if (!file.exists(png_path)) {
    stop(
      "Panel ", panel_letter,
      " PNG missing after source script: ",
      png_path
    )
  }

  img <- png::readPNG(png_path)

  list(
    panel = panel_letter,
    path = png_path,
    image = img
  )
}

trim_white_png_array <- function(
  img,
  threshold = 0.992,
  padding = 18
) {

  dd <- dim(img)

  if (length(dd) != 3 || dd[3] < 3) {
    stop("Unexpected PNG array dimensions.")
  }

  rgb <- img[, , 1:3, drop = FALSE]

  content_mask <- apply(
    rgb,
    c(1, 2),
    function(z) any(z < threshold)
  )

  idx <- which(content_mask, arr.ind = TRUE)

  if (nrow(idx) == 0) {
    return(
      list(
        image = img,
        crop = c(
          top = 1,
          bottom = dd[1],
          left = 1,
          right = dd[2]
        )
      )
    )
  }

  r1 <- max(1, min(idx[, 1]) - padding)
  r2 <- min(dd[1], max(idx[, 1]) + padding)
  c1 <- max(1, min(idx[, 2]) - padding)
  c2 <- min(dd[2], max(idx[, 2]) + padding)

  list(
    image = img[r1:r2, c1:c2, , drop = FALSE],
    crop = c(
      top = r1,
      bottom = r2,
      left = c1,
      right = c2
    )
  )
}

panel_pngs <- lapply(
  names(panel_defs),
  function(pp) read_panel_png(pp, panel_defs[[pp]])
)
names(panel_pngs) <- names(panel_defs)

trimmed_panels <- list()
trim_manifest <- list()

for (pp in names(panel_pngs)) {

  rr <- panel_pngs[[pp]]

  tt <- trim_white_png_array(
    rr$image,
    threshold = trim_white_threshold,
    padding = trim_padding_px
  )

  # Panel B needs a little extra room above the plot area so the upper end of
  # its vertical y-axis title is not clipped after trimming.
  if (pp == "B") {
    d_full <- dim(rr$image)
    new_top <- max(1, tt$crop["top"] - b_extra_top_padding_px)
    tt$image <- rr$image[
      new_top:tt$crop["bottom"],
      tt$crop["left"]:tt$crop["right"],
      ,
      drop = FALSE
    ]
    tt$crop["top"] <- new_top
  }

  trimmed_panels[[pp]] <- tt$image

  d0 <- dim(rr$image)
  d1 <- dim(tt$image)

  trim_manifest[[pp]] <- data.frame(
    panel = pp,
    source_png = rr$path,
    source_width_px = d0[2],
    source_height_px = d0[1],
    trimmed_width_px = d1[2],
    trimmed_height_px = d1[1],
    crop_top_px = tt$crop["top"],
    crop_bottom_px = tt$crop["bottom"],
    crop_left_px = tt$crop["left"],
    crop_right_px = tt$crop["right"],
    stringsAsFactors = FALSE
  )
}

trim_manifest_df <- do.call(rbind, trim_manifest)

write.csv(
  trim_manifest_df,
  file.path(out_dir, "03_trim_manifest.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 6. MASTER LAYOUT COORDINATES
# ------------------------------------------------------------------------------
inner_x0 <- margin_left_mm

row3_y0 <- margin_bottom_mm
row2_y0 <- row3_y0 + row3_h_mm + gap_y_mm
row1_y0 <- row2_y0 + row2_h_mm + gap_y_mm

layout_df <- data.frame(
  panel = c("A", "B", "C", "D", "E", "F"),

  x_mm = c(
    inner_x0,
    inner_x0 + row1_A_w_mm + gap_x_mm,
    inner_x0,
    inner_x0 + row2_C_w_mm + gap_x_mm,
    inner_x0,
    inner_x0 + row3_E_w_mm + gap_x_mm
  ),

  y_mm = c(
    row1_y0,
    row1_y0,
    row2_y0,
    row2_y0,
    row3_y0,
    row3_y0
  ),

  width_mm = c(
    row1_A_w_mm,
    row1_B_w_mm,
    row2_C_w_mm,
    row2_D_w_mm,
    row3_E_w_mm,
    row3_F_w_mm
  ),

  height_mm = c(
    row1_h_mm,
    row1_h_mm,
    row2_h_mm,
    row2_h_mm,
    row3_h_mm,
    row3_h_mm
  ),

  stringsAsFactors = FALSE
)

write.csv(
  layout_df,
  file.path(out_dir, "02_master_layout_manifest.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 7. FIT A RASTER PANEL INTO A FIXED CELL WITHOUT DISTORTION
# ------------------------------------------------------------------------------
draw_panel_raster <- function(
  img,
  x_mm,
  y_mm,
  width_mm,
  height_mm
) {

  dd <- dim(img)
  img_h <- dd[1]
  img_w <- dd[2]

  img_aspect <- img_w / img_h
  box_aspect <- width_mm / height_mm

  if (img_aspect >= box_aspect) {
    draw_w <- width_mm
    draw_h <- width_mm / img_aspect
  } else {
    draw_h <- height_mm
    draw_w <- height_mm * img_aspect
  }

  pushViewport(
    viewport(
      x = unit(x_mm, "mm"),
      y = unit(y_mm, "mm"),
      width = unit(width_mm, "mm"),
      height = unit(height_mm, "mm"),
      just = c("left", "bottom"),
      clip = "on"
    )
  )

  grid.raster(
    img,
    x = unit((width_mm - draw_w) / 2, "mm"),
    y = unit((height_mm - draw_h) / 2, "mm"),
    width = unit(draw_w, "mm"),
    height = unit(draw_h, "mm"),
    just = c("left", "bottom"),
    interpolate = TRUE
  )

  popViewport()
}

# ------------------------------------------------------------------------------
# 8. PLACEMENT + EFFECTIVE RESOLUTION PREFLIGHT
# ------------------------------------------------------------------------------
calc_fit <- function(img, width_mm, height_mm) {

  dd <- dim(img)
  img_h <- dd[1]
  img_w <- dd[2]

  img_aspect <- img_w / img_h
  box_aspect <- width_mm / height_mm

  if (img_aspect >= box_aspect) {
    draw_w <- width_mm
    draw_h <- width_mm / img_aspect
  } else {
    draw_h <- height_mm
    draw_w <- height_mm * img_aspect
  }

  effective_dpi_x <- img_w / (draw_w / 25.4)
  effective_dpi_y <- img_h / (draw_h / 25.4)

  data.frame(
    draw_width_mm = draw_w,
    draw_height_mm = draw_h,
    effective_dpi_x = effective_dpi_x,
    effective_dpi_y = effective_dpi_y,
    effective_dpi_min = min(effective_dpi_x, effective_dpi_y)
  )
}

placement_audit <- do.call(
  rbind,
  lapply(seq_len(nrow(layout_df)), function(ii) {

    pp <- layout_df$panel[ii]

    fit <- calc_fit(
      trimmed_panels[[pp]],
      layout_df$width_mm[ii],
      layout_df$height_mm[ii]
    )

    data.frame(
      panel = pp,
      allocated_width_mm = layout_df$width_mm[ii],
      allocated_height_mm = layout_df$height_mm[ii],
      fit,
      stringsAsFactors = FALSE
    )
  })
)

write.csv(
  placement_audit,
  file.path(out_dir, "04_panel_placement_resolution.csv"),
  row.names = FALSE
)

if (any(placement_audit$effective_dpi_min < 300)) {
  warning(
    "At least one panel has effective source resolution below 300 dpi. ",
    "See 04_panel_placement_resolution.csv."
  )
}

# ------------------------------------------------------------------------------
# 9. DRAW COMPLETE FIGURE 2
# ------------------------------------------------------------------------------
draw_master_figure <- function() {

  grid.newpage()

  grid.rect(
    gp = gpar(
      fill = "white",
      col = NA
    )
  )

  for (ii in seq_len(nrow(layout_df))) {

    pp <- layout_df$panel[ii]
    img_i <- trimmed_panels[[pp]]

    draw_panel_raster(
      img = img_i,
      x_mm = layout_df$x_mm[ii],
      y_mm = layout_df$y_mm[ii],
      width_mm = layout_df$width_mm[ii],
      height_mm = layout_df$height_mm[ii]
    )

    # Locate the actual aspect-ratio-preserved raster within its panel box.
    fit_i <- calc_fit(
      img_i,
      layout_df$width_mm[ii],
      layout_df$height_mm[ii]
    )
    draw_x_i <- layout_df$x_mm[ii] +
      (layout_df$width_mm[ii] - fit_i$draw_width_mm) / 2
    draw_y_i <- layout_df$y_mm[ii] +
      (layout_df$height_mm[ii] - fit_i$draw_height_mm) / 2
    draw_top_i <- draw_y_i + fit_i$draw_height_mm

    # Remove the small source-script panel letter only.
    # Panel B v3 is already tag-free, so do NOT paint a white mask over B.
    # The old mask overlapped the upper part of B's vertical y-axis title.
    if (pp != "B") {
      grid.rect(
        x = unit(draw_x_i, "mm"),
        y = unit(draw_top_i, "mm"),
        width = unit(8.5, "mm"),
        height = unit(7.5, "mm"),
        just = c("left", "top"),
        gp = gpar(fill = "white", col = NA)
      )
    }

    # Uniform manuscript-level panel letter, matching Figures 1 and 4.
    grid.text(
      pp,
      x = unit(layout_df$x_mm[ii] + 0.7, "mm"),
      y = unit(layout_df$y_mm[ii] + layout_df$height_mm[ii] - 0.7, "mm"),
      just = c("left", "top"),
      gp = gpar(
        fontfamily = "Arial",
        fontface = "bold",
        fontsize = 16,
        col = "black"
      )
    )
  }
}

# ------------------------------------------------------------------------------
# 10. SAVE MASTER OUTPUTS
# ------------------------------------------------------------------------------
master_stem <- "Figure2_complete_A-F_v10_CDsplit_Ev6"

if (!capabilities("cairo")) {
  stop(
    "This R installation does not report Cairo support. ",
    "The submission master requires cairo_pdf()."
  )
}

grDevices::cairo_pdf(
  filename = file.path(out_dir, paste0(master_stem, ".pdf")),
  width = figure_width_mm / 25.4,
  height = figure_height_mm / 25.4,
  onefile = FALSE,
  bg = "white"
)
draw_master_figure()
dev.off()

svglite::svglite(
  file = file.path(out_dir, paste0(master_stem, ".svg")),
  width = figure_width_mm / 25.4,
  height = figure_height_mm / 25.4,
  bg = "white"
)
draw_master_figure()
dev.off()

ragg::agg_tiff(
  filename = file.path(out_dir, paste0(master_stem, ".tiff")),
  width = figure_width_mm,
  height = figure_height_mm,
  units = "mm",
  res = 600,
  compression = "lzw",
  background = "white"
)
draw_master_figure()
dev.off()

ragg::agg_png(
  filename = file.path(out_dir, paste0(master_stem, ".png")),
  width = figure_width_mm,
  height = figure_height_mm,
  units = "mm",
  res = 300,
  background = "white"
)
draw_master_figure()
dev.off()

# ------------------------------------------------------------------------------
# 11. SUBMISSION PREFLIGHT
# ------------------------------------------------------------------------------
tiff_px_width  <- round(figure_width_mm  / 25.4 * 600)
tiff_px_height <- round(figure_height_mm / 25.4 * 600)

png_px_width  <- round(figure_width_mm  / 25.4 * 300)
png_px_height <- round(figure_height_mm / 25.4 * 300)

preflight <- data.frame(
  item = c(
    "master_width_mm",
    "master_height_mm",
    "layout_row1_height_mm",
    "layout_row2_height_mm",
    "layout_row3_height_mm",
    "TIFF_dpi",
    "TIFF_width_px",
    "TIFF_height_px",
    "PNG_dpi",
    "PNG_width_px",
    "PNG_height_px",
    "panel_overlap_possible",
    "assembly_source",
    "panel_C_structure",
    "panel_D_structure",
    "panel_E_module_labels"
  ),
  value = c(
    figure_width_mm,
    figure_height_mm,
    row1_h_mm,
    row2_h_mm,
    row3_h_mm,
    600,
    tiff_px_width,
    tiff_px_height,
    300,
    png_px_width,
    png_px_height,
    "FALSE - each panel uses a clipped fixed viewport",
    "standalone panel PNGs + vector 16-pt panel-letter overlay",
    "three side-by-side maps: 1 dpi / 3 dpi / 7 dpi",
    "corrected Step5BR2 pseudotime only",
    "Mod1-Mod5 with increased left spacing"
  ),
  stringsAsFactors = FALSE
)

write.csv(
  preflight,
  file.path(out_dir, "05_submission_preflight.csv"),
  row.names = FALSE
)

capture.output(
  sessionInfo(),
  file = file.path(out_dir, "sessionInfo.txt")
)

message("")
message("============================================================")
message("FIGURE 2 A-F GLIA-UNIFIED MASTER v10 (revised C/D + E v6) COMPLETE")
message("============================================================")
message("Master size: ", figure_width_mm, " x ", figure_height_mm, " mm")
message("")
message("Complete figure PDF:")
message(file.path(out_dir, paste0(master_stem, ".pdf")))
message("")
message("Complete figure TIFF (600 dpi):")
message(file.path(out_dir, paste0(master_stem, ".tiff")))
message("")
message("Standalone approved A-F outputs:")
message(panel_copy_dir)
message("")
message("No panel can draw outside its assigned viewport.")
message("============================================================")
