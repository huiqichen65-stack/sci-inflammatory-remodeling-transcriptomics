# ============================================================
# Figure 4A | Global 7 dpi CellChat interaction-count network
# v8 visual-strength refinement
#
# Changes vs v7:
# - SAME displayed 10 compartments
# - SAME frozen 7 dpi full-cell nboot=100 significant communications
# - SAME source × target interaction-count matrix logic
# - NO CellChat analysis rerun
#
# Visual strengthening only:
#   * darker / more saturated node colors
#   * thicker and less transparent edges
#   * slightly larger nodes overall
#   * slightly tighter visual circle impression
#   * title kept, panel letter retained
#
# Output:
#   <PROJECT_ROOT>/09_final_figures/GLIA/Figure4/
#   single_panels/Figure4A/
# ============================================================

options(stringsAsFactors = FALSE)

# ------------------------------------------------------------
# 0. Packages
# ------------------------------------------------------------
required_pkgs <- c("data.table", "dplyr", "CellChat", "scales")

missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_pkgs) > 0) {
  stop(
    "Missing package(s): ",
    paste(missing_pkgs, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(CellChat)
  library(scales)
})

# ------------------------------------------------------------
# 1. Paths
# ------------------------------------------------------------
ROOT <- Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
STEP8B_CONFIRM <- file.path(
  ROOT,
  "04_results",
  "step8B_CONFIRM_GSE162610_7dpi_targeted_fullcell_nboot100"
)

COMM_FILE <- file.path(
  STEP8B_CONFIRM,
  "tables",
  "03_fullcell_targeted_all_communications.tsv"
)

COUNT_FILE <- file.path(
  STEP8B_CONFIRM,
  "tables",
  "01_full_7dpi_cell_counts_by_group_sample.csv"
)

OUTDIR <- file.path(
  ROOT,
  "09_final_figures",
  "GLIA",
  "Figure4",
  "single_panels",
  "Figure4A"
)

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

PNG_FILE <- file.path(
  OUTDIR,
  "Figure4A_global_interaction_count_network_v8_review.png"
)

PDF_FILE <- file.path(
  OUTDIR,
  "Figure4A_global_interaction_count_network_v8.pdf"
)

SVG_FILE <- file.path(
  OUTDIR,
  "Figure4A_global_interaction_count_network_v8.svg"
)

if (file.exists(PNG_FILE)) unlink(PNG_FILE)
if (file.exists(PDF_FILE)) unlink(PDF_FILE)
if (file.exists(SVG_FILE)) unlink(SVG_FILE)

# ------------------------------------------------------------
# 2. Input audit
# ------------------------------------------------------------
input_audit <- data.frame(
  input = c(
    "03_fullcell_targeted_all_communications.tsv",
    "01_full_7dpi_cell_counts_by_group_sample.csv"
  ),
  path = c(COMM_FILE, COUNT_FILE),
  exists = file.exists(c(COMM_FILE, COUNT_FILE)),
  stringsAsFactors = FALSE
)

write.csv(
  input_audit,
  file.path(OUTDIR, "00_input_audit_v8.csv"),
  row.names = FALSE
)

if (any(!input_audit$exists)) {
  print(input_audit)
  stop("Required frozen Figure4A input file is missing.")
}

# ------------------------------------------------------------
# 3. Read frozen results
# ------------------------------------------------------------
comm <- fread(
  COMM_FILE,
  sep = "\t",
  data.table = FALSE
)

cell_counts <- fread(
  COUNT_FILE,
  data.table = FALSE
)

required_comm_cols <- c(
  "timepoint",
  "source",
  "target",
  "interaction_name",
  "probability",
  "p_value"
)

required_count_cols <- c(
  "cellchat_group",
  "n_cells"
)

miss_comm <- setdiff(required_comm_cols, colnames(comm))
miss_count <- setdiff(required_count_cols, colnames(cell_counts))

if (length(miss_comm) > 0) {
  stop(
    "Communication table missing columns: ",
    paste(miss_comm, collapse = ", ")
  )
}

if (length(miss_count) > 0) {
  stop(
    "Cell-count table missing columns: ",
    paste(miss_count, collapse = ", ")
  )
}

# ------------------------------------------------------------
# 4. Use ONLY frozen significant 7 dpi communications
# ------------------------------------------------------------
comm_7dpi <- comm %>%
  filter(
    timepoint == "7dpi",
    is.finite(p_value),
    p_value <= 0.05
  )

if (nrow(comm_7dpi) == 0) {
  stop("No significant 7 dpi communications found.")
}

# ------------------------------------------------------------
# 5. Exclude only the historical Other branch from display
# ------------------------------------------------------------
OLD_OTHER_BRANCH <- "Microglia_Arm_1"

comm_plot <- comm_7dpi %>%
  filter(
    source != OLD_OTHER_BRANCH,
    target != OLD_OTHER_BRANCH
  )

# ------------------------------------------------------------
# 6. Fixed display order (same 10 compartments)
# ------------------------------------------------------------
preferred_raw_order <- c(
  "Microglia_Arm_2",  # Target program
  "Other_myeloid",
  "Microglia_other",
  "OPC",
  "Pericyte_VSMC",
  "Fibroblast",
  "Endothelial",
  "Astrocyte",
  "T_cell",
  "Neutrophil"
)

present_groups <- union(
  unique(comm_plot$source),
  unique(comm_plot$target)
)

node_raw <- preferred_raw_order[
  preferred_raw_order %in% present_groups
]

if (length(node_raw) < 2) {
  stop("Too few groups available for Figure4A.")
}

# ------------------------------------------------------------
# 7. Display labels
# ------------------------------------------------------------
pretty_group <- function(x) {
  out <- x
  out[out == "Microglia_Arm_2"] <- "Target program"
  out[out == "Microglia_other"] <- "Microglia other"
  out[out == "Other_myeloid"] <- "Other myeloid"
  out[out == "T_cell"] <- "T cell"
  out[out == "Pericyte_VSMC"] <- "Pericyte/VSMC"
  out
}

node_display <- pretty_group(node_raw)

# ------------------------------------------------------------
# 8. Build FULL source × target interaction-count matrix
# ------------------------------------------------------------
pair_counts <- comm_plot %>%
  filter(
    source %in% node_raw,
    target %in% node_raw
  ) %>%
  count(
    source,
    target,
    name = "n_interactions"
  )

mat <- matrix(
  0,
  nrow = length(node_raw),
  ncol = length(node_raw),
  dimnames = list(node_display, node_display)
)

for (i in seq_len(nrow(pair_counts))) {
  src_display <- pretty_group(pair_counts$source[i])
  tgt_display <- pretty_group(pair_counts$target[i])
  mat[src_display, tgt_display] <- pair_counts$n_interactions[i]
}

if (sum(mat) <= 0) {
  stop("The global interaction-count matrix is empty.")
}

# ------------------------------------------------------------
# 9. Real vertex weights, compressed for readability
# ------------------------------------------------------------
group_sizes <- cell_counts %>%
  group_by(cellchat_group) %>%
  summarise(
    n_cells = sum(n_cells),
    .groups = "drop"
  )

vertex_ncells <- group_sizes$n_cells[
  match(node_raw, group_sizes$cellchat_group)
]

vertex_ncells[is.na(vertex_ncells)] <- 1
names(vertex_ncells) <- node_display

vertex_weight <- sqrt(vertex_ncells)
vertex_weight <- scales::rescale(
  vertex_weight,
  to = c(7.2, 16.5)
)
names(vertex_weight) <- node_display

# Mild extra emphasis for the core node without changing data logic
vertex_weight["Target program"] <- max(
  vertex_weight["Target program"],
  12.5
)

# ------------------------------------------------------------
# 10. Darker, more saturated color vector
# ------------------------------------------------------------
COL_TARGET <- "#B13A36"

all_colors <- c(
  "Target program"  = COL_TARGET,
  "Other myeloid"   = "#D95F02",
  "Microglia other" = "#7663B0",
  "OPC"             = "#3C6EA8",
  "Pericyte/VSMC"   = "#5E9E1B",
  "Fibroblast"      = "#6B61B8",
  "Endothelial"     = "#159570",
  "Astrocyte"       = "#D81B80",
  "T cell"          = "#149A99",
  "Neutrophil"      = "#D9A300"
)

missing_colors <- setdiff(node_display, names(all_colors))
if (length(missing_colors) > 0) {
  stop(
    "Missing display color for: ",
    paste(missing_colors, collapse = ", ")
  )
}

color_use <- all_colors[node_display]

if (
  is.null(names(color_use)) ||
  !identical(names(color_use), colnames(mat))
) {
  stop("Named color vector does not match network matrix nodes.")
}

# ------------------------------------------------------------
# 11. Export exact data
# ------------------------------------------------------------
write.csv(
  as.data.frame(mat),
  file.path(OUTDIR, "01_global_interaction_count_matrix_v8.csv")
)

write.csv(
  pair_counts,
  file.path(OUTDIR, "02_global_source_target_interaction_counts_v8.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(
    raw_group = node_raw,
    display_group = node_display,
    n_cells = as.numeric(vertex_ncells),
    plot_vertex_weight = as.numeric(vertex_weight),
    color = unname(color_use),
    stringsAsFactors = FALSE
  ),
  file.path(OUTDIR, "03_node_and_color_audit_v8.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 12. Font
# ------------------------------------------------------------
FONT_FAMILY <- "Arial"

if (.Platform$OS.type == "windows") {
  try(
    grDevices::windowsFonts(
      Arial = grDevices::windowsFont("Arial")
    ),
    silent = TRUE
  )
}

# ------------------------------------------------------------
# 13. Draw Figure4A
# ------------------------------------------------------------
draw_Figure4A <- function() {
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar), add = TRUE)

  par(
    mar = c(0.22, 0.22, 2.2, 0.22),
    oma = c(0.02, 0.02, 0.02, 0.02),
    xpd = TRUE,
    family = FONT_FAMILY
  )

  CellChat::netVisual_circle(
    mat,
    color.use = color_use,
    vertex.weight = as.numeric(vertex_weight),
    vertex.weight.max = max(vertex_weight),
    vertex.size.max = 17.5,
    weight.scale = TRUE,
    edge.weight.max = max(mat),
    edge.width.max = 7.2,
    alpha.edge = 0.70,
    label.edge = FALSE,
    vertex.label.cex = 0.82,
    vertex.label.color = "black",
    title.name = "Number of interactions – 7 dpi SCI"
  )

  mtext(
    "A",
    side = 3,
    adj = -0.006,
    line = -0.10,
    cex = 1.95,
    font = 2,
    family = FONT_FAMILY
  )
}

# ------------------------------------------------------------
# 14. Save PNG first
# ------------------------------------------------------------
grDevices::png(
  filename = PNG_FILE,
  width = 2350,
  height = 1800,
  res = 320,
  bg = "white"
)

tryCatch(
  draw_Figure4A(),
  error = function(e) {
    grDevices::dev.off()
    if (file.exists(PNG_FILE)) unlink(PNG_FILE)
    stop(e)
  }
)

grDevices::dev.off()

if (
  !file.exists(PNG_FILE) ||
  is.na(file.info(PNG_FILE)$size) ||
  file.info(PNG_FILE)$size < 10000
) {
  stop("PNG output was not written correctly.")
}

# ------------------------------------------------------------
# 15. Save PDF
# ------------------------------------------------------------
if (capabilities("cairo")) {
  grDevices::cairo_pdf(
    filename = PDF_FILE,
    width = 7.6,
    height = 5.9,
    family = FONT_FAMILY,
    onefile = TRUE,
    bg = "white"
  )
} else {
  grDevices::pdf(
    file = PDF_FILE,
    width = 7.6,
    height = 5.9,
    family = "Helvetica",
    useDingbats = FALSE,
    bg = "white"
  )
}

tryCatch(
  draw_Figure4A(),
  error = function(e) {
    grDevices::dev.off()
    if (file.exists(PDF_FILE)) unlink(PDF_FILE)
    stop(e)
  }
)

grDevices::dev.off()

# ------------------------------------------------------------
# 16. Save SVG
# ------------------------------------------------------------
grDevices::svg(
  filename = SVG_FILE,
  width = 7.6,
  height = 5.9,
  family = FONT_FAMILY,
  onefile = TRUE,
  bg = "white"
)

tryCatch(
  draw_Figure4A(),
  error = function(e) {
    grDevices::dev.off()
    if (file.exists(SVG_FILE)) unlink(SVG_FILE)
    stop(e)
  }
)

grDevices::dev.off()

# ------------------------------------------------------------
# 17. Final audit
# ------------------------------------------------------------
outfiles <- c(PNG_FILE, PDF_FILE, SVG_FILE)

output_audit <- data.frame(
  file = outfiles,
  exists = file.exists(outfiles),
  size_bytes = file.info(outfiles)$size,
  stringsAsFactors = FALSE
)

write.csv(
  output_audit,
  file.path(OUTDIR, "04_output_audit_v8.csv"),
  row.names = FALSE
)

capture.output(
  sessionInfo(),
  file = file.path(OUTDIR, "sessionInfo_v8.txt")
)

cat("\n============================================================\n")
cat("Figure 4A GLOBAL network v8 completed.\n")
cat("Only visual strength was increased; displayed compartments and matrix logic are unchanged.\n")
cat("No CellChat analysis was rerun.\n\n")
cat("Open this PNG first:\n")
cat("  ", PNG_FILE, "\n", sep = "")
cat("\nPDF:\n")
cat("  ", PDF_FILE, "\n", sep = "")
cat("\nSVG:\n")
cat("  ", SVG_FILE, "\n", sep = "")
cat("============================================================\n")
