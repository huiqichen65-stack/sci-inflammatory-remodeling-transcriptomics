# ==============================================================================
# 23_step5BR2_GSE162610_single_principal_root_branch_topology_audit.R
# ==============================================================================
# Project:
#   Adult SCI GSE162610 activated/DAM-like microglia fate-bifurcation project
#
# STEP 5B-R2
#   Single-principal-root orientation + corrected branch-topology audit
#
# WHY THIS STEP EXISTS
#   Step5B-R established a useful injury-response continuum:
#       1 dpi -> 3 dpi -> 7 dpi
#   with:
#       cell-level true-time Spearman ~0.65
#       sample-level Spearman ~0.79
#       monotonic 1 -> 3 -> 7 dpi medians
#
#   However:
#     1) 225 primary 1 dpi root cells mapped to 16 principal vertices;
#     2) four degree>=3 graph vertices were detected;
#     3) the previous branch audit reported downstream_components = 0 for all
#        four branch vertices, which is topologically suspicious.
#
# THIS SCRIPT DOES TWO THINGS ONLY:
#
#   A. Re-orient the EXISTING Step5B-R graph with ONE objective principal root:
#        root_pr_nodes = the principal vertex containing the largest number of
#        primary 1 dpi root cells.
#
#      In the current Step5B-R report this was Y_19, but the script detects it
#      from the stored object rather than hard-coding the answer.
#
#   B. Re-audit each degree>=3 branch point using an explicit mapping:
#
#        cell barcode
#          -> closest principal graph vertex
#          -> graph component after removing branch vertex
#
#      This fixes the previous ambiguity where vertex-indexed names could be
#      mistakenly re-indexed by cell barcode.
#
# IMPORTANT:
#   - The principal graph is NOT relearned.
#   - UMAP is NOT recalculated.
#   - No Monocle clustering is rerun.
#   - No repair/scar score is used to orient pseudotime or choose a branchpoint.
#   - Repair/scar programs are projected only AFTER objective branch topology
#     has been established.
#
# PRIMARY BIFURCATION CRITERION:
#   A branch point must yield at least TWO downstream cell-supported components,
#   each with:
#       >=100 total cells
#       >=100 cells at 7 dpi
#
# BRANCHPOINT RANKING AMONG VALID CANDIDATES:
#   1) greater nearby high-confidence 3 dpi transition support;
#   2) larger minimum 7 dpi support across the two largest downstream arms;
#   3) smaller graph distance from the single root anchor.
#
# INPUTS:
#   Step5B-R Monocle3 CDS
#   Step5B-R full Seurat result object
#
# OUTPUT:
#   04_results/step5BR2_GSE162610_single_root_branch_audit/
# ==============================================================================

options(stringsAsFactors = FALSE)
set.seed(20260731)

# ==============================================================================
# 0. PATHS / FIXED DEFINITIONS
# ==============================================================================

PROJECT_ROOT <- Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
STEP5BR_DIR <- file.path(
  PROJECT_ROOT,
  "04_results",
  "step5BR_GSE162610_injury_response_trajectory_refinement"
)

INPUT_CDS <- file.path(
  STEP5BR_DIR,
  "GSE162610_step5BR_injury_monocle3_cds.rds"
)

INPUT_SEURAT <- file.path(
  STEP5BR_DIR,
  "GSE162610_step5BR_injury_trajectory_seurat.rds"
)

OUTDIR <- file.path(
  PROJECT_ROOT,
  "04_results",
  "step5BR2_GSE162610_single_root_branch_audit"
)
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

OUT_FIG1 <- file.path(
  OUTDIR,
  "Fig5BR2_01_single_principal_root_pseudotime.pdf"
)

OUT_FIG2 <- file.path(
  OUTDIR,
  "Fig5BR2_02_corrected_branch_topology_audit.pdf"
)

OUT_FIG3 <- file.path(
  OUTDIR,
  "Fig5BR2_03_objective_branch_arms.pdf"
)

OUT_FIG4 <- file.path(
  OUTDIR,
  "Fig5BR2_04_posthoc_repair_scar_by_arm.pdf"
)

OUT_FIG5 <- file.path(
  OUTDIR,
  "Fig5BR2_05_old_vs_single_root_orientation.pdf"
)

OUT_XLSX <- file.path(
  OUTDIR,
  "Table5BR2_single_root_branch_audit.xlsx"
)

OUT_CDS <- file.path(
  OUTDIR,
  "GSE162610_step5BR2_single_root_cds.rds"
)

OUT_RDS <- file.path(
  OUTDIR,
  "GSE162610_step5BR2_single_root_seurat.rds"
)

OUT_REP <- file.path(
  OUTDIR,
  "Step5BR2_report.txt"
)

OUT_SESS <- file.path(
  OUTDIR,
  "Step5BR2_sessionInfo.txt"
)

TIME_COL    <- "timepoint"
SAMPLE_COL  <- "sample_id"
CLUSTER_COL <- "step3A_cluster_r04"

PRIMARY_ROOT_COL <- "step5BR_primary_root"

REPAIR_SCORE_COL <- "repair_candidate_score"
SCAR_SCORE_COL   <- "scar_candidate_score"

MIN_TOTAL_CELLS_PER_DOWNSTREAM_ARM <- 100
MIN_7DPI_CELLS_PER_DOWNSTREAM_ARM  <- 100

BRANCHPOINT_GRAPH_RADIUS <- 2

INJURY_TIME_LEVELS <- c(
  "1 dpi",
  "3 dpi",
  "7 dpi"
)

TRUE_DAY <- c(
  "1 dpi" = 1,
  "3 dpi" = 3,
  "7 dpi" = 7
)

# ==============================================================================
# 1. PACKAGES
# ==============================================================================

required_pkgs <- c(
  "Seurat",
  "SeuratObject",
  "monocle3",
  "SingleCellExperiment",
  "igraph",
  "dplyr",
  "tidyr",
  "tibble",
  "ggplot2",
  "patchwork",
  "openxlsx"
)

missing_pkgs <- required_pkgs[
  !vapply(
    required_pkgs,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_pkgs) > 0) {
  stop(
    paste0(
      "\n缺少以下R包：\n  ",
      paste(missing_pkgs, collapse = ", "),
      "\n请在RStudio中安装后重新source。"
    ),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(monocle3)
  library(SingleCellExperiment)
  library(igraph)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
})

# ==============================================================================
# 2. COLORS
# ==============================================================================

TIME_COLORS <- c(
  "1 dpi" = "#FC8D62",
  "3 dpi" = "#8DA0CB",
  "7 dpi" = "#E78AC3"
)

ROOT_COL <- "#D73027"
BG_COL <- "#D9D9D9"

# ==============================================================================
# 3. HELPERS
# ==============================================================================

msg <- function(...) {
  cat(
    sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")),
    ...,
    "\n",
    sep = ""
  )
}

capture_table <- function(x) {
  paste(
    capture.output(
      print(
        tibble::as_tibble(x),
        n = Inf,
        width = Inf
      )
    ),
    collapse = "\n"
  )
}

get_pseudotime_safe <- function(cds) {

  pt <- tryCatch(
    monocle3::pseudotime(
      cds,
      reduction_method = "UMAP"
    ),
    error = function(e) NULL
  )

  if (!is.null(pt)) return(pt)

  aux <- cds@principal_graph_aux[["UMAP"]]

  if (
    !is.null(aux) &&
    "pseudotime" %in% names(aux)
  ) {
    return(aux$pseudotime)
  }

  stop(
    "Could not extract Monocle3 pseudotime.",
    call. = FALSE
  )
}

get_closest_vertex <- function(cds, pg) {

  aux <- cds@principal_graph_aux[["UMAP"]]
  raw <- aux$pr_graph_cell_proj_closest_vertex

  if (is.null(raw)) {
    stop(
      "Monocle3 closest-principal-vertex mapping not found.",
      call. = FALSE
    )
  }

  if (is.matrix(raw) || is.data.frame(raw)) {
    vv <- raw[, 1]
    names(vv) <- rownames(raw)
  } else {
    vv <- raw
  }

  if (is.null(names(vv))) {
    names(vv) <- colnames(cds)
  }

  vv_chr <- as.character(vv)
  pg_names <- igraph::V(pg)$name

  if (all(vv_chr %in% pg_names)) {
    return(
      stats::setNames(
        vv_chr,
        names(vv)
      )
    )
  }

  y_names <- paste0(
    "Y_",
    vv_chr
  )

  if (all(y_names %in% pg_names)) {
    return(
      stats::setNames(
        y_names,
        names(vv)
      )
    )
  }

  idx <- suppressWarnings(
    as.integer(vv_chr)
  )

  if (
    all(!is.na(idx)) &&
    all(idx >= 1) &&
    all(idx <= igraph::vcount(pg))
  ) {
    return(
      stats::setNames(
        pg_names[idx],
        names(vv)
      )
    )
  }

  stop(
    "Could not map closest-principal-vertex IDs to graph vertex names.",
    call. = FALSE
  )
}

cells_near_graph_vertex <- function(
    pg,
    closest_vertex,
    vertex_name,
    radius = 2) {

  dd <- igraph::distances(
    pg,
    v = vertex_name,
    to = igraph::V(pg),
    mode = "all"
  )

  near_vertices <- colnames(dd)[
    as.numeric(dd[1, ]) <= radius
  ]

  names(closest_vertex)[
    closest_vertex %in%
      near_vertices
  ]
}

safe_add_sheet <- function(
    wb,
    sheet_name,
    df) {

  openxlsx::addWorksheet(
    wb,
    sheetName = sheet_name
  )

  ddf <- tryCatch(
    as.data.frame(df),
    error = function(e) NULL
  )

  if (
    is.null(ddf) ||
    nrow(ddf) == 0 ||
    ncol(ddf) == 0
  ) {

    openxlsx::writeData(
      wb,
      sheet = sheet_name,
      x = data.frame(
        note =
          "No results generated under the predefined Step5B-R2 criteria.",
        stringsAsFactors = FALSE
      )
    )

  } else {

    openxlsx::writeDataTable(
      wb,
      sheet = sheet_name,
      x = ddf,
      tableStyle = "TableStyleMedium2"
    )
  }

  openxlsx::freezePane(
    wb,
    sheet = sheet_name,
    firstRow = TRUE
  )
}

# Correct cell -> graph component mapper.
#
# IMPORTANT:
# membership is named by PRINCIPAL VERTEX.
# closest_vertex is named by CELL.
#
# We explicitly preserve cell names after indexing membership by each cell's
# closest principal vertex.
map_cells_to_components <- function(
    membership,
    closest_vertex,
    cell_names) {

  cv <- closest_vertex[cell_names]

  mapped <- membership[
    as.character(cv)
  ]

  out <- as.integer(mapped)

  names(out) <- cell_names

  out
}

# ==============================================================================
# 4. READ EXISTING STEP5B-R GRAPH / SEURAT OBJECT
# ==============================================================================

msg("====================================================")
msg("GSE162610 Step 5B-R2 | single principal root + branch audit")
msg("====================================================")

if (!file.exists(INPUT_CDS)) {
  stop(
    paste0(
      "找不到Step5B-R CDS：\n",
      INPUT_CDS
    ),
    call. = FALSE
  )
}

if (!file.exists(INPUT_SEURAT)) {
  stop(
    paste0(
      "找不到Step5B-R Seurat结果：\n",
      INPUT_SEURAT
    ),
    call. = FALSE
  )
}

cds <- readRDS(INPUT_CDS)
obj <- readRDS(INPUT_SEURAT)

if (!inherits(obj, "Seurat")) {
  stop(
    "Step5B-R Seurat input is invalid.",
    call. = FALSE
  )
}

pg <- monocle3::principal_graph(
  cds
)[["UMAP"]]

if (is.null(pg)) {
  stop(
    "Existing Step5B-R CDS does not contain a UMAP principal graph.",
    call. = FALSE
  )
}

msg(
  "Existing injury graph: ",
  igraph::vcount(pg),
  " vertices / ",
  igraph::ecount(pg),
  " edges."
)

msg(
  "Graph connected components: ",
  igraph::components(pg)$no
)

# ==============================================================================
# 5. EXTRACT PREVIOUS ROOT / PSEUDOTIME / CELL METADATA
# ==============================================================================

required_meta <- c(
  TIME_COL,
  SAMPLE_COL,
  CLUSTER_COL,
  PRIMARY_ROOT_COL,
  "shared_transition_3dpi_highconf",
  "transition_status",
  REPAIR_SCORE_COL,
  SCAR_SCORE_COL
)

missing_meta <- required_meta[
  !required_meta %in%
    colnames(obj[[]])
]

if (length(missing_meta) > 0) {
  stop(
    paste0(
      "Step5B-R Seurat结果缺少metadata：",
      paste(missing_meta, collapse = ", ")
    ),
    call. = FALSE
  )
}

cds_cells <- colnames(cds)

if (!all(cds_cells %in% colnames(obj))) {
  stop(
    "Some CDS cells are absent from Step5B-R Seurat object.",
    call. = FALSE
  )
}

primary_root_cells <- colnames(obj)[
  as.logical(
    obj[[PRIMARY_ROOT_COL, drop = TRUE]]
  )
]

primary_root_cells <- intersect(
  primary_root_cells,
  cds_cells
)

if (length(primary_root_cells) == 0) {
  stop(
    "No Step5B-R primary root cells found in CDS.",
    call. = FALSE
  )
}

old_pt <- NULL

if (
  "pseudotime_primary" %in%
    colnames(
      SummarizedExperiment::colData(cds)
    )
) {
  old_pt <- as.numeric(
    SummarizedExperiment::colData(cds)$pseudotime_primary
  )
  names(old_pt) <- cds_cells
} else if (
  "step5BR_pseudotime" %in%
    colnames(obj[[]])
) {
  old_pt <- as.numeric(
    obj$step5BR_pseudotime[
      match(
        cds_cells,
        colnames(obj)
      )
    ]
  )
  names(old_pt) <- cds_cells
}

# ==============================================================================
# 6. CLOSEST PRINCIPAL VERTEX + SINGLE ROOT ANCHOR
# ==============================================================================

closest_vertex <- get_closest_vertex(
  cds,
  pg
)

closest_vertex <- closest_vertex[
  cds_cells
]

root_vertex_counts <- tibble::tibble(
  cell = primary_root_cells,
  principal_vertex =
    as.character(
      closest_vertex[
        primary_root_cells
      ]
    )
) |>
  dplyr::count(
    principal_vertex,
    name = "n_primary_root_cells"
  ) |>
  dplyr::arrange(
    dplyr::desc(
      n_primary_root_cells
    )
  )

if (nrow(root_vertex_counts) == 0) {
  stop(
    "Primary root cells could not be mapped to principal vertices.",
    call. = FALSE
  )
}

ROOT_ANCHOR_VERTEX <-
  root_vertex_counts$principal_vertex[1]

msg(
  "Primary root cells: ",
  length(primary_root_cells)
)

msg(
  "Primary root cells map to ",
  nrow(root_vertex_counts),
  " principal vertices."
)

msg(
  "Single objective root principal node: ",
  ROOT_ANCHOR_VERTEX,
  " (",
  root_vertex_counts$n_primary_root_cells[1],
  " primary root cells)"
)

# ==============================================================================
# 7. RE-ORDER CELLS ON SAME GRAPH USING ONE root_pr_node
# ==============================================================================

msg(
  "Re-ordering SAME principal graph with root_pr_nodes = ",
  ROOT_ANCHOR_VERTEX
)

cds_single <- monocle3::order_cells(
  cds,
  reduction_method = "UMAP",
  root_pr_nodes = ROOT_ANCHOR_VERTEX
)

pt_single <- get_pseudotime_safe(
  cds_single
)

pt_single <- pt_single[
  cds_cells
]

finite_fraction <- mean(
  is.finite(
    pt_single
  )
)

msg(
  "Finite single-root pseudotime fraction: ",
  round(
    100 *
      finite_fraction,
    2
  ),
  "%"
)

# ==============================================================================
# 8. SINGLE-ROOT TRUE-TIME VALIDATION
# ==============================================================================

cell_df <- tibble::tibble(
  cell = cds_cells,

  timepoint =
    as.character(
      obj[[TIME_COL, drop = TRUE]][
        match(
          cds_cells,
          colnames(obj)
        )
      ]
    ),

  sample =
    as.character(
      obj[[SAMPLE_COL, drop = TRUE]][
        match(
          cds_cells,
          colnames(obj)
        )
      ]
    ),

  cluster =
    as.character(
      obj[[CLUSTER_COL, drop = TRUE]][
        match(
          cds_cells,
          colnames(obj)
        )
      ]
    ),

  shared_transition_3dpi_highconf =
    as.logical(
      obj$shared_transition_3dpi_highconf[
        match(
          cds_cells,
          colnames(obj)
        )
      ]
    ),

  transition_status =
    as.character(
      obj$transition_status[
        match(
          cds_cells,
          colnames(obj)
        )
      ]
    ),

  repair_candidate_score =
    as.numeric(
      obj[[REPAIR_SCORE_COL, drop = TRUE]][
        match(
          cds_cells,
          colnames(obj)
        )
      ]
    ),

  scar_candidate_score =
    as.numeric(
      obj[[SCAR_SCORE_COL, drop = TRUE]][
        match(
          cds_cells,
          colnames(obj)
        )
      ]
    ),

  old_pseudotime =
    if (
      is.null(old_pt)
    ) {
      NA_real_
    } else {
      as.numeric(
        old_pt[
          cds_cells
        ]
      )
    },

  pseudotime_single_root =
    as.numeric(
      pt_single
    ),

  principal_vertex =
    as.character(
      closest_vertex[
        cds_cells
      ]
    ),

  primary_root_cell =
    cds_cells %in%
      primary_root_cells
)

cell_df$repair_minus_scar <-
  cell_df$repair_candidate_score -
  cell_df$scar_candidate_score

cell_df$true_day <-
  unname(
    TRUE_DAY[
      cell_df$timepoint
    ]
  )

finite_idx <- is.finite(
  cell_df$pseudotime_single_root
)

cell_level_time_spearman <-
  suppressWarnings(
    stats::cor(
      cell_df$true_day[
        finite_idx
      ],
      cell_df$pseudotime_single_root[
        finite_idx
      ],
      method = "spearman"
    )
  )

time_summary <- cell_df |>
  dplyr::filter(
    is.finite(
      pseudotime_single_root
    )
  ) |>
  dplyr::mutate(
    timepoint =
      factor(
        timepoint,
        levels =
          INJURY_TIME_LEVELS,
        ordered = TRUE
      )
  ) |>
  dplyr::group_by(
    timepoint
  ) |>
  dplyr::summarise(
    n_cells = dplyr::n(),
    pseudotime_mean =
      mean(
        pseudotime_single_root
      ),
    pseudotime_median =
      median(
        pseudotime_single_root
      ),
    pseudotime_q25 =
      stats::quantile(
        pseudotime_single_root,
        0.25
      ),
    pseudotime_q75 =
      stats::quantile(
        pseudotime_single_root,
        0.75
      ),
    .groups = "drop"
  )

time_med <- time_summary$pseudotime_median
names(time_med) <-
  as.character(
    time_summary$timepoint
  )

time_monotonic <- all(
  diff(
    time_med[
      INJURY_TIME_LEVELS
    ]
  ) > 0
)

sample_summary <- cell_df |>
  dplyr::filter(
    is.finite(
      pseudotime_single_root
    )
  ) |>
  dplyr::group_by(
    timepoint,
    sample
  ) |>
  dplyr::summarise(
    n_cells = dplyr::n(),
    median_pseudotime =
      median(
        pseudotime_single_root
      ),
    mean_pseudotime =
      mean(
        pseudotime_single_root
      ),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    true_day =
      unname(
        TRUE_DAY[
          timepoint
        ]
      )
  )

sample_level_time_spearman <-
  suppressWarnings(
    stats::cor(
      sample_summary$true_day,
      sample_summary$median_pseudotime,
      method = "spearman"
    )
  )

reference_df <- cell_df |>
  dplyr::mutate(
    reference_state =
      dplyr::case_when(
        primary_root_cell ~
          "1dpi_primary_root_cells",

        timepoint == "1 dpi" &
          cluster %in%
            c("2", "5", "10") ~
          "1dpi_early_activated",

        timepoint == "3 dpi" &
          shared_transition_3dpi_highconf ~
          "3dpi_shared_transition",

        timepoint == "7 dpi" &
          cluster == "1" &
          transition_status ==
            "Mature_DAM_like" ~
          "7dpi_established_DAM",

        TRUE ~
          "Other"
      )
  )

reference_summary <- reference_df |>
  dplyr::filter(
    reference_state != "Other",
    is.finite(
      pseudotime_single_root
    )
  ) |>
  dplyr::group_by(
    reference_state
  ) |>
  dplyr::summarise(
    n_cells = dplyr::n(),
    pseudotime_mean =
      mean(
        pseudotime_single_root
      ),
    pseudotime_median =
      median(
        pseudotime_single_root
      ),
    pseudotime_q25 =
      stats::quantile(
        pseudotime_single_root,
        0.25
      ),
    pseudotime_q75 =
      stats::quantile(
        pseudotime_single_root,
        0.75
      ),
    .groups = "drop"
  )

# Compare previous cell-root orientation with single-principal-node orientation.
old_vs_single_spearman <- NA_real_

if (!all(is.na(cell_df$old_pseudotime))) {

  ok <- is.finite(
    cell_df$old_pseudotime
  ) &
    is.finite(
      cell_df$pseudotime_single_root
    )

  old_vs_single_spearman <-
    suppressWarnings(
      stats::cor(
        cell_df$old_pseudotime[
          ok
        ],
        cell_df$pseudotime_single_root[
          ok
        ],
        method = "spearman"
      )
    )
}

msg(
  "Cell-level true-time Spearman: ",
  round(
    cell_level_time_spearman,
    4
  )
)

msg(
  "Sample-level true-time Spearman: ",
  round(
    sample_level_time_spearman,
    4
  )
)

msg(
  "1 -> 3 -> 7 dpi medians monotonic: ",
  time_monotonic
)

msg(
  "Old cell-root vs single-node pseudotime Spearman: ",
  round(
    old_vs_single_spearman,
    4
  )
)

# ==============================================================================
# 9. PURE GRAPH TOPOLOGY
# ==============================================================================

degree_vec <- igraph::degree(
  pg
)

branch_vertices <- names(
  degree_vec[
    degree_vec >= 3
  ]
)

leaf_vertices <- names(
  degree_vec[
    degree_vec == 1
  ]
)

root_dist <- igraph::distances(
  pg,
  v = ROOT_ANCHOR_VERTEX,
  to = igraph::V(pg)
)

graph_vertex_table <- tibble::tibble(
  vertex = names(
    degree_vec
  ),
  degree =
    as.integer(
      degree_vec
    ),
  is_root_anchor =
    names(degree_vec) ==
      ROOT_ANCHOR_VERTEX,
  is_branch =
    names(degree_vec) %in%
      branch_vertices,
  is_leaf =
    names(degree_vec) %in%
      leaf_vertices,
  graph_distance_from_root =
    as.numeric(
      root_dist[
        1,
        names(
          degree_vec
        )
      ]
    )
)

msg(
  "Degree>=3 branch vertices: ",
  length(
    branch_vertices
  )
)

# ==============================================================================
# 10. CORRECTED BRANCH COMPONENT AUDIT
# ==============================================================================

branchpoint_rows <- list()
branch_component_rows <- list()
branch_cache <- list()

if (length(branch_vertices) > 0) {

  for (bp in branch_vertices) {

    # If branchpoint itself is the root anchor, there is no unique upstream side.
    if (bp == ROOT_ANCHOR_VERTEX) {

      branchpoint_rows[[bp]] <-
        tibble::tibble(
          branch_vertex = bp,
          degree =
            as.integer(
              igraph::degree(
                pg,
                v = bp
              )
            ),
          topology_components_after_deletion =
            NA_integer_,
          topology_downstream_components =
            NA_integer_,
          cell_supported_downstream_components =
            NA_integer_,
          downstream_major_arms =
            NA_integer_,
          min_top2_7dpi_cells =
            NA_integer_,
          top2_7dpi_sum =
            NA_integer_,
          near_cells =
            NA_integer_,
          near_3dpi_cells =
            NA_integer_,
          near_transition_cells =
            NA_integer_,
          near_transition_fraction =
            NA_real_,
          note =
            "Branch vertex equals root anchor; not eligible as downstream bifurcation."
        )

      next
    }

    g2 <- igraph::delete_vertices(
      pg,
      bp
    )

    comp <- igraph::components(
      g2
    )

    membership <- comp$membership
    names(membership) <-
      igraph::V(g2)$name

    topology_components <-
      comp$no

    root_side_component <-
      membership[
        ROOT_ANCHOR_VERTEX
      ]

    if (
      length(root_side_component) == 0 ||
      is.na(root_side_component)
    ) {

      branchpoint_rows[[bp]] <-
        tibble::tibble(
          branch_vertex = bp,
          degree =
            as.integer(
              igraph::degree(
                pg,
                v = bp
              )
            ),
          topology_components_after_deletion =
            topology_components,
          topology_downstream_components =
            NA_integer_,
          cell_supported_downstream_components =
            NA_integer_,
          downstream_major_arms =
            NA_integer_,
          min_top2_7dpi_cells =
            NA_integer_,
          top2_7dpi_sum =
            NA_integer_,
          near_cells =
            NA_integer_,
          near_3dpi_cells =
            NA_integer_,
          near_transition_cells =
            NA_integer_,
          near_transition_fraction =
            NA_real_,
          note =
            "Root anchor absent after branch deletion; unexpected graph condition."
        )

      next
    }

    topology_downstream_components <-
      topology_components - 1L

    # --------------------------------------------------------------------------
    # CRITICAL FIX:
    # Explicit CELL-NAMED component vector.
    # --------------------------------------------------------------------------
    cell_component <- map_cells_to_components(
      membership = membership,
      closest_vertex = closest_vertex,
      cell_names = cell_df$cell
    )

    tmp <- cell_df |>
      dplyr::mutate(
        graph_component =
          as.integer(
            cell_component[
              cell
            ]
          )
      )

    # Cells whose closest vertex is exactly the removed BP naturally become NA.
    # Those are recorded separately and are NOT silently assigned downstream.
    n_cells_at_bp <-
      sum(
        tmp$principal_vertex ==
          bp
      )

    component_summary <- tmp |>
      dplyr::filter(
        !is.na(
          graph_component
        )
      ) |>
      dplyr::group_by(
        graph_component
      ) |>
      dplyr::summarise(
        n_cells = dplyr::n(),
        n_1dpi =
          sum(
            timepoint ==
              "1 dpi"
          ),
        n_3dpi =
          sum(
            timepoint ==
              "3 dpi"
          ),
        n_7dpi =
          sum(
            timepoint ==
              "7 dpi"
          ),
        n_transition =
          sum(
            shared_transition_3dpi_highconf,
            na.rm = TRUE
          ),
        pseudotime_median =
          median(
            pseudotime_single_root,
            na.rm = TRUE
          ),
        .groups = "drop"
      ) |>
      dplyr::mutate(
        branch_vertex = bp,
        is_root_side =
          graph_component ==
            root_side_component
      )

    downstream <- component_summary |>
      dplyr::filter(
        !is_root_side
      ) |>
      dplyr::arrange(
        dplyr::desc(
          n_7dpi
        ),
        dplyr::desc(
          n_cells
        )
      )

    major <- downstream |>
      dplyr::filter(
        n_cells >=
          MIN_TOTAL_CELLS_PER_DOWNSTREAM_ARM,
        n_7dpi >=
          MIN_7DPI_CELLS_PER_DOWNSTREAM_ARM
      )

    top2 <- head(
      major$n_7dpi,
      2
    )

    near_cells <- cells_near_graph_vertex(
      pg,
      closest_vertex,
      bp,
      radius =
        BRANCHPOINT_GRAPH_RADIUS
    )

    near_df <- cell_df |>
      dplyr::filter(
        cell %in%
          near_cells
      )

    branchpoint_rows[[bp]] <-
      tibble::tibble(
        branch_vertex = bp,

        degree =
          as.integer(
            igraph::degree(
              pg,
              v = bp
            )
          ),

        topology_components_after_deletion =
          topology_components,

        topology_downstream_components =
          topology_downstream_components,

        cell_supported_downstream_components =
          nrow(
            downstream
          ),

        downstream_major_arms =
          nrow(
            major
          ),

        min_top2_7dpi_cells =
          ifelse(
            length(top2) >= 2,
            min(top2),
            0
          ),

        top2_7dpi_sum =
          sum(top2),

        cells_mapped_exactly_to_branch_vertex =
          n_cells_at_bp,

        near_cells =
          nrow(
            near_df
          ),

        near_3dpi_cells =
          sum(
            near_df$timepoint ==
              "3 dpi"
          ),

        near_transition_cells =
          sum(
            near_df$shared_transition_3dpi_highconf,
            na.rm = TRUE
          ),

        near_transition_fraction =
          ifelse(
            nrow(near_df) > 0,
            mean(
              near_df$shared_transition_3dpi_highconf,
              na.rm = TRUE
            ),
            NA_real_
          ),

        note = "OK"
      )

    if (nrow(component_summary) > 0) {

      branch_component_rows[[bp]] <-
        component_summary |>
        dplyr::select(
          branch_vertex,
          dplyr::everything()
        )
    }

    branch_cache[[bp]] <- list(
      membership =
        membership,
      root_side_component =
        as.integer(
          root_side_component
        ),
      component_summary =
        component_summary,
      downstream =
        downstream,
      major =
        major,
      cell_component =
        cell_component
    )
  }
}

branchpoint_summary <- if (
  length(
    branchpoint_rows
  ) > 0
) {
  dplyr::bind_rows(
    branchpoint_rows
  ) |>
    dplyr::left_join(
      graph_vertex_table |>
        dplyr::select(
          vertex,
          graph_distance_from_root
        ),
      by = c(
        "branch_vertex" =
          "vertex"
      )
    ) |>
    dplyr::mutate(
      candidate_bifurcation =
        !is.na(
          downstream_major_arms
        ) &
        downstream_major_arms >= 2
    ) |>
    dplyr::arrange(
      dplyr::desc(
        candidate_bifurcation
      ),
      dplyr::desc(
        near_transition_cells
      ),
      dplyr::desc(
        min_top2_7dpi_cells
      ),
      graph_distance_from_root
    )
} else {
  tibble::tibble()
}

branch_component_summary <- if (
  length(
    branch_component_rows
  ) > 0
) {
  dplyr::bind_rows(
    branch_component_rows
  )
} else {
  tibble::tibble()
}

primary_branchpoint <- NA_character_

if (
  nrow(
    branchpoint_summary
  ) > 0 &&
  any(
    branchpoint_summary$candidate_bifurcation,
    na.rm = TRUE
  )
) {

  candidate_rows <- branchpoint_summary |>
    dplyr::filter(
      candidate_bifurcation
    ) |>
    dplyr::arrange(
      dplyr::desc(
        near_transition_cells
      ),
      dplyr::desc(
        min_top2_7dpi_cells
      ),
      graph_distance_from_root
    )

  primary_branchpoint <-
    candidate_rows$branch_vertex[1]

  msg(
    "Objective primary bifurcation candidate: ",
    primary_branchpoint
  )

} else {

  msg(
    "No branch point satisfies the predefined two-major-7dpi-arm criterion."
  )
}

# ==============================================================================
# 11. FORMAL ARM ASSIGNMENT IF CANDIDATE EXISTS
# ==============================================================================

cell_df$branch_arm <- NA_character_
cell_df$major_late_arm <- FALSE

arm_summary <- tibble::tibble()
branch_program_summary <- tibble::tibble()
branch_program_sample_summary <- tibble::tibble()

if (!is.na(primary_branchpoint)) {

  cache <-
    branch_cache[[
      primary_branchpoint
    ]]

  root_side_component <-
    cache$root_side_component

  cell_component <-
    cache$cell_component

  downstream <- cache$downstream |>
    dplyr::arrange(
      dplyr::desc(
        n_7dpi
      ),
      dplyr::desc(
        n_cells
      )
    ) |>
    dplyr::mutate(
      arm_rank =
        dplyr::row_number(),
      branch_arm =
        paste0(
          "Arm_",
          arm_rank
        ),
      major_late_arm =
        n_cells >=
          MIN_TOTAL_CELLS_PER_DOWNSTREAM_ARM &
        n_7dpi >=
          MIN_7DPI_CELLS_PER_DOWNSTREAM_ARM
    )

  arm_summary <- downstream

  arm_map <- stats::setNames(
    downstream$branch_arm,
    downstream$graph_component
  )

  major_map <- stats::setNames(
    downstream$major_late_arm,
    downstream$graph_component
  )

  cell_df$branch_component <-
    as.integer(
      cell_component[
        cell_df$cell
      ]
    )

  root_side_idx <- which(
    !is.na(
      cell_df$branch_component
    ) &
      cell_df$branch_component ==
        root_side_component
  )

  cell_df$branch_arm[
    root_side_idx
  ] <- "Pre_branch_root_side"

  downstream_idx <- which(
    !is.na(
      cell_df$branch_component
    ) &
      cell_df$branch_component !=
        root_side_component
  )

  downstream_comp_chr <-
    as.character(
      cell_df$branch_component[
        downstream_idx
      ]
    )

  cell_df$branch_arm[
    downstream_idx
  ] <- unname(
    arm_map[
      downstream_comp_chr
    ]
  )

  cell_df$major_late_arm[
    downstream_idx
  ] <- unname(
    major_map[
      downstream_comp_chr
    ]
  )

  cell_df$branch_arm[
    cell_df$principal_vertex ==
      primary_branchpoint
  ] <- "Branchpoint_vertex"

  # Post-hoc program summaries ONLY for 7 dpi arm cells.
  branch_program_summary <- cell_df |>
    dplyr::filter(
      timepoint == "7 dpi",
      !is.na(
        branch_arm
      ),
      grepl(
        "^Arm_",
        branch_arm
      )
    ) |>
    dplyr::group_by(
      branch_arm
    ) |>
    dplyr::summarise(
      n_cells = dplyr::n(),

      repair_mean =
        mean(
          repair_candidate_score,
          na.rm = TRUE
        ),

      repair_median =
        median(
          repair_candidate_score,
          na.rm = TRUE
        ),

      scar_mean =
        mean(
          scar_candidate_score,
          na.rm = TRUE
        ),

      scar_median =
        median(
          scar_candidate_score,
          na.rm = TRUE
        ),

      repair_minus_scar_mean =
        mean(
          repair_minus_scar,
          na.rm = TRUE
        ),

      repair_minus_scar_median =
        median(
          repair_minus_scar,
          na.rm = TRUE
        ),

      pseudotime_median =
        median(
          pseudotime_single_root,
          na.rm = TRUE
        ),

      .groups = "drop"
    )

  branch_program_sample_summary <- cell_df |>
    dplyr::filter(
      timepoint == "7 dpi",
      !is.na(
        branch_arm
      ),
      grepl(
        "^Arm_",
        branch_arm
      )
    ) |>
    dplyr::group_by(
      sample,
      branch_arm
    ) |>
    dplyr::summarise(
      n_cells = dplyr::n(),

      repair_median =
        median(
          repair_candidate_score,
          na.rm = TRUE
        ),

      scar_median =
        median(
          scar_candidate_score,
          na.rm = TRUE
        ),

      repair_minus_scar_median =
        median(
          repair_minus_scar,
          na.rm = TRUE
        ),

      .groups = "drop"
    )
}

# ==============================================================================
# 12. ADD RESULTS TO CDS / SEURAT
# ==============================================================================

SummarizedExperiment::colData(
  cds_single
)$step5BR2_pseudotime <-
  as.numeric(
    pt_single
  )

SummarizedExperiment::colData(
  cds_single
)$step5BR2_branch_arm <-
  cell_df$branch_arm

obj$step5BR2_single_root_pseudotime <-
  NA_real_

obj$step5BR2_principal_root_node <-
  ROOT_ANCHOR_VERTEX

obj$step5BR2_branch_arm <-
  NA_character_

obj$step5BR2_major_late_arm <-
  FALSE

idx_obj <- match(
  cell_df$cell,
  colnames(obj)
)

obj$step5BR2_single_root_pseudotime[
  idx_obj
] <- cell_df$pseudotime_single_root

obj$step5BR2_branch_arm[
  idx_obj
] <- cell_df$branch_arm

obj$step5BR2_major_late_arm[
  idx_obj
] <- cell_df$major_late_arm

obj$step5BR2_primary_branchpoint <-
  ifelse(
    is.na(
      primary_branchpoint
    ),
    NA_character_,
    primary_branchpoint
  )

# ==============================================================================
# 13. FIGURE 1 -- SINGLE ROOT PSEUDOTIME
# ==============================================================================

msg(
  "Creating Fig5BR2_01_single_principal_root_pseudotime.pdf..."
)

p1a <- monocle3::plot_cells(
  cds_single,
  reduction_method = "UMAP",
  color_cells_by = TIME_COL,
  show_trajectory_graph = TRUE,
  label_cell_groups = FALSE,
  label_leaves = TRUE,
  label_branch_points = TRUE,
  graph_label_size = 3,
  cell_size = 0.45
) +
  ggplot2::scale_color_manual(
    values = TIME_COLORS
  ) +
  ggplot2::ggtitle(
    paste0(
      "Same graph oriented from single principal root ",
      ROOT_ANCHOR_VERTEX
    )
  )

p1b <- monocle3::plot_cells(
  cds_single,
  reduction_method = "UMAP",
  color_cells_by = "pseudotime",
  show_trajectory_graph = TRUE,
  label_cell_groups = FALSE,
  label_leaves = TRUE,
  label_branch_points = TRUE,
  graph_label_size = 3,
  cell_size = 0.45
) +
  ggplot2::ggtitle(
    "Single-principal-node pseudotime"
  )

pt_plot <- cell_df |>
  dplyr::filter(
    is.finite(
      pseudotime_single_root
    )
  ) |>
  dplyr::mutate(
    timepoint =
      factor(
        timepoint,
        levels =
          INJURY_TIME_LEVELS,
        ordered = TRUE
      )
  )

p1c <- ggplot2::ggplot(
  pt_plot,
  ggplot2::aes(
    x = timepoint,
    y =
      pseudotime_single_root,
    fill = timepoint
  )
) +
  ggplot2::geom_violin(
    trim = TRUE,
    scale = "width",
    alpha = 0.75,
    linewidth = 0.25
  ) +
  ggplot2::geom_boxplot(
    width = 0.12,
    outlier.shape = NA,
    fill = "#FFFFFF",
    linewidth = 0.25
  ) +
  ggplot2::scale_fill_manual(
    values = TIME_COLORS,
    guide = "none"
  ) +
  ggplot2::labs(
    title = "Single-root pseudotime versus true injury time",
    subtitle = paste0(
      "Cell Spearman = ",
      round(
        cell_level_time_spearman,
        3
      ),
      " | Sample Spearman = ",
      round(
        sample_level_time_spearman,
        3
      ),
      " | Monotonic = ",
      time_monotonic
    ),
    x = "True time",
    y = "Pseudotime"
  ) +
  ggplot2::theme_classic(
    base_size = 11
  )

reference_order <- c(
  "1dpi_primary_root_cells",
  "1dpi_early_activated",
  "3dpi_shared_transition",
  "7dpi_established_DAM"
)

ref_plot <- reference_df |>
  dplyr::mutate(
    reference_state =
      factor(
        reference_state,
        levels =
          reference_order,
        ordered = TRUE
      )
  ) |>
  dplyr::filter(
    reference_state != "Other",
    is.finite(
      pseudotime_single_root
    )
  )

p1d <- ggplot2::ggplot(
  ref_plot,
  ggplot2::aes(
    x = reference_state,
    y =
      pseudotime_single_root
  )
) +
  ggplot2::geom_violin(
    fill = "#9ECAE1",
    trim = TRUE,
    scale = "width",
    linewidth = 0.25
  ) +
  ggplot2::geom_boxplot(
    width = 0.12,
    outlier.shape = NA,
    fill = "#FFFFFF",
    linewidth = 0.25
  ) +
  ggplot2::labs(
    title = "Formal biological references",
    x = NULL,
    y = "Single-root pseudotime"
  ) +
  ggplot2::theme_classic(
    base_size = 10
  ) +
  ggplot2::theme(
    axis.text.x =
      ggplot2::element_text(
        angle = 25,
        hjust = 1
      )
  )

grDevices::pdf(
  OUT_FIG1,
  width = 14,
  height = 10,
  useDingbats = FALSE
)

print(
  (p1a | p1b) /
    (p1c | p1d) +
    patchwork::plot_annotation(
      title =
        "Step 5B-R2 | Single-principal-root orientation"
    )
)

grDevices::dev.off()

# ==============================================================================
# 14. FIGURE 2 -- CORRECTED BRANCH TOPOLOGY AUDIT
# ==============================================================================

msg(
  "Creating Fig5BR2_02_corrected_branch_topology_audit.pdf..."
)

if (nrow(branchpoint_summary) > 0) {

  p2a_df <- branchpoint_summary |>
    dplyr::select(
      branch_vertex,
      topology_downstream_components,
      cell_supported_downstream_components
    ) |>
    tidyr::pivot_longer(
      cols = c(
        topology_downstream_components,
        cell_supported_downstream_components
      ),
      names_to = "component_type",
      values_to = "n_components"
    ) |>
    dplyr::mutate(
      component_type =
        dplyr::recode(
          component_type,
          topology_downstream_components =
            "Topology components",
          cell_supported_downstream_components =
            "Cell-supported components"
        )
    )

  p2a <- ggplot2::ggplot(
    p2a_df,
    ggplot2::aes(
      x = branch_vertex,
      y = n_components,
      fill = component_type
    )
  ) +
    ggplot2::geom_col(
      position = "dodge",
      width = 0.70
    ) +
    ggplot2::labs(
      title =
        "Branch deletion: topology versus actual cell support",
      x =
        "Principal-graph branch vertex",
      y =
        "Downstream components",
      fill = NULL
    ) +
    ggplot2::theme_classic(
      base_size = 11
    )

  p2b <- ggplot2::ggplot(
    branchpoint_summary,
    ggplot2::aes(
      x = branch_vertex,
      y = min_top2_7dpi_cells,
      fill = candidate_bifurcation
    )
  ) +
    ggplot2::geom_col(
      width = 0.70
    ) +
    ggplot2::geom_hline(
      yintercept =
        MIN_7DPI_CELLS_PER_DOWNSTREAM_ARM,
      linetype = 2,
      linewidth = 0.45,
      color = "#555555"
    ) +
    ggplot2::scale_fill_manual(
      values = c(
        "FALSE" = "#BDBDBD",
        "TRUE" = "#2C7BB6"
      ),
      name =
        "Objective\nbifurcation"
    ) +
    ggplot2::labs(
      title =
        "Minimum 7 dpi support across two largest major arms",
      subtitle = paste0(
        "Dashed line = ",
        MIN_7DPI_CELLS_PER_DOWNSTREAM_ARM,
        " cells/arm"
      ),
      x =
        "Principal-graph branch vertex",
      y =
        "Minimum top-two 7 dpi cells"
    ) +
    ggplot2::theme_classic(
      base_size = 11
    )

  p2c <- ggplot2::ggplot(
    branchpoint_summary,
    ggplot2::aes(
      x = graph_distance_from_root,
      y = near_transition_cells,
      label = branch_vertex,
      size = min_top2_7dpi_cells
    )
  ) +
    ggplot2::geom_point(
      alpha = 0.8
    ) +
    ggplot2::geom_text(
      vjust = -0.7,
      size = 3.5
    ) +
    ggplot2::labs(
      title =
        "Branch points relative to root and 3 dpi transition support",
      x =
        "Graph distance from single root anchor",
      y =
        "Nearby high-confidence 3 dpi transition cells",
      size =
        "Min top-two\n7 dpi cells"
    ) +
    ggplot2::theme_classic(
      base_size = 11
    )

  if (
    nrow(
      branch_component_summary
    ) > 0
  ) {

    comp_plot_df <-
      branch_component_summary |>
      dplyr::filter(
        !is_root_side
      )

    p2d <- ggplot2::ggplot(
      comp_plot_df,
      ggplot2::aes(
        x = factor(
          graph_component
        ),
        y = n_cells,
        fill = branch_vertex
      )
    ) +
      ggplot2::geom_col(
        position = "dodge"
      ) +
      ggplot2::facet_wrap(
        ~branch_vertex,
        scales = "free_x"
      ) +
      ggplot2::labs(
        title =
          "Actual cell-supported downstream components",
        x =
          "Component after branch-vertex deletion",
        y =
          "Cells"
      ) +
      ggplot2::theme_classic(
        base_size = 9
      ) +
      ggplot2::theme(
        legend.position = "none"
      )

  } else {

    p2d <- ggplot2::ggplot(
      tibble::tibble(
        x = 1,
        y = 1
      ),
      ggplot2::aes(x, y)
    ) +
      ggplot2::annotate(
        "text",
        x = 1,
        y = 1,
        label =
          "No cell-supported downstream components.",
        size = 5
      ) +
      ggplot2::theme_void()
  }

} else {

  placeholder <- tibble::tibble(
    x = 1,
    y = 1
  )

  p2a <- ggplot2::ggplot(
    placeholder,
    ggplot2::aes(x, y)
  ) +
    ggplot2::annotate(
      "text",
      x = 1,
      y = 1,
      label =
        "No degree>=3 branch vertices.",
      size = 5
    ) +
    ggplot2::theme_void()

  p2b <- p2a
  p2c <- p2a
  p2d <- p2a
}

grDevices::pdf(
  OUT_FIG2,
  width = 14,
  height = 10,
  useDingbats = FALSE
)

print(
  (p2a | p2b) /
    (p2c | p2d) +
    patchwork::plot_annotation(
      title =
        "Step 5B-R2 | Corrected branch-topology audit"
    )
)

grDevices::dev.off()

# ==============================================================================
# 15. FIGURE 3 -- OBJECTIVE BRANCH ARMS
# ==============================================================================

msg(
  "Creating Fig5BR2_03_objective_branch_arms.pdf..."
)

if (!is.na(primary_branchpoint)) {

  SummarizedExperiment::colData(
    cds_single
  )$step5BR2_branch_arm_display <-
    ifelse(
      is.na(
        cell_df$branch_arm
      ),
      "Unassigned",
      cell_df$branch_arm
    )

  p3a <- monocle3::plot_cells(
    cds_single,
    reduction_method = "UMAP",
    color_cells_by =
      "step5BR2_branch_arm_display",
    show_trajectory_graph = TRUE,
    label_cell_groups = FALSE,
    label_leaves = TRUE,
    label_branch_points = TRUE,
    graph_label_size = 3,
    cell_size = 0.45
  ) +
    ggplot2::ggtitle(
      paste0(
        "Objective branch arms from ",
        primary_branchpoint
      )
    )

  arm_time_df <- cell_df |>
    dplyr::filter(
      !is.na(
        branch_arm
      ),
      grepl(
        "^Arm_",
        branch_arm
      )
    ) |>
    dplyr::count(
      branch_arm,
      timepoint,
      name = "n_cells"
    ) |>
    dplyr::group_by(
      branch_arm
    ) |>
    dplyr::mutate(
      fraction =
        n_cells /
        sum(n_cells)
    ) |>
    dplyr::ungroup()

  p3b <- ggplot2::ggplot(
    arm_time_df,
    ggplot2::aes(
      x = branch_arm,
      y = fraction,
      fill = timepoint
    )
  ) +
    ggplot2::geom_col(
      width = 0.70
    ) +
    ggplot2::scale_fill_manual(
      values = TIME_COLORS,
      name = "True time"
    ) +
    ggplot2::labs(
      title =
        "True-time composition of downstream arms",
      x = "Graph arm",
      y = "Fraction"
    ) +
    ggplot2::theme_classic(
      base_size = 11
    )

  near_bp_cells <- cells_near_graph_vertex(
    pg,
    closest_vertex,
    primary_branchpoint,
    radius =
      BRANCHPOINT_GRAPH_RADIUS
  )

  near_bp_df <- cell_df |>
    dplyr::filter(
      cell %in%
        near_bp_cells
    )

  p3c <- ggplot2::ggplot(
    near_bp_df,
    ggplot2::aes(
      x = factor(
        timepoint,
        levels =
          INJURY_TIME_LEVELS
      ),
      fill =
        shared_transition_3dpi_highconf
    )
  ) +
    ggplot2::geom_bar() +
    ggplot2::scale_fill_manual(
      values = c(
        "FALSE" = "#D9D9D9",
        "TRUE" = "#D73027"
      ),
      name =
        "High-conf\n3 dpi transition"
    ) +
    ggplot2::labs(
      title = paste0(
        "Cells near ",
        primary_branchpoint
      ),
      x = "True time",
      y = "Cells"
    ) +
    ggplot2::theme_classic(
      base_size = 11
    )

  p3d <- ggplot2::ggplot(
    arm_summary,
    ggplot2::aes(
      x = branch_arm,
      y = n_7dpi,
      fill = major_late_arm
    )
  ) +
    ggplot2::geom_col(
      width = 0.70
    ) +
    ggplot2::geom_hline(
      yintercept =
        MIN_7DPI_CELLS_PER_DOWNSTREAM_ARM,
      linetype = 2,
      linewidth = 0.45
    ) +
    ggplot2::scale_fill_manual(
      values = c(
        "FALSE" = "#BDBDBD",
        "TRUE" = "#2C7BB6"
      )
    ) +
    ggplot2::labs(
      title =
        "7 dpi support of objective downstream arms",
      x = "Graph arm",
      y = "7 dpi cells",
      fill = "Major late arm"
    ) +
    ggplot2::theme_classic(
      base_size = 11
    )

} else {

  p3a <- monocle3::plot_cells(
    cds_single,
    reduction_method = "UMAP",
    color_cells_by = "pseudotime",
    show_trajectory_graph = TRUE,
    label_cell_groups = FALSE,
    label_leaves = TRUE,
    label_branch_points = TRUE,
    graph_label_size = 3,
    cell_size = 0.45
  ) +
    ggplot2::ggtitle(
      "No objective two-major-7dpi-arm bifurcation"
    )

  placeholder <- tibble::tibble(
    x = 1,
    y = 1
  )

  p3b <- ggplot2::ggplot(
    placeholder,
    ggplot2::aes(x, y)
  ) +
    ggplot2::annotate(
      "text",
      x = 1,
      y = 1,
      label =
        "No branch-arm assignment was forced.",
      size = 5
    ) +
    ggplot2::theme_void()

  p3c <- p3b
  p3d <- p3b
}

grDevices::pdf(
  OUT_FIG3,
  width = 14,
  height = 10,
  useDingbats = FALSE
)

print(
  (p3a | p3b) /
    (p3c | p3d) +
    patchwork::plot_annotation(
      title =
        "Step 5B-R2 | Objective downstream branch arms"
    )
)

grDevices::dev.off()

# ==============================================================================
# 16. FIGURE 4 -- POST-HOC REPAIR / SCAR
# ==============================================================================

msg(
  "Creating Fig5BR2_04_posthoc_repair_scar_by_arm.pdf..."
)

SummarizedExperiment::colData(
  cds_single
)$repair_minus_scar <-
  cell_df$repair_minus_scar

p4a <- monocle3::plot_cells(
  cds_single,
  reduction_method = "UMAP",
  color_cells_by =
    REPAIR_SCORE_COL,
  show_trajectory_graph = TRUE,
  label_cell_groups = FALSE,
  label_leaves = FALSE,
  label_branch_points = TRUE,
  graph_label_size = 3,
  cell_size = 0.45
) +
  ggplot2::ggtitle(
    "Post-hoc repair-related program"
  )

p4b <- monocle3::plot_cells(
  cds_single,
  reduction_method = "UMAP",
  color_cells_by =
    SCAR_SCORE_COL,
  show_trajectory_graph = TRUE,
  label_cell_groups = FALSE,
  label_leaves = FALSE,
  label_branch_points = TRUE,
  graph_label_size = 3,
  cell_size = 0.45
) +
  ggplot2::ggtitle(
    "Post-hoc scar-related program"
  )

p4c <- monocle3::plot_cells(
  cds_single,
  reduction_method = "UMAP",
  color_cells_by =
    "repair_minus_scar",
  show_trajectory_graph = TRUE,
  label_cell_groups = FALSE,
  label_leaves = FALSE,
  label_branch_points = TRUE,
  graph_label_size = 3,
  cell_size = 0.45
) +
  ggplot2::ggtitle(
    "Post-hoc repair minus scar"
  )

if (
  !is.na(
    primary_branchpoint
  ) &&
  nrow(
    branch_program_sample_summary
  ) > 0
) {

  branch_long <-
    branch_program_sample_summary |>
    dplyr::select(
      sample,
      branch_arm,
      repair_median,
      scar_median
    ) |>
    tidyr::pivot_longer(
      cols = c(
        repair_median,
        scar_median
      ),
      names_to = "program",
      values_to =
        "sample_median_score"
    ) |>
    dplyr::mutate(
      program =
        dplyr::recode(
          program,
          repair_median =
            "Repair-related",
          scar_median =
            "Scar-related"
        )
    )

  p4d <- ggplot2::ggplot(
    branch_long,
    ggplot2::aes(
      x = branch_arm,
      y =
        sample_median_score,
      shape = sample
    )
  ) +
    ggplot2::geom_point(
      size = 3
    ) +
    ggplot2::facet_wrap(
      ~program,
      nrow = 1,
      scales = "free_y"
    ) +
    ggplot2::labs(
      title =
        "7 dpi arm-level programs by biological sample",
      subtitle =
        "Post-hoc only; graph topology fixed first",
      x = "Graph arm",
      y =
        "Sample median program score"
    ) +
    ggplot2::theme_classic(
      base_size = 10
    )

} else {

  p4d <- ggplot2::ggplot(
    tibble::tibble(
      x = 1,
      y = 1
    ),
    ggplot2::aes(x, y)
  ) +
    ggplot2::annotate(
      "text",
      x = 1,
      y = 1,
      label =
        "No arm-level repair/scar comparison:\nno objective bifurcation was forced.",
      size = 4.5
    ) +
    ggplot2::theme_void()
}

grDevices::pdf(
  OUT_FIG4,
  width = 14,
  height = 10,
  useDingbats = FALSE
)

print(
  (p4a | p4b) /
    (p4c | p4d) +
    patchwork::plot_annotation(
      title =
        "Step 5B-R2 | Step4 programs after corrected topology audit"
    )
)

grDevices::dev.off()

# ==============================================================================
# 17. FIGURE 5 -- OLD VS SINGLE ROOT ORIENTATION
# ==============================================================================

msg(
  "Creating Fig5BR2_05_old_vs_single_root_orientation.pdf..."
)

if (!all(is.na(cell_df$old_pseudotime))) {

  compare_df <- cell_df |>
    dplyr::filter(
      is.finite(
        old_pseudotime
      ),
      is.finite(
        pseudotime_single_root
      )
    )

  p5a <- ggplot2::ggplot(
    compare_df,
    ggplot2::aes(
      x = old_pseudotime,
      y =
        pseudotime_single_root
    )
  ) +
    ggplot2::geom_point(
      size = 0.45,
      alpha = 0.25
    ) +
    ggplot2::labs(
      title =
        "Old 225-cell-root versus single-principal-node pseudotime",
      subtitle = paste0(
        "Spearman = ",
        round(
          old_vs_single_spearman,
          4
        )
      ),
      x =
        "Previous Step5B-R pseudotime",
      y =
        "Single-root-node pseudotime"
    ) +
    ggplot2::theme_classic(
      base_size = 11
    )

  compare_time <- compare_df |>
    dplyr::group_by(
      timepoint
    ) |>
    dplyr::summarise(
      old_median =
        median(
          old_pseudotime
        ),
      single_root_median =
        median(
          pseudotime_single_root
        ),
      .groups = "drop"
    ) |>
    tidyr::pivot_longer(
      cols = c(
        old_median,
        single_root_median
      ),
      names_to =
        "orientation",
      values_to =
        "median_pseudotime"
    ) |>
    dplyr::mutate(
      timepoint =
        factor(
          timepoint,
          levels =
            INJURY_TIME_LEVELS,
          ordered = TRUE
        ),
      orientation =
        dplyr::recode(
          orientation,
          old_median =
            "Previous 225-cell root",
          single_root_median =
            "Single principal root node"
        )
    )

  p5b <- ggplot2::ggplot(
    compare_time,
    ggplot2::aes(
      x = timepoint,
      y = median_pseudotime,
      group = orientation,
      linetype = orientation
    )
  ) +
    ggplot2::geom_line(
      linewidth = 0.8
    ) +
    ggplot2::geom_point(
      size = 2.5
    ) +
    ggplot2::labs(
      title =
        "True-time ordering under two orientations",
      x = "True time",
      y = "Median pseudotime",
      linetype = "Orientation"
    ) +
    ggplot2::theme_classic(
      base_size = 11
    )

} else {

  placeholder <- tibble::tibble(
    x = 1,
    y = 1
  )

  p5a <- ggplot2::ggplot(
    placeholder,
    ggplot2::aes(x, y)
  ) +
    ggplot2::annotate(
      "text",
      x = 1,
      y = 1,
      label =
        "Previous Step5B-R pseudotime was not stored.",
      size = 5
    ) +
    ggplot2::theme_void()

  p5b <- p5a
}

root_vertex_plot <- root_vertex_counts |>
  dplyr::mutate(
    is_anchor =
      principal_vertex ==
        ROOT_ANCHOR_VERTEX
  )

p5c <- ggplot2::ggplot(
  root_vertex_plot,
  ggplot2::aes(
    x = reorder(
      principal_vertex,
      n_primary_root_cells
    ),
    y =
      n_primary_root_cells,
    fill =
      is_anchor
  )
) +
  ggplot2::geom_col() +
  ggplot2::coord_flip() +
  ggplot2::scale_fill_manual(
    values = c(
      "FALSE" = "#BDBDBD",
      "TRUE" = ROOT_COL
    ),
    guide = "none"
  ) +
  ggplot2::labs(
    title =
      "Distribution of the original 225 root cells across graph vertices",
    subtitle = paste0(
      "Single anchor = ",
      ROOT_ANCHOR_VERTEX
    ),
    x =
      "Principal graph vertex",
    y =
      "Primary root cells"
  ) +
  ggplot2::theme_classic(
    base_size = 10
  )

p5d <- ggplot2::ggplot(
  sample_summary,
  ggplot2::aes(
    x = true_day,
    y = median_pseudotime,
    label = sample
  )
) +
  ggplot2::geom_point(
    size = 2.7
  ) +
  ggplot2::geom_smooth(
    method = "lm",
    se = FALSE,
    linewidth = 0.7
  ) +
  ggplot2::geom_text(
    size = 3,
    vjust = -0.6,
    check_overlap = TRUE
  ) +
  ggplot2::labs(
    title =
      "Biological-sample median pseudotime",
    subtitle = paste0(
      "Spearman = ",
      round(
        sample_level_time_spearman,
        3
      )
    ),
    x =
      "True day after injury",
    y =
      "Median single-root pseudotime"
  ) +
  ggplot2::theme_classic(
    base_size = 11
  )

grDevices::pdf(
  OUT_FIG5,
  width = 14,
  height = 10,
  useDingbats = FALSE
)

print(
  (p5a | p5b) /
    (p5c | p5d) +
    patchwork::plot_annotation(
      title =
        "Step 5B-R2 | Orientation and root-topology audit"
    )
)

grDevices::dev.off()

# ==============================================================================
# 18. EXCEL
# ==============================================================================

msg(
  "Writing Table5BR2_single_root_branch_audit.xlsx..."
)

parameter_tbl <- tibble::tibble(
  parameter = c(
    "Existing CDS cells",
    "Existing graph vertices",
    "Existing graph edges",
    "Existing graph connected components",
    "Original primary root cells",
    "Original root principal vertices",
    "Single root principal node",
    "Single-root finite pseudotime fraction",
    "Single-root cell-level true-time Spearman",
    "Single-root sample-level true-time Spearman",
    "Single-root 1->3->7 medians monotonic",
    "Old-vs-single-root pseudotime Spearman",
    "Degree>=3 branch vertices",
    "Objective primary bifurcation"
  ),
  value = c(
    ncol(cds_single),
    igraph::vcount(pg),
    igraph::ecount(pg),
    igraph::components(pg)$no,
    length(primary_root_cells),
    nrow(root_vertex_counts),
    ROOT_ANCHOR_VERTEX,
    finite_fraction,
    cell_level_time_spearman,
    sample_level_time_spearman,
    time_monotonic,
    old_vs_single_spearman,
    length(branch_vertices),
    ifelse(
      is.na(
        primary_branchpoint
      ),
      "None",
      primary_branchpoint
    )
  )
)

wb <- openxlsx::createWorkbook()

safe_add_sheet(
  wb,
  "Parameters",
  parameter_tbl
)

safe_add_sheet(
  wb,
  "Root_vertex_counts",
  root_vertex_counts
)

safe_add_sheet(
  wb,
  "Time_pseudotime",
  time_summary
)

safe_add_sheet(
  wb,
  "Sample_pseudotime",
  sample_summary
)

safe_add_sheet(
  wb,
  "Reference_states",
  reference_summary
)

safe_add_sheet(
  wb,
  "Graph_vertices",
  graph_vertex_table
)

safe_add_sheet(
  wb,
  "Branchpoint_summary",
  branchpoint_summary
)

safe_add_sheet(
  wb,
  "Branch_component_summary",
  branch_component_summary
)

safe_add_sheet(
  wb,
  "Arm_summary",
  arm_summary
)

safe_add_sheet(
  wb,
  "Branch_program_7dpi",
  branch_program_summary
)

safe_add_sheet(
  wb,
  "Branch_program_sample",
  branch_program_sample_summary
)

safe_add_sheet(
  wb,
  "Cell_results",
  cell_df
)

openxlsx::saveWorkbook(
  wb,
  OUT_XLSX,
  overwrite = TRUE
)

# ==============================================================================
# 19. SAVE OBJECTS
# ==============================================================================

saveRDS(
  cds_single,
  OUT_CDS,
  compress = TRUE
)

saveRDS(
  obj,
  OUT_RDS,
  compress = TRUE
)

# ==============================================================================
# 20. REPORT
# ==============================================================================

msg(
  "Writing Step5BR2_report.txt..."
)

report_lines <- c(
  "======================================================================",
  "GSE162610 | STEP 5B-R2 REPORT",
  "Single-principal-root orientation + corrected branch-topology audit",
  "======================================================================",
  "",
  "SCOPE:",
  "The Step5B-R principal graph was reused WITHOUT relearning.",
  "Only pseudotime orientation and branch-component mapping were re-evaluated.",
  "",
  "-------------------------- EXISTING GRAPH ------------------------------",
  paste0(
    "Cells: ",
    ncol(cds_single)
  ),
  paste0(
    "Graph vertices: ",
    igraph::vcount(pg)
  ),
  paste0(
    "Graph edges: ",
    igraph::ecount(pg)
  ),
  paste0(
    "Connected components: ",
    igraph::components(pg)$no
  ),
  "",
  "---------------------- SINGLE PRINCIPAL ROOT ---------------------------",
  paste0(
    "Original Step5B-R primary root cells: ",
    length(primary_root_cells)
  ),
  paste0(
    "Principal vertices occupied by those root cells: ",
    nrow(root_vertex_counts)
  ),
  paste0(
    "Single objective root principal node: ",
    ROOT_ANCHOR_VERTEX
  ),
  "",
  "Original root-cell distribution across principal vertices:",
  capture_table(
    root_vertex_counts
  ),
  "",
  "---------------------- SINGLE-ROOT TIME QC -----------------------------",
  paste0(
    "Finite pseudotime fraction: ",
    round(
      100 *
        finite_fraction,
      2
    ),
    "%"
  ),
  paste0(
    "Cell-level true-time/pseudotime Spearman: ",
    round(
      cell_level_time_spearman,
      4
    )
  ),
  paste0(
    "Sample-level true-time/median-pseudotime Spearman: ",
    round(
      sample_level_time_spearman,
      4
    )
  ),
  paste0(
    "1 -> 3 -> 7 dpi medians strictly monotonic: ",
    time_monotonic
  ),
  paste0(
    "Previous cell-root vs single-node pseudotime Spearman: ",
    round(
      old_vs_single_spearman,
      4
    )
  ),
  "",
  "Pseudotime by true time:",
  capture_table(
    time_summary
  ),
  "",
  "Formal biological references:",
  capture_table(
    reference_summary
  ),
  "",
  "-------------------- CORRECTED BRANCH AUDIT ---------------------------",
  "Critical implementation rule:",
  "  cell barcode -> closest principal vertex -> component membership",
  "is now mapped explicitly with the CELL barcode retained as the vector name.",
  "",
  "This distinguishes:",
  "  1. graph-topology components after deleting a degree>=3 vertex;",
  "  2. components that actually contain mapped biological cells;",
  "  3. major late arms containing >=100 total and >=100 7 dpi cells.",
  "",
  if (
    nrow(
      branchpoint_summary
    ) > 0
  ) {
    paste0(
      "Branch-point summary:\n",
      capture_table(
        branchpoint_summary
      )
    )
  } else {
    "No degree>=3 branch vertices detected."
  },
  "",
  if (
    nrow(
      branch_component_summary
    ) > 0
  ) {
    paste0(
      "Per-component cell support:\n",
      capture_table(
        branch_component_summary
      )
    )
  } else {
    "No downstream component table generated."
  },
  "",
  paste0(
    "Objective primary bifurcation: ",
    ifelse(
      is.na(
        primary_branchpoint
      ),
      "None satisfying predefined two-major-7dpi-arm criterion",
      primary_branchpoint
    )
  ),
  "",
  "--------------------- POST-HOC REPAIR / SCAR --------------------------",
  "Repair/scar scores were not used to orient the graph or select branchpoints.",
  if (
    nrow(
      branch_program_summary
    ) > 0
  ) {
    capture_table(
      branch_program_summary
    )
  } else {
    "No arm-level program comparison because no objective bifurcation was assigned."
  },
  "",
  "--------------------------- DECISION RULE ------------------------------",
  "Interpretation after Step5B-R2:",
  "  A. If the corrected audit now reveals >=2 large cell-supported late arms,",
  "     inspect whether the branch point is biologically positioned near the",
  "     3 dpi shared-transition state before any repair/scar naming.",
  "  B. If graph topology has side components but they contain few/no cells,",
  "     they are graph stubs and should not be called biological fate branches.",
  "  C. If no branch point has two substantial 7 dpi arms even after corrected",
  "     mapping, the current dataset does not support a discrete two-arm",
  "     repair/scar bifurcation under this Monocle topology.",
  "",
  "Formal outputs:",
  paste0("  ", OUT_FIG1),
  paste0("  ", OUT_FIG2),
  paste0("  ", OUT_FIG3),
  paste0("  ", OUT_FIG4),
  paste0("  ", OUT_FIG5),
  paste0("  ", OUT_XLSX),
  paste0("  ", OUT_CDS),
  paste0("  ", OUT_RDS),
  paste0("  ", OUT_REP),
  paste0("  ", OUT_SESS)
)

writeLines(
  report_lines,
  OUT_REP,
  useBytes = TRUE
)

writeLines(
  capture.output(
    sessionInfo()
  ),
  OUT_SESS,
  useBytes = TRUE
)

# ==============================================================================
# 21. FINAL CONSOLE SUMMARY
# ==============================================================================

msg("====================================================")
msg("STEP 5B-R2 COMPLETE")
msg("====================================================")

msg(
  "Single principal root node: ",
  ROOT_ANCHOR_VERTEX
)

msg(
  "Cell-level true-time Spearman: ",
  round(
    cell_level_time_spearman,
    4
  )
)

msg(
  "Sample-level true-time Spearman: ",
  round(
    sample_level_time_spearman,
    4
  )
)

msg(
  "1->3->7 dpi medians monotonic: ",
  time_monotonic
)

msg(
  "Old vs single-root pseudotime Spearman: ",
  round(
    old_vs_single_spearman,
    4
  )
)

msg(
  "Degree>=3 branch vertices: ",
  length(
    branch_vertices
  )
)

msg(
  "Objective primary bifurcation: ",
  ifelse(
    is.na(
      primary_branchpoint
    ),
    "NONE",
    primary_branchpoint
  )
)

msg("")
msg("请上传以下6个结果用于最终判断：")
msg("  1) Fig5BR2_01_single_principal_root_pseudotime.pdf")
msg("  2) Fig5BR2_02_corrected_branch_topology_audit.pdf")
msg("  3) Fig5BR2_03_objective_branch_arms.pdf")
msg("  4) Fig5BR2_04_posthoc_repair_scar_by_arm.pdf")
msg("  5) Fig5BR2_05_old_vs_single_root_orientation.pdf")
msg("  6) Step5BR2_report.txt")
msg("")
msg("这一轮先判定真实拓扑，不进入branch-specific DEG。")
