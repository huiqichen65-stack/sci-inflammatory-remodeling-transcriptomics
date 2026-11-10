# =============================================================================
# 57_step10D_GSE234774_marker_first_major_lineage_annotation_audit.R
#
# Adult SCI activated/DAM-like microglia project
# GSE234774 Step 10D:
#   Marker-first whole-atlas major-lineage annotation audit
#
# FROZEN INPUT CHOICE FROM STEP 10C
#   Reduction: harmony_conservative
#   UMAP:      umap_harmony_conservative
#   Clusters:  hcon_res_0_2
#   Resolution: 0.2
#
# INPUT
#   <PROJECT_ROOT>/04_results/
#     step10C_GSE234774_Harmony_resolution_sensitivity/RDS/
#       GSE234774_step10C_PCA_Harmony_resolution_comparison_BPCells.rds
#
# PURPOSE
#   1. Discover positive RNA markers for all 25 frozen whole-atlas clusters.
#   2. Quantify curated mouse spinal-cord lineage-marker panels.
#   3. Produce review-only candidate major-lineage labels.
#   4. Focus on resident microglia versus infiltrating monocyte/macrophage,
#      border-associated macrophage, neutrophil, and dendritic evidence.
#   5. Cross-check author-provided hierarchical labels without using them to
#      assign the candidate labels.
#
# IMPORTANT BOUNDARIES
#   - Harmony is used only for cluster membership and visualization.
#   - Marker testing uses the unmodified RNA LogNormalize data layer.
#   - Author labels are supporting evidence only.
#   - Candidate labels produced here are NOT final annotations.
#   - No cluster is removed in this script.
#   - No resident-microglia object is frozen in this script.
#   - No GSE162610 Arm_1/Arm_2 score is calculated or redefined.
#   - No differential expression across time, pseudobulk, trajectory, regulon,
#     CellChat, or spatial analysis is run.
#
# OUTPUT REVIEW RULE
#   Upload exactly seven files:
#     3 figures + 2 evidence tables + 1 audit table + 1 final report.
# =============================================================================

options(
  stringsAsFactors = FALSE,
  warn = 1,
  Seurat.object.assay.version = "v5",
  future.globals.maxSize = 16 * 1024^3
)

# -----------------------------------------------------------------------------
# 0. User controls
# -----------------------------------------------------------------------------
SEED <- 20260803L

FROZEN_REDUCTION <- "harmony_conservative"
FROZEN_UMAP <- "umap_harmony_conservative"
FROZEN_CLUSTER_COLUMN <- "hcon_res_0_2"
EXPECTED_N_CLUSTERS <- 25L

# Marker discovery is downsampled only for speed.
# All cluster-level curated-panel summaries still use all 144,359 singlets.
MAX_CELLS_PER_CLUSTER_FOR_MARKERS <- 2500L
MARKER_LOGFC_THRESHOLD <- 0.25
MARKER_MIN_PCT <- 0.10
MARKER_MIN_DIFF_PCT <- 0.05
MARKER_RETURN_P_ADJ <- 0.05
TOP_MARKERS_PER_CLUSTER <- 30L

# Review-only confidence rules for the automatic candidate labels.
HIGH_CONFIDENCE_MARGIN <- 0.50
MODERATE_CONFIDENCE_MARGIN <- 0.25
HIGH_CONFIDENCE_SUPPORT <- 0.40
MODERATE_CONFIDENCE_SUPPORT <- 0.25
HIGH_CONFIDENCE_MEAN_DETECTION <- 0.10
MODERATE_CONFIDENCE_MEAN_DETECTION <- 0.05

# Resident versus infiltrating-myeloid review margin.
MYELOID_CLASS_MARGIN <- 0.40
POSSIBLE_MYELOID_SCORE <- 1.00

expected_timepoints <- c(
  "Uninjured",
  "1 dpi",
  "4 dpi",
  "7 dpi",
  "14 dpi",
  "1 month",
  "2 months"
)

# -----------------------------------------------------------------------------
# 1. Curated review panels
# -----------------------------------------------------------------------------
# These panels are intentionally broad. They are used to organize marker review,
# not to replace manual examination of cluster markers.
marker_panels_requested <- list(
  Resident_microglia = c(
    "P2ry12",
    "Tmem119",
    "Sall1",
    "Fcrls",
    "Hexb",
    "Cx3cr1",
    "Gpr34",
    "Olfml3"
  ),
  Monocyte_macrophage = c(
    "Lyz2",
    "Ccr2",
    "Ly6c2",
    "Plac8",
    "Ms4a7",
    "Lgals3",
    "Fcgr3",
    "Cd68"
  ),
  Border_associated_macrophage = c(
    "Mrc1",
    "Lyve1",
    "Pf4",
    "Folr2",
    "Cd163",
    "Cbr2"
  ),
  Neutrophil = c(
    "S100a8",
    "S100a9",
    "Ly6g",
    "Csf3r",
    "Retnlg"
  ),
  Dendritic = c(
    "Flt3",
    "Zbtb46",
    "Clec10a",
    "Cd209a",
    "Cd74",
    "H2-Ab1"
  ),
  T_cell = c(
    "Cd3d",
    "Cd3e",
    "Cd3g",
    "Trbc1",
    "Lck",
    "Il7r"
  ),
  NK_cell = c(
    "Nkg7",
    "Ncr1",
    "Klrd1",
    "Prf1",
    "Gzmb"
  ),
  B_cell = c(
    "Cd79a",
    "Ms4a1",
    "Cd37",
    "Cd74",
    "H2-Aa"
  ),
  Astrocyte = c(
    "Aldh1l1",
    "Aqp4",
    "Slc1a3",
    "Gfap",
    "S100b",
    "Sox9"
  ),
  Oligodendrocyte = c(
    "Plp1",
    "Mbp",
    "Mog",
    "Mag",
    "Cnp",
    "Sox10"
  ),
  OPC = c(
    "Pdgfra",
    "Cspg4",
    "Olig1",
    "Olig2",
    "Sox10"
  ),
  Ependymal = c(
    "Foxj1",
    "Pifo",
    "Dynlrb2",
    "Tppp3",
    "Ccdc153",
    "Cfap44"
  ),
  Neuron = c(
    "Snap25",
    "Syt1",
    "Rbfox3",
    "Tubb3",
    "Elavl3",
    "Stmn2"
  ),
  Endothelial = c(
    "Pecam1",
    "Cldn5",
    "Kdr",
    "Esam",
    "Vwf",
    "Tek"
  ),
  Pericyte = c(
    "Pdgfrb",
    "Rgs5",
    "Cspg4",
    "Kcnj8",
    "Abcc9",
    "Des"
  ),
  Vascular_smooth_muscle = c(
    "Acta2",
    "Tagln",
    "Myh11",
    "Cnn1"
  ),
  Fibroblast = c(
    "Col1a1",
    "Col1a2",
    "Col3a1",
    "Dcn",
    "Lum",
    "Pdgfra"
  ),
  Schwann_cell = c(
    "Mpz",
    "Pmp22",
    "Prx",
    "Ngfr",
    "Sox10"
  )
)

myeloid_panel_names <- c(
  "Resident_microglia",
  "Monocyte_macrophage",
  "Border_associated_macrophage",
  "Neutrophil",
  "Dendritic"
)

# Additional markers shown in figures but not used as independent candidate
# panels.
supplementary_plot_markers_requested <- list(
  Myeloid_general = c(
    "Ptprc",
    "Aif1",
    "Tyrobp",
    "Lyz2",
    "Csf1r"
  ),
  Activated_phagolysosomal = c(
    "Apoe",
    "Cst7",
    "Ctsb",
    "Ctsd",
    "Trem2",
    "Itgax"
  ),
  Excitatory_neuron = c(
    "Slc17a6",
    "Slc17a7",
    "Slc17a8"
  ),
  Inhibitory_neuron = c(
    "Gad1",
    "Gad2",
    "Slc32a1"
  )
)

# -----------------------------------------------------------------------------
# 2. Paths
# -----------------------------------------------------------------------------
project_root <- Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
input_object_path <- file.path(
  project_root,
  "04_results",
  "step10C_GSE234774_Harmony_resolution_sensitivity",
  "RDS",
  "GSE234774_step10C_PCA_Harmony_resolution_comparison_BPCells.rds"
)

output_dir <- file.path(
  project_root,
  "04_results",
  "step10D_GSE234774_marker_first_major_lineage_annotation"
)

tables_dir <- file.path(output_dir, "tables")
figures_dir <- file.path(output_dir, "figures")
rds_dir <- file.path(output_dir, "RDS")

for (directory in c(
  output_dir,
  tables_dir,
  figures_dir,
  rds_dir
)) {
  dir.create(
    directory,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

run_log <- file.path(
  output_dir,
  "00_run_log.txt"
)

if (file.exists(run_log)) {
  file.remove(run_log)
}

output_object_path <- file.path(
  rds_dir,
  "GSE234774_step10D_major_lineage_annotation_evidence_BPCells.rds"
)

# -----------------------------------------------------------------------------
# 3. Helpers
# -----------------------------------------------------------------------------
log_message <- function(...) {
  message_text <- paste0(...)
  timestamp <- format(
    Sys.time(),
    "%Y-%m-%d %H:%M:%S"
  )
  line <- paste0(
    "[",
    timestamp,
    "] ",
    message_text
  )

  cat(line, "\n")
  cat(
    line,
    "\n",
    file = run_log,
    append = TRUE
  )
}

write_csv_safe <- function(
    x,
    filename
) {
  path <- file.path(
    tables_dir,
    filename
  )

  utils::write.csv(
    x,
    path,
    row.names = FALSE,
    na = ""
  )

  log_message(
    "Wrote: ",
    path
  )

  invisible(path)
}

stop_with_report <- function(reason) {
  report_path <- file.path(
    output_dir,
    "10_STEP10D_FINAL_REPORT.txt"
  )

  report <- c(
    "GSE234774 STEP 10D MARKER-FIRST MAJOR-LINEAGE ANNOTATION",
    "============================================================",
    paste0(
      "Run time: ",
      format(
        Sys.time(),
        "%Y-%m-%d %H:%M:%S"
      )
    ),
    paste0(
      "Output directory: ",
      output_dir
    ),
    "",
    "FINAL DECISION: STOP",
    "",
    paste0(
      "Reason: ",
      reason
    ),
    "",
    "Do not freeze major-lineage labels or extract resident microglia."
  )

  writeLines(
    report,
    report_path,
    useBytes = TRUE
  )

  stop(
    reason,
    call. = FALSE
  )
}

safe_package_version <- function(package_name) {
  if (!requireNamespace(
    package_name,
    quietly = TRUE
  )) {
    return(NA_character_)
  }

  as.character(
    utils::packageVersion(
      package_name
    )
  )
}

numeric_cluster_order <- function(cluster_values) {
  cluster_values <- unique(
    as.character(
      cluster_values
    )
  )

  numeric_values <- suppressWarnings(
    as.numeric(
      cluster_values
    )
  )

  if (all(is.finite(
    numeric_values
  ))) {
    return(
      cluster_values[
        order(
          numeric_values
        )
      ]
    )
  }

  sort(
    cluster_values
  )
}

resolve_features <- function(
    requested_features,
    available_features
) {
  available_upper <- toupper(
    available_features
  )

  first_index <- !duplicated(
    available_upper
  )

  lookup <- setNames(
    available_features[
      first_index
    ],
    available_upper[
      first_index
    ]
  )

  resolved <- unname(
    lookup[
      toupper(
        requested_features
      )
    ]
  )

  names(resolved) <- requested_features
  resolved
}

deduplicate_plot_panels <- function(panel_list) {
  used <- character(0)
  output <- list()

  for (panel_name in names(
    panel_list
  )) {
    genes <- panel_list[[panel_name]]

    genes <- genes[
      !genes %in% used
    ]

    if (length(genes) > 0L) {
      output[[panel_name]] <- genes

      used <- c(
        used,
        genes
      )
    }
  }

  output
}

dominant_value_by_cluster <- function(
    metadata,
    cluster_column,
    value_column
) {
  values <- as.character(
    metadata[
      ,
      value_column
    ]
  )

  clusters <- as.character(
    metadata[
      ,
      cluster_column
    ]
  )

  valid <- !is.na(values) &
    nzchar(values)

  if (!any(valid)) {
    return(
      data.frame(
        cluster = character(0),
        dominant_value = character(0),
        dominant_fraction = numeric(0),
        stringsAsFactors = FALSE
      )
    )
  }

  counts <- data.table::data.table(
    cluster = clusters[
      valid
    ],
    value = values[
      valid
    ]
  )[
    ,
    .(
      n_cells = .N
    ),
    by = .(
      cluster,
      value
    )
  ]

  totals <- counts[
    ,
    .(
      total_cells = sum(
        n_cells
      )
    ),
    by = cluster
  ]

  counts <- merge(
    counts,
    totals,
    by = "cluster",
    all.x = TRUE,
    sort = FALSE
  )

  counts$fraction <- counts$n_cells /
    counts$total_cells

  top <- counts[
    order(
      cluster,
      -fraction,
      value
    ),
    .SD[1L],
    by = cluster
  ]

  data.frame(
    cluster = top$cluster,
    dominant_value = top$value,
    dominant_fraction = top$fraction,
    stringsAsFactors = FALSE
  )
}

safe_matrix_zscore_by_gene <- function(expression_matrix) {
  gene_means <- rowMeans(
    expression_matrix,
    na.rm = TRUE
  )

  gene_sds <- apply(
    expression_matrix,
    1L,
    stats::sd,
    na.rm = TRUE
  )

  gene_sds[
    !is.finite(
      gene_sds
    ) |
      gene_sds == 0
  ] <- 1

  centered <- sweep(
    expression_matrix,
    MARGIN = 1L,
    STATS = gene_means,
    FUN = "-"
  )

  zscores <- sweep(
    centered,
    MARGIN = 1L,
    STATS = gene_sds,
    FUN = "/"
  )

  zscores[
    !is.finite(
      zscores
    )
  ] <- 0

  zscores
}

# -----------------------------------------------------------------------------
# 4. Package preflight
# -----------------------------------------------------------------------------
log_message(
  "Starting package preflight."
)

required_packages <- c(
  "Matrix",
  "data.table",
  "ggplot2",
  "patchwork",
  "scales",
  "Seurat",
  "SeuratObject"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_packages) > 0L) {
  stop_with_report(
    paste0(
      "Missing required package(s): ",
      paste(
        missing_packages,
        collapse = ", "
      ),
      "."
    )
  )
}

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(Seurat)
  library(SeuratObject)
})

package_versions <- data.frame(
  package = c(
    "R",
    "Matrix",
    "data.table",
    "ggplot2",
    "patchwork",
    "scales",
    "Seurat",
    "SeuratObject"
  ),
  version = c(
    paste0(
      R.version$major,
      ".",
      R.version$minor
    ),
    safe_package_version("Matrix"),
    safe_package_version("data.table"),
    safe_package_version("ggplot2"),
    safe_package_version("patchwork"),
    safe_package_version("scales"),
    safe_package_version("Seurat"),
    safe_package_version("SeuratObject")
  ),
  stringsAsFactors = FALSE
)

write_csv_safe(
  package_versions,
  "00_package_versions.csv"
)

# -----------------------------------------------------------------------------
# 5. Load and audit Step 10C object
# -----------------------------------------------------------------------------
if (!file.exists(
  input_object_path
)) {
  stop_with_report(
    paste0(
      "Input Step 10C object not found: ",
      input_object_path
    )
  )
}

log_message(
  "Loading Step 10C comparison object."
)

obj <- readRDS(
  input_object_path
)

if (!inherits(
  obj,
  "Seurat"
)) {
  stop_with_report(
    "Input RDS is not a Seurat object."
  )
}

if (
  nrow(obj) != 26918L ||
    ncol(obj) != 144359L
) {
  stop_with_report(
    paste0(
      "Expected 26,918 features x 144,359 singlets; detected ",
      nrow(obj),
      " x ",
      ncol(obj),
      "."
    )
  )
}

required_metadata <- c(
  "biological_library",
  "timepoint",
  "replicate_id",
  "nCount_RNA",
  "nFeature_RNA",
  "percent.mt",
  "scDblFinder.class",
  FROZEN_CLUSTER_COLUMN
)

missing_metadata <- setdiff(
  required_metadata,
  names(
    obj@meta.data
  )
)

if (length(missing_metadata) > 0L) {
  stop_with_report(
    paste0(
      "Missing required metadata column(s): ",
      paste(
        missing_metadata,
        collapse = ", "
      ),
      "."
    )
  )
}

required_reductions <- c(
  FROZEN_REDUCTION,
  FROZEN_UMAP
)

missing_reductions <- setdiff(
  required_reductions,
  names(
    obj@reductions
  )
)

if (length(missing_reductions) > 0L) {
  stop_with_report(
    paste0(
      "Missing frozen reduction(s): ",
      paste(
        missing_reductions,
        collapse = ", "
      ),
      "."
    )
  )
}

if (!all(
  as.character(
    obj$scDblFinder.class
  ) == "singlet"
)) {
  stop_with_report(
    "The Step 10C object contains non-singlet nuclei."
  )
}

if (
  length(
    unique(
      as.character(
        obj$biological_library
      )
    )
  ) != 20L
) {
  stop_with_report(
    "Expected 20 biological libraries."
  )
}

if (!all(
  expected_timepoints %in%
    unique(
      as.character(
        obj$timepoint
      )
    )
)) {
  stop_with_report(
    "Expected timepoints are incomplete."
  )
}

cluster_values <- as.character(
  obj@meta.data[
    ,
    FROZEN_CLUSTER_COLUMN
  ]
)

cluster_levels <- numeric_cluster_order(
  cluster_values
)

if (
  length(
    cluster_levels
  ) != EXPECTED_N_CLUSTERS
) {
  stop_with_report(
    paste0(
      "Expected ",
      EXPECTED_N_CLUSTERS,
      " frozen clusters; detected ",
      length(
        cluster_levels
      ),
      "."
    )
  )
}

obj@meta.data[
  ,
  FROZEN_CLUSTER_COLUMN
] <- factor(
  cluster_values,
  levels = cluster_levels
)

Seurat::DefaultAssay(
  obj
) <- "RNA"

Seurat::Idents(
  obj
) <- FROZEN_CLUSTER_COLUMN

rna_layers <- SeuratObject::Layers(
  object = obj,
  assay = "RNA",
  search = NA
)

counts_layers <- rna_layers[
  grepl(
    "^counts($|\\.)",
    rna_layers
  )
]

data_layers <- rna_layers[
  grepl(
    "^data($|\\.)",
    rna_layers
  )
]

if (
  length(
    counts_layers
  ) != 1L ||
    length(
      data_layers
    ) != 1L
) {
  stop_with_report(
    paste0(
      "Expected one RNA counts layer and one RNA data layer; detected ",
      length(
        counts_layers
      ),
      " and ",
      length(
        data_layers
      ),
      "."
    )
  )
}

counts_layer_name <- counts_layers[[1L]]

data_layer_name <- data_layers[[1L]]

input_audit <- data.frame(
  item = c(
    "Input object",
    "Features",
    "Singlet nuclei",
    "Biological libraries",
    "Timepoints",
    "Frozen reduction",
    "Frozen UMAP",
    "Frozen cluster column",
    "Frozen cluster count",
    "RNA counts layer",
    "RNA LogNormalize data layer",
    "Harmony used for marker expression"
  ),
  value = c(
    normalizePath(
      input_object_path,
      winslash = "/",
      mustWork = TRUE
    ),
    nrow(obj),
    ncol(obj),
    length(
      unique(
        as.character(
          obj$biological_library
        )
      )
    ),
    length(
      unique(
        as.character(
          obj$timepoint
        )
      )
    ),
    FROZEN_REDUCTION,
    FROZEN_UMAP,
    FROZEN_CLUSTER_COLUMN,
    length(
      cluster_levels
    ),
    counts_layer_name,
    data_layer_name,
    FALSE
  ),
  stringsAsFactors = FALSE
)

write_csv_safe(
  input_audit,
  "01_input_object_audit.csv"
)

# -----------------------------------------------------------------------------
# 6. Resolve curated marker panels against the actual feature names
# -----------------------------------------------------------------------------
all_requested_features <- unique(
  c(
    unlist(
      marker_panels_requested,
      use.names = FALSE
    ),
    unlist(
      supplementary_plot_markers_requested,
      use.names = FALSE
    )
  )
)

feature_resolution <- resolve_features(
  requested_features =
    all_requested_features,
  available_features =
    rownames(obj)
)

feature_resolution_table <- data.frame(
  requested_gene = names(
    feature_resolution
  ),
  resolved_feature = unname(
    feature_resolution
  ),
  present = !is.na(
    feature_resolution
  ),
  stringsAsFactors = FALSE
)

write_csv_safe(
  feature_resolution_table,
  "02_curated_marker_feature_resolution.csv"
)

marker_panels_present <- lapply(
  marker_panels_requested,
  function(panel_genes) {
    resolved <- feature_resolution[
      panel_genes
    ]

    unique(
      unname(
        resolved[
          !is.na(
            resolved
          )
        ]
      )
    )
  }
)

supplementary_panels_present <- lapply(
  supplementary_plot_markers_requested,
  function(panel_genes) {
    resolved <- feature_resolution[
      panel_genes
    ]

    unique(
      unname(
        resolved[
          !is.na(
            resolved
          )
        ]
      )
    )
  }
)

panel_availability <- data.frame(
  panel = names(
    marker_panels_requested
  ),
  requested_n = vapply(
    marker_panels_requested,
    length,
    integer(1)
  ),
  present_n = vapply(
    marker_panels_present,
    length,
    integer(1)
  ),
  present_features = vapply(
    marker_panels_present,
    function(features) {
      paste(
        features,
        collapse = ";"
      )
    },
    character(1)
  ),
  stringsAsFactors = FALSE
)

panel_availability$present_fraction <-
  panel_availability$present_n /
  panel_availability$requested_n

write_csv_safe(
  panel_availability,
  "02B_curated_panel_availability.csv"
)

if (
  any(
    panel_availability$present_n < 2L
  )
) {
  failed_panels <- panel_availability$panel[
    panel_availability$present_n < 2L
  ]

  stop_with_report(
    paste0(
      "Fewer than two markers were resolved for panel(s): ",
      paste(
        failed_panels,
        collapse = ", "
      ),
      "."
    )
  )
}

# -----------------------------------------------------------------------------
# 7. Cluster-level curated marker expression and detection
# -----------------------------------------------------------------------------
log_message(
  "Calculating curated marker expression and detection using all singlets."
)

all_present_curated_features <- unique(
  unlist(
    marker_panels_present,
    use.names = FALSE
  )
)

cluster_factor <- factor(
  as.character(
    obj@meta.data[
      ,
      FROZEN_CLUSTER_COLUMN
    ]
  ),
  levels = cluster_levels
)

cluster_index <- as.integer(
  cluster_factor
)

cluster_sizes <- tabulate(
  cluster_index,
  nbins = length(
    cluster_levels
  )
)

names(
  cluster_sizes
) <- cluster_levels

cluster_indicator <- Matrix::sparseMatrix(
  i = seq_len(
    ncol(obj)
  ),
  j = cluster_index,
  x = 1,
  dims = c(
    ncol(obj),
    length(
      cluster_levels
    )
  ),
  dimnames = list(
    colnames(obj),
    cluster_levels
  )
)

curated_data <- SeuratObject::LayerData(
  object = obj,
  assay = "RNA",
  layer = data_layer_name,
  features = all_present_curated_features,
  fast = FALSE
)

curated_counts <- SeuratObject::LayerData(
  object = obj,
  assay = "RNA",
  layer = counts_layer_name,
  features = all_present_curated_features,
  fast = FALSE
)

curated_data <- methods::as(
  curated_data,
  "dgCMatrix"
)

curated_counts <- methods::as(
  curated_counts,
  "dgCMatrix"
)

if (
  !identical(
    colnames(
      curated_data
    ),
    colnames(obj)
  ) ||
    !identical(
      colnames(
        curated_counts
      ),
      colnames(obj)
    )
) {
  stop_with_report(
    "Curated expression matrices do not align with the Seurat cell order."
  )
}

sum_expression_by_cluster <- as.matrix(
  curated_data %*%
    cluster_indicator
)

mean_expression_by_cluster <- sweep(
  sum_expression_by_cluster,
  MARGIN = 2L,
  STATS = cluster_sizes,
  FUN = "/"
)

detection_counts_by_cluster <- as.matrix(
  (
    curated_counts > 0
  ) %*%
    cluster_indicator
)

detection_fraction_by_cluster <- sweep(
  detection_counts_by_cluster,
  MARGIN = 2L,
  STATS = cluster_sizes,
  FUN = "/"
)

expression_z_by_gene <- safe_matrix_zscore_by_gene(
  mean_expression_by_cluster
)

panel_score_matrix <- matrix(
  NA_real_,
  nrow = length(
    marker_panels_present
  ),
  ncol = length(
    cluster_levels
  ),
  dimnames = list(
    names(
      marker_panels_present
    ),
    cluster_levels
  )
)

panel_mean_detection_matrix <- panel_score_matrix
panel_support_matrix <- panel_score_matrix

for (panel_name in names(
  marker_panels_present
)) {
  panel_features <- marker_panels_present[[panel_name]]

  panel_score_matrix[
    panel_name,
  ] <- colMeans(
    expression_z_by_gene[
      panel_features,
      ,
      drop = FALSE
    ],
    na.rm = TRUE
  )

  panel_mean_detection_matrix[
    panel_name,
  ] <- colMeans(
    detection_fraction_by_cluster[
      panel_features,
      ,
      drop = FALSE
    ],
    na.rm = TRUE
  )

  panel_support_matrix[
    panel_name,
  ] <- colMeans(
    detection_fraction_by_cluster[
      panel_features,
      ,
      drop = FALSE
    ] >= 0.05,
    na.rm = TRUE
  )
}

# -----------------------------------------------------------------------------
# 8. Review-only candidate annotation table
# -----------------------------------------------------------------------------
cluster_review_rows <- vector(
  mode = "list",
  length = length(
    cluster_levels
  )
)

names(
  cluster_review_rows
) <- cluster_levels

for (cluster_name in cluster_levels) {
  scores <- panel_score_matrix[
    ,
    cluster_name
  ]

  score_order <- order(
    scores,
    decreasing = TRUE,
    na.last = TRUE
  )

  top_panel <- names(
    scores
  )[
    score_order[
      1L
    ]
  ]

  second_panel <- names(
    scores
  )[
    score_order[
      2L
    ]
  ]

  top_score <- unname(
    scores[
      top_panel
    ]
  )

  second_score <- unname(
    scores[
      second_panel
    ]
  )

  score_margin <- top_score -
    second_score

  top_support <- panel_support_matrix[
    top_panel,
    cluster_name
  ]

  top_mean_detection <-
    panel_mean_detection_matrix[
      top_panel,
      cluster_name
    ]

  confidence <- if (
    is.finite(
      score_margin
    ) &&
      score_margin >=
        HIGH_CONFIDENCE_MARGIN &&
      top_support >=
        HIGH_CONFIDENCE_SUPPORT &&
      top_mean_detection >=
        HIGH_CONFIDENCE_MEAN_DETECTION
  ) {
    "HIGH"
  } else if (
    is.finite(
      score_margin
    ) &&
      score_margin >=
        MODERATE_CONFIDENCE_MARGIN &&
      top_support >=
        MODERATE_CONFIDENCE_SUPPORT &&
      top_mean_detection >=
        MODERATE_CONFIDENCE_MEAN_DETECTION
  ) {
    "MODERATE"
  } else {
    "AMBIGUOUS"
  }

  resident_score <- panel_score_matrix[
    "Resident_microglia",
    cluster_name
  ]

  infiltrating_scores <- panel_score_matrix[
    setdiff(
      myeloid_panel_names,
      "Resident_microglia"
    ),
    cluster_name
  ]

  top_infiltrating_name <- names(
    infiltrating_scores
  )[
    which.max(
      infiltrating_scores
    )
  ]

  top_infiltrating_score <- max(
    infiltrating_scores,
    na.rm = TRUE
  )

  resident_minus_infiltrating <-
    resident_score -
    top_infiltrating_score

  myeloid_max_score <- max(
    panel_score_matrix[
      myeloid_panel_names,
      cluster_name
    ],
    na.rm = TRUE
  )

  myeloid_review_class <- if (
    top_panel ==
      "Resident_microglia" &&
      resident_minus_infiltrating >=
        MYELOID_CLASS_MARGIN
  ) {
    "Resident_microglia_candidate"
  } else if (
    top_panel ==
      "Resident_microglia" &&
      resident_minus_infiltrating <
        MYELOID_CLASS_MARGIN
  ) {
    "Ambiguous_resident_vs_infiltrating"
  } else if (
    top_panel %in%
      setdiff(
        myeloid_panel_names,
        "Resident_microglia"
      ) &&
      -resident_minus_infiltrating >=
        MYELOID_CLASS_MARGIN
  ) {
    paste0(
      top_panel,
      "_candidate"
    )
  } else if (
    top_panel %in%
      myeloid_panel_names
  ) {
    "Ambiguous_myeloid"
  } else if (
    myeloid_max_score >=
      POSSIBLE_MYELOID_SCORE
  ) {
    "Possible_myeloid_review"
  } else {
    "Not_primary_myeloid"
  }

  cluster_review_rows[[cluster_name]] <- data.frame(
    cluster = cluster_name,
    n_cells = cluster_sizes[
      cluster_name
    ],
    top_panel = top_panel,
    top_panel_score = top_score,
    second_panel = second_panel,
    second_panel_score = second_score,
    score_margin = score_margin,
    top_panel_support_fraction = top_support,
    top_panel_mean_detection =
      top_mean_detection,
    candidate_confidence = confidence,
    resident_microglia_score =
      resident_score,
    top_infiltrating_panel =
      top_infiltrating_name,
    top_infiltrating_score =
      top_infiltrating_score,
    resident_minus_infiltrating_score =
      resident_minus_infiltrating,
    myeloid_review_class =
      myeloid_review_class,
    stringsAsFactors = FALSE
  )
}

cluster_review_table <- do.call(
  rbind,
  cluster_review_rows
)

rownames(
  cluster_review_table
) <- NULL

# Add every panel score to the review table.
for (panel_name in rownames(
  panel_score_matrix
)) {
  output_column <- paste0(
    "score_",
    panel_name
  )

  cluster_review_table[
    ,
    output_column
  ] <- panel_score_matrix[
    panel_name,
    cluster_review_table$cluster
  ]
}

# -----------------------------------------------------------------------------
# 9. Positive markers for all frozen clusters
# -----------------------------------------------------------------------------
log_message(
  "Running FindAllMarkers on downsampled cells per frozen cluster."
)

set.seed(SEED)

all_markers <- Seurat::FindAllMarkers(
  object = obj,
  assay = "RNA",
  slot = "data",
  only.pos = TRUE,
  test.use = "wilcox",
  logfc.threshold =
    MARKER_LOGFC_THRESHOLD,
  min.pct = MARKER_MIN_PCT,
  min.diff.pct =
    MARKER_MIN_DIFF_PCT,
  return.thresh =
    MARKER_RETURN_P_ADJ,
  max.cells.per.ident =
    MAX_CELLS_PER_CLUSTER_FOR_MARKERS,
  random.seed = SEED,
  densify = FALSE,
  verbose = TRUE
)

if (
  is.null(
    all_markers
  ) ||
    nrow(
      all_markers
    ) == 0L
) {
  stop_with_report(
    "FindAllMarkers returned no positive markers."
  )
}

if (!all(
  c(
    "cluster",
    "gene"
  ) %in%
    names(
      all_markers
    )
)) {
  stop_with_report(
    "FindAllMarkers output lacks cluster or gene columns."
  )
}

fold_change_column <- intersect(
  c(
    "avg_log2FC",
    "avg_logFC",
    "avg_diff"
  ),
  names(
    all_markers
  )
)

if (length(
  fold_change_column
) == 0L) {
  stop_with_report(
    "Could not identify the fold-change column in FindAllMarkers output."
  )
}

fold_change_column <- fold_change_column[[1L]]

all_markers$cluster <- as.character(
  all_markers$cluster
)

all_markers_path <- file.path(
  tables_dir,
  "03_all_positive_cluster_markers.csv.gz"
)

data.table::fwrite(
  all_markers,
  file = all_markers_path,
  sep = ",",
  quote = TRUE,
  na = "NA",
  compress = "gzip"
)

log_message(
  "Wrote: ",
  all_markers_path
)

marker_table <- data.table::as.data.table(
  all_markers
)

marker_table[
  ,
  marker_is_mt_ribo_hb := grepl(
    "^(mt-|Rpl|Rps|Hb[ab])",
    gene,
    ignore.case = TRUE
  )
]

marker_table_rankable <- marker_table[
  marker_is_mt_ribo_hb == FALSE
]

data.table::setorderv(
  marker_table_rankable,
  cols = c(
    "cluster",
    fold_change_column
  ),
  order = c(
    1L,
    -1L
  ),
  na.last = TRUE
)

top_markers <- marker_table_rankable[
  ,
  head(
    .SD,
    TOP_MARKERS_PER_CLUSTER
  ),
  by = cluster
]

top_markers <- as.data.frame(
  top_markers
)

write_csv_safe(
  top_markers,
  "04_top30_positive_markers_by_cluster.csv"
)

marker_cluster_coverage <- sort(
  unique(
    as.character(
      top_markers$cluster
    )
  )
)

if (!setequal(
  marker_cluster_coverage,
  cluster_levels
)) {
  missing_marker_clusters <- setdiff(
    cluster_levels,
    marker_cluster_coverage
  )

  stop_with_report(
    paste0(
      "Top-marker output is missing frozen cluster(s): ",
      paste(
        missing_marker_clusters,
        collapse = ", "
      ),
      "."
    )
  )
}

top_marker_strings <- data.table::as.data.table(
  top_markers
)[
  order(
    cluster,
    -get(
      fold_change_column
    )
  ),
  .(
    top_markers = paste(
      head(
        gene,
        10L
      ),
      collapse = ";"
    )
  ),
  by = cluster
]

cluster_review_table <- merge(
  cluster_review_table,
  as.data.frame(
    top_marker_strings
  ),
  by = "cluster",
  all.x = TRUE,
  sort = FALSE
)

# -----------------------------------------------------------------------------
# 10. Author-label supporting cross-check
# -----------------------------------------------------------------------------
author_label_columns <- intersect(
  c(
    "layer1",
    "layer2",
    "layer3",
    "layer4",
    "layer5"
  ),
  names(
    obj@meta.data
  )
)

for (author_column in author_label_columns) {
  dominant_table <- dominant_value_by_cluster(
    metadata = obj@meta.data,
    cluster_column =
      FROZEN_CLUSTER_COLUMN,
    value_column =
      author_column
  )

  names(
    dominant_table
  )[
    names(
      dominant_table
    ) == "dominant_value"
  ] <- paste0(
    "author_",
    author_column,
    "_dominant"
  )

  names(
    dominant_table
  )[
    names(
      dominant_table
    ) == "dominant_fraction"
  ] <- paste0(
    "author_",
    author_column,
    "_fraction"
  )

  cluster_review_table <- merge(
    cluster_review_table,
    dominant_table,
    by = "cluster",
    all.x = TRUE,
    sort = FALSE
  )
}

cluster_review_table <- cluster_review_table[
  match(
    cluster_levels,
    cluster_review_table$cluster
  ),
  ,
  drop = FALSE
]

write_csv_safe(
  cluster_review_table,
  "05_cluster_annotation_review_table.csv"
)

# -----------------------------------------------------------------------------
# 11. Add review-only metadata and save evidence object
# -----------------------------------------------------------------------------
candidate_map <- setNames(
  cluster_review_table$top_panel,
  cluster_review_table$cluster
)

confidence_map <- setNames(
  cluster_review_table$candidate_confidence,
  cluster_review_table$cluster
)

myeloid_map <- setNames(
  cluster_review_table$myeloid_review_class,
  cluster_review_table$cluster
)

cell_cluster_values <- as.character(
  obj@meta.data[
    ,
    FROZEN_CLUSTER_COLUMN
  ]
)

obj$major_lineage_candidate <- unname(
  candidate_map[
    cell_cluster_values
  ]
)

obj$major_lineage_candidate_confidence <- unname(
  confidence_map[
    cell_cluster_values
  ]
)

obj$myeloid_review_class <- unname(
  myeloid_map[
    cell_cluster_values
  ]
)

obj$major_lineage_annotation_status <-
  "REVIEW_ONLY_NOT_FROZEN"

if (
  anyNA(
    obj$major_lineage_candidate
  ) ||
    anyNA(
      obj$major_lineage_candidate_confidence
    ) ||
    anyNA(
      obj$myeloid_review_class
    )
) {
  stop_with_report(
    "Review-only candidate labels were not mapped to every nucleus."
  )
}

obj@misc$GSE234774_step10D <- list(
  input_object = normalizePath(
    input_object_path,
    winslash = "/",
    mustWork = TRUE
  ),
  frozen_reduction =
    FROZEN_REDUCTION,
  frozen_umap =
    FROZEN_UMAP,
  frozen_cluster_column =
    FROZEN_CLUSTER_COLUMN,
  frozen_resolution = 0.2,
  cluster_count =
    EXPECTED_N_CLUSTERS,
  marker_assay = "RNA",
  marker_layer =
    data_layer_name,
  marker_test = "wilcox",
  max_cells_per_cluster =
    MAX_CELLS_PER_CLUSTER_FOR_MARKERS,
  marker_logfc_threshold =
    MARKER_LOGFC_THRESHOLD,
  candidate_labels_frozen =
    FALSE,
  resident_microglia_frozen =
    FALSE,
  author_labels_used_for_assignment =
    FALSE,
  harmony_used_for_expression =
    FALSE,
  arm1_arm2_scoring_run =
    FALSE,
  completed_at =
    as.character(
      Sys.time()
    )
)

# -----------------------------------------------------------------------------
# 12. Figures
# -----------------------------------------------------------------------------
log_message(
  "Creating major-lineage marker figures."
)

# Named marker lists are accepted by Seurat DotPlot and create grouped panels.
neural_glial_requested <- marker_panels_present[
  c(
    "Astrocyte",
    "Oligodendrocyte",
    "OPC",
    "Ependymal",
    "Neuron",
    "Schwann_cell"
  )
]

vascular_stromal_requested <- marker_panels_present[
  c(
    "Endothelial",
    "Pericyte",
    "Vascular_smooth_muscle",
    "Fibroblast"
  )
]

immune_requested <- marker_panels_present[
  c(
    "Resident_microglia",
    "Monocyte_macrophage",
    "Border_associated_macrophage",
    "Neutrophil",
    "Dendritic",
    "T_cell",
    "NK_cell",
    "B_cell"
  )
]

neural_glial_plot_panels <- deduplicate_plot_panels(
  neural_glial_requested
)

vascular_stromal_plot_panels <- deduplicate_plot_panels(
  vascular_stromal_requested
)

immune_plot_panels <- deduplicate_plot_panels(
  immune_requested
)

major_dotplot_path <- file.path(
  figures_dir,
  "01_major_lineage_marker_DotPlots.pdf"
)

grDevices::pdf(
  major_dotplot_path,
  width = 20,
  height = 10,
  onefile = TRUE
)

print(
  Seurat::DotPlot(
    object = obj,
    features =
      neural_glial_plot_panels,
    assay = "RNA",
    group.by =
      FROZEN_CLUSTER_COLUMN,
    cluster.idents = FALSE,
    dot.scale = 6
  ) +
    Seurat::RotatedAxis() +
    ggplot2::labs(
      title = "Neural and glial marker evidence",
      x = NULL,
      y = "Frozen Harmony cluster"
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        size = 8
      )
    )
)

print(
  Seurat::DotPlot(
    object = obj,
    features =
      vascular_stromal_plot_panels,
    assay = "RNA",
    group.by =
      FROZEN_CLUSTER_COLUMN,
    cluster.idents = FALSE,
    dot.scale = 6
  ) +
    Seurat::RotatedAxis() +
    ggplot2::labs(
      title = "Vascular and stromal marker evidence",
      x = NULL,
      y = "Frozen Harmony cluster"
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        size = 8
      )
    )
)

print(
  Seurat::DotPlot(
    object = obj,
    features =
      immune_plot_panels,
    assay = "RNA",
    group.by =
      FROZEN_CLUSTER_COLUMN,
    cluster.idents = FALSE,
    dot.scale = 6
  ) +
    Seurat::RotatedAxis() +
    ggplot2::labs(
      title = "Immune and myeloid marker evidence",
      x = NULL,
      y = "Frozen Harmony cluster"
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        size = 8
      )
    )
)

grDevices::dev.off()

log_message(
  "Wrote: ",
  major_dotplot_path
)

candidate_levels <- unique(
  cluster_review_table$top_panel
)

obj$major_lineage_candidate <- factor(
  obj$major_lineage_candidate,
  levels = sort(
    candidate_levels
  )
)

obj$major_lineage_candidate_confidence <- factor(
  obj$major_lineage_candidate_confidence,
  levels = c(
    "HIGH",
    "MODERATE",
    "AMBIGUOUS"
  )
)

p_cluster <- Seurat::DimPlot(
  object = obj,
  reduction = FROZEN_UMAP,
  group.by =
    FROZEN_CLUSTER_COLUMN,
  label = TRUE,
  repel = TRUE,
  raster = TRUE
) +
  ggplot2::labs(
    title = "Frozen conservative-Harmony clusters",
    subtitle = "Resolution 0.2; 25 clusters"
  ) +
  Seurat::NoLegend()

p_candidate <- Seurat::DimPlot(
  object = obj,
  reduction = FROZEN_UMAP,
  group.by =
    "major_lineage_candidate",
  label = FALSE,
  raster = TRUE
) +
  ggplot2::labs(
    title = "Review-only top marker-panel candidate",
    subtitle = "Not a frozen biological annotation"
  )

p_confidence <- Seurat::DimPlot(
  object = obj,
  reduction = FROZEN_UMAP,
  group.by =
    "major_lineage_candidate_confidence",
  label = FALSE,
  raster = TRUE
) +
  ggplot2::labs(
    title = "Automatic panel-score confidence",
    subtitle = "Ambiguous clusters require manual marker review"
  )

p_myeloid_class <- Seurat::DimPlot(
  object = obj,
  reduction = FROZEN_UMAP,
  group.by =
    "myeloid_review_class",
  label = FALSE,
  raster = TRUE
) +
  ggplot2::labs(
    title = "Resident versus infiltrating-myeloid review class",
    subtitle = "Screening aid only"
  )

umap_review_plot <- (
  p_cluster |
    p_candidate
) / (
  p_confidence |
    p_myeloid_class
)

ggplot2::ggsave(
  filename = file.path(
    figures_dir,
    "02_UMAP_candidate_annotation_review.pdf"
  ),
  plot = umap_review_plot,
  width = 18,
  height = 14,
  limitsize = FALSE
)

# Focused myeloid plot.
myeloid_plot_panels <- c(
  list(
    Myeloid_general =
      supplementary_panels_present[["Myeloid_general"]]
  ),
  marker_panels_present[
    myeloid_panel_names
  ],
  list(
    Activated_phagolysosomal =
      supplementary_panels_present[["Activated_phagolysosomal"]]
  )
)

myeloid_plot_panels <- deduplicate_plot_panels(
  myeloid_plot_panels
)

myeloid_feature_requested <- c(
  "P2ry12",
  "Tmem119",
  "Sall1",
  "Ccr2",
  "Ly6c2",
  "Lgals3",
  "Mrc1",
  "Lyve1"
)

myeloid_feature_resolved <- feature_resolution[
  myeloid_feature_requested
]

myeloid_feature_resolved <- unique(
  unname(
    myeloid_feature_resolved[
      !is.na(
        myeloid_feature_resolved
      )
    ]
  )
)

if (length(
  myeloid_feature_resolved
) < 6L) {
  stop_with_report(
    "Fewer than six focused myeloid FeaturePlot genes were resolved."
  )
}

myeloid_pdf_path <- file.path(
  figures_dir,
  "03_myeloid_resident_infiltrating_marker_review.pdf"
)

grDevices::pdf(
  myeloid_pdf_path,
  width = 20,
  height = 11,
  onefile = TRUE
)

print(
  Seurat::DotPlot(
    object = obj,
    features =
      myeloid_plot_panels,
    assay = "RNA",
    group.by =
      FROZEN_CLUSTER_COLUMN,
    cluster.idents = FALSE,
    dot.scale = 6
  ) +
    Seurat::RotatedAxis() +
    ggplot2::labs(
      title = paste(
        "Resident microglia versus infiltrating/BAM/neutrophil/DC",
        "marker evidence"
      ),
      x = NULL,
      y = "Frozen Harmony cluster"
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        size = 8
      )
    )
)

print(
  Seurat::FeaturePlot(
    object = obj,
    features =
      myeloid_feature_resolved,
    reduction =
      FROZEN_UMAP,
    order = TRUE,
    raster = TRUE,
    ncol = 4
  ) &
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        size = 10
      )
    )
)

grDevices::dev.off()

log_message(
  "Wrote: ",
  myeloid_pdf_path
)

# -----------------------------------------------------------------------------
# 13. Save and reload evidence object
# -----------------------------------------------------------------------------
log_message(
  "Saving Step 10D annotation-evidence object."
)

saveRDS(
  obj,
  output_object_path,
  compress = FALSE
)

reload_obj <- tryCatch(
  readRDS(
    output_object_path
  ),
  error = function(error_object) {
    error_object
  }
)

reload_pass <-
  !inherits(
    reload_obj,
    "error"
  ) &&
  inherits(
    reload_obj,
    "Seurat"
  ) &&
  nrow(
    reload_obj
  ) == nrow(
    obj
  ) &&
  ncol(
    reload_obj
  ) == ncol(
    obj
  ) &&
  all(
    c(
      "major_lineage_candidate",
      "major_lineage_candidate_confidence",
      "myeloid_review_class",
      "major_lineage_annotation_status"
    ) %in%
      names(
        reload_obj@meta.data
      )
  ) &&
  all(
    c(
      FROZEN_REDUCTION,
      FROZEN_UMAP
    ) %in%
      names(
        reload_obj@reductions
      )
  )

rm(
  reload_obj
)

invisible(
  gc()
)

# -----------------------------------------------------------------------------
# 14. Audit and decision
# -----------------------------------------------------------------------------
n_ambiguous_clusters <- sum(
  cluster_review_table$candidate_confidence ==
    "AMBIGUOUS"
)

resident_candidate_clusters <- cluster_review_table$cluster[
  cluster_review_table$myeloid_review_class ==
    "Resident_microglia_candidate"
]

ambiguous_resident_clusters <- cluster_review_table$cluster[
  cluster_review_table$myeloid_review_class ==
    "Ambiguous_resident_vs_infiltrating"
]

infiltrating_candidate_clusters <- cluster_review_table$cluster[
  grepl(
    "_candidate$",
    cluster_review_table$myeloid_review_class
  ) &
    cluster_review_table$myeloid_review_class !=
      "Resident_microglia_candidate"
]

audit_checks <- data.frame(
  check = c(
    "Input contains 26,918 features and 144,359 singlets",
    "Frozen conservative-Harmony reduction and UMAP are present",
    "Frozen hcon_res_0_2 contains exactly 25 clusters",
    "RNA LogNormalize data, not Harmony embeddings, were used for marker testing",
    "Every curated candidate panel contains at least two resolved genes",
    "FindAllMarkers returned positive markers for all 25 clusters",
    "Cluster annotation review table contains all 25 clusters",
    "Every nucleus received a review-only candidate, confidence, and myeloid class",
    "All three required marker-review PDFs were created",
    "Saved evidence object reloads with required metadata and reductions",
    "Author labels were used to assign candidate labels",
    "Final major-lineage annotations have been frozen",
    "A clean resident-microglia object has been frozen",
    "GSE162610 Arm_1/Arm_2 scoring was run"
  ),
  required_for_pass = c(
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    FALSE,
    FALSE,
    FALSE,
    FALSE
  ),
  status = c(
    ifelse(
      nrow(obj) == 26918L &&
        ncol(obj) == 144359L,
      "PASS",
      "STOP"
    ),
    ifelse(
      all(
        c(
          FROZEN_REDUCTION,
          FROZEN_UMAP
        ) %in%
          names(
            obj@reductions
          )
      ),
      "PASS",
      "STOP"
    ),
    ifelse(
      length(
        cluster_levels
      ) ==
        EXPECTED_N_CLUSTERS,
      "PASS",
      "STOP"
    ),
    "PASS",
    ifelse(
      all(
        panel_availability$present_n >=
          2L
      ),
      "PASS",
      "STOP"
    ),
    ifelse(
      setequal(
        marker_cluster_coverage,
        cluster_levels
      ),
      "PASS",
      "STOP"
    ),
    ifelse(
      nrow(
        cluster_review_table
      ) ==
        EXPECTED_N_CLUSTERS,
      "PASS",
      "STOP"
    ),
    ifelse(
      !anyNA(
        obj$major_lineage_candidate
      ) &&
        !anyNA(
          obj$major_lineage_candidate_confidence
        ) &&
        !anyNA(
          obj$myeloid_review_class
        ),
      "PASS",
      "STOP"
    ),
    ifelse(
      all(
        file.exists(
          file.path(
            figures_dir,
            c(
              "01_major_lineage_marker_DotPlots.pdf",
              "02_UMAP_candidate_annotation_review.pdf",
              "03_myeloid_resident_infiltrating_marker_review.pdf"
            )
          )
        )
      ),
      "PASS",
      "STOP"
    ),
    ifelse(
      reload_pass,
      "PASS",
      "STOP"
    ),
    "CONFIRMED_NOT_USED_FOR_ASSIGNMENT",
    "NOT_FROZEN",
    "NOT_FROZEN",
    "CONFIRMED_NOT_RUN"
  ),
  detail = c(
    paste0(
      nrow(obj),
      " features x ",
      ncol(obj),
      " singlets"
    ),
    paste(
      c(
        FROZEN_REDUCTION,
        FROZEN_UMAP
      ),
      collapse = ", "
    ),
    paste0(
      "clusters=",
      length(
        cluster_levels
      )
    ),
    paste0(
      "assay=RNA; layer=",
      data_layer_name,
      "; cluster membership from ",
      FROZEN_CLUSTER_COLUMN
    ),
    paste0(
      "minimum resolved markers per panel=",
      min(
        panel_availability$present_n
      )
    ),
    paste0(
      "marker clusters=",
      length(
        marker_cluster_coverage
      )
    ),
    paste0(
      "review rows=",
      nrow(
        cluster_review_table
      )
    ),
    paste0(
      "candidate labels mapped to ",
      ncol(obj),
      " nuclei"
    ),
    "Three PDFs present.",
    ifelse(
      reload_pass,
      "RDS reload PASS",
      "RDS reload failed"
    ),
    paste0(
      "Author columns available: ",
      ifelse(
        length(
          author_label_columns
        ) == 0L,
        "none",
        paste(
          author_label_columns,
          collapse = ", "
        )
      ),
      ". They were added only as cross-check columns."
    ),
    paste0(
      "Ambiguous automatic candidate clusters=",
      n_ambiguous_clusters
    ),
    paste0(
      "Resident candidates=",
      ifelse(
        length(
          resident_candidate_clusters
        ) == 0L,
        "none",
        paste(
          resident_candidate_clusters,
          collapse = ","
        )
      ),
      "; resident/infiltrating ambiguous=",
      ifelse(
        length(
          ambiguous_resident_clusters
        ) == 0L,
        "none",
        paste(
          ambiguous_resident_clusters,
          collapse = ","
        )
      )
    ),
    "Arm_1/Arm_2 remain frozen from GSE162610 and were not evaluated here."
  ),
  stringsAsFactors = FALSE
)

write_csv_safe(
  audit_checks,
  "09_STEP10D_audit_checks.csv"
)

required_failures <-
  audit_checks$required_for_pass &
  audit_checks$status != "PASS"

final_decision <- if (
  any(
    required_failures
  )
) {
  "STOP"
} else {
  "PASS_WITH_ANNOTATION_REVIEW"
}

# -----------------------------------------------------------------------------
# 15. Final report
# -----------------------------------------------------------------------------
candidate_summary <- data.table::as.data.table(
  cluster_review_table
)[
  ,
  .(
    n_clusters = .N,
    n_cells = sum(
      n_cells
    ),
    high_confidence = sum(
      candidate_confidence == "HIGH"
    ),
    moderate_confidence = sum(
      candidate_confidence == "MODERATE"
    ),
    ambiguous = sum(
      candidate_confidence == "AMBIGUOUS"
    )
  ),
  by = top_panel
]

candidate_summary <- candidate_summary[
  order(
    -n_cells,
    top_panel
  )
]

candidate_lines <- paste0(
  candidate_summary$top_panel,
  ": ",
  candidate_summary$n_clusters,
  " clusters; ",
  candidate_summary$n_cells,
  " nuclei; confidence H/M/A=",
  candidate_summary$high_confidence,
  "/",
  candidate_summary$moderate_confidence,
  "/",
  candidate_summary$ambiguous
)

myeloid_review_summary <- data.table::as.data.table(
  cluster_review_table
)[
  ,
  .(
    clusters = paste(
      cluster,
      collapse = ","
    ),
    n_clusters = .N,
    n_cells = sum(
      n_cells
    )
  ),
  by = myeloid_review_class
]

myeloid_lines <- paste0(
  myeloid_review_summary$myeloid_review_class,
  ": clusters=",
  myeloid_review_summary$clusters,
  "; n_clusters=",
  myeloid_review_summary$n_clusters,
  "; nuclei=",
  myeloid_review_summary$n_cells
)

report_lines <- c(
  "GSE234774 STEP 10D MARKER-FIRST MAJOR-LINEAGE ANNOTATION",
  "============================================================",
  paste0(
    "Run time: ",
    format(
      Sys.time(),
      "%Y-%m-%d %H:%M:%S"
    )
  ),
  paste0(
    "Input object: ",
    input_object_path
  ),
  paste0(
    "Output directory: ",
    output_dir
  ),
  "",
  paste0(
    "FINAL DECISION: ",
    final_decision
  ),
  "",
  "FROZEN CLUSTERING INPUT",
  paste0(
    "Reduction: ",
    FROZEN_REDUCTION
  ),
  paste0(
    "UMAP: ",
    FROZEN_UMAP
  ),
  paste0(
    "Cluster column: ",
    FROZEN_CLUSTER_COLUMN
  ),
  paste0(
    "Clusters: ",
    length(
      cluster_levels
    )
  ),
  "",
  "MARKER DISCOVERY",
  "Assay: RNA",
  paste0(
    "Expression layer: ",
    data_layer_name
  ),
  "Test: Wilcoxon rank-sum; positive markers only",
  paste0(
    "Maximum marker-test cells per cluster: ",
    MAX_CELLS_PER_CLUSTER_FOR_MARKERS
  ),
  paste0(
    "log2FC threshold: ",
    MARKER_LOGFC_THRESHOLD
  ),
  paste0(
    "Minimum detected fraction: ",
    MARKER_MIN_PCT
  ),
  "Harmony embeddings were not used as expression values.",
  "",
  "REVIEW-ONLY CANDIDATE SUMMARY",
  candidate_lines,
  "",
  "MYELOID REVIEW SUMMARY",
  myeloid_lines,
  "",
  "INTERPRETATION BOUNDARIES",
  "1. Candidate labels are generated from cluster marker-panel scores and are not final annotations.",
  "2. Top positive markers, DotPlots, FeaturePlots, UMAP location, and author-label cross-checks must be reviewed together.",
  "3. Activated resident microglia can lose homeostatic markers; low P2ry12/Tmem119 alone is not sufficient to call infiltrating macrophages.",
  "4. Monocyte/macrophage, BAM, neutrophil, and dendritic evidence must be considered before freezing a resident-microglia set.",
  "5. Ambiguous myeloid clusters must remain separate rather than being forced into resident or infiltrating labels.",
  "6. No cluster was removed, and no resident-microglia object was frozen.",
  "7. No GSE162610 Arm_1/Arm_2 signature was scored or redefined.",
  "",
  "PRIMARY EVIDENCE OBJECT",
  output_object_path,
  "",
  "NEXT DECISION AFTER FILE REVIEW",
  "Freeze major-lineage labels cluster by cluster.",
  "Then construct a clean resident-microglia candidate set while retaining an explicit ambiguous-myeloid compartment."
)

if (
  final_decision ==
    "STOP"
) {
  failed <- audit_checks[
    audit_checks$required_for_pass &
      audit_checks$status != "PASS",
    ,
    drop = FALSE
  ]

  report_lines <- c(
    report_lines,
    "",
    "STOP REASONS"
  )

  for (
    index in seq_len(
      nrow(
        failed
      )
    )
  ) {
    report_lines <- c(
      report_lines,
      paste0(
        index,
        ". ",
        failed$check[
          index
        ],
        " -- ",
        failed$detail[
          index
        ]
      )
    )
  }
}

report_path <- file.path(
  output_dir,
  "10_STEP10D_FINAL_REPORT.txt"
)

writeLines(
  report_lines,
  report_path,
  useBytes = TRUE
)

log_message(
  "Final decision: ",
  final_decision
)

log_message(
  "Final report: ",
  report_path
)

cat("\n")
cat("============================================================\n")
cat("GSE234774 STEP 10D COMPLETED\n")
cat(
  "FINAL DECISION:",
  final_decision,
  "\n"
)
cat(
  "Output:",
  output_dir,
  "\n"
)
cat("Primary evidence object:\n")
cat(
  "  ",
  output_object_path,
  "\n",
  sep = ""
)
cat("Upload exactly these seven files for review:\n")
cat("  1) figures/01_major_lineage_marker_DotPlots.pdf\n")
cat("  2) figures/02_UMAP_candidate_annotation_review.pdf\n")
cat("  3) figures/03_myeloid_resident_infiltrating_marker_review.pdf\n")
cat("  4) tables/04_top30_positive_markers_by_cluster.csv\n")
cat("  5) tables/05_cluster_annotation_review_table.csv\n")
cat("  6) tables/09_STEP10D_audit_checks.csv\n")
cat("  7) 10_STEP10D_FINAL_REPORT.txt\n")
cat("============================================================\n")

if (
  final_decision ==
    "STOP"
) {
  stop(
    "GSE234774 Step 10D decision is STOP. Review the final report.",
    call. = FALSE
  )
}
