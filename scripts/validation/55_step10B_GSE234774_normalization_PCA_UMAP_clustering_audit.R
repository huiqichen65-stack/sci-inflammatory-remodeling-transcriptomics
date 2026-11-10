# =============================================================================
# 55_step10B_GSE234774_normalization_PCA_UMAP_clustering_audit.R
#
# Adult SCI activated/DAM-like microglia project
# GSE234774 Step 10B:
#   LogNormalize + HVG + PCA + unintegrated UMAP + clustering audit
#
# INPUT
#   <PROJECT_ROOT>/04_results/
#     step9C_GSE234774_samplewise_scDblFinder/RDS/
#       GSE234774_step9C_QC_singlets_BPCells.rds
#
# OUTPUT
#   <PROJECT_ROOT>/04_results/
#     step10B_GSE234774_normalization_PCA_UMAP_clustering/
#
# FROZEN DESIGN RULES
#   - LogNormalize is the primary normalization.
#   - SCTransform is not run in the primary workflow.
#   - Cell-cycle scores are not regressed.
#   - percent.mt, nCount_RNA, timepoint, and biological_library are not regressed.
#   - Harmony/integration is not run automatically because biological_library
#     is nested within timepoint and correction could remove real time biology.
#   - PCA/UMAP/clustering are used only to construct and audit the whole atlas.
#   - No cell type is named in this script.
#   - Author-provided labels are retained only as supporting metadata.
#   - GSE234774 does not redefine GSE162610 Arm_1 or Arm_2.
#
# NEXT STEP AFTER REVIEW
#   Marker-first separation of resident microglia from infiltrating
#   monocyte/macrophage, with ambiguous myeloid nuclei retained separately.
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

N_VARIABLE_FEATURES <- 3000L
N_PCS_COMPUTE <- 50L
PCS_USE <- 1:40

NORMALIZATION_SCALE_FACTOR <- 10000
NORMALIZATION_BLOCK_SIZE <- 10000L

K_PARAM <- 30L
N_TREES <- 50L
SNN_PRUNE <- 1 / 15

CLUSTER_RESOLUTIONS <- c(
  0.2,
  0.4,
  0.6,
  0.8,
  1.0
)

PRIMARY_RESOLUTION <- 0.6

UMAP_N_NEIGHBORS <- 30L
UMAP_MIN_DIST <- 0.3

# Set TRUE only when the saved normalized BPCells directory is incomplete
# or must be regenerated.
FORCE_REWRITE_NORMALIZED_DATA <- FALSE

# Review flags only; they do not automatically exclude a cluster.
LIBRARY_DOMINANCE_REVIEW_THRESHOLD <- 0.80
TIMEPOINT_DOMINANCE_REVIEW_THRESHOLD <- 0.95
MIN_CLUSTER_SIZE_FOR_DOMINANCE_REVIEW <- 500L

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
# 1. Paths
# -----------------------------------------------------------------------------
project_root <- Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
input_object_path <- file.path(
  project_root,
  "04_results",
  "step9C_GSE234774_samplewise_scDblFinder",
  "RDS",
  "GSE234774_step9C_QC_singlets_BPCells.rds"
)

output_dir <- file.path(
  project_root,
  "04_results",
  "step10B_GSE234774_normalization_PCA_UMAP_clustering"
)

tables_dir <- file.path(output_dir, "tables")
figures_dir <- file.path(output_dir, "figures")
rds_dir <- file.path(output_dir, "RDS")
bpcells_dir <- file.path(output_dir, "BPCells")

for (d in c(
  output_dir,
  tables_dir,
  figures_dir,
  rds_dir,
  bpcells_dir
)) {
  dir.create(
    d,
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

normalized_data_dir <- file.path(
  bpcells_dir,
  "GSE234774_step10B_LogNormalize_data"
)

output_object_path <- file.path(
  rds_dir,
  "GSE234774_step10B_LogNormalize_PCA_UMAP_clusters_BPCells.rds"
)

# -----------------------------------------------------------------------------
# 2. Helpers
# -----------------------------------------------------------------------------
log_message <- function(...) {
  msg <- paste0(...)
  stamp <- format(
    Sys.time(),
    "%Y-%m-%d %H:%M:%S"
  )
  line <- paste0(
    "[",
    stamp,
    "] ",
    msg
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
    "12_STEP10B_FINAL_REPORT.txt"
  )

  report <- c(
    "GSE234774 STEP 10B NORMALIZATION / PCA / UMAP / CLUSTERING",
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
    "Do not start resident-microglia identification until this error is resolved."
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

safe_package_version <- function(pkg) {
  if (!requireNamespace(
    pkg,
    quietly = TRUE
  )) {
    return(NA_character_)
  }

  as.character(
    utils::packageVersion(pkg)
  )
}

open_bpcells_safe <- function(path) {
  tryCatch(
    BPCells::open_matrix_dir(path),
    error = function(e) e
  )
}

resolution_column_name <- function(resolution) {
  paste0(
    "cluster_res_",
    gsub(
      "\\.",
      "_",
      format(
        resolution,
        scientific = FALSE,
        trim = TRUE
      )
    )
  )
}

safe_spearman <- function(
    x,
    y
) {
  ok <- is.finite(x) &
    is.finite(y)

  if (sum(ok) < 3L) {
    return(NA_real_)
  }

  suppressWarnings(
    stats::cor(
      x[ok],
      y[ok],
      method = "spearman"
    )
  )
}

collapse_unique <- function(
    x,
    max_values = 20L
) {
  x <- unique(
    as.character(x)
  )

  x <- x[
    !is.na(x) &
      nzchar(x)
  ]

  if (length(x) == 0L) {
    return(NA_character_)
  }

  if (length(x) > max_values) {
    x <- c(
      x[seq_len(max_values)],
      paste0(
        "... +",
        length(x) - max_values,
        " more"
      )
    )
  }

  paste(
    x,
    collapse = " | "
  )
}

# -----------------------------------------------------------------------------
# 3. Package preflight
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
  "SeuratObject",
  "BPCells"
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

if (
  utils::packageVersion("Seurat") <
    "5.2.0"
) {
  stop_with_report(
    paste0(
      "Seurat >=5.2.0 is required for the improved BPCells PCA path. ",
      "Detected Seurat ",
      utils::packageVersion("Seurat"),
      ". Update Seurat and rerun."
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
  library(BPCells)
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
    "SeuratObject",
    "BPCells"
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
    safe_package_version("SeuratObject"),
    safe_package_version("BPCells")
  ),
  stringsAsFactors = FALSE
)

write_csv_safe(
  package_versions,
  "00_package_versions.csv"
)

# -----------------------------------------------------------------------------
# 4. Load and audit the singlet object
# -----------------------------------------------------------------------------
if (!file.exists(input_object_path)) {
  stop_with_report(
    paste0(
      "Input singlet object not found: ",
      input_object_path
    )
  )
}

log_message(
  "Loading Step 9C singlet-only object."
)

obj <- readRDS(
  input_object_path
)

if (!inherits(
  obj,
  "Seurat"
)) {
  stop_with_report(
    "The input RDS is not a Seurat object."
  )
}

if (!"RNA" %in% names(obj@assays)) {
  stop_with_report(
    "RNA assay is absent."
  )
}

Seurat::DefaultAssay(obj) <- "RNA"

required_metadata <- c(
  "biological_library",
  "timepoint",
  "replicate_id",
  "nCount_RNA",
  "nFeature_RNA",
  "percent.mt",
  "percent.HB",
  "scDblFinder.score",
  "scDblFinder.class"
)

missing_metadata <- setdiff(
  required_metadata,
  names(obj@meta.data)
)

if (length(missing_metadata) > 0L) {
  stop_with_report(
    paste0(
      "Input object is missing metadata column(s): ",
      paste(
        missing_metadata,
        collapse = ", "
      ),
      "."
    )
  )
}

if (
  ncol(obj) != 144359L
) {
  stop_with_report(
    paste0(
      "Expected 144,359 singlets, detected ",
      ncol(obj),
      "."
    )
  )
}

if (
  nrow(obj) != 26918L
) {
  stop_with_report(
    paste0(
      "Expected 26,918 features, detected ",
      nrow(obj),
      "."
    )
  )
}

if (
  !all(
    as.character(
      obj$scDblFinder.class
    ) == "singlet"
  )
) {
  stop_with_report(
    "The input object contains non-singlet nuclei."
  )
}

n_libraries <- length(
  unique(
    as.character(
      obj$biological_library
    )
  )
)

if (n_libraries != 20L) {
  stop_with_report(
    paste0(
      "Expected 20 biological libraries, detected ",
      n_libraries,
      "."
    )
  )
}

observed_timepoints <- unique(
  as.character(
    obj$timepoint
  )
)

if (
  !all(
    expected_timepoints %in%
      observed_timepoints
  )
) {
  stop_with_report(
    "The input object does not contain all expected timepoints."
  )
}

layer_names_before <- SeuratObject::Layers(
  object = obj,
  assay = "RNA",
  search = NA
)

counts_layers <- layer_names_before[
  grepl(
    "^counts($|\\.)",
    layer_names_before
  )
]

if (length(counts_layers) != 1L) {
  stop_with_report(
    paste0(
      "Expected exactly one RNA counts layer, detected: ",
      paste(
        counts_layers,
        collapse = ", "
      )
    )
  )
}

counts_layer <- counts_layers[[1L]]

counts_mat <- SeuratObject::LayerData(
  object = obj,
  assay = "RNA",
  layer = counts_layer,
  fast = FALSE
)

if (
  !identical(
    colnames(counts_mat),
    colnames(obj)
  )
) {
  stop_with_report(
    "Counts-layer barcodes do not match the Seurat object."
  )
}

input_audit <- data.frame(
  item = c(
    "Input object",
    "Features",
    "Singlet nuclei",
    "Biological libraries",
    "Timepoints",
    "Counts layer",
    "All nuclei are singlets",
    "Normalization present before Step 10B",
    "PCA present before Step 10B",
    "Harmony present before Step 10B"
  ),
  value = c(
    normalizePath(
      input_object_path,
      winslash = "/",
      mustWork = TRUE
    ),
    nrow(obj),
    ncol(obj),
    n_libraries,
    length(observed_timepoints),
    counts_layer,
    all(
      as.character(
        obj$scDblFinder.class
      ) == "singlet"
    ),
    any(
      grepl(
        "^data($|\\.)",
        layer_names_before
      )
    ),
    "pca" %in% names(
      obj@reductions
    ),
    "harmony" %in% names(
      obj@reductions
    )
  ),
  stringsAsFactors = FALSE
)

write_csv_safe(
  input_audit,
  "01_input_object_audit.csv"
)

# -----------------------------------------------------------------------------
# 5. LogNormalize
# -----------------------------------------------------------------------------
log_message(
  "Running LogNormalize on the RNA assay."
)

set.seed(SEED)

obj <- Seurat::NormalizeData(
  object = obj,
  assay = "RNA",
  normalization.method = "LogNormalize",
  scale.factor = NORMALIZATION_SCALE_FACTOR,
  block.size = NORMALIZATION_BLOCK_SIZE,
  verbose = TRUE
)

data_layers <- SeuratObject::Layers(
  object = obj,
  assay = "RNA",
  search = "^data"
)

if (length(data_layers) != 1L) {
  stop_with_report(
    paste0(
      "Expected exactly one normalized RNA data layer after NormalizeData, detected: ",
      paste(
        data_layers,
        collapse = ", "
      )
    )
  )
}

normalized_layer_name <- data_layers[[1L]]

normalized_layer <- SeuratObject::LayerData(
  object = obj,
  assay = "RNA",
  layer = normalized_layer_name,
  fast = FALSE
)

if (
  isTRUE(
    FORCE_REWRITE_NORMALIZED_DATA
  ) &&
    dir.exists(normalized_data_dir)
) {
  log_message(
    "Deleting the existing normalized BPCells directory."
  )

  unlink(
    normalized_data_dir,
    recursive = TRUE,
    force = TRUE
  )
}

if (!dir.exists(normalized_data_dir)) {
  log_message(
    "Writing LogNormalize data to a stable BPCells directory."
  )

  BPCells::write_matrix_dir(
    mat = normalized_layer,
    dir = normalized_data_dir,
    compress = TRUE,
    overwrite = FALSE
  )

  log_message(
    "Completed normalized BPCells directory."
  )
} else {
  log_message(
    "Reusing existing normalized BPCells directory: ",
    normalized_data_dir
  )
}

normalized_disk <- open_bpcells_safe(
  normalized_data_dir
)

if (inherits(
  normalized_disk,
  "error"
)) {
  stop_with_report(
    paste0(
      "The normalized BPCells directory cannot be opened: ",
      conditionMessage(normalized_disk),
      ". Set FORCE_REWRITE_NORMALIZED_DATA <- TRUE and rerun."
    )
  )
}

if (
  !identical(
    as.integer(
      dim(normalized_disk)
    ),
    as.integer(
      dim(obj)
    )
  ) ||
    !identical(
      rownames(normalized_disk),
      rownames(obj)
    ) ||
    !identical(
      colnames(normalized_disk),
      colnames(obj)
    )
) {
  stop_with_report(
    paste(
      "The normalized BPCells matrix does not match the current singlet object.",
      "Set FORCE_REWRITE_NORMALIZED_DATA <- TRUE and rerun."
    )
  )
}

# Replace the normalized layer with its stable on-disk copy.
obj[["RNA"]]$data <- normalized_disk

rm(
  normalized_layer
)

invisible(gc())

normalization_audit <- data.frame(
  item = c(
    "Normalization method",
    "Scale factor",
    "Block size",
    "Normalized data storage",
    "Normalized features",
    "Normalized nuclei",
    "SCTransform run",
    "Cell-cycle regression run",
    "Technical-covariate regression run"
  ),
  value = c(
    "LogNormalize",
    NORMALIZATION_SCALE_FACTOR,
    NORMALIZATION_BLOCK_SIZE,
    normalizePath(
      normalized_data_dir,
      winslash = "/",
      mustWork = TRUE
    ),
    nrow(normalized_disk),
    ncol(normalized_disk),
    FALSE,
    FALSE,
    FALSE
  ),
  stringsAsFactors = FALSE
)

write_csv_safe(
  normalization_audit,
  "02_normalization_audit.csv"
)

# -----------------------------------------------------------------------------
# 6. Variable features
# -----------------------------------------------------------------------------
log_message(
  "Selecting variable features."
)

obj <- Seurat::FindVariableFeatures(
  object = obj,
  assay = "RNA",
  selection.method = "vst",
  nfeatures = N_VARIABLE_FEATURES,
  verbose = TRUE
)

variable_features <- Seurat::VariableFeatures(
  object = obj,
  assay = "RNA"
)

if (
  length(variable_features) !=
    N_VARIABLE_FEATURES
) {
  stop_with_report(
    paste0(
      "Expected ",
      N_VARIABLE_FEATURES,
      " variable features, detected ",
      length(variable_features),
      "."
    )
  )
}

mt_hvg <- variable_features[
  grepl(
    "^mt-",
    variable_features,
    ignore.case = TRUE
  )
]

hb_hvg <- variable_features[
  grepl(
    "^hb[ab]",
    variable_features,
    ignore.case = TRUE
  )
]

pca_features <- setdiff(
  variable_features,
  union(
    mt_hvg,
    hb_hvg
  )
)

if (length(pca_features) < 2000L) {
  stop_with_report(
    paste0(
      "Only ",
      length(pca_features),
      " PCA features remain after excluding mitochondrial and hemoglobin HVGs."
    )
  )
}

variable_feature_table <- data.frame(
  rank = seq_along(
    variable_features
  ),
  gene = variable_features,
  excluded_from_PCA = variable_features %in%
    union(
      mt_hvg,
      hb_hvg
    ),
  exclusion_reason = ifelse(
    variable_features %in% mt_hvg,
    "mitochondrial",
    ifelse(
      variable_features %in% hb_hvg,
      "hemoglobin",
      ""
    )
  ),
  stringsAsFactors = FALSE
)

write_csv_safe(
  variable_feature_table,
  "03_variable_features.csv"
)

# -----------------------------------------------------------------------------
# 7. Scale PCA features and run PCA
# -----------------------------------------------------------------------------
log_message(
  "Scaling PCA features without regressing cell cycle, timepoint, library, or QC covariates."
)

obj <- Seurat::ScaleData(
  object = obj,
  assay = "RNA",
  features = pca_features,
  vars.to.regress = NULL,
  do.center = TRUE,
  do.scale = TRUE,
  block.size = 1000L,
  verbose = TRUE
)

log_message(
  "Running PCA."
)

set.seed(SEED)

obj <- Seurat::RunPCA(
  object = obj,
  assay = "RNA",
  features = pca_features,
  npcs = N_PCS_COMPUTE,
  reduction.name = "pca",
  reduction.key = "PC_",
  approx = TRUE,
  seed.use = SEED,
  verbose = TRUE
)

if (!"pca" %in% names(
  obj@reductions
)) {
  stop_with_report(
    "PCA reduction was not created."
  )
}

pca_embeddings <- Seurat::Embeddings(
  object = obj,
  reduction = "pca"
)

if (
  nrow(pca_embeddings) != ncol(obj) ||
    ncol(pca_embeddings) <
    max(PCS_USE)
) {
  stop_with_report(
    "PCA embedding dimensions are insufficient."
  )
}

pca_stdev <- Seurat::Stdev(
  object = obj,
  reduction = "pca"
)

pca_variance <- pca_stdev^2
pca_percent <- 100 *
  pca_variance /
  sum(pca_variance)

pca_table <- data.frame(
  PC = seq_along(
    pca_stdev
  ),
  stdev = pca_stdev,
  variance = pca_variance,
  percent_variance = pca_percent,
  cumulative_percent_variance = cumsum(
    pca_percent
  ),
  used_for_neighbors_UMAP = seq_along(
    pca_stdev
  ) %in% PCS_USE,
  stringsAsFactors = FALSE
)

write_csv_safe(
  pca_table,
  "04_PCA_variance.csv"
)

loadings <- Seurat::Loadings(
  object = obj,
  reduction = "pca"
)

loading_rows <- list()

for (pc_index in seq_len(
  min(
    10L,
    ncol(loadings)
  )
)) {
  pc_name <- colnames(
    loadings
  )[[pc_index]]

  loading_vector <- loadings[
    ,
    pc_index
  ]

  positive_genes <- names(
    sort(
      loading_vector,
      decreasing = TRUE
    )
  )[seq_len(
    min(
      30L,
      length(loading_vector)
    )
  )]

  negative_genes <- names(
    sort(
      loading_vector,
      decreasing = FALSE
    )
  )[seq_len(
    min(
      30L,
      length(loading_vector)
    )
  )]

  loading_rows[[length(
    loading_rows
  ) + 1L]] <- data.frame(
    PC = pc_name,
    direction = "positive",
    rank = seq_along(
      positive_genes
    ),
    gene = positive_genes,
    loading = loading_vector[
      positive_genes
    ],
    stringsAsFactors = FALSE
  )

  loading_rows[[length(
    loading_rows
  ) + 1L]] <- data.frame(
    PC = pc_name,
    direction = "negative",
    rank = seq_along(
      negative_genes
    ),
    gene = negative_genes,
    loading = loading_vector[
      negative_genes
    ],
    stringsAsFactors = FALSE
  )
}

pca_loading_table <- do.call(
  rbind,
  loading_rows
)

write_csv_safe(
  pca_loading_table,
  "05_top_PCA_loadings.csv"
)

# -----------------------------------------------------------------------------
# 8. PCA technical-covariate audit
# -----------------------------------------------------------------------------
technical_variables <- c(
  "nCount_RNA",
  "nFeature_RNA",
  "percent.mt",
  "percent.HB",
  "scDblFinder.score"
)

technical_variables <- technical_variables[
  technical_variables %in%
    names(obj@meta.data)
]

technical_rows <- list()

for (pc_name in colnames(
  pca_embeddings
)[PCS_USE]) {
  pc_vector <- pca_embeddings[
    ,
    pc_name
  ]

  for (variable in technical_variables) {
    technical_rows[[length(
      technical_rows
    ) + 1L]] <- data.frame(
      PC = pc_name,
      variable = variable,
      spearman_r = safe_spearman(
        pc_vector,
        as.numeric(
          obj@meta.data[
            ,
            variable
          ]
        )
      ),
      stringsAsFactors = FALSE
    )
  }
}

pc_technical_correlations <- do.call(
  rbind,
  technical_rows
)

pc_technical_correlations$absolute_r <- abs(
  pc_technical_correlations$spearman_r
)

pc_technical_correlations$review_abs_r_ge_0_50 <-
  pc_technical_correlations$absolute_r >=
  0.50

write_csv_safe(
  pc_technical_correlations,
  "06_PC_technical_correlations.csv"
)

# -----------------------------------------------------------------------------
# 9. Neighbor graph and resolution sweep
# -----------------------------------------------------------------------------
log_message(
  "Constructing PCA neighbor and SNN graphs."
)

obj <- Seurat::FindNeighbors(
  object = obj,
  reduction = "pca",
  dims = PCS_USE,
  k.param = K_PARAM,
  compute.SNN = TRUE,
  prune.SNN = SNN_PRUNE,
  nn.method = "annoy",
  n.trees = N_TREES,
  graph.name = c(
    "RNA_pca_nn",
    "RNA_pca_snn"
  ),
  verbose = TRUE
)

if (
  !all(
    c(
      "RNA_pca_nn",
      "RNA_pca_snn"
    ) %in% names(
      obj@graphs
    )
  )
) {
  stop_with_report(
    "PCA neighbor/SNN graphs were not created."
  )
}

resolution_rows <- list()

for (i in seq_along(
  CLUSTER_RESOLUTIONS
)) {
  resolution <- CLUSTER_RESOLUTIONS[[i]]
  cluster_column <- resolution_column_name(
    resolution
  )

  log_message(
    "Clustering resolution ",
    resolution,
    "."
  )

  obj <- Seurat::FindClusters(
    object = obj,
    graph.name = "RNA_pca_snn",
    cluster.name = cluster_column,
    resolution = resolution,
    algorithm = 1,
    n.start = 10L,
    n.iter = 10L,
    random.seed = SEED + i,
    group.singletons = TRUE,
    verbose = TRUE
  )

  if (!cluster_column %in%
      names(obj@meta.data)) {
    stop_with_report(
      paste0(
        "Cluster column was not created: ",
        cluster_column
      )
    )
  }

  cluster_sizes_i <- table(
    obj@meta.data[
      ,
      cluster_column
    ]
  )

  resolution_rows[[i]] <- data.frame(
    resolution = resolution,
    cluster_column = cluster_column,
    n_clusters = length(
      cluster_sizes_i
    ),
    minimum_cluster_size = min(
      cluster_sizes_i
    ),
    median_cluster_size = stats::median(
      as.numeric(
        cluster_sizes_i
      )
    ),
    maximum_cluster_size = max(
      cluster_sizes_i
    ),
    stringsAsFactors = FALSE
  )
}

resolution_summary <- do.call(
  rbind,
  resolution_rows
)

write_csv_safe(
  resolution_summary,
  "07_clustering_resolution_summary.csv"
)

primary_cluster_column <- resolution_column_name(
  PRIMARY_RESOLUTION
)

if (!primary_cluster_column %in%
    names(obj@meta.data)) {
  stop_with_report(
    paste0(
      "Primary cluster column is absent: ",
      primary_cluster_column
    )
  )
}

obj$seurat_clusters <- factor(
  obj@meta.data[
    ,
    primary_cluster_column
  ]
)

Seurat::Idents(
  object = obj
) <- primary_cluster_column

# -----------------------------------------------------------------------------
# 10. Unintegrated UMAP
# -----------------------------------------------------------------------------
log_message(
  "Running unintegrated PCA UMAP."
)

set.seed(SEED)

obj <- Seurat::RunUMAP(
  object = obj,
  reduction = "pca",
  dims = PCS_USE,
  reduction.name = "umap_pca_unintegrated",
  reduction.key = "UMAPpca_",
  umap.method = "uwot",
  n.neighbors = UMAP_N_NEIGHBORS,
  n.components = 2L,
  metric = "cosine",
  min.dist = UMAP_MIN_DIST,
  uwot.sgd = TRUE,
  seed.use = SEED,
  verbose = TRUE
)

if (!"umap_pca_unintegrated" %in%
    names(obj@reductions)) {
  stop_with_report(
    "Unintegrated PCA UMAP was not created."
  )
}

umap_embeddings <- Seurat::Embeddings(
  object = obj,
  reduction = "umap_pca_unintegrated"
)

if (
  nrow(umap_embeddings) != ncol(obj) ||
    ncol(umap_embeddings) != 2L
) {
  stop_with_report(
    "UMAP dimensions are incorrect."
  )
}

# -----------------------------------------------------------------------------
# 11. Primary-cluster composition audit
# -----------------------------------------------------------------------------
cluster_meta <- obj@meta.data
cluster_meta$barcode <- rownames(
  cluster_meta
)
cluster_meta$primary_cluster <- as.character(
  cluster_meta[
    ,
    primary_cluster_column
  ]
)

library_composition <- data.table::as.data.table(
  cluster_meta
)[
  ,
  .(
    n_cells = .N
  ),
  by = .(
    primary_cluster,
    biological_library,
    timepoint
  )
]

cluster_totals <- data.table::as.data.table(
  cluster_meta
)[
  ,
  .(
    cluster_n_cells = .N
  ),
  by = primary_cluster
]

library_composition <- merge(
  library_composition,
  cluster_totals,
  by = "primary_cluster",
  all.x = TRUE,
  sort = FALSE
)

library_composition$within_cluster_fraction <-
  library_composition$n_cells /
  library_composition$cluster_n_cells

write_csv_safe(
  as.data.frame(
    library_composition
  ),
  "08_primary_cluster_library_composition.csv"
)

timepoint_composition <- data.table::as.data.table(
  cluster_meta
)[
  ,
  .(
    n_cells = .N
  ),
  by = .(
    primary_cluster,
    timepoint
  )
]

timepoint_composition <- merge(
  timepoint_composition,
  cluster_totals,
  by = "primary_cluster",
  all.x = TRUE,
  sort = FALSE
)

timepoint_composition$within_cluster_fraction <-
  timepoint_composition$n_cells /
  timepoint_composition$cluster_n_cells

write_csv_safe(
  as.data.frame(
    timepoint_composition
  ),
  "08B_primary_cluster_timepoint_composition.csv"
)

library_dominance <- data.table::as.data.table(
  library_composition
)[
  order(
    primary_cluster,
    -within_cluster_fraction,
    biological_library
  ),
  .SD[1L],
  by = primary_cluster
]

timepoint_dominance <- data.table::as.data.table(
  timepoint_composition
)[
  order(
    primary_cluster,
    -within_cluster_fraction,
    timepoint
  ),
  .SD[1L],
  by = primary_cluster
]

cluster_summary <- data.table::as.data.table(
  cluster_meta
)[
  ,
  .(
    n_cells = .N,
    n_libraries = uniqueN(
      biological_library
    ),
    n_timepoints = uniqueN(
      timepoint
    ),
    median_nFeature_RNA = median(
      nFeature_RNA
    ),
    median_nCount_RNA = median(
      nCount_RNA
    ),
    median_percent_mt = median(
      percent.mt
    )
  ),
  by = primary_cluster
]

cluster_summary <- merge(
  cluster_summary,
  library_dominance[
    ,
    c(
      "primary_cluster",
      "biological_library",
      "within_cluster_fraction"
    ),
    with = FALSE
  ],
  by = "primary_cluster",
  all.x = TRUE,
  sort = FALSE
)

names(cluster_summary)[
  names(cluster_summary) ==
    "biological_library"
] <- "dominant_library"

names(cluster_summary)[
  names(cluster_summary) ==
    "within_cluster_fraction"
] <- "maximum_library_fraction"

cluster_summary <- merge(
  cluster_summary,
  timepoint_dominance[
    ,
    c(
      "primary_cluster",
      "timepoint",
      "within_cluster_fraction"
    ),
    with = FALSE
  ],
  by = "primary_cluster",
  all.x = TRUE,
  sort = FALSE
)

names(cluster_summary)[
  names(cluster_summary) ==
    "timepoint"
] <- "dominant_timepoint"

names(cluster_summary)[
  names(cluster_summary) ==
    "within_cluster_fraction"
] <- "maximum_timepoint_fraction"

cluster_summary$library_dominance_review <-
  cluster_summary$n_cells >=
  MIN_CLUSTER_SIZE_FOR_DOMINANCE_REVIEW &
  cluster_summary$maximum_library_fraction >=
  LIBRARY_DOMINANCE_REVIEW_THRESHOLD

cluster_summary$timepoint_dominance_review <-
  cluster_summary$n_cells >=
  MIN_CLUSTER_SIZE_FOR_DOMINANCE_REVIEW &
  cluster_summary$maximum_timepoint_fraction >=
  TIMEPOINT_DOMINANCE_REVIEW_THRESHOLD

cluster_summary <- as.data.frame(
  cluster_summary
)

cluster_summary <- cluster_summary[
  order(
    suppressWarnings(
      as.numeric(
        cluster_summary$primary_cluster
      )
    ),
    cluster_summary$primary_cluster
  ),
  ,
  drop = FALSE
]

write_csv_safe(
  cluster_summary,
  "09_primary_cluster_summary.csv"
)

# Author labels are retained only as a descriptive cross-check.
author_label_columns <- intersect(
  c(
    "layer1",
    "layer2",
    "layer3",
    "layer4",
    "layer5",
    "name_RNA_snn_res.0.2",
    "name_RNA_snn_res.0.5"
  ),
  names(cluster_meta)
)

author_cross_rows <- list()

for (author_column in author_label_columns) {
  author_values <- as.character(
    cluster_meta[
      ,
      author_column
    ]
  )

  valid <- !is.na(author_values) &
    nzchar(author_values)

  if (!any(valid)) {
    next
  }

  tab <- data.table::data.table(
    primary_cluster = cluster_meta$primary_cluster[
      valid
    ],
    author_label = author_values[
      valid
    ]
  )[
    ,
    .(
      n_cells = .N
    ),
    by = .(
      primary_cluster,
      author_label
    )
  ]

  tab$author_label_column <- author_column

  author_cross_rows[[length(
    author_cross_rows
  ) + 1L]] <- as.data.frame(
    tab
  )
}

if (length(author_cross_rows) > 0L) {
  author_cross_table <- do.call(
    rbind,
    author_cross_rows
  )
} else {
  author_cross_table <- data.frame(
    primary_cluster = character(0),
    author_label = character(0),
    n_cells = integer(0),
    author_label_column = character(0),
    stringsAsFactors = FALSE
  )
}

write_csv_safe(
  author_cross_table,
  "10_primary_cluster_vs_author_labels.csv"
)

# -----------------------------------------------------------------------------
# 12. Figures
# -----------------------------------------------------------------------------
log_message(
  "Creating PCA/UMAP technical-audit figures."
)

variable_feature_plot <- Seurat::VariableFeaturePlot(
  object = obj,
  assay = "RNA"
)

top_variable_features <- head(
  variable_features,
  20L
)

variable_feature_plot_labeled <- Seurat::LabelPoints(
  plot = variable_feature_plot,
  points = top_variable_features,
  repel = TRUE
)

ggplot2::ggsave(
  filename = file.path(
    figures_dir,
    "01_variable_features.pdf"
  ),
  plot = variable_feature_plot_labeled,
  width = 9,
  height = 7,
  limitsize = FALSE
)

p_elbow <- Seurat::ElbowPlot(
  object = obj,
  reduction = "pca",
  ndims = N_PCS_COMPUTE
) +
  ggplot2::geom_vline(
    xintercept = max(PCS_USE),
    linetype = 2
  ) +
  ggplot2::labs(
    title = "PCA elbow",
    subtitle = paste0(
      "PC1-PC",
      max(PCS_USE),
      " used for neighbors and UMAP"
    )
  )

p_variance <- ggplot2::ggplot(
  pca_table,
  ggplot2::aes(
    x = PC,
    y = cumulative_percent_variance
  )
) +
  ggplot2::geom_line() +
  ggplot2::geom_point(
    size = 1
  ) +
  ggplot2::geom_vline(
    xintercept = max(PCS_USE),
    linetype = 2
  ) +
  ggplot2::labs(
    x = "PC",
    y = "Cumulative variance (%)",
    title = "PCA cumulative variance",
    subtitle = "Variance uses stdev squared"
  ) +
  ggplot2::theme_bw()

ggplot2::ggsave(
  filename = file.path(
    figures_dir,
    "02_PCA_elbow_and_variance.pdf"
  ),
  plot = p_elbow + p_variance,
  width = 13,
  height = 6,
  limitsize = FALSE
)

p_umap_cluster <- Seurat::DimPlot(
  object = obj,
  reduction = "umap_pca_unintegrated",
  group.by = primary_cluster_column,
  label = TRUE,
  repel = TRUE,
  raster = TRUE
) +
  ggplot2::labs(
    title = paste0(
      "Unintegrated PCA UMAP: resolution ",
      PRIMARY_RESOLUTION
    )
  ) +
  Seurat::NoLegend()

p_umap_timepoint <- Seurat::DimPlot(
  object = obj,
  reduction = "umap_pca_unintegrated",
  group.by = "timepoint",
  raster = TRUE
) +
  ggplot2::labs(
    title = "Unintegrated PCA UMAP by timepoint"
  )

p_umap_library <- Seurat::DimPlot(
  object = obj,
  reduction = "umap_pca_unintegrated",
  group.by = "biological_library",
  raster = TRUE
) +
  ggplot2::labs(
    title = "Unintegrated PCA UMAP by biological library"
  )

ggplot2::ggsave(
  filename = file.path(
    figures_dir,
    "03_UMAP_cluster_timepoint_library.pdf"
  ),
  plot = (
    p_umap_cluster |
      p_umap_timepoint
  ) / p_umap_library,
  width = 16,
  height = 13,
  limitsize = FALSE
)

p_cluster_size <- ggplot2::ggplot(
  cluster_summary,
  ggplot2::aes(
    x = reorder(
      primary_cluster,
      n_cells
    ),
    y = n_cells
  )
) +
  ggplot2::geom_col() +
  ggplot2::coord_flip() +
  ggplot2::scale_y_continuous(
    labels = scales::comma
  ) +
  ggplot2::labs(
    x = paste0(
      "Primary cluster (resolution ",
      PRIMARY_RESOLUTION,
      ")"
    ),
    y = "Nuclei",
    title = "Primary cluster sizes"
  ) +
  ggplot2::theme_bw()

ggplot2::ggsave(
  filename = file.path(
    figures_dir,
    "04_primary_cluster_sizes.pdf"
  ),
  plot = p_cluster_size,
  width = 9,
  height = 7,
  limitsize = FALSE
)

# -----------------------------------------------------------------------------
# 13. Remove scale.data after PCA and save
# -----------------------------------------------------------------------------
# PCA, graphs, clusters, and UMAP are already stored.
# Removing scale.data avoids saving a large dense 3000 x 144359 matrix.
obj[["RNA"]]$scale.data <- NULL

obj@misc$GSE234774_step10B <- list(
  input_object = normalizePath(
    input_object_path,
    winslash = "/",
    mustWork = TRUE
  ),
  normalization = "LogNormalize",
  normalization_scale_factor =
    NORMALIZATION_SCALE_FACTOR,
  variable_features =
    N_VARIABLE_FEATURES,
  pca_features =
    length(pca_features),
  pca_dimensions_computed =
    N_PCS_COMPUTE,
  pca_dimensions_used =
    PCS_USE,
  primary_reduction =
    "pca",
  umap_reduction =
    "umap_pca_unintegrated",
  cluster_resolutions =
    CLUSTER_RESOLUTIONS,
  primary_resolution =
    PRIMARY_RESOLUTION,
  primary_cluster_column =
    primary_cluster_column,
  harmony_run = FALSE,
  sctransform_run = FALSE,
  cell_cycle_regression_run = FALSE,
  technical_covariate_regression_run = FALSE,
  biological_unit =
    "biological_library",
  age_status =
    "METHODS_REVIEW_REQUIRED",
  completed_at =
    as.character(
      Sys.time()
    ),
  interpretation = paste(
    "Whole-atlas technical clustering only.",
    "No cell-type annotation and no Arm_1/Arm_2 scoring."
  )
)

log_message(
  "Saving the Step 10B object."
)

saveRDS(
  obj,
  output_object_path,
  compress = FALSE
)

log_message(
  "Saved object: ",
  output_object_path
)

# -----------------------------------------------------------------------------
# 14. Reload and final audit
# -----------------------------------------------------------------------------
reload_obj <- tryCatch(
  readRDS(
    output_object_path
  ),
  error = function(e) e
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
  nrow(reload_obj) == 26918L &&
  ncol(reload_obj) == 144359L &&
  "pca" %in%
    names(
      reload_obj@reductions
    ) &&
  "umap_pca_unintegrated" %in%
    names(
      reload_obj@reductions
    ) &&
  primary_cluster_column %in%
    names(
      reload_obj@meta.data
    )

rm(
  reload_obj
)

invisible(gc())

n_primary_clusters <- length(
  unique(
    as.character(
      obj@meta.data[
        ,
        primary_cluster_column
      ]
    )
  )
)

maximum_abs_technical_r <- max(
  pc_technical_correlations$absolute_r,
  na.rm = TRUE
)

n_library_dominant_clusters <- sum(
  cluster_summary$library_dominance_review,
  na.rm = TRUE
)

n_timepoint_dominant_clusters <- sum(
  cluster_summary$timepoint_dominance_review,
  na.rm = TRUE
)

audit_checks <- data.frame(
  check = c(
    "Input contains 26,918 features and 144,359 singlet nuclei",
    "Input contains 20 biological libraries and seven expected timepoints",
    "LogNormalize data are stored in a readable BPCells directory",
    "Exactly 3,000 variable features were selected",
    "At least 2,000 non-mitochondrial/non-hemoglobin PCA features remain",
    "Fifty PCs were computed",
    "PC1-PC40 are available for neighbors and UMAP",
    "PCA neighbor and SNN graphs were created",
    "All requested clustering resolutions were created",
    "Primary clustering contains at least two clusters",
    "UMAP contains all singlet nuclei",
    "Saved object reloads with PCA, UMAP, and primary clusters",
    "Harmony/integration was not run",
    "SCTransform was not run",
    "Cell-cycle and technical covariates were not regressed",
    "Cluster library/timepoint dominance has been manually reviewed",
    "Cohort age has been confirmed from the study methods"
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
    TRUE,
    TRUE,
    FALSE,
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
      n_libraries == 20L &&
        all(
          expected_timepoints %in%
            observed_timepoints
        ),
      "PASS",
      "STOP"
    ),
    ifelse(
      nrow(normalized_disk) ==
        nrow(obj) &&
        ncol(normalized_disk) ==
        ncol(obj),
      "PASS",
      "STOP"
    ),
    ifelse(
      length(variable_features) ==
        N_VARIABLE_FEATURES,
      "PASS",
      "STOP"
    ),
    ifelse(
      length(pca_features) >=
        2000L,
      "PASS",
      "STOP"
    ),
    ifelse(
      ncol(pca_embeddings) >=
        N_PCS_COMPUTE,
      "PASS",
      "STOP"
    ),
    ifelse(
      max(PCS_USE) <=
        ncol(pca_embeddings),
      "PASS",
      "STOP"
    ),
    ifelse(
      all(
        c(
          "RNA_pca_nn",
          "RNA_pca_snn"
        ) %in% names(
          obj@graphs
        )
      ),
      "PASS",
      "STOP"
    ),
    ifelse(
      all(
        vapply(
          CLUSTER_RESOLUTIONS,
          resolution_column_name,
          character(1)
        ) %in%
          names(
            obj@meta.data
          )
      ),
      "PASS",
      "STOP"
    ),
    ifelse(
      n_primary_clusters >= 2L,
      "PASS",
      "STOP"
    ),
    ifelse(
      nrow(umap_embeddings) ==
        ncol(obj),
      "PASS",
      "STOP"
    ),
    ifelse(
      reload_pass,
      "PASS",
      "STOP"
    ),
    "CONFIRMED_NOT_RUN",
    "CONFIRMED_NOT_RUN",
    "CONFIRMED_NOT_RUN",
    "REVIEW",
    "REVIEW"
  ),
  detail = c(
    paste0(
      nrow(obj),
      " features x ",
      ncol(obj),
      " singlets"
    ),
    paste0(
      n_libraries,
      " libraries; timepoints=",
      paste(
        expected_timepoints,
        collapse = ", "
      )
    ),
    normalizePath(
      normalized_data_dir,
      winslash = "/",
      mustWork = TRUE
    ),
    paste0(
      "selected=",
      length(variable_features)
    ),
    paste0(
      "PCA features=",
      length(pca_features),
      "; excluded mt=",
      length(mt_hvg),
      "; excluded HB=",
      length(hb_hvg)
    ),
    paste0(
      "PCs=",
      ncol(pca_embeddings)
    ),
    paste0(
      "used dimensions=",
      min(PCS_USE),
      "-",
      max(PCS_USE)
    ),
    paste(
      names(
        obj@graphs
      ),
      collapse = ", "
    ),
    paste(
      vapply(
        CLUSTER_RESOLUTIONS,
        resolution_column_name,
        character(1)
      ),
      collapse = ", "
    ),
    paste0(
      "primary clusters=",
      n_primary_clusters,
      " at resolution ",
      PRIMARY_RESOLUTION
    ),
    paste0(
      "UMAP rows=",
      nrow(umap_embeddings)
    ),
    ifelse(
      reload_pass,
      "RDS reload PASS",
      "RDS reload failed"
    ),
    "No harmony or integrated reduction exists in this step.",
    "Primary normalization is LogNormalize.",
    "No vars.to.regress were supplied.",
    paste0(
      "maximum absolute PC technical correlation=",
      round(
        maximum_abs_technical_r,
        3
      ),
      "; library-dominant cluster flags=",
      n_library_dominant_clusters,
      "; timepoint-dominant cluster flags=",
      n_timepoint_dominant_clusters
    ),
    "Age remains METHODS_REVIEW_REQUIRED."
  ),
  stringsAsFactors = FALSE
)

write_csv_safe(
  audit_checks,
  "11_STEP10B_audit_checks.csv"
)

required_fail <-
  audit_checks$required_for_pass &
  audit_checks$status != "PASS"

final_decision <- if (
  any(
    required_fail
  )
) {
  "STOP"
} else {
  "PASS_WITH_CLUSTER_REVIEW"
}

# -----------------------------------------------------------------------------
# 15. Final report
# -----------------------------------------------------------------------------
resolution_lines <- paste0(
  "resolution ",
  resolution_summary$resolution,
  ": ",
  resolution_summary$n_clusters,
  " clusters; min/median/max size=",
  resolution_summary$minimum_cluster_size,
  "/",
  round(
    resolution_summary$median_cluster_size,
    1
  ),
  "/",
  resolution_summary$maximum_cluster_size
)

report_lines <- c(
  "GSE234774 STEP 10B NORMALIZATION / PCA / UMAP / CLUSTERING",
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
  "INPUT",
  paste0(
    "Features: ",
    nrow(obj)
  ),
  paste0(
    "QC singlets: ",
    ncol(obj)
  ),
  paste0(
    "Biological libraries: ",
    n_libraries
  ),
  paste0(
    "Timepoints: ",
    paste(
      expected_timepoints,
      collapse = ", "
    )
  ),
  "",
  "NORMALIZATION",
  "Primary method: LogNormalize",
  paste0(
    "Scale factor: ",
    NORMALIZATION_SCALE_FACTOR
  ),
  paste0(
    "Variable features: ",
    length(variable_features)
  ),
  paste0(
    "PCA features after mt/HB exclusion: ",
    length(pca_features)
  ),
  "SCTransform: NOT RUN",
  "Cell-cycle regression: NOT RUN",
  "Technical-covariate regression: NOT RUN",
  "",
  "PCA / GRAPH / UMAP",
  paste0(
    "PCs computed: ",
    ncol(pca_embeddings)
  ),
  paste0(
    "PCs used: ",
    min(PCS_USE),
    "-",
    max(PCS_USE)
  ),
  paste0(
    "Maximum absolute PC-versus-technical Spearman correlation: ",
    round(
      maximum_abs_technical_r,
      3
    )
  ),
  paste0(
    "Neighbor k: ",
    K_PARAM
  ),
  "Primary dimensional reduction: unintegrated PCA",
  "UMAP: umap_pca_unintegrated",
  "Harmony/integration: NOT RUN",
  "",
  "CLUSTERING RESOLUTION SWEEP",
  resolution_lines,
  "",
  paste0(
    "Primary resolution: ",
    PRIMARY_RESOLUTION
  ),
  paste0(
    "Primary cluster column: ",
    primary_cluster_column
  ),
  paste0(
    "Primary clusters: ",
    n_primary_clusters
  ),
  paste0(
    "Clusters flagged for >=80% single-library dominance: ",
    n_library_dominant_clusters
  ),
  paste0(
    "Clusters flagged for >=95% single-timepoint dominance: ",
    n_timepoint_dominant_clusters
  ),
  "",
  "INTERPRETATION",
  "1. This is a whole-atlas technical clustering, not a final biological annotation.",
  "2. biological_library is nested within timepoint; automatic global Harmony could remove true longitudinal biology.",
  "3. A library- or timepoint-dominant cluster is a review flag, not automatic evidence of artifact.",
  "4. Author labels are retained only as supporting cross-checks.",
  "5. No cluster is called resident microglia, macrophage, Arm_1, or Arm_2 in this step.",
  "6. GSE234774 validates frozen GSE162610 programs and does not redefine them.",
  "",
  "PRIMARY OBJECT",
  output_object_path,
  "",
  "NEXT STEP",
  "Review PCA technical correlations, UMAP structure, and cluster composition.",
  "Then perform marker-first resident microglia versus monocyte/macrophage identification.",
  "Do not start frozen Arm_1/Arm_2 signature transfer before the clean resident-microglia object is frozen."
)

if (final_decision == "STOP") {
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

  for (i in seq_len(
    nrow(failed)
  )) {
    report_lines <- c(
      report_lines,
      paste0(
        i,
        ". ",
        failed$check[[i]],
        " -- ",
        failed$detail[[i]]
      )
    )
  }
}

report_path <- file.path(
  output_dir,
  "12_STEP10B_FINAL_REPORT.txt"
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
cat("GSE234774 STEP 10B COMPLETED\n")
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
cat("Primary object:\n")
cat(
  "  ",
  output_object_path,
  "\n",
  sep = ""
)
cat("Please upload these five files for review:\n")
cat("  1) tables/04_PCA_variance.csv\n")
cat("  2) tables/06_PC_technical_correlations.csv\n")
cat("  3) tables/09_primary_cluster_summary.csv\n")
cat("  4) tables/11_STEP10B_audit_checks.csv\n")
cat("  5) 12_STEP10B_FINAL_REPORT.txt\n")
cat("============================================================\n")

if (final_decision == "STOP") {
  stop(
    "GSE234774 Step 10B decision is STOP. Review the final report.",
    call. = FALSE
  )
}
