# ============================================================
# 13_step3A_GSE162610_microglia_reclustering.R
#
# Project:
#   Adult SCI activated/DAM-like microglia fate bifurcation
#
# Step 3A purpose:
#   De novo reclustering of the FINAL clean resident-derived microglia
#
# Formal input:
#   <PROJECT_ROOT>/04_results/
#   step2C_GSE162610_final_microglia_extraction/
#   GSE162610_microglia_clean.rds
#
# Core workflow:
#   RNA raw counts
#   -> LogNormalize
#   -> 3000 microglia-only HVGs
#   -> ScaleData
#   -> PCA 50 PCs
#   -> Harmony(sample_id ONLY)
#   -> Harmony-based neighbors
#   -> multi-resolution clustering
#   -> Harmony-based UMAP
#   -> exploratory state marker audit
#
# IMPORTANT:
#   - Do NOT regress cell cycle.
#   - Do NOT regress timepoint.
#   - Do NOT use SCT as the primary workflow.
#   - Do NOT define DAM/repair/scar yet.
#   - This step is STATE DISCOVERY, not final state annotation.
#
# Directory convention:
#   Root scripts: <PROJECT_ROOT>/NN_stepX_*.R
#   Raw data:     <PROJECT_ROOT>/01_raw_data
#   Results/RDS:  <PROJECT_ROOT>/04_results/stepX_*
# ============================================================

options(stringsAsFactors = FALSE)
set.seed(1234)

# ------------------------------------------------------------
# 0. Packages
# ------------------------------------------------------------
required_pkgs <- c(
  "Seurat",
  "SeuratObject",
  "dplyr",
  "ggplot2",
  "patchwork",
  "harmony",
  "clustree",
  "pheatmap",
  "Matrix"
)

missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_pkgs) > 0) {
  stop(
    "缺少以下 R 包，请先安装后再 Source：\n",
    paste(missing_pkgs, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(harmony)
  library(clustree)
  library(pheatmap)
  library(Matrix)
})

# ------------------------------------------------------------
# 1. Unified project paths
# ------------------------------------------------------------
project_dir <- Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
input_rds <- file.path(
  project_dir,
  "04_results",
  "step2C_GSE162610_final_microglia_extraction",
  "GSE162610_microglia_clean.rds"
)

out_dir <- file.path(
  project_dir,
  "04_results",
  "step3A_GSE162610_microglia_reclustering"
)

fig_dir <- file.path(out_dir, "figures")
tab_dir <- file.path(out_dir, "tables")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input_rds)) {
  stop(
    "找不到 Step 3A 正式输入对象：\n",
    input_rds
  )
}

# ------------------------------------------------------------
# 2. Read formal clean microglia object
# ------------------------------------------------------------
message("====================================================")
message("GSE162610 Step 3A | clean microglia de novo reclustering")
message("====================================================")
message("Reading:")
message(input_rds)

obj <- readRDS(input_rds)

if (!inherits(obj, "Seurat")) {
  stop("输入 RDS 不是 Seurat 对象。")
}

if (!"RNA" %in% names(obj@assays)) {
  stop("输入对象中没有 RNA assay。")
}

DefaultAssay(obj) <- "RNA"

message("Input cells: ", ncol(obj))
message("Input features: ", nrow(obj))

# expected formal Step 2C cell number
expected_n <- 20747

if (ncol(obj) != expected_n) {
  warning(
    "Step 3A 输入细胞数与 Step 2C 预期不一致：",
    " expected = ", expected_n,
    ", observed = ", ncol(obj)
  )
}

# ------------------------------------------------------------
# 3. Metadata audit
# ------------------------------------------------------------
meta <- obj[[]]

required_meta <- c(
  "sample_id",
  "timepoint",
  "step2_final_cluster"
)

missing_meta <- required_meta[
  !required_meta %in% colnames(meta)
]

if (length(missing_meta) > 0) {
  stop(
    "缺少必要 metadata：",
    paste(missing_meta, collapse = ", ")
  )
}

# keep old Step 2 provenance as character
obj$step2_final_cluster <- as.character(obj$step2_final_cluster)

message("Sample counts:")
print(table(obj$sample_id))

message("Timepoint counts:")
print(table(obj$timepoint))

message("Step 2 provenance counts:")
print(table(obj$step2_final_cluster))

# ------------------------------------------------------------
# 4. Seurat v5 layer handling
# ------------------------------------------------------------
# JoinLayers if available; failure should not stop if the object is already joined.
joinlayers_available <- exists(
  "JoinLayers",
  where = asNamespace("SeuratObject"),
  mode = "function",
  inherits = FALSE
)

if (isTRUE(joinlayers_available)) {
  obj <- tryCatch(
    {
      message("Trying JoinLayers(RNA)...")
      SeuratObject::JoinLayers(
        obj,
        assay = "RNA"
      )
    },
    error = function(e) {
      message(
        "JoinLayers skipped safely: ",
        e$message
      )
      obj
    }
  )
}

DefaultAssay(obj) <- "RNA"

# ------------------------------------------------------------
# 5. Preserve Step 2 provenance and remove dependence on old reductions
# ------------------------------------------------------------
# We deliberately create NEW reduction names in Step 3A.
# Old reductions remain only as provenance and are never used for Step 3A.

old_reductions <- Reductions(obj)

writeLines(
  old_reductions,
  file.path(
    tab_dir,
    "00_input_old_reductions.txt"
  )
)

# ------------------------------------------------------------
# 6. De novo LogNormalize
# ------------------------------------------------------------
message("Step 3A normalization: LogNormalize")

obj <- NormalizeData(
  obj,
  assay = "RNA",
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = TRUE
)

# ------------------------------------------------------------
# 7. Microglia-only variable features
# ------------------------------------------------------------
message("Finding 3000 microglia-only HVGs")

obj <- FindVariableFeatures(
  obj,
  assay = "RNA",
  selection.method = "vst",
  nfeatures = 3000,
  verbose = TRUE
)

hvg <- VariableFeatures(obj)

writeLines(
  hvg,
  file.path(
    tab_dir,
    "01_microglia_HVG_3000.txt"
  )
)

# Variable-feature figure
p_vf <- VariableFeaturePlot(obj)

top10_hvg <- head(hvg, 10)

p_vf_label <- LabelPoints(
  plot = p_vf,
  points = top10_hvg,
  repel = TRUE
)

ggsave(
  file.path(
    fig_dir,
    "01_microglia_variable_features.pdf"
  ),
  p_vf_label,
  width = 8,
  height = 6
)

# ------------------------------------------------------------
# 8. ScaleData
# ------------------------------------------------------------
# IMPORTANT:
# No cell-cycle regression.
# No sample regression.
# No timepoint regression.
message("ScaleData WITHOUT cell-cycle/timepoint/sample regression")

obj <- ScaleData(
  obj,
  assay = "RNA",
  features = hvg,
  verbose = TRUE
)

# ------------------------------------------------------------
# 9. PCA
# ------------------------------------------------------------
message("Running PCA: 50 PCs")

obj <- RunPCA(
  obj,
  assay = "RNA",
  features = hvg,
  npcs = 50,
  reduction.name = "pca.step3A",
  reduction.key = "PC3A_",
  verbose = TRUE
)

pca_stdev <- Stdev(
  obj,
  reduction = "pca.step3A"
)

pca_var <- pca_stdev^2 / sum(pca_stdev^2)
pca_cum <- cumsum(pca_var)

pca_tbl <- data.frame(
  PC = seq_along(pca_var),
  variance_fraction = pca_var,
  cumulative_variance = pca_cum,
  stringsAsFactors = FALSE
)

write.csv(
  pca_tbl,
  file.path(
    tab_dir,
    "02_microglia_PCA_variance.csv"
  ),
  row.names = FALSE
)

# Candidate PC summary
pc80 <- which(pca_cum >= 0.80)[1]
pc90 <- which(pca_cum >= 0.90)[1]

pc_summary <- data.frame(
  criterion = c(
    "80_percent_cumulative_variance",
    "90_percent_cumulative_variance",
    "formal_Step3A_dims"
  ),
  PC = c(
    pc80,
    pc90,
    30
  )
)

write.csv(
  pc_summary,
  file.path(
    tab_dir,
    "03_candidate_PC_summary.csv"
  ),
  row.names = FALSE
)

# Formal Step 3A choice:
# 30 PCs to retain subtle within-microglia state biology.
dims_use <- 1:30

# PCA elbow
pdf(
  file.path(
    fig_dir,
    "02_microglia_PCA_ElbowPlot_50PCs.pdf"
  ),
  width = 7,
  height = 6
)

print(
  ElbowPlot(
    obj,
    ndims = 50,
    reduction = "pca.step3A"
  )
)

dev.off()

# PCA by sample/timepoint
p_pca_sample <- DimPlot(
  obj,
  reduction = "pca.step3A",
  group.by = "sample_id"
) +
  ggtitle("Microglia PCA | biological sample")

p_pca_time <- DimPlot(
  obj,
  reduction = "pca.step3A",
  group.by = "timepoint"
) +
  ggtitle("Microglia PCA | SCI timepoint")

ggsave(
  file.path(
    fig_dir,
    "03_microglia_PCA_sample_timepoint.pdf"
  ),
  p_pca_sample | p_pca_time,
  width = 14,
  height = 6
)

# ------------------------------------------------------------
# 10. Pre-Harmony UMAP audit
# ------------------------------------------------------------
# This is only for before/after comparison.
# Formal clustering will NOT use this reduction.
obj <- RunUMAP(
  obj,
  reduction = "pca.step3A",
  dims = dims_use,
  reduction.name = "umap.preHarmony.step3A",
  reduction.key = "preH3A_",
  seed.use = 1234,
  verbose = TRUE
)

p_pre_sample <- DimPlot(
  obj,
  reduction = "umap.preHarmony.step3A",
  group.by = "sample_id"
) +
  ggtitle("Pre-Harmony | sample")

p_pre_time <- DimPlot(
  obj,
  reduction = "umap.preHarmony.step3A",
  group.by = "timepoint"
) +
  ggtitle("Pre-Harmony | timepoint")

ggsave(
  file.path(
    fig_dir,
    "04_preHarmony_microglia_UMAP.pdf"
  ),
  p_pre_sample | p_pre_time,
  width = 14,
  height = 6
)

# ------------------------------------------------------------
# 11. Harmony
# ------------------------------------------------------------
# Correct approach for this project:
#   group.by.vars = sample_id
#   DO NOT regress/correct timepoint
#   Formal neighbors/UMAP/clustering MUST use Harmony reduction.
message("Running Harmony by sample_id ONLY")

obj <- tryCatch(
  {
    RunHarmony(
      object = obj,
      group.by.vars = "sample_id",
      reduction.use = "pca.step3A",
      dims.use = dims_use,
      theta = 1,
      max.iter.harmony = 20,
      reduction.save = "harmony.step3A",
      verbose = TRUE
    )
  },
  error = function(e1) {
    message(
      "Harmony signature 1 failed: ",
      e1$message
    )
    message("Trying alternate Harmony syntax...")

    RunHarmony(
      object = obj,
      group.by.vars = "sample_id",
      reduction = "pca.step3A",
      dims.use = dims_use,
      theta = 1,
      max.iter.harmony = 20,
      reduction.save = "harmony.step3A",
      verbose = TRUE
    )
  }
)

if (!"harmony.step3A" %in% Reductions(obj)) {
  # Some harmony versions ignore reduction.save and create "harmony".
  if ("harmony" %in% Reductions(obj)) {
    obj[["harmony.step3A"]] <- obj[["harmony"]]
  } else {
    stop("Harmony reduction 未成功生成。")
  }
}

# ------------------------------------------------------------
# 12. Harmony-based UMAP
# ------------------------------------------------------------
obj <- RunUMAP(
  obj,
  reduction = "harmony.step3A",
  dims = dims_use,
  reduction.name = "umap.microglia.step3A",
  reduction.key = "mg3A_",
  seed.use = 1234,
  verbose = TRUE
)

p_post_sample <- DimPlot(
  obj,
  reduction = "umap.microglia.step3A",
  group.by = "sample_id"
) +
  ggtitle("Post-Harmony | biological sample")

p_post_time <- DimPlot(
  obj,
  reduction = "umap.microglia.step3A",
  group.by = "timepoint"
) +
  ggtitle("Post-Harmony | SCI timepoint")

p_post_origin <- DimPlot(
  obj,
  reduction = "umap.microglia.step3A",
  group.by = "step2_final_cluster",
  label = TRUE,
  repel = TRUE
) +
  ggtitle("Step 2 lineage provenance")

ggsave(
  file.path(
    fig_dir,
    "05_postHarmony_microglia_UMAP_overview.pdf"
  ),
  (p_post_sample | p_post_time) / p_post_origin,
  width = 15,
  height = 11
)

# Before vs after comparison
ggsave(
  file.path(
    fig_dir,
    "06_pre_vs_post_Harmony_microglia.pdf"
  ),
  (p_pre_sample | p_post_sample) /
    (p_pre_time | p_post_time),
  width = 15,
  height = 11
)

# ------------------------------------------------------------
# 13. Harmony-based neighbors
# ------------------------------------------------------------
# IMPORTANT:
# The old reference script used PCA for neighbors/UMAP after Harmony.
# Here we explicitly use harmony.step3A.

obj <- FindNeighbors(
  obj,
  reduction = "harmony.step3A",
  dims = dims_use,
  graph.name = c(
    "microglia_step3A_nn",
    "microglia_step3A_snn"
  ),
  verbose = TRUE
)

# ------------------------------------------------------------
# 14. Multi-resolution clustering
# ------------------------------------------------------------
resolutions <- seq(
  0.1,
  1.0,
  by = 0.1
)

cluster_number_list <- list()

for (res in resolutions) {

  cl_name <- paste0(
    "step3A_res.",
    format(
      res,
      nsmall = 1,
      trim = TRUE,
      scientific = FALSE
    )
  )

  obj <- FindClusters(
    obj,
    graph.name = "microglia_step3A_snn",
    resolution = res,
    algorithm = 1,
    random.seed = 1234,
    cluster.name = cl_name,
    verbose = FALSE
  )

  cluster_number_list[[cl_name]] <- data.frame(
    resolution = res,
    cluster_column = cl_name,
    n_clusters = length(
      unique(
        as.character(
          obj[[cl_name]][, 1]
        )
      )
    )
  )
}

cluster_number_tbl <- bind_rows(
  cluster_number_list
)

write.csv(
  cluster_number_tbl,
  file.path(
    tab_dir,
    "04_cluster_number_by_resolution.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 15. Correct Step 3A-only clustree
# ------------------------------------------------------------
tree_meta <- obj[[]][
  ,
  paste0(
    "step3A_res.",
    format(
      resolutions,
      nsmall = 1,
      trim = TRUE,
      scientific = FALSE
    )
  ),
  drop = FALSE
]

p_tree <- clustree(
  tree_meta,
  prefix = "step3A_res."
)

ggsave(
  file.path(
    fig_dir,
    "07_microglia_clustree_resolution_0.1_to_1.0.pdf"
  ),
  p_tree,
  width = 12,
  height = 9
)

# ------------------------------------------------------------
# 16. Provisional anchor resolution = 0.4
# ------------------------------------------------------------
# This is ONLY an anchor for marker/state audit.
# It is not considered final until results are reviewed.

anchor_col <- "step3A_res.0.4"

if (!anchor_col %in% colnames(obj[[]])) {
  stop(
    "未找到 provisional resolution 0.4 metadata。"
  )
}

obj$step3A_cluster_r04 <- as.character(
  obj[[anchor_col]][, 1]
)

Idents(obj) <- "step3A_cluster_r04"

message("Provisional Step 3A resolution 0.4 cluster counts:")
print(table(obj$step3A_cluster_r04))

# ------------------------------------------------------------
# 17. Provisional UMAP at r = 0.4
# ------------------------------------------------------------
p_cluster <- DimPlot(
  obj,
  reduction = "umap.microglia.step3A",
  group.by = "step3A_cluster_r04",
  label = TRUE,
  repel = TRUE
) +
  ggtitle(
    "Clean microglia de novo reclustering | provisional resolution 0.4"
  )

p_cluster_time <- DimPlot(
  obj,
  reduction = "umap.microglia.step3A",
  group.by = "timepoint"
) +
  ggtitle("SCI timepoint")

p_cluster_sample <- DimPlot(
  obj,
  reduction = "umap.microglia.step3A",
  group.by = "sample_id"
) +
  ggtitle("Biological sample")

p_cluster_origin <- DimPlot(
  obj,
  reduction = "umap.microglia.step3A",
  group.by = "step2_final_cluster",
  label = TRUE,
  repel = TRUE
) +
  ggtitle("Step 2 cluster provenance")

ggsave(
  file.path(
    fig_dir,
    "08_provisional_r04_microglia_UMAP.pdf"
  ),
  (p_cluster | p_cluster_time) /
    (p_cluster_sample | p_cluster_origin),
  width = 15,
  height = 11
)

# ------------------------------------------------------------
# 18. Exploratory state marker panel
# ------------------------------------------------------------
# These are exploratory programs only.
# They DO NOT create state labels in this step.

marker_programs <- list(

  Homeostatic = c(
    "P2ry12",
    "Tmem119",
    "Sall1",
    "Fcrls",
    "Hexb",
    "Gpr34",
    "P2ry13",
    "Siglech"
  ),

  Activated_injury = c(
    "Apoe",
    "Lpl",
    "Cst7",
    "Cd68",
    "Ctsb",
    "Ctsd",
    "Cd9",
    "Tyrobp"
  ),

  DAM_like_candidate = c(
    "Trem2",
    "Apoe",
    "Lpl",
    "Cst7",
    "Axl",
    "Clec7a",
    "Itgax",
    "Gpnmb"
  ),

  Phagocytic_lipid = c(
    "Lpl",
    "Apoe",
    "Cd68",
    "Ctsd",
    "Lgals3",
    "Fabp5",
    "Cst7"
  ),

  IFN_response = c(
    "Ifit1",
    "Ifit2",
    "Ifit3",
    "Isg15",
    "Mx1",
    "Rsad2",
    "Usp18",
    "Oasl1"
  ),

  Cycling_S = c(
    "Mcm2",
    "Mcm5",
    "Mcm6",
    "Pcna",
    "Ccne2",
    "Uhrf1",
    "Cdc6"
  ),

  Cycling_G2M = c(
    "Mki67",
    "Top2a",
    "Ube2c",
    "Birc5",
    "Cenpf",
    "Ccna2",
    "Aurkb"
  )
)

rna_genes <- rownames(obj[["RNA"]])

marker_programs_use <- lapply(
  marker_programs,
  function(gs) intersect(gs, rna_genes)
)

marker_availability <- do.call(
  rbind,
  lapply(
    names(marker_programs),
    function(nm) {
      data.frame(
        program = nm,
        n_requested = length(marker_programs[[nm]]),
        n_available = length(marker_programs_use[[nm]]),
        available_genes = paste(
          marker_programs_use[[nm]],
          collapse = ";"
        ),
        stringsAsFactors = FALSE
      )
    }
  )
)

write.csv(
  marker_availability,
  file.path(
    tab_dir,
    "05_state_marker_availability.csv"
  ),
  row.names = FALSE
)

dot_features <- unique(
  unlist(
    marker_programs_use,
    use.names = FALSE
  )
)

p_dot <- DotPlot(
  obj,
  features = dot_features,
  group.by = "step3A_cluster_r04",
  assay = "RNA"
) +
  RotatedAxis() +
  ggtitle(
    "Exploratory microglial state markers | provisional r=0.4"
  )

ggsave(
  file.path(
    fig_dir,
    "09_exploratory_state_marker_DotPlot_r04.pdf"
  ),
  p_dot,
  width = 20,
  height = 8
)

# ------------------------------------------------------------
# 19. FindAllMarkers at provisional r = 0.4
# ------------------------------------------------------------
DefaultAssay(obj) <- "RNA"
Idents(obj) <- "step3A_cluster_r04"

message("Running FindAllMarkers at provisional resolution 0.4")

markers <- FindAllMarkers(
  obj,
  assay = "RNA",
  only.pos = TRUE,
  min.pct = 0.10,
  logfc.threshold = 0.25,
  test.use = "wilcox",
  verbose = TRUE
)

write.csv(
  markers,
  file.path(
    tab_dir,
    "06_all_cluster_markers_provisional_r04.csv"
  ),
  row.names = FALSE
)

fc_col <- intersect(
  c(
    "avg_log2FC",
    "avg_logFC"
  ),
  colnames(markers)
)[1]

if (!is.na(fc_col)) {

  top20 <- markers %>%
    group_by(cluster) %>%
    arrange(
      desc(.data[[fc_col]]),
      .by_group = TRUE
    ) %>%
    slice_head(n = 20) %>%
    ungroup()

} else {

  top20 <- markers %>%
    group_by(cluster) %>%
    slice_head(n = 20) %>%
    ungroup()
}

write.csv(
  top20,
  file.path(
    tab_dir,
    "07_top20_markers_per_cluster_provisional_r04.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 20. New cluster x timepoint
# ------------------------------------------------------------
cluster_time_counts <- as.data.frame(
  table(
    cluster = obj$step3A_cluster_r04,
    timepoint = obj$timepoint
  )
)

cluster_time_counts <- cluster_time_counts %>%
  group_by(cluster) %>%
  mutate(
    cluster_total = sum(Freq),
    proportion_within_cluster =
      ifelse(
        cluster_total > 0,
        Freq / cluster_total,
        NA_real_
      )
  ) %>%
  ungroup()

write.csv(
  cluster_time_counts,
  file.path(
    tab_dir,
    "08_cluster_by_timepoint_provisional_r04.csv"
  ),
  row.names = FALSE
)

time_mat <- xtabs(
  Freq ~ cluster + timepoint,
  data = cluster_time_counts
)

time_prop <- prop.table(
  time_mat,
  margin = 1
)

pdf(
  file.path(
    fig_dir,
    "10_cluster_timepoint_heatmap_provisional_r04.pdf"
  ),
  width = 8,
  height = max(
    6,
    0.45 * nrow(time_prop) + 3
  )
)

pheatmap(
  time_prop,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  main = "Provisional microglia clusters × timepoint | row proportion"
)

dev.off()

# ------------------------------------------------------------
# 21. New cluster x sample
# ------------------------------------------------------------
cluster_sample_counts <- as.data.frame(
  table(
    cluster = obj$step3A_cluster_r04,
    sample_id = obj$sample_id
  )
)

cluster_sample_counts <- cluster_sample_counts %>%
  group_by(cluster) %>%
  mutate(
    cluster_total = sum(Freq),
    proportion_within_cluster =
      ifelse(
        cluster_total > 0,
        Freq / cluster_total,
        NA_real_
      )
  ) %>%
  ungroup()

write.csv(
  cluster_sample_counts,
  file.path(
    tab_dir,
    "09_cluster_by_sample_provisional_r04.csv"
  ),
  row.names = FALSE
)

sample_mat <- xtabs(
  Freq ~ cluster + sample_id,
  data = cluster_sample_counts
)

sample_prop <- prop.table(
  sample_mat,
  margin = 1
)

pdf(
  file.path(
    fig_dir,
    "11_cluster_sample_heatmap_provisional_r04.pdf"
  ),
  width = 11,
  height = max(
    6,
    0.45 * nrow(sample_prop) + 3
  )
)

pheatmap(
  sample_prop,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  main = "Provisional microglia clusters × biological sample | row proportion"
)

dev.off()

# ------------------------------------------------------------
# 22. New cluster x Step 2 provenance
# ------------------------------------------------------------
cluster_origin_counts <- as.data.frame(
  table(
    step3A_cluster = obj$step3A_cluster_r04,
    step2_cluster = obj$step2_final_cluster
  )
)

cluster_origin_counts <- cluster_origin_counts %>%
  group_by(step3A_cluster) %>%
  mutate(
    cluster_total = sum(Freq),
    proportion_within_step3A_cluster =
      ifelse(
        cluster_total > 0,
        Freq / cluster_total,
        NA_real_
      )
  ) %>%
  ungroup()

write.csv(
  cluster_origin_counts,
  file.path(
    tab_dir,
    "10_step3A_cluster_vs_step2_provenance.csv"
  ),
  row.names = FALSE
)

origin_mat <- xtabs(
  Freq ~ step3A_cluster + step2_cluster,
  data = cluster_origin_counts
)

origin_prop <- prop.table(
  origin_mat,
  margin = 1
)

pdf(
  file.path(
    fig_dir,
    "12_step3A_cluster_vs_step2_provenance_heatmap.pdf"
  ),
  width = 8,
  height = max(
    6,
    0.45 * nrow(origin_prop) + 3
  )
)

pheatmap(
  origin_prop,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  main = "New Step 3A cluster × Step 2 provenance | row proportion"
)

dev.off()

# ------------------------------------------------------------
# 23. Provisional cluster summary
# ------------------------------------------------------------
cluster_summary <- data.frame(
  cluster = sort(
    unique(
      as.character(
        obj$step3A_cluster_r04
      )
    )
  ),
  stringsAsFactors = FALSE
)

cluster_summary$n_cells <- vapply(
  cluster_summary$cluster,
  function(cl) {
    sum(
      as.character(
        obj$step3A_cluster_r04
      ) == cl
    )
  },
  integer(1)
)

cluster_summary$proportion <- (
  cluster_summary$n_cells /
    ncol(obj)
)

write.csv(
  cluster_summary,
  file.path(
    tab_dir,
    "11_provisional_cluster_summary_r04.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 24. Save Step 3A object
# ------------------------------------------------------------
output_rds <- file.path(
  out_dir,
  "GSE162610_step3A_microglia_reclustered_multires.rds"
)

saveRDS(
  obj,
  output_rds
)

# ------------------------------------------------------------
# 25. Formal Step 3A report
# ------------------------------------------------------------
report_file <- file.path(
  out_dir,
  "12_STEP3A_microglia_reclustering_report.txt"
)

sink(report_file)

cat("GSE162610 Step 3A | clean microglia de novo reclustering\n")
cat("====================================================\n\n")

cat("Project directory:\n")
cat(project_dir, "\n\n")

cat("Input RDS:\n")
cat(input_rds, "\n\n")

cat("Input cells: ", ncol(obj), "\n", sep = "")
cat("Input genes: ", nrow(obj), "\n", sep = "")

cat("\nFormal preprocessing:\n")
cat("- LogNormalize scale.factor = 10000\n")
cat("- HVG = 3000\n")
cat("- ScaleData = no regression\n")
cat("- Cell cycle NOT regressed\n")
cat("- timepoint NOT regressed\n")
cat("- sample_id NOT regressed in ScaleData\n")
cat("- PCA = 50 PCs\n")
cat("- Harmony dims = 1:30\n")
cat("- Harmony group = sample_id ONLY\n")
cat("- Harmony theta = 1\n")
cat("- Harmony max.iter = 20\n")
cat("- neighbors = harmony.step3A\n")
cat("- UMAP = harmony.step3A\n")
cat("- clustering graph = microglia_step3A_snn\n")

cat("\nPC summary:\n")
print(pc_summary)

cat("\nResolution scan:\n")
print(cluster_number_tbl)

cat("\nProvisional anchor resolution = 0.4\n")

cat("\nProvisional r=0.4 cluster counts:\n")
print(
  table(
    obj$step3A_cluster_r04
  )
)

cat("\nStep 2 provenance:\n")
print(
  table(
    obj$step2_final_cluster
  )
)

cat("\nTimepoint distribution:\n")
print(
  table(
    obj$timepoint
  )
)

cat("\nBiological sample distribution:\n")
print(
  table(
    obj$sample_id
  )
)

cat("\nExploratory marker availability:\n")
print(marker_availability)

cat("\nIMPORTANT INTERPRETATION RULES\n")
cat("----------------------------------------\n")
cat("1. Step 3A performs de novo within-microglia state discovery.\n")
cat("2. The old Step 2 clusters are provenance labels only.\n")
cat("3. Provisional resolution 0.4 is NOT final until clustree/markers/time/sample are reviewed.\n")
cat("4. No cluster should yet be formally named DAM.\n")
cat("5. No cluster should yet be named repair-associated or scar-associated.\n")
cat("6. Cycling populations are intentionally preserved; cell cycle was not regressed.\n")
cat("7. Any later fate trajectory must be re-checked with and without cycling cells.\n")
cat("8. Harmony corrects sample_id only and must preserve timepoint biology.\n")

cat("\nSaved Step 3A object:\n")
cat(output_rds, "\n")

cat("\nSTEP 3A COMPLETE = TRUE\n")

sink()

# ------------------------------------------------------------
# 26. Compact Step 3B handoff
# ------------------------------------------------------------
handoff_file <- file.path(
  out_dir,
  "13_STEP3B_handoff.txt"
)

sink(handoff_file)

cat("GSE162610 | Step 3B handoff\n")
cat("============================\n\n")

cat("Input for Step 3B after review:\n")
cat(output_rds, "\n\n")

cat("Step 3A purpose completed:\n")
cat("- clean microglia de novo normalization\n")
cat("- microglia-only HVG/PCA\n")
cat("- sample-level Harmony\n")
cat("- multi-resolution clustering\n")
cat("- preliminary state marker audit\n\n")

cat("Step 3B should ONLY begin after resolution/state structure is reviewed.\n")
cat("Step 3B will perform program scoring and formal state annotation.\n")

sink()

# ------------------------------------------------------------
# 27. Console output
# ------------------------------------------------------------
message("====================================================")
message("STEP 3A COMPLETE")
message("Results:")
message(out_dir)
message("")
message("Saved object:")
message(output_rds)
message("")
message("请优先上传以下 7 个文件用于 Step 3A 判断：")
message("1) figures/05_postHarmony_microglia_UMAP_overview.pdf")
message("2) figures/07_microglia_clustree_resolution_0.1_to_1.0.pdf")
message("3) figures/08_provisional_r04_microglia_UMAP.pdf")
message("4) figures/09_exploratory_state_marker_DotPlot_r04.pdf")
message("5) figures/10_cluster_timepoint_heatmap_provisional_r04.pdf")
message("6) tables/07_top20_markers_per_cluster_provisional_r04.csv")
message("7) 12_STEP3A_microglia_reclustering_report.txt")
message("")
message("Step 3A 仍然不要正式命名 DAM / repair / scar。")
message("====================================================")
