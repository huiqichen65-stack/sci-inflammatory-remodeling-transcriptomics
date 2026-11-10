# ============================================================
# Step 1D｜GSE162610：正式 LogNormalize + HVG + ScaleData + PCA（Harmony 前）
# 课题：成人 SCI 后 activated/DAM-like microglia 的 repair/scar 状态分叉
#
# 基于你上传的 7Normalize 模板改写，但做以下调整：
#   1) 不做 CellCycle 回归（本课题前面已决定跳过 cell cycle）
#   2) 不使用 SCTransform，主线采用 LogNormalize
#   3) 使用 3000 个高变基因，保留更多损伤后细胞状态差异
#   4) ScaleData 不回归 S.Score/G2M.Score，也不回归 timepoint
#   5) 正式 RunPCA 50 PCs
#   6) 输出 PCA 方差解释、ElbowPlot
#   7) 额外生成 Harmony 前 UMAP（仅用于观察 sample batch / timepoint 分布）
#   8) 本步骤不做 Harmony、不做正式聚类
#
# 输入：
#   Desktop/xjz/02_objects/step1C_GSE162610_scDblFinder/
#     GSE162610_QC_singlets.rds
#
# 输出：
#   Desktop/xjz/02_objects/step1D_GSE162610_LogNorm_PCA_preHarmony/
#
# 运行方式：
#   RStudio 打开本脚本 -> Source
#   不需要终端
# ============================================================

options(stringsAsFactors = FALSE)
options(timeout = 1200)

set.seed(20260730)

# ------------------------------------------------------------
# 0. 路径
# ------------------------------------------------------------
project_dir <- Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
in_rds <- file.path(
  project_dir,
  "02_objects",
  "step1C_GSE162610_scDblFinder",
  "GSE162610_QC_singlets.rds"
)

out_dir <- file.path(
  project_dir,
  "02_objects",
  "step1D_GSE162610_LogNorm_PCA_preHarmony"
)

fig_dir <- file.path(out_dir, "figures")
tab_dir <- file.path(out_dir, "tables")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(out_dir, "00_run_log.txt")

log_msg <- function(...) {
  msg <- paste0(...)
  cat(msg, "\n")
  cat(
    paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", msg, "\n"),
    file = log_file,
    append = TRUE
  )
}

writeLines(
  paste0("GSE162610 Step 1D started: ", Sys.time()),
  log_file
)

if (!file.exists(in_rds)) {
  stop("找不到 Step 1C singlet 对象：\n", in_rds)
}

# ------------------------------------------------------------
# 1. R 包
# ------------------------------------------------------------
pkgs <- c(
  "Seurat",
  "data.table",
  "ggplot2",
  "patchwork"
)

missing_pkgs <- pkgs[
  !vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_pkgs) > 0) {
  stop(
    "缺少 R 包：",
    paste(missing_pkgs, collapse = ", "),
    "\n请先安装后重新 Source。"
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

log_msg("============================================================")
log_msg("Step 1D | Formal LogNormalize + PCA before Harmony")
log_msg("============================================================")

# ------------------------------------------------------------
# 2. 读取 singlet 对象
# ------------------------------------------------------------
log_msg("Reading singlet Seurat object ...")

pbmc <- readRDS(in_rds)

log_msg(
  "Input object: ",
  nrow(pbmc),
  " genes x ",
  ncol(pbmc),
  " singlets"
)

required_meta <- c(
  "sample_id",
  "timepoint",
  "nFeature_RNA",
  "nCount_RNA",
  "percent.mt",
  "percent.HB"
)

if (!all(required_meta %in% colnames(pbmc@meta.data))) {
  stop(
    "对象缺少 metadata：",
    paste(
      setdiff(required_meta, colnames(pbmc@meta.data)),
      collapse = ", "
    )
  )
}

# 确保使用 RNA assay
DefaultAssay(pbmc) <- "RNA"

# ------------------------------------------------------------
# 3. 正式 LogNormalize
# ------------------------------------------------------------
log_msg("Running NormalizeData(LogNormalize, scale.factor=10000) ...")

t0 <- Sys.time()

pbmc <- NormalizeData(
  pbmc,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

log_msg("NormalizeData finished.")

# ------------------------------------------------------------
# 4. 高变基因：3000
# ------------------------------------------------------------
log_msg("Finding 3000 highly variable genes ...")

pbmc <- FindVariableFeatures(
  pbmc,
  selection.method = "vst",
  nfeatures = 3000,
  verbose = FALSE
)

hvg <- VariableFeatures(pbmc)

if (length(hvg) != 3000) {
  warning(
    "实际识别到的 HVG 数为 ",
    length(hvg),
    "，不是预期的 3000。"
  )
}

writeLines(
  hvg,
  file.path(
    tab_dir,
    "01_HVG_3000_genes.txt"
  )
)

# HVG 信息表
vf_info <- HVFInfo(pbmc)

vf_info$gene <- rownames(vf_info)
vf_info$is_HVG <- vf_info$gene %in% hvg

fwrite(
  vf_info,
  file.path(
    tab_dir,
    "02_variable_feature_statistics.csv"
  )
)

log_msg("HVG finished: ", length(hvg), " genes.")

# ------------------------------------------------------------
# 5. ScaleData
# ------------------------------------------------------------
# 不回归：
#   S.Score / G2M.Score
#   timepoint
#   sample_id
#   nCount_RNA
#
# 理由：
#   timepoint 是我们真正需要保留的 SCI 生物学信号；
#   sample batch 下一步由 Harmony 处理；
#   不希望把 activated/DAM-like 状态相关变化一起回归掉。
#
# 仅缩放 3000 HVGs，节省内存和时间。
# ------------------------------------------------------------
log_msg("Scaling 3000 HVGs without cell-cycle regression ...")

pbmc <- ScaleData(
  pbmc,
  features = hvg,
  verbose = FALSE
)

log_msg("ScaleData finished.")

# ------------------------------------------------------------
# 6. PCA：正式计算 50 PCs
# ------------------------------------------------------------
log_msg("Running PCA with 50 PCs ...")

pbmc <- RunPCA(
  pbmc,
  features = hvg,
  npcs = 50,
  verbose = FALSE
)

log_msg("PCA finished.")

if (!"pca" %in% Reductions(pbmc)) {
  stop("PCA reduction 未成功生成。")
}

# ------------------------------------------------------------
# 7. PCA 方差解释
# ------------------------------------------------------------
stdev <- Stdev(pbmc, reduction = "pca")

variance_pct <- 100 * (stdev^2) / sum(stdev^2)
cumulative_pct <- cumsum(variance_pct)

pca_variance <- data.frame(
  PC = seq_along(stdev),
  stdev = stdev,
  variance_pct = variance_pct,
  cumulative_variance_pct = cumulative_pct
)

fwrite(
  pca_variance,
  file.path(
    tab_dir,
    "03_PCA_variance_explained.csv"
  )
)

# 自动给一个“候选 PC 数”，仅供参考，不直接替代人工看 Elbow
# 规则：
#   取累计解释达到 80% 的 PC；
#   限制在 20–40 PCs 范围内。
pc80 <- which(cumulative_pct >= 80)[1]

if (is.na(pc80)) {
  pc80 <- 30
}

candidate_pc <- max(20, min(40, pc80))

writeLines(
  c(
    paste0("PC reaching 80% cumulative variance: ", pc80),
    paste0("Candidate PCs for Harmony/UMAP: 1:", candidate_pc),
    "",
    "IMPORTANT:",
    "This is only an automatic reference.",
    "Final PC number should be decided after reviewing ElbowPlot and biological separation."
  ),
  file.path(
    tab_dir,
    "04_candidate_PC_number.txt"
  )
)

log_msg(
  "Automatic candidate PC number = ",
  candidate_pc
)

# ------------------------------------------------------------
# 8. ElbowPlot：独立 PDF
# ------------------------------------------------------------
p_elbow <- ElbowPlot(
  pbmc,
  ndims = 50
) +
  ggtitle(
    "GSE162610 PCA elbow plot (pre-Harmony)"
  )

ggsave(
  file.path(
    fig_dir,
    "01_PCA_ElbowPlot_50PCs.pdf"
  ),
  plot = p_elbow,
  width = 7,
  height = 5
)

# ------------------------------------------------------------
# 9. PCA scatter：PC1 vs PC2
# ------------------------------------------------------------
p_pca_sample <- DimPlot(
  pbmc,
  reduction = "pca",
  group.by = "sample_id",
  dims = c(1, 2),
  pt.size = 0.15
) +
  ggtitle(
    "PCA before Harmony - colored by sample"
  )

ggsave(
  file.path(
    fig_dir,
    "02_PCA_PC1_PC2_by_sample.pdf"
  ),
  plot = p_pca_sample,
  width = 10,
  height = 7
)

p_pca_time <- DimPlot(
  pbmc,
  reduction = "pca",
  group.by = "timepoint",
  dims = c(1, 2),
  pt.size = 0.15
) +
  ggtitle(
    "PCA before Harmony - colored by timepoint"
  )

ggsave(
  file.path(
    fig_dir,
    "03_PCA_PC1_PC2_by_timepoint.pdf"
  ),
  plot = p_pca_time,
  width = 9,
  height = 7
)

# ------------------------------------------------------------
# 10. Harmony 前诊断 UMAP
# ------------------------------------------------------------
# 注意：
#   这不是最终 UMAP。
#   只是看看 sample batch 在校正前有多明显。
#   下一步 Harmony 后会重新跑正式 UMAP。
# ------------------------------------------------------------
dims_preHarmony <- 1:candidate_pc

log_msg(
  "Running diagnostic pre-Harmony UMAP using PCs 1:",
  candidate_pc,
  " ..."
)

pbmc <- RunUMAP(
  pbmc,
  reduction = "pca",
  dims = dims_preHarmony,
  reduction.name = "umap.preharmony",
  reduction.key = "UMAPPRE_",
  seed.use = 20260730,
  verbose = FALSE
)

p_umap_sample <- DimPlot(
  pbmc,
  reduction = "umap.preharmony",
  group.by = "sample_id",
  pt.size = 0.10
) +
  ggtitle(
    "Pre-Harmony UMAP - sample_id"
  )

ggsave(
  file.path(
    fig_dir,
    "04_preHarmony_UMAP_by_sample.pdf"
  ),
  plot = p_umap_sample,
  width = 10,
  height = 7
)

p_umap_time <- DimPlot(
  pbmc,
  reduction = "umap.preharmony",
  group.by = "timepoint",
  pt.size = 0.10
) +
  ggtitle(
    "Pre-Harmony UMAP - timepoint"
  )

ggsave(
  file.path(
    fig_dir,
    "05_preHarmony_UMAP_by_timepoint.pdf"
  ),
  plot = p_umap_time,
  width = 9,
  height = 7
)

# ------------------------------------------------------------
# 11. Top PCA loadings
# ------------------------------------------------------------
loadings <- Loadings(
  pbmc,
  reduction = "pca"
)

n_top <- 30
pcs_to_report <- 1:min(20, ncol(loadings))

top_loadings_list <- lapply(
  pcs_to_report,
  function(pc) {

    vals <- loadings[, pc]

    ord <- order(
      abs(vals),
      decreasing = TRUE
    )[1:min(n_top, length(vals))]

    data.frame(
      PC = paste0("PC", pc),
      gene = rownames(loadings)[ord],
      loading = vals[ord],
      abs_loading = abs(vals[ord]),
      stringsAsFactors = FALSE
    )
  }
)

top_loadings <- do.call(
  rbind,
  top_loadings_list
)

fwrite(
  top_loadings,
  file.path(
    tab_dir,
    "05_top30_gene_loadings_PC1_to_PC20.csv"
  )
)

# ------------------------------------------------------------
# 12. 样本/时间点细胞数确认
# ------------------------------------------------------------
sample_counts <- as.data.frame(
  table(pbmc$sample_id),
  stringsAsFactors = FALSE
)

colnames(sample_counts) <- c(
  "sample_id",
  "n_cells"
)

fwrite(
  sample_counts,
  file.path(
    tab_dir,
    "06_cell_counts_by_sample_preHarmony.csv"
  )
)

time_counts <- as.data.frame(
  table(pbmc$timepoint),
  stringsAsFactors = FALSE
)

colnames(time_counts) <- c(
  "timepoint",
  "n_cells"
)

fwrite(
  time_counts,
  file.path(
    tab_dir,
    "07_cell_counts_by_timepoint_preHarmony.csv"
  )
)

# ------------------------------------------------------------
# 13. 保存正式 pre-Harmony 对象
# ------------------------------------------------------------
saveRDS(
  pbmc,
  file.path(
    out_dir,
    "GSE162610_LogNorm_PCA_preHarmony.rds"
  ),
  compress = FALSE
)

# ------------------------------------------------------------
# 14. 最终检查
# ------------------------------------------------------------
all_samples_preserved <- length(
  unique(as.character(pbmc$sample_id))
) == 10

all_timepoints_preserved <- length(
  unique(as.character(pbmc$timepoint))
) == 4

pca_ok <- (
  "pca" %in% Reductions(pbmc) &&
    ncol(Embeddings(pbmc, "pca")) >= 30
)

umap_ok <- (
  "umap.preharmony" %in% Reductions(pbmc)
)

step_pass <- all(
  all_samples_preserved,
  all_timepoints_preserved,
  length(hvg) >= 2000,
  pca_ok,
  umap_ok
)

runtime_min <- as.numeric(
  difftime(
    Sys.time(),
    t0,
    units = "mins"
  )
)

# ------------------------------------------------------------
# 15. 报告
# ------------------------------------------------------------
report <- c(
  "============================================================",
  "GSE162610 STEP 1D | LOGNORMALIZE + PCA PRE-HARMONY REPORT",
  "============================================================",
  paste0("Date: ", Sys.time()),
  paste0("Runtime (min): ", round(runtime_min, 2)),
  "",
  paste0(
    "Input singlets: ",
    ncol(pbmc)
  ),
  paste0(
    "Genes: ",
    nrow(pbmc)
  ),
  "",
  "Normalization:",
  "  LogNormalize",
  "  scale.factor = 10000",
  "",
  paste0(
    "Variable features: ",
    length(hvg)
  ),
  "Cell-cycle regression: NO",
  "S.Score/G2M.Score regression: NO",
  "Timepoint regression: NO",
  "sample_id regression: NO",
  "",
  paste0(
    "PCA components calculated: ",
    ncol(
      Embeddings(
        pbmc,
        "pca"
      )
    )
  ),
  paste0(
    "PC reaching 80% cumulative variance: ",
    pc80
  ),
  paste0(
    "Automatic candidate PCs for next step: 1:",
    candidate_pc
  ),
  "",
  paste0(
    "All 10 samples preserved: ",
    all_samples_preserved
  ),
  paste0(
    "All 4 timepoints preserved: ",
    all_timepoints_preserved
  ),
  paste0(
    "PCA successful: ",
    pca_ok
  ),
  paste0(
    "Diagnostic pre-Harmony UMAP successful: ",
    umap_ok
  ),
  "",
  "IMPORTANT:",
  "The pre-Harmony UMAP is diagnostic only.",
  "Do not use it as the final clustering/annotation embedding.",
  "Final Harmony dimensionality should be chosen after reviewing ElbowPlot and PCA variance table.",
  "",
  paste0(
    "STEP 1D COMPLETE = ",
    step_pass
  ),
  "",
  "NEXT:",
  "Step 1E | Harmony batch correction by sample_id + formal UMAP/neighbor graph/clustering.",
  "============================================================"
)

writeLines(
  report,
  file.path(
    out_dir,
    "08_STEP1D_LogNorm_PCA_report.txt"
  )
)

log_msg(
  "Step 1D complete = ",
  step_pass
)

cat("\n")
cat(paste(report, collapse = "\n"))
cat("\n\n")

if (step_pass) {
  cat(
    "✅ Step 1D 完成。请先把 PCA/UMAP 结果发给我判断，再进入 Harmony。\n"
  )
} else {
  cat(
    "⚠ Step 1D 有检查项未通过，请先查看报告。\n"
  )
}
