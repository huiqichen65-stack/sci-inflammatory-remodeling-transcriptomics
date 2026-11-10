# ============================================================
# Step 1B｜GSE162610 主数据集：保守 QC 过滤
# 课题：成人 SCI 后 activated/DAM-like microglia 的 repair/scar 状态分叉
#
# 基于 Step 1A 的真实 QC 分布，本步骤采用保守阈值：
#   nFeature_RNA >= 500
#   nCount_RNA   >= 1000
#   percent.mt   < 15
#   percent.HB   < 1
#
# 重要：
#   1) 暂时不设置 nFeature_RNA / nCount_RNA 的统一上限
#   2) 高 feature / 高 UMI 细胞先保留，下一步交给 doublet 检测处理
#   3) 不做 Normalize / Scale / PCA / Harmony
#   4) 本步骤只完成 QC 过滤并保存干净对象
#
# 运行方式：
#   RStudio 打开本脚本 -> Source
#   不需要终端
# ============================================================

options(stringsAsFactors = FALSE)
options(timeout = 600)

# ------------------------------------------------------------
# 0. 路径
# ------------------------------------------------------------
project_dir <- Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
in_rds <- file.path(
  project_dir,
  "02_objects",
  "step1A_GSE162610_QC_audit",
  "GSE162610_raw_with_QCmetrics.rds"
)

out_dir <- file.path(
  project_dir,
  "02_objects",
  "step1B_GSE162610_QC_filtered"
)

fig_dir <- file.path(out_dir, "figures")
tab_dir <- file.path(out_dir, "tables")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(in_rds)) {
  stop("找不到 Step 1A Seurat 对象：\n", in_rds)
}

# ------------------------------------------------------------
# 1. R 包
# ------------------------------------------------------------
pkgs <- c("Seurat", "data.table", "ggplot2", "patchwork")

miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]

if (length(miss) > 0) {
  stop(
    "缺少 R 包：", paste(miss, collapse = ", "),
    "\n请先安装后重新 Source。"
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

cat("============================================================\n")
cat("Step 1B | GSE162610 conservative QC filtering\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# 2. 读取 Step 1A 对象
# ------------------------------------------------------------
pbmc <- readRDS(in_rds)

cat("QC 前对象：", nrow(pbmc), "genes x", ncol(pbmc), "cells\n\n")

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
    "Step 1A 对象缺少以下 metadata：",
    paste(setdiff(required_meta, colnames(pbmc@meta.data)), collapse = ", ")
  )
}

# ------------------------------------------------------------
# 3. QC 阈值
# ------------------------------------------------------------
minFeature <- 500
minCount   <- 1000
maxMT      <- 15
maxHB      <- 1

thresholds <- data.frame(
  metric = c(
    "nFeature_RNA",
    "nCount_RNA",
    "percent.mt",
    "percent.HB"
  ),
  rule = c(
    paste0(">=", minFeature),
    paste0(">=", minCount),
    paste0("<", maxMT),
    paste0("<", maxHB)
  ),
  rationale = c(
    "去除极低复杂度/低质量细胞",
    "去除极低 UMI 细胞",
    "保守排除高线粒体比例细胞；避免过度删除 SCI 激活细胞",
    "去除明显血红蛋白污染；本数据绝大多数细胞接近 0"
  ),
  stringsAsFactors = FALSE
)

fwrite(
  thresholds,
  file.path(tab_dir, "01_QC_thresholds.csv")
)

cat("本步骤 QC 阈值：\n")
cat("  nFeature_RNA >=", minFeature, "\n")
cat("  nCount_RNA   >=", minCount, "\n")
cat("  percent.mt   <", maxMT, "\n")
cat("  percent.HB   <", maxHB, "\n")
cat("  不设置 nFeature / nCount 上限\n\n")

# ------------------------------------------------------------
# 4. 对每个细胞标记 QC 原因
# ------------------------------------------------------------
md <- pbmc@meta.data
md$cell <- rownames(md)

md$fail_low_feature <- md$nFeature_RNA < minFeature
md$fail_low_count   <- md$nCount_RNA < minCount
md$fail_high_mt     <- md$percent.mt >= maxMT
md$fail_high_HB     <- md$percent.HB >= maxHB

md$QC_pass <- !(
  md$fail_low_feature |
  md$fail_low_count |
  md$fail_high_mt |
  md$fail_high_HB
)

# 写回对象，便于追踪
pbmc$QC_pass <- md$QC_pass

# ------------------------------------------------------------
# 5. QC 前后：每个样本保留率
# ------------------------------------------------------------
before_sample <- as.data.frame(
  table(md$sample_id),
  stringsAsFactors = FALSE
)
colnames(before_sample) <- c("sample_id", "n_before")

after_sample <- as.data.frame(
  table(md$sample_id[md$QC_pass]),
  stringsAsFactors = FALSE
)
colnames(after_sample) <- c("sample_id", "n_after")

retention_sample <- merge(
  before_sample,
  after_sample,
  by = "sample_id",
  all.x = TRUE
)

retention_sample$n_after[is.na(retention_sample$n_after)] <- 0
retention_sample$n_removed <- (
  retention_sample$n_before - retention_sample$n_after
)
retention_sample$retention_pct <- round(
  100 * retention_sample$n_after / retention_sample$n_before,
  2
)

# 加入 timepoint
sample_time <- unique(
  md[, c("sample_id", "timepoint"), drop = FALSE]
)
sample_time$timepoint <- as.character(sample_time$timepoint)

retention_sample <- merge(
  retention_sample,
  sample_time,
  by = "sample_id",
  all.x = TRUE
)

time_levels <- c("Uninjured", "1 dpi", "3 dpi", "7 dpi")
retention_sample$time_order <- match(
  retention_sample$timepoint,
  time_levels
)
retention_sample <- retention_sample[
  order(retention_sample$time_order, retention_sample$sample_id),
]
retention_sample$time_order <- NULL

fwrite(
  retention_sample,
  file.path(tab_dir, "02_QC_retention_by_sample.csv")
)

# ------------------------------------------------------------
# 6. 每个时间点保留率
# ------------------------------------------------------------
before_time <- as.data.frame(
  table(as.character(md$timepoint)),
  stringsAsFactors = FALSE
)
colnames(before_time) <- c("timepoint", "n_before")

after_time <- as.data.frame(
  table(as.character(md$timepoint[md$QC_pass])),
  stringsAsFactors = FALSE
)
colnames(after_time) <- c("timepoint", "n_after")

retention_time <- merge(
  before_time,
  after_time,
  by = "timepoint",
  all.x = TRUE
)

retention_time$n_after[is.na(retention_time$n_after)] <- 0
retention_time$n_removed <- (
  retention_time$n_before - retention_time$n_after
)
retention_time$retention_pct <- round(
  100 * retention_time$n_after / retention_time$n_before,
  2
)

retention_time$time_order <- match(
  retention_time$timepoint,
  time_levels
)

retention_time <- retention_time[
  order(retention_time$time_order),
]
retention_time$time_order <- NULL

fwrite(
  retention_time,
  file.path(tab_dir, "03_QC_retention_by_timepoint.csv")
)

# ------------------------------------------------------------
# 7. 统计不同 QC 失败原因
#    注意：同一个细胞可以同时失败多个条件，因此各原因之和可 > 总删除数
# ------------------------------------------------------------
fail_summary <- do.call(
  rbind,
  lapply(split(md, md$sample_id), function(x) {
    data.frame(
      sample_id = unique(as.character(x$sample_id))[1],
      timepoint = unique(as.character(x$timepoint))[1],
      n_cells = nrow(x),
      fail_low_feature = sum(x$fail_low_feature),
      fail_low_count = sum(x$fail_low_count),
      fail_high_mt = sum(x$fail_high_mt),
      fail_high_HB = sum(x$fail_high_HB),
      fail_any = sum(!x$QC_pass),
      pass = sum(x$QC_pass),
      stringsAsFactors = FALSE
    )
  })
)

rownames(fail_summary) <- NULL

fail_summary$time_order <- match(
  fail_summary$timepoint,
  time_levels
)

fail_summary <- fail_summary[
  order(fail_summary$time_order, fail_summary$sample_id),
]
fail_summary$time_order <- NULL

fwrite(
  fail_summary,
  file.path(tab_dir, "04_QC_failure_reasons_by_sample.csv")
)

# ------------------------------------------------------------
# 8. 保存被删除细胞清单
# ------------------------------------------------------------
removed_cells <- md[
  !md$QC_pass,
  c(
    "cell",
    "sample_id",
    "timepoint",
    "nFeature_RNA",
    "nCount_RNA",
    "percent.mt",
    "percent.HB",
    "fail_low_feature",
    "fail_low_count",
    "fail_high_mt",
    "fail_high_HB"
  ),
  drop = FALSE
]

fwrite(
  removed_cells,
  file.path(tab_dir, "05_removed_cells_QC.csv")
)

# ------------------------------------------------------------
# 9. 正式过滤
# ------------------------------------------------------------
pbmc_qc <- subset(
  pbmc,
  cells = rownames(md)[md$QC_pass]
)

# QC_pass 已经全部 TRUE，没有必要保留在最终对象
pbmc_qc$QC_pass <- NULL

cat("QC 后对象：", nrow(pbmc_qc), "genes x", ncol(pbmc_qc), "cells\n")
cat(
  "总保留率：",
  round(100 * ncol(pbmc_qc) / ncol(pbmc), 2),
  "%\n\n"
)

# ------------------------------------------------------------
# 10. 过滤后的 QC summary
# ------------------------------------------------------------
md_after <- pbmc_qc@meta.data

qc_after_summary <- do.call(
  rbind,
  lapply(split(md_after, md_after$sample_id), function(x) {
    data.frame(
      sample_id = unique(as.character(x$sample_id))[1],
      timepoint = unique(as.character(x$timepoint))[1],
      n_cells = nrow(x),
      median_nFeature = median(x$nFeature_RNA, na.rm = TRUE),
      median_nCount = median(x$nCount_RNA, na.rm = TRUE),
      median_percent_mt = median(x$percent.mt, na.rm = TRUE),
      median_percent_HB = median(x$percent.HB, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
)

rownames(qc_after_summary) <- NULL
qc_after_summary$time_order <- match(
  qc_after_summary$timepoint,
  time_levels
)
qc_after_summary <- qc_after_summary[
  order(qc_after_summary$time_order, qc_after_summary$sample_id),
]
qc_after_summary$time_order <- NULL

fwrite(
  qc_after_summary,
  file.path(tab_dir, "06_QC_summary_by_sample_AFTER_filter.csv")
)

# ------------------------------------------------------------
# 11. QC 后小提琴图：按样本
# ------------------------------------------------------------
features_qc <- c(
  "nFeature_RNA",
  "nCount_RNA",
  "percent.mt",
  "percent.HB"
)

theme_sample <- theme(
  axis.title.x = element_blank(),
  axis.text.x = element_text(
    angle = 45,
    hjust = 1,
    size = 7
  ),
  legend.position = "none"
)

plots_sample <- lapply(features_qc, function(f) {
  VlnPlot(
    pbmc_qc,
    features = f,
    group.by = "sample_id",
    pt.size = 0
  ) + theme_sample
})

vln_sample <- wrap_plots(
  plots_sample,
  ncol = 2
)

ggsave(
  file.path(
    fig_dir,
    "01_QC_violin_by_sample_AFTER_filter.pdf"
  ),
  plot = vln_sample,
  width = 15,
  height = 9
)

# ------------------------------------------------------------
# 12. QC 后小提琴图：按时间点
# ------------------------------------------------------------
plots_time <- lapply(features_qc, function(f) {
  VlnPlot(
    pbmc_qc,
    features = f,
    group.by = "timepoint",
    pt.size = 0
  ) +
    theme(
      axis.title.x = element_blank(),
      legend.position = "none"
    )
})

vln_time <- wrap_plots(
  plots_time,
  ncol = 2
)

ggsave(
  file.path(
    fig_dir,
    "02_QC_violin_by_timepoint_AFTER_filter.pdf"
  ),
  plot = vln_time,
  width = 12,
  height = 8
)

# ------------------------------------------------------------
# 13. QC 后散点图
# ------------------------------------------------------------
p1 <- FeatureScatter(
  pbmc_qc,
  feature1 = "nCount_RNA",
  feature2 = "nFeature_RNA",
  group.by = "timepoint"
)

p2 <- FeatureScatter(
  pbmc_qc,
  feature1 = "nCount_RNA",
  feature2 = "percent.mt",
  group.by = "timepoint"
)

scatter_after <- p1 + p2

ggsave(
  file.path(
    fig_dir,
    "03_QC_feature_scatter_AFTER_filter.pdf"
  ),
  plot = scatter_after,
  width = 14,
  height = 6
)

# ------------------------------------------------------------
# 14. 每个样本保留率柱状图
# ------------------------------------------------------------
p_ret <- ggplot(
  retention_sample,
  aes(
    x = reorder(sample_id, retention_pct),
    y = retention_pct
  )
) +
  geom_col() +
  coord_flip() +
  geom_text(
    aes(label = paste0(retention_pct, "%")),
    hjust = -0.05,
    size = 3
  ) +
  ylim(
    0,
    min(
      105,
      max(retention_sample$retention_pct, na.rm = TRUE) + 5
    )
  ) +
  labs(
    x = NULL,
    y = "Cells retained (%)",
    title = "GSE162610 QC retention by sample"
  ) +
  theme_classic()

ggsave(
  file.path(
    fig_dir,
    "04_QC_retention_by_sample.pdf"
  ),
  plot = p_ret,
  width = 8,
  height = 6
)

# ------------------------------------------------------------
# 15. 保存 QC 后对象
# ------------------------------------------------------------
saveRDS(
  pbmc_qc,
  file.path(
    out_dir,
    "GSE162610_QC_filtered.rds"
  ),
  compress = FALSE
)

# ------------------------------------------------------------
# 16. 安全检查
# ------------------------------------------------------------
all_samples_preserved <- (
  length(unique(pbmc_qc$sample_id)) ==
    length(unique(pbmc$sample_id))
)

all_timepoints_preserved <- (
  length(unique(as.character(pbmc_qc$timepoint))) ==
    length(unique(as.character(pbmc$timepoint)))
)

min_retention <- min(
  retention_sample$retention_pct,
  na.rm = TRUE
)

# 如果单一样本被删太多，提示人工复核
retention_warning <- min_retention < 70

# ------------------------------------------------------------
# 17. 报告
# ------------------------------------------------------------
report <- c(
  "============================================================",
  "GSE162610 STEP 1B | CONSERVATIVE QC REPORT",
  "============================================================",
  paste0("Date: ", Sys.time()),
  "",
  "QC thresholds:",
  paste0("  nFeature_RNA >= ", minFeature),
  paste0("  nCount_RNA   >= ", minCount),
  paste0("  percent.mt   < ", maxMT),
  paste0("  percent.HB   < ", maxHB),
  "  No upper cutoff for nFeature_RNA or nCount_RNA.",
  "",
  paste0("Before QC: ", ncol(pbmc), " cells"),
  paste0("After QC:  ", ncol(pbmc_qc), " cells"),
  paste0(
    "Removed:   ",
    ncol(pbmc) - ncol(pbmc_qc),
    " cells"
  ),
  paste0(
    "Overall retention: ",
    round(100 * ncol(pbmc_qc) / ncol(pbmc), 2),
    "%"
  ),
  "",
  "Retention by sample:",
  paste0(
    "  ",
    retention_sample$sample_id,
    " | ",
    retention_sample$timepoint,
    " | ",
    retention_sample$n_before,
    " -> ",
    retention_sample$n_after,
    " | retained ",
    retention_sample$retention_pct,
    "%"
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
    "Lowest sample retention: ",
    round(min_retention, 2),
    "%"
  ),
  paste0(
    "Retention warning (<70% in any sample): ",
    retention_warning
  ),
  "",
  "IMPORTANT:",
  "High nFeature/high nCount tails were intentionally retained.",
  "They will be evaluated in the next doublet-detection step.",
  "No normalization, PCA, clustering, or Harmony was performed here.",
  "",
  paste0(
    "STEP 1B PASS = ",
    all_samples_preserved &&
      all_timepoints_preserved &&
      !retention_warning
  ),
  "",
  "NEXT:",
  "Step 1C | Doublet detection / audit.",
  "============================================================"
)

writeLines(
  report,
  file.path(
    out_dir,
    "07_STEP1B_QC_report.txt"
  )
)

cat(paste(report, collapse = "\n"))

cat("\n\n")
if (
  all_samples_preserved &&
  all_timepoints_preserved &&
  !retention_warning
) {
  cat("✅ Step 1B 完成，可以进入 doublet 检测。\n")
} else {
  cat("⚠ Step 1B 完成，但请先检查样本保留率再继续。\n")
}
