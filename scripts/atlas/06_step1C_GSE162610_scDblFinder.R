# ============================================================
# Step 1C-v2｜GSE162610：用 scDblFinder 按 sample_id 检测双细胞
#
# 说明：
#   上一版 DoubletFinder 在 ScaleData 后自行结束，截图中没有实际 Error，
#   只有“包由 R 4.5.3 构建”的版本警告；这类警告本身不会导致停止。
#
#   为避免继续耗时排查 paramSweep / DoubletFinder 版本兼容问题，
#   本版改用 scDblFinder：
#     1) 直接使用原始 counts
#     2) 按 sample_id 作为独立 capture 分别估计 doublet
#     3) 不需要 Normalize / Scale / PCA / pK sweep
#     4) 更快、更稳，且不影响后续正式 Seurat 标准化流程
#
# 输入：
#   Desktop/xjz/02_objects/step1B_GSE162610_QC_filtered/
#     GSE162610_QC_filtered.rds
#
# 输出：
#   Desktop/xjz/02_objects/step1C_GSE162610_scDblFinder/
#
# 运行：
#   RStudio 打开本脚本 -> Source
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
  "step1B_GSE162610_QC_filtered",
  "GSE162610_QC_filtered.rds"
)

out_dir <- file.path(
  project_dir,
  "02_objects",
  "step1C_GSE162610_scDblFinder"
)

fig_dir <- file.path(out_dir, "figures")
tab_dir <- file.path(out_dir, "tables")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(out_dir, "00_run_log.txt")

log_msg <- function(...) {
  x <- paste0(...)
  cat(x, "\n")
  cat(
    paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", x, "\n"),
    file = log_file,
    append = TRUE
  )
}

writeLines(
  paste0("GSE162610 Step 1C-v2 run started: ", Sys.time()),
  log_file
)

if (!file.exists(in_rds)) {
  stop("找不到 Step 1B 对象：\n", in_rds)
}

# ------------------------------------------------------------
# 1. 安装/加载包
# ------------------------------------------------------------
cran_pkgs <- c("Seurat", "data.table", "ggplot2")
missing_cran <- cran_pkgs[
  !vapply(cran_pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_cran) > 0) {
  install.packages(
    missing_cran,
    repos = "https://cloud.r-project.org"
  )
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages(
    "BiocManager",
    repos = "https://cloud.r-project.org"
  )
}

bioc_pkgs <- c(
  "SingleCellExperiment",
  "scDblFinder",
  "SummarizedExperiment"
)

missing_bioc <- bioc_pkgs[
  !vapply(bioc_pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_bioc) > 0) {
  BiocManager::install(
    missing_bioc,
    ask = FALSE,
    update = FALSE
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(ggplot2)
  library(SingleCellExperiment)
  library(scDblFinder)
  library(SummarizedExperiment)
})

log_msg("============================================================")
log_msg("Step 1C-v2 | GSE162610 scDblFinder")
log_msg("============================================================")

# ------------------------------------------------------------
# 2. 读取 Step 1B 对象
# ------------------------------------------------------------
log_msg("Reading Seurat object ...")
pbmc <- readRDS(in_rds)

log_msg(
  "Input object: ",
  nrow(pbmc), " genes x ",
  ncol(pbmc), " cells"
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

sample_ids <- sort(unique(as.character(pbmc$sample_id)))

log_msg("Detected biological samples: ", length(sample_ids))
log_msg(paste0("  ", sample_ids, collapse = "\n"))

if (length(sample_ids) != 10) {
  warning("当前不是预期的 10 个 sample_id，请检查。")
}

# ------------------------------------------------------------
# 3. 从 Seurat 提取原始 counts
# ------------------------------------------------------------
log_msg("Extracting RNA counts ...")

counts <- GetAssayData(
  object = pbmc,
  assay = "RNA",
  layer = "counts"
)

log_msg(
  "Counts matrix: ",
  nrow(counts), " genes x ",
  ncol(counts), " cells"
)

if (!identical(colnames(counts), colnames(pbmc))) {
  stop("Counts barcode 顺序与 Seurat 对象不一致。")
}

# ------------------------------------------------------------
# 4. 创建 SingleCellExperiment
# ------------------------------------------------------------
log_msg("Building SingleCellExperiment ...")

sce <- SingleCellExperiment(
  assays = list(counts = counts)
)

colData(sce)$sample_id <- as.character(pbmc$sample_id)
colData(sce)$timepoint <- as.character(pbmc$timepoint)

# 释放一个 counts 引用
rm(counts)
gc()

# ------------------------------------------------------------
# 5. scDblFinder
# ------------------------------------------------------------
# samples="sample_id"：
# 把每个生物学样本视为独立 capture，分别估计 doublet。
#
# clusters=TRUE：
# scDblFinder 内部先快速聚类，再结合 cluster-aware 方法；
# 对组织 scRNA-seq 通常比完全随机模式更稳。
#
# dbr=NULL：
# 让 scDblFinder 根据每个 capture 的细胞数自动估计 doublet rate，
# 不手工写死 7%。
# ------------------------------------------------------------
log_msg("Running scDblFinder by sample_id ...")
log_msg("This is the main computation step; please do not press Stop.")

t0 <- Sys.time()

sce <- scDblFinder(
  sce,
  samples = "sample_id",
  clusters = TRUE,
  dbr = NULL,
  verbose = TRUE
)

t1 <- Sys.time()

runtime_min <- as.numeric(
  difftime(t1, t0, units = "mins")
)

log_msg(
  "scDblFinder completed. Runtime = ",
  round(runtime_min, 2),
  " min"
)

# ------------------------------------------------------------
# 6. 检查输出列
# ------------------------------------------------------------
needed_cols <- c(
  "scDblFinder.score",
  "scDblFinder.class"
)

if (!all(needed_cols %in% colnames(colData(sce)))) {
  stop(
    "scDblFinder 完成，但没有找到预期输出列：",
    paste(
      setdiff(
        needed_cols,
        colnames(colData(sce))
      ),
      collapse = ", "
    )
  )
}

df_score <- colData(sce)$scDblFinder.score
df_class <- as.character(
  colData(sce)$scDblFinder.class
)

# scDblFinder 常见标签为 singlet / doublet
table_class <- table(df_class)

log_msg(
  "Classification: ",
  paste(
    names(table_class),
    as.integer(table_class),
    collapse = " | "
  )
)

# ------------------------------------------------------------
# 7. 映射回 Seurat
# ------------------------------------------------------------
if (!identical(colnames(sce), colnames(pbmc))) {
  stop("SCE barcode 顺序与 Seurat 对象不一致。")
}

pbmc$Double_score <- as.numeric(df_score)
pbmc$Is_Double <- ifelse(
  tolower(df_class) == "doublet",
  "Doublet",
  "Singlet"
)

# ------------------------------------------------------------
# 8. cell-level 结果
# ------------------------------------------------------------
cell_result <- data.frame(
  barcode = colnames(pbmc),
  sample_id = as.character(pbmc$sample_id),
  timepoint = as.character(pbmc$timepoint),
  nFeature_RNA = pbmc$nFeature_RNA,
  nCount_RNA = pbmc$nCount_RNA,
  percent.mt = pbmc$percent.mt,
  Double_score = pbmc$Double_score,
  Is_Double = pbmc$Is_Double,
  stringsAsFactors = FALSE
)

fwrite(
  cell_result,
  file.path(
    tab_dir,
    "01_scDblFinder_cell_level_results.csv"
  )
)

# ------------------------------------------------------------
# 9. 按样本汇总
# ------------------------------------------------------------
summary_sample <- do.call(
  rbind,
  lapply(
    split(cell_result, cell_result$sample_id),
    function(x) {

      n_total <- nrow(x)
      n_doublet <- sum(
        x$Is_Double == "Doublet",
        na.rm = TRUE
      )
      n_singlet <- sum(
        x$Is_Double == "Singlet",
        na.rm = TRUE
      )

      data.frame(
        sample_id = unique(x$sample_id)[1],
        timepoint = unique(x$timepoint)[1],
        n_total = n_total,
        n_singlet = n_singlet,
        n_doublet = n_doublet,
        detected_doublet_pct = round(
          100 * n_doublet / n_total,
          2
        ),
        median_score_singlet = median(
          x$Double_score[
            x$Is_Double == "Singlet"
          ],
          na.rm = TRUE
        ),
        median_score_doublet = median(
          x$Double_score[
            x$Is_Double == "Doublet"
          ],
          na.rm = TRUE
        ),
        median_nCount_singlet = median(
          x$nCount_RNA[
            x$Is_Double == "Singlet"
          ],
          na.rm = TRUE
        ),
        median_nCount_doublet = median(
          x$nCount_RNA[
            x$Is_Double == "Doublet"
          ],
          na.rm = TRUE
        ),
        median_nFeature_singlet = median(
          x$nFeature_RNA[
            x$Is_Double == "Singlet"
          ],
          na.rm = TRUE
        ),
        median_nFeature_doublet = median(
          x$nFeature_RNA[
            x$Is_Double == "Doublet"
          ],
          na.rm = TRUE
        ),
        stringsAsFactors = FALSE
      )
    }
  )
)

rownames(summary_sample) <- NULL

time_levels <- c(
  "Uninjured",
  "1 dpi",
  "3 dpi",
  "7 dpi"
)

summary_sample$time_order <- match(
  summary_sample$timepoint,
  time_levels
)

summary_sample <- summary_sample[
  order(
    summary_sample$time_order,
    summary_sample$sample_id
  ),
]

summary_sample$time_order <- NULL

fwrite(
  summary_sample,
  file.path(
    tab_dir,
    "02_scDblFinder_summary_by_sample.csv"
  )
)

# ------------------------------------------------------------
# 10. 按时间点汇总
# ------------------------------------------------------------
summary_time <- aggregate(
  cbind(
    n_total,
    n_singlet,
    n_doublet
  ) ~ timepoint,
  data = summary_sample,
  FUN = sum
)

summary_time$detected_doublet_pct <- round(
  100 *
    summary_time$n_doublet /
    summary_time$n_total,
  2
)

summary_time$time_order <- match(
  summary_time$timepoint,
  time_levels
)

summary_time <- summary_time[
  order(summary_time$time_order),
]

summary_time$time_order <- NULL

fwrite(
  summary_time,
  file.path(
    tab_dir,
    "03_scDblFinder_summary_by_timepoint.csv"
  )
)

# ------------------------------------------------------------
# 11. 图 1：每个样本 doublet rate
# ------------------------------------------------------------
p_rate <- ggplot(
  summary_sample,
  aes(
    x = reorder(
      sample_id,
      detected_doublet_pct
    ),
    y = detected_doublet_pct
  )
) +
  geom_col() +
  coord_flip() +
  geom_text(
    aes(
      label = paste0(
        detected_doublet_pct,
        "%"
      )
    ),
    hjust = -0.08,
    size = 3
  ) +
  expand_limits(
    y = max(
      summary_sample$detected_doublet_pct
    ) + 1
  ) +
  labs(
    x = NULL,
    y = "Detected doublets (%)",
    title = "GSE162610 scDblFinder rate by sample"
  ) +
  theme_classic()

ggsave(
  file.path(
    fig_dir,
    "01_doublet_rate_by_sample.pdf"
  ),
  p_rate,
  width = 8,
  height = 6
)

# ------------------------------------------------------------
# 12. 图 2：Singlet vs Doublet QC
# ------------------------------------------------------------
plot_qc <- data.frame(
  Is_Double = pbmc$Is_Double,
  nCount_RNA = pbmc$nCount_RNA,
  nFeature_RNA = pbmc$nFeature_RNA,
  stringsAsFactors = FALSE
)

p_count <- ggplot(
  plot_qc,
  aes(
    x = Is_Double,
    y = nCount_RNA
  )
) +
  geom_violin(
    scale = "width",
    trim = TRUE
  ) +
  labs(
    x = NULL,
    y = "nCount_RNA",
    title = "UMI counts"
  ) +
  theme_classic()

p_feature <- ggplot(
  plot_qc,
  aes(
    x = Is_Double,
    y = nFeature_RNA
  )
) +
  geom_violin(
    scale = "width",
    trim = TRUE
  ) +
  labs(
    x = NULL,
    y = "nFeature_RNA",
    title = "Detected genes"
  ) +
  theme_classic()

if (requireNamespace("patchwork", quietly = TRUE)) {
  p_qc <- p_count + p_feature

  ggsave(
    file.path(
      fig_dir,
      "02_doublet_vs_singlet_QC_violin.pdf"
    ),
    p_qc,
    width = 10,
    height = 5
  )
} else {
  ggsave(
    file.path(
      fig_dir,
      "02A_doublet_vs_singlet_nCount.pdf"
    ),
    p_count,
    width = 5,
    height = 5
  )

  ggsave(
    file.path(
      fig_dir,
      "02B_doublet_vs_singlet_nFeature.pdf"
    ),
    p_feature,
    width = 5,
    height = 5
  )
}

# ------------------------------------------------------------
# 13. 图 3：score distribution by sample
# ------------------------------------------------------------
p_score <- ggplot(
  cell_result,
  aes(
    x = Double_score,
    fill = Is_Double
  )
) +
  geom_histogram(
    bins = 70,
    position = "identity",
    alpha = 0.55
  ) +
  facet_wrap(
    ~ sample_id,
    scales = "free_y"
  ) +
  labs(
    x = "scDblFinder score",
    y = "Cell count",
    title = "scDblFinder score distribution by sample",
    fill = NULL
  ) +
  theme_classic()

ggsave(
  file.path(
    fig_dir,
    "03_doublet_score_distribution_by_sample.pdf"
  ),
  p_score,
  width = 14,
  height = 9
)

# ------------------------------------------------------------
# 14. 保存带标签的完整对象
# ------------------------------------------------------------
saveRDS(
  pbmc,
  file.path(
    out_dir,
    "GSE162610_QC_with_scDblFinder_labels.rds"
  ),
  compress = FALSE
)

# ------------------------------------------------------------
# 15. 去除 doublet
# ------------------------------------------------------------
pbmc_singlet <- subset(
  pbmc,
  subset = Is_Double == "Singlet"
)

saveRDS(
  pbmc_singlet,
  file.path(
    out_dir,
    "GSE162610_QC_singlets.rds"
  ),
  compress = FALSE
)

# ------------------------------------------------------------
# 16. 最终检查
# ------------------------------------------------------------
all_samples_preserved <- (
  length(unique(pbmc_singlet$sample_id)) ==
    length(unique(pbmc$sample_id))
)

all_timepoints_preserved <- (
  length(
    unique(as.character(pbmc_singlet$timepoint))
  ) ==
    length(
      unique(as.character(pbmc$timepoint))
    )
)

overall_doublet_pct <- round(
  100 *
    (ncol(pbmc) - ncol(pbmc_singlet)) /
    ncol(pbmc),
  2
)

# 宽松审计警戒线
rate_warning <- any(
  summary_sample$detected_doublet_pct > 15 |
    summary_sample$detected_doublet_pct < 0.05
)

# ------------------------------------------------------------
# 17. 报告
# ------------------------------------------------------------
report <- c(
  "============================================================",
  "GSE162610 STEP 1C-v2 | scDblFinder REPORT",
  "============================================================",
  paste0("Date: ", Sys.time()),
  paste0("Runtime (min): ", round(runtime_min, 2)),
  "",
  paste0("Input cells: ", ncol(pbmc)),
  paste0("Singlets retained: ", ncol(pbmc_singlet)),
  paste0(
    "Doublets removed: ",
    ncol(pbmc) - ncol(pbmc_singlet)
  ),
  paste0(
    "Overall detected doublet rate: ",
    overall_doublet_pct,
    "%"
  ),
  "",
  "Doublet rate by sample:",
  paste0(
    "  ",
    summary_sample$sample_id,
    " | ",
    summary_sample$timepoint,
    " | total=",
    summary_sample$n_total,
    " | doublet=",
    summary_sample$n_doublet,
    " | rate=",
    summary_sample$detected_doublet_pct,
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
    "Rate warning (>15% or <0.05% in any sample): ",
    rate_warning
  ),
  "",
  "IMPORTANT:",
  "scDblFinder was run with sample_id as independent capture.",
  "No formal downstream normalization/Harmony has been performed.",
  "",
  paste0(
    "STEP 1C PASS = ",
    all_samples_preserved &&
      all_timepoints_preserved &&
      !rate_warning
  ),
  "",
  "NEXT:",
  "Step 1D | Formal LogNormalize + HVG + PCA before Harmony.",
  "============================================================"
)

writeLines(
  report,
  file.path(
    out_dir,
    "04_STEP1C_scDblFinder_report.txt"
  )
)

log_msg("Final report written.")
log_msg(
  "Singlets retained: ",
  ncol(pbmc_singlet)
)
log_msg(
  "Overall doublet rate: ",
  overall_doublet_pct,
  "%"
)

cat("\n")
cat(paste(report, collapse = "\n"))
cat("\n\n")

if (
  all_samples_preserved &&
  all_timepoints_preserved &&
  !rate_warning
) {
  cat("✅ Step 1C-v2 完成，请把结果发给我判断。\n")
} else {
  cat("⚠ Step 1C-v2 完成，但有结果需要人工复核。\n")
}
