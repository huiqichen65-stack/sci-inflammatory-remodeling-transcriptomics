# ============================================================
# Step 0A｜GSE162610 主数据集审计 + sample manifest
# 课题：成人 SCI 后 activated/DAM-like microglia 的 repair/scar 状态分叉
# 目标：
#   1) 检查 GSE162610 5 个核心文件是否齐全
#   2) 不加载完整表达矩阵，仅读取 MatrixMarket 头部检查矩阵维度
#   3) 读取 barcode / gene / barcode_metadata / gene_metadata
#   4) 自动识别样本与时间信息，生成统一 sample manifest
#   5) 输出 Step 0 审计报告，为下一步正式读入 Seurat 做准备
#
# 运行方式：
#   直接在 RStudio 打开本脚本，点击 Source 即可
#   不需要终端
# ============================================================

options(stringsAsFactors = FALSE)
options(timeout = 600)

# -----------------------------
# 0. 项目路径
# -----------------------------
project_dir <- Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
raw_dir     <- file.path(project_dir, "01_raw_data", "GSE162610")
admin_dir   <- file.path(project_dir, "00_admin", "step0_GSE162610_audit")
dir.create(admin_dir, recursive = TRUE, showWarnings = FALSE)

cat("============================================================\n")
cat("Step 0A | GSE162610 主数据集审计\n")
cat("项目目录：", project_dir, "\n")
cat("原始数据目录：", raw_dir, "\n")
cat("输出目录：", admin_dir, "\n")
cat("============================================================\n\n")

# -----------------------------
# 1. 依赖
# -----------------------------
need_pkgs <- c("data.table")
to_install <- need_pkgs[!vapply(need_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(to_install) > 0) {
  install.packages(to_install, repos = "https://cloud.r-project.org")
}
suppressPackageStartupMessages(library(data.table))

# -----------------------------
# 2. 检查核心文件
# -----------------------------
files <- c(
  barcode_metadata = file.path(raw_dir, "GSE162610_barcode_metadata.tsv.gz"),
  barcodes         = file.path(raw_dir, "GSE162610_barcodes.tsv.gz"),
  gene_metadata    = file.path(raw_dir, "GSE162610_gene_metadata.tsv.gz"),
  genes            = file.path(raw_dir, "GSE162610_genes.tsv.gz"),
  matrix           = file.path(raw_dir, "GSE162610_sci_mat.mtx.gz")
)

file_check <- data.frame(
  item = names(files),
  path = unname(files),
  exists = file.exists(files),
  size_MB = ifelse(file.exists(files), round(file.info(files)$size / 1024^2, 3), NA_real_)
)

fwrite(file_check, file.path(admin_dir, "01_file_check.csv"))

if (!all(file_check$exists)) {
  stop(
    "以下核心文件缺失：\n",
    paste(file_check$path[!file_check$exists], collapse = "\n")
  )
}

cat("5 个核心文件均存在。\n\n")

# -----------------------------
# 3. 安全读取 gzipped tsv
# -----------------------------
read_tsv_gz <- function(path, header = TRUE) {
  data.table::fread(
    cmd = paste("gzip -dc", shQuote(path)),
    sep = "\t",
    header = header,
    data.table = FALSE,
    showProgress = FALSE
  )
}

# Windows 上如果 gzip 命令不可用，则 fallback 到 gzfile + read.delim
read_tsv_gz_safe <- function(path, header = TRUE) {
  out <- try(
    read_tsv_gz(path, header = header),
    silent = TRUE
  )
  if (!inherits(out, "try-error")) return(out)

  con <- gzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)
  read.delim(
    con,
    header = header,
    sep = "\t",
    quote = "",
    comment.char = "",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

# -----------------------------
# 4. 读取 barcodes / genes
# -----------------------------
barcodes <- read_tsv_gz_safe(files["barcodes"], header = FALSE)
genes    <- read_tsv_gz_safe(files["genes"], header = FALSE)

colnames(barcodes)[1] <- "barcode"
if (ncol(genes) == 1) {
  colnames(genes) <- "gene"
} else {
  colnames(genes)[1:2] <- c("gene_id", "gene_symbol")
}

cat("barcodes 数：", nrow(barcodes), "\n")
cat("genes 数：", nrow(genes), "\n\n")

# -----------------------------
# 5. 读取 metadata
# -----------------------------
barcode_meta <- read_tsv_gz_safe(files["barcode_metadata"], header = TRUE)
gene_meta    <- read_tsv_gz_safe(files["gene_metadata"], header = TRUE)

cat("barcode metadata：", nrow(barcode_meta), "行 ×", ncol(barcode_meta), "列\n")
cat("gene metadata：", nrow(gene_meta), "行 ×", ncol(gene_meta), "列\n\n")

# 保存原始字段名与前几行，方便人工检查
writeLines(colnames(barcode_meta), file.path(admin_dir, "02_barcode_metadata_columns.txt"))
writeLines(colnames(gene_meta),    file.path(admin_dir, "03_gene_metadata_columns.txt"))
fwrite(head(barcode_meta, 30), file.path(admin_dir, "04_barcode_metadata_head30.csv"))
fwrite(head(gene_meta, 30),    file.path(admin_dir, "05_gene_metadata_head30.csv"))

# -----------------------------
# 6. 读取 MatrixMarket 头部，不加载整个表达矩阵
# -----------------------------
read_mtx_dims_gz <- function(path) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con))

  first <- readLines(con, n = 1)
  if (!grepl("^%%MatrixMarket", first)) {
    stop("矩阵文件不是标准 MatrixMarket 格式：", path)
  }

  repeat {
    x <- readLines(con, n = 1)
    if (length(x) == 0) stop("未找到 MatrixMarket 维度行。")
    if (!grepl("^%", x)) break
  }

  vals <- strsplit(trimws(x), "\\s+")[[1]]
  if (length(vals) < 3) stop("MatrixMarket 维度行格式异常：", x)

  c(
    n_rows = as.integer(vals[1]),
    n_cols = as.integer(vals[2]),
    n_nonzero = as.double(vals[3])
  )
}

mtx_dims <- read_mtx_dims_gz(files["matrix"])

dim_check <- data.frame(
  metric = c("matrix_n_rows", "matrix_n_cols", "matrix_n_nonzero",
             "genes_n", "barcodes_n", "barcode_metadata_n", "gene_metadata_n"),
  value = c(
    mtx_dims["n_rows"],
    mtx_dims["n_cols"],
    mtx_dims["n_nonzero"],
    nrow(genes),
    nrow(barcodes),
    nrow(barcode_meta),
    nrow(gene_meta)
  )
)
fwrite(dim_check, file.path(admin_dir, "06_dimension_check.csv"))

cat("MatrixMarket 维度：",
    mtx_dims["n_rows"], "genes ×",
    mtx_dims["n_cols"], "cells；non-zero =",
    format(mtx_dims["n_nonzero"], scientific = FALSE), "\n\n")

# -----------------------------
# 7. 基础维度一致性
# -----------------------------
genes_match    <- as.integer(mtx_dims["n_rows"]) == nrow(genes)
barcodes_match <- as.integer(mtx_dims["n_cols"]) == nrow(barcodes)

# barcode metadata 可能是一行一个细胞，也可能有额外格式，因此这里只先记录
meta_same_as_cells <- nrow(barcode_meta) == nrow(barcodes)
gene_meta_same     <- nrow(gene_meta) == nrow(genes)

basic_check <- data.frame(
  check = c(
    "matrix rows == genes",
    "matrix columns == barcodes",
    "barcode_metadata rows == barcodes",
    "gene_metadata rows == genes"
  ),
  pass = c(
    genes_match,
    barcodes_match,
    meta_same_as_cells,
    gene_meta_same
  )
)
fwrite(basic_check, file.path(admin_dir, "07_basic_consistency_check.csv"))

# -----------------------------
# 8. 自动寻找 barcode 列
# -----------------------------
norm_name <- function(x) {
  tolower(gsub("[^a-z0-9]+", "_", x))
}
meta_cols_norm <- norm_name(colnames(barcode_meta))

barcode_candidates <- which(
  meta_cols_norm %in% c(
    "barcode", "barcodes", "cell", "cell_id", "cellid",
    "cell_barcode", "cellbarcode", "cell_name", "cellname"
  ) |
    grepl("barcode", meta_cols_norm)
)

barcode_col <- NA_character_
barcode_overlap <- NA_integer_

if (length(barcode_candidates) > 0) {
  overlaps <- vapply(barcode_candidates, function(i) {
    sum(as.character(barcode_meta[[i]]) %in% as.character(barcodes$barcode))
  }, integer(1))
  best <- barcode_candidates[which.max(overlaps)]
  barcode_col <- colnames(barcode_meta)[best]
  barcode_overlap <- max(overlaps)
}

# 如果找不到明显 barcode 列，再扫描前若干字符列
if (is.na(barcode_col)) {
  char_cols <- which(vapply(barcode_meta, function(x) is.character(x) || is.factor(x), logical(1)))
  if (length(char_cols) > 0) {
    overlaps <- vapply(char_cols, function(i) {
      sum(as.character(barcode_meta[[i]]) %in% as.character(barcodes$barcode))
    }, integer(1))
    if (max(overlaps) > 0) {
      best <- char_cols[which.max(overlaps)]
      barcode_col <- colnames(barcode_meta)[best]
      barcode_overlap <- max(overlaps)
    }
  }
}

# -----------------------------
# 9. 合并 cell metadata
# -----------------------------
if (!is.na(barcode_col) && barcode_overlap > 0) {
  cell_meta <- merge(
    data.frame(barcode = as.character(barcodes$barcode), stringsAsFactors = FALSE),
    barcode_meta,
    by.x = "barcode",
    by.y = barcode_col,
    all.x = TRUE,
    sort = FALSE
  )
} else if (nrow(barcode_meta) == nrow(barcodes)) {
  cell_meta <- cbind(
    data.frame(barcode = as.character(barcodes$barcode), stringsAsFactors = FALSE),
    barcode_meta
  )
  barcode_col <- "[row-order matched]"
  barcode_overlap <- nrow(barcodes)
} else {
  cell_meta <- data.frame(barcode = as.character(barcodes$barcode), stringsAsFactors = FALSE)
  warning(
    "无法可靠把 barcode_metadata 与 barcodes 对齐。\n",
    "请把 04_barcode_metadata_head30.csv 发给我，我会按真实格式改脚本。"
  )
}

fwrite(cell_meta, file.path(admin_dir, "08_cell_metadata_merged.csv"))

# -----------------------------
# 10. 自动寻找样本/时间相关列
# -----------------------------
nm <- norm_name(colnames(cell_meta))

sample_idx <- which(
  nm %in% c("sample", "sample_id", "sampleid", "orig_ident", "library", "replicate") |
    grepl("sample|library|replicate|orig_ident", nm)
)

time_idx <- which(
  nm %in% c("time", "timepoint", "time_point", "dpi", "day", "days", "condition", "group") |
    grepl("time|dpi|day|condition|group|injury", nm)
)

possible_sample_cols <- colnames(cell_meta)[sample_idx]
possible_time_cols   <- colnames(cell_meta)[time_idx]

writeLines(
  c(
    paste0("barcode column: ", barcode_col),
    paste0("barcode overlap: ", barcode_overlap),
    "",
    "Possible sample columns:",
    if (length(possible_sample_cols)) possible_sample_cols else "[none]",
    "",
    "Possible time/group columns:",
    if (length(possible_time_cols)) possible_time_cols else "[none]"
  ),
  file.path(admin_dir, "09_detected_metadata_fields.txt")
)

# -----------------------------
# 11. 自动构建 sample manifest
# -----------------------------
# 优先选最像 sample 的列
choose_best_col <- function(df, cols, role = "sample") {
  if (length(cols) == 0) return(NA_character_)

  # 倾向选择非唯一、但有多个类别且不是每个 cell 唯一的列
  scores <- sapply(cols, function(cc) {
    x <- df[[cc]]
    nuniq <- length(unique(na.omit(x)))
    nonmiss <- sum(!is.na(x))
    if (nonmiss == 0) return(-Inf)

    score <- 0
    if (role == "sample") {
      if (nuniq >= 2 && nuniq <= 100) score <- score + 5
      if (grepl("sample|orig|library|replicate", norm_name(cc))) score <- score + 3
    } else {
      if (nuniq >= 2 && nuniq <= 30) score <- score + 5
      if (grepl("dpi|time|day|condition|group|injury", norm_name(cc))) score <- score + 3
    }
    score
  })

  if (all(!is.finite(scores))) return(NA_character_)
  cols[which.max(scores)]
}

sample_col <- choose_best_col(cell_meta, possible_sample_cols, "sample")
time_col   <- choose_best_col(cell_meta, possible_time_cols, "time")

if (!is.na(sample_col)) {
  cell_meta$sample_auto <- as.character(cell_meta[[sample_col]])
} else {
  cell_meta$sample_auto <- "UNRESOLVED"
}

if (!is.na(time_col)) {
  cell_meta$time_auto <- as.character(cell_meta[[time_col]])
} else {
  cell_meta$time_auto <- NA_character_
}

# 如果 time 列未找到，尝试从 sample / metadata 字符串中抽取 0/1/3/7 dpi
parse_timepoint <- function(x) {
  x0 <- tolower(as.character(x))
  out <- rep(NA_character_, length(x0))

  out[grepl("uninj|naive|sham|control|ctrl|0[ _-]*dpi|day[ _-]*0|0d", x0)] <- "Uninjured"
  out[grepl("(^|[^0-9])1[ _-]*(dpi|d|day)|day[ _-]*1([^0-9]|$)", x0)] <- "1 dpi"
  out[grepl("(^|[^0-9])3[ _-]*(dpi|d|day)|day[ _-]*3([^0-9]|$)", x0)] <- "3 dpi"
  out[grepl("(^|[^0-9])7[ _-]*(dpi|d|day)|day[ _-]*7([^0-9]|$)", x0)] <- "7 dpi"
  out
}

if (all(is.na(cell_meta$time_auto)) || all(cell_meta$time_auto == "")) {
  scan_cols <- unique(c(sample_col, possible_time_cols, possible_sample_cols))
  scan_cols <- scan_cols[!is.na(scan_cols)]

  parsed <- rep(NA_character_, nrow(cell_meta))
  for (cc in scan_cols) {
    tmp <- parse_timepoint(cell_meta[[cc]])
    parsed[is.na(parsed) & !is.na(tmp)] <- tmp[is.na(parsed) & !is.na(tmp)]
  }
  cell_meta$time_auto <- parsed
}

# 标准化常见时间标签
normalize_timepoint <- function(x) {
  z <- tolower(trimws(as.character(x)))
  out <- as.character(x)

  out[grepl("uninj|naive|sham|control|ctrl|^0$|0.*dpi|day.*0", z)] <- "Uninjured"
  out[grepl("(^|[^0-9])1([^0-9]|$).*dpi|1.*day|day.*1|^1$", z)] <- "1 dpi"
  out[grepl("(^|[^0-9])3([^0-9]|$).*dpi|3.*day|day.*3|^3$", z)] <- "3 dpi"
  out[grepl("(^|[^0-9])7([^0-9]|$).*dpi|7.*day|day.*7|^7$", z)] <- "7 dpi"
  out
}
cell_meta$time_auto <- normalize_timepoint(cell_meta$time_auto)

# sample-level manifest
sample_manifest <- unique(
  cell_meta[, c("sample_auto", "time_auto"), drop = FALSE]
)
colnames(sample_manifest) <- c("sample_id", "timepoint")

# 每个 sample 的 cell 数
cell_counts <- as.data.frame(table(cell_meta$sample_auto), stringsAsFactors = FALSE)
colnames(cell_counts) <- c("sample_id", "n_cells")

sample_manifest <- merge(sample_manifest, cell_counts, by = "sample_id", all.x = TRUE)
sample_manifest$dataset <- "GSE162610"
sample_manifest$modality <- "scRNA-seq"
sample_manifest$role <- "Discovery"
sample_manifest$species <- "mouse"

# 排序：Uninjured → 1 → 3 → 7
time_levels <- c("Uninjured", "1 dpi", "3 dpi", "7 dpi")
sample_manifest$time_order <- match(sample_manifest$timepoint, time_levels)
sample_manifest <- sample_manifest[
  order(sample_manifest$time_order, sample_manifest$sample_id),
]
sample_manifest$time_order <- NULL

fwrite(sample_manifest, file.path(admin_dir, "10_GSE162610_sample_manifest_auto.csv"))
fwrite(cell_meta, file.path(admin_dir, "11_GSE162610_cell_metadata_with_auto_fields.csv"))

# -----------------------------
# 12. 各时间点细胞数
# -----------------------------
time_counts <- as.data.frame(table(cell_meta$time_auto, useNA = "ifany"), stringsAsFactors = FALSE)
colnames(time_counts) <- c("timepoint", "n_cells")
fwrite(time_counts, file.path(admin_dir, "12_cell_counts_by_timepoint.csv"))

# -----------------------------
# 13. 自动检查主时间点是否齐全
# -----------------------------
detected_times <- unique(na.omit(cell_meta$time_auto))
expected_times <- c("Uninjured", "1 dpi", "3 dpi", "7 dpi")
timepoint_check <- data.frame(
  expected_timepoint = expected_times,
  detected = expected_times %in% detected_times
)
fwrite(timepoint_check, file.path(admin_dir, "13_expected_timepoint_check.csv"))

# -----------------------------
# 14. Step 0 完成判定
# -----------------------------
# barcode metadata 不要求一定同 barcodes 行数，只要最终能对齐即可
cell_meta_match_ok <- nrow(cell_meta) == nrow(barcodes)

step0_pass <- all(c(
  genes_match,
  barcodes_match,
  cell_meta_match_ok,
  all(timepoint_check$detected)
))

# -----------------------------
# 15. 输出审计报告
# -----------------------------
report <- c(
  "============================================================",
  "GSE162610 Step 0A Audit Report",
  "============================================================",
  paste0("Date: ", Sys.time()),
  "",
  paste0("Matrix genes: ", mtx_dims["n_rows"]),
  paste0("Matrix cells: ", mtx_dims["n_cols"]),
  paste0("Matrix non-zero entries: ", format(mtx_dims["n_nonzero"], scientific = FALSE)),
  paste0("Genes file rows: ", nrow(genes)),
  paste0("Barcode file rows: ", nrow(barcodes)),
  paste0("Barcode metadata rows: ", nrow(barcode_meta)),
  paste0("Gene metadata rows: ", nrow(gene_meta)),
  "",
  paste0("Detected barcode metadata column: ", barcode_col),
  paste0("Barcode overlap / matched rows: ", barcode_overlap),
  paste0("Auto-selected sample column: ", ifelse(is.na(sample_col), "NOT FOUND", sample_col)),
  paste0("Auto-selected time/group column: ", ifelse(is.na(time_col), "NOT FOUND", time_col)),
  "",
  "Detected timepoints:",
  paste0("  - ", detected_times),
  "",
  "Expected timepoints:",
  paste0("  - ", expected_times, ": ", ifelse(timepoint_check$detected, "PASS", "MISSING")),
  "",
  "Basic dimension checks:",
  paste0("  - ", basic_check$check, ": ", ifelse(basic_check$pass, "PASS", "CHECK")),
  "",
  paste0("STEP 0A PASS = ", step0_pass),
  "",
  "Next step:",
  if (step0_pass) {
    "Proceed to Step 1A: create Seurat object, sample-level QC, and first-pass whole-cell atlas."
  } else {
    paste0(
      "Do NOT proceed yet. Please send back: ",
      "10_GSE162610_sample_manifest_auto.csv, ",
      "09_detected_metadata_fields.txt, ",
      "04_barcode_metadata_head30.csv, ",
      "14_step0A_report.txt"
    )
  },
  "============================================================"
)

writeLines(report, file.path(admin_dir, "14_step0A_report.txt"))

cat(paste(report, collapse = "\n"))
cat("\n\n输出文件夹：", admin_dir, "\n")

if (!step0_pass) {
  cat("\n⚠ Step 0A 尚未完全通过。请不要急着进入 Seurat QC。\n")
} else {
  cat("\n✅ Step 0A 已通过，可以进入 Step 1A。\n")
}
