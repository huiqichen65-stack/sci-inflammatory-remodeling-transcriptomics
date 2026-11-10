# ============================================================
# Step 0B｜GSE162610 元数据格式修正 + 最终 sample manifest
#
# 为什么需要本脚本：
#   GSE162610_barcode_metadata.tsv.gz 和 gene_metadata.tsv.gz 实际为“无表头文件”。
#   上一步按 header=TRUE 读取，导致第一行数据被误当作列名，
#   因而少 1 个 cell 和少 1 个 gene，并把 sample/time 字段识别错位。
#
# 本脚本将：
#   1) 用 header=FALSE 重新读取 metadata
#   2) 固定 V1=barcode, V2=sample, V12=timepoint
#   3) 检查 33177 genes × 66178 cells 是否完全一致
#   4) 生成最终 sample manifest
#   5) 输出 PASS/FAIL
# ============================================================

options(stringsAsFactors = FALSE)

project_dir <- Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
raw_dir <- file.path(project_dir, "01_raw_data", "GSE162610")
out_dir <- file.path(project_dir, "00_admin", "step0_GSE162610_final")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

need <- "data.table"
if (!requireNamespace(need, quietly = TRUE)) {
  install.packages(need, repos = "https://cloud.r-project.org")
}
suppressPackageStartupMessages(library(data.table))

files <- list(
  barcode_metadata = file.path(raw_dir, "GSE162610_barcode_metadata.tsv.gz"),
  barcodes = file.path(raw_dir, "GSE162610_barcodes.tsv.gz"),
  gene_metadata = file.path(raw_dir, "GSE162610_gene_metadata.tsv.gz"),
  genes = file.path(raw_dir, "GSE162610_genes.tsv.gz"),
  matrix = file.path(raw_dir, "GSE162610_sci_mat.mtx.gz")
)

stopifnot(all(file.exists(unlist(files))))

read_gz_tab <- function(path, header = FALSE) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con))
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

read_mtx_dims <- function(path) {
  con <- gzfile(path, "rt")
  on.exit(close(con))
  x <- readLines(con, n = 1)
  if (!grepl("^%%MatrixMarket", x)) stop("Not a MatrixMarket file.")
  repeat {
    z <- readLines(con, n = 1)
    if (length(z) == 0) stop("Matrix dimension line not found.")
    if (!grepl("^%", z)) break
  }
  v <- strsplit(trimws(z), "\\s+")[[1]]
  c(n_genes = as.integer(v[1]),
    n_cells = as.integer(v[2]),
    n_nonzero = as.double(v[3]))
}

cat("============================================================\n")
cat("Step 0B | GSE162610 metadata 修正\n")
cat("============================================================\n\n")

# 1. Headerless read
barcodes <- read_gz_tab(files$barcodes, header = FALSE)
genes <- read_gz_tab(files$genes, header = FALSE)
barcode_meta <- read_gz_tab(files$barcode_metadata, header = FALSE)
gene_meta <- read_gz_tab(files$gene_metadata, header = FALSE)
dims <- read_mtx_dims(files$matrix)

colnames(barcodes)[1] <- "barcode"

# 2. Basic dimensions
checks <- data.frame(
  check = c(
    "matrix genes == genes.tsv",
    "matrix cells == barcodes.tsv",
    "barcode metadata rows == cells",
    "gene metadata rows == genes"
  ),
  observed = c(
    paste0(dims["n_genes"], " vs ", nrow(genes)),
    paste0(dims["n_cells"], " vs ", nrow(barcodes)),
    paste0(nrow(barcode_meta), " vs ", dims["n_cells"]),
    paste0(nrow(gene_meta), " vs ", dims["n_genes"])
  ),
  pass = c(
    dims["n_genes"] == nrow(genes),
    dims["n_cells"] == nrow(barcodes),
    nrow(barcode_meta) == dims["n_cells"],
    nrow(gene_meta) == dims["n_genes"]
  )
)
fwrite(checks, file.path(out_dir, "01_dimension_checks.csv"))

# 3. Known structure from actual file inspection
if (ncol(barcode_meta) < 12) {
  stop("barcode metadata 列数不足 12，无法按 GSE162610 当前格式解析。")
}

colnames(barcode_meta) <- paste0("V", seq_len(ncol(barcode_meta)))
barcode_meta$barcode <- as.character(barcode_meta$V1)
barcode_meta$sample_id <- as.character(barcode_meta$V2)
barcode_meta$timepoint_raw <- as.character(barcode_meta$V12)

normalize_time <- function(x) {
  z <- tolower(trimws(as.character(x)))
  out <- rep(NA_character_, length(z))
  out[grepl("uninj|naive|control|ctrl|sham|^0dpi$", z)] <- "Uninjured"
  out[grepl("^1dpi$|^1 dpi$|^day ?1$", z)] <- "1 dpi"
  out[grepl("^3dpi$|^3 dpi$|^day ?3$", z)] <- "3 dpi"
  out[grepl("^7dpi$|^7 dpi$|^day ?7$", z)] <- "7 dpi"
  out
}

barcode_meta$timepoint <- normalize_time(barcode_meta$timepoint_raw)

# 4. Barcode exact-order check
barcode_order_match <- identical(
  as.character(barcodes$barcode),
  as.character(barcode_meta$barcode)
)

barcode_set_match <- setequal(
  as.character(barcodes$barcode),
  as.character(barcode_meta$barcode)
)

# 5. Build final sample manifest
cell_counts <- as.data.frame(
  table(barcode_meta$sample_id),
  stringsAsFactors = FALSE
)
colnames(cell_counts) <- c("sample_id", "n_cells")

sample_time <- unique(
  barcode_meta[, c("sample_id", "timepoint_raw", "timepoint")]
)

manifest <- merge(sample_time, cell_counts, by = "sample_id", all.x = TRUE)
manifest$dataset <- "GSE162610"
manifest$modality <- "scRNA-seq"
manifest$role <- "Discovery"
manifest$species <- "mouse"

levels_time <- c("Uninjured", "1 dpi", "3 dpi", "7 dpi")
manifest$time_order <- match(manifest$timepoint, levels_time)
manifest <- manifest[order(manifest$time_order, manifest$sample_id), ]
manifest$time_order <- NULL

fwrite(manifest, file.path(out_dir, "02_GSE162610_sample_manifest_FINAL.csv"))

# 6. Sample/time consistency
sample_n_times <- aggregate(
  timepoint ~ sample_id,
  data = barcode_meta,
  FUN = function(x) length(unique(na.omit(x)))
)
colnames(sample_n_times)[2] <- "n_timepoints_per_sample"

sample_time_consistency <- all(sample_n_times$n_timepoints_per_sample == 1)

expected_samples <- c(
  "uninj_sample1","uninj_sample2","uninj_sample3",
  "1dpi_sample1","1dpi_sample2","1dpi_sample3",
  "3dpi_sample1","3dpi_sample2",
  "7dpi_sample1","7dpi_sample2"
)

expected_times <- c("Uninjured", "1 dpi", "3 dpi", "7 dpi")

sample_check <- data.frame(
  expected_sample = expected_samples,
  detected = expected_samples %in% manifest$sample_id
)
time_check <- data.frame(
  expected_timepoint = expected_times,
  detected = expected_times %in% manifest$timepoint
)

fwrite(sample_check, file.path(out_dir, "03_expected_sample_check.csv"))
fwrite(time_check, file.path(out_dir, "04_expected_timepoint_check.csv"))

# 7. Cell count by time
time_counts <- aggregate(
  barcode ~ timepoint,
  data = barcode_meta,
  FUN = length
)
colnames(time_counts) <- c("timepoint", "n_cells")
time_counts$time_order <- match(time_counts$timepoint, levels_time)
time_counts <- time_counts[order(time_counts$time_order), c("timepoint","n_cells")]
fwrite(time_counts, file.path(out_dir, "05_cell_counts_by_timepoint.csv"))

# 8. Preserve useful metadata preview
preview_cols <- unique(c(
  "barcode","sample_id","timepoint_raw","timepoint",
  if (ncol(barcode_meta) >= 30) "V30" else NULL
))
fwrite(
  head(barcode_meta[, preview_cols, drop = FALSE], 50),
  file.path(out_dir, "06_metadata_preview50.csv")
)

# If V30 looks like author broad annotation, summarize it
if (ncol(barcode_meta) >= 30) {
  anno <- as.character(barcode_meta$V30)
  if (length(unique(anno)) > 1 && length(unique(anno)) < 100) {
    anno_tab <- as.data.frame(table(anno), stringsAsFactors = FALSE)
    colnames(anno_tab) <- c("V30_value", "n_cells")
    anno_tab <- anno_tab[order(-anno_tab$n_cells), ]
    fwrite(anno_tab, file.path(out_dir, "07_V30_annotation_summary.csv"))
  }
}

# 9. Final pass
pass <- all(checks$pass) &&
        barcode_set_match &&
        all(sample_check$detected) &&
        all(time_check$detected) &&
        sample_time_consistency &&
        sum(manifest$n_cells) == dims["n_cells"]

report <- c(
  "============================================================",
  "GSE162610 Step 0B FINAL REPORT",
  "============================================================",
  paste0("Matrix: ", dims["n_genes"], " genes x ", dims["n_cells"], " cells"),
  paste0("Non-zero entries: ", format(dims["n_nonzero"], scientific = FALSE)),
  "",
  paste0("barcode metadata rows: ", nrow(barcode_meta)),
  paste0("gene metadata rows: ", nrow(gene_meta)),
  paste0("barcode order identical: ", barcode_order_match),
  paste0("barcode set identical: ", barcode_set_match),
  "",
  "Samples:",
  paste0("  ", manifest$sample_id, " | ", manifest$timepoint,
         " | n=", manifest$n_cells),
  "",
  "Timepoint totals:",
  paste0("  ", time_counts$timepoint, " | n=", time_counts$n_cells),
  "",
  paste0("All 10 expected samples detected: ", all(sample_check$detected)),
  paste0("All 4 expected timepoints detected: ", all(time_check$detected)),
  paste0("Each sample maps to one timepoint: ", sample_time_consistency),
  "",
  paste0("STEP 0 FINAL PASS = ", pass),
  "",
  if (pass) {
    "NEXT: Step 1A | Construct Seurat object and perform sample-level QC."
  } else {
    "STOP: Send this report plus 02_GSE162610_sample_manifest_FINAL.csv back for review."
  },
  "============================================================"
)

writeLines(report, file.path(out_dir, "08_step0_FINAL_report.txt"))
cat(paste(report, collapse = "\n"))

if (pass) {
  cat("\n\n✅ Step 0 最终通过。\n")
} else {
  cat("\n\n⚠ Step 0 仍未通过，请先检查输出。\n")
}
