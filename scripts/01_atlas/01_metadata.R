options(stringsAsFactors = FALSE)
project_dir <- Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
raw_dir <- file.path(project_dir, "01_raw_data", "GSE162610")
out_dir <- file.path(project_dir, "00_admin", "step0_GSE162610_confirmed")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
if (!requireNamespace("data.table", quietly = TRUE)) {
  install.packages("data.table", repos = "https://cloud.r-project.org")
}
suppressPackageStartupMessages(library(data.table))
f_barcodes <- file.path(raw_dir, "GSE162610_barcodes.tsv.gz")
f_bmeta    <- file.path(raw_dir, "GSE162610_barcode_metadata.tsv.gz")
f_genes    <- file.path(raw_dir, "GSE162610_genes.tsv.gz")
f_gmeta    <- file.path(raw_dir, "GSE162610_gene_metadata.tsv.gz")
f_matrix   <- file.path(raw_dir, "GSE162610_sci_mat.mtx.gz")
stopifnot(all(file.exists(c(f_barcodes, f_bmeta, f_genes, f_gmeta, f_matrix))))
read_tab <- function(path, header = FALSE) {
  con <- gzfile(path, "rt")
  on.exit(close(con))
  read.delim( con, header = header, sep = "\t", quote = "", comment.char = "", check.names = FALSE, stringsAsFactors = FALSE )
}
clean_chr <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub('^"|"$', '', x)
  x
}
read_mtx_dims <- function(path) {
  con <- gzfile(path, "rt")
  on.exit(close(con))
  first <- readLines(con, n = 1)
  if (!grepl("^%%MatrixMarket", first)) stop("Matrix file is not MatrixMarket.")
  repeat {
    z <- readLines(con, n = 1)
    if (length(z) == 0) stop("Matrix dimensions not found.")
    if (!grepl("^%", z)) break
  }
  v <- strsplit(trimws(z), "\\s+")[[1]]
  c(n_genes = as.integer(v[1]), n_cells = as.integer(v[2]), n_nonzero = as.double(v[3]))
}
cat("============================================================\n")
cat("Step 0C | GSE162610 metadata 自适应解析\n")
cat("============================================================\n\n")
barcodes_df <- read_tab(f_barcodes, header = FALSE)
genes_df    <- read_tab(f_genes, header = FALSE)
bmeta_raw   <- read_tab(f_bmeta, header = FALSE)
gmeta_raw   <- read_tab(f_gmeta, header = FALSE)
dims        <- read_mtx_dims(f_matrix)
barcodes <- clean_chr(barcodes_df[[1]])
barcodes <- barcodes[nzchar(barcodes)]
cat("Matrix:", dims["n_genes"], "genes x", dims["n_cells"], "cells\n")
cat("barcodes.tsv:", length(barcodes), "\n")
cat("barcode metadata raw rows:", nrow(bmeta_raw), "\n")
cat("gene metadata raw rows:", nrow(gmeta_raw), "\n\n")
for (j in seq_len(ncol(bmeta_raw))) {
  bmeta_raw[[j]] <- clean_chr(bmeta_raw[[j]])
}
for (j in seq_len(ncol(gmeta_raw))) {
  gmeta_raw[[j]] <- clean_chr(gmeta_raw[[j]])
}
barcode_overlap <- sapply(seq_len(ncol(bmeta_raw)), function(j) {
  sum(bmeta_raw[[j]] %in% barcodes)
})
barcode_col_idx <- which.max(barcode_overlap)
barcode_overlap_max <- barcode_overlap[barcode_col_idx]
if (barcode_overlap_max < 0.95 * length(barcodes)) {
  stop( "无法可靠识别 barcode 列。最大重叠仅 ", barcode_overlap_max, " / ", length(barcodes) )
}
real_row <- bmeta_raw[[barcode_col_idx]] %in% barcodes
bmeta <- bmeta_raw[real_row, , drop = FALSE]
cat("自动识别 barcode 列: V", barcode_col_idx, "\n", sep = "")
cat("真实 barcode 行:", nrow(bmeta), "\n")
cat("自动删除 header/非细胞行:", nrow(bmeta_raw) - nrow(bmeta), "\n\n")
idx <- match(barcodes, bmeta[[barcode_col_idx]])
if (anyNA(idx)) {
  stop("metadata 中仍有 barcode 缺失，缺失数 = ", sum(is.na(idx)))
}
bmeta <- bmeta[idx, , drop = FALSE]
barcode_order_ok <- identical( barcodes, as.character(bmeta[[barcode_col_idx]]) )
sample_regex <- '^(uninj|1dpi|3dpi|7dpi)_sample[0-9]+$'
sample_scores <- sapply(seq_len(ncol(bmeta)), function(j) {
  x <- tolower(clean_chr(bmeta[[j]]))
  mean(grepl(sample_regex, x))
})
sample_col_idx <- which.max(sample_scores)
if (sample_scores[sample_col_idx] < 0.90) {
  stop( "无法可靠识别 sample 列。最佳列匹配比例 = ", round(sample_scores[sample_col_idx], 4) )
}
sample_id <- clean_chr(bmeta[[sample_col_idx]])
normalize_time <- function(x) {
  z <- tolower(clean_chr(x))
  out <- rep(NA_character_, length(z))
  out[z %in% c("uninjured","uninj","naive","control","ctrl","0dpi","0 dpi")] <- "Uninjured"
  out[z %in% c("1dpi","1 dpi","day1","day 1","1d")] <- "1 dpi"
  out[z %in% c("3dpi","3 dpi","day3","day 3","3d")] <- "3 dpi"
  out[z %in% c("7dpi","7 dpi","day7","day 7","7d")] <- "7 dpi"
  out
}
time_scores <- sapply(seq_len(ncol(bmeta)), function(j) {
  mean(!is.na(normalize_time(bmeta[[j]])))
})
time_col_idx <- which.max(time_scores)
timepoint_from_meta <- normalize_time(bmeta[[time_col_idx]])
if (time_scores[time_col_idx] >= 0.90) {
  timepoint <- timepoint_from_meta
  time_source <- paste0("metadata V", time_col_idx)
} else {
  sid <- tolower(sample_id)
  timepoint <- rep(NA_character_, length(sid))
  timepoint[grepl("^uninj_", sid)] <- "Uninjured"
  timepoint[grepl("^1dpi_", sid)]  <- "1 dpi"
  timepoint[grepl("^3dpi_", sid)]  <- "3 dpi"
  timepoint[grepl("^7dpi_", sid)]  <- "7 dpi"
  time_source <- "derived from sample_id"
}
cell_metadata <- data.frame( barcode = barcodes, sample_id = sample_id, timepoint = timepoint, dataset = "GSE162610", modality = "scRNA-seq", role = "Discovery", stringsAsFactors = FALSE )
orig <- as.data.frame(bmeta, stringsAsFactors = FALSE)
colnames(orig) <- paste0("original_V", seq_len(ncol(orig)))
cell_metadata <- cbind(cell_metadata, orig)
manifest <- aggregate( barcode ~ sample_id + timepoint, data = cell_metadata, FUN = length )
colnames(manifest)[3] <- "n_cells"
time_levels <- c("Uninjured","1 dpi","3 dpi","7 dpi")
manifest$time_order <- match(manifest$timepoint, time_levels)
manifest <- manifest[order(manifest$time_order, manifest$sample_id), ]
manifest$time_order <- NULL
manifest$dataset <- "GSE162610"
manifest$role <- "Discovery"
expected_samples <- c( "uninj_sample1","uninj_sample2","uninj_sample3", "1dpi_sample1","1dpi_sample2","1dpi_sample3", "3dpi_sample1","3dpi_sample2", "7dpi_sample1","7dpi_sample2" )
expected_times <- c("Uninjured","1 dpi","3 dpi","7 dpi")
sample_check <- setequal(sort(unique(manifest$sample_id)), sort(expected_samples))
time_check <- setequal(sort(unique(manifest$timepoint)), sort(expected_times))
time_counts <- aggregate( n_cells ~ timepoint, data = manifest, FUN = sum )
time_counts$time_order <- match(time_counts$timepoint, time_levels)
time_counts <- time_counts[order(time_counts$time_order), c("timepoint","n_cells")]
if (nrow(gmeta_raw) == as.integer(dims["n_genes"]) + 1) {
  gmeta_data <- gmeta_raw[-1, , drop = FALSE]
  gene_header_removed <- TRUE
} else if (nrow(gmeta_raw) == as.integer(dims["n_genes"])) {
  gmeta_data <- gmeta_raw
  gene_header_removed <- FALSE
} else {
  gmeta_data <- gmeta_raw
  gene_header_removed <- NA
}
gene_meta_count_ok <- nrow(gmeta_data) == as.integer(dims["n_genes"])
pass <- all(c( as.integer(dims["n_cells"]) == length(barcodes), nrow(cell_metadata) == as.integer(dims["n_cells"]), barcode_overlap_max == length(barcodes), barcode_order_ok, sum(manifest$n_cells) == as.integer(dims["n_cells"]), sample_check, time_check, all(!is.na(cell_metadata$timepoint)), as.integer(dims["n_genes"]) == nrow(genes_df), gene_meta_count_ok ))
fwrite(cell_metadata, file.path(out_dir, "01_GSE162610_cell_metadata_CONFIRMED.csv"))
fwrite(manifest, file.path(out_dir, "02_GSE162610_sample_manifest_CONFIRMED.csv"))
fwrite(time_counts, file.path(out_dir, "03_GSE162610_timepoint_counts.csv"))
detect_info <- data.frame( item = c( "barcode_col", "barcode_overlap", "sample_col", "sample_match_fraction", "time_col_candidate", "time_match_fraction", "time_source", "gene_header_removed" ), value = c( paste0("V", barcode_col_idx), barcode_overlap_max, paste0("V", sample_col_idx), round(sample_scores[sample_col_idx], 6), paste0("V", time_col_idx), round(time_scores[time_col_idx], 6), time_source, gene_header_removed ) )
fwrite(detect_info, file.path(out_dir, "04_auto_detection_info.csv"))
report <- c( "============================================================", "GSE162610 STEP 0C CONFIRMED REPORT", "============================================================", paste0("Matrix: ", dims["n_genes"], " genes x ", dims["n_cells"], " cells"), paste0("barcodes.tsv: ", length(barcodes)), paste0("barcode metadata raw rows: ", nrow(bmeta_raw)), paste0("real cell metadata rows after barcode filtering: ", nrow(cell_metadata)), paste0("barcode column: V", barcode_col_idx), paste0("barcode exact overlap: ", barcode_overlap_max, " / ", length(barcodes)), paste0("barcode order identical after reordering: ", barcode_order_ok), paste0("sample column: V", sample_col_idx), paste0("time source: ", time_source), "", "Samples:", paste0("  ", manifest$sample_id, " | ", manifest$timepoint, " | n=", manifest$n_cells), "", "Timepoint totals:", paste0("  ", time_counts$timepoint, " | n=", time_counts$n_cells), "", paste0("All expected 10 samples detected exactly: ", sample_check), paste0("All expected 4 timepoints detected exactly: ", time_check), paste0("All cells have timepoint: ", all(!is.na(cell_metadata$timepoint))), paste0("Gene metadata corrected rows: ", nrow(gmeta_data)), paste0("Gene metadata count matches matrix genes: ", gene_meta_count_ok), "", paste0("STEP 0 FINAL PASS = ", pass), "", if (pass) {
    "NEXT: Step 1A | Build Seurat object and perform sample-level QC."
  } else {
    "STOP: Send this report and 04_auto_detection_info.csv for review."
  },
  "============================================================"
)
writeLines(report, file.path(out_dir, "05_STEP0_CONFIRMED_report.txt"))
cat(paste(report, collapse = "\n"))
if (pass) {
  cat("\n\nStep 0 status: PASS\n")
} else {
  cat("\n\nWarning: Step 0 仍有自动解析项未通过。\n")
}
