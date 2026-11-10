options(stringsAsFactors = FALSE)
options(timeout = 600)
project_dir <- Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
raw_dir <- file.path(project_dir, "01_raw_data", "GSE162610")
meta_dir <- file.path(project_dir, "00_admin", "step0_GSE162610_confirmed")
out_dir <- file.path(project_dir, "02_objects", "step1A_GSE162610_QC_audit")
fig_dir <- file.path(out_dir, "figures")
tab_dir <- file.path(out_dir, "tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
f_matrix <- file.path(raw_dir, "GSE162610_sci_mat.mtx.gz")
f_barcodes <- file.path(raw_dir, "GSE162610_barcodes.tsv.gz")
f_genes <- file.path(raw_dir, "GSE162610_genes.tsv.gz")
f_cellmeta <- file.path(meta_dir, "01_GSE162610_cell_metadata_CONFIRMED.csv")
f_manifest <- file.path(meta_dir, "02_GSE162610_sample_manifest_CONFIRMED.csv")
required_files <- c(f_matrix, f_barcodes, f_genes, f_cellmeta, f_manifest)
if (!all(file.exists(required_files))) {
  stop( "Required input not found: ", paste(required_files[!file.exists(required_files)], collapse = "\n") )
}
pkgs <- c("Matrix", "Seurat", "data.table", "ggplot2", "patchwork")
miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss) > 0) {
  stop("Missing R packages: ", paste(miss, collapse = ", "))
}
suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
  library(data.table)
  library(ggplot2)
  library(patchwork)
})
read_gz_tab <- function(path) {
  con <- gzfile(path, "rt")
  on.exit(close(con))
  read.delim( con, header = FALSE, sep = "\t", quote = "", comment.char = "", check.names = FALSE, stringsAsFactors = FALSE )
}
clean_chr <- function(x) {
  x <- trimws(as.character(x))
  gsub('^"|"$', '', x)
}
barcodes_df <- read_gz_tab(f_barcodes)
genes_df <- read_gz_tab(f_genes)
barcodes <- clean_chr(barcodes_df[[1]])
if (ncol(genes_df) >= 2) {
  gene_names <- clean_chr(genes_df[[2]])
  if (mean(nzchar(gene_names)) < 0.90) {
    gene_names <- clean_chr(genes_df[[1]])
  }
} else {
  gene_names <- clean_chr(genes_df[[1]])
}
if (mean(grepl("^ENSMUSG", gene_names, ignore.case = TRUE)) > 0.50 && ncol(genes_df) >= 2) {
  alt <- clean_chr(genes_df[[1]])
  if (mean(grepl("^ENSMUSG", alt, ignore.case = TRUE)) < 0.50) {
    gene_names <- alt
  }
}
gene_names <- make.unique(gene_names)
t0 <- Sys.time()
con <- gzfile(f_matrix, open = "rt")
counts <- Matrix::readMM(con)
close(con)
counts <- as(counts, "CsparseMatrix")
rownames(counts) <- gene_names
colnames(counts) <- barcodes
t1 <- Sys.time()
load_minutes <- as.numeric(difftime(t1, t0, units = "mins"))
if (nrow(counts) != length(gene_names)) {
  stop("Consistency check failed.")
}
if (ncol(counts) != length(barcodes)) {
  stop("Consistency check failed.")
}
pbmc <- CreateSeuratObject( counts = counts, project = "GSE162610", min.cells = 0, min.features = 0 )
rm(counts)
gc()
cell_meta <- data.table::fread(f_cellmeta, data.table = FALSE)
if (!all(c("barcode", "sample_id", "timepoint") %in% colnames(cell_meta))) {
  stop( "Missing metadata: ", "Check", f_cellmeta )
}
cell_meta$barcode <- clean_chr(cell_meta$barcode)
cell_meta$sample_id <- clean_chr(cell_meta$sample_id)
cell_meta$timepoint <- clean_chr(cell_meta$timepoint)
idx <- match(colnames(pbmc), cell_meta$barcode)
if (anyNA(idx)) {
  stop("Validation check.", sum(is.na(idx)), "Validation failed.")
}
cell_meta <- cell_meta[idx, , drop = FALSE]
rownames(cell_meta) <- cell_meta$barcode
core_meta <- cell_meta[, c("sample_id", "timepoint", "dataset", "modality", "role"), drop = FALSE]
pbmc <- AddMetaData(pbmc, metadata = core_meta)
pbmc$orig.ident <- pbmc$sample_id
pbmc$timepoint <- factor( pbmc$timepoint, levels = c("Uninjured", "1 dpi", "3 dpi", "7 dpi"), ordered = TRUE )
mt_genes <- grep("^mt-", rownames(pbmc), value = TRUE, ignore.case = TRUE)
if (length(mt_genes) > 0) {
  pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, features = mt_genes)
} else {
  warning("Validation check.")
  pbmc$percent.mt <- 0
}
hb_genes <- grep("^Hb[ab]", rownames(pbmc), value = TRUE, ignore.case = TRUE)
if (length(hb_genes) > 0) {
  pbmc[["percent.HB"]] <- PercentageFeatureSet(pbmc, features = hb_genes)
} else {
  warning("Validation check.")
  pbmc$percent.HB <- 0
}
writeLines(mt_genes, file.path(tab_dir, "01_detected_mitochondrial_genes.txt"))
writeLines(hb_genes, file.path(tab_dir, "02_detected_hemoglobin_genes.txt"))
sample_counts <- as.data.frame(table(pbmc$sample_id), stringsAsFactors = FALSE)
colnames(sample_counts) <- c("sample_id", "n_cells")
sample_counts <- sample_counts[order(sample_counts$sample_id), ]
time_counts <- as.data.frame(table(pbmc$timepoint), stringsAsFactors = FALSE)
colnames(time_counts) <- c("timepoint", "n_cells")
fwrite(sample_counts, file.path(tab_dir, "03_cell_counts_by_sample_before_QC.csv"))
fwrite(time_counts, file.path(tab_dir, "04_cell_counts_by_time_before_QC.csv"))
md <- pbmc@meta.data
md$cell <- rownames(md)
qc_summary <- do.call( rbind, lapply(split(md, md$sample_id), function(x) {
    data.frame( sample_id = unique(as.character(x$sample_id))[1], timepoint = unique(as.character(x$timepoint))[1], n_cells = nrow(x), median_nFeature = median(x$nFeature_RNA, na.rm = TRUE), median_nCount = median(x$nCount_RNA, na.rm = TRUE), median_percent_mt = median(x$percent.mt, na.rm = TRUE), median_percent_HB = median(x$percent.HB, na.rm = TRUE), mean_percent_mt = mean(x$percent.mt, na.rm = TRUE), mean_percent_HB = mean(x$percent.HB, na.rm = TRUE) )
  })
)
rownames(qc_summary) <- NULL
fwrite(qc_summary, file.path(tab_dir, "05_QC_summary_by_sample.csv"))
quant_probs <- c(0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99)
qc_quantiles <- do.call( rbind, lapply(split(md, md$sample_id), function(x) {
    metrics <- c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.HB")
    do.call( rbind, lapply(metrics, function(m) {
        qs <- quantile(x[[m]], probs = quant_probs, na.rm = TRUE, names = FALSE)
        data.frame( sample_id = unique(as.character(x$sample_id))[1], timepoint = unique(as.character(x$timepoint))[1], metric = m, quantile = paste0(quant_probs * 100, "%"), value = as.numeric(qs) )
      })
    )
  })
)
rownames(qc_quantiles) <- NULL
fwrite(qc_quantiles, file.path(tab_dir, "06_QC_quantiles_by_sample.csv"))
global_qc <- data.frame( metric = c( "n_cells", "n_genes", "median_nFeature_RNA", "median_nCount_RNA", "median_percent_mt", "median_percent_HB" ), value = c( ncol(pbmc), nrow(pbmc), median(pbmc$nFeature_RNA), median(pbmc$nCount_RNA), median(pbmc$percent.mt), median(pbmc$percent.HB) ) )
fwrite(global_qc, file.path(tab_dir, "07_global_QC_summary.csv"))
theme_qc <- theme( axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, size = 7), legend.position = "none" )
features_qc <- c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.HB")
plots_sample <- lapply(features_qc, function(f) {
  VlnPlot( pbmc, features = f, group.by = "sample_id", pt.size = 0 ) + theme_qc
})
vln_sample <- wrap_plots(plots_sample, ncol = 2)
ggsave( file.path(fig_dir, "01_QC_violin_by_sample_BEFORE_filter.pdf"), plot = vln_sample, width = 15, height = 9 )
plots_time <- lapply(features_qc, function(f) {
  VlnPlot( pbmc, features = f, group.by = "timepoint", pt.size = 0 ) + theme( axis.title.x = element_blank(), legend.position = "none" )
})
vln_time <- wrap_plots(plots_time, ncol = 2)
ggsave( file.path(fig_dir, "02_QC_violin_by_timepoint_BEFORE_filter.pdf"), plot = vln_time, width = 12, height = 8 )
p_scatter1 <- FeatureScatter( pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "timepoint" )
p_scatter2 <- FeatureScatter( pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "timepoint" )
scatter_all <- p_scatter1 + p_scatter2
ggsave( file.path(fig_dir, "03_QC_feature_scatter_BEFORE_filter.pdf"), plot = scatter_all, width = 14, height = 6 )
saveRDS( pbmc, file.path(out_dir, "GSE162610_raw_with_QCmetrics.rds"), compress = FALSE )
manifest <- fread(f_manifest, data.table = FALSE)
sample_counts_ok <- ( ncol(pbmc) == sum(manifest$n_cells) && ncol(pbmc) == 66178 )
