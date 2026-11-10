options(stringsAsFactors = FALSE)
options(timeout = 1200)
set.seed(20260730)
project_dir <- Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
in_rds <- file.path( project_dir, "02_objects", "step1B_GSE162610_QC_filtered", "GSE162610_QC_filtered.rds" )
out_dir <- file.path( project_dir, "02_objects", "step1C_GSE162610_scDblFinder" )
fig_dir <- file.path(out_dir, "figures")
tab_dir <- file.path(out_dir, "tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(out_dir, "00_run_log.txt")
writeLines( paste0("GSE162610 Step 1C-v2 run started: ", Sys.time()), log_file )
if (!file.exists(in_rds)) {
  stop("Required input not found: ", in_rds)
}
cran_pkgs <- c("Seurat", "data.table", "ggplot2")
missing_cran <- cran_pkgs[ !vapply(cran_pkgs, requireNamespace, logical(1), quietly = TRUE) ]
if (length(missing_cran) > 0) {
  install.packages( missing_cran, repos = "https://cloud.r-project.org" )
}
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages( "BiocManager", repos = "https://cloud.r-project.org" )
}
bioc_pkgs <- c( "SingleCellExperiment", "scDblFinder", "SummarizedExperiment" )
missing_bioc <- bioc_pkgs[ !vapply(bioc_pkgs, requireNamespace, logical(1), quietly = TRUE) ]
if (length(missing_bioc) > 0) {
  BiocManager::install( missing_bioc, ask = FALSE, update = FALSE )
}
suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(ggplot2)
  library(SingleCellExperiment)
  library(scDblFinder)
  library(SummarizedExperiment)
})
pbmc <- readRDS(in_rds)
required_meta <- c( "sample_id", "timepoint", "nFeature_RNA", "nCount_RNA", "percent.mt", "percent.HB" )
if (!all(required_meta %in% colnames(pbmc@meta.data))) {
  stop( "Missing metadata: ", paste( setdiff(required_meta, colnames(pbmc@meta.data)), collapse = ", " ) )
}
sample_ids <- sort(unique(as.character(pbmc$sample_id)))
if (length(sample_ids) != 10) {
  warning("Expected-value check failed.")
}
counts <- GetAssayData( object = pbmc, assay = "RNA", layer = "counts" )
if (!identical(colnames(counts), colnames(pbmc))) {
  stop("Input object is not a Seurat object.")
}
sce <- SingleCellExperiment( assays = list(counts = counts) )
colData(sce)$sample_id <- as.character(pbmc$sample_id)
colData(sce)$timepoint <- as.character(pbmc$timepoint)
rm(counts)
gc()
t0 <- Sys.time()
sce <- scDblFinder( sce, samples = "sample_id", clusters = TRUE, dbr = NULL, verbose = TRUE )
t1 <- Sys.time()
runtime_min <- as.numeric( difftime(t1, t0, units = "mins") )
needed_cols <- c( "scDblFinder.score", "scDblFinder.class" )
if (!all(needed_cols %in% colnames(colData(sce)))) {
  stop( "Required result not found.", paste( setdiff( needed_cols, colnames(colData(sce)) ), collapse = ", " ) )
}
df_score <- colData(sce)$scDblFinder.score
df_class <- as.character( colData(sce)$scDblFinder.class )
table_class <- table(df_class)
if (!identical(colnames(sce), colnames(pbmc))) {
  stop("Input object is not a Seurat object.")
}
pbmc$Double_score <- as.numeric(df_score)
pbmc$Is_Double <- ifelse( tolower(df_class) == "doublet", "Doublet", "Singlet" )
cell_result <- data.frame( barcode = colnames(pbmc), sample_id = as.character(pbmc$sample_id), timepoint = as.character(pbmc$timepoint), nFeature_RNA = pbmc$nFeature_RNA, nCount_RNA = pbmc$nCount_RNA, percent.mt = pbmc$percent.mt, Double_score = pbmc$Double_score, Is_Double = pbmc$Is_Double, stringsAsFactors = FALSE )
fwrite( cell_result, file.path( tab_dir, "01_scDblFinder_cell_level_results.csv" ) )
summary_sample <- do.call( rbind, lapply( split(cell_result, cell_result$sample_id), function(x) {
      n_total <- nrow(x)
      n_doublet <- sum( x$Is_Double == "Doublet", na.rm = TRUE )
      n_singlet <- sum( x$Is_Double == "Singlet", na.rm = TRUE )
      data.frame( sample_id = unique(x$sample_id)[1], timepoint = unique(x$timepoint)[1], n_total = n_total, n_singlet = n_singlet, n_doublet = n_doublet, detected_doublet_pct = round( 100 * n_doublet / n_total, 2 ), median_score_singlet = median( x$Double_score[ x$Is_Double == "Singlet" ], na.rm = TRUE ), median_score_doublet = median( x$Double_score[ x$Is_Double == "Doublet" ], na.rm = TRUE ), median_nCount_singlet = median( x$nCount_RNA[ x$Is_Double == "Singlet" ], na.rm = TRUE ), median_nCount_doublet = median( x$nCount_RNA[ x$Is_Double == "Doublet" ], na.rm = TRUE ), median_nFeature_singlet = median( x$nFeature_RNA[ x$Is_Double == "Singlet" ], na.rm = TRUE ), median_nFeature_doublet = median( x$nFeature_RNA[ x$Is_Double == "Doublet" ], na.rm = TRUE ), stringsAsFactors = FALSE )
    }
  )
)
rownames(summary_sample) <- NULL
time_levels <- c( "Uninjured", "1 dpi", "3 dpi", "7 dpi" )
summary_sample$time_order <- match( summary_sample$timepoint, time_levels )
summary_sample <- summary_sample[ order( summary_sample$time_order, summary_sample$sample_id ), ]
summary_sample$time_order <- NULL
fwrite( summary_sample, file.path( tab_dir, "02_scDblFinder_summary_by_sample.csv" ) )
summary_time <- aggregate( cbind( n_total, n_singlet, n_doublet ) ~ timepoint, data = summary_sample, FUN = sum )
summary_time$detected_doublet_pct <- round( 100 * summary_time$n_doublet / summary_time$n_total, 2 )
summary_time$time_order <- match( summary_time$timepoint, time_levels )
summary_time <- summary_time[ order(summary_time$time_order), ]
summary_time$time_order <- NULL
fwrite( summary_time, file.path( tab_dir, "03_scDblFinder_summary_by_timepoint.csv" ) )
p_rate <- ggplot( summary_sample, aes( x = reorder( sample_id, detected_doublet_pct ), y = detected_doublet_pct ) ) + geom_col() + coord_flip() + geom_text( aes( label = paste0( detected_doublet_pct, "%" ) ), hjust = -0.08, size = 3 ) + expand_limits( y = max( summary_sample$detected_doublet_pct ) + 1 ) + labs( x = NULL, y = "Detected doublets (%)", title = "GSE162610 scDblFinder rate by sample" ) + theme_classic()
ggsave( file.path( fig_dir, "01_doublet_rate_by_sample.pdf" ), p_rate, width = 8, height = 6 )
plot_qc <- data.frame( Is_Double = pbmc$Is_Double, nCount_RNA = pbmc$nCount_RNA, nFeature_RNA = pbmc$nFeature_RNA, stringsAsFactors = FALSE )
p_count <- ggplot( plot_qc, aes( x = Is_Double, y = nCount_RNA ) ) + geom_violin( scale = "width", trim = TRUE ) + labs( x = NULL, y = "nCount_RNA", title = "UMI counts" ) + theme_classic()
p_feature <- ggplot( plot_qc, aes( x = Is_Double, y = nFeature_RNA ) ) + geom_violin( scale = "width", trim = TRUE ) + labs( x = NULL, y = "nFeature_RNA", title = "Detected genes" ) + theme_classic()
if (requireNamespace("patchwork", quietly = TRUE)) {
  p_qc <- p_count + p_feature
  ggsave( file.path( fig_dir, "02_doublet_vs_singlet_QC_violin.pdf" ), p_qc, width = 10, height = 5 )
} else {
  ggsave( file.path( fig_dir, "02A_doublet_vs_singlet_nCount.pdf" ), p_count, width = 5, height = 5 )
  ggsave( file.path( fig_dir, "02B_doublet_vs_singlet_nFeature.pdf" ), p_feature, width = 5, height = 5 )
}
p_score <- ggplot( cell_result, aes( x = Double_score, fill = Is_Double ) ) + geom_histogram( bins = 70, position = "identity", alpha = 0.55 ) + facet_wrap( ~ sample_id, scales = "free_y" ) + labs( x = "scDblFinder score", y = "Cell count", title = "scDblFinder score distribution by sample", fill = NULL ) + theme_classic()
ggsave( file.path( fig_dir, "03_doublet_score_distribution_by_sample.pdf" ), p_score, width = 14, height = 9 )
saveRDS( pbmc, file.path( out_dir, "GSE162610_QC_with_scDblFinder_labels.rds" ), compress = FALSE )
pbmc_singlet <- subset( pbmc, subset = Is_Double == "Singlet" )
saveRDS( pbmc_singlet, file.path( out_dir, "GSE162610_QC_singlets.rds" ), compress = FALSE )
all_samples_preserved <- ( length(unique(pbmc_singlet$sample_id)) == length(unique(pbmc$sample_id)) )
all_timepoints_preserved <- ( length( unique(as.character(pbmc_singlet$timepoint)) ) == length( unique(as.character(pbmc$timepoint)) ) )
overall_doublet_pct <- round( 100 * (ncol(pbmc) - ncol(pbmc_singlet)) / ncol(pbmc), 2 )
rate_warning <- any( summary_sample$detected_doublet_pct > 15 | summary_sample$detected_doublet_pct < 0.05 )
