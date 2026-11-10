options(stringsAsFactors = FALSE)
options(timeout = 600)
project_dir <- Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
in_rds <- file.path( project_dir, "02_objects", "step1A_GSE162610_QC_audit", "GSE162610_raw_with_QCmetrics.rds" )
out_dir <- file.path( project_dir, "02_objects", "step1B_GSE162610_QC_filtered" )
fig_dir <- file.path(out_dir, "figures")
tab_dir <- file.path(out_dir, "tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(in_rds)) {
  stop("Input object is not a Seurat object.", in_rds)
}
pkgs <- c("Seurat", "data.table", "ggplot2", "patchwork")
miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss) > 0) {
  stop("Missing R packages: ", paste(miss, collapse = ", "))
}
suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(ggplot2)
  library(patchwork)
})
pbmc <- readRDS(in_rds)
required_meta <- c( "sample_id", "timepoint", "nFeature_RNA", "nCount_RNA", "percent.mt", "percent.HB" )
if (!all(required_meta %in% colnames(pbmc@meta.data))) {
  stop( "Missing metadata: ", paste(setdiff(required_meta, colnames(pbmc@meta.data)), collapse = ", ") )
}
minFeature <- 500
minCount   <- 1000
maxMT      <- 15
maxHB      <- 1
thresholds <- data.frame( metric = c( "nFeature_RNA", "nCount_RNA", "percent.mt", "percent.HB" ), rule = c( paste0(">=", minFeature), paste0(">=", minCount), paste0("<", maxMT), paste0("<", maxHB) ), rationale = c( "QC filter.", "QC filter.", "QC filter.", "QC filter." ), stringsAsFactors = FALSE )
fwrite( thresholds, file.path(tab_dir, "01_QC_thresholds.csv") )
md <- pbmc@meta.data
md$cell <- rownames(md)
md$fail_low_feature <- md$nFeature_RNA < minFeature
md$fail_low_count   <- md$nCount_RNA < minCount
md$fail_high_mt     <- md$percent.mt >= maxMT
md$fail_high_HB     <- md$percent.HB >= maxHB
md$QC_pass <- !( md$fail_low_feature | md$fail_low_count | md$fail_high_mt | md$fail_high_HB )
pbmc$QC_pass <- md$QC_pass
before_sample <- as.data.frame( table(md$sample_id), stringsAsFactors = FALSE )
colnames(before_sample) <- c("sample_id", "n_before")
after_sample <- as.data.frame( table(md$sample_id[md$QC_pass]), stringsAsFactors = FALSE )
colnames(after_sample) <- c("sample_id", "n_after")
retention_sample <- merge( before_sample, after_sample, by = "sample_id", all.x = TRUE )
retention_sample$n_after[is.na(retention_sample$n_after)] <- 0
retention_sample$n_removed <- ( retention_sample$n_before - retention_sample$n_after )
retention_sample$retention_pct <- round( 100 * retention_sample$n_after / retention_sample$n_before, 2 )
sample_time <- unique( md[, c("sample_id", "timepoint"), drop = FALSE] )
sample_time$timepoint <- as.character(sample_time$timepoint)
retention_sample <- merge( retention_sample, sample_time, by = "sample_id", all.x = TRUE )
time_levels <- c("Uninjured", "1 dpi", "3 dpi", "7 dpi")
retention_sample$time_order <- match( retention_sample$timepoint, time_levels )
retention_sample <- retention_sample[ order(retention_sample$time_order, retention_sample$sample_id), ]
retention_sample$time_order <- NULL
fwrite( retention_sample, file.path(tab_dir, "02_QC_retention_by_sample.csv") )
before_time <- as.data.frame( table(as.character(md$timepoint)), stringsAsFactors = FALSE )
colnames(before_time) <- c("timepoint", "n_before")
after_time <- as.data.frame( table(as.character(md$timepoint[md$QC_pass])), stringsAsFactors = FALSE )
colnames(after_time) <- c("timepoint", "n_after")
retention_time <- merge( before_time, after_time, by = "timepoint", all.x = TRUE )
retention_time$n_after[is.na(retention_time$n_after)] <- 0
retention_time$n_removed <- ( retention_time$n_before - retention_time$n_after )
retention_time$retention_pct <- round( 100 * retention_time$n_after / retention_time$n_before, 2 )
retention_time$time_order <- match( retention_time$timepoint, time_levels )
retention_time <- retention_time[ order(retention_time$time_order), ]
retention_time$time_order <- NULL
fwrite( retention_time, file.path(tab_dir, "03_QC_retention_by_timepoint.csv") )
fail_summary <- do.call( rbind, lapply(split(md, md$sample_id), function(x) {
    data.frame( sample_id = unique(as.character(x$sample_id))[1], timepoint = unique(as.character(x$timepoint))[1], n_cells = nrow(x), fail_low_feature = sum(x$fail_low_feature), fail_low_count = sum(x$fail_low_count), fail_high_mt = sum(x$fail_high_mt), fail_high_HB = sum(x$fail_high_HB), fail_any = sum(!x$QC_pass), pass = sum(x$QC_pass), stringsAsFactors = FALSE )
  })
)
rownames(fail_summary) <- NULL
fail_summary$time_order <- match( fail_summary$timepoint, time_levels )
fail_summary <- fail_summary[ order(fail_summary$time_order, fail_summary$sample_id), ]
fail_summary$time_order <- NULL
fwrite( fail_summary, file.path(tab_dir, "04_QC_failure_reasons_by_sample.csv") )
removed_cells <- md[ !md$QC_pass, c( "cell", "sample_id", "timepoint", "nFeature_RNA", "nCount_RNA", "percent.mt", "percent.HB", "fail_low_feature", "fail_low_count", "fail_high_mt", "fail_high_HB" ), drop = FALSE ]
fwrite( removed_cells, file.path(tab_dir, "05_removed_cells_QC.csv") )
pbmc_qc <- subset( pbmc, cells = rownames(md)[md$QC_pass] )
pbmc_qc$QC_pass <- NULL
md_after <- pbmc_qc@meta.data
qc_after_summary <- do.call( rbind, lapply(split(md_after, md_after$sample_id), function(x) {
    data.frame( sample_id = unique(as.character(x$sample_id))[1], timepoint = unique(as.character(x$timepoint))[1], n_cells = nrow(x), median_nFeature = median(x$nFeature_RNA, na.rm = TRUE), median_nCount = median(x$nCount_RNA, na.rm = TRUE), median_percent_mt = median(x$percent.mt, na.rm = TRUE), median_percent_HB = median(x$percent.HB, na.rm = TRUE), stringsAsFactors = FALSE )
  })
)
rownames(qc_after_summary) <- NULL
qc_after_summary$time_order <- match( qc_after_summary$timepoint, time_levels )
qc_after_summary <- qc_after_summary[ order(qc_after_summary$time_order, qc_after_summary$sample_id), ]
qc_after_summary$time_order <- NULL
fwrite( qc_after_summary, file.path(tab_dir, "06_QC_summary_by_sample_AFTER_filter.csv") )
features_qc <- c( "nFeature_RNA", "nCount_RNA", "percent.mt", "percent.HB" )
theme_sample <- theme( axis.title.x = element_blank(), axis.text.x = element_text( angle = 45, hjust = 1, size = 7 ), legend.position = "none" )
plots_sample <- lapply(features_qc, function(f) {
  VlnPlot( pbmc_qc, features = f, group.by = "sample_id", pt.size = 0 ) + theme_sample
})
vln_sample <- wrap_plots( plots_sample, ncol = 2 )
ggsave( file.path( fig_dir, "01_QC_violin_by_sample_AFTER_filter.pdf" ), plot = vln_sample, width = 15, height = 9 )
plots_time <- lapply(features_qc, function(f) {
  VlnPlot( pbmc_qc, features = f, group.by = "timepoint", pt.size = 0 ) + theme( axis.title.x = element_blank(), legend.position = "none" )
})
vln_time <- wrap_plots( plots_time, ncol = 2 )
ggsave( file.path( fig_dir, "02_QC_violin_by_timepoint_AFTER_filter.pdf" ), plot = vln_time, width = 12, height = 8 )
p1 <- FeatureScatter( pbmc_qc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "timepoint" )
p2 <- FeatureScatter( pbmc_qc, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "timepoint" )
scatter_after <- p1 + p2
ggsave( file.path( fig_dir, "03_QC_feature_scatter_AFTER_filter.pdf" ), plot = scatter_after, width = 14, height = 6 )
p_ret <- ggplot( retention_sample, aes( x = reorder(sample_id, retention_pct), y = retention_pct ) ) + geom_col() + coord_flip() + geom_text( aes(label = paste0(retention_pct, "%")), hjust = -0.05, size = 3 ) + ylim( 0, min( 105, max(retention_sample$retention_pct, na.rm = TRUE) + 5 ) ) + labs( x = NULL, y = "Cells retained (%)", title = "GSE162610 QC retention by sample" ) + theme_classic()
ggsave( file.path( fig_dir, "04_QC_retention_by_sample.pdf" ), plot = p_ret, width = 8, height = 6 )
saveRDS( pbmc_qc, file.path( out_dir, "GSE162610_QC_filtered.rds" ), compress = FALSE )
all_samples_preserved <- ( length(unique(pbmc_qc$sample_id)) == length(unique(pbmc$sample_id)) )
all_timepoints_preserved <- ( length(unique(as.character(pbmc_qc$timepoint))) == length(unique(as.character(pbmc$timepoint))) )
min_retention <- min( retention_sample$retention_pct, na.rm = TRUE )
retention_warning <- min_retention < 70
