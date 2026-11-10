options(stringsAsFactors = FALSE)
set.seed(1234)
required_pkgs <- c( "Seurat", "SeuratObject", "dplyr", "ggplot2", "patchwork", "pheatmap" )
missing_pkgs <- required_pkgs[ !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE) ]
if (length(missing_pkgs) > 0) {
  stop("Missing R packages: ", paste(missing_pkgs, collapse = ", "))
}
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(pheatmap)
})
project_dir <- Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
step2B_dir <- file.path( project_dir, "04_results", "step2B_GSE162610_lineage_adjudication" )
input_rds <- file.path( step2B_dir, "GSE162610_step2B_lineage_adjudication.rds" )
out_dir <- file.path( project_dir, "04_results", "step2C_GSE162610_final_microglia_extraction" )
fig_dir <- file.path(out_dir, "figures")
tab_dir <- file.path(out_dir, "tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(input_rds)) {
  stop( "Required input not found: ", input_rds )
}
obj <- readRDS(input_rds)
if (!inherits(obj, "Seurat")) {
  stop("Input object is not a Seurat object.")
}
DefaultAssay(obj) <- "RNA"
meta <- obj[[]]
required_meta <- c( "myeloid_cluster_r04", "sample_id", "timepoint" )
missing_meta <- required_meta[ !required_meta %in% colnames(meta) ]
if (length(missing_meta) > 0) {
  stop( "Missing metadata: ", paste(missing_meta, collapse = ", ") )
}
obj$myeloid_cluster_r04 <- as.character(obj$myeloid_cluster_r04)
print(table(obj$myeloid_cluster_r04))
clean_clusters <- c("0", "3", "4", "5", "7", "9")
excluded_clusters <- c("1", "2", "6", "8", "10", "11")
observed_clusters <- sort( unique(as.character(obj$myeloid_cluster_r04)) )
if (!all(clean_clusters %in% observed_clusters)) {
  stop( "Validation check.", paste( clean_clusters[ !clean_clusters %in% observed_clusters ], collapse = ", " ) )
}
if (!all(excluded_clusters %in% observed_clusters)) {
  warning( "Expected-value check failed.", paste( excluded_clusters[ !excluded_clusters %in% observed_clusters ], collapse = ", " ) )
}
final_lineage_map <- data.frame( myeloid_cluster_r04 = c( "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11" ), final_lineage_class = c( "canonical_resident_microglia", "infiltrating_macrophage", "infiltrating_monocyte", "activated_resident_microglia", "activated_resident_microglia", "cycling_S_microglia", "macrophage_biased_infiltrating_myeloid", "cycling_G2M_microglia", "BAM_like_macrophage", "IFN_responsive_resident_microglia", "IFN_responsive_infiltrating_mono_macro", "MHCII_APC_like_infiltrating_macrophage" ), final_step2_status = c( "KEEP", "EXCLUDE", "EXCLUDE", "KEEP", "KEEP", "KEEP", "EXCLUDE", "KEEP", "EXCLUDE", "KEEP", "EXCLUDE", "EXCLUDE" ), enter_step3 = c( TRUE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE ), stringsAsFactors = FALSE )
write.csv( final_lineage_map, file.path( tab_dir, "01_final_step2_lineage_map.csv" ), row.names = FALSE )
map_class <- setNames( final_lineage_map$final_lineage_class, final_lineage_map$myeloid_cluster_r04 )
map_status <- setNames( final_lineage_map$final_step2_status, final_lineage_map$myeloid_cluster_r04 )
obj$final_step2_lineage <- unname( map_class[obj$myeloid_cluster_r04] )
obj$final_step2_status <- unname( map_status[obj$myeloid_cluster_r04] )
full_counts <- as.data.frame( table( cluster = obj$myeloid_cluster_r04, lineage = obj$final_step2_lineage, status = obj$final_step2_status ) )
write.csv( full_counts, file.path( tab_dir, "02_full_myeloid_final_lineage_counts.csv" ), row.names = FALSE )
clean_cells <- colnames(obj)[ obj$myeloid_cluster_r04 %in% clean_clusters ]
if (length(clean_cells) == 0) {
  stop("Required result not found.")
}
microglia <- subset( obj, cells = clean_cells )
DefaultAssay(microglia) <- "RNA"
microglia$step2_final_cluster <- factor( as.character(microglia$myeloid_cluster_r04), levels = c("0", "3", "4", "5", "7", "9") )
microglia$step2_final_state_label <- unname( map_class[ as.character(microglia$myeloid_cluster_r04) ] )
Idents(microglia) <- "step2_final_cluster"
print(table(microglia$step2_final_cluster))
expected_total <- sum( c( 11808, 2797, 2196, 1864, 1515, 567 ) )
if (ncol(microglia) != expected_total) {
  warning( "Consistency check failed.", "Expected = ", expected_total, " | Observed = ", ncol(microglia) )
}
cluster_counts <- as.data.frame( table( cluster = microglia$step2_final_cluster ) )
cluster_counts$proportion <- ( cluster_counts$Freq / sum(cluster_counts$Freq) )
write.csv( cluster_counts, file.path( tab_dir, "03_clean_microglia_cluster_counts.csv" ), row.names = FALSE )
sample_counts <- as.data.frame( table( cluster = microglia$step2_final_cluster, sample_id = microglia$sample_id ) )
sample_counts <- sample_counts %>% group_by(cluster) %>% mutate( cluster_total = sum(Freq), proportion_within_cluster = ifelse(cluster_total > 0, Freq / cluster_total, NA_real_) ) %>% ungroup()
write.csv( sample_counts, file.path( tab_dir, "04_clean_microglia_cluster_by_sample.csv" ), row.names = FALSE )
sample_retention <- data.frame( sample_id = unique(as.character(obj$sample_id)), stringsAsFactors = FALSE )
sample_retention$all_step2B_myeloid <- vapply( sample_retention$sample_id, function(s) {
    sum(as.character(obj$sample_id) == s)
  },
  numeric(1)
)
sample_retention$clean_microglia <- vapply( sample_retention$sample_id, function(s) {
    sum(as.character(microglia$sample_id) == s)
  },
  numeric(1)
)
sample_retention$retention_fraction <- ( sample_retention$clean_microglia / sample_retention$all_step2B_myeloid )
sample_retention <- sample_retention[ order(sample_retention$sample_id), ]
write.csv( sample_retention, file.path( tab_dir, "05_clean_microglia_sample_retention.csv" ), row.names = FALSE )
time_counts <- as.data.frame( table( cluster = microglia$step2_final_cluster, timepoint = microglia$timepoint ) )
time_counts <- time_counts %>% group_by(cluster) %>% mutate( cluster_total = sum(Freq), proportion_within_cluster = ifelse(cluster_total > 0, Freq / cluster_total, NA_real_) ) %>% ungroup()
write.csv( time_counts, file.path( tab_dir, "06_clean_microglia_cluster_by_timepoint.csv" ), row.names = FALSE )
time_retention <- data.frame( timepoint = unique(as.character(obj$timepoint)), stringsAsFactors = FALSE )
time_retention$all_step2B_myeloid <- vapply( time_retention$timepoint, function(tp) {
    sum(as.character(obj$timepoint) == tp)
  },
  numeric(1)
)
time_retention$clean_microglia <- vapply( time_retention$timepoint, function(tp) {
    sum(as.character(microglia$timepoint) == tp)
  },
  numeric(1)
)
time_retention$retention_fraction <- ( time_retention$clean_microglia / time_retention$all_step2B_myeloid )
write.csv( time_retention, file.path( tab_dir, "07_clean_microglia_timepoint_retention.csv" ), row.names = FALSE )
author_col <- "original_V30_step1F"
if (author_col %in% colnames(microglia[[]])) {
  author_cross <- as.data.frame( table( cluster = microglia$step2_final_cluster, author_annotation = microglia[[author_col]][, 1], useNA = "ifany" ) )
  author_cross <- author_cross %>% group_by(cluster) %>% mutate( cluster_total = sum(Freq), proportion_within_cluster = ifelse(cluster_total > 0, Freq / cluster_total, NA_real_) ) %>% ungroup()
  write.csv( author_cross, file.path( tab_dir, "08_clean_microglia_vs_author_annotation.csv" ), row.names = FALSE )
}
if (!"umap.myeloid" %in% Reductions(microglia)) {
  stop( "Required input not found: ", "Validation check." )
}
p1 <- DimPlot( microglia, reduction = "umap.myeloid", group.by = "step2_final_cluster", label = TRUE, repel = TRUE ) + ggtitle("Final clean resident-derived microglia | Step 2C")
p2 <- DimPlot( microglia, reduction = "umap.myeloid", group.by = "step2_final_state_label", label = TRUE, repel = TRUE ) + ggtitle("Final Step 2 lineage labels")
p3 <- DimPlot( microglia, reduction = "umap.myeloid", group.by = "sample_id" ) + ggtitle("Clean microglia by biological sample")
p4 <- DimPlot( microglia, reduction = "umap.myeloid", group.by = "timepoint" ) + ggtitle("Clean microglia by SCI timepoint")
p_all <- (p1 | p2) / (p3 | p4)
ggsave( file.path( fig_dir, "01_final_clean_microglia_UMAP_overview.pdf" ), p_all, width = 16, height = 11 )
cluster_prop_plot <- cluster_counts
p_cluster_prop <- ggplot( cluster_prop_plot, aes( x = cluster, y = proportion ) ) + geom_col() + theme_bw() + labs( title = "Final clean microglia cluster composition", x = "Step 2 final cluster", y = "Proportion" )
ggsave( file.path( fig_dir, "02_clean_microglia_cluster_proportion.pdf" ), p_cluster_prop, width = 7, height = 5 )
time_mat <- xtabs( Freq ~ cluster + timepoint, data = time_counts )
time_prop_mat <- prop.table( time_mat, margin = 1 )
pdf( file.path( fig_dir, "03_clean_microglia_cluster_timepoint_heatmap.pdf" ), width = 8, height = 6 )
pheatmap( time_prop_mat, cluster_rows = FALSE, cluster_cols = FALSE, main = "Clean microglia cluster × timepoint | row proportion" )
dev.off()
sample_mat <- xtabs( Freq ~ cluster + sample_id, data = sample_counts )
sample_prop_mat <- prop.table( sample_mat, margin = 1 )
pdf( file.path( fig_dir, "04_clean_microglia_cluster_sample_heatmap.pdf" ), width = 10, height = 6 )
pheatmap( sample_prop_mat, cluster_rows = FALSE, cluster_cols = FALSE, main = "Clean microglia cluster × biological sample | row proportion" )
dev.off()
metadata_keep_check <- c( "myeloid_cluster_r04", "broad_origin", "sample_id", "timepoint", "original_V30_step1F", "final_step2_lineage", "final_step2_status", "step2_final_cluster", "step2_final_state_label" )
metadata_status <- data.frame( metadata_field = metadata_keep_check, present = metadata_keep_check %in% colnames(microglia[[]]), stringsAsFactors = FALSE )
write.csv( metadata_status, file.path( tab_dir, "09_final_object_metadata_check.csv" ), row.names = FALSE )
final_rds <- file.path( out_dir, "GSE162610_microglia_clean.rds" )
saveRDS( microglia, final_rds )
full_final_rds <- file.path( out_dir, "GSE162610_step2_full_myeloid_with_final_lineage.rds" )
saveRDS( obj, full_final_rds )
report_file <- file.path( out_dir, "10_STEP2C_final_microglia_extraction_report.txt" )
sink(report_file)
print(table(microglia$step2_final_cluster))
print( table( microglia$step2_final_cluster, microglia$step2_final_state_label ) )
print(sample_retention)
print(time_retention)
if (author_col %in% colnames(microglia[[]])) {
  print( sort( table( microglia[[author_col]][, 1], useNA = "ifany" ), decreasing = TRUE ) )
}
if (ncol(microglia) == expected_total) {
} else {
}
sink()
handoff_file <- file.path( out_dir, "11_STEP3_handoff.txt" )
sink(handoff_file)
sink()
