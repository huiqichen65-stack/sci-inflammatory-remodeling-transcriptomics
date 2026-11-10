options(stringsAsFactors = FALSE)
set.seed(1234)
required_pkgs <- c( "Seurat", "SeuratObject", "dplyr", "ggplot2", "patchwork", "AUCell", "GSVA", "pheatmap", "Matrix" )
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
  library(AUCell)
  library(GSVA)
  library(pheatmap)
  library(Matrix)
})
project_dir <- Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
step3A_dir <- file.path( project_dir, "04_results", "step3A_GSE162610_microglia_reclustering" )
input_rds <- file.path( step3A_dir, "GSE162610_step3A_microglia_reclustered_multires.rds" )
step3A_top20_file <- file.path( step3A_dir, "tables", "07_top20_markers_per_cluster_provisional_r04.csv" )
out_dir <- file.path( project_dir, "04_results", "step3B_GSE162610_state_scoring" )
fig_dir <- file.path(out_dir, "figures")
tab_dir <- file.path(out_dir, "tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(input_rds)) {
  stop( "Required input not found: ", input_rds )
}
obj_all <- readRDS(input_rds)
if (!inherits(obj_all, "Seurat")) {
  stop("Input object is not a Seurat object.")
}
DefaultAssay(obj_all) <- "RNA"
required_meta <- c( "step3A_cluster_r04", "sample_id", "timepoint", "step2_final_cluster" )
missing_meta <- required_meta[ !required_meta %in% colnames(obj_all[[]]) ]
if (length(missing_meta) > 0) {
  stop( "Missing metadata: ", paste(missing_meta, collapse = ", ") )
}
obj_all$step3A_cluster_r04 <- as.character( obj_all$step3A_cluster_r04 )
print(table(obj_all$step3A_cluster_r04))
remove_clusters <- c("7", "11")
keep_clusters <- c( "0", "1", "2", "3", "4", "5", "6", "8", "9", "10" )
observed_clusters <- sort( unique( as.character( obj_all$step3A_cluster_r04 ) ) )
if (!all(remove_clusters %in% observed_clusters)) {
  stop( "Expected-value check failed.", paste( remove_clusters[ !remove_clusters %in% observed_clusters ], collapse = ", " ) )
}
if (!all(keep_clusters %in% observed_clusters)) {
  stop( "Expected-value check failed.", paste( keep_clusters[ !keep_clusters %in% observed_clusters ], collapse = ", " ) )
}
cleanup_map <- data.frame( step3A_cluster_r04 = observed_clusters, step3B_status = ifelse( observed_clusters %in% remove_clusters, "REMOVE", "KEEP" ), step3B_reason = ifelse( observed_clusters == "7", "residual_monocyte_like_contamination", ifelse( observed_clusters == "11", "residual_astrocyte_glial_contamination", "retained_microglial_state_candidate" ) ), stringsAsFactors = FALSE )
write.csv( cleanup_map, file.path( tab_dir, "01_step3B_residual_cleanup_map.csv" ), row.names = FALSE )
obj_all$step3B_status <- ifelse( obj_all$step3A_cluster_r04 %in% remove_clusters, "REMOVE", "KEEP" )
obj_all$step3B_cleanup_reason <- ifelse( obj_all$step3A_cluster_r04 == "7", "residual_monocyte_like_contamination", ifelse( obj_all$step3A_cluster_r04 == "11", "residual_astrocyte_glial_contamination", "retained_microglial_state_candidate" ) )
removed_counts <- as.data.frame( table( cluster = obj_all$step3A_cluster_r04[ obj_all$step3B_status == "REMOVE" ] ) )
write.csv( removed_counts, file.path( tab_dir, "02_removed_residual_contaminant_counts.csv" ), row.names = FALSE )
keep_cells <- colnames(obj_all)[ obj_all$step3B_status == "KEEP" ]
obj <- subset( obj_all, cells = keep_cells )
DefaultAssay(obj) <- "RNA"
obj$step3A_cluster_r04 <- factor( as.character(obj$step3A_cluster_r04), levels = keep_clusters )
Idents(obj) <- "step3A_cluster_r04"
expected_clean_n <- 20747 - 620 - 169
if (ncol(obj) != expected_clean_n) {
  warning( "Step 3B clean cell number differs from expected value: ", expected_clean_n, " | observed = ", ncol(obj) )
}
if (!"umap.microglia.step3A" %in% Reductions(obj)) {
  stop( "Step 3A formal UMAP reduction umap.microglia.step3A is missing." )
}
p_cluster <- DimPlot( obj, reduction = "umap.microglia.step3A", group.by = "step3A_cluster_r04", label = TRUE, repel = TRUE ) + ggtitle( "Step 3B retained microglial state candidates" )
p_time <- DimPlot( obj, reduction = "umap.microglia.step3A", group.by = "timepoint" ) + ggtitle("SCI timepoint")
p_sample <- DimPlot( obj, reduction = "umap.microglia.step3A", group.by = "sample_id" ) + ggtitle("Biological sample")
p_origin <- DimPlot( obj, reduction = "umap.microglia.step3A", group.by = "step2_final_cluster", label = TRUE, repel = TRUE ) + ggtitle("Step 2 lineage provenance")
ggsave( file.path( fig_dir, "01_step3B_cleaned_microglia_UMAP_overview.pdf" ), (p_cluster | p_time) / (p_sample | p_origin), width = 15, height = 11 )
gene_sets <- list( Homeostatic = c( "P2ry12", "Tmem119", "Sall1", "Fcrls", "Hexb", "Gpr34", "P2ry13", "Siglech", "Cx3cr1" ), Activated_injury = c( "Apoe", "Lpl", "Cst7", "Cd68", "Ctsb", "Ctsd", "Cd9", "Tyrobp", "Trem2" ), DAM_like_candidate = c( "Trem2", "Apoe", "Lpl", "Cst7", "Axl", "Clec7a", "Itgax", "Gpnmb", "Spp1", "Lgals3", "Cd9", "Tyrobp" ), Phagocytic_lipid = c( "Lpl", "Apoe", "Cd68", "Ctsd", "Lgals3", "Fabp5", "Gpnmb", "Spp1", "Msr1", "Cd9" ), Early_inflammatory = c( "Il1b", "Il1a", "Il6", "Ccl7", "Cxcl3", "Acod1", "Ptges", "Fosl1", "Nfkbia", "Junb" ), IFN_response = c( "Ifit1", "Ifit2", "Ifit3", "Isg15", "Mx1", "Rsad2", "Usp18", "Oasl1" ), Cycling_S = c( "Mcm2", "Mcm5", "Mcm6", "Pcna", "Ccne2", "Uhrf1", "Cdc6", "Dtl", "Rad51" ), Cycling_G2M = c( "Mki67", "Top2a", "Ube2c", "Birc5", "Cenpf", "Ccna2", "Aurkb", "Cdc20", "Ccnb1" ) )
rna_genes <- rownames(obj[["RNA"]])
gene_sets_use <- lapply( gene_sets, function(gs) {
    intersect(gs, rna_genes)
  }
)
gene_set_availability <- do.call( rbind, lapply( names(gene_sets), function(nm) {
      data.frame( program = nm, n_requested = length(gene_sets[[nm]]), n_available = length(gene_sets_use[[nm]]), available_genes = paste( gene_sets_use[[nm]], collapse = ";" ), stringsAsFactors = FALSE )
    }
  )
)
write.csv( gene_set_availability, file.path( tab_dir, "03_state_program_gene_set_availability.csv" ), row.names = FALSE )
bad_sets <- names(gene_sets_use)[ vapply( gene_sets_use, length, integer(1) ) < 3 ]
if (length(bad_sets) > 0) {
  stop( "Threshold check failed.", paste(bad_sets, collapse = ", ") )
}
dot_features <- unique( unlist( gene_sets_use, use.names = FALSE ) )
p_dot <- DotPlot( obj, features = dot_features, group.by = "step3A_cluster_r04", assay = "RNA" ) + RotatedAxis() + ggtitle( "Formal microglial state-program markers | Step 3B" )
ggsave( file.path( fig_dir, "02_formal_state_program_marker_DotPlot.pdf" ), p_dot, width = 22, height = 8 )
get_rna_matrix <- function(seurat_obj) {
  x <- tryCatch( GetAssayData( seurat_obj, assay = "RNA", layer = "counts" ), error = function(e) NULL )
  if (is.null(x)) {
    x <- tryCatch( GetAssayData( seurat_obj, assay = "RNA", slot = "counts" ), error = function(e) NULL )
  }
  if (is.null(x)) {
    x <- tryCatch( GetAssayData( seurat_obj, assay = "RNA", layer = "data" ), error = function(e) NULL )
  }
  if (is.null(x)) {
    stop("Validation failed.")
  }
  x
}
expr_auc <- get_rna_matrix(obj)
if (!inherits(expr_auc, "dgCMatrix")) {
  expr_auc <- as( expr_auc, "dgCMatrix" )
}
keep_gene_auc <- Matrix::rowSums( expr_auc ) > 0
expr_auc <- expr_auc[ keep_gene_auc, , drop = FALSE ]
cells_rankings <- AUCell_buildRankings( expr_auc, plotStats = FALSE, nCores = 1, verbose = TRUE )
auc_max_rank <- max( 50, ceiling( 0.05 * nrow(expr_auc) ) )
cells_AUC <- tryCatch( AUCell_calcAUC( gene_sets_use, cells_rankings, aucMaxRank = auc_max_rank, nCores = 1 ), error = function(e1) {
    AUCell_calcAUC( gene_sets_use, cells_rankings, aucMaxRank = auc_max_rank )
  }
)
auc_mat <- as.matrix( getAUC(cells_AUC) )
auc_mat <- auc_mat[ , colnames(obj), drop = FALSE ]
for (nm in rownames(auc_mat)) {
  obj[[paste0("AUCell_", nm)]] <- as.numeric( auc_mat[nm, ] )
}
cluster_levels <- keep_clusters
auc_summary_list <- lapply( cluster_levels, function(cl) {
    idx <- which( as.character( obj$step3A_cluster_r04 ) == cl )
    out <- data.frame( cluster = cl, n_cells = length(idx), stringsAsFactors = FALSE )
    for (nm in rownames(auc_mat)) {
      vals <- auc_mat[ nm, idx ]
      out[[paste0(nm, "_mean")]] <- mean( vals, na.rm = TRUE )
      out[[paste0(nm, "_median")]] <- median( vals, na.rm = TRUE )
    }
    out
  }
)
auc_cluster_summary <- bind_rows( auc_summary_list )
write.csv( auc_cluster_summary, file.path( tab_dir, "04_AUCell_state_program_cluster_summary.csv" ), row.names = FALSE )
auc_cluster_median <- sapply( cluster_levels, function(cl) {
    idx <- which( as.character( obj$step3A_cluster_r04 ) == cl )
    apply( auc_mat[ , idx, drop = FALSE ], 1, median, na.rm = TRUE )
  }
)
colnames(auc_cluster_median) <- cluster_levels
row_z <- function(x) {
  sx <- sd( x, na.rm = TRUE )
  if ( !is.finite(sx) || sx == 0 ) {
    return( rep( 0, length(x) ) )
  }
  ( x - mean( x, na.rm = TRUE ) ) / sx
}
auc_cluster_z <- t( apply( auc_cluster_median, 1, row_z ) )
write.csv( data.frame( program = rownames(auc_cluster_z), auc_cluster_z, check.names = FALSE ), file.path( tab_dir, "05_AUCell_cluster_Zscores.csv" ), row.names = FALSE )
pdf( file.path( fig_dir, "03_AUCell_state_program_heatmap.pdf" ), width = 11, height = 7 )
pheatmap( auc_cluster_z, cluster_rows = FALSE, cluster_cols = FALSE, main = "AUCell state programs | cluster median Z-score" )
dev.off()
auc_features <- paste0( "AUCell_", rownames(auc_mat) )
p_auc <- FeaturePlot( obj, features = auc_features, reduction = "umap.microglia.step3A", ncol = 4, order = TRUE )
ggsave( file.path( fig_dir, "04_AUCell_state_program_UMAP.pdf" ), p_auc, width = 17, height = 10 )
avg_expr <- AverageExpression( obj, assays = "RNA", group.by = "step3A_cluster_r04", layer = "data", verbose = FALSE )[["RNA"]]
avg_expr <- as.matrix( avg_expr )
run_ssgsea_compat <- function( expr_matrix, gene_sets_list ) {
  ns_exports <- getNamespaceExports( "GSVA" )
  if ( any( ns_exports == "ssgseaParam" ) ) {
    param <- tryCatch( GSVA::ssgseaParam( exprData = expr_matrix, geneSets = gene_sets_list, minSize = 2, maxSize = Inf ), error = function(e1) {
        GSVA::ssgseaParam( expr_matrix, gene_sets_list, minSize = 2, maxSize = Inf )
      }
    )
    out <- GSVA::gsva( param, verbose = FALSE )
  } else {
    out <- GSVA::gsva( expr_matrix, gene_sets_list, method = "ssgsea", kcdf = "Gaussian", abs.ranking = TRUE, verbose = FALSE )
  }
  as.matrix(out)
}
ssgsea_mat <- run_ssgsea_compat( avg_expr, gene_sets_use )
write.csv( data.frame( program = rownames(ssgsea_mat), ssgsea_mat, check.names = FALSE ), file.path( tab_dir, "06_cluster_average_ssGSEA_scores.csv" ), row.names = FALSE )
ssgsea_z <- t( apply( ssgsea_mat, 1, row_z ) )
write.csv( data.frame( program = rownames(ssgsea_z), ssgsea_z, check.names = FALSE ), file.path( tab_dir, "07_ssGSEA_cluster_Zscores.csv" ), row.names = FALSE )
pdf( file.path( fig_dir, "05_ssGSEA_state_program_heatmap.pdf" ), width = 11, height = 7 )
pheatmap( ssgsea_z, cluster_rows = FALSE, cluster_cols = FALSE, main = "Cluster-average ssGSEA state programs | row Z-score" )
dev.off()
cluster_time <- as.data.frame( table( cluster = obj$step3A_cluster_r04, timepoint = obj$timepoint ) )
cluster_time <- cluster_time %>% group_by(cluster) %>% mutate( cluster_total = sum(Freq), proportion_within_cluster = ifelse( cluster_total > 0, Freq / cluster_total, NA_real_ ) ) %>% ungroup()
write.csv( cluster_time, file.path( tab_dir, "08_cluster_by_timepoint_after_cleanup.csv" ), row.names = FALSE )
time_mat <- xtabs( Freq ~ cluster + timepoint, data = cluster_time )
time_prop <- prop.table( time_mat, margin = 1 )
pdf( file.path( fig_dir, "06_cluster_timepoint_heatmap_after_cleanup.pdf" ), width = 8, height = 7 )
pheatmap( time_prop, cluster_rows = FALSE, cluster_cols = FALSE, main = "Step 3B retained clusters x timepoint | row proportion" )
dev.off()
cluster_sample <- as.data.frame( table( cluster = obj$step3A_cluster_r04, sample_id = obj$sample_id ) )
cluster_sample <- cluster_sample %>% group_by(cluster) %>% mutate( cluster_total = sum(Freq), proportion_within_cluster = ifelse( cluster_total > 0, Freq / cluster_total, NA_real_ ) ) %>% ungroup()
write.csv( cluster_sample, file.path( tab_dir, "09_cluster_by_sample_after_cleanup.csv" ), row.names = FALSE )
cluster_provenance <- as.data.frame( table( step3A_cluster = obj$step3A_cluster_r04, step2_cluster = obj$step2_final_cluster ) )
cluster_provenance <- cluster_provenance %>% group_by(step3A_cluster) %>% mutate( cluster_total = sum(Freq), proportion_within_cluster = ifelse( cluster_total > 0, Freq / cluster_total, NA_real_ ) ) %>% ungroup()
write.csv( cluster_provenance, file.path( tab_dir, "10_cluster_vs_step2_provenance_after_cleanup.csv" ), row.names = FALSE )
if (file.exists(step3A_top20_file)) {
  top20_step3A <- read.csv( step3A_top20_file, stringsAsFactors = FALSE, check.names = FALSE )
  cluster_col_candidates <- intersect( c( "cluster", "Cluster", "ident" ), colnames(top20_step3A) )
  if (length(cluster_col_candidates) > 0) {
    marker_cluster_col <- cluster_col_candidates[1]
    top20_keep <- top20_step3A[ as.character( top20_step3A[[marker_cluster_col]] ) %in% keep_clusters, , drop = FALSE ]
    write.csv( top20_keep, file.path( tab_dir, "11_top20_markers_retained_clusters.csv" ), row.names = FALSE )
  }
}
get_top_timepoint <- function(cl) {
  x <- cluster_time[ as.character( cluster_time$cluster ) == cl, , drop = FALSE ]
  if (nrow(x) == 0) {
    return(NA_character_)
  }
  as.character( x$timepoint[ which.max( x$proportion_within_cluster ) ] )
}
get_top_provenance <- function(cl) {
  x <- cluster_provenance[ as.character( cluster_provenance$step3A_cluster ) == cl, , drop = FALSE ]
  if (nrow(x) == 0) {
    return(NA_character_)
  }
  as.character( x$step2_cluster[ which.max( x$proportion_within_cluster ) ] )
}
get_top_program <- function( zmat, cl ) {
  candidate_cols <- c( cl, paste0( "g", cl ) )
  hit <- candidate_cols[ candidate_cols %in% colnames(zmat) ]
  if (length(hit) == 0) {
    return(NA_character_)
  }
  vals <- zmat[ , hit[1] ]
  names(vals)[ which.max(vals) ]
}
get_zscore <- function( zmat, program, cl ) {
  if ( !program %in% rownames(zmat) ) {
    return(NA_real_)
  }
  candidates <- c( cl, paste0( "g", cl ) )
  hit <- candidates[ candidates %in% colnames(zmat) ]
  if (length(hit) == 0) {
    return(NA_real_)
  }
  as.numeric( zmat[ program, hit[1] ] )
}
evidence_tbl <- data.frame( cluster = keep_clusters, n_cells = vapply( keep_clusters, function(cl) {
      sum( as.character( obj$step3A_cluster_r04 ) == cl )
    },
    integer(1)
  ), top_timepoint = vapply(
    keep_clusters, function(cl) as.character(get_top_timepoint(cl)), character(1)
  ), dominant_step2_provenance = vapply(
    keep_clusters, function(cl) as.character(get_top_provenance(cl)), character(1)
  ), AUCell_top_program = vapply(
    keep_clusters, function(cl) {
      get_top_program( auc_cluster_z, cl )
    },
    character(1)
  ), ssGSEA_top_program = vapply(
    keep_clusters, function(cl) {
      get_top_program( ssgsea_z, cl )
    },
    character(1)
  ), stringsAsFactors = FALSE
)
for (program in names(gene_sets_use)) {
  evidence_tbl[[ paste0( "AUCellZ_", program ) ]] <- vapply( keep_clusters, function(cl) {
      get_zscore( auc_cluster_z, program, cl )
    },
    numeric(1)
  )
  evidence_tbl[[ paste0( "ssGSEAZ_", program ) ]] <- vapply( keep_clusters, function(cl) {
      get_zscore( ssgsea_z, program, cl )
    },
    numeric(1)
  )
}
evidence_tbl$final_state_annotation <- NA_character_
evidence_tbl$reviewer_note <- NA_character_
write.csv( evidence_tbl, file.path( tab_dir, "12_cluster_state_evidence_matrix.csv" ), row.names = FALSE )
output_rds <- file.path( out_dir, "GSE162610_step3B_clean_state_scored.rds" )
saveRDS( obj, output_rds )
full_output_rds <- file.path( out_dir, "GSE162610_step3B_full_with_cleanup_flags.rds" )
saveRDS( obj_all, full_output_rds )
report_file <- file.path( out_dir, "13_STEP3B_state_scoring_report.txt" )
sink(report_file)
print(removed_counts)
print( table( obj$step3A_cluster_r04 ) )
print(gene_set_availability)
print(evidence_tbl)
if (ncol(obj) == expected_clean_n) {
} else {
}
sink()
handoff_file <- file.path( out_dir, "14_STEP3C_handoff.txt" )
sink(handoff_file)
sink()
