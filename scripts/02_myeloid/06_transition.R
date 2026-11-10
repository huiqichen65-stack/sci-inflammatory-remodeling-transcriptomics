options(stringsAsFactors = FALSE)
set.seed(1234)
pkgs <- c("Seurat","SeuratObject","dplyr","ggplot2","patchwork", "AUCell","Matrix","pheatmap")
miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss) > 0) stop("缺少R包: ", paste(miss, collapse=", "))
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(AUCell)
  library(Matrix)
  library(pheatmap)
})
project_dir <- Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
input_rds <- file.path( project_dir, "04_results", "step3B_GSE162610_state_scoring", "GSE162610_step3B_clean_state_scored.rds" )
out_dir <- file.path( project_dir, "04_results", "step3C_GSE162610_activated_DAM_transition_mapping" )
fig_dir <- file.path(out_dir, "figures")
tab_dir <- file.path(out_dir, "tables")
dir.create(fig_dir, recursive=TRUE, showWarnings=FALSE)
dir.create(tab_dir, recursive=TRUE, showWarnings=FALSE)
if (!file.exists(input_rds)) stop("找不到输入RDS:\n", input_rds)
message("====================================================")
message("GSE162610 Step 3C | activated-DAM transition mapping")
message("====================================================")
obj <- readRDS(input_rds)
if (!inherits(obj, "Seurat")) stop("输入对象不是Seurat对象")
DefaultAssay(obj) <- "RNA"
req_meta <- c("step3A_cluster_r04","sample_id","timepoint","step2_final_cluster")
miss_meta <- req_meta[!req_meta %in% colnames(obj[[]])]
if (length(miss_meta)>0) stop("缺少metadata: ", paste(miss_meta, collapse=", "))
obj$step3A_cluster_r04 <- as.character(obj$step3A_cluster_r04)
message("Input cells: ", ncol(obj))
print(table(obj$step3A_cluster_r04))
print(table(obj$timepoint))
print(table(obj$sample_id))
if (ncol(obj) != 19958) {
  warning("Step3C输入细胞数不是预期19958，observed=", ncol(obj))
}
if (!"umap.microglia.step3A" %in% Reductions(obj)) {
  stop("缺少 umap.microglia.step3A reduction")
}
gene_sets <- list( Homeostatic = c( "P2ry12","Tmem119","Sall1","Fcrls","Hexb", "Gpr34","P2ry13","Siglech","Cx3cr1" ), Activated_injury_core = c( "Apoe","Cst7","Cd68","Ctsb","Ctsd","Cd9","Tyrobp","Trem2" ), Early_DAM_transition = c( "Trem2","Apoe","Lpl","Cst7","Cd9","Tyrobp","Axl" ), Mature_DAM_phagocytic = c( "Gpnmb","Spp1","Lgals3","Fabp5","Cd68","Ctsd", "Lpl","Itgax","Clec7a" ) )
rna_genes <- rownames(obj[["RNA"]])
gene_sets_use <- lapply(gene_sets, intersect, y=rna_genes)
avail <- bind_rows(lapply(names(gene_sets), function(nm) {
  data.frame( program=nm, n_requested=length(gene_sets[[nm]]), n_available=length(gene_sets_use[[nm]]), available_genes=paste(gene_sets_use[[nm]], collapse=";") )
}))
write.csv(avail, file.path(tab_dir,"01_transition_program_gene_availability.csv"), row.names=FALSE)
bad <- names(gene_sets_use)[vapply(gene_sets_use,length,integer(1)) < 3]
if (length(bad)>0) stop("以下gene set可用基因<3: ", paste(bad,collapse=", "))
message("Starting AUCell transition scoring...")
get_rna <- function(x) {
  m <- tryCatch(GetAssayData(x, assay="RNA", layer="counts"), error=function(e) NULL)
  if (is.null(m)) {
    m <- tryCatch(GetAssayData(x, assay="RNA", slot="counts"), error=function(e) NULL)
  }
  if (is.null(m)) {
    m <- tryCatch(GetAssayData(x, assay="RNA", layer="data"), error=function(e) NULL)
  }
  if (is.null(m)) stop("无法取得RNA counts/data")
  m
}
expr <- get_rna(obj)
if (!inherits(expr,"dgCMatrix")) expr <- as(expr,"dgCMatrix")
expr <- expr[Matrix::rowSums(expr)>0,,drop=FALSE]
message("AUCell sparse matrix: ", nrow(expr), " genes x ", ncol(expr), " cells")
rankings <- AUCell_buildRankings(expr, plotStats=FALSE, verbose=TRUE)
auc_max_rank <- max(50, ceiling(0.05*nrow(expr)))
cells_auc <- AUCell_calcAUC(gene_sets_use, rankings, aucMaxRank=auc_max_rank)
auc <- as.matrix(getAUC(cells_auc))
auc <- auc[,colnames(obj),drop=FALSE]
for (nm in rownames(auc)) {
  obj[[paste0("TransAUC_",nm)]] <- as.numeric(auc[nm,])
}
homeo_col <- "TransAUC_Homeostatic"
act_col <- "TransAUC_Activated_injury_core"
early_col <- "TransAUC_Early_DAM_transition"
mature_col <- "TransAUC_Mature_DAM_phagocytic"
md <- obj[[]]
homeo_idx <- as.character(md$step3A_cluster_r04)=="3"
mature_idx <- as.character(md$step3A_cluster_r04)=="1"
if (sum(homeo_idx)<100) stop("cluster3 reference cells <100")
if (sum(mature_idx)<100) stop("cluster1 reference cells <100")
qsafe <- function(x,p) as.numeric(quantile(x,p,na.rm=TRUE,names=FALSE,type=7))
act_thr <- qsafe(md[[act_col]][homeo_idx],0.90)
early_thr <- qsafe(md[[early_col]][homeo_idx],0.90)
mature_thr <- qsafe(md[[mature_col]][mature_idx],0.25)
homeo_thr <- qsafe(md[[homeo_col]][homeo_idx],0.50)
thr_tbl <- data.frame( threshold=c( "Activated_high_q90_cluster3", "EarlyDAM_high_q90_cluster3", "MatureDAM_established_q25_cluster1", "Homeostatic_median_cluster3" ), reference_cluster=c("3","3","1","3"), value=c(act_thr,early_thr,mature_thr,homeo_thr) )
write.csv(thr_tbl,file.path(tab_dir,"02_main_transition_thresholds.csv"), row.names=FALSE)
obj$Activated_high <- obj[[act_col]][,1] > act_thr
obj$EarlyDAM_high <- obj[[early_col]][,1] > early_thr
obj$MatureDAM_established <- obj[[mature_col]][,1] >= mature_thr
obj$Homeostatic_high <- obj[[homeo_col]][,1] >= homeo_thr
side_clusters <- c("6","8","9")
obj$transition_core_eligible <- !obj$step3A_cluster_r04 %in% side_clusters
status <- rep("Other",ncol(obj))
names(status) <- colnames(obj)
A <- obj$Activated_high
E <- obj$EarlyDAM_high
M <- obj$MatureDAM_established
H <- obj$Homeostatic_high
status[E & M] <- "Mature_DAM_like"
status[A & E & !M] <- "Activated_EarlyDAM_transition"
status[A & !E & !M] <- "Activated_only"
status[!A & E & !M] <- "EarlyDAM_only"
status[H & !A & !E & !M] <- "Homeostatic_like"
obj$transition_status <- factor( status[colnames(obj)], levels=c( "Homeostatic_like","Activated_only","EarlyDAM_only", "Activated_EarlyDAM_transition","Mature_DAM_like","Other" ) )
obj$shared_transition_3dpi_broad <- ( as.character(obj$transition_status)=="Activated_EarlyDAM_transition" & as.character(obj$timepoint)=="3 dpi" )
obj$shared_transition_3dpi_highconf <- ( obj$shared_transition_3dpi_broad & obj$transition_core_eligible )
status_counts <- as.data.frame(table(transition_status=obj$transition_status))
status_counts$proportion <- status_counts$Freq/sum(status_counts$Freq)
write.csv(status_counts,file.path(tab_dir,"03_transition_status_counts.csv"), row.names=FALSE)
status_time <- as.data.frame(table( timepoint=obj$timepoint, transition_status=obj$transition_status )) %>% group_by(timepoint) %>% mutate(timepoint_total=sum(Freq), proportion_within_timepoint=Freq/timepoint_total) %>% ungroup()
write.csv(status_time,file.path(tab_dir,"04_transition_status_by_timepoint.csv"), row.names=FALSE)
status_cluster <- as.data.frame(table( cluster=obj$step3A_cluster_r04, transition_status=obj$transition_status )) %>% group_by(cluster) %>% mutate(cluster_total=sum(Freq), proportion_within_cluster=Freq/cluster_total) %>% ungroup()
write.csv(status_cluster,file.path(tab_dir,"05_transition_status_by_cluster.csv"), row.names=FALSE)
high_cells <- colnames(obj)[obj$shared_transition_3dpi_highconf]
high_n <- length(high_cells)
message("High-confidence 3 dpi shared-transition candidates: ", high_n)
high_cluster <- as.data.frame(table( cluster=obj$step3A_cluster_r04[obj$shared_transition_3dpi_highconf] ))
if (sum(high_cluster$Freq)>0) high_cluster$proportion <- high_cluster$Freq/sum(high_cluster$Freq)
write.csv(high_cluster,file.path(tab_dir,"06_highconf_3dpi_transition_by_cluster.csv"), row.names=FALSE)
high_sample <- as.data.frame(table( sample_id=obj$sample_id[obj$shared_transition_3dpi_highconf] ))
if (sum(high_sample$Freq)>0) high_sample$proportion <- high_sample$Freq/sum(high_sample$Freq)
write.csv(high_sample,file.path(tab_dir,"07_highconf_3dpi_transition_by_sample.csv"), row.names=FALSE)
hm <- obj[[]][obj$shared_transition_3dpi_highconf,,drop=FALSE]
cols_keep <- intersect(c( "sample_id","timepoint","step3A_cluster_r04","step2_final_cluster", homeo_col,act_col,early_col,mature_col,"transition_status" ),colnames(hm))
high_export <- data.frame( cell_barcode=rownames(hm), hm[,cols_keep,drop=FALSE], check.names=FALSE )
write.csv(high_export,file.path(tab_dir,"08_highconf_3dpi_transition_cells.csv"), row.names=FALSE)
time_summary <- obj[[]] %>% mutate(timepoint_chr=as.character(timepoint)) %>% group_by(timepoint_chr) %>% summarise( n_cells=n(), Homeostatic_median=median(.data[[homeo_col]],na.rm=TRUE), Activated_median=median(.data[[act_col]],na.rm=TRUE), EarlyDAM_median=median(.data[[early_col]],na.rm=TRUE), MatureDAM_median=median(.data[[mature_col]],na.rm=TRUE), activated_high_fraction=mean(Activated_high,na.rm=TRUE), earlyDAM_high_fraction=mean(EarlyDAM_high,na.rm=TRUE), mature_DAM_fraction=mean(MatureDAM_established,na.rm=TRUE), transition_fraction=mean( as.character(transition_status)=="Activated_EarlyDAM_transition", na.rm=TRUE ), .groups="drop" )
write.csv(time_summary,file.path(tab_dir,"09_transition_score_summary_by_timepoint.csv"), row.names=FALSE)
three <- obj[[]] %>% filter(as.character(timepoint)=="3 dpi") %>% group_by(sample_id) %>% summarise( n_cells=n(), n_transition_broad=sum(shared_transition_3dpi_broad), fraction_transition_broad=mean(shared_transition_3dpi_broad), n_transition_highconf=sum(shared_transition_3dpi_highconf), fraction_transition_highconf=mean(shared_transition_3dpi_highconf), Activated_median=median(.data[[act_col]],na.rm=TRUE), EarlyDAM_median=median(.data[[early_col]],na.rm=TRUE), MatureDAM_median=median(.data[[mature_col]],na.rm=TRUE), .groups="drop" )
write.csv(three,file.path(tab_dir,"10_3dpi_transition_by_biological_sample.csv"), row.names=FALSE)
cluster_summary <- obj[[]] %>% mutate(cluster_chr=as.character(step3A_cluster_r04)) %>% group_by(cluster_chr) %>% summarise( n_cells=n(), Homeostatic_median=median(.data[[homeo_col]],na.rm=TRUE), Activated_median=median(.data[[act_col]],na.rm=TRUE), EarlyDAM_median=median(.data[[early_col]],na.rm=TRUE), MatureDAM_median=median(.data[[mature_col]],na.rm=TRUE), transition_fraction=mean( as.character(transition_status)=="Activated_EarlyDAM_transition" ), mature_DAM_fraction=mean( as.character(transition_status)=="Mature_DAM_like" ), .groups="drop" )
write.csv(cluster_summary,file.path(tab_dir,"11_transition_score_summary_by_cluster.csv"), row.names=FALSE)
sens <- bind_rows(lapply(c(0.85,0.90,0.95), function(pp) {
  at <- qsafe(md[[act_col]][homeo_idx],pp)
  et <- qsafe(md[[early_col]][homeo_idx],pp)
  af <- obj[[act_col]][,1] > at
  ef <- obj[[early_col]][,1] > et
  mf <- obj[[mature_col]][,1] >= mature_thr
  tf <- af & ef & !mf
  data.frame( timepoint=as.character(obj$timepoint), transition=tf ) %>% group_by(timepoint) %>% summarise( n_cells=n(), n_transition=sum(transition), transition_fraction=mean(transition), .groups="drop" ) %>% mutate( activated_early_reference_quantile=pp, activated_threshold=at, earlyDAM_threshold=et, mature_threshold_cluster1_q25=mature_thr )
}))
write.csv(sens,file.path(tab_dir,"12_transition_threshold_sensitivity.csv"), row.names=FALSE)
if (high_n >= 30) {
  obj$transition_marker_group <- ifelse( obj$shared_transition_3dpi_highconf, "3dpi_transition_highconf", "other_microglia" )
  Idents(obj) <- "transition_marker_group"
  mk <- FindMarkers( obj, ident.1="3dpi_transition_highconf", ident.2="other_microglia", assay="RNA", only.pos=TRUE, min.pct=0.10, logfc.threshold=0.25, test.use="wilcox", verbose=FALSE )
  mk$gene <- rownames(mk)
  mk <- mk[order(mk$p_val_adj,-mk$avg_log2FC),,drop=FALSE]
  write.csv(mk,file.path(tab_dir,"13_highconf_transition_marker_audit.csv"), row.names=FALSE)
  write.csv(head(mk,30),file.path(tab_dir,"14_top30_highconf_transition_markers.csv"), row.names=FALSE)
}
Idents(obj) <- "step3A_cluster_r04"
p_scores <- FeaturePlot( obj, features=c(homeo_col,act_col,early_col,mature_col), reduction="umap.microglia.step3A", ncol=2, order=TRUE )
ggsave(file.path(fig_dir,"01_transition_program_score_UMAPs.pdf"), p_scores,width=13,height=10)
p_status <- DimPlot( obj, reduction="umap.microglia.step3A", group.by="transition_status" ) + ggtitle("Activated-to-DAM transition status")
ggsave(file.path(fig_dir,"02_transition_status_UMAP.pdf"), p_status,width=8,height=7)
cells3 <- colnames(obj)[as.character(obj$timepoint)=="3 dpi"]
obj3 <- subset(obj,cells=cells3)
obj3$transition_highlight <- ifelse( obj3$shared_transition_3dpi_highconf, "Shared_transition_highconf", "Other_3dpi_microglia" )
p3a <- DimPlot( obj3,reduction="umap.microglia.step3A", group.by="transition_highlight" ) + ggtitle("3 dpi | high-confidence activated/Early-DAM transition")
p3b <- DimPlot( obj3,reduction="umap.microglia.step3A", group.by="step3A_cluster_r04", label=TRUE,repel=TRUE ) + ggtitle("3 dpi | Step 3 cluster context")
ggsave(file.path(fig_dir,"03_3dpi_shared_transition_UMAP.pdf"), p3a|p3b,width=14,height=6)
scatter_df <- obj[[]]
p_scatter <- ggplot(scatter_df,aes( x=.data[[act_col]],y=.data[[early_col]] )) + geom_point(aes(shape=as.character(timepoint)),alpha=0.25,size=0.7) + geom_vline(xintercept=act_thr,linetype="dashed") + geom_hline(yintercept=early_thr,linetype="dashed") + theme_bw() + labs( title="Activated vs Early-DAM programs", subtitle="Dashed lines = q90 of cluster-3 homeostatic reference", x="Activated-injury AUCell AUC", y="Early-DAM AUCell AUC", shape="Timepoint" )
ggsave(file.path(fig_dir,"04_Activated_vs_EarlyDAM_scatter.pdf"), p_scatter,width=8,height=7)
time_levels <- c("Uninjured","1 dpi","3 dpi","7 dpi")
pd <- obj[[]]
pd$timepoint_plot <- factor(as.character(pd$timepoint),levels=time_levels)
long <- bind_rows( data.frame(timepoint=pd$timepoint_plot,program="Homeostatic",score=pd[[homeo_col]]), data.frame(timepoint=pd$timepoint_plot,program="Activated_injury",score=pd[[act_col]]), data.frame(timepoint=pd$timepoint_plot,program="Early_DAM",score=pd[[early_col]]), data.frame(timepoint=pd$timepoint_plot,program="Mature_DAM",score=pd[[mature_col]]) )
p_v <- ggplot(long,aes(x=timepoint,y=score)) + geom_violin(scale="width",trim=TRUE) + geom_boxplot(width=0.12,outlier.shape=NA) + facet_wrap(~program,scales="free_y",ncol=2) + theme_bw() + labs(title="Transition-stage program activity across SCI time", x="SCI timepoint",y="AUCell AUC")
ggsave(file.path(fig_dir,"05_transition_programs_by_timepoint.pdf"), p_v,width=12,height=8)
stp <- status_time
stp$timepoint <- factor(as.character(stp$timepoint),levels=time_levels)
p_st <- ggplot(stp,aes( x=timepoint,y=proportion_within_timepoint,fill=transition_status )) + geom_col() + theme_bw() + labs(title="Transition-state composition across SCI time", x="SCI timepoint",y="Proportion within timepoint", fill="Transition status")
ggsave(file.path(fig_dir,"06_transition_status_proportion_by_timepoint.pdf"), p_st,width=9,height=6)
mat <- xtabs(Freq ~ cluster + transition_status,data=status_cluster)
prop <- prop.table(mat,margin=1)
pdf(file.path(fig_dir,"07_cluster_transition_status_heatmap.pdf"), width=10,height=7)
pheatmap(prop,cluster_rows=FALSE,cluster_cols=FALSE, main="Step 3 clusters x activated-DAM transition status | row proportion")
dev.off()
if (nrow(three)>0) {
  p_s <- ggplot(three,aes(x=sample_id,y=fraction_transition_highconf)) + geom_col() + theme_bw() + labs( title="High-confidence shared transition in 3 dpi biological replicates", x="3 dpi biological sample", y="Fraction high-confidence transition" ) + theme(axis.text.x=element_text(angle=45,hjust=1))
  ggsave(file.path(fig_dir,"08_3dpi_transition_by_sample.pdf"), p_s,width=8,height=6)
}
tm <- as.matrix(time_summary[,c( "Homeostatic_median","Activated_median", "EarlyDAM_median","MatureDAM_median" )])
rownames(tm) <- time_summary$timepoint_chr
tm <- t(tm)
row_z <- function(x) {
  s <- sd(x,na.rm=TRUE)
  if (!is.finite(s) || s==0) return(rep(0,length(x)))
  (x-mean(x,na.rm=TRUE))/s
}
tmz <- t(apply(tm,1,row_z))
pdf(file.path(fig_dir,"09_transition_program_time_heatmap.pdf"), width=8,height=5)
pheatmap(tmz,cluster_rows=FALSE,cluster_cols=FALSE, main="Median transition-program activity across SCI time | row Z-score")
dev.off()
output_rds <- file.path( out_dir, "GSE162610_step3C_transition_mapped.rds" )
saveRDS(obj,output_rds)
report_file <- file.path(out_dir,"15_STEP3C_transition_mapping_report.txt")
sink(report_file)
cat("GSE162610 Step 3C | activated-DAM transition mapping\n")
cat("====================================================\n\n")
cat("Input RDS:\n",input_rds,"\n\n",sep="")
cat("Input cells: ",ncol(obj),"\n",sep="")
cat("\nBiological references:\n")
cat("cluster 3 = canonical homeostatic reference\n")
cat("cluster 1 = established DAM-like/phagocytic reference\n")
cat("\nProgram availability:\n")
print(avail)
cat("\nMain data-driven thresholds:\n")
print(thr_tbl)
cat("\nTransition definition:\n")
cat("Activated_high > q90 of cluster3 Activated score\n")
cat("EarlyDAM_high > q90 of cluster3 Early-DAM score\n")
cat("MatureDAM_established >= q25 of cluster1 mature-DAM score\n")
cat("Transition = Activated_high + EarlyDAM_high + NOT MatureDAM_established\n")
cat("High-confidence shared transition = Transition + 3 dpi + not clusters 6/8/9\n")
cat("\nTransition status counts:\n")
print(status_counts)
cat("\nTransition status by timepoint:\n")
print(status_time)
cat("\nHigh-confidence 3 dpi transition cells: ",high_n,"\n",sep="")
cat("\nHigh-confidence 3 dpi transition by cluster:\n")
print(high_cluster)
cat("\nHigh-confidence 3 dpi transition by sample:\n")
print(high_sample)
cat("\n3 dpi biological replicate validation:\n")
print(three)
cat("\nScore summary by timepoint:\n")
print(time_summary)
cat("\nCluster score summary:\n")
print(cluster_summary)
cat("\nThreshold sensitivity q85/q90/q95:\n")
print(sens)
cat("\nIMPORTANT\n")
cat("1. Step3C does NOT define repair/scar.\n")
cat("2. Step3C does NOT run Monocle.\n")
cat("3. Shared precursor should preferentially appear around 3 dpi.\n")
cat("4. Stronger evidence requires both 3 dpi biological replicates.\n")
cat("5. Stronger evidence requires robustness across q85/q90/q95.\n")
cat("6. Clusters 6/8/9 remain in object but are excluded from high-confidence precursor call.\n")
cat("7. Full cluster/time DEG and GO/KEGG/GSEA belong to Step3D.\n")
cat("\nSaved object:\n",output_rds,"\n",sep="")
cat("\nSTEP 3C TECHNICAL COMPLETE = TRUE\n")
sink()
handoff <- file.path(out_dir,"16_STEP3D_handoff.txt")
sink(handoff)
cat("GSE162610 | Step 3D handoff\n")
cat("============================\n\n")
cat("Formal input:\n",output_rds,"\n\n",sep="")
cat("Step3C completed:\n")
cat("- Homeostatic / Activated / Early-DAM / Mature-DAM per-cell mapping\n")
cat("- reference-based thresholds\n")
cat("- 3 dpi shared-transition identification\n")
cat("- biological-replicate validation\n")
cat("- q85/q90/q95 sensitivity analysis\n\n")
cat("Step3D goal:\n")
cat("cluster/time DEG + GO/KEGG/GSEA confirmation of\n")
cat("1 dpi activation -> 3 dpi shared transition -> 7 dpi mature DAM remodeling.\n")
cat("Do NOT define repair/scar until Step4.\n")
sink()
message("====================================================")
message("STEP 3C COMPLETE")
message("Results: ",out_dir)
message("High-confidence 3 dpi transition candidates: ",high_n)
message("")
message("Key outputs:")
message("1) figures/01_transition_program_score_UMAPs.pdf")
message("2) figures/03_3dpi_shared_transition_UMAP.pdf")
message("3) figures/05_transition_programs_by_timepoint.pdf")
message("4) figures/06_transition_status_proportion_by_timepoint.pdf")
message("5) tables/10_3dpi_transition_by_biological_sample.csv")
message("6) tables/12_transition_threshold_sensitivity.csv")
message("7) 15_STEP3C_transition_mapping_report.txt")
message("====================================================")
