# =============================================================================
# 79_step6G_GSE162610_R2MC6_trajectory_topology_sensitivity_v2_fix.R
# GLIA | Reviewer 2 | Major Comment 6
#
# v2 fix: reset Monocle3 S4 slots with S4Vectors::SimpleList(), not list().
# Objective: independently RELEARN Monocle3 principal graphs under stratified
# subsampling and graph parameters, then audit two 7-dpi sample-specific arm
# recoveries and 1-dpi high-pseudotime cells. Frozen original Step5BR2 labels
# are used ONLY for POST-HOC recovery, NEVER for choosing the new branchpoint.
#
# Input (READ ONLY; NOT a ZIP):
#   xjz/04_results/step5BR2_GSE162610_single_root_branch_audit/
#       GSE162610_step5BR2_single_root_cds.rds
#       GSE162610_step5BR2_single_root_seurat.rds (root-cell flags/labels)
# Optional context, read only:
#   xjz/04_results/step6FRH_GSE162610_Harmony_matched_trajectory_sensitivity/
#       31N_source_matched_Marsh_paired_trajectory/
#       31N_STEP6FRH_SOURCE_MATCHED_MARSH_PAIRED_FINAL_REPORT.txt
# Output (all NEW files):
#   xjz/04_results/step6G_GSE162610_trajectory_topology_sensitivity/
#       00_run_log.txt, 00_sessionInfo.txt, tables/, figures/, scenarios/,
#       10_STEP6G_FINAL_REPORT.txt
#
# Analysis deliberately keeps the original 2D UMAP FIXED; this tests stability
# of GRAPH LEARNING conditional on that embedding, not full-pipeline robustness.
# It does NOT prove cell fate transitions, longitudinal lineage, or causality.
# No previous scripts/results are overwritten. No uploaded ZIP is accessed.
# =============================================================================

options(stringsAsFactors = FALSE)
PROJECT_ROOT <- Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
RES <- file.path(PROJECT_ROOT, "04_results")
IN_CDS <- file.path(RES, "step5BR2_GSE162610_single_root_branch_audit",
                    "GSE162610_step5BR2_single_root_cds.rds")
IN_SEURAT <- file.path(RES, "step5BR2_GSE162610_single_root_branch_audit",
                       "GSE162610_step5BR2_single_root_seurat.rds")
PRIOR_MARSH <- file.path(RES,
   "step6FRH_GSE162610_Harmony_matched_trajectory_sensitivity",
   "31N_source_matched_Marsh_paired_trajectory",
   "31N_STEP6FRH_SOURCE_MATCHED_MARSH_PAIRED_FINAL_REPORT.txt")
OUT <- file.path(RES, "step6G_GSE162610_trajectory_topology_sensitivity")
TBL <- file.path(OUT, "tables")
FIG <- file.path(OUT, "figures")
CASES <- file.path(OUT, "scenarios")
for (d in c(OUT, TBL, FIG, CASES)) dir.create(d, recursive=TRUE, showWarnings=FALSE)

# Set TRUE if intentionally discarding/recomputing previously completed scenarios.
FORCE_RERUN <- FALSE
# Set FALSE for a quick baseline-only audit before the full sensitivity run.
RUN_SENSITIVITY <- TRUE
# Predefined descriptive threshold; not a biological-replicate significance test.
RECALL_THRESHOLD <- 0.95
HIGH_1DPI_PT <- 15
ORIGINAL_SEED <- 20260731L

# Parameters preserve the original full-data controls (k=30, ratio=0.8), while
# varying k and graph pruning; subsamples are stratified by sample x timepoint.
scenarios <- data.frame(
  scenario_id = c("00_full_original_parameters",
    "01_full_k20", "02_full_k40", "03_full_ratio065", "04_full_ratio095",
    "05_sub80_seed1", "06_sub80_seed2", "07_sub80_seed3",
    "08_sub65_seed1", "09_sub65_seed2"),
  fraction = c(1,1,1,1,1, .80,.80,.80,.65,.65),
  seed = c(ORIGINAL_SEED, ORIGINAL_SEED, ORIGINAL_SEED,
           ORIGINAL_SEED, ORIGINAL_SEED, 20260801L,20260802L,
           20260803L,20260811L,20260812L),
  cluster_k = c(30,20,40,30,30,30,30,30,30,30),
  euclidean_distance_ratio = c(.8,.8,.8,.65,.95,.8,.8,.8,.8,.8),
  stringsAsFactors=FALSE
)
if (!RUN_SENSITIVITY) scenarios <- scenarios[1,,drop=FALSE]

write_csv <- function(x, path) {
  utils::write.csv(x, file=path, row.names=FALSE, na="")
}
write_gz_csv <- function(x, path) {
  co <- gzfile(path, "wt", compression=6)
  on.exit(close(co))
  utils::write.csv(x, file=co, row.names=FALSE, na="")
}
log_file <- file.path(OUT, "00_run_log.txt")
log_msg <- function(...) {
  msg <- paste0(..., collapse="")
  line <- paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "|", msg)
  cat(line,"\n")
  cat(line,"\n",file=log_file,append=TRUE)
}
for (input in c(IN_CDS, IN_SEURAT))
  if (!file.exists(input)) stop("找不到原始 Step5BR2 对象: ", input)
for (pkg in c("monocle3", "igraph", "S4Vectors",
              "SingleCellExperiment", "SummarizedExperiment"))
  if (!requireNamespace(pkg, quietly=TRUE))
    stop("缺少R包: ", pkg, ". 请使用原项目的R环境运行。")
writeLines(capture.output(sessionInfo()), file.path(OUT,"00_sessionInfo.txt"))
write_csv(scenarios,file.path(TBL,"01_scenario_design.csv"))
log_msg("Reviewer 2 MC6 topology sensitivity | input ", IN_CDS)
log_msg("Reading Step5BR2 Seurat METADATA first (read-only).")
if(!requireNamespace("SeuratObject",quietly=TRUE))
  stop("缺少SeuratObject包；无法读取冻结的root/branch标签")
ref_seurat <- readRDS(IN_SEURAT)
ref_seurat_meta <- ref_seurat[[]]
rm(ref_seurat); gc(verbose=FALSE)
log_msg("Reading FROZEN Step5BR2 Monocle3 CDS (read-only).")
ref <- readRDS(IN_CDS)
cell_ids <- colnames(ref)
if (length(cell_ids)!=14005L || anyDuplicated(cell_ids))
  stop("冻结的轨迹细胞数/唯一性与14005细胞不符，请检查输入。")
meta <- as.data.frame(SummarizedExperiment::colData(ref))
if(!all(cell_ids %in% rownames(meta)) ||
   !all(cell_ids %in% rownames(ref_seurat_meta)))
   stop("CDS/Seurat 的barcode不能一一对齐，已停止以防错位")
meta <- meta[cell_ids,,drop=FALSE]
# Original Step5BR root flag was added to the Seurat object AFTER constructing
# the Step5BR CDS; do not assume it exists in colData of the CDS.
meta$step5BR_primary_root <- as.logical(
      ref_seurat_meta[cell_ids,"step5BR_primary_root"])
# Match the original frozen branch labels against saved Seurat labels.
if(!"step5BR2_branch_arm" %in% names(ref_seurat_meta))
  stop("参考Seurat缺少冻结的branch arm标签")
if("step5BR2_branch_arm" %in% names(meta)) {
  a <- as.character(meta$step5BR2_branch_arm)
  b <- as.character(ref_seurat_meta[cell_ids,"step5BR2_branch_arm"])
  if(!identical(a,b)) stop("参考CDS和Seurat的分支标签不一致")
}
meta$step5BR2_branch_arm <- as.character(
      ref_seurat_meta[cell_ids,"step5BR2_branch_arm"])
rm(ref_seurat_meta);gc(verbose=FALSE)
needed <- c("sample_id", "timepoint", "step5BR_primary_root",
            "shared_transition_3dpi_highconf", "step5BR2_branch_arm",
            "step3A_cluster_r04")
miss <- setdiff(needed, names(meta))
if (length(miss)) stop("Step5BR2 CDS 缺少字段: ", paste(miss,collapse=", "))
sample_id <- as.character(meta$sample_id)
timepoint <- as.character(meta$timepoint)
root_flag <- as.logical(meta$step5BR_primary_root)
if (anyNA(root_flag) || sum(root_flag)!=225L)
  stop("原始 root flag 不是预期225个细胞，禁止自动推测根细胞。")
transition <- as.logical(meta$shared_transition_3dpi_highconf)
transition[is.na(transition)] <- FALSE
old_arm <- as.character(meta$step5BR2_branch_arm)
if (!all(timepoint %in% c("1 dpi","3 dpi","7 dpi")))
  stop("timepoint 标签异常；请先核查源数据。")
orig_7 <- sum(timepoint=="7 dpi" & old_arm %in% c("Arm_1","Arm_2"))
if (orig_7 != 2541L) stop("Step5BR2 7dpi两条原始分支应为2541细胞，实际: ",orig_7)
if (length(unique(sample_id[timepoint=="7 dpi"]))!=2L)
  stop("预期两个7dpi生物学样本，当前不符。")
ref_umap <- as.matrix(SingleCellExperiment::reducedDim(ref, "UMAP"))
if (is.null(ref_umap) || ncol(ref_umap)<2L) stop("参考CDS缺少UMAP")
if (is.null(rownames(ref_umap))) rownames(ref_umap) <- cell_ids
ref_umap <- ref_umap[cell_ids,1:2,drop=FALSE]
ref_graph <- monocle3::principal_graph(ref)[["UMAP"]]
if (is.null(ref_graph)) stop("参考CDS中缺少Monocle principal graph")
ref_pt <- tryCatch(monocle3::pseudotime(ref,reduction_method="UMAP"),error=function(e)NULL)
if (is.null(ref_pt)) stop("Step5BR2没有原始pseudotime")
if (is.null(names(ref_pt))) names(ref_pt) <- cell_ids
ref_pt <- as.numeric(ref_pt[cell_ids])
if (any(!is.finite(ref_pt))) stop("原始pseudotime含非有限值")

input_audit <- data.frame(
  key=c("analysis_script_version","input_file","input_reference_seurat","n_cells","n_primary_root_cells","n_7dpi_old_arm_cells",
        "ref_graph_nodes","ref_graph_edges","ref_graph_components",
        "n_7dpi_biological_samples","embedding_scope","original_parameters",
        "prior_Marsh_report_exists"),
  value=c("Step6G_v2_SimpleList_S4_fix",IN_CDS,IN_SEURAT,length(cell_ids),sum(root_flag),orig_7,
          igraph::vcount(ref_graph),igraph::ecount(ref_graph),
          igraph::components(ref_graph)$no,2,
          "Frozen Step5BR2 UMAP; graph relearned in every scenario",
          "Leiden; k=30; ratio=0.8; partition=FALSE; close_loop=FALSE",
          file.exists(PRIOR_MARSH)))
write_csv(input_audit,file.path(TBL,"00_input_and_parameter_audit.csv"))

# Explicitly audit the Figure 2D 1-dpi-high-pseudotime issue in the ORIGINAL
# saved graph, rather than explaining it from a graph reconstructed later.
ref_table <- data.frame(cell=cell_ids,sample_id=sample_id,timepoint=timepoint,
  original_arm=old_arm,step3A_cluster=as.character(meta$step3A_cluster_r04),
  primary_root_cell=root_flag,shared_transition=transition,
  original_pseudotime=ref_pt,UMAP_1=ref_umap[,1],UMAP_2=ref_umap[,2])
write_gz_csv(ref_table, file.path(TBL,"02_reference_cell_pseudotime_and_labels.csv.gz"))
high_1dpi <- ref_table[ref_table$timepoint=="1 dpi" &
                           ref_table$original_pseudotime >= HIGH_1DPI_PT,,drop=FALSE]
write_gz_csv(high_1dpi,file.path(TBL,"03_reference_1dpi_high_pseudotime_cells.csv.gz"))
ref_groups <- split(ref_table,interaction(ref_table$sample_id,
               ref_table$timepoint,drop=TRUE))
ref_summary <- do.call(rbind,lapply(ref_groups,function(d) data.frame(
  sample_id=d$sample_id[1], timepoint=d$timepoint[1],n_cells=nrow(d),
  median_pt=stats::median(d$original_pseudotime),
  q75_pt=as.numeric(stats::quantile(d$original_pseudotime,.75)),
  n_high_ge15=sum(d$original_pseudotime>=HIGH_1DPI_PT),
  fraction_high_ge15=mean(d$original_pseudotime>=HIGH_1DPI_PT),
  n_primary_root=sum(d$primary_root_cell))))
write_csv(ref_summary,file.path(TBL,"04_reference_pseudotime_by_sample.csv"))
high_clu <- as.data.frame(table(sample_id=high_1dpi$sample_id,
       step3A_cluster=high_1dpi$step3A_cluster,
       primary_root=high_1dpi$primary_root_cell))
high_clu <- high_clu[high_clu$Freq>0,,drop=FALSE]
write_csv(high_clu,file.path(TBL,"05_reference_high_1dpi_sample_cluster_root.csv"))

pdf(file.path(FIG,"01_reference_1dpi_high_pseudotime_audit.pdf"),width=12,height=5.4)
par(mfrow=c(1,2),mar=c(4,4,3,1))
cols <- ifelse(timepoint=="1 dpi" & ref_pt>=HIGH_1DPI_PT,"#D64F3C",
          ifelse(timepoint=="1 dpi","#3E83AD","#BDBDBD"))
plot(ref_umap[,1],ref_umap[,2],col=cols,pch=16,cex=.2,
     xlab="UMAP 1",ylab="UMAP 2",main="Original frozen graph: high-PT 1 dpi")
legend("topleft",c("1 dpi, PT >=15","1 dpi, PT <15","3/7 dpi"),
       pch=16,col=c("#D64F3C","#3E83AD","#BDBDBD"),bty="n",cex=.8)
boxplot(split(ref_pt,timepoint),ylab="Original pseudotime",
        main="Original pseudotime by actual timepoint")
abline(h=HIGH_1DPI_PT,lty=2,col="gray40")
dev.off()

# Monocle projection returns graph-local vertex ID; safely map back by BARCODE.
closest_graph_vertices <- function(cds,g) {
  v <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  if (is.null(v)) stop("New graph has no cell-to-closest-vertex map")
  if (is.matrix(v) || is.data.frame(v)) {
    nn <- rownames(v); v <- v[,1L]; names(v) <- nn
  }
  if (is.null(names(v))) {
    if (length(v)!=ncol(cds)) stop("Unlabeled closest-vertex vector length mismatch")
    names(v) <- colnames(cds)
  }
  cc <- colnames(cds)
  if (!all(cc %in% names(v))) stop("Cell barcodes missing from closest-vertex map")
  z <- as.character(v[cc]); gn <- igraph::V(g)$name
  if (all(z %in% gn)) return(stats::setNames(z,cc))
  if (all(paste0("Y_",z) %in% gn)) return(stats::setNames(paste0("Y_",z),cc))
  zz <- suppressWarnings(as.integer(z))
  if (all(!is.na(zz) & zz>=1L & zz<=length(gn)))
    return(stats::setNames(gn[zz],cc))
  stop("Cannot map projection IDs to named Monocle graph vertices")
}

sample_stratified_cells <- function(frac,seed) {
  if (frac>=1) return(cell_ids)
  set.seed(seed)
  strat <- split(cell_ids,interaction(sample_id,timepoint,drop=TRUE))
  sort(unlist(lapply(strat,function(x) sample(x,size=max(1L,round(length(x)*frac)),
                                          replace=FALSE)),use.names=FALSE))
}

# New branchpoint ranking exactly mirrors original Step5BR2/Step6FRH objective:
#  two downstream cell-supported arms, each with >=100 total and 7dpi cells
#  on full data; proportional threshold in subsamples, then by 3dpi high-confidence
#  transition proximity within graph radius 2, min 7dpi arm size, root distance.
choose_branchpoint <- function(g, closest,cc,case_time,case_trans,root_anchor,frac) {
  root_distance <- as.numeric(igraph::distances(g,v=root_anchor,
                                  to=igraph::V(g))[1,])
  names(root_distance) <- igraph::V(g)$name
  bp <- names(which(igraph::degree(g)>=3))
  rows <- list(); caches <- list()
  nmin <- max(40L, ceiling(100*frac))
  for (v in bp) {
    if (identical(v,root_anchor)) next
    gdel <- igraph::delete_vertices(g,v)
    membership <- igraph::components(gdel)$membership
    names(membership) <- igraph::V(gdel)$name
    root_side <- membership[[root_anchor]]
    if (is.null(root_side) || is.na(root_side)) next
    member <- unname(membership[unname(closest[cc])]); names(member) <- cc
    downstream <- unique(member[!is.na(member) & member!=root_side])
    downstream <- downstream[!is.na(downstream)]
    if (!length(downstream)) next
    arms <- do.call(rbind,lapply(downstream,function(z) data.frame(
       component=as.integer(z),
       n_cells=sum(member==z,na.rm=TRUE),
       n_7dpi=sum(member==z & case_time=="7 dpi",na.rm=TRUE))))
    arms <- arms[order(-arms$n_7dpi,-arms$n_cells),,drop=FALSE]
    major <- arms[arms$n_cells>=nmin & arms$n_7dpi>=nmin,,drop=FALSE]
    distv <- as.numeric(igraph::distances(g,v=v,to=igraph::V(g))[1,])
    nearby <- igraph::V(g)$name[distv<=2]
    nnear <- sum(unname(closest[cc]) %in% nearby &
                   case_time=="3 dpi" & case_trans)
    rows[[length(rows)+1L]] <- data.frame(
       branchpoint=v,degree=as.integer(igraph::degree(g,v=v)),
       downstream_supported=nrow(arms),downstream_major=nrow(major),
       min_top2_7dpi=if(nrow(major)>=2) min(major$n_7dpi[1:2]) else 0L,
       near_3dpi_highconf=nnear,dist_from_root=root_distance[[v]],
       candidate_bifurcation=nrow(major)>=2,
       min_support_threshold=nmin)
    caches[[v]] <- list(member=member,root_side=root_side,arms=arms,major=major)
  }
  blank <- data.frame(branchpoint=character(),degree=integer(),
    downstream_supported=integer(),downstream_major=integer(),
    min_top2_7dpi=integer(),near_3dpi_highconf=integer(),
    dist_from_root=numeric(),candidate_bifurcation=logical(),
    min_support_threshold=integer())
  if (!length(rows)) return(list(audit=blank,selected=NA_character_,cache=NULL))
  dat <- do.call(rbind,rows)
  dat <- dat[order(-as.integer(dat$candidate_bifurcation),
                    -dat$near_3dpi_highconf,-dat$min_top2_7dpi,
                    dat$dist_from_root,dat$branchpoint),,drop=FALSE]
  sel <- if(isTRUE(dat$candidate_bifurcation[1])) dat$branchpoint[1] else NA_character_
  list(audit=dat,selected=sel,
       cache=if(!is.na(sel)) caches[[sel]] else NULL)
}

scenario_dir <- function(id) file.path(CASES,id)
case_cols <- c("scenario_id","fraction","seed","cluster_k","euclidean_distance_ratio")
write_csv(scenarios, file.path(TBL,"01_scenario_design.csv"))

run_case <- function(param) {
  id <- param$scenario_id
  dd <- scenario_dir(id)
  dir.create(dd,recursive=TRUE,showWarnings=FALSE)
  cfg_path <- file.path(dd,"00_case_parameters.csv")
  success_path <- file.path(dd,"SUCCESS.flag")
  current <- as.data.frame(param[,case_cols,drop=FALSE])
  if (file.exists(success_path) && !FORCE_RERUN) {
    prev <- tryCatch(read.csv(cfg_path,stringsAsFactors=FALSE),error=function(e)NULL)
    same_cfg <- !is.null(prev) && identical(names(prev),names(current)) &&
      nrow(prev)==1L && all(vapply(names(current), function(n)
        identical(as.character(prev[[n]][1]),as.character(current[[n]][1])),
        logical(1)))
    if (!same_cfg)
      stop("Existing result has different/unknown parameters: ",dd,
           " ; set FORCE_RERUN=TRUE to replace it")
    log_msg("SKIP completed scenario: ",id)
    return(invisible(TRUE))
  }
  if(file.exists(success_path)) unlink(success_path)
  if(file.exists(file.path(dd,"FAILURE.txt"))) unlink(file.path(dd,"FAILURE.txt"))
  write_csv(current,cfg_path)
  log_msg("RUN ",id," | fraction=",param$fraction," k=",param$cluster_k,
          " ratio=",param$euclidean_distance_ratio)
  stage <- "subsample_and_clear_graph_slots"
  tryCatch({
    keep <- sample_stratified_cells(param$fraction,param$seed)
    obj <- ref[,keep]  # deep subset of ORIGINAL dataset; NEVER alter ref on disk
    # All three slots in monocle3::cell_data_set are S4Vectors SimpleList.
    # Base list() assignment throws a class-validation error before graph learning.
    obj@principal_graph <- S4Vectors::SimpleList()
    obj@principal_graph_aux <- S4Vectors::SimpleList()
    obj@clusters <- S4Vectors::SimpleList()
    is_valid <- methods::validObject(obj, test=TRUE)
    if (!identical(is_valid, TRUE))
      stop("Cell-data-set slot validation failed: ", paste(is_valid, collapse="; "))
    orig_emb <- ref_umap[keep,,drop=FALSE]
    now <- as.matrix(SingleCellExperiment::reducedDim(obj,"UMAP"))
    if(is.null(rownames(now))) rownames(now)<-colnames(obj)
    if(max(abs(now[keep,1:2,drop=FALSE]-orig_emb))>1e-8)
       stop("Frozen original UMAP changed during subsetting")
    stage <- "monocle3_cluster_cells"
    set.seed(as.integer(param$seed))
    obj <- monocle3::cluster_cells(obj,reduction_method="UMAP",
       k=as.integer(param$cluster_k),cluster_method="leiden",
       random_seed=as.integer(param$seed),verbose=FALSE)
    stage <- "monocle3_learn_graph"
    obj <- monocle3::learn_graph(obj,use_partition=FALSE,close_loop=FALSE,
      learn_graph_control=list(euclidean_distance_ratio=
                                  param$euclidean_distance_ratio))
    graph <- monocle3::principal_graph(obj)[["UMAP"]]
    if(is.null(graph) || igraph::components(graph)$no!=1L)
       stop("Relearned principal graph is absent/disconnected")
    if(max(abs(as.matrix(SingleCellExperiment::reducedDim(obj,"UMAP"))[
                         keep,1:2,drop=FALSE]-orig_emb))>1e-8)
       stop("UMAP changed during graph learning")
    stage <- "cell_projection_and_root_selection"
    closest <- closest_graph_vertices(obj,graph)
    cmeta <- meta[keep,,drop=FALSE]
    ctime <- as.character(cmeta$timepoint)
    csample <- as.character(cmeta$sample_id)
    ctrans <- as.logical(cmeta$shared_transition_3dpi_highconf)
    ctrans[is.na(ctrans)] <- FALSE
    roots <- keep[as.logical(cmeta$step5BR_primary_root)]
    if(length(roots)<20L) stop("Too few retained primary roots")
    rv <- sort(table(closest[roots]),decreasing=TRUE)
    root_anchor <- names(rv)[1L]
    stage <- "monocle3_order_cells"
    obj <- monocle3::order_cells(obj,reduction_method="UMAP",
                                root_pr_nodes=root_anchor)
    pt <- monocle3::pseudotime(obj,reduction_method="UMAP")
    if (is.null(names(pt))) names(pt) <- colnames(obj)
    pt <- as.numeric(pt[keep])
    if(any(!is.finite(pt))) stop("Nonfinite Monocle pseudotime")
    stage <- "objective_branchpoint_selection"
    rule <- choose_branchpoint(graph,closest,keep,ctime,ctrans,
                               root_anchor,param$fraction)
    new_arm <- setNames(rep("Unassigned",length(keep)),keep)
    if(!is.na(rule$selected)) {
      arm <- rule$cache
      new_arm[!is.na(arm$member) & arm$member==arm$root_side] <- "Pre_branch_root_side"
      top <- head(arm$major$component,2)
      new_arm[!is.na(arm$member) & arm$member==top[1]] <- "New_Arm_1"
      new_arm[!is.na(arm$member) & arm$member==top[2]] <- "New_Arm_2"
      other <- !is.na(arm$member) & arm$member!=arm$root_side &
               new_arm=="Unassigned"
      new_arm[other] <- "Other_downstream"
    }
    old <- as.character(cmeta$step5BR2_branch_arm)
    eval <- ctime=="7 dpi" & old %in% c("Arm_1","Arm_2")
    # Global flip only, never per sample; labels are used only AFTER graph selection.
    aligned <- new_arm
    if (any(eval)) {
      same <- sum(eval & ((old=="Arm_1" & new_arm=="New_Arm_1") |
                          (old=="Arm_2" & new_arm=="New_Arm_2")))
      flip <- sum(eval & ((old=="Arm_1" & new_arm=="New_Arm_2") |
                          (old=="Arm_2" & new_arm=="New_Arm_1")))
      if(flip>same) {
        aligned[new_arm=="New_Arm_1"] <- "Mapped_Arm_2"
        aligned[new_arm=="New_Arm_2"] <- "Mapped_Arm_1"
        arm_orientation <- "global_flip"
      } else {
        aligned[new_arm=="New_Arm_1"] <- "Mapped_Arm_1"
        aligned[new_arm=="New_Arm_2"] <- "Mapped_Arm_2"
        arm_orientation <- "direct"
      }
    } else arm_orientation <- "no_7dpi_reference"
    stage <- "save_correspondence_and_sample_audit"
    corr <- data.frame(cell=keep,sample_id=csample,timepoint=ctime,
      original_arm=old,relearned_arm=unname(new_arm[keep]),
      globally_aligned_arm=unname(aligned[keep]),
      new_pseudotime=pt,original_pseudotime=ref_pt[match(keep,cell_ids)],
      UMAP_1=orig_emb[,1],UMAP_2=orig_emb[,2],
      primary_root=as.logical(cmeta$step5BR_primary_root),
      stringsAsFactors=FALSE)
    write_gz_csv(corr,file.path(dd,"03_cell_correspondence.csv.gz"))
    rule$audit$scenario_id <- rep(id,nrow(rule$audit))
    write_csv(rule$audit,file.path(dd,"02_branchpoint_candidates.csv"))
    # Per-sample recall is conditional on ORIGINAL labeled arm cells RETAINED
    # in that subsample; report retention denominator separately.
    sample_levels <- sort(unique(sample_id[timepoint=="7 dpi"]))
    rr <- do.call(rbind,lapply(sample_levels,function(s) do.call(rbind,
      lapply(c("Arm_1","Arm_2"),function(a) {
        orig_n <- sum(timepoint=="7 dpi" & sample_id==s & old_arm==a)
        take <- ctime=="7 dpi" & csample==s & old==a
        nret <- sum(take)
        nmatch <- sum(take & aligned==paste0("Mapped_",a))
        data.frame(scenario_id=id,sample_id=s,reference_arm=a,
          n_original=orig_n,n_retained=nret,retained_fraction=nret/orig_n,
          n_recovered=nmatch,recall=if(nret) nmatch/nret else NA_real_)
      })) ))
    write_csv(rr,file.path(dd,"04_7dpi_sample_arm_recovery.csv"))
    # Independently verify BOTH new arms contain adequate 7dpi cells in BOTH samples.
    armcounts <- as.data.frame(table(sample_id=csample[ctime=="7 dpi"],
                                     relearned_arm=new_arm[ctime=="7 dpi"]),
                              stringsAsFactors=FALSE)
    armcounts <- armcounts[armcounts$relearned_arm %in% c("New_Arm_1","New_Arm_2"),,
                           drop=FALSE]
    write_csv(armcounts,file.path(dd,"05_new_arm_counts_by_7dpi_sample.csv"))
    independent_sample_support <- (nrow(armcounts)==4L &&
                                    all(armcounts$Freq>=50L))
    tp <- do.call(rbind,lapply(c("1 dpi","3 dpi","7 dpi"),function(t) {
      values <- pt[ctime==t]
      data.frame(scenario_id=id,timepoint=t,n_cells=length(values),
         median_pt=median(values),q25=as.numeric(quantile(values,.25)),
         q75=as.numeric(quantile(values,.75)),
         fraction_in_within_case_top_quartile=mean(values >=
                                as.numeric(quantile(pt,.75))),
         n_pseudotime_ge15=sum(values>=HIGH_1DPI_PT),
         fraction_pseudotime_ge15=mean(values>=HIGH_1DPI_PT))
    }))
    write_csv(tp,file.path(dd,"06_pseudotime_by_timepoint.csv"))
    sd <- data.frame(scenario_id=id,fraction=param$fraction,seed=param$seed,
      cluster_k=param$cluster_k,
      euclidean_distance_ratio=param$euclidean_distance_ratio,
      n_cells=length(keep),n_root_cells=length(roots),root_anchor=root_anchor,
      graph_nodes=igraph::vcount(graph),graph_edges=igraph::ecount(graph),
      graph_size_matches_reference=(igraph::vcount(graph)==igraph::vcount(ref_graph) &&
                                    igraph::ecount(graph)==igraph::ecount(ref_graph)),
      graph_components=igraph::components(graph)$no,
      n_degree_ge3=sum(igraph::degree(graph)>=3L),
      n_valid_branchpoints=sum(rule$audit$candidate_bifurcation),
      selected_branchpoint=rule$selected,
      has_two_major_late_arms=!is.na(rule$selected),
      both_7dpi_samples_have_two_arms=independent_sample_support,
      arm_orientation=arm_orientation,
      min_sample_arm_recall=if(all(is.finite(rr$recall))) min(rr$recall) else NA_real_,
      passes_descriptive_95pct_recall=all(is.finite(rr$recall)) &&
                                           all(rr$recall>=RECALL_THRESHOLD),
      cell_true_day_spearman=suppressWarnings(cor(
          c("1 dpi"=1,"3 dpi"=3,"7 dpi"=7)[ctime],pt,method="spearman")),
      status="SUCCESS",error="",stringsAsFactors=FALSE)
    write_csv(sd,file.path(dd,"01_scenario_topology_summary.csv"))
    # Compact graphical comparison; plotted in one FIXED source UMAP space.
    pal <- c(Arm_1="#B65437",Arm_2="#327CAC",Mapped_Arm_1="#B65437",
             Mapped_Arm_2="#327CAC",Pre_branch_root_side="#BDBDBD",
             Other_downstream="#7E7E7E",Unassigned="#DDDDDD")
    keep7 <- which(ctime=="7 dpi")
    pdf(file.path(dd,"07_7dpi_arm_correspondence.pdf"),width=11,height=5)
    par(mfrow=c(1,2),mar=c(4,4,3,1))
    plot(orig_emb[keep7,1],orig_emb[keep7,2],pch=16,cex=.26,
      col=pal[old[keep7]],xlab="UMAP 1",ylab="UMAP 2",
      main=paste(id,"| original labels"))
    plot(orig_emb[keep7,1],orig_emb[keep7,2],pch=16,cex=.26,
      col=pal[aligned[keep7]],xlab="UMAP 1",ylab="UMAP 2",
      main=paste(id,"| relearned, posthoc aligned"))
    dev.off()
    writeLines("SUCCESS",success_path)
    log_msg("DONE ",id," | selected=",rule$selected,
      " | min two-sample recall=",format(sd$min_sample_arm_recall,digits=3),
      " | both new 7dpi arms by sample=",independent_sample_support)
    rm(obj,graph,corr);gc(verbose=FALSE)
    invisible(TRUE)
  },error=function(e) {
    err <- paste0("stage=",stage," | ",conditionMessage(e))
    log_msg("FAILED ",id," | ",err)
    write_csv(data.frame(scenario_id=id,fraction=param$fraction,seed=param$seed,
      cluster_k=param$cluster_k,
      euclidean_distance_ratio=param$euclidean_distance_ratio,
      n_cells=NA_integer_,n_root_cells=NA_integer_,root_anchor=NA_character_,
      graph_nodes=NA_integer_,graph_edges=NA_integer_,graph_size_matches_reference=FALSE,
      graph_components=NA_integer_,
      n_degree_ge3=NA_integer_,n_valid_branchpoints=NA_integer_,
      selected_branchpoint=NA_character_,has_two_major_late_arms=FALSE,
      both_7dpi_samples_have_two_arms=FALSE,arm_orientation=NA_character_,
      min_sample_arm_recall=NA_real_,passes_descriptive_95pct_recall=FALSE,
      cell_true_day_spearman=NA_real_,status="FAILED",error=err),
      file.path(dd,"01_scenario_topology_summary.csv"))
    writeLines(err,file.path(dd,"FAILURE.txt"))
    gc(verbose=FALSE)
    invisible(FALSE)
  })
}

for (i in seq_len(nrow(scenarios))) {
  run_case(scenarios[i,,drop=FALSE])
  # Without a successful full-data control, variants cannot establish robustness.
  # Stop early instead of repeating the same structural failure 10 times.
  if (i==1L && !file.exists(file.path(scenario_dir(scenarios$scenario_id[1]),
                                 "SUCCESS.flag"))) {
    log_msg("BASELINE FAILED: remaining scenarios intentionally not run. ",
            "Review scenarios/00_full_original_parameters/FAILURE.txt")
    break
  }
}

# Aggregate scenario outputs using ONLY completed analysis files (not uploaded ZIP).
collect <- function(filename) {
  xs <- lapply(scenarios$scenario_id,function(id) {
    path <- file.path(scenario_dir(id),filename)
    if(!file.exists(path)) return(NULL)
    tryCatch(read.csv(path,stringsAsFactors=FALSE),error=function(e)NULL)
  })
  xs <- Filter(function(z) !is.null(z) && nrow(z)>0,xs)
  if(!length(xs)) return(NULL)
  do.call(rbind,xs)
}
agg <- collect("01_scenario_topology_summary.csv")
arm <- collect("04_7dpi_sample_arm_recovery.csv")
branches <- collect("02_branchpoint_candidates.csv")
temporal <- collect("06_pseudotime_by_timepoint.csv")
if(!is.null(agg)) write_csv(agg,file.path(TBL,"06_all_scenario_topology_summary.csv"))
if(!is.null(arm)) write_csv(arm,file.path(TBL,"07_all_scenario_7dpi_sample_arm_recovery.csv"))
if(!is.null(branches)) write_csv(branches,file.path(TBL,"08_all_scenario_branchpoint_candidates.csv"))
if(!is.null(temporal)) write_csv(temporal,file.path(TBL,"09_all_scenario_temporal_pseudotime.csv"))

checks <- data.frame(check=c("input_14005_cells","original_7dpi_2541_cells",
   "primary_root_225_cells","two_7dpi_biological_samples",
   "full_baseline_executed","all_scenarios_successful",
   "baseline_95pct_recovery_and_two_sample_support",
   "baseline_graph_node_edge_counts_match_original"),
 status=c("PASS","PASS","PASS","PASS",
   if(!is.null(agg) && any(agg$scenario_id==scenarios$scenario_id[1] &
                               agg$status=="SUCCESS")) "PASS" else "FAIL",
   if(!is.null(agg) && nrow(agg)==nrow(scenarios) &&
                          all(agg$status=="SUCCESS")) "PASS" else "FAIL",
   if(!is.null(agg) && any(agg$scenario_id==scenarios$scenario_id[1] &
                    agg$status=="SUCCESS" &
                    agg$passes_descriptive_95pct_recall &
                    agg$both_7dpi_samples_have_two_arms)) "PASS" else "FAIL",
    if(!is.null(agg) && any(agg$scenario_id==scenarios$scenario_id[1] &
                    agg$status=="SUCCESS" &
                    agg$graph_size_matches_reference)) "PASS" else "FAIL"))
write_csv(checks,file.path(TBL,"10_STEP6G_audit_checks.csv"))

# Show balanced per-sample arm recovery without falsely treating cells as independent samples.
if(!is.null(arm) && nrow(arm)) {
  arm <- arm[order(arm$scenario_id,arm$sample_id,arm$reference_arm),,drop=FALSE]
  pdf(file.path(FIG,"02_7dpi_two_sample_arm_recovery.pdf"),width=12,height=6)
  vals <- arm$recall
  plot(seq_along(vals),vals,pch=16,cex=.8,ylim=c(0,1.05),xaxt="n",
       xlab="Scenario / sample / original arm",ylab="Recovery among retained cells",
       main="Post-hoc branch correspondence (not fate validation)")
  abline(h=RECALL_THRESHOLD,lty=2,col="gray40")
  axis(1,at=seq_along(vals),labels=paste(arm$scenario_id,arm$sample_id,
      arm$reference_arm,sep=" / "),las=2,cex.axis=.38)
  dev.off()
}

baseline <- if(!is.null(agg)) agg[agg$scenario_id==scenarios$scenario_id[1],,drop=FALSE] else NULL
n_pass <- if(!is.null(agg)) sum(agg$status=="SUCCESS") else 0L
n_stable <- if(!is.null(agg)) sum(agg$status=="SUCCESS" &
   agg$both_7dpi_samples_have_two_arms & agg$passes_descriptive_95pct_recall,
   na.rm=TRUE) else 0L
report <- c(
 "GLIA | Reviewer 2 Major Comment 6 | Step6G trajectory topology sensitivity",
 paste("Run:",as.character(Sys.time())),
 paste("Input (read-only):",IN_CDS),
 paste("Output:",OUT),
 paste("Reference graph nodes/edges:",igraph::vcount(ref_graph),igraph::ecount(ref_graph)),
 paste("Reference root cells:",sum(root_flag),"; original 7dpi Arm1+Arm2:",orig_7),
 paste("Reference 1dpi pseudotime >=",HIGH_1DPI_PT,":",nrow(high_1dpi),
       "/",sum(timepoint=="1 dpi"),
       sprintf("(%.1f%%)",100*nrow(high_1dpi)/sum(timepoint=="1 dpi"))),
 paste("Scenario count:",nrow(scenarios),"; successful:",n_pass),
 paste("Scenarios with both-sample two-arm support AND all arm recall>=",
       RECALL_THRESHOLD,":",n_stable),
 "Descriptive recovery threshold is an INTERNAL sensitivity criterion, not a P-value.",
 "Relearned full-control must recover original graph/arms before variant comparisons are interpreted.",
 "The root/branch Y_ vertex identifiers are local to each graph, not stable identifiers.",
 "Pseudotime units vary across reconstructed graphs; raw values and the >=15 threshold are NOT cross-scenario comparable. Within-case top-quartile fraction is provided.",
 "Original arm labels were used ONLY after objective graph/branch selection, for correspondence.",
 "The UMAP was kept FROZEN throughout: this is NOT an end-to-end PCA/Harmony/UMAP robustness test.",
 "Cross-sectional scRNA-seq cannot establish cell-fate bifurcation or lineage transitions.",
 "Changing root on an unchanged graph tests pseudotime orientation, NOT topology stability.",
 "If baseline failed, subsequent scenarios were intentionally skipped; investigate the baseline FAILURE.txt.",
 "See all tables in tables/; per-scenario cell-level correspondence in scenarios/.",
 if(file.exists(PRIOR_MARSH)) paste("Existing source-matched Marsh sensitivity report:",
                                   PRIOR_MARSH) else
   "Existing Marsh source-matched report not found; no effect on this analysis.",
 "", "ALL-SCENARIO STATUS:",
 if(!is.null(agg)) capture.output(print(agg[,c("scenario_id","status","n_cells",
    "graph_nodes","selected_branchpoint","n_valid_branchpoints",
    "both_7dpi_samples_have_two_arms","min_sample_arm_recall",
    "passes_descriptive_95pct_recall")],row.names=FALSE))
 else "No readable scenario outputs.",
 "", "AUDIT:",capture.output(print(checks,row.names=FALSE)))
writeLines(report,file.path(OUT,"10_STEP6G_FINAL_REPORT.txt"))
log_msg("FINISHED. Review 10_STEP6G_FINAL_REPORT.txt and tables/06-10.")
if(any(checks$status=="FAIL"))
  warning("One or more Step6G checks FAILED. Do not claim robust bifurcation without reviewing the report.")
