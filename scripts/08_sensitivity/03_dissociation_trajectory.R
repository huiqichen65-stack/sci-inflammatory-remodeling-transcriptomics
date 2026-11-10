options(stringsAsFactors=FALSE, width=180)
ROOT <- Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
RES <- file.path(ROOT,"04_results")
BASE <- file.path(RES,"step6FRH_GSE162610_Harmony_matched_trajectory_sensitivity")
OUT <- file.path(BASE,"31N_source_matched_Marsh_paired_trajectory")
TDIR <- file.path(OUT,"tables")
FDIR <- file.path(OUT,"figures")
RDIR <- file.path(OUT,"RDS")
LIB124 <- file.path(BASE,"RDS","isolated_Harmony_1.2.4_library")
IN2C <- file.path(RES,"step2C_GSE162610_final_microglia_extraction","GSE162610_microglia_clean.rds")
IN3A <- file.path(RES,"step3A_GSE162610_microglia_reclustering","GSE162610_step3A_microglia_reclustered_multires.rds")
IN5BR <- file.path(RES,"step5BR_GSE162610_injury_response_trajectory_refinement","GSE162610_step5BR_injury_monocle3_cds.rds")
IN5BRSEU <- file.path(RES,"step5BR_GSE162610_injury_response_trajectory_refinement","GSE162610_step5BR_injury_trajectory_seurat.rds")
IN5BR2 <- file.path(RES,"step5BR2_GSE162610_single_root_branch_audit","GSE162610_step5BR2_single_root_cds.rds")
IN5BR2SEU <- file.path(RES,"step5BR2_GSE162610_single_root_branch_audit","GSE162610_step5BR2_single_root_seurat.rds")
MARSH_XLSX <- file.path(ROOT,"01_raw_data","reference_gene_sets","Marsh2022_Supplementary_Table_4.xlsx")
REPORT <- file.path(OUT,"31N_STEP6FRH_SOURCE_MATCHED_MARSH_PAIRED_FINAL_REPORT.txt")
LOG <- file.path(OUT,"31N_worker_log.txt")
for (d in c(OUT,TDIR,FDIR,RDIR)) dir.create(d,recursive=TRUE,showWarnings=FALSE)
fail <- function(...) stop(paste0(...,collapse=""),call.=FALSE)
write_table <- function(x,n) utils::write.csv(x,file.path(TDIR,n),row.names=FALSE)
needed <- c(IN2C,IN3A,IN5BR,IN5BRSEU,IN5BR2,IN5BR2SEU,MARSH_XLSX)
get_version <- function(p) if(requireNamespace(p,quietly=TRUE)) as.character(utils::packageVersion(p)) else "not_installed"
private_version <- function() {
  q <- file.path(LIB124,"harmony","DESCRIPTION")
  if(!file.exists(q)) return(NA_character_)
  d <- tryCatch(read.dcf(q,fields="Version"),error=function(e) NA)
  if(length(d)==1 && is.na(d)) return(NA_character_)
  as.character(d[1,1])
}
self_path <- function() {
  a <- grep("^--file=",commandArgs(FALSE),value=TRUE)
  if(length(a)) {
    p <- sub("^--file=","",a[1]);if(file.exists(p)) return(normalizePath(p,winslash="/"))
  }
  for(i in rev(seq_along(sys.frames()))) {
    p <- sys.frames()[[i]]$ofile
    if(is.character(p)&&length(p)==1L&&!is.na(p)&&file.exists(p))
      return(normalizePath(p,winslash="/"))
  }
  fail("Cannot locate running script. Use source() with an absolute .R path.")
}
get_counts <- function(x) {
  z <- tryCatch(SeuratObject::GetAssayData(x,assay="RNA",layer="counts"),error=function(e) NULL)
  if(is.null(z)) fail("RNA counts layer not accessible; source input may contain unjoined layers")
  z
}
get_vertex <- function(cds,g) {
  v <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  if(is.null(v)) fail("Monocle3 closest-vertex data absent")
  if(is.matrix(v)||is.data.frame(v)) { nm<-rownames(v);v<-v[,1L];names(v)<-nm }
  if(is.null(names(v))) names(v)<-colnames(cds)
  nm <- names(v)
  v <- as.character(v);gn <- igraph::V(g)$name
  if(all(v%in%gn)) return(stats::setNames(v,nm))
  vv <- paste0("Y_",v)
  if(all(vv%in%gn)) return(stats::setNames(vv,nm))
  ind <- suppressWarnings(as.integer(v))
  if(all(!is.na(ind)&ind>=1L&ind<=length(gn))) return(stats::setNames(gn[ind],nm))
  fail("Unable to translate Monocle3 vertex mapping")
}
get_umap <- function(cds) {
  a <- as.matrix(SingleCellExperiment::reducedDim(cds,"UMAP"))
  if(is.null(rownames(a))) rownames(a)<-colnames(cds)
  if(ncol(a)<2 || !setequal(rownames(a),colnames(cds))) fail("CDS UMAP invalid")
  a[colnames(cds),1:2,drop=FALSE]
}
procrustes <- function(a,b) {
  a <- scale(a,center=TRUE,scale=FALSE);b<-scale(b,center=TRUE,scale=FALSE)
  a <- a/sqrt(sum(a*a));b<-b/sqrt(sum(b*b))
  min(1,max(0,sum(svd(crossprod(b,a),nu=0,nv=0)$d)))
}
if("--31N-worker" %in% commandArgs(trailingOnly=TRUE)) {
  if(!identical(private_version(),"1.2.4")) fail("Private Harmony 1.2.4 library missing: ",LIB124)
  .libPaths(unique(c(normalizePath(LIB124,winslash="/",mustWork=TRUE),.libPaths())))
  if("harmony" %in% loadedNamespaces()) fail("Harmony unexpectedly loaded before isolation")
  suppressPackageStartupMessages(library(harmony,lib.loc=LIB124))
  if(as.character(utils::packageVersion("harmony"))!="1.2.4") fail("Wrong isolated Harmony version")
  pkgs <- c("Seurat","SeuratObject","Matrix","readxl","ggplot2","monocle3", "SingleCellExperiment","SummarizedExperiment","igraph")
  for(p in pkgs) if(!requireNamespace(p,quietly=TRUE)) fail("Required R package missing: ",p)
  for(p in needed) if(!file.exists(p)) fail("Input file missing: ",p)
  input_md5 <- unname(tools::md5sum(needed))
  write_table(data.frame(path=needed,md5=input_md5,bytes=file.info(needed)$size),"01_input_MD5_read_only.csv")
  if(unname(tools::md5sum(IN3A))!="8911e6daad5a4aa94b0b99cabd0f9d7d")
    fail("Original Step3A RDS MD5 changed from previously audited frozen reference")
  write_table(data.frame(package=c("R","harmony",pkgs,"uwot","RANN","leidenbase"), version=vapply(c("R","harmony",pkgs,"uwot","RANN","leidenbase"),function(p) if(p=="R") as.character(getRversion()) else get_version(p),character(1))),"02_package_versions.csv")
  writeLines(utils::capture.output(sessionInfo()),file.path(OUT,"31N_sessionInfo.txt"))
  raw <- as.data.frame(readxl::read_excel(MARSH_XLSX,sheet=1,col_names=FALSE,.name_repair="minimal"))
  headers <- trimws(as.character(unlist(raw[1,,drop=TRUE],use.names=FALSE)))
  extract_genes <- function(label) {
    col <- which(headers==label)
    if(length(col)!=1L) fail("Marsh workbook missing/ambiguous column: ",label)
    z <- trimws(as.character(raw[-1,col,drop=TRUE]))
    unique(z[!is.na(z)&nzchar(z)&toupper(z)!="NA"])
  }
  identity <- extract_genes("Microglial Identity Score")
  cns <- extract_genes("All CNS Cells Act. Score")
  myeloid <- extract_genes("Micro/Myeloid Shared Act. Score")
  if(!identical(as.integer(c(length(identity),length(cns),length(myeloid))),c(18L,19L,25L)))
    fail("Marsh Table 4 unexpected gene numbers")
  marsh <- sort(unique(c(cns,myeloid)))
  if(length(marsh)!=34L) fail("Marsh activation union must contain 34 unique genes")
  ref3 <- readRDS(IN3A)
  src_cds <- readRDS(IN5BR)
  ref_cds <- readRDS(IN5BR2)
  src_obj <- readRDS(IN5BRSEU)
  ref_obj <- readRDS(IN5BR2SEU)
  if(!inherits(ref3,"Seurat")) fail("Original Step3A reference not Seurat")
  full_cells <- colnames(ref3); full_features <- rownames(ref3)
  if(length(full_cells)!=20747L||length(full_features)!=33177L||anyDuplicated(full_cells))
    fail("Step3A benchmark 20,747 cells / 33,177 features invalid")
  rc <- colnames(src_cds)
  if(length(rc)!=14005L||anyDuplicated(rc)||!setequal(rc,colnames(ref_cds)))
    fail("Source CDS and frozen Step5BR2 cell cohort disagreement")
  if(!all(rc%in%full_cells)||!all(rc%in%colnames(src_obj))||!all(rc%in%colnames(ref_obj)))
    fail("Frozen trajectory barcode missing from parent Step3A or Step5BR/BR2 Seurat")
  md <- src_obj[[]][match(rc,rownames(src_obj[[]])),,drop=FALSE]
  labels <- ref_obj[[]][match(rc,rownames(ref_obj[[]])),,drop=FALSE]
  if(!identical(rownames(md),rc)||!identical(rownames(labels),rc)) fail("Metadata barcode match failure")
  req <- c("step5BR_primary_root","sample_id","timepoint","shared_transition_3dpi_highconf")
  if(!all(req%in%colnames(md))||!"step5BR2_branch_arm"%in%colnames(labels))
    fail("Missing frozen Step5BR/BR2 metadata fields")
  roots <- rc[as.logical(md$step5BR_primary_root)%in%TRUE]
  if(length(roots)!=225L) fail("Expected 225 original frozen root cells")
  tp <- as.character(md$timepoint);samp <- as.character(md$sample_id)
  trans <- as.logical(md$shared_transition_3dpi_highconf)
  trans[is.na(trans)] <- FALSE
  original_arm <- as.character(labels$step5BR2_branch_arm)
  if(sum(tp=="7 dpi"&original_arm=="Arm_1")!=1399L|| sum(tp=="7 dpi"&original_arm=="Arm_2")!=1142L) fail("Frozen 7dpi arm identities changed")
  samples7 <- sort(unique(samp[tp=="7 dpi"]))
  if(length(samples7)!=2L) fail("Expected two 7dpi biological samples")
  src_count <- SummarizedExperiment::assay(src_cds,"counts")
  if(!identical(colnames(src_count),rc)) src_count<-src_count[,rc,drop=FALSE]
  present <- intersect(marsh,rownames(src_count))
  if(length(present)!=26L) fail("Detected Marsh activation genes !=26; verify original workbook")
  write_table(data.frame(gene=marsh,in_source_counts=marsh%in%rownames(src_count), All_CNS_Activation=marsh%in%cns,Micro_Myeloid_Activation=marsh%in%myeloid), "03_Marsh_gene_exclusion_manifest.csv")
  write_table(data.frame(barcode=rc,sample_id=samp,timepoint=tp, original_arm=original_arm,primary_root=rc%in%roots),"04_frozen_cell_manifest.csv")
  ref_pca <- SeuratObject::Embeddings(ref3[["pca.step3A"]])[full_cells,1:50,drop=FALSE]
  ref_h <- SeuratObject::Embeddings(ref3[["harmony.step3A"]])[full_cells,1:30,drop=FALSE]
  ref_umap <- SeuratObject::Embeddings(ref3[["umap.microglia.step3A"]])[full_cells,1:2,drop=FALSE]
  ref_hvg <- SeuratObject::VariableFeatures(ref3)
  ref_graph <- monocle3::principal_graph(ref_cds)[["UMAP"]]
  if(is.null(ref_graph)) fail("Frozen reference Monocle graph absent")
  compute_step3A <- function(condition,excluded) {
    x <- readRDS(IN2C)
    if(!inherits(x,"Seurat")||!identical(colnames(x),full_cells)|| !identical(rownames(x),full_features))
      fail(condition,": Step2C cell/feature ORDER differs from historical Step3A")
    SeuratObject::DefaultAssay(x)<-"RNA"
    if(exists("JoinLayers",envir=asNamespace("SeuratObject"),mode="function",inherits=FALSE))
      x <- tryCatch(SeuratObject::JoinLayers(x,assay="RNA"),error=function(e){
        x})
    ct <- get_counts(x)
    if(!identical(colnames(ct),full_cells)||!identical(rownames(ct),full_features))
      fail(condition,": RNA count order mismatch")
    mask_present <- intersect(excluded,rownames(ct))
    if(length(excluded) && !setequal(mask_present,present)) fail("Marsh masked gene set mismatch")
    if(length(mask_present)) {
      x <- subset(x,features=setdiff(rownames(x),mask_present))
      SeuratObject::DefaultAssay(x)<-"RNA"
      if(any(mask_present%in%rownames(x))||ncol(x)!=20747L || !identical(colnames(x),full_cells) || !identical(rownames(x),setdiff(full_features,mask_present)))
        fail("Mask did not preserve exact feature/cell order after removing 26 genes")
    }
    rm(ct);invisible(gc())
    set.seed(1234L)
    x <- Seurat::NormalizeData(x,assay="RNA",normalization.method="LogNormalize", scale.factor=10000,verbose=FALSE)
    x <- Seurat::FindVariableFeatures(x,assay="RNA",selection.method="vst", nfeatures=3000,verbose=FALSE)
    hv <- SeuratObject::VariableFeatures(x)
    if(length(hv)!=3000L||any(hv%in%mask_present)) fail("HVG computation invalid")
    if(!length(excluded) && !identical(hv,ref_hvg))
      fail("FULL CONTROL did not recover original 3000 HVGs. STOP before masking.")
    p <- Seurat::VariableFeaturePlot(x)
    p <- Seurat::LabelPoints(plot=p,points=head(hv,10),repel=TRUE)
    tmp <- tempfile(fileext=".pdf")
    ggplot2::ggsave(tmp,p,width=8,height=6);unlink(tmp)
    x <- Seurat::ScaleData(x,assay="RNA",features=hv,verbose=FALSE)
    x <- Seurat::RunPCA(x,assay="RNA",features=hv,npcs=50, reduction.name="pca.step3A",reduction.key="PC3A_",verbose=FALSE)
    pc <- SeuratObject::Embeddings(x[["pca.step3A"]])[full_cells,1:50,drop=FALSE]
    if(!length(excluded)) {
      abs_cor <- vapply(1:50,function(j)abs(stats::cor(pc[,j],ref_pca[,j])),numeric(1))
      pc_sim <- procrustes(ref_pca,pc)
      write_table(data.frame(stage="PCA50",min_abs_component_pearson=min(abs_cor), procrustes_similarity=pc_sim),"05_full_control_PCA_reproduction.csv")
      if(min(abs_cor)<0.999999 || pc_sim<0.999999)
        fail("Full-gene PCA is not faithfully reproduced. STOP before masking.")
    }
    set.seed(1234L)
    x <- Seurat::RunUMAP(x,reduction="pca.step3A",dims=1:30, reduction.name="umap.preHarmony.step3A",reduction.key="preH3A_", seed.use=1234L,verbose=FALSE)
    x <- harmony::RunHarmony(object=x,group.by.vars="sample_id", reduction.use="pca.step3A",dims.use=1:30,theta=1, max.iter.harmony=20L,reduction.save="harmony.step3A",verbose=FALSE)
    hh <- SeuratObject::Embeddings(x[["harmony.step3A"]])[full_cells,1:30,drop=FALSE]
    if(!length(excluded)) {
      e <- max(abs(hh-ref_h))
      write_table(data.frame(stage="Harmony30",max_abs_coordinate_delta=e, procrustes_similarity=procrustes(ref_h,hh)),"06_full_control_Harmony_reproduction.csv")
      if(!is.finite(e)||e>1e-3) fail("Full control Harmony did not match source reference; STOP before masking; max diff=",signif(e,6))
    }
    x <- Seurat::RunUMAP(x,reduction="harmony.step3A",dims=1:30, reduction.name="umap.microglia.step3A",reduction.key="mg3A_", seed.use=1234L,verbose=FALSE)
    um <- SeuratObject::Embeddings(x[["umap.microglia.step3A"]])[full_cells,1:2,drop=FALSE]
    if(any(!is.finite(um))) fail("Nonfinite UMAP coords")
    if(!length(excluded)) {
      e <- max(abs(um-ref_umap))
      write_table(data.frame(stage="UMAP2",max_abs_coordinate_delta=e, procrustes_similarity=procrustes(ref_umap,um)),"07_full_control_UMAP_reproduction.csv")
      if(!is.finite(e)||e>1e-3) fail("Full control UMAP did not match source reference; STOP before masking; max diff=",signif(e,6))
    }
    write_table(data.frame(condition=condition,n_cells=ncol(x),n_features=nrow(x), detected_removed_genes=length(mask_present),n_HVG=length(hv), overlap_HVG_with_Marsh=sum(hv%in%marsh)),paste0("08_",condition,"_Step3A_summary.csv"))
    saveRDS(list(condition=condition,umap=um,hvg=hv,removed=mask_present, harmony_version=as.character(utils::packageVersion("harmony"))), file.path(RDIR,paste0("31N_",condition,"_Step3A_UMAP.rds")))
    rm(x);invisible(gc());um
  }
  infer_graph <- function(condition,um,excluded) {
    count_mat <- src_count
    if(length(excluded)) count_mat<-count_mat[!rownames(count_mat)%in%excluded,,drop=FALSE]
    if(!identical(colnames(count_mat),rc)||any(rownames(count_mat)%in%excluded))
      fail("Masked Monocle counts not aligned or excluded genes remain")
    cell_md <- as.data.frame(SummarizedExperiment::colData(src_cds))
    cell_md <- cell_md[rc,,drop=FALSE]
    gene_md <- data.frame(gene_short_name=rownames(count_mat),row.names=rownames(count_mat))
    cds <- monocle3::new_cell_data_set(count_mat,cell_metadata=cell_md,gene_metadata=gene_md)
    set.seed(20260731L)
    cds <- monocle3::preprocess_cds(cds,num_dim=50L,method="PCA")
    u <- um[rc,1:2,drop=FALSE];colnames(u)<-c("UMAP_1","UMAP_2")
    SingleCellExperiment::reducedDim(cds,"UMAP")<-u
    if(max(abs(get_umap(cds)-u))>1e-8) fail("Monocle imported UMAP changed")
    set.seed(20260731L)
    cds <- monocle3::cluster_cells(cds,reduction_method="UMAP",k=30L, cluster_method="leiden",random_seed=20260731L,verbose=FALSE)
    cds <- monocle3::learn_graph(cds,use_partition=FALSE,close_loop=FALSE, learn_graph_control=list(euclidean_distance_ratio=0.8))
    if(max(abs(get_umap(cds)-u))>1e-8) fail("UMAP changed while learning graph")
    g <- monocle3::principal_graph(cds)[["UMAP"]]
    if(is.null(g)||igraph::components(g)$no!=1L) fail("No connected Monocle graph")
    cv <- get_vertex(cds,g)[rc]
    if(anyNA(cv)) fail("Missing cell to vertex projection")
    rv <- sort(table(cv[roots]),decreasing=TRUE)
    if(!length(rv)) fail("Cannot define original 225-cell root anchor")
    root_anchor <- names(rv)[1]
    cds <- monocle3::order_cells(cds,reduction_method="UMAP",root_pr_nodes=root_anchor)
    pt <- as.numeric(monocle3::pseudotime(cds,reduction_method="UMAP")[rc])
    graph_dist <- as.numeric(igraph::distances(g,v=root_anchor,to=igraph::V(g))[1,])
    names(graph_dist)<-igraph::V(g)$name
    deg <- igraph::degree(g)
    br <- names(deg[deg>=3L]);rows<-list();cache<-list()
    for(bp in br) {
      if(identical(bp,root_anchor)) next
      cut <- igraph::delete_vertices(g,bp)
      memb <- igraph::components(cut)$membership
      names(memb)<-igraph::V(cut)$name
      root_side <- unname(memb[root_anchor])
      if(length(root_side)!=1L||is.na(root_side)) next
      comp <- as.integer(memb[cv]);names(comp)<-rc
      valid <- !is.na(comp)
      groups <- sort(unique(comp[valid]));groups<-groups[groups!=root_side]
      if(!length(groups)) next
      d <- data.frame(component=groups,n_cells=integer(length(groups)), n_7dpi=integer(length(groups)))
      for(i in seq_along(groups)) {
        flag<-valid&comp==groups[i]
        d$n_cells[i]<-sum(flag);d$n_7dpi[i]<-sum(flag&tp=="7 dpi")
      }
      d<-d[order(-d$n_7dpi,-d$n_cells),,drop=FALSE]
      major<-d[d$n_cells>=100L&d$n_7dpi>=100L,,drop=FALSE]
      distance_bp<-as.numeric(igraph::distances(g,v=bp,to=igraph::V(g))[1,])
      near_vertices<-igraph::V(g)$name[distance_bp<=2]
      near_transition<-sum(cv%in%near_vertices&tp=="3 dpi"&trans)
      rows[[length(rows)+1L]]<-data.frame(condition=condition,branch_vertex=bp, degree=as.integer(deg[bp]),downstream_components=nrow(d), qualifying_downstream_arms=nrow(major), min_top2_7dpi=if(nrow(major)>=2L)min(major$n_7dpi[1:2])else 0L, near_transition_cells=near_transition, distance_from_root=graph_dist[[bp]],candidate=nrow(major)>=2L)
      cache[[bp]]<-list(member=comp,major=major,root_side=root_side)
    }
    audit <- if(length(rows)) do.call(rbind,rows) else data.frame( condition=character(),branch_vertex=character(),degree=integer(), downstream_components=integer(),qualifying_downstream_arms=integer(), min_top2_7dpi=integer(),near_transition_cells=integer(), distance_from_root=numeric(),candidate=logical())
    if(nrow(audit)) audit<-audit[order(-as.integer(audit$candidate), -audit$near_transition_cells,-audit$min_top2_7dpi,audit$distance_from_root),,drop=FALSE]
    write_table(audit,paste0("09_",condition,"_objective_branchpoint_audit.csv"))
    bp <- if(nrow(audit)&&isTRUE(audit$candidate[1])) as.character(audit$branch_vertex[1]) else NA_character_
    fresh <- stats::setNames(rep("Unassigned",length(rc)),rc)
    if(!is.na(bp)) {
      z <- cache[[bp]]
      fresh[!is.na(z$member)&z$member==z$root_side] <- "Root_side"
      top <- head(z$major$component,2)
      if(length(top)==2L) for(j in 1:2)
        fresh[!is.na(z$member)&z$member==top[j]] <- paste0("New_Arm_",j)
      other <- !is.na(z$member)&z$member!=z$root_side&fresh=="Unassigned"
      fresh[other]<-"Other_downstream"
    }
    ix <- tp=="7 dpi"&original_arm%in%c("Arm_1","Arm_2")
    choices <- list(c("New_Arm_1","New_Arm_2"),c("New_Arm_2","New_Arm_1"))
    hits <- vapply(choices,function(z)sum((original_arm[ix]=="Arm_1"&fresh[rc[ix]]==z[1])| (original_arm[ix]=="Arm_2"&fresh[rc[ix]]==z[2])),integer(1))
    orient <- choices[[which.max(hits)]]
    aligned <- stats::setNames(rep("Other",length(rc)),rc)
    aligned[fresh==orient[1]]<-"Aligned_Alternative"
    aligned[fresh==orient[2]]<-"Aligned_Target"
    aligned[fresh=="Root_side"]<-"Root_side"
    eval <- list()
    for(ss in c(samples7,"ALL_7dpi")) for(j in 1:2) {
      flag <- ix&original_arm==paste0("Arm_",j)&(ss=="ALL_7dpi"|samp==ss)
      n<-sum(flag);found<-sum(fresh[rc[flag]]==orient[j])
      eval[[length(eval)+1L]]<-data.frame(condition=condition,sample_id=ss, original_arm=paste0("Arm_",j),mapped_new_arm=orient[j], n_original=n,n_recovered=found,recall=if(n>0)found/n else NA_real_)
    }
    eval <- do.call(rbind,eval)
    write_table(eval,paste0("10_",condition,"_7dpi_arm_recovery.csv"))
    cell_tbl <- data.frame(barcode=rc,sample_id=samp,timepoint=tp, original_arm=original_arm,new_arm=unname(fresh[rc]), aligned_biological_arm=unname(aligned[rc]),pseudotime=pt, UMAP_1=u[,1],UMAP_2=u[,2])
    write_table(cell_tbl,paste0("11_",condition,"_cell_level_correspondence.csv"))
    if(!all(is.finite(pt))) fail("Some Monocle pseudotime values are not finite in ",condition)
    stats <- data.frame(condition=condition,n_cells=length(rc),root_cells=length(roots), nodes=igraph::vcount(g),edges=igraph::ecount(g), components=igraph::components(g)$no,degree_ge3=sum(deg>=3L), root_anchor=root_anchor,objective_branchpoint=bp, finite_pseudotime_fraction=mean(is.finite(pt)))
    write_table(stats,paste0("12_",condition,"_graph_summary.csv"))
    saveRDS(cds,file.path(RDIR,paste0("31N_",condition,"_relearned_Monocle_CDS.rds")))
    list(stats=stats,evaluation=eval,cells=cell_tbl,branchpoint=bp,passed=(!is.na(bp)&& all(eval$recall[eval$sample_id!="ALL_7dpi"]>=0.95,na.rm=FALSE)), umap=u)
  }
  um_control <- compute_step3A("all_genes_control",character())
  control <- infer_graph("all_genes_control",um_control,character())
  gate <- data.frame(condition="all_genes_control",candidate_present=!is.na(control$branchpoint), recovered_all_four_sample_arms=isTRUE(control$passed), minimum_sample_arm_recall=min(control$evaluation$recall[control$evaluation$sample_id!="ALL_7dpi"]), passes_control_gate=isTRUE(control$passed))
  write_table(gate,"13_control_validity_gate.csv")
  if(!isTRUE(control$passed)) {
    writeLines(c("31N PAIRED ANALYSIS: CONTROL VALIDITY GATE FAILED", "No Marsh exclusion was executed or interpreted.", "The full-gene control did not recover both original biological arms in BOTH 7-dpi samples.", utils::capture.output(print(gate,row.names=FALSE)), "Inspect 05/06/07 full-gene stage metrics and 09/10/12 control graph tables."), REPORT,useBytes=TRUE)
    fail("Full-gene control did not reproduce frozen arms. STOP; see ",REPORT)
  }
  um_mask <- compute_step3A("Marsh_activation_masked",present)
  masked <- infer_graph("Marsh_activation_masked",um_mask,present)
  joint <- rbind(control$evaluation,masked$evaluation)
  write_table(joint,"14_paired_sample_level_recovery_comparison.csv")
  js <- merge(control$cells[,c("barcode","new_arm","aligned_biological_arm")], masked$cells[,c("barcode","new_arm","aligned_biological_arm")],by="barcode", suffixes=c("_control","_masked"),sort=FALSE)
  js<-js[match(rc,js$barcode),,drop=FALSE]
  js$sample_id<-samp;js$timepoint<-tp;js$original_arm<-original_arm
  write_table(js,"15_paired_cell_level_comparison.csv")
  paired_umap_similarity <- procrustes(um_control,um_mask)
  paired_agree <- mean(js$aligned_biological_arm_control==js$aligned_biological_arm_masked)
  ix_labeled <- js$timepoint=="7 dpi" & js$original_arm%in%c("Arm_1","Arm_2")
  paired_agree_labeled <- mean(js$aligned_biological_arm_control[ix_labeled]== js$aligned_biological_arm_masked[ix_labeled])
  summary <- data.frame(condition=c("all_genes_control","Marsh_activation_masked"), branchpoint=c(control$branchpoint,masked$branchpoint), branch_supported_both_samples=c(control$passed,masked$passed), graph_nodes=c(control$stats$nodes,masked$stats$nodes), graph_edges=c(control$stats$edges,masked$stats$edges), min_7dpi_sample_arm_recall=c( min(control$evaluation$recall[control$evaluation$sample_id!="ALL_7dpi"]), min(masked$evaluation$recall[masked$evaluation$sample_id!="ALL_7dpi"])))
  write_table(summary,"16_paired_final_topology_summary.csv")
  plot_three <- function(x,colour,title) {
    ii<-tp=="7 dpi";code<-c(Aligned_Alternative="#3A78B0",Aligned_Target="#C76049", Root_side="#BBBBBB",Other="#D2D2D2")
    cc<-unname(code[as.character(colour[ii])]);cc[is.na(cc)]<-"#D2D2D2"
    plot(x[ii,1],x[ii,2],col=cc,pch=16,cex=0.24, xlab="UMAP 1",ylab="UMAP 2",main=title)
  }
  grDevices::pdf(file.path(FDIR,"31N_original_fullgene_masked_7dpi_arms.pdf"),width=16,height=5.4)
  op<-par(mfrow=c(1,3),mar=c(4,4,3,1))
  plot_three(ref_umap[rc,,drop=FALSE],ifelse(original_arm=="Arm_1","Aligned_Alternative", ifelse(original_arm=="Arm_2","Aligned_Target","Other")),"Original Step5BR2")
  plot_three(control$umap,control$cells$aligned_biological_arm,"Full gene control")
  plot_three(masked$umap,masked$cells$aligned_biological_arm,"Marsh activation masked")
  par(op);grDevices::dev.off()
  full_ok <- isTRUE(control$passed) && isTRUE(masked$passed)
  text <- c( "GLIA Reviewer 2 Major Comment 2 | Source-matched paired Marsh trajectory sensitivity", paste("Completed:",as.character(Sys.time())), "Input objects READ ONLY. Harmony 1.2.4 installed only in existing private library.", "Historical Harmony version UNKNOWN; candidate 1.2.4 reproduces 31M original Harmony under original source order.", paste("Original Step3A MD5:",unname(tools::md5sum(IN3A))), paste("Cells Step3A:",length(full_cells),"; frozen trajectory:",length(rc),"; roots:",length(roots)), paste("Marsh activation union:",length(marsh),"; detected/excluded:",length(present)), "Marsh identity signature genes are not additionally removed.", "Both arms are chosen from objective branch rules; original arm labels are used ONLY for post hoc barcode evaluation.", paste("Control full-gene validity gate:",if(control$passed)"PASS" else "FAIL"), paste("Masked source-matched two-sample topology:",if(masked$passed)"PASS" else "NOT SUPPORTED"), paste("Paired UMAP Procrustes similarity:",signif(paired_umap_similarity,8)), "Predeclared technical gate: each 7dpi sample must recover >=95% of BOTH frozen Arm_1 and Arm_2 barcode sets.", paste("Paired mapped label agreement (all 14,005 cells; descriptive only):",signif(paired_agree,8)), paste("Paired mapped label agreement (originally labeled 7 dpi cells):",signif(paired_agree_labeled,8)), "This is an internal sensitivity analysis, not independent validation; no causal lineage transition is inferred.", "PASS applies to a predefined operational cell-barcode correspondence test, not proof of biological causality.", "If masking is not supported, report loss of support; claim that the original branch is artefactual.", "-- Per-condition summary --",utils::capture.output(print(summary,row.names=FALSE)), "-- Per-sample and pooled 7 dpi recall --",utils::capture.output(print(joint,row.names=FALSE)), paste("Outputs:",OUT), paste("Figure:",file.path(FDIR,"31N_original_fullgene_masked_7dpi_arms.pdf")) )
  writeLines(text,REPORT,useBytes=TRUE)
  quit(save="no",status=0L,runLast=FALSE)
}
if(!identical(private_version(),"1.2.4"))
  fail("Private Harmony 1.2.4 from 31K is missing: ",LIB124, "\nalter globally installed Harmony 2.0.5.")
for(p in needed) if(!file.exists(p)) fail("Required source input missing: ",p)
self<-self_path()
commands<-unique(c(Sys.which("Rscript"),file.path(R.home("bin"),"Rscript.exe"), file.path(R.home(),"bin","Rscript.exe")))
commands<-commands[nzchar(commands)&file.exists(commands)]
if(!length(commands)) fail("Rscript executable not found")
if(file.exists(REPORT)) unlink(REPORT)
status<-tryCatch(system2(commands[1],args=c("--vanilla",shQuote(self),"--31N-worker"), stdout=LOG,stderr=LOG,wait=TRUE), error=function(e){writeLines(conditionMessage(e),LOG);127L})
if(!identical(as.integer(status),0L)||!file.exists(REPORT))
  fail("31N stopped or failed. Read ",LOG, " and any report in ",OUT, ". In particular, baseline control gate MUST pass before masked interpretation.")
