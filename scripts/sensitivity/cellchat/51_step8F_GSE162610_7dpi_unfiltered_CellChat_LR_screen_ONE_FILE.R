# ==============================================================================
# 51_step8F_GSE162610_7dpi_unfiltered_CellChat_LR_screen_ONE_FILE.R
# GLIA Reviewer 1 | Step 8F | Full-database, full-cell, 7-dpi LR sensitivity
# ==============================================================================
# IMPORTANT: The previous Step8B full-cell nboot=100 run was restricted to 20
# preview-prioritized LR interactions. The older all-database Step8B PREVIEW
# used downsampling and nboot=20. Neither is a full-cell all-database analysis.
# This script estimates the 7-dpi network de novo using all eligible 7-dpi
# cells and CellChatDB.mouse WITHOUT Step7G or the earlier 18/20-axis list.
# Standard expression-based CellChat screening is used, not prior target-program
# evidence. Rankings are exploratory within the detected LR interactions.
# It does NOT test downstream biological responses for all recipient pathways,
# and does NOT constitute independent biological validation.
#
# Run in RStudio:
# source("<PROJECT_ROOT>/51_step8F_GSE162610_7dpi_unfiltered_CellChat_LR_screen_ONE_FILE.R", encoding="UTF-8")
# This is substantially heavier than the earlier 20-interaction confirmation.
# ==============================================================================

options(stringsAsFactors=FALSE, width=240, timeout=86400,
        future.globals.maxSize=40*1024^3)
set.seed(20260803L)
ROOT <- Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
OUT <- file.path(ROOT,"04_results","step8F_GSE162610_7dpi_unfiltered_CellChat_LR_screen")
TAB <- file.path(OUT,"tables")
RDS <- file.path(OUT,"RDS")
REP <- file.path(OUT,"reports")
for(p in c(OUT,TAB,RDS,REP)) dir.create(p,recursive=TRUE,showWarnings=FALSE)
ATLAS <- file.path(ROOT,"04_results","step1F_GSE162610_broad_annotation",
                   "GSE162610_broad_annotation_audit.rds")
MANIFEST <- file.path(ROOT,"04_results","step8A_GSE162610_CellChat_input_preflight",
                      "RDS","GSE162610_step8A_CellChat_input_manifest.rds")
OLD18 <- file.path(ROOT,"04_results","step8B_CONFIRM_GSE162610_7dpi_targeted_fullcell_nboot100",
                   "tables","07_confirmed_core_axes.tsv")
CHECKPOINT <- file.path(RDS,"GSE162610_step8F_7dpi_full_database_nboot100_CellChat.rds")
RUN_REPORT <- file.path(REP,"Step8F_unfiltered_LR_screen_report.txt")
ERROR_FILE <- file.path(REP,"Step8F_ERROR.txt")
NBOOT <- 100L
MIN_GROUP_CELLS <- 20L
MIN_SAMPLE_GROUP_CELLS <- 20L
MIN_SENDER_DETECTION <- 0.01  # uniformly required in EACH sample
MIN_RECEPTOR_DETECTION <- 0.05 # uniformly required in EACH sample
CELLCHAT_P_THRESHOLD <- 0.05
FORMAL_SENDER <- "Microglia_Arm_2"
COMPARATOR <- "Microglia_Arm_1"
EXCLUDE_RECIPIENTS <- c("Microglia_Arm_1","Microglia_Arm_2","Microglia_Shared_prebranch")

fail <- function(...) {
  message <- paste0(...,collapse="")
  writeLines(message,ERROR_FILE,useBytes=TRUE)
  stop(message,call.=FALSE)
}
read_tsv <- function(x) data.table::fread(x,data.table=FALSE,showProgress=FALSE)
write_tsv <- function(x,name) data.table::fwrite(x,file.path(TAB,name),sep="\t",na="NA")
call_supported <- function(fn, args) {
  f <- names(formals(fn))
  if(!"..." %in% f) args <- args[names(args) %in% f]
  do.call(fn,args)
}
as_flag <- function(x) {
  if(is.logical(x)) return(x)
  toupper(trimws(as.character(x))) %in% c("TRUE","T","1","YES","PASS")
}
read_rna <- function(x) {
  out <- tryCatch(SeuratObject::LayerData(x,assay="RNA",layer="data"),error=function(e) NULL)
  if(is.null(out)) out <- tryCatch(Seurat::GetAssayData(x,assay="RNA",slot="data"),error=function(e) NULL)
  if(is.null(out)) fail("No RNA normalized layer in frozen atlas")
  out
}
load_db <- function() {
  db <- tryCatch(get("CellChatDB.mouse",envir=asNamespace("CellChat"),inherits=FALSE),error=function(e) NULL)
  if(!is.null(db)) return(db)
  e <- new.env(parent=emptyenv())
  suppressWarnings(try(utils::data("CellChatDB.mouse",package="CellChat",envir=e),silent=TRUE))
  if(!exists("CellChatDB.mouse",envir=e,inherits=FALSE)) fail("Cannot load CellChatDB.mouse")
  get("CellChatDB.mouse",envir=e)
}
run <- function() {
  required <- c("CellChat","Seurat","SeuratObject","Matrix","data.table","future")
  missing_packages <- required[!vapply(required,requireNamespace,logical(1),quietly=TRUE)]
  if(length(missing_packages)) fail("Missing R packages: ",paste(missing_packages,collapse=", "))
  if(!file.exists(ATLAS) || !file.exists(MANIFEST)) fail("Missing frozen atlas or Step8A manifest; check project paths")
  manifest <- readRDS(MANIFEST)
  if(!identical(manifest$status,"PASS")) fail("Step8A manifest status is not PASS")
  m <- as.data.frame(manifest$cell_metadata,stringsAsFactors=FALSE)
  need <- c("cell","sample_id","timepoint_standard","cellchat_group","include_step8_timecourse")
  if(length(setdiff(need,names(m)))) fail("Missing Step8A metadata columns: ",paste(setdiff(need,names(m)),collapse=", "))
  eligible <- as.character(manifest$eligible_groups_by_time[["7dpi"]])
  m <- m[as_flag(m$include_step8_timecourse) & as.character(m$timepoint_standard)=="7dpi" &
         as.character(m$cellchat_group) %in% eligible,,drop=FALSE]
  m[] <- lapply(m,function(x) if(is.factor(x)) as.character(x) else x)
  m <- m[order(m$cell),,drop=FALSE]
  if(anyDuplicated(m$cell) || !nrow(m)) fail("Invalid 7dpi cell identifiers")
  if(length(unique(m$sample_id))!=2L) fail("Expected two 7dpi biological samples; got ",length(unique(m$sample_id)))
  groups <- sort(unique(as.character(m$cellchat_group)))
  if(!all(c(FORMAL_SENDER,COMPARATOR) %in% groups)) fail("Missing frozen microglial sender groups")
  if(any(table(m$cellchat_group) < MIN_GROUP_CELLS)) fail("Some groups have <20 cells; audit Step8A")
  if(nrow(m)!=14505L) fail("Expected 14,505 eligible 7dpi cells; found ",nrow(m),"; STOP for audit")
  count_table <- as.data.frame(table(sample=m$sample_id,group=m$cellchat_group),stringsAsFactors=FALSE)
  names(count_table)[3] <- "n_cells"
  write_tsv(count_table,"00_7dpi_cell_counts_by_sample_group.tsv")
  message("Step8F: input audit passed; loading full 7dpi RNA data")
  atlas <- readRDS(ATLAS)
  if(!inherits(atlas,"Seurat")) fail("Frozen atlas is not a Seurat object")
  Seurat::DefaultAssay(atlas) <- "RNA"
  expr <- read_rna(atlas)
  if(length(setdiff(m$cell,colnames(expr)))) fail("Some 7dpi cells missing from RNA layer")
  expr <- expr[,m$cell,drop=FALSE]
  rm(atlas); invisible(gc())
  db <- load_db()
  interaction_db <- as.data.frame(db$interaction,stringsAsFactors=FALSE)
  if(nrow(interaction_db) < 2000L) fail("Unexpectedly restricted CellChatDB.mouse; expected full database")
  if(!"interaction_name" %in% names(interaction_db)) interaction_db$interaction_name <- rownames(interaction_db)
  expected_interactions <- sort(unique(as.character(interaction_db$interaction_name)))
  if(anyNA(expected_interactions)) fail("CellChat database interaction identifiers missing")
  parms <- list(nboot=NBOOT,type="triMean",trim=0.10,raw.use=TRUE,
                population.size=FALSE,distance.use=FALSE,seed.use=20260803L)
  cc <- NULL
  if(file.exists(CHECKPOINT)) {
    cp <- tryCatch(readRDS(CHECKPOINT),error=function(e)NULL)
    if(!is.null(cp) && inherits(cp$cellchat,"CellChat") &&
       identical(cp$cells,as.character(m$cell)) &&
       identical(cp$db_interactions,expected_interactions) &&
       identical(cp$parameters,parms)) {
      cc <- cp$cellchat
      message("Reusing validated full-database Step8F checkpoint")
    } else fail("Step8F checkpoint exists but parameters/cells/DB differ; manually inspect it before retrying")
  }
  if(is.null(cc)) {
    future::plan("sequential")
    meta <- data.frame(cellchat_group=factor(m$cellchat_group,levels=groups),
                       samples=as.character(m$sample_id),row.names=m$cell)
    cc <- CellChat::createCellChat(object=expr,meta=meta,group.by="cellchat_group")
    cc@DB <- db  # FULL database; no Step7G, no 18/20 LR shortlist
    cc <- CellChat::subsetData(cc)
    cc <- CellChat::identifyOverExpressedGenes(cc)
    cc <- CellChat::identifyOverExpressedInteractions(cc)
    n_lr <- nrow(cc@LR$LRsig)
    if(!length(n_lr) || is.na(n_lr) || n_lr<1L) fail("Expression-only LR screening selected zero interactions")
    message("Unrestricted expression-based CellChat LR count: ",n_lr,
            "; computing on 14,505 cells, nboot=100")
    cc <- call_supported(CellChat::computeCommunProb,
                         c(list(object=cc),parms))
    cc <- call_supported(CellChat::filterCommunication,list(object=cc,min.cells=MIN_GROUP_CELLS))
    cc <- call_supported(CellChat::computeCommunProbPathway,list(object=cc,thresh=CELLCHAT_P_THRESHOLD))
    cc <- call_supported(CellChat::aggregateNet,list(object=cc,thresh=CELLCHAT_P_THRESHOLD,return.object=TRUE))
    saveRDS(list(cellchat=cc,cells=as.character(m$cell),
                 db_interactions=expected_interactions,parameters=parms),
            CHECKPOINT,compress=FALSE)
    message("Full-database CellChat checkpoint saved")
  }
  raw <- call_supported(CellChat::subsetCommunication,
                        list(object=cc,slot.name="net",thresh=1))
  raw <- as.data.frame(raw,stringsAsFactors=FALSE)
  if(!nrow(raw)) fail("No CellChat communication rows returned")
  get_col <- function(z,nms,default=NA_character_) {
    hit <- intersect(nms,names(z))
    if(length(hit)) z[[hit[1L]]] else rep(default,nrow(z))
  }
  all_lr <- data.frame(source=as.character(get_col(raw,"source")),
                       target=as.character(get_col(raw,"target")),
                       ligand=as.character(get_col(raw,"ligand")),
                       receptor=as.character(get_col(raw,"receptor")),
                       interaction_name=as.character(get_col(raw,c("interaction_name","interaction_name_2"))),
                       pathway_name=as.character(get_col(raw,c("pathway_name","pathway"))),
                       probability=suppressWarnings(as.numeric(get_col(raw,c("prob","probability")))),
                       p_value=suppressWarnings(as.numeric(get_col(raw,c("pval","p_value")))),
                       stringsAsFactors=FALSE)
  if(all(is.na(all_lr$interaction_name))) all_lr$interaction_name <- paste(all_lr$ligand,all_lr$receptor,sep="_")
  if(anyNA(all_lr$source) || anyNA(all_lr$target) || all(is.na(all_lr$probability))) fail("CellChat output columns incomplete")
  write_tsv(all_lr,"01_full_7dpi_all_cellchat_interactions.tsv")
  # Compare frozen Arm_2 against Arm_1 without selecting ligands from Step7G.
  recipient_groups <- setdiff(groups,EXCLUDE_RECIPIENTS)
  a <- all_lr[all_lr$source==FORMAL_SENDER & all_lr$target %in% recipient_groups,,drop=FALSE]
  b <- all_lr[all_lr$source==COMPARATOR & all_lr$target %in% recipient_groups,,drop=FALSE]
  if(!nrow(a)) fail("No Arm_2 -> recipient LR rows; inspect all-interactions output")
  k <- c("target","ligand","receptor","interaction_name","pathway_name")
  dta <- data.table::as.data.table(a)
  dtb <- data.table::as.data.table(b)
  A <- dta[,.(probability_Arm_2=max(probability,na.rm=TRUE),
              p_value_Arm_2=min(p_value,na.rm=TRUE)),by=k]
  B <- dtb[,.(probability_Arm_1=max(probability,na.rm=TRUE)),by=k]
  C <- merge(A,B,by=k,all.x=TRUE,sort=FALSE)
  C[is.na(probability_Arm_1),probability_Arm_1:=0]
  C[,delta_Arm2_minus_Arm1:=probability_Arm_2-probability_Arm_1]
  C[,branch_ratio:=(probability_Arm_2+1e-9)/(probability_Arm_1+1e-9)]
  C[,cellchat_nominal_support:=is.finite(p_value_Arm_2) & p_value_Arm_2<=CELLCHAT_P_THRESHOLD &
      is.finite(probability_Arm_2) & probability_Arm_2>0]
  # Biological-sample audit: uniform ligand and receptor detection in each 7dpi
  # sample. Complexes require all component genes to be detected. Unknown
  # complex names -> NA -> cannot qualify as two-sample supported.
  gene_names <- rownames(expr)
  gn_low <- tolower(gene_names)
  subunits <- function(symbol) {
    symbols <- unique(trimws(unlist(strsplit(as.character(symbol),"[_+|]"))))
    match_idx <- match(tolower(symbols),gn_low)
    if(!length(symbols) || anyNA(match_idx)) return(character(0))
    gene_names[match_idx]
  }
  samples <- sort(unique(m$sample_id))
  combos <- unique(m[c("sample_id","cellchat_group")])
  keys <- paste(combos$sample_id,combos$cellchat_group,sep="___")
  fractions <- vector("list",length(keys));names(fractions)<-keys
  n_by_key <- setNames(integer(length(keys)),keys)
  for(i in seq_len(nrow(combos))) {
    ids <- m$cell[m$sample_id==combos$sample_id[i] & m$cellchat_group==combos$cellchat_group[i]]
    n_by_key[i] <- length(ids)
    fractions[[i]] <- if(length(ids)) Matrix::rowSums(expr[,ids,drop=FALSE]>0)/length(ids) else numeric(0)
    names(fractions[[i]]) <- gene_names
  }
  gene_detection <- function(symbol,sample,group) {
    genes <- subunits(symbol)
    key <- paste(sample,group,sep="___")
    if(!length(genes) || !key %in% names(fractions) || n_by_key[[key]] < MIN_SAMPLE_GROUP_CELLS) return(NA_real_)
    min(as.numeric(fractions[[key]][genes]))
  }
  C <- as.data.frame(C,stringsAsFactors=FALSE)
  for(s in samples) {
    suffix <- gsub("[^A-Za-z0-9]","_",s)
    C[[paste0("sender_ligand_fraction_",suffix)]] <- vapply(C$ligand,gene_detection,numeric(1),sample=s,group=FORMAL_SENDER)
    C[[paste0("recipient_receptor_fraction_",suffix)]] <- mapply(gene_detection,C$receptor,MoreArgs=list(sample=s),group=C$target)
    # mapply above sends each target with its corresponding receptor.
    C[[paste0("sample_support_",suffix)]] <- !is.na(C[[paste0("sender_ligand_fraction_",suffix)]]) &
      !is.na(C[[paste0("recipient_receptor_fraction_",suffix)]]) &
      C[[paste0("sender_ligand_fraction_",suffix)]]>=MIN_SENDER_DETECTION &
      C[[paste0("recipient_receptor_fraction_",suffix)]]>=MIN_RECEPTOR_DETECTION
  }
  support_cols <- paste0("sample_support_",gsub("[^A-Za-z0-9]","_",samples))
  C$n_sample_supported <- rowSums(as.data.frame(C[support_cols]))
  # Candidate gate: CellChat P values are descriptive, not replicate-level P.
  C$unfiltered_candidate <- C$cellchat_nominal_support & C$n_sample_supported==2L
  # An explanatory/lexicographic rank, NOT mechanistic validation:
  # first the number of supported biological samples, then CellChat
  # significance status, then model probability (raw scale).
  C <- C[order(-as.integer(C$unfiltered_candidate),-C$n_sample_supported,
               -as.integer(C$cellchat_nominal_support),-C$probability_Arm_2,
               C$interaction_name,C$target),,drop=FALSE]
  C$display_rank <- seq_len(nrow(C))
  C$rank_among_two_sample_candidates <- NA_integer_
  candidate_indices <- which(C$unfiltered_candidate)
  C$rank_among_two_sample_candidates[candidate_indices] <- seq_along(candidate_indices)
  if(file.exists(OLD18)) {
    old <- read_tsv(OLD18)
    kk <- c("target","ligand","receptor","interaction_name")
    if(all(kk %in% names(old))) {
      old$old_key <- do.call(paste,c(old[kk],sep="___"))
      C$in_old_18 <- do.call(paste,c(C[kk],sep="___")) %in% old$old_key
    } else C$in_old_18 <- NA
  } else C$in_old_18 <- NA
  write_tsv(C,"02_all_Arm2_to_recipients_unfiltered_ranked.tsv")
  tnf <- tolower(C$ligand)=="tnf" & grepl("tnfrsf1a",C$receptor,ignore.case=TRUE)
  write_tsv(C[tnf,,drop=FALSE],"03_Tnf_Tnfrsf1a_7dpi_unfiltered_details.tsv")
  report <- c(
    "GSE162610 | Step8F | Complete 7dpi full-cell all-database LR re-screen",
    paste("CellChat version:",as.character(utils::packageVersion("CellChat"))),
    paste("Eligible 7dpi cells:",nrow(m),"(all; no downsampling)"),
    paste("Biological samples:",paste(samples,collapse=", ")),
    paste("CellChatDB interactions:",length(expected_interactions)),
    paste("CellChat tested expression-qualified LR:",nrow(cc@LR$LRsig)),
    paste("Computational bootstrap count:",NBOOT),
    paste("Arm2->recipient rows (all CellChat LR):",nrow(C)),
    paste("Arm2->recipient CellChat nominal support:",sum(C$cellchat_nominal_support)),
    paste("2-sample ligand+receptor detection with CellChat support:",sum(C$unfiltered_candidate)),
    paste("Tnf-Tnfrsf1a rows:",sum(tnf)),
    if(any(tnf)) paste("Best Tnf-Tnfrsf1a candidate rank:",
                       if(any(C$unfiltered_candidate[tnf])) min(C$rank_among_two_sample_candidates[tnf],na.rm=TRUE) else "none qualifies") else "Tnf-Tnfrsf1a not detected",
    "",
    "The analysis did NOT use Step7G TF/regulon/target-program evidence or the original 18/20-axis list for selection.",
    "It does use an expression-based CellChat standard filter and two 7dpi samples to evaluate detection consistency.",
    "Preassigned PRIMARY/CONTEXT labels and manually curated recipient response programs are NOT in the rank.",
    "Comparisons of CellChat communication probability are exploratory; they do not establish biologically independent signaling.",
    "CellChat P values are computational bootstrap statistics, not animal-level P values.",
    "No downstream response-program validation for the full discovery set was performed.",
    "Failure to meet expression detection may reflect transcript dropout, and receptor complex decomposition is conservative.",
    "Please report the FULL candidate distribution, not only Tnf-Tnfrsf1a.",
    paste("Output:",OUT))
  writeLines(report,RUN_REPORT,useBytes=TRUE)
  cat(paste(report,collapse="\n"),"\n")
  capture.output(sessionInfo(),file=file.path(REP,"Step8F_sessionInfo.txt"))
  invisible(C)
}
tryCatch(run(),error=function(e) {
  if(!file.exists(ERROR_FILE)) writeLines(conditionMessage(e),ERROR_FILE)
  stop(e)
})
