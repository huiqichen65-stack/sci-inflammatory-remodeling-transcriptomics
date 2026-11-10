rm(list = ls())
gc()
options(stringsAsFactors = FALSE, width = 220, warn = 1, timeout = 86400)
set.seed(20260801)
PROJECT_ROOT <- Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
CONDA_ENV_NAME <- "pyscenic_env"
NUM_WORKERS <- 10L
GRN_SEED <- 20260801L
RUN_SMOKE_TEST <- TRUE
FORCE_REBUILD_LOOM <- FALSE
FORCE_RERUN_GRN <- FALSE
EXPECTED_CELLS <- 4203L
EXPECTED_GENES <- 7897L
EXPECTED_STATE_COUNTS <- c(Shared_prebranch = 116L, Arm_1 = 2289L, Arm_2 = 1798L)
MIN_TF_OVERLAP <- 900L
MIN_FEATHER_OVERLAP <- 0.90
MIN_EDGES <- 10000L
MIN_TF_COVERAGE <- 0.80
MIN_TARGET_COVERAGE <- 0.70
KEY_TFS <- c("Egr1","Junb","Jun","Fos","Fosb","Atf3","Nfkb1","Nfkb2", "Rela","Relb","Stat1","Stat3","Irf1","Irf7","Spi1","Cebpb", "Mafb","Tfeb","Mitf","Nr1h3","Srebf1","Pparg")
STEP7A_INPUT <- file.path(PROJECT_ROOT, "04_results", "step7A_GSE162610_pySCENIC_input_environment_audit", "input")
MTX_FILE <- file.path(STEP7A_INPUT, "GSE162610_step7A_frozen_4203_raw_counts_genes_by_cells.mtx")
GENE_FILE <- file.path(STEP7A_INPUT, "GSE162610_step7A_frozen_4203_genes.tsv")
BARCODE_FILE <- file.path(STEP7A_INPUT, "GSE162610_step7A_frozen_4203_barcodes.tsv")
METADATA_FILE <- file.path(STEP7A_INPUT, "GSE162610_step7A_frozen_4203_cell_metadata.csv")
RESOURCE_DIR <- file.path(PROJECT_ROOT, "01_raw_data", "pySCENIC_resources")
TF_FILE <- file.path(RESOURCE_DIR, "mm_mgi_tfs.txt")
FEATHER_FILE <- file.path(RESOURCE_DIR, "mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather")
OUTDIR <- file.path(PROJECT_ROOT, "04_results", "step7B_GSE162610_GRNBoost2_network_inference")
INPUTDIR <- file.path(OUTDIR, "input")
NETWORKDIR <- file.path(OUTDIR, "network")
TABDIR <- file.path(OUTDIR, "tables")
LOGDIR <- file.path(OUTDIR, "logs")
REPORTDIR <- file.path(OUTDIR, "reports")
PYDIR <- file.path(OUTDIR, "python")
TMPDIR_GRN <- file.path(OUTDIR, "tmp")
for (d in c(OUTDIR, INPUTDIR, NETWORKDIR, TABDIR, LOGDIR, REPORTDIR, PYDIR, TMPDIR_GRN))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
FULL_LOOM <- file.path(INPUTDIR, "GSE162610_step7B_frozen_4203_raw_counts.loom")
SMOKE_LOOM <- file.path(INPUTDIR, "GSE162610_step7B_smoke_test.loom")
SMOKE_TFS <- file.path(INPUTDIR, "GSE162610_step7B_smoke_test_TFs.txt")
SMOKE_ADJ <- file.path(NETWORKDIR, "GSE162610_step7B_smoke_GRNBoost2_adjacencies.tsv")
ADJ_FILE <- file.path(NETWORKDIR, "GSE162610_step7B_GRNBoost2_adjacencies.tsv")
PRECHECK_FILE <- file.path(TABDIR, "01_step7B_preflight_audit.tsv")
REPORT_FILE <- file.path(REPORTDIR, "Step7B_report.txt")
ERROR_FILE <- file.path(REPORTDIR, "Step7B_ERROR.txt")
SESSION_FILE <- file.path(REPORTDIR, "Step7B_sessionInfo.txt")
if (!requireNamespace("data.table", quietly = TRUE))
  stop("Missing R package: data.table", call. = FALSE)
ts <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")
qcmd <- function(x) shQuote(as.character(x), type = "cmd")
safe <- function(x) normalizePath(x, winslash = "/", mustWork = FALSE)
fmt_time <- function(s) {
  s <- as.numeric(s); if (!is.finite(s)) return(NA_character_)
  if (s < 60) paste0(round(s, 1), " sec") else if (s < 3600) paste0(round(s/60, 1), " min") else paste0(round(s/3600, 2), " h")
}
fmt_size <- function(x) {
  x <- as.numeric(x); u <- c("B","KB","MB","GB","TB"); i <- 1L
  while (is.finite(x) && x >= 1024 && i < length(u)) { x <- x/1024; i <- i+1L }
  paste0(round(x, 3), " ", u[i])
}
write_csv <- function(x, p) utils::write.csv(x, p, row.names = FALSE, fileEncoding = "UTF-8", na = "")
write_tsv <- function(x, p) data.table::fwrite(x, p, sep = "\t", quote = FALSE, na = "")
log_tail <- function(p, n = 120L) {
  if (!file.exists(p)) return("<log not found>")
  z <- readLines(p, warn = FALSE, encoding = "UTF-8"); paste(utils::tail(z, n), collapse = "\n")
}
abort_step <- function(x) {
  writeLines(c("GSE162610 STEP 7B ERROR", paste0("Time: ", ts()), "", x), ERROR_FILE, useBytes = TRUE)
  stop(x, call. = FALSE)
}
run_cmd <- function(command, args, log_file, env = character()) {
  if (file.exists(log_file)) file.remove(log_file)
  env_names <- character(0)
  old_values <- character(0)
  if (length(env) > 0L) {
    env_names <- sub("=.*$", "", env)
    env_values <- sub("^[^=]*=", "", env)
    old_values <- Sys.getenv(env_names, unset = NA_character_)
    do.call(Sys.setenv, as.list(stats::setNames(env_values, env_names)))
  }
  on.exit({
    if (length(env_names) > 0L) {
      unset_names <- env_names[is.na(old_values)]
      if (length(unset_names) > 0L) Sys.unsetenv(unset_names)
      restore_index <- !is.na(old_values)
      if (any(restore_index)) {
        do.call( Sys.setenv, as.list(stats::setNames(old_values[restore_index], env_names[restore_index])) )
      }
    }
  }, add = TRUE)
  st <- Sys.time()
  code <- tryCatch( system2(command, args = args, stdout = log_file, stderr = log_file), error = function(e) {
      999L
    }
  )
  if (is.null(code)) code <- 0L
  list(status = as.integer(code), elapsed = as.numeric(difftime(Sys.time(), st, units = "secs")), log = log_file)
}
valid_adj <- function(p) {
  if (!file.exists(p) || file.info(p)$size < 1024) return(FALSE)
  h <- tryCatch(data.table::fread(p, nrows = 5L, showProgress = FALSE), error = function(e) NULL)
  !is.null(h) && all(c("TF","target","importance") %in% names(h))
}
start_all <- Sys.time()
required <- c(MTX_FILE, GENE_FILE, BARCODE_FILE, METADATA_FILE, TF_FILE, FEATHER_FILE)
if (any(!file.exists(required))) abort_step(paste("Missing files:", paste(required[!file.exists(required)], collapse = "\n")))
python_candidates <- c( file.path(Sys.getenv("USERPROFILE"), "anaconda3", "envs", CONDA_ENV_NAME, "python.exe"), file.path(Sys.getenv("USERPROFILE"), "miniconda3", "envs", CONDA_ENV_NAME, "python.exe") )
ENV_PYTHON <- python_candidates[file.exists(python_candidates)][1]
if (is.na(ENV_PYTHON)) abort_step("Cannot find the validated pyscenic_env Python executable.")
PYSCENIC_EXE <- file.path(dirname(ENV_PYTHON), "Scripts", "pyscenic.exe")
if (!file.exists(PYSCENIC_EXE)) abort_step(paste("Cannot find pySCENIC CLI:", PYSCENIC_EXE))
ENV_PYTHON <- safe(ENV_PYTHON); PYSCENIC_EXE <- safe(PYSCENIC_EXE)
grn_env <- c(paste0("TEMP=", safe(TMPDIR_GRN)), paste0("TMP=", safe(TMPDIR_GRN)), paste0("TMPDIR=", safe(TMPDIR_GRN)), "PYTHONUTF8=1", "OMP_NUM_THREADS=1", "OPENBLAS_NUM_THREADS=1", "MKL_NUM_THREADS=1")
PREP_PY <- file.path(PYDIR, "35_prepare_GRNBoost2_input.py")
py <- c( "import sys", "from pathlib import Path", "import numpy as np", "import pandas as pd", "import loompy", "from scipy import sparse", "from scipy.io import mmread", "from ctxcore.rnkdb import FeatherRankingDatabase", "mtx,gf,bf,mf,tf,ff,loom,smoke,smoketf,audit = map(Path, sys.argv[1:11])", "force = sys.argv[11].lower() == 'true'", "eg=int(sys.argv[12]); ec=int(sys.argv[13])", "def onecol(p):", " return [x.strip().split('\\t')[0] for x in p.open(encoding='utf-8-sig') if x.strip()]", "genes=onecol(gf); cells=onecol(bf); tfs=onecol(tf); md=pd.read_csv(mf, low_memory=False)", "req=['cell','step5CR_state','sample_id','timepoint','retained_step5BR_pseudotime','Y29_signed_progress']", "miss=[x for x in req if x not in md.columns]", "assert not miss, f'Missing metadata: {miss}'", "assert len(genes)==eg and len(cells)==ec", "assert len(set(genes))==len(genes) and len(set(cells))==len(cells)", "assert md.shape[0]==len(cells) and md['cell'].astype(str).tolist()==cells", "states=md['step5CR_state'].astype(str).value_counts().to_dict()", "assert states=={'Shared_prebranch':116,'Arm_1':2289,'Arm_2':1798}, states", "x=mmread(str(mtx)); x=sparse.csr_matrix(x, dtype=np.float32)", "assert x.shape==(len(genes),len(cells)) and x.nnz>0", "assert np.all(np.isfinite(x.data)) and np.all(x.data>=0)", "tfset=set(tfs); overlap=[g for g in genes if g in tfset]", "db=FeatherRankingDatabase(fname=str(ff),name='mm10_10kb'); dbg=set(db.genes); dbo=[g for g in genes if g in dbg]", "rebuild=force or not loom.exists()", "if not rebuild:", " try:", " with loompy.connect(str(loom),'r') as ds:", " rebuild=tuple(ds.shape)!=tuple(x.shape) or list(ds.ra.Gene.astype(str))!=genes or list(ds.ca.CellID.astype(str))!=cells", " except Exception:", " rebuild=True", "if rebuild:", " if loom.exists(): loom.unlink()", " ca={'CellID':np.asarray(cells,dtype=object),'BranchState':np.asarray(md['step5CR_state'].astype(str),dtype=object),", " 'SampleID':np.asarray(md['sample_id'].astype(str),dtype=object),'Timepoint':np.asarray(md['timepoint'].astype(str),dtype=object),", " 'Pseudotime':pd.to_numeric(md['retained_step5BR_pseudotime'],errors='coerce').to_numpy(dtype=float),", " 'Y29SignedProgress':pd.to_numeric(md['Y29_signed_progress'],errors='coerce').to_numpy(dtype=float)}", " loompy.create(str(loom),x,{'Gene':np.asarray(genes,dtype=object)},ca)", "# smoke subset: all prebranch + 100 cells from each arm; 1000 variable genes + 60 variable TFs", "idx=[]", "for s,n in [('Shared_prebranch',116),('Arm_1',100),('Arm_2',100)]: idx += np.flatnonzero(md['step5CR_state'].astype(str).to_numpy()==s)[:n].tolist()", "idx=np.asarray(sorted(set(idx)),dtype=int)", "mu=np.asarray(x.mean(axis=1)).ravel(); va=np.maximum(np.asarray(x.power(2).mean(axis=1)).ravel()-mu**2,0)", "top=np.argsort(va)[-1000:]; tfi=np.asarray([i for i,g in enumerate(genes) if g in tfset],dtype=int); stfi=tfi[np.argsort(va[tfi])[-60:]]", "gi=np.asarray(sorted(set(top.tolist()+stfi.tolist())),dtype=int); sx=x[gi,:][:,idx].tocsr(); sg=[genes[i] for i in gi]; sc=[cells[i] for i in idx]; st=[g for g in sg if g in tfset]", "assert len(st)>=20", "if smoke.exists(): smoke.unlink()", "loompy.create(str(smoke),sx,{'Gene':np.asarray(sg,dtype=object)},{'CellID':np.asarray(sc,dtype=object)})", "smoketf.write_text('\\n'.join(st)+'\\n',encoding='utf-8')", "metrics=[('matrix_genes',len(genes)),('matrix_cells',len(cells)),('matrix_nnz',x.nnz),('TF_overlap',len(overlap)),", "('feather_overlap',len(dbo)),('feather_overlap_fraction',len(dbo)/len(genes)),('loom_rebuilt',rebuild),", "('smoke_genes',len(sg)),('smoke_cells',len(sc)),('smoke_TFs',len(st)),", "('state_Shared_prebranch',states['Shared_prebranch']),('state_Arm_1',states['Arm_1']),('state_Arm_2',states['Arm_2'])]", "with audit.open('w',encoding='utf-8') as h:", " h.write('metric\\tvalue\\n')", " [h.write(f'{k}\\t{v}\\n') for k,v in metrics]", "print('STEP7B_PREFLIGHT_PASS'); [print(k,v) for k,v in metrics]" )
writeLines(py, PREP_PY, useBytes = TRUE)
prep_log <- file.path(LOGDIR, "01_prepare_sparse_loom_and_preflight.log")
prep <- run_cmd(ENV_PYTHON, c(qcmd(PREP_PY), qcmd(MTX_FILE), qcmd(GENE_FILE), qcmd(BARCODE_FILE), qcmd(METADATA_FILE), qcmd(TF_FILE), qcmd(FEATHER_FILE), qcmd(FULL_LOOM), qcmd(SMOKE_LOOM), qcmd(SMOKE_TFS), qcmd(PRECHECK_FILE), tolower(as.character(FORCE_REBUILD_LOOM)), EXPECTED_GENES, EXPECTED_CELLS), prep_log, grn_env)
if (prep$status != 0L) abort_step(paste("Preflight failed:\n", log_tail(prep_log, 160)))
pre <- data.table::fread(PRECHECK_FILE, data.table = FALSE)
pm <- setNames(as.character(pre$value), pre$metric)
pre_pass <- as.integer(pm['matrix_genes'])==EXPECTED_GENES && as.integer(pm['matrix_cells'])==EXPECTED_CELLS && as.integer(pm['TF_overlap'])>=MIN_TF_OVERLAP && as.numeric(pm['feather_overlap_fraction'])>=MIN_FEATHER_OVERLAP && as.integer(pm['state_Shared_prebranch'])==116L && as.integer(pm['state_Arm_1'])==2289L && as.integer(pm['state_Arm_2'])==1798L
if (!pre_pass) abort_step("Step 7B preflight thresholds were not met.")
write_csv(pre, file.path(TABDIR, "01_step7B_preflight_audit.csv"))
smoke_pass <- TRUE; smoke_edges <- NA_integer_; smoke_tfs_n <- NA_integer_; smoke_elapsed <- NA_real_
if (RUN_SMOKE_TEST) {
  if (file.exists(SMOKE_ADJ)) file.remove(SMOKE_ADJ)
  smoke_log <- file.path(LOGDIR, "02_GRNBoost2_smoke_test.log")
  sr <- run_cmd(PYSCENIC_EXE, c("grn","--num_workers","2","--client_or_address","local","--sparse","--method","grnboost2", "--seed",GRN_SEED,"--cell_id_attribute","CellID","--gene_attribute","Gene","--output",qcmd(SMOKE_ADJ), qcmd(SMOKE_LOOM),qcmd(SMOKE_TFS)), smoke_log, grn_env)
  smoke_elapsed <- sr$elapsed
  if (sr$status==0L && valid_adj(SMOKE_ADJ)) {
    sa <- data.table::fread(SMOKE_ADJ, showProgress = FALSE); smoke_edges <- nrow(sa); smoke_tfs_n <- data.table::uniqueN(sa$TF)
    smoke_pass <- smoke_edges>=1000L && smoke_tfs_n>=20L; rm(sa); gc()
  } else smoke_pass <- FALSE
  write_csv(data.frame(metric=c("status","elapsed_seconds","edges","unique_TFs","pass"), value=c(sr$status,sr$elapsed,smoke_edges,smoke_tfs_n,smoke_pass)), file.path(TABDIR,"02_smoke_test_audit.csv"))
  if (!smoke_pass) abort_step(paste("Smoke test failed; full run was not started.\n",log_tail(smoke_log,160)))
}
gr_status <- 0L; gr_elapsed <- 0; reused <- FALSE
if (!FORCE_RERUN_GRN && valid_adj(ADJ_FILE)) {
  reused <- TRUE; msg("Reusing the existing structurally valid full adjacency file.")
} else {
  if (file.exists(ADJ_FILE)) file.remove(ADJ_FILE)
  grn_log <- file.path(LOGDIR, "03_GRNBoost2_full_run.log")
  gr <- run_cmd(PYSCENIC_EXE, c("grn","--num_workers",NUM_WORKERS,"--client_or_address","local","--sparse","--method","grnboost2", "--seed",GRN_SEED,"--cell_id_attribute","CellID","--gene_attribute","Gene","--output",qcmd(ADJ_FILE), qcmd(FULL_LOOM),qcmd(TF_FILE)), grn_log, grn_env)
  gr_status <- gr$status; gr_elapsed <- gr$elapsed
  if (gr_status!=0L || !valid_adj(ADJ_FILE)) abort_step(paste("Full GRNBoost2 failed.\n",log_tail(grn_log,220)))
}
adj <- data.table::fread(ADJ_FILE, showProgress = TRUE)
if (!all(c("TF","target","importance") %in% names(adj))) abort_step("Adjacency columns are invalid.")
adj <- adj[, .(TF=as.character(TF), target=as.character(target), importance=as.numeric(importance))]
data.table::setorder(adj, -importance)
n_edges <- nrow(adj); n_tfs <- data.table::uniqueN(adj$TF); n_targets <- data.table::uniqueN(adj$target)
input_tfs <- as.integer(pm['TF_overlap']); input_genes <- as.integer(pm['matrix_genes'])
tf_cov <- n_tfs/input_tfs; target_cov <- n_targets/input_genes
bad_imp <- sum(!is.finite(adj$importance) | adj$importance<=0)
dups <- nrow(adj[, .N, by=.(TF,target)][N>1L]); self_edges <- nrow(adj[TF==target])
qs <- stats::quantile(adj$importance, c(0,.25,.5,.75,.9,.95,.99,.999,1), na.rm=TRUE)
summary <- data.frame(metric=c("method","seed","workers","reused_existing","status","elapsed_seconds","elapsed_human", "adjacency_size_bytes","adjacency_size_human","input_genes","input_cells","input_TFs","network_edges","unique_TFs", "unique_targets","TF_coverage_fraction","target_coverage_fraction","bad_importance_values","duplicate_pairs","self_edges", "median_importance","p99_importance","maximum_importance"), value=c("GRNBoost2",GRN_SEED,NUM_WORKERS,reused,gr_status,gr_elapsed,fmt_time(gr_elapsed),file.info(ADJ_FILE)$size, fmt_size(file.info(ADJ_FILE)$size),input_genes,EXPECTED_CELLS,input_tfs,n_edges,n_tfs,n_targets,tf_cov,target_cov,bad_imp, dups,self_edges,unname(qs['50%']),unname(qs['99%']),unname(qs['100%'])), stringsAsFactors=FALSE)
write_csv(summary, file.path(TABDIR,"03_GRNBoost2_network_summary.csv"))
write_csv(data.frame(quantile=names(qs),importance=as.numeric(qs)), file.path(TABDIR,"04_importance_quantiles.csv"))
tf_summary <- adj[, .(n_targets=.N,total_importance=sum(importance),mean_importance=mean(importance), median_importance=median(importance),maximum_importance=max(importance),top_target=target[1], top_target_importance=importance[1]), by=TF]
data.table::setorder(tf_summary,-total_importance,-maximum_importance)
tf_summary[, total_importance_rank := seq_len(.N)]
write_tsv(tf_summary,file.path(TABDIR,"05_TF_network_summary.tsv"))
write_tsv(head(adj,2000L),file.path(TABDIR,"06_top2000_GRNBoost2_edges.tsv"))
key <- merge(data.table::data.table(TF=KEY_TFS),tf_summary,by="TF",all.x=TRUE,sort=FALSE)
key[,present_in_network:=!is.na(n_targets)]
data.table::setcolorder(key,c("TF","present_in_network",setdiff(names(key),c("TF","present_in_network"))))
write_tsv(key,file.path(TABDIR,"07_key_TF_network_audit.tsv"))
key_targets <- adj[TF %in% KEY_TFS][order(TF,-importance),head(.SD,10L),by=TF]
write_tsv(key_targets,file.path(TABDIR,"08_key_TF_top10_targets.tsv"))
write_tsv(head(tf_summary,100L),file.path(TABDIR,"09_top100_TFs_by_total_importance.tsv"))
network_pass <- gr_status==0L && n_edges>=MIN_EDGES && tf_cov>=MIN_TF_COVERAGE && target_cov>=MIN_TARGET_COVERAGE && bad_imp==0L && dups==0L
decision <- data.frame(domain=c("Input files","Sparse loom preflight","Smoke test","Full GRNBoost2","Minimum edge count", "TF coverage","Target coverage","Positive finite importance","No duplicate TF-target pairs","Overall Step 7B"), pass=c(all(file.exists(required)),pre_pass,smoke_pass,gr_status==0L && valid_adj(ADJ_FILE),n_edges>=MIN_EDGES, tf_cov>=MIN_TF_COVERAGE,target_cov>=MIN_TARGET_COVERAGE,bad_imp==0L,dups==0L,pre_pass&&smoke_pass&&network_pass), stringsAsFactors=FALSE)
decision$status <- ifelse(decision$pass,"PASS","FAIL")
write_csv(decision,file.path(TABDIR,"10_final_decision.csv"))
FINAL_STATUS <- if (tail(decision$pass,1)) "PASS" else "FAIL"
result <- list(step="Step7B",generated_at=ts(),method="GRNBoost2",seed=GRN_SEED,workers=NUM_WORKERS, input_genes=input_genes,input_cells=EXPECTED_CELLS,input_TFs=input_tfs,full_loom=safe(FULL_LOOM), adjacency_file=safe(ADJ_FILE),network_summary=summary,TF_summary=as.data.frame(tf_summary), key_TF_summary=as.data.frame(key),key_TF_top_targets=as.data.frame(key_targets),preflight=pre,decision=decision)
saveRDS(result,file.path(OUTDIR,"GSE162610_step7B_GRNBoost2_results.rds"),compress=TRUE)
elapsed_all <- as.numeric(difftime(Sys.time(),start_all,units="secs"))
writeLines(capture.output(sessionInfo()),SESSION_FILE,useBytes=TRUE)
rm(adj,tf_summary,key,key_targets); gc()
