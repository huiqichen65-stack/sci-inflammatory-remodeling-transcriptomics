options(stringsAsFactors = FALSE, warn = 1)
PROJECT_ROOT <- Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
STEP11B <- file.path(PROJECT_ROOT, "04_results", "step11B_GSE234774_refined_main_Arm_specific_shared_AUCell_ssGSEA")
INPUT_RESULTS <- file.path(STEP11B, "RDS", "GSE234774_step11B_refined_Arm_AUCell_ssGSEA_results.rds")
INPUT_SCORED <- file.path(STEP11B, "RDS", "GSE234774_step11B_refined_main_Arm_scored_BPCells.rds")
OUT <- file.path(PROJECT_ROOT, "04_results", "step11D2_GSE234774_target_category_matched_null_sensitivity")
PREVIOUS_STEP11D_SUMMARY <- file.path(PROJECT_ROOT, "04_results", "step11D_GSE234774_target_size_expression_matched_null_benchmark", "tables", "06_category_matched_target_vs_null_timepoint_summary.csv")
TBL <- file.path(OUT, "tables")
FIG <- file.path(OUT, "figures")
dir.create(TBL, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
LOG <- file.path(OUT, "00_run_log.txt")
if (file.exists(LOG)) file.remove(LOG)
logmsg <- function(...) {
  x <- paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(..., collapse = " "))
  cat(x, "\n")
  cat(x, "\n", file = LOG, append = TRUE)
}
stopf <- function(...) stop(paste0(..., collapse = ""), call. = FALSE)
savecsv <- function(x, basename) {
  utils::write.csv(x, file.path(TBL, basename), row.names = FALSE, na = "")
}
SEED_MATCHING <- 20261028L
SEED_AUC_RANKINGS <- 20260803L
N_NULL <- 500L
NEAREST_POOL <- 500L
RANDOM_TOP_K <- 100L
EXPECTED_NUCLEI <- 19617L
EXPECTED_LIBRARIES <- 13L
EXPECTED_TARGET_GENES <- 496L
MAIN_TIMEPOINTS <- c("Uninjured", "7 dpi", "14 dpi", "1 month", "2 months")
PROGRAM <- "Arm_2_specific"
needed <- c("SeuratObject", "Matrix", "AUCell", "BiocParallel")
missing_packages <- needed[ !vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stopf("Required packages unavailable: ", paste(missing_packages, collapse=", "))
}
for (p in c(INPUT_RESULTS, INPUT_SCORED)) {
  if (!file.exists(p)) stopf("Missing frozen Step11B input: ", p)
}
logmsg("Reading Step11B results and scored Seurat object (read only).")
result <- readRDS(INPUT_RESULTS)
required_slots <- c("scored_gene_sets_present", "AUCell_aucMaxRank", "cell_AUCell")
if (!all(required_slots %in% names(result)))
  stopf("Step11B results RDS missing fields: ", paste(setdiff(required_slots, names(result)), collapse = ", "))
obj <- readRDS(INPUT_SCORED)
if (!inherits(obj, "Seurat")) stopf("Scored object is not a Seurat object.")
md <- obj@meta.data
if (!all(c("biological_library", "timepoint") %in% colnames(md)))
  stopf("Seurat metadata lacks biological_library or timepoint.")
cell_ids <- colnames(obj)
if (length(cell_ids) != EXPECTED_NUCLEI || nrow(md) != EXPECTED_NUCLEI)
  stopf("Expected exactly 19,617 nuclei; found ", length(cell_ids), ".")
if (!setequal(rownames(md), cell_ids)) stopf("Metadata/cell barcode mismatch.")
md <- md[cell_ids, , drop = FALSE]
lib <- as.character(md$biological_library)
tp <- as.character(md$timepoint)
if (anyNA(lib) || anyNA(tp) || !all(tp %in% MAIN_TIMEPOINTS))
  stopf("Missing/unexpected library or timepoint metadata.")
library_to_timepoint <- tapply(tp, lib, function(x) {
  z <- unique(x)
  if (length(z) != 1L) stopf("Mixed timepoints in one library.")
  z
})
if (length(library_to_timepoint) != EXPECTED_LIBRARIES)
  stopf("Expected exactly 13 libraries; found ", length(library_to_timepoint), ".")
expected_libraries_by_time <- c("Uninjured"=2L, "7 dpi"=2L, "14 dpi"=3L, "1 month"=3L, "2 months"=3L)
actual_libraries_by_time <- table(factor(unname(library_to_timepoint), levels = MAIN_TIMEPOINTS))
if (!identical(as.integer(actual_libraries_by_time), as.integer(expected_libraries_by_time)))
  stopf("Unexpected libraries per timepoint.")
libraries <- names(library_to_timepoint)
uninj_libs <- libraries[unname(library_to_timepoint) == "Uninjured"]
uninj_cells <- cell_ids[lib %in% uninj_libs]
gene_sets <- result$scored_gene_sets_present
if (!PROGRAM %in% names(gene_sets)) stopf("Target gene set not available.")
target <- unique(as.character(gene_sets[[PROGRAM]]))
if (length(target) != EXPECTED_TARGET_GENES) {
  stopf("Expected 496 mapped target genes, observed ", length(target), ".")
}
if (anyNA(target) || any(!nzchar(target)))
  stopf("Invalid target gene IDs.")
rna_layers <- SeuratObject::Layers(obj[["RNA"]])
count_layers <- rna_layers[grepl("^counts($|\\.)", rna_layers)]
if (length(count_layers) != 1L)
  stopf("Expected exactly one RNA counts layer, found: ", paste(count_layers, collapse = ", "))
logmsg("Reading sparse RNA counts layer: ", count_layers[[1L]])
raw_matrix <- SeuratObject::LayerData( object=obj, assay="RNA", layer=count_layers[[1L]], fast=FALSE)
if (inherits(raw_matrix, "dgCMatrix")) {
  counts <- raw_matrix
} else {
  counts <- tryCatch(as(raw_matrix, "dgCMatrix"), error=function(e) stopf( "Cannot coerce counts to dgCMatrix. ", "Inspect BPCells paths and available RAM: ", conditionMessage(e)))
}
rm(raw_matrix); invisible(gc())
if (!identical(colnames(counts), cell_ids))
  stopf("Counts columns do not match scored-object cell order.")
if (anyDuplicated(rownames(counts)) || anyDuplicated(colnames(counts)))
  stopf("Duplicate gene or cell IDs.")
features <- rownames(counts)
if (!all(target %in% features))
  stopf("Frozen target genes not found in counts: ", paste(head(setdiff(target, features), 10), collapse = ", "))
if (any(Matrix::colSums(counts) == 0))
  stopf("At least one retained nucleus has zero total counts.")
logmsg("Computing baseline expression/detection directly from RNA counts in ", length(uninj_libs), " uninjured biological libraries.")
baseline_by_library <- lapply(uninj_libs, function(l) {
  selected_cells <- which(lib == l)
  if (!length(selected_cells))
    stopf("Missing uninjured library in counts: ", l)
  subcounts <- counts[, selected_cells, drop=FALSE]
  cell_depth <- as.numeric(Matrix::colSums(subcounts))
  if (any(!is.finite(cell_depth)) || any(cell_depth <= 0))
    stopf("Nonpositive/nonfinite RNA library size in uninjured library: ", l)
  normalized <- subcounts %*%
    Matrix::Diagonal(x = 10000 / cell_depth)
  expression <- log1p(as.numeric(Matrix::rowMeans(normalized)))
  detection <- as.numeric(Matrix::rowMeans(subcounts > 0))
  rm(subcounts, normalized)
  list(library=l, expression=expression, detection=detection, n_nuclei=length(selected_cells))
})
baseline_expression <- rowMeans(do.call(cbind, lapply(baseline_by_library, function(z) z$expression)))
baseline_detection <- rowMeans(do.call(cbind, lapply(baseline_by_library, function(z) z$detection)))
names(baseline_expression) <- features
names(baseline_detection) <- features
if (any(!is.finite(baseline_expression)) || any(!is.finite(baseline_detection)) || any(baseline_detection < 0 | baseline_detection > 1))
  stopf("Invalid baseline expression/detection metrics.")
savecsv(data.frame( biological_library=vapply(baseline_by_library, `[[`, "", "library"), n_nuclei=vapply(baseline_by_library, `[[`, integer(1), "n_nuclei"), timepoint="Uninjured" ), "01_uninjured_library_matching_audit.csv")
logmsg("Baseline matching PASS: ", length(features), " genes; ", paste(vapply(baseline_by_library, function(z) paste0(z$library, "=", z$n_nuclei), ""), collapse=", "), " uninjured nuclei per library.")
rm(baseline_by_library)
invisible(gc())
n_genes <- length(features)
scaled_rank <- function(v) (rank(v, ties.method="average") - 0.5) / length(v)
x <- scaled_rank(baseline_expression)
y <- scaled_rank(baseline_detection)
names(x) <- names(y) <- features
category_levels <- c( "Mitochondrial", "Ribosomal", "Hemoglobin", "Gm_numeric", "Other" )
classify_symbol <- function(v) {
  z <- rep("Other", length(v))
  z[grepl("^Gm[0-9]+", v, ignore.case=TRUE)] <- "Gm_numeric"
  z[grepl("^Hb[ab]", v, ignore.case=TRUE)] <- "Hemoglobin"
  z[grepl("^(Rpl|Rps)", v, ignore.case=TRUE)] <- "Ribosomal"
  z[grepl("^mt-", v, ignore.case=TRUE)] <- "Mitochondrial"
  z
}
background <- setdiff(features, target)
if (length(background) < length(target)*2L) {
  stopf("Too few background genes after target exclusion.")
}
target_category <- stats::setNames(classify_symbol(target), target)
background_category <- stats::setNames(classify_symbol(background), background)
target_class_count <- table(factor(unname(target_category), levels=category_levels))
background_class_count <- table(factor(unname(background_category), levels=category_levels))
availability <- data.frame( category=category_levels, n_target=as.integer(target_class_count), n_available_background=as.integer(background_class_count), adequate_background=as.integer(background_class_count) >= as.integer(target_class_count), stringsAsFactors=FALSE )
savecsv(availability, "01B_target_background_category_availability.csv")
if (!all(availability$adequate_background)) {
  insufficient <- availability$category[!availability$adequate_background]
  stopf("Cannot perform EXACT category matching in: ", paste(insufficient, collapse=", "), ". Inspect 01B_target_background_category_availability.csv. ", "Do NOT silently relax the required category counts.")
}
logmsg("Frozen target category counts: ", paste(availability$category, availability$n_target, sep="=", collapse=", "))
n_neighbor <- min(NEAREST_POOL, length(background))
candidates <- vector("list", length(target))
names(candidates) <- target
for (i in seq_along(target)) {
  g <- target[[i]]
  eligible <- background[background_category == target_category[[g]]]
  d2 <- (x[eligible] - x[[g]])^2 + (y[eligible] - y[[g]])^2
  ids <- order(d2, seq_along(d2))[seq_len(min(n_neighbor, length(eligible)))]
  candidates[[i]] <- eligible[ids]
}
if (any(vapply(candidates, length, integer(1)) == 0L))
  stopf("At least one target gene has no within-category match.")
set.seed(SEED_MATCHING)
null_sets <- vector("list", N_NULL)
names(null_sets) <- sprintf("CategoryMatchedNull_%03d", seq_len(N_NULL))
match_err <- numeric(N_NULL)
match_p90 <- numeric(N_NULL)
match_expression <- numeric(N_NULL)
match_detection <- numeric(N_NULL)
null_category_rows <- vector("list", N_NULL)
for (j in seq_len(N_NULL)) {
  used <- stats::setNames(rep(FALSE, length(background)), background)
  picked <- character(length(target))
  order_target <- order( vapply(candidates, length, integer(1)), stats::runif(length(target)) )
  for (i in order_target) {
    cand <- candidates[[i]]
    remaining <- cand[!used[cand]]
    if (!length(remaining)) {
      g <- target[[i]]
      pool <- background[background_category == target_category[[g]]]
      pool <- pool[!used[pool]]
      if (!length(pool)) stopf( "No remaining within-category control for ", g, " in random set ", j, ". Exact matching impossible.")
      d2 <- (x[pool]-x[[g]])^2 + (y[pool]-y[[g]])^2
      remaining <- pool[order(d2, seq_along(d2))]
    }
    choice <- sample(head(remaining, min(RANDOM_TOP_K, length(remaining))), 1L)
    picked[[i]] <- choice
    used[[choice]] <- TRUE
  }
  if (anyDuplicated(picked) || length(intersect(picked,target)) || length(picked) != EXPECTED_TARGET_GENES) {
    stopf("Null set size/uniqueness/target-exclusion audit FAILED: ", j)
  }
  this_class_count <- table(factor(classify_symbol(picked), levels=category_levels))
  if (!identical(as.integer(this_class_count), as.integer(target_class_count))) {
    stopf("Exact prefix-category balance FAILED for null set ", j)
  }
  null_sets[[j]] <- picked
  distances <- sqrt((x[picked]-x[target])^2 + (y[picked]-y[target])^2)
  match_err[[j]] <- stats::median(distances)
  match_p90[[j]] <- unname(stats::quantile(distances, 0.90))
  match_expression[[j]] <- stats::median(abs(x[picked]-x[target]))
  match_detection[[j]] <- stats::median(abs(y[picked]-y[target]))
  null_category_rows[[j]] <- data.frame( null_set=names(null_sets)[[j]], category=category_levels, target_count=as.integer(target_class_count), matched_count=as.integer(this_class_count), stringsAsFactors=FALSE )
}
sets <- c(list(Target_496=target), null_sets)
category_audit <- do.call(rbind, null_category_rows)
if (!all(category_audit$target_count == category_audit$matched_count))
  stopf("Postsampling exact-category audit FAILED.")
savecsv(category_audit, "03B_exact_category_balance_500_sets.csv")
paired_genes <- data.frame( null_set=rep(names(null_sets), each=length(target)), target_gene=rep(target, times=N_NULL), target_category=rep(unname(target_category), times=N_NULL), matched_gene=unlist(null_sets, use.names=FALSE), matched_category=classify_symbol(unlist(null_sets, use.names=FALSE)), stringsAsFactors=FALSE )
if (!all(paired_genes$target_category == paired_genes$matched_category))
  stopf("Per-gene strict category matching FAILED.")
con <- gzfile(file.path(TBL, "02_category_matched_null_genes_by_target.csv.gz"), "wt")
utils::write.csv(paired_genes, con, row.names=FALSE)
close(con)
savecsv(data.frame( null_set=names(null_sets), n_genes=vapply(null_sets, length, integer(1)), median_rank_distance=match_err, p90_rank_distance=match_p90, median_abs_expression_rank_difference=match_expression, median_abs_detection_rank_difference=match_detection, stringsAsFactors=FALSE ), "03_category_matched_null_expression_detection_QC.csv")
logmsg("Generated ", N_NULL, " size-, prefix-class-, baseline-expression/detection-matched controls.")
logmsg("Building cell rankings using the original Step11B AUCell settings.")
auc_rank <- as.integer(result$AUCell_aucMaxRank)
if (length(auc_rank) != 1L || !is.finite(auc_rank) || auc_rank != ceiling(0.05*nrow(counts)))
  stopf("Unexpected Step11B AUCell aucMaxRank: ", auc_rank)
set.seed(SEED_AUC_RANKINGS)
ranking <- AUCell::AUCell_buildRankings( exprMat=counts, featureType="genes", plotStats=FALSE, splitByBlocks=TRUE, BPPARAM=BiocParallel::SerialParam(), keepZeroesAsNA=FALSE, verbose=TRUE)
rm(counts); invisible(gc())
logmsg("Scoring target plus ", N_NULL, " matched-null signatures.")
auc_obj <- AUCell::AUCell_calcAUC( geneSets=sets, rankings=ranking, nCores=1L, normAUC=TRUE, aucMaxRank=auc_rank, verbose=TRUE)
auc <- as.matrix(AUCell::getAUC(auc_obj))
rm(auc_obj, ranking); invisible(gc())
if (!identical(colnames(auc), cell_ids) || !setequal(rownames(auc), names(sets)) || any(!is.finite(auc)))
  stopf("AUCell score matrix failed name/finite-value audit.")
original_auc <- result$cell_AUCell
if (!PROGRAM %in% rownames(original_auc) || !all(cell_ids %in% colnames(original_auc)))
  stopf("Saved Step11B original target AUCell row/cell IDs unavailable.")
target_difference <- max(abs( as.numeric(auc["Target_496",cell_ids]) - as.numeric(original_auc[PROGRAM,cell_ids])))
if (!is.finite(target_difference) || target_difference > 1e-6)
  stopf("Reconstructed target AUCell is not identical to frozen Step11B: ", "max_abs_difference=", format(target_difference, digits=6), ". Do not interpret the null benchmark before resolving this.")
logmsg("Original AUCell reproduction PASS: maximum absolute difference = ", format(target_difference, digits=7))
logmsg("Summarizing median AUCell per gene set and biological library.")
library_scores <- do.call(rbind, lapply(libraries, function(l) {
  cell_subset <- which(lib == l)
  med <- apply(auc[, cell_subset, drop=FALSE], 1L, stats::median)
  data.frame( gene_set=names(med), biological_library=l, timepoint=unname(library_to_timepoint[[l]]), n_nuclei=length(cell_subset), AUCell_median=as.numeric(med), stringsAsFactors=FALSE )
}))
library_scores$is_target <- library_scores$gene_set == "Target_496"
savecsv(library_scores, "04_category_matched_target_and_null_library_AUCell.csv")
rm(auc); invisible(gc())
contrasts <- do.call(rbind, lapply(names(sets), function(gs) {
  df <- library_scores[library_scores$gene_set == gs, , drop=FALSE]
  base <- mean(df$AUCell_median[df$timepoint == "Uninjured"])
  do.call(rbind, lapply(MAIN_TIMEPOINTS[-1L], function(t) {
    data.frame(gene_set=gs, timepoint=t, n_injured_libraries=sum(df$timepoint==t), n_uninjured_libraries=sum(df$timepoint=="Uninjured"), injured_mean_library_median=mean(df$AUCell_median[df$timepoint==t]), uninjured_mean_library_median=base, injury_minus_uninjured= mean(df$AUCell_median[df$timepoint==t]) - base)
  }))
}))
savecsv(contrasts, "05_category_matched_all_gene_set_timepoint_contrasts.csv")
summary_table <- do.call(rbind, lapply(MAIN_TIMEPOINTS[-1L], function(t) {
  dt <- contrasts[contrasts$timepoint==t, , drop=FALSE]
  obs <- dt$injury_minus_uninjured[dt$gene_set=="Target_496"]
  nul <- dt$injury_minus_uninjured[dt$gene_set!="Target_496"]
  data.frame( timepoint=t, n_injured_libraries=dt$n_injured_libraries[1L], n_uninjured_libraries=dt$n_uninjured_libraries[1L], target_delta_AUCell=obs, null_median_delta_AUCell=stats::median(nul), null_2.5pct=unname(stats::quantile(nul, 0.025)), null_97.5pct=unname(stats::quantile(nul, 0.975)), target_percentile_vs_null=100*mean(nul <= obs), empirical_upper_tail_fraction=(1+sum(nul>=obs))/(N_NULL+1), interpretation_note="Gene-set randomization only; not biological replicate P" )
}))
savecsv(summary_table, "06_category_matched_target_vs_null_timepoint_summary.csv")
pdf(file.path(FIG,"01_category_matched_target_vs_null_library_level_contrasts.pdf"), width=9, height=6.5, useDingbats=FALSE)
op <- par(mar=c(4.5,5.0,2.5,1.5))
null_delta <- lapply(MAIN_TIMEPOINTS[-1L], function(t) {
  v <- contrasts$injury_minus_uninjured[contrasts$timepoint==t & contrasts$gene_set!="Target_496"]
  v[is.finite(v)]
})
names(null_delta) <- MAIN_TIMEPOINTS[-1L]
plot_all <- c(unlist(null_delta, use.names=FALSE), summary_table$target_delta_AUCell)
plot_limits <- range(plot_all[is.finite(plot_all)])
plot_span <- max(diff(plot_limits), 0.002)
plot_limits <- plot_limits + c(-0.09, 0.13)*plot_span
boxplot(null_delta, outline=FALSE, col="grey85", border="grey35", ylim=plot_limits, ylab="Injured minus uninjured (mean of library-median AUCell)", xlab="Time after spinal cord injury", main="Frozen target versus category- and baseline-matched null sets")
abline(h=0, col="grey55", lty=3)
points(seq_along(null_delta), summary_table$target_delta_AUCell, pch=19, col="black", cex=1.5)
legend("topleft",legend=c("Category- and baseline-matched null", "Frozen target program"), fill=c("grey85", NA), border=c("grey35", NA), pch=c(NA,19), bty="n", cex=.85)
par(op); dev.off()
if (file.exists(PREVIOUS_STEP11D_SUMMARY)) {
  previous <- utils::read.csv(PREVIOUS_STEP11D_SUMMARY, stringsAsFactors=FALSE, check.names=FALSE)
  needed_previous <- c("timepoint", "target_delta_AUCell", "null_median_delta_AUCell", "null_97.5pct")
  if (!all(needed_previous %in% colnames(previous)))
    stopf("Old Step11D summary columns have changed; manual audit required.")
  if (anyDuplicated(previous$timepoint) || !setequal(previous$timepoint, MAIN_TIMEPOINTS[-1L]))
    stopf("Previous Step11D summary timepoints mismatch.")
  previous <- previous[match(summary_table$timepoint,previous$timepoint),]
  prior_target_delta <- as.numeric(previous$target_delta_AUCell)
  max_target_delta_difference <- max(abs(prior_target_delta - summary_table$target_delta_AUCell))
  if (!is.finite(max_target_delta_difference) || max_target_delta_difference > 1e-8) {
    stopf("Old and new frozen target timepoint deltas disagree: ", max_target_delta_difference)
  }
  compare_table <- data.frame( timepoint=summary_table$timepoint, frozen_target_delta_AUCell=summary_table$target_delta_AUCell, original_null_median=as.numeric(previous$null_median_delta_AUCell), category_matched_null_median=summary_table$null_median_delta_AUCell, original_null_97.5pct=as.numeric(previous$null_97.5pct), category_matched_null_97.5pct=summary_table$null_97.5pct, old_new_target_delta_max_abs_diff=max_target_delta_difference )
  savecsv(compare_table, "07_original_vs_category_matched_sensitivity.csv")
  logmsg("Original-vs-strict target AUCell delta comparison PASS; max diff=", format(max_target_delta_difference, digits=7))
} else {
  logmsg("Optional old Step11D summary not found; cross-run comparison skipped.")
  max_target_delta_difference <- NA_real_
}
checks <- data.frame( check=c( "frozen_step11B_results_read_only", "nuclei_19617", "biological_libraries_13", "timepoint_library_composition", "frozen_target_gene_count_496", "background_excludes_target", "null_sets_500_size_496_and_distinct_genes", "exact_prefix_category_counts_all_null_sets", "target_matched_gene_pair_category_consistency", "same_AUCell_rank_cutoff_as_Step11B", "original_target_AUCell_reproduced", "library_level_summary_units", "timepoint_contrast_summary_complete" ), status=rep("PASS", 13L), stringsAsFactors=FALSE )
savecsv(checks, "09_STEP11D2_audit_checks.csv")
report_lines <- c( "Step11D2 | Reviewer 2 optional suggestion 1 | Category-matched sensitivity", paste0("Input frozen results: ", INPUT_RESULTS), paste0("Input frozen scored object: ", INPUT_SCORED), paste0("Target program: ", PROGRAM, ", ", length(target), " mapped genes"), paste0("Retained cohort: ", length(libraries), " biological libraries, ", length(cell_ids), " nuclei"), paste0("Nulls: ", N_NULL, " sets of exactly ", length(target), " distinct genes, each excluding all target genes"), paste0("Target symbol category counts: ", paste(availability$category, availability$n_target, sep="=", collapse=", ")), "All 500 sets have EXACT per-category equality to frozen target.", "Matching uses within-category baseline expression/detection of two", "uninjured biological libraries, averaged with equal library weight.", "Categories are gene-symbol prefix strata; not assertions of biology.", paste0("AUC normalized = TRUE; aucMaxRank = ", auc_rank, "; frozen target max_abs_reproduction_diff = ", format(target_difference,digits=8)), "Main descriptive metric = injured minus uninjured difference", "in the mean of library-median AUCell scores, at each injured timepoint.", "Gene-set randomization tail fraction is NOT a biological-replicate P value.", "These findings do NOT establish microglial exclusivity, cell fate", "bifurcation, same-cell persistence, or protein-level mechanism.", paste0("Matching rank-distance median across 500 nulls = ", format(stats::median(match_err),digits=6)), paste0("Matching rank-distance 90th percentile across null medians = ", format(unname(stats::quantile(match_err, 0.90)),digits=6)), paste0("Old-vs-new frozen target delta max absolute diff = ", format(max_target_delta_difference, digits=6)), "", "Timepoint results:", paste(capture.output(print(summary_table, row.names=FALSE)), collapse="\n") )
writeLines(report_lines, file.path(OUT,"10_STEP11D2_FINAL_REPORT.txt"))
capture.output(sessionInfo(), file=file.path(OUT,"00_sessionInfo.txt"))
logmsg("Step11D2 finished successfully. Original Step11B and Step11D unchanged.")
