options(stringsAsFactors = FALSE, warn = 1, width = 160)
PROJECT_ROOT <- Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
if (!dir.exists(PROJECT_ROOT)) {
  stop("Project root was not found: ", PROJECT_ROOT)
}
if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("Required R package missing: data.table. Install it first.")
}
STEP12D_DIR <- file.path(PROJECT_ROOT, "04_results", "step12D_GSE234774_full_mode_RCTD_per_capture_area")
STEP12F_DIR <- file.path(PROJECT_ROOT, "04_results", "step12F_GSE234774_independent_lesion_zones_target_program_validation")
rctd_file <- file.path(STEP12D_DIR, "tables", "03_all_capture_areas_normalized_weights.csv.gz")
zone_file <- file.path(STEP12F_DIR, "tables", "02_spot_independent_lesion_zones_and_scores.csv.gz")
original_file <- file.path(STEP12F_DIR, "tables", "05_animal_zone_vs_remote_contrasts.csv")
output_dir <- file.path(PROJECT_ROOT, "04_results", "step12H_GSE234774_spatial_cell_source_sensitivity")
tables_dir <- file.path(output_dir, "tables")
figures_dir <- file.path(output_dir, "figures")
for (d in c(output_dir, tables_dir, figures_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}
run_log <- file.path(output_dir, "00_run_log.txt")
writeLines(c("Step 12H: Reviewer 2 Major Comment 4 spatial sensitivity", paste("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), paste("Project root:", PROJECT_ROOT)), run_log)
input_files <- c( step12D_normalized_RCTD_weights = rctd_file, step12F_spot_zones_scores = zone_file, step12F_original_animal_contrasts = original_file )
missing_inputs <- input_files[!file.exists(input_files)]
if (length(missing_inputs)) {
  stop("Required input files missing:\n", paste(missing_inputs, collapse = "\n"))
}
library(data.table)
parameter_audit <- data.table( item = c("project_root", "analysis", "unit_of_inference", "zone_definition", "high_microglia_thresholds", "BAM_covariate", "BAM_model", "expected_total_spots", "expected_injured_spots", "expected_injured_sections", "expected_injured_animals", paste0("file_", names(input_files)), paste0("md5_", names(input_files))), value = c(PROJECT_ROOT, "Step12H R2MC4 spatial source sensitivity", "biological animal", "existing Step12F lesion zones, unchanged", "within-section Microglia RCTD 50th and 75th percentiles", "BAM_macrophage normalized RCTD weight", "target score ~ BAM_macrophage + section fixed effects", "33941", "23791", "24", "6", unname(input_files), unname(as.character(tools::md5sum(input_files)))) )
fwrite(parameter_audit, file.path(tables_dir, "01_input_and_parameter_audit.csv"))
read_gz <- function(path) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)
  as.data.table(read.csv(con, check.names = FALSE, stringsAsFactors = FALSE))
}
z <- read_gz(zone_file)
w <- read_gz(rctd_file)
stopifnot(nrow(z) == 33941L, nrow(w) == 33941L, uniqueN(z$barcode) == 33941L, uniqueN(w$barcode) == 33941L)
needed <- c("barcode", "section_id", "timepoint", "biological_replicate_id", "lesion_zone", "Microglia_RCTD_weight", "target_program_adjusted_z_within_section", "shared_prebranch_adjusted_z_within_section")
stopifnot(all(needed %in% names(z)), all(c("barcode", "BAM_macrophage", "Microglia") %in% names(w)))
x <- merge(z, w[, .(barcode, BAM_macrophage, microglia_check = Microglia)], by = "barcode", all.x = TRUE, sort = FALSE)
stopifnot(nrow(x) == 33941L, !anyNA(x$BAM_macrophage), isTRUE(all.equal(x$Microglia_RCTD_weight, x$microglia_check, tolerance = 1e-10, check.attributes = FALSE)))
x <- x[timepoint %in% c("7days", "2months")]
stopifnot(nrow(x) == 23791L, uniqueN(x$section_id) == 24L, uniqueN(x$biological_replicate_id) == 6L, all(x$lesion_zone %in% c("core", "border", "perilesional", "remote")))
stopifnot(all(is.finite(x$BAM_macrophage)), all(is.finite(x$Microglia_RCTD_weight)), all(is.finite(x$target_program_adjusted_z_within_section)))
fit_pooled <- lm(target_program_adjusted_z_within_section ~ BAM_macrophage + factor(section_id), data = x)
x[, target_BAM_residual_pooled := as.numeric(residuals(fit_pooled))]
model_record <- data.table(model = "injured_pooled_section_fixed_effects", n_spots = nrow(x), bam_coefficient = unname(coef(fit_pooled)["BAM_macrophage"]))
x[, target_BAM_residual_by_timepoint := NA_real_]
for (tp in c("7days", "2months")) {
  ind <- which(x$timepoint == tp)
  df <- x[ind]
  fit <- lm(target_program_adjusted_z_within_section ~ BAM_macrophage + factor(section_id), data = df)
  x[ind, target_BAM_residual_by_timepoint := as.numeric(residuals(fit))]
  model_record <- rbind(model_record, data.table( model = paste0("timepoint_specific_", tp), n_spots = length(ind), bam_coefficient = unname(coef(fit)["BAM_macrophage"])))
}
x[, microglia_q50 := quantile(Microglia_RCTD_weight, 0.50, na.rm = TRUE), by = section_id]
x[, microglia_q75 := quantile(Microglia_RCTD_weight, 0.75, na.rm = TRUE), by = section_id]
scenarios <- list( list(id = "unrestricted_original_target", df = x, metric = "target_program_adjusted_z_within_section"), list(id = "unrestricted_original_shared", df = x, metric = "shared_prebranch_adjusted_z_within_section"), list(id = "microglia_top50_within_section", df = x[Microglia_RCTD_weight >= microglia_q50], metric = "target_program_adjusted_z_within_section"), list(id = "microglia_top25_within_section", df = x[Microglia_RCTD_weight >= microglia_q75], metric = "target_program_adjusted_z_within_section"), list(id = "BAM_adjusted_pooled", df = x, metric = "target_BAM_residual_pooled"), list(id = "BAM_adjusted_by_timepoint", df = x, metric = "target_BAM_residual_by_timepoint"), list(id = "BAM_adjusted_pooled_microglia_top50", df = x[Microglia_RCTD_weight >= microglia_q50], metric = "target_BAM_residual_pooled") )
section_results <- list()
counts <- list()
for (sc in scenarios) {
  dt <- sc$df
  count_table <- dt[, .(n_spots = .N), by = .(section_id, timepoint, biological_replicate_id, lesion_zone)]
  count_table[, analysis := sc$id]
  counts[[length(counts) + 1L]] <- count_table
  means <- dt[, .(mean_score = mean(get(sc$metric))), by = .(section_id, timepoint, biological_replicate_id, lesion_zone)]
  remote <- means[lesion_zone == "remote", .(section_id, remote_mean = mean_score)]
  for (zone in c("core", "border", "perilesional")) {
    sel <- means[lesion_zone == zone]
    v <- merge(sel, remote, by = "section_id", all = FALSE)
    if (nrow(v) != 24L) {
      stop("Missing zone or remote comparison: analysis=", sc$id, ", zone=", zone, ", n_sections=", nrow(v), " (expected 24).")
    }
    v[, `:=`(analysis = sc$id, metric = sc$metric, contrast = paste0(zone, "_minus_remote"), section_contrast = mean_score - remote_mean)]
    section_results[[length(section_results) + 1L]] <- v
  }
}
section_results <- rbindlist(section_results, fill = TRUE)
counts <- rbindlist(counts, fill = TRUE)
animal <- section_results[, .( animal_mean_contrast = mean(section_contrast), n_sections = .N), by = .(analysis, metric, timepoint, biological_replicate_id, contrast)]
stopifnot(all(animal$n_sections == 4L))
summary <- animal[, .( n_animals = .N, mean_animal_effect = mean(animal_mean_contrast), min_animal_effect = min(animal_mean_contrast), max_animal_effect = max(animal_mean_contrast), positive_animals = sum(animal_mean_contrast > 0), negative_animals = sum(animal_mean_contrast < 0)), by = .(analysis, metric, timepoint, contrast)]
summary[, two_sided_sign_p := mapply(function(k, n) {
  binom.test(k, n, p = 0.5)$p.value
}, positive_animals, n_animals)]
orig <- fread(original_file)
ref <- orig[metric == "target_adjusted_within_section_z", .(timepoint, biological_replicate_id, contrast, original_animal_contrast = animal_mean_contrast)]
ours <- animal[analysis == "unrestricted_original_target", .(timepoint, biological_replicate_id, contrast, recalculated_animal_contrast = animal_mean_contrast)]
check <- merge(ours, ref, by = c("timepoint", "biological_replicate_id", "contrast"))
check[, absolute_difference := abs(recalculated_animal_contrast - original_animal_contrast)]
stopifnot(nrow(check) == 18L, max(check$absolute_difference) < 1e-6)
fwrite(counts, file.path(tables_dir, "02_section_zone_spot_counts.csv"))
fwrite(section_results, file.path(tables_dir, "03_section_zone_vs_remote_contrasts.csv"))
fwrite(animal, file.path(tables_dir, "04_animal_zone_vs_remote_contrasts.csv"))
fwrite(summary, file.path(tables_dir, "05_animal_sensitivity_summary.csv"))
fwrite(model_record, file.path(tables_dir, "06_BAM_adjustment_model_audit.csv"))
fwrite(check, file.path(tables_dir, "07_original_step12F_reproduction_audit.csv"))
qc <- data.table( check = c("all_33941_spots_unique_barcode", "injured_23791_spots", "injured_24_sections", "injured_6_animals", "Microglia_weight_matches_Step12D", "all_3_zones_and_remote_in_every_scenario", "each_animal_contains_4_sections", "original_Step12F_18_contrasts_reproduced", "BAM_regressions_have_finite_coefficients", "all_sign_tests_use_3_animals"), passed = c(uniqueN(z$barcode) == 33941L && uniqueN(w$barcode) == 33941L, nrow(x) == 23791L, uniqueN(x$section_id) == 24L, uniqueN(x$biological_replicate_id) == 6L, isTRUE(all.equal(x$Microglia_RCTD_weight, x$microglia_check, tolerance = 1e-10, check.attributes = FALSE)), all(section_results[, .N, by = .(analysis, contrast)]$N == 24L), all(animal$n_sections == 4L), nrow(check) == 18L && max(check$absolute_difference) < 1e-6, all(is.finite(model_record$bam_coefficient)), all(summary$n_animals == 3L)) )
fwrite(qc, file.path(tables_dir, "09_STEP12H_audit_checks.csv"))
if (!all(qc$passed)) {
  stop("Step12H QC failure; inspect tables/09_STEP12H_audit_checks.csv")
}
plot_scenarios <- c( "unrestricted_original_target", "microglia_top50_within_section", "microglia_top25_within_section", "BAM_adjusted_pooled", "BAM_adjusted_by_timepoint" )
plot_names <- c( "Original", "Microglia top 50%", "Microglia top 25%", "BAM-adjusted pooled", "BAM-adjusted by time" )
plot_file <- file.path(figures_dir, "01_exploratory_core_minus_remote_sensitivity.pdf")
pdf(plot_file, width = 12, height = 5.5, useDingbats = FALSE)
tryCatch({
  par(mfrow = c(1, 2), mar = c(9, 4.3, 3, 1.2), oma = c(0, 0, 1, 0))
  all_core <- animal[contrast == "core_minus_remote" & analysis %in% plot_scenarios]
  ylim_all <- range(c(0, all_core$animal_mean_contrast), finite = TRUE)
  ylim_all <- ylim_all + c(-0.08, 0.08) * max(1, diff(ylim_all))
  for (tp in c("7days", "2months")) {
    aa <- all_core[timepoint == tp]
    plot(NA_real_, xlim = c(0.5, length(plot_scenarios) + 0.5), ylim = ylim_all, xaxt = "n", xlab = "", ylab = "Core minus remote score", main = if (tp == "7days") "7 dpi (3 animals)" else "2 months (3 animals)")
    axis(1, at = seq_along(plot_scenarios), labels = plot_names, las = 2, cex.axis = 0.75)
    abline(h = 0, col = "gray50", lty = 2)
    for (i in seq_along(plot_scenarios)) {
      v <- aa[analysis == plot_scenarios[i]]$animal_mean_contrast
      if (length(v) != 3L) stop("Plot expected exactly 3 animals per scenario")
      points(i + c(-0.1, 0, 0.1), sort(v), pch = 1, cex = 1.05, col = "#4C78A8")
      points(i, mean(v), pch = 18, cex = 1.25, col = "#B64B45")
    }
  }
  mtext("Exploratory only; animal-level effects, no cell-of-origin inference", outer = TRUE, line = -0.3, cex = 0.85)
}, finally = dev.off())
core_summary <- summary[contrast == "core_minus_remote", .(analysis, timepoint, mean_animal_effect, positive_animals, n_animals, two_sided_sign_p)]
setorder(core_summary, timepoint, analysis)
capture.output(sessionInfo(), file = file.path(output_dir, "00_sessionInfo.txt"))
print(core_summary)
