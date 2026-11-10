# ============================================================
# 17_step3E_GSE162610_noncycling_transition_refinement_COLORFIX_v2.R
#
# Step 3E purpose:
#   After Step 3D, refine the activated/DAM-like transition by
#   removing the dominant proliferation signal and comparing
#   biologically defined non-cycling state groups.
#
#   This version explicitly sets COLORS for every figure.
#
# Main biological states:
#   Homeostatic_root        = Uninjured + Step3 cluster 3
#   Early_activated_1dpi    = 1 dpi + Step3 clusters 2/5/10
#   Shared_transition_3dpi  = Step3C shared_transition_3dpi_highconf == TRUE
#   Established_DAM_7dpi    = 7 dpi + Step3 cluster 1 + Mature_DAM_like
#
# Main pairwise comparisons:
#   C1 Early_activated_1dpi   vs Homeostatic_root
#   C2 Shared_transition_3dpi vs Early_activated_1dpi
#   C3 Established_DAM_7dpi   vs Shared_transition_3dpi
#
# IMPORTANT:
#   - no repair/scar labels here
#   - no Monocle here
#   - all figures use explicit color/fill scales
# ============================================================

options(stringsAsFactors = FALSE)
set.seed(1234)

# ------------------------------------------------------------
# 0. Packages
# ------------------------------------------------------------
required_pkgs <- c(
  "Seurat", "SeuratObject", "dplyr", "tidyr", "ggplot2",
  "ggrepel", "patchwork", "Matrix", "pheatmap",
  "clusterProfiler", "org.Mm.eg.db", "enrichplot"
)

missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_pkgs) > 0) {
  stop("缺少以下R包，请先安装：\n", paste(missing_pkgs, collapse = ", "))
}

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(Matrix)
  library(pheatmap)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(enrichplot)
})

if (requireNamespace("R.utils", quietly = TRUE)) {
  R.utils::setOption("clusterProfiler.download.method", "auto")
}

# ------------------------------------------------------------
# 1. Paths
# ------------------------------------------------------------
project_dir <- Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
input_rds <- file.path(
  project_dir,
  "04_results",
  "step3D_GSE162610_time_DEG_pathway_confirmation",
  "GSE162610_step3D_time_DEG_pathway_confirmed.rds"
)

out_dir <- file.path(
  project_dir,
  "04_results",
  "step3E_GSE162610_noncycling_transition_refinement"
)

fig_dir  <- file.path(out_dir, "figures")
tab_dir  <- file.path(out_dir, "tables")
go_dir   <- file.path(out_dir, "GO_KEGG")
gsea_dir <- file.path(out_dir, "GSEA")
pb_dir   <- file.path(out_dir, "pseudobulk")
rds_dir  <- file.path(out_dir, "RDS")

for (dd in c(fig_dir, tab_dir, go_dir, gsea_dir, pb_dir, rds_dir)) {
  dir.create(dd, recursive = TRUE, showWarnings = FALSE)
}

if (!file.exists(input_rds)) {
  stop("找不到输入RDS:\n", input_rds)
}

# ------------------------------------------------------------
# 2. Global plotting colors
# ------------------------------------------------------------
# explicit colors for all figures
sig_cols <- c(
  "Sig" = "#D73027",      # red
  "Not Sig" = "#1F9BCF"   # blue
)

state_cols <- c(
  "Homeostatic_root"       = "#66C2A5",
  "Early_activated_1dpi"   = "#FC8D62",
  "Shared_transition_3dpi" = "#8DA0CB",
  "Established_DAM_7dpi"   = "#E78AC3"
)

comparison_tile_cols <- c(
  "Early 1 dpi vs Homeostatic"              = "#8DD3C7",
  "Transition 3 dpi vs Early 1 dpi"         = "#FDB462",
  "Established DAM 7 dpi vs Transition 3 dpi" = "#B3DE69"
)

stage_cols <- c(
  "Early_activation" = "#FC8D62",
  "Transition_3dpi"  = "#8DA0CB",
  "Late_DAM_7dpi"    = "#E78AC3"
)

comparison_point_cols <- c(
  "C1_Early_vs_Homeostatic" = "#D95F02",
  "C2_Transition_vs_Early"  = "#7570B3",
  "C3_DAM_vs_Transition"    = "#1B9E77"
)

hm_colors <- colorRampPalette(
  c("#2166AC", "#F7F7F7", "#B2182B")
)(101)

# ------------------------------------------------------------
# 3. Read object
# ------------------------------------------------------------
message("====================================================")
message("GSE162610 Step 3E | non-cycling transition refinement")
message("====================================================")
message("Reading: ", input_rds)

obj <- readRDS(input_rds)
if (!inherits(obj, "Seurat")) stop("输入对象不是 Seurat object。")

DefaultAssay(obj) <- "RNA"

required_meta <- c(
  "sample_id",
  "timepoint",
  "step3A_cluster_r04",
  "transition_status",
  "shared_transition_3dpi_highconf"
)

missing_meta <- required_meta[!required_meta %in% colnames(obj[[]])]
if (length(missing_meta) > 0) {
  stop("缺少metadata: ", paste(missing_meta, collapse = ", "))
}

obj$step3A_cluster_r04 <- as.character(obj$step3A_cluster_r04)

# ------------------------------------------------------------
# 4. Define Step 3E states
# ------------------------------------------------------------
state_vec <- rep(NA_character_, ncol(obj))
names(state_vec) <- colnames(obj)

state_vec[
  as.character(obj$timepoint) == "Uninjured" &
    obj$step3A_cluster_r04 == "3"
] <- "Homeostatic_root"

state_vec[
  as.character(obj$timepoint) == "1 dpi" &
    obj$step3A_cluster_r04 %in% c("2", "5", "10")
] <- "Early_activated_1dpi"

state_vec[
  as.logical(obj$shared_transition_3dpi_highconf)
] <- "Shared_transition_3dpi"

state_vec[
  as.character(obj$timepoint) == "7 dpi" &
    obj$step3A_cluster_r04 == "1" &
    as.character(obj$transition_status) == "Mature_DAM_like"
] <- "Established_DAM_7dpi"

obj$step3E_state <- factor(
  state_vec[colnames(obj)],
  levels = names(state_cols)
)

state_summary <- obj[[]] %>%
  dplyr::filter(!is.na(step3E_state)) %>%
  dplyr::group_by(step3E_state) %>%
  dplyr::summarise(
    n_cells = dplyr::n(),
    n_samples = dplyr::n_distinct(sample_id),
    samples = paste(sort(unique(as.character(sample_id))), collapse = ";"),
    .groups = "drop"
  )

write.csv(
  state_summary,
  file.path(tab_dir, "01_step3E_state_summary.csv"),
  row.names = FALSE
)

state_audit <- obj[[]] %>%
  dplyr::filter(!is.na(step3E_state)) %>%
  dplyr::count(step3E_state, timepoint, sample_id, step3A_cluster_r04, name = "n_cells")

write.csv(
  state_audit,
  file.path(tab_dir, "02_step3E_state_cell_audit.csv"),
  row.names = FALSE
)

if (any(state_summary$n_cells < 100)) {
  stop("至少一个Step3E state细胞数 <100，请查看 01_step3E_state_summary.csv")
}
if (any(state_summary$n_samples < 2)) {
  stop("至少一个Step3E state来自 <2个biological samples，不宜继续比较。")
}

stage_cells <- colnames(obj)[!is.na(obj$step3E_state)]
obj_stage <- subset(obj, cells = stage_cells)
obj_stage$step3E_state <- factor(as.character(obj_stage$step3E_state), levels = names(state_cols))

saveRDS(
  obj_stage,
  file.path(rds_dir, "GSE162610_step3E_trajectory_candidate_cells.rds")
)

# ------------------------------------------------------------
# 5. State UMAP + sample composition plots (colored)
# ------------------------------------------------------------
if ("umap.microglia.step3A" %in% Reductions(obj_stage)) {
  p_umap_state <- DimPlot(
    obj_stage,
    reduction = "umap.microglia.step3A",
    group.by = "step3E_state",
    cols = state_cols
  ) +
    ggtitle("Step 3E refined biological states") +
    theme_bw()

  p_umap_time <- DimPlot(
    obj_stage,
    reduction = "umap.microglia.step3A",
    group.by = "timepoint"
  ) +
    ggtitle("Step 3E cells by timepoint") +
    theme_bw()

  ggsave(
    file.path(fig_dir, "00_step3E_state_UMAPs.pdf"),
    p_umap_state | p_umap_time,
    width = 13, height = 6
  )
}

state_sample_plot <- obj_stage[[]] %>%
  dplyr::count(step3E_state, sample_id, name = "n_cells") %>%
  dplyr::group_by(step3E_state) %>%
  dplyr::mutate(frac_within_state = n_cells / sum(n_cells)) %>%
  dplyr::ungroup()

p_state_sample <- ggplot(
  state_sample_plot,
  aes(x = sample_id, y = frac_within_state, fill = step3E_state)
) +
  geom_col(position = "dodge") +
  facet_wrap(~step3E_state, scales = "free_x", ncol = 2) +
  scale_fill_manual(values = state_cols, drop = FALSE) +
  theme_bw() +
  labs(
    title = "Biological-sample contribution to each Step 3E state",
    x = "Sample",
    y = "Fraction within state",
    fill = "State"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(
  file.path(fig_dir, "01_step3E_state_sample_composition.pdf"),
  p_state_sample,
  width = 12, height = 8
)

# ------------------------------------------------------------
# 6. DEG helpers
# ------------------------------------------------------------
get_fc_col <- function(df) {
  hit <- intersect(c("avg_log2FC", "avg_logFC"), colnames(df))
  if (length(hit) == 0) stop("FindMarkers输出中无 avg_log2FC/avg_logFC")
  hit[1]
}

run_deg <- function(sobj, case, ref) {
  Idents(sobj) <- "step3E_state"

  x <- FindMarkers(
    sobj,
    ident.1 = case,
    ident.2 = ref,
    assay = "RNA",
    test.use = "wilcox",
    logfc.threshold = 0,
    min.pct = 0.01,
    only.pos = FALSE,
    verbose = FALSE
  )

  x$gene <- rownames(x)
  fc_col <- get_fc_col(x)
  x$log2FC <- as.numeric(x[[fc_col]])
  x$Sig <- ifelse(
    !is.na(x$p_val_adj) & x$p_val_adj < 0.05 & abs(x$log2FC) >= 0.25,
    "Sig", "Not Sig"
  )
  x$DEG_class <- "NS"
  x$DEG_class[x$Sig == "Sig" & x$log2FC >= 0.25] <- "UP"
  x$DEG_class[x$Sig == "Sig" & x$log2FC <= -0.25] <- "DOWN"
  x <- x[order(x$p_val_adj, -abs(x$log2FC)), , drop = FALSE]
  rownames(x) <- NULL
  x
}

# ------------------------------------------------------------
# 7. Pairwise DEGs
# ------------------------------------------------------------
comparisons <- list(
  C1_Early_vs_Homeostatic = c(case = "Early_activated_1dpi", ref = "Homeostatic_root"),
  C2_Transition_vs_Early  = c(case = "Shared_transition_3dpi", ref = "Early_activated_1dpi"),
  C3_DAM_vs_Transition    = c(case = "Established_DAM_7dpi", ref = "Shared_transition_3dpi")
)

comparison_display <- c(
  C1_Early_vs_Homeostatic = "Early 1 dpi vs Homeostatic",
  C2_Transition_vs_Early  = "Transition 3 dpi vs Early 1 dpi",
  C3_DAM_vs_Transition    = "Established DAM 7 dpi vs Transition 3 dpi"
)

deg_list <- list()

for (nm in names(comparisons)) {
  case <- unname(comparisons[[nm]]["case"])
  ref  <- unname(comparisons[[nm]]["ref"])
  message("DEG: ", comparison_display[[nm]])

  dd <- run_deg(obj_stage, case = case, ref = ref)
  dd$comparison <- nm
  dd$comparison_display <- comparison_display[[nm]]
  dd$case <- case
  dd$reference <- ref
  deg_list[[nm]] <- dd

  write.csv(
    dd,
    file.path(tab_dir, paste0("03_DEG_", nm, ".csv")),
    row.names = FALSE
  )
}

deg_all <- dplyr::bind_rows(deg_list)

write.csv(
  deg_all,
  file.path(tab_dir, "04_all_step3E_DEGs_combined.csv"),
  row.names = FALSE
)

deg_summary <- deg_all %>%
  dplyr::group_by(comparison, comparison_display, DEG_class) %>%
  dplyr::summarise(n_genes = dplyr::n(), .groups = "drop")

write.csv(
  deg_summary,
  file.path(tab_dir, "05_step3E_DEG_count_summary.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 8. Multi-volcano with explicit colors (RED/BLUE)
# ------------------------------------------------------------
# User requested colored multi-volcano, not black points.
volcano_df <- deg_all %>%
  dplyr::mutate(
    neg_abs_fc = ifelse(log2FC < 0, -abs(log2FC), abs(log2FC)),
    abs_fc = abs(log2FC)
  )

label_df <- volcano_df %>%
  dplyr::filter(Sig == "Sig") %>%
  dplyr::group_by(comparison_display) %>%
  dplyr::arrange(p_val_adj, dplyr::desc(abs_fc), .by_group = TRUE) %>%
  dplyr::slice_head(n = 12) %>%
  dplyr::ungroup()

# background rectangles per facet
bg_df <- data.frame(
  comparison_display = unname(comparison_display),
  xmin = -Inf,
  xmax = Inf,
  ymin = -Inf,
  ymax = Inf
)

# tile strip inside each panel
tile_df <- data.frame(
  comparison_display = unname(comparison_display),
  label = unname(comparison_display),
  x = 0,
  y = 0,
  stringsAsFactors = FALSE
)

p_volcano <- ggplot(volcano_df, aes(x = log2FC, y = log2FC)) +
  geom_rect(
    data = bg_df,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = "grey92",
    color = NA
  ) +
  geom_hline(yintercept = 0, color = "grey25", linewidth = 0.4) +
  geom_point(
    aes(color = Sig),
    alpha = 0.9,
    size = 1.35
  ) +
  geom_label(
    data = tile_df,
    aes(x = x, y = y, label = label, fill = comparison_display),
    inherit.aes = FALSE,
    color = "black",
    fontface = "bold",
    linewidth = 0.4,
    label.r = unit(0.08, "lines"),
    label.padding = unit(0.18, "lines"),
    size = 4.1
  ) +
  ggrepel::geom_text_repel(
    data = label_df,
    aes(label = gene),
    color = "black",
    size = 3.4,
    max.overlaps = 50,
    segment.color = "black",
    min.segment.length = 0
  ) +
  facet_wrap(~comparison_display, nrow = 1, scales = "free_x") +
  scale_color_manual(values = sig_cols, drop = FALSE) +
  scale_fill_manual(values = comparison_tile_cols, drop = FALSE) +
  theme_bw() +
  labs(
    title = "Step 3E refined-state differential expression",
    subtitle = "Red = significant; Blue = not significant",
    x = "Genes",
    y = "log2FC",
    color = NULL,
    fill = NULL
  ) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "white"),
    strip.text = element_blank(),
    legend.position = "right",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

ggsave(
  file.path(fig_dir, "02_step3E_multiVolcano_colored.pdf"),
  p_volcano,
  width = 18, height = 7
)

# ------------------------------------------------------------
# 9. Data-driven signatures
# ------------------------------------------------------------
c1 <- deg_list[["C1_Early_vs_Homeostatic"]]
c2 <- deg_list[["C2_Transition_vs_Early"]]
c3 <- deg_list[["C3_DAM_vs_Transition"]]

early_sig <- c1 %>%
  dplyr::filter(p_val_adj < 0.05, log2FC >= 0.25) %>%
  dplyr::arrange(dplyr::desc(log2FC))

transition_base <- c2 %>%
  dplyr::filter(p_val_adj < 0.05, log2FC >= 0.25) %>%
  dplyr::select(gene, p_val_adj_transition = p_val_adj, log2FC_transition = log2FC)

late_join <- c3 %>%
  dplyr::select(gene, p_val_adj_late = p_val_adj, log2FC_late = log2FC)

transition_sig <- transition_base %>%
  dplyr::left_join(late_join, by = "gene") %>%
  dplyr::mutate(
    p_val_adj_late = ifelse(is.na(p_val_adj_late), 1, p_val_adj_late),
    log2FC_late = ifelse(is.na(log2FC_late), 0, log2FC_late)
  ) %>%
  dplyr::filter(log2FC_late <= 0.15 | p_val_adj_late >= 0.05) %>%
  dplyr::arrange(dplyr::desc(log2FC_transition))

late_sig <- c3 %>%
  dplyr::filter(p_val_adj < 0.05, log2FC >= 0.25) %>%
  dplyr::arrange(dplyr::desc(log2FC))

write.csv(early_sig, file.path(tab_dir, "06_signature_Early_activation.csv"), row.names = FALSE)
write.csv(transition_sig, file.path(tab_dir, "07_signature_Transition_3dpi.csv"), row.names = FALSE)
write.csv(late_sig, file.path(tab_dir, "08_signature_Late_DAM_7dpi.csv"), row.names = FALSE)

sig_summary <- data.frame(
  stage = c("Early_activation", "Transition_3dpi", "Late_DAM_7dpi"),
  n_genes = c(nrow(early_sig), nrow(transition_sig), nrow(late_sig)),
  definition = c(
    "UP in Early_activated_1dpi vs Homeostatic_root",
    "UP in Shared_transition_3dpi vs Early_activated_1dpi and not further strongly increased in Established_DAM_7dpi",
    "UP in Established_DAM_7dpi vs Shared_transition_3dpi"
  ),
  stringsAsFactors = FALSE
)

write.csv(
  sig_summary,
  file.path(tab_dir, "09_step3E_signature_summary.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 10. Average-expression heatmap with explicit colors
# ------------------------------------------------------------
# Use Seurat-safe labels containing only letters/numbers.
# Seurat may sanitize underscores/hyphens in identity names.
safe_state <- c(
  "Homeostatic_root"       = "State1Homeostatic",
  "Early_activated_1dpi"   = "State2Early",
  "Shared_transition_3dpi" = "State3Transition",
  "Established_DAM_7dpi"   = "State4DAM"
)

safe_levels <- unname(
  safe_state[names(state_cols)]
)

obj_stage$step3E_state_safe <- unname(
  safe_state[
    as.character(obj_stage$step3E_state)
  ]
)

obj_stage$step3E_state_safe <- factor(
  obj_stage$step3E_state_safe,
  levels = safe_levels
)

if (any(is.na(obj_stage$step3E_state_safe))) {
  stop("step3E_state_safe mapping failed.")
}

avg_expr <- AverageExpression(
  obj_stage,
  assays = "RNA",
  group.by = "step3E_state_safe",
  layer = "data",
  verbose = FALSE
)[["RNA"]]

avg_expr <- as.matrix(avg_expr)

observed_avg_cols <- colnames(avg_expr)

if (!all(safe_levels %in% observed_avg_cols)) {
  stop(
    paste0(
      "AverageExpression state columns missing. Expected: ",
      paste(safe_levels, collapse = ", "),
      " | Observed: ",
      paste(observed_avg_cols, collapse = ", ")
    )
  )
}

avg_expr <- avg_expr[
  ,
  safe_levels,
  drop = FALSE
]

# Restore biological state labels after aggregation
colnames(avg_expr) <- names(state_cols)

top_n <- 15
heat_genes <- unique(c(
  head(early_sig$gene, top_n),
  head(transition_sig$gene, top_n),
  head(late_sig$gene, top_n)
))
heat_genes <- intersect(heat_genes, rownames(avg_expr))

row_z <- function(x) {
  sx <- sd(x, na.rm = TRUE)
  if (!is.finite(sx) || sx == 0) return(rep(0, length(x)))
  (x - mean(x, na.rm = TRUE)) / sx
}

heat_z <- t(apply(avg_expr[heat_genes, , drop = FALSE], 1, row_z))

row_stage <- rep("Other", nrow(heat_z))
names(row_stage) <- rownames(heat_z)
row_stage[rownames(heat_z) %in% head(early_sig$gene, top_n)] <- "Early_activation"
row_stage[rownames(heat_z) %in% head(transition_sig$gene, top_n)] <- "Transition_3dpi"
row_stage[rownames(heat_z) %in% head(late_sig$gene, top_n)] <- "Late_DAM_7dpi"

row_anno <- data.frame(
  Signature = row_stage[rownames(heat_z)],
  row.names = rownames(heat_z)
)

col_anno <- data.frame(
  State = colnames(heat_z),
  row.names = colnames(heat_z)
)

pdf(
  file.path(fig_dir, "03_step3E_signature_heatmap_colored.pdf"),
  width = 8.5,
  height = max(9, 0.2 * nrow(heat_z) + 3)
)

pheatmap(
  heat_z,
  color = hm_colors,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_row = row_anno,
  annotation_col = col_anno,
  annotation_colors = list(
    Signature = stage_cols,
    State = state_cols
  ),
  border_color = NA,
  main = "Step 3E data-driven signatures | row Z-score"
)

dev.off()

# ------------------------------------------------------------
# 11. GO/KEGG enrichment
# ------------------------------------------------------------
symbol_to_entrez <- function(symbols) {
  symbols <- unique(symbols[!is.na(symbols) & symbols != ""])
  if (length(symbols) == 0) return(data.frame())
  suppressMessages(
    bitr(
      symbols,
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = org.Mm.eg.db
    )
  )
}

stage_gene_lists <- list(
  Early_activation = early_sig$gene,
  Transition_3dpi  = transition_sig$gene,
  Late_DAM_7dpi    = late_sig$gene
)

go_all <- list()
kegg_all <- list()

for (stage_nm in names(stage_gene_lists)) {
  message("GO/KEGG ORA: ", stage_nm)

  id_map <- symbol_to_entrez(stage_gene_lists[[stage_nm]])
  if (nrow(id_map) < 10) {
    warning(stage_nm, ": mapped genes <10; ORA skipped.")
    next
  }

  entrez <- unique(id_map$ENTREZID)

  ego <- tryCatch(
    enrichGO(
      gene = entrez,
      OrgDb = org.Mm.eg.db,
      keyType = "ENTREZID",
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.20,
      readable = TRUE
    ),
    error = function(e) NULL
  )

  if (!is.null(ego)) {
    ego_df <- as.data.frame(ego)
    if (nrow(ego_df) > 0) {
      ego_df$stage <- stage_nm
      go_all[[stage_nm]] <- ego_df
      write.csv(ego_df, file.path(go_dir, paste0(stage_nm, "_GO_BP.csv")), row.names = FALSE)
    }
  }

  ekegg <- tryCatch(
    enrichKEGG(
      gene = entrez,
      organism = "mmu",
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      qvalueCutoff = 0.20
    ),
    error = function(e) NULL
  )

  if (!is.null(ekegg)) {
    ekegg_df <- as.data.frame(ekegg)
    if (nrow(ekegg_df) > 0) {
      ekegg_df$stage <- stage_nm
      kegg_all[[stage_nm]] <- ekegg_df
      write.csv(ekegg_df, file.path(go_dir, paste0(stage_nm, "_KEGG.csv")), row.names = FALSE)
    }
  }
}

go_combined <- if (length(go_all) > 0) dplyr::bind_rows(go_all) else data.frame()
kegg_combined <- if (length(kegg_all) > 0) dplyr::bind_rows(kegg_all) else data.frame()

write.csv(go_combined, file.path(tab_dir, "10_GO_BP_combined.csv"), row.names = FALSE)
write.csv(kegg_combined, file.path(tab_dir, "11_KEGG_combined.csv"), row.names = FALSE)

make_enrich_plot <- function(df, title, outfile) {
  if (is.null(df) || nrow(df) == 0) return(invisible(NULL))

  dd <- df %>%
    dplyr::filter(!is.na(p.adjust)) %>%
    dplyr::group_by(stage) %>%
    dplyr::arrange(p.adjust, .by_group = TRUE) %>%
    dplyr::slice_head(n = 8) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      minus_log10_padj = -log10(pmax(p.adjust, 1e-300))
    )

  p <- ggplot(
    dd,
    aes(
      x = minus_log10_padj,
      y = reorder(Description, minus_log10_padj),
      color = stage,
      size = Count
    )
  ) +
    geom_point(alpha = 0.9) +
    facet_wrap(~stage, scales = "free_y", ncol = 1) +
    scale_color_manual(values = stage_cols, drop = FALSE) +
    theme_bw() +
    labs(
      title = title,
      x = "-log10 adjusted P",
      y = NULL,
      color = "Stage",
      size = "Gene count"
    )

  ggsave(outfile, p, width = 11, height = max(8, 0.33 * nrow(dd)))
}

make_enrich_plot(
  go_combined,
  "Step 3E GO Biological Process overview",
  file.path(fig_dir, "04_step3E_GO_BP_overview_colored.pdf")
)

make_enrich_plot(
  kegg_combined,
  "Step 3E KEGG overview",
  file.path(fig_dir, "05_step3E_KEGG_overview_colored.pdf")
)

# ------------------------------------------------------------
# 12. GSEA
# ------------------------------------------------------------
make_ranked_entrez <- function(deg_df) {
  df <- deg_df %>%
    dplyr::filter(!is.na(gene), !is.na(log2FC)) %>%
    dplyr::select(gene, log2FC)

  id_map <- symbol_to_entrez(df$gene)
  if (nrow(id_map) == 0) return(numeric(0))

  xx <- df %>%
    dplyr::inner_join(id_map, by = c("gene" = "SYMBOL")) %>%
    dplyr::arrange(dplyr::desc(abs(log2FC))) %>%
    dplyr::distinct(ENTREZID, .keep_all = TRUE)

  rr <- xx$log2FC
  names(rr) <- xx$ENTREZID
  sort(rr, decreasing = TRUE)
}

gsea_kegg_all <- list()
gsea_go_all <- list()

for (nm in names(deg_list)) {
  message("GSEA: ", comparison_display[[nm]])

  rank_list <- make_ranked_entrez(deg_list[[nm]])
  if (length(rank_list) < 100) {
    warning(nm, ": rank list <100; GSEA skipped.")
    next
  }

  ggo <- tryCatch(
    gseGO(
      geneList = rank_list,
      OrgDb = org.Mm.eg.db,
      keyType = "ENTREZID",
      ont = "BP",
      minGSSize = 10,
      maxGSSize = 500,
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      verbose = FALSE
    ),
    error = function(e) NULL
  )

  if (!is.null(ggo)) {
    ggo_df <- as.data.frame(ggo)
    if (nrow(ggo_df) > 0) {
      ggo_df$comparison <- nm
      gsea_go_all[[nm]] <- ggo_df
      write.csv(ggo_df, file.path(gsea_dir, paste0(nm, "_GSEA_GO_BP.csv")), row.names = FALSE)
    }
  }

  gkegg <- tryCatch(
    gseKEGG(
      geneList = rank_list,
      organism = "mmu",
      minGSSize = 10,
      maxGSSize = 500,
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      verbose = FALSE
    ),
    error = function(e) NULL
  )

  if (!is.null(gkegg)) {
    gkegg_df <- as.data.frame(gkegg)
    if (nrow(gkegg_df) > 0) {
      gkegg_df$comparison <- nm
      gsea_kegg_all[[nm]] <- gkegg_df
      write.csv(gkegg_df, file.path(gsea_dir, paste0(nm, "_GSEA_KEGG.csv")), row.names = FALSE)
    }
  }
}

gsea_go_combined <- if (length(gsea_go_all) > 0) dplyr::bind_rows(gsea_go_all) else data.frame()
gsea_kegg_combined <- if (length(gsea_kegg_all) > 0) dplyr::bind_rows(gsea_kegg_all) else data.frame()

write.csv(gsea_go_combined, file.path(tab_dir, "12_GSEA_GO_BP_combined.csv"), row.names = FALSE)
write.csv(gsea_kegg_combined, file.path(tab_dir, "13_GSEA_KEGG_combined.csv"), row.names = FALSE)

make_gsea_plot <- function(df, title, outfile) {
  if (is.null(df) || nrow(df) == 0) return(invisible(NULL))

  dd <- df %>%
    dplyr::filter(!is.na(NES), !is.na(p.adjust)) %>%
    dplyr::group_by(comparison) %>%
    dplyr::arrange(p.adjust, dplyr::desc(abs(NES)), .by_group = TRUE) %>%
    dplyr::slice_head(n = 8) %>%
    dplyr::ungroup()

  dd$comparison_display <- comparison_display[dd$comparison]

  p <- ggplot(
    dd,
    aes(
      x = NES,
      y = reorder(Description, NES),
      color = comparison,
      size = -log10(pmax(p.adjust, 1e-300))
    )
  ) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    geom_point(alpha = 0.9) +
    facet_wrap(~comparison_display, scales = "free_y", ncol = 1) +
    scale_color_manual(values = comparison_point_cols, drop = FALSE) +
    theme_bw() +
    labs(
      title = title,
      x = "NES",
      y = NULL,
      color = "Comparison",
      size = "-log10 adj.P"
    )

  ggsave(outfile, p, width = 11, height = max(8, 0.35 * nrow(dd)))
}

make_gsea_plot(
  gsea_kegg_combined,
  "Step 3E GSEA KEGG overview",
  file.path(fig_dir, "06_step3E_GSEA_KEGG_overview_colored.pdf")
)

make_gsea_plot(
  gsea_go_combined,
  "Step 3E GSEA GO-BP overview",
  file.path(fig_dir, "07_step3E_GSEA_GO_BP_overview_colored.pdf")
)

# ------------------------------------------------------------
# 13. Sample-level pseudobulk edgeR
# ------------------------------------------------------------
pseudobulk_run <- FALSE
pb_deg_list <- list()

if (requireNamespace("edgeR", quietly = TRUE)) {
  message("edgeR detected: running pseudobulk validation...")

  get_counts <- function(sobj) {
    x <- tryCatch(GetAssayData(sobj, assay = "RNA", layer = "counts"), error = function(e) NULL)
    if (is.null(x)) x <- tryCatch(GetAssayData(sobj, assay = "RNA", slot = "counts"), error = function(e) NULL)
    if (is.null(x)) stop("无法取得RNA counts。")
    x
  }

  counts <- get_counts(obj_stage)

  if (!inherits(counts, "dgCMatrix")) {
    counts <- as(counts, "dgCMatrix")
  }

  # Ensure exact cell order between expression matrix and metadata.
  if (!identical(colnames(counts), colnames(obj_stage))) {
    counts <- counts[
      ,
      colnames(obj_stage),
      drop = FALSE
    ]
  }

  meta_all <- obj_stage[[]]
  meta_all <- meta_all[
    colnames(obj_stage),
    ,
    drop = FALSE
  ]

  sample_meta <- data.frame(
    cell = colnames(obj_stage),
    sample_id = as.character(meta_all$sample_id),
    step3E_state = as.character(meta_all$step3E_state),
    stringsAsFactors = FALSE
  )

  if (any(is.na(sample_meta$step3E_state))) {
    stop("Unexpected NA step3E_state inside obj_stage.")
  }

  sample_meta$sample_state <- paste(
    sample_meta$step3E_state,
    sample_meta$sample_id,
    sep = "__"
  )

  pb_groups <- factor(
    sample_meta$sample_state,
    levels = unique(sample_meta$sample_state)
  )

  # cells x sample-state design matrix
  mm <- Matrix::sparse.model.matrix(
    ~ 0 + pb_groups
  )

  if (nrow(mm) != ncol(counts)) {
    stop(
      paste0(
        "Pseudobulk design mismatch: design rows = ",
        nrow(mm),
        ", count cells = ",
        ncol(counts)
      )
    )
  }

  # genes x cells  %*%  cells x sample-state
  # = genes x sample-state pseudobulk count matrix
  pb_counts <- counts %*% mm

  pb_colnames <- colnames(mm)
  pb_colnames <- sub(
    "^pb_groups",
    "",
    pb_colnames
  )

  colnames(pb_counts) <- pb_colnames

  pb_sample_state <- data.frame(
    sample_state = colnames(pb_counts),
    stringsAsFactors = FALSE
  )
  tmp <- strsplit(
    pb_sample_state$sample_state,
    "__",
    fixed = TRUE
  )

  if (any(lengths(tmp) != 2)) {
    stop(
      "Unexpected sample_state parsing failure: ",
      paste(
        pb_sample_state$sample_state[lengths(tmp) != 2],
        collapse = ", "
      )
    )
  }

  pb_sample_state$state <- vapply(
    tmp,
    function(z) z[[1]],
    character(1)
  )

  pb_sample_state$sample_id <- vapply(
    tmp,
    function(z) z[[2]],
    character(1)
  )

  write.csv(pb_sample_state, file.path(pb_dir, "00_pseudobulk_sample_state_metadata.csv"), row.names = FALSE)

  run_pb <- function(case, ref, nm) {
    keep <- pb_sample_state$state %in% c(ref, case)
    cc <- pb_counts[, keep, drop = FALSE]
    meta <- pb_sample_state[keep, , drop = FALSE]
    group <- factor(meta$state, levels = c(ref, case))

    y <- edgeR::DGEList(counts = cc, group = group)
    keep_gene <- edgeR::filterByExpr(y, group = group)
    y <- y[keep_gene, , keep.lib.sizes = FALSE]
    y <- edgeR::calcNormFactors(y)
    design <- model.matrix(~group)
    y <- edgeR::estimateDisp(y, design, robust = TRUE)
    fit <- edgeR::glmQLFit(y, design, robust = TRUE)
    qlf <- edgeR::glmQLFTest(fit, coef = 2)
    tt <- edgeR::topTags(qlf, n = Inf, sort.by = "PValue")$table
    tt$gene <- rownames(tt)
    tt$comparison <- nm
    rownames(tt) <- NULL
    tt
  }

  for (nm in names(comparisons)) {
    case <- unname(comparisons[[nm]]["case"])
    ref  <- unname(comparisons[[nm]]["ref"])
    pb <- tryCatch(
      run_pb(case, ref, nm),
      error = function(e) {
        message(
          "Pseudobulk edgeR failed for ",
          nm,
          ": ",
          conditionMessage(e)
        )
        NULL
      }
    )
    if (!is.null(pb)) {
      pb_deg_list[[nm]] <- pb
      write.csv(pb, file.path(pb_dir, paste0("01_pseudobulk_", nm, ".csv")), row.names = FALSE)
    }
  }

  if (length(pb_deg_list) == length(comparisons)) pseudobulk_run <- TRUE
} else {
  message("edgeR not installed: pseudobulk skipped safely.")
}

# ------------------------------------------------------------
# 14. Cell vs pseudobulk concordance figure (colored)
# ------------------------------------------------------------
pb_concordance <- data.frame()

if (length(pb_deg_list) > 0) {
  cc_list <- list()

  for (nm in names(pb_deg_list)) {
    cell_df <- deg_list[[nm]] %>%
      dplyr::select(gene, cell_log2FC = log2FC, cell_padj = p_val_adj, Sig)

    pb_df <- pb_deg_list[[nm]] %>%
      dplyr::select(gene, pb_log2FC = logFC, pb_FDR = FDR)

    xx <- dplyr::inner_join(cell_df, pb_df, by = "gene")
    xx$comparison <- nm
    xx$direction_concordant <- sign(xx$cell_log2FC) == sign(xx$pb_log2FC)
    cc_list[[nm]] <- xx
  }

  pb_concordance <- dplyr::bind_rows(cc_list)

  write.csv(
    pb_concordance,
    file.path(tab_dir, "14_cell_vs_pseudobulk_concordance.csv"),
    row.names = FALSE
  )

  p_pb <- ggplot(
    pb_concordance,
    aes(
      x = cell_log2FC,
      y = pb_log2FC,
      color = comparison
    )
  ) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    geom_point(alpha = 0.35, size = 0.9) +
    facet_wrap(~comparison, scales = "free", nrow = 1) +
    scale_color_manual(values = comparison_point_cols, drop = FALSE) +
    theme_bw() +
    labs(
      title = "Step 3E cell-level vs sample-level pseudobulk concordance",
      x = "Cell-level log2FC",
      y = "Pseudobulk log2FC",
      color = "Comparison"
    )

  ggsave(
    file.path(fig_dir, "08_step3E_cell_vs_pseudobulk_concordance_colored.pdf"),
    p_pb,
    width = 15, height = 5
  )
}

# ------------------------------------------------------------
# 15. Save final object
# ------------------------------------------------------------
saveRDS(
  obj_stage,
  file.path(out_dir, "GSE162610_step3E_noncycling_transition_refined.rds")
)

# ------------------------------------------------------------
# 16. Report
# ------------------------------------------------------------
report_file <- file.path(out_dir, "17_STEP3E_noncycling_transition_refinement_report.txt")
sink(report_file)

cat("GSE162610 Step 3E | non-cycling transition refinement\n")
cat("====================================================\n\n")
cat("Input RDS:\n", input_rds, "\n\n", sep = "")

cat("Selected states:\n")
print(state_summary)

cat("\nState cell audit (top rows):\n")
print(utils::head(state_audit, 20))

cat("\nPairwise comparisons:\n")
cat("C1 Early_activated_1dpi vs Homeostatic_root\n")
cat("C2 Shared_transition_3dpi vs Early_activated_1dpi\n")
cat("C3 Established_DAM_7dpi vs Shared_transition_3dpi\n")

cat("\nDEG count summary:\n")
print(deg_summary)

cat("\nSignature summary:\n")
print(sig_summary)

cat("\nGO rows: ", nrow(go_combined), "\n", sep = "")
cat("KEGG rows: ", nrow(kegg_combined), "\n", sep = "")
cat("GSEA GO rows: ", nrow(gsea_go_combined), "\n", sep = "")
cat("GSEA KEGG rows: ", nrow(gsea_kegg_combined), "\n", sep = "")
cat("Pseudobulk run: ", pseudobulk_run, "\n", sep = "")

if (nrow(pb_concordance) > 0) {
  pb_summary <- pb_concordance %>%
    dplyr::group_by(comparison) %>%
    dplyr::summarise(
      n_genes = dplyr::n(),
      direction_concordance = mean(direction_concordant, na.rm = TRUE),
      .groups = "drop"
    )
  cat("\nCell vs pseudobulk direction concordance:\n")
  print(pb_summary)
}

cat("\nIMPORTANT\n")
cat("1. This step refines the shared activated/DAM-like transition after removing the dominant whole-timepoint cell-cycle effect.\n")
cat("2. No repair/scar labels are assigned here.\n")
cat("3. No Monocle is run here.\n")
cat("4. All figures in this script use explicit colors (not black defaults).\n")
cat("5. Review C2 and C3 carefully before Step 4.\n")

sink()

# ------------------------------------------------------------
# 17. Console output
# ------------------------------------------------------------
message("====================================================")
message("STEP 3E COMPLETE")
message("Results: ", out_dir)
message("")
message("ALL FIGURES IN THIS VERSION USE EXPLICIT COLORS.")
message("Multi-volcano red/blue fix applied.")
message("")
message("请优先上传以下 7 个文件：")
message("1) figures/02_step3E_multiVolcano_colored.pdf")
message("2) figures/03_step3E_signature_heatmap_colored.pdf")
message("3) figures/04_step3E_GO_BP_overview_colored.pdf （若生成）")
message("4) figures/06_step3E_GSEA_KEGG_overview_colored.pdf （若生成）")
message("5) tables/09_step3E_signature_summary.csv")
message("6) figures/08_step3E_cell_vs_pseudobulk_concordance_colored.pdf （若edgeR已运行）")
message("7) 17_STEP3E_noncycling_transition_refinement_report.txt")
message("====================================================")
