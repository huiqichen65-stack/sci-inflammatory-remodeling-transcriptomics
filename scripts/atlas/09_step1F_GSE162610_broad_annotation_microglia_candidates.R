# ============================================================
# Step 1F｜GSE162610：19 clusters marker 鉴定 + broad lineage 初筛
# 课题：成人 SCI 后 activated/DAM-like microglia 的 repair/scar 状态分叉
#
# 基于你上传的 10annotation 模板改写。
#
# 本步骤目标：
#   1) 固定使用 Step 1E 的 broad atlas resolution = 0.4（19 clusters）
#   2) 找每个 cluster 的 markers
#   3) 用小鼠脊髓 major-lineage marker panel 做 DotPlot/FeaturePlot
#   4) 计算 cluster-level lineage program score
#   5) 生成“broad_type_suggestion”，但暂不把它当最终人工注释
#   6) 重点区分：
#        resident microglia
#        vs infiltrating monocyte/macrophage
#        vs neutrophil
#   7) 自动给出 strict / extended candidate microglia clusters
#   8) 如 Step 0 原始 metadata 中存在作者 cell-type annotation，
#      自动识别并作为辅助证据，不直接照搬
#
# 重要：
#   - 本步骤不直接把 broad_type_suggestion 当最终 celltype。
#   - 先把结果发给我判断后，再决定哪些 cluster 是真正 microglia。
#   - 后续会提取 clean microglia，再单独重新 HVG/PCA/Harmony/聚类。
#
# 输入：
#   Desktop/xjz/02_objects/step1E_GSE162610_Harmony_cluster/
#     GSE162610_Harmony_multires_clustered.rds
#
# 输出：
#   Desktop/xjz/02_objects/step1F_GSE162610_broad_annotation/
#
# 运行方式：
#   RStudio 打开本脚本 -> Source
# ============================================================

options(stringsAsFactors = FALSE)
options(timeout = 1800)

set.seed(20260730)

# ------------------------------------------------------------
# 0. 路径
# ------------------------------------------------------------
project_dir <- Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
in_rds <- file.path(
  project_dir,
  "02_objects",
  "step1E_GSE162610_Harmony_cluster",
  "GSE162610_Harmony_multires_clustered.rds"
)

step0_meta_file <- file.path(
  project_dir,
  "00_admin",
  "step0_GSE162610_confirmed",
  "01_GSE162610_cell_metadata_CONFIRMED.csv"
)

out_dir <- file.path(
  project_dir,
  "02_objects",
  "step1F_GSE162610_broad_annotation"
)

fig_dir <- file.path(out_dir, "figures")
tab_dir <- file.path(out_dir, "tables")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(out_dir, "00_run_log.txt")

log_msg <- function(...) {
  msg <- paste0(...)
  cat(msg, "\n")
  cat(
    paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", msg, "\n"),
    file = log_file,
    append = TRUE
  )
}

writeLines(
  paste0("GSE162610 Step 1F started: ", Sys.time()),
  log_file
)

if (!file.exists(in_rds)) {
  stop("找不到 Step 1E 对象：\n", in_rds)
}

# ------------------------------------------------------------
# 1. R 包
# ------------------------------------------------------------
pkgs <- c(
  "Seurat",
  "data.table",
  "dplyr",
  "ggplot2",
  "patchwork"
)

missing_pkgs <- pkgs[
  !vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_pkgs) > 0) {
  stop(
    "缺少 R 包：",
    paste(missing_pkgs, collapse = ", "),
    "\n请先安装后重新 Source。"
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

log_msg("============================================================")
log_msg("Step 1F | Cluster markers + broad lineage annotation audit")
log_msg("============================================================")

# ------------------------------------------------------------
# 2. 读取对象并固定 broad clustering
# ------------------------------------------------------------
pbmc <- readRDS(in_rds)

DefaultAssay(pbmc) <- "RNA"

required_meta <- c(
  "sample_id",
  "timepoint",
  "cluster_provisional"
)

if (!all(required_meta %in% colnames(pbmc@meta.data))) {
  stop(
    "对象缺少 metadata：",
    paste(
      setdiff(required_meta, colnames(pbmc@meta.data)),
      collapse = ", "
    )
  )
}

if (!"umap.harmony" %in% Reductions(pbmc)) {
  stop("对象中没有 umap.harmony reduction。")
}

pbmc$cluster_provisional <- as.character(pbmc$cluster_provisional)
Idents(pbmc) <- "cluster_provisional"

cluster_ids <- sort(
  unique(as.character(pbmc$cluster_provisional))
)

log_msg(
  "Input object: ",
  nrow(pbmc), " genes x ",
  ncol(pbmc), " cells"
)

log_msg(
  "Broad clusters: ",
  length(cluster_ids),
  " (resolution 0.4)"
)

if (length(cluster_ids) != 19) {
  warning(
    "当前 cluster 数不是预期的 19，而是 ",
    length(cluster_ids)
  )
}

# ------------------------------------------------------------
# 3. Marker gene sets（mouse spinal cord）
# ------------------------------------------------------------
marker_sets_raw <- list(

  Microglia = c(
    "P2ry12", "Tmem119", "Fcrls", "Hexb",
    "Cx3cr1", "Sall1", "Siglech", "Csf1r"
  ),

  Myeloid_general = c(
    "Aif1", "Tyrobp", "Lyz2", "Csf1r",
    "Fcerg1", "Ctss", "Itgam", "Spi1"
  ),

  Macrophage_Monocyte = c(
    "Ccr2", "Ly6c2", "Ms4a7", "Lgals3",
    "Fcgr1", "Fn1", "S100a8", "S100a9"
  ),

  Neutrophil = c(
    "S100a8", "S100a9", "Retnlg", "Csf3r",
    "Ly6g", "Mpo", "Ngp", "Camp"
  ),

  Astrocyte = c(
    "Aqp4", "Gfap", "Aldoc", "Slc1a3",
    "S100b", "Sox9", "Fabp7", "Glul"
  ),

  Oligodendrocyte = c(
    "Mbp", "Plp1", "Mog", "Mag",
    "Cnp", "Mobp", "Mal", "Opalin"
  ),

  OPC = c(
    "Pdgfra", "Cspg4", "Olig1", "Olig2",
    "Sox10", "Ptprz1", "Vcan", "Nkx2-2"
  ),

  Endothelial = c(
    "Pecam1", "Cdh5", "Kdr", "Vwf",
    "Emcn", "Esam", "Klf2", "Rbp7"
  ),

  Fibroblast = c(
    "Col1a1", "Col1a2", "Col3a1", "Dcn",
    "Lum", "Pdgfra", "Col5a1", "Fn1"
  ),

  Pericyte_VSMC = c(
    "Pdgfrb", "Rgs5", "Mcam", "Cspg4",
    "Des", "Acta2", "Tagln", "Cnn1"
  ),

  Neuron = c(
    "Snap25", "Syt1", "Rbfox3", "Tubb3",
    "Map2", "Slc17a7", "Gad1", "Gad2"
  ),

  T_cell = c(
    "Cd3d", "Cd3e", "Cd3g", "Cd247",
    "Trbc1", "Trbc2", "Lck", "Il7r"
  ),

  B_cell = c(
    "Cd79a", "Cd79b", "Ms4a1", "Cd37",
    "Cd74", "H2-Aa", "Cd22", "Cd19"
  ),

  Mast = c(
    "Kit", "Ms4a2", "Cpa3", "Tpsb2",
    "Fcer1a", "Hdc", "Mcpt4", "Gata2"
  ),

  Ependymal = c(
    "Foxj1", "Pifo", "Tppp3", "Rsph1",
    "Cfap44", "Dynlrb2", "Fam183b", "Ccdc153"
  )
)

# cell-cycle 仅用于识别 cycling cluster，不做回归
cycling_markers_raw <- c(
  "Mki67", "Top2a", "Ube2c", "Birc5",
  "Cenpf", "Ccnb1", "Ccnb2", "Tpx2"
)

# ------------------------------------------------------------
# 4. Case-insensitive 匹配基因名
# ------------------------------------------------------------
match_genes <- function(x, universe) {

  uni_lower <- tolower(universe)
  idx <- match(tolower(x), uni_lower)

  unique(
    universe[
      idx[!is.na(idx)]
    ]
  )
}

marker_sets <- lapply(
  marker_sets_raw,
  match_genes,
  universe = rownames(pbmc)
)

cycling_markers <- match_genes(
  cycling_markers_raw,
  rownames(pbmc)
)

marker_availability <- do.call(
  rbind,
  lapply(
    names(marker_sets_raw),
    function(nm) {

      original <- marker_sets_raw[[nm]]
      matched <- marker_sets[[nm]]

      data.frame(
        program = nm,
        n_requested = length(original),
        n_found = length(matched),
        found_genes = paste(matched, collapse = ";"),
        missing_genes = paste(
          setdiff(
            original,
            original[
              tolower(original) %in%
                tolower(matched)
            ]
          ),
          collapse = ";"
        ),
        stringsAsFactors = FALSE
      )
    }
  )
)

fwrite(
  marker_availability,
  file.path(
    tab_dir,
    "00_marker_gene_availability.csv"
  )
)

# 去掉找到少于 2 个基因的 program
marker_sets <- marker_sets[
  vapply(marker_sets, length, integer(1)) >= 2
]

log_msg(
  "Usable lineage programs: ",
  paste(names(marker_sets), collapse = ", ")
)

# ------------------------------------------------------------
# 5. FindAllMarkers
# ------------------------------------------------------------
# 与上传模板相比：
#   logFC 阈值从 1 调低到 0.25，
#   这样不会漏掉 broad annotation 所需的中等强度 marker。
#
# only.pos = TRUE
# min.pct = 0.1
# return.thresh = 0.01
# ------------------------------------------------------------
log_msg("Running FindAllMarkers for 19 broad clusters ...")
log_msg("This may take several minutes. Please do not press Stop.")

t0 <- Sys.time()

pbmc_markers <- FindAllMarkers(
  object = pbmc,
  assay = "RNA",
  only.pos = TRUE,
  min.pct = 0.10,
  logfc.threshold = 0.25,
  test.use = "wilcox",
  return.thresh = 0.01,
  random.seed = 20260730,
  verbose = TRUE
)

marker_runtime <- as.numeric(
  difftime(
    Sys.time(),
    t0,
    units = "mins"
  )
)

log_msg(
  "FindAllMarkers finished in ",
  round(marker_runtime, 2),
  " min"
)

if (nrow(pbmc_markers) == 0) {
  stop("FindAllMarkers 没有返回任何 marker。")
}

fwrite(
  pbmc_markers,
  file.path(
    tab_dir,
    "01_all_cluster_markers.csv"
  )
)

# ------------------------------------------------------------
# 6. Top20 markers / cluster
# ------------------------------------------------------------
if (!"avg_log2FC" %in% colnames(pbmc_markers)) {
  stop("marker 表中没有 avg_log2FC 列。")
}

top20_markers <- pbmc_markers %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  slice_head(n = 20) %>%
  ungroup()

fwrite(
  top20_markers,
  file.path(
    tab_dir,
    "02_top20_markers_per_cluster.csv"
  )
)

top50_markers <- pbmc_markers %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  slice_head(n = 50) %>%
  ungroup()

# ------------------------------------------------------------
# 7. Broad lineage DotPlot
# ------------------------------------------------------------
broad_panel <- c(
  marker_sets$Microglia,
  marker_sets$Macrophage_Monocyte,
  marker_sets$Astrocyte,
  marker_sets$Oligodendrocyte,
  marker_sets$OPC,
  marker_sets$Endothelial,
  marker_sets$Fibroblast,
  marker_sets$Pericyte_VSMC,
  marker_sets$Neuron,
  marker_sets$T_cell,
  marker_sets$B_cell,
  marker_sets$Mast,
  marker_sets$Ependymal
)

broad_panel <- unique(broad_panel)

p_broad_dot <- DotPlot(
  pbmc,
  features = broad_panel,
  group.by = "cluster_provisional",
  assay = "RNA"
) +
  RotatedAxis() +
  ggtitle(
    "Broad lineage markers across 19 clusters"
  ) +
  theme(
    axis.text.x = element_text(
      angle = 60,
      hjust = 1,
      vjust = 1
    )
  )

ggsave(
  file.path(
    fig_dir,
    "01_broad_lineage_marker_DotPlot.pdf"
  ),
  p_broad_dot,
  width = 24,
  height = 9
)

# ------------------------------------------------------------
# 8. Myeloid discrimination DotPlot
# ------------------------------------------------------------
myeloid_panel <- unique(c(
  marker_sets$Microglia,
  marker_sets$Myeloid_general,
  marker_sets$Macrophage_Monocyte,
  marker_sets$Neutrophil
))

p_myeloid_dot <- DotPlot(
  pbmc,
  features = myeloid_panel,
  group.by = "cluster_provisional",
  assay = "RNA"
) +
  RotatedAxis() +
  ggtitle(
    "Microglia vs infiltrating myeloid discrimination"
  ) +
  theme(
    axis.text.x = element_text(
      angle = 60,
      hjust = 1,
      vjust = 1
    )
  )

ggsave(
  file.path(
    fig_dir,
    "02_myeloid_discrimination_DotPlot.pdf"
  ),
  p_myeloid_dot,
  width = 18,
  height = 8
)

# ------------------------------------------------------------
# 9. Cycling marker DotPlot
# ------------------------------------------------------------
if (length(cycling_markers) >= 2) {

  p_cycle <- DotPlot(
    pbmc,
    features = cycling_markers,
    group.by = "cluster_provisional",
    assay = "RNA"
  ) +
    RotatedAxis() +
    ggtitle(
      "Cycling markers (diagnostic only; not regressed)"
    )

  ggsave(
    file.path(
      fig_dir,
      "03_cycling_marker_DotPlot.pdf"
    ),
    p_cycle,
    width = 10,
    height = 7
  )
}

# ------------------------------------------------------------
# 10. Representative FeaturePlots
# ------------------------------------------------------------
representative_raw <- c(
  "P2ry12",
  "Tmem119",
  "Ccr2",
  "S100a8",
  "Aqp4",
  "Mbp",
  "Pdgfra",
  "Pecam1",
  "Col1a1",
  "Pdgfrb",
  "Snap25",
  "Cd3d"
)

representative <- match_genes(
  representative_raw,
  rownames(pbmc)
)

if (length(representative) > 0) {

  p_feat <- FeaturePlot(
    pbmc,
    features = representative,
    reduction = "umap.harmony",
    ncol = 4,
    pt.size = 0.08,
    order = TRUE
  )

  ggsave(
    file.path(
      fig_dir,
      "04_representative_marker_FeaturePlots.pdf"
    ),
    p_feat,
    width = 16,
    height = 12
  )
}

# ------------------------------------------------------------
# 11. Cluster UMAP
# ------------------------------------------------------------
p_cluster <- DimPlot(
  pbmc,
  reduction = "umap.harmony",
  group.by = "cluster_provisional",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.08
) +
  ggtitle(
    "GSE162610 broad clusters (resolution 0.4)"
  )

ggsave(
  file.path(
    fig_dir,
    "05_cluster_UMAP_res0.4.pdf"
  ),
  p_cluster,
  width = 10,
  height = 8
)

# ------------------------------------------------------------
# 12. Cluster-level lineage program scores
# ------------------------------------------------------------
all_score_genes <- unique(
  unlist(marker_sets)
)

# 平均 log-normalized expression
avg_expr <- AverageExpression(
  pbmc,
  assays = "RNA",
  features = all_score_genes,
  group.by = "cluster_provisional",
  slot = "data",
  verbose = FALSE
)$RNA

# 保证列顺序可追溯
avg_expr <- as.matrix(avg_expr)

# 每个 gene 在 cluster 之间做 z-score
z_by_gene <- t(
  apply(
    avg_expr,
    1,
    function(x) {
      s <- sd(x)
      if (is.na(s) || s == 0) {
        rep(0, length(x))
      } else {
        (x - mean(x)) / s
      }
    }
  )
)

colnames(z_by_gene) <- colnames(avg_expr)
rownames(z_by_gene) <- rownames(avg_expr)

program_score_mat <- sapply(
  names(marker_sets),
  function(program) {

    gg <- intersect(
      marker_sets[[program]],
      rownames(z_by_gene)
    )

    if (length(gg) == 1) {
      as.numeric(z_by_gene[gg, ])
    } else {
      colMeans(
        z_by_gene[gg, , drop = FALSE],
        na.rm = TRUE
      )
    }
  }
)

program_score_mat <- as.matrix(program_score_mat)

rownames(program_score_mat) <- colnames(avg_expr)
colnames(program_score_mat) <- names(marker_sets)

program_score_df <- data.frame(
  cluster = rownames(program_score_mat),
  program_score_mat,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

fwrite(
  program_score_df,
  file.path(
    tab_dir,
    "03_cluster_lineage_program_scores.csv"
  )
)

# ------------------------------------------------------------
# 13. broad_type_suggestion
# ------------------------------------------------------------
score_only <- program_score_df[
  ,
  setdiff(
    colnames(program_score_df),
    "cluster"
  ),
  drop = FALSE
]

broad_type_suggestion <- colnames(score_only)[
  max.col(
    as.matrix(score_only),
    ties.method = "first"
  )
]

broad_suggestion_df <- data.frame(
  cluster = program_score_df$cluster,
  broad_type_suggestion = broad_type_suggestion,
  top_program_score = apply(
    as.matrix(score_only),
    1,
    max
  ),
  second_program_score = apply(
    as.matrix(score_only),
    1,
    function(x) {
      sort(x, decreasing = TRUE)[2]
    }
  ),
  stringsAsFactors = FALSE
)

broad_suggestion_df$score_margin <- (
  broad_suggestion_df$top_program_score -
    broad_suggestion_df$second_program_score
)

# 标记低置信度
broad_suggestion_df$suggestion_confidence <- ifelse(
  broad_suggestion_df$score_margin >= 0.5,
  "Higher",
  "Needs_manual_review"
)

fwrite(
  broad_suggestion_df,
  file.path(
    tab_dir,
    "04_cluster_broad_type_suggestions.csv"
  )
)

# 映射回每个 cell
cluster_to_type <- setNames(
  broad_suggestion_df$broad_type_suggestion,
  broad_suggestion_df$cluster
)

pbmc$broad_type_suggestion <- unname(
  cluster_to_type[
    as.character(pbmc$cluster_provisional)
  ]
)

# ------------------------------------------------------------
# 14. broad type suggestion UMAP
# ------------------------------------------------------------
p_broad_umap <- DimPlot(
  pbmc,
  reduction = "umap.harmony",
  group.by = "broad_type_suggestion",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.08
) +
  ggtitle(
    "Broad lineage suggestion (algorithm-assisted; not final)"
  )

ggsave(
  file.path(
    fig_dir,
    "06_broad_type_suggestion_UMAP.pdf"
  ),
  p_broad_umap,
  width = 11,
  height = 8
)

# ------------------------------------------------------------
# 15. Lineage score heatmap
# ------------------------------------------------------------
score_long <- data.table::melt(
  data.table::as.data.table(program_score_df),
  id.vars = "cluster",
  variable.name = "program",
  value.name = "score"
)

p_score_heat <- ggplot(
  score_long,
  aes(
    x = program,
    y = factor(
      cluster,
      levels = rev(cluster_ids)
    ),
    fill = score
  )
) +
  geom_tile() +
  labs(
    x = NULL,
    y = "Cluster",
    fill = "Program\nz-score",
    title = "Cluster lineage-program score heatmap"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 55,
      hjust = 1
    )
  )

ggsave(
  file.path(
    fig_dir,
    "07_cluster_lineage_score_heatmap.pdf"
  ),
  p_score_heat,
  width = 12,
  height = 9
)

# ------------------------------------------------------------
# 16. 尝试读取作者原始 annotation
# ------------------------------------------------------------
author_annotation_col <- NA_character_
author_annotation_available <- FALSE

if (file.exists(step0_meta_file)) {

  log_msg(
    "Step 0 confirmed metadata found; checking for author annotation ..."
  )

  raw_meta <- fread(
    step0_meta_file,
    data.table = FALSE
  )

  if ("barcode" %in% colnames(raw_meta)) {

    original_cols <- grep(
      "^original_V",
      colnames(raw_meta),
      value = TRUE
    )

    keyword_pattern <- paste(
      c(
        "microglia",
        "myeloid",
        "macroph",
        "mono",
        "neutro",
        "astro",
        "oligo",
        "endothel",
        "fibro",
        "pericy",
        "neuron",
        "epend",
        "mast",
        "t.cell",
        "b.cell",
        "lymph"
      ),
      collapse = "|"
    )

    if (length(original_cols) > 0) {

      candidate_scores <- sapply(
        original_cols,
        function(cc) {

          vals <- unique(
            na.omit(
              as.character(
                raw_meta[[cc]]
              )
            )
          )

          vals <- vals[
            nzchar(vals)
          ]

          nuniq <- length(vals)

          if (nuniq < 2 || nuniq > 100) {
            return(-Inf)
          }

          sum(
            grepl(
              keyword_pattern,
              vals,
              ignore.case = TRUE
            )
          )
        }
      )

      if (
        length(candidate_scores) > 0 &&
        max(candidate_scores, na.rm = TRUE) >= 2
      ) {

        author_annotation_col <- names(
          candidate_scores
        )[
          which.max(candidate_scores)
        ]

        author_annotation_available <- TRUE

        log_msg(
          "Detected probable author annotation column: ",
          author_annotation_col
        )

        idx_raw <- match(
          colnames(pbmc),
          raw_meta$barcode
        )

        if (anyNA(idx_raw)) {

          warning(
            "作者 metadata 有 ",
            sum(is.na(idx_raw)),
            " 个 barcode 无法映射。"
          )

        } else {

          pbmc$author_annotation <- as.character(
            raw_meta[[author_annotation_col]][idx_raw]
          )

          author_ct <- as.data.frame(
            table(
              cluster = pbmc$cluster_provisional,
              author_annotation = pbmc$author_annotation
            ),
            stringsAsFactors = FALSE
          )

          colnames(author_ct)[3] <- "n_cells"

          fwrite(
            author_ct,
            file.path(
              tab_dir,
              "05_author_annotation_by_cluster.csv"
            )
          )

          # 作者 annotation UMAP
          p_author <- DimPlot(
            pbmc,
            reduction = "umap.harmony",
            group.by = "author_annotation",
            label = TRUE,
            repel = TRUE,
            pt.size = 0.08
          ) +
            ggtitle(
              paste0(
                "Original author annotation (",
                author_annotation_col,
                ")"
              )
            )

          ggsave(
            file.path(
              fig_dir,
              "08_author_annotation_UMAP.pdf"
            ),
            p_author,
            width = 12,
            height = 8
          )

          # cluster 内 author annotation 比例
          author_prop <- author_ct

          cluster_totals_author <- aggregate(
            n_cells ~ cluster,
            data = author_prop,
            FUN = sum
          )

          colnames(cluster_totals_author)[2] <- "cluster_total"

          author_prop <- merge(
            author_prop,
            cluster_totals_author,
            by = "cluster",
            all.x = TRUE
          )

          author_prop$prop_within_cluster <- (
            author_prop$n_cells /
              author_prop$cluster_total
          )

          p_author_heat <- ggplot(
            author_prop,
            aes(
              x = author_annotation,
              y = cluster,
              fill = prop_within_cluster
            )
          ) +
            geom_tile() +
            labs(
              x = NULL,
              y = "Cluster",
              fill = "Fraction",
              title = "Current clusters vs original author annotation"
            ) +
            theme_classic() +
            theme(
              axis.text.x = element_text(
                angle = 55,
                hjust = 1
              )
            )

          ggsave(
            file.path(
              fig_dir,
              "09_cluster_vs_author_annotation_heatmap.pdf"
            ),
            p_author_heat,
            width = 13,
            height = 9
          )
        }
      }
    }
  }
}

# ------------------------------------------------------------
# 17. Candidate microglia evidence
# ------------------------------------------------------------
microglia_markers <- marker_sets$Microglia
macrophage_markers <- marker_sets$Macrophage_Monocyte
neutrophil_markers <- marker_sets$Neutrophil

micro_score <- setNames(
  program_score_df$Microglia,
  program_score_df$cluster
)

macro_score <- setNames(
  program_score_df$Macrophage_Monocyte,
  program_score_df$cluster
)

neut_score <- setNames(
  program_score_df$Neutrophil,
  program_score_df$cluster
)

micro_hits <- sapply(
  cluster_ids,
  function(cl) {

    genes_cl <- top50_markers$gene[
      as.character(top50_markers$cluster) == cl
    ]

    sum(
      tolower(genes_cl) %in%
        tolower(microglia_markers)
    )
  }
)

macro_hits <- sapply(
  cluster_ids,
  function(cl) {

    genes_cl <- top50_markers$gene[
      as.character(top50_markers$cluster) == cl
    ]

    sum(
      tolower(genes_cl) %in%
        tolower(macrophage_markers)
    )
  }
)

neut_hits <- sapply(
  cluster_ids,
  function(cl) {

    genes_cl <- top50_markers$gene[
      as.character(top50_markers$cluster) == cl
    ]

    sum(
      tolower(genes_cl) %in%
        tolower(neutrophil_markers)
    )
  }
)

candidate_df <- data.frame(
  cluster = cluster_ids,
  broad_type_suggestion = cluster_to_type[cluster_ids],
  Microglia_score = as.numeric(
    micro_score[cluster_ids]
  ),
  Macrophage_Monocyte_score = as.numeric(
    macro_score[cluster_ids]
  ),
  Neutrophil_score = as.numeric(
    neut_score[cluster_ids]
  ),
  microglia_markers_in_top50 = as.integer(
    micro_hits[cluster_ids]
  ),
  macrophage_markers_in_top50 = as.integer(
    macro_hits[cluster_ids]
  ),
  neutrophil_markers_in_top50 = as.integer(
    neut_hits[cluster_ids]
  ),
  stringsAsFactors = FALSE
)

# strict：
#   Microglia program 为最高 broad suggestion
#   + microglia marker 至少 2 个进入 top50
candidate_df$strict_microglia_candidate <- (
  candidate_df$broad_type_suggestion == "Microglia" &
    candidate_df$microglia_markers_in_top50 >= 2
)

# extended：
#   保留可能因激活而丢失 homeostatic markers 的 myeloid clusters
#   条件：
#     microglia/myeloid/macrophage 任一方向较强，
#     但不能是明显 neutrophil
candidate_df$extended_myeloid_candidate <- (
  (
    candidate_df$broad_type_suggestion %in%
      c(
        "Microglia",
        "Myeloid_general",
        "Macrophage_Monocyte"
      ) |
      candidate_df$Microglia_score > 0 |
      candidate_df$microglia_markers_in_top50 >= 1
  ) &
    candidate_df$broad_type_suggestion != "Neutrophil" &
    candidate_df$Neutrophil_score <
      max(
        candidate_df$Neutrophil_score,
        na.rm = TRUE
      )
)

candidate_df$candidate_level <- "No"

candidate_df$candidate_level[
  candidate_df$extended_myeloid_candidate
] <- "Extended_myeloid_review"

candidate_df$candidate_level[
  candidate_df$strict_microglia_candidate
] <- "Strict_microglia"

# 若存在作者 annotation，则加入 cluster 中作者 Microglia 占比
if ("author_annotation" %in% colnames(pbmc@meta.data)) {

  author_micro_prop <- sapply(
    cluster_ids,
    function(cl) {

      vals <- as.character(
        pbmc$author_annotation[
          pbmc$cluster_provisional == cl
        ]
      )

      if (length(vals) == 0) {
        return(NA_real_)
      }

      mean(
        grepl(
          "microglia",
          vals,
          ignore.case = TRUE
        ),
        na.rm = TRUE
      )
    }
  )

  candidate_df$author_microglia_fraction <- as.numeric(
    author_micro_prop[cluster_ids]
  )

  # 作者 annotation 高度支持时，至少进入 extended review
  upgrade <- (
    !is.na(candidate_df$author_microglia_fraction) &
      candidate_df$author_microglia_fraction >= 0.20 &
      candidate_df$candidate_level == "No"
  )

  candidate_df$candidate_level[upgrade] <- (
    "Extended_myeloid_review"
  )
}

fwrite(
  candidate_df,
  file.path(
    tab_dir,
    "06_candidate_microglia_cluster_evidence.csv"
  )
)

strict_clusters <- candidate_df$cluster[
  candidate_df$candidate_level ==
    "Strict_microglia"
]

extended_clusters <- candidate_df$cluster[
  candidate_df$candidate_level ==
    "Extended_myeloid_review"
]

writeLines(
  c(
    "STRICT MICROGLIA CANDIDATES:",
    if (length(strict_clusters) > 0) {
      paste(strict_clusters, collapse = ", ")
    } else {
      "[none]"
    },
    "",
    "EXTENDED MYELOID CLUSTERS TO REVIEW:",
    if (length(extended_clusters) > 0) {
      paste(extended_clusters, collapse = ", ")
    } else {
      "[none]"
    },
    "",
    "IMPORTANT:",
    "Do not extract microglia automatically from this file.",
    "Final selection must be made after reviewing DotPlot, top markers, UMAP, and author annotation."
  ),
  file.path(
    tab_dir,
    "07_candidate_microglia_clusters.txt"
  )
)

# ------------------------------------------------------------
# 18. Candidate microglia / myeloid by sample
# ------------------------------------------------------------
candidate_keep <- unique(
  c(
    strict_clusters,
    extended_clusters
  )
)

if (length(candidate_keep) > 0) {

  candidate_cells <- colnames(pbmc)[
    pbmc$cluster_provisional %in%
      candidate_keep
  ]

  tmp <- pbmc@meta.data[
    candidate_cells,
    ,
    drop = FALSE
  ]

  candidate_by_sample <- as.data.frame(
    table(
      cluster = tmp$cluster_provisional,
      sample_id = tmp$sample_id
    ),
    stringsAsFactors = FALSE
  )

  colnames(candidate_by_sample)[3] <- "n_cells"

  fwrite(
    candidate_by_sample,
    file.path(
      tab_dir,
      "08_candidate_microglia_by_sample.csv"
    )
  )

  # candidate level 映射到细胞
  cluster_to_candidate_level <- setNames(
    candidate_df$candidate_level,
    candidate_df$cluster
  )

  pbmc$microglia_candidate_level <- unname(
    cluster_to_candidate_level[
      pbmc$cluster_provisional
    ]
  )

  p_candidate <- DimPlot(
    pbmc,
    reduction = "umap.harmony",
    group.by = "microglia_candidate_level",
    pt.size = 0.08
  ) +
    ggtitle(
      "Microglia / myeloid candidate clusters for manual review"
    )

  ggsave(
    file.path(
      fig_dir,
      "11_candidate_microglia_UMAP.pdf"
    ),
    p_candidate,
    width = 10,
    height = 8
  )
}

# ------------------------------------------------------------
# 19. Top5 marker heatmap（每 cluster 最多抽 250 cells）
# ------------------------------------------------------------
top5_markers <- pbmc_markers %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  slice_head(n = 5) %>%
  ungroup()

heat_genes <- unique(top5_markers$gene)
heat_genes <- intersect(
  heat_genes,
  rownames(pbmc)
)

if (length(heat_genes) > 1) {

  heat_cells <- unlist(
    lapply(
      cluster_ids,
      function(cl) {

        cells_cl <- colnames(pbmc)[
          pbmc$cluster_provisional == cl
        ]

        if (length(cells_cl) <= 250) {
          cells_cl
        } else {
          sample(
            cells_cl,
            250
          )
        }
      }
    ),
    use.names = FALSE
  )

  pbmc_heat <- subset(
    pbmc,
    cells = heat_cells
  )

  Idents(pbmc_heat) <- "cluster_provisional"

  pbmc_heat <- ScaleData(
    pbmc_heat,
    features = heat_genes,
    verbose = FALSE
  )

  p_heat <- DoHeatmap(
    pbmc_heat,
    features = heat_genes,
    group.by = "cluster_provisional",
    raster = TRUE
  ) +
    ggtitle(
      "Top 5 markers per broad cluster"
    )

  ggsave(
    file.path(
      fig_dir,
      "10_top5_cluster_marker_heatmap.pdf"
    ),
    p_heat,
    width = 15,
    height = 12
  )

  rm(pbmc_heat)
  gc()
}

# ------------------------------------------------------------
# 20. 保存带 suggestion 的对象
# ------------------------------------------------------------
saveRDS(
  pbmc,
  file.path(
    out_dir,
    "GSE162610_broad_annotation_audit.rds"
  ),
  compress = FALSE
)

# ------------------------------------------------------------
# 21. 报告
# ------------------------------------------------------------
n_markers <- nrow(pbmc_markers)
n_strict <- length(strict_clusters)
n_extended <- length(extended_clusters)

report <- c(
  "============================================================",
  "GSE162610 STEP 1F | BROAD ANNOTATION AUDIT REPORT",
  "============================================================",
  paste0("Date: ", Sys.time()),
  "",
  paste0("Input cells: ", ncol(pbmc)),
  paste0("Genes: ", nrow(pbmc)),
  paste0("Broad clusters (resolution 0.4): ", length(cluster_ids)),
  paste0("FindAllMarkers rows: ", n_markers),
  paste0("Marker runtime (min): ", round(marker_runtime, 2)),
  "",
  "Broad-lineage suggestion by cluster:",
  paste0(
    "  cluster ",
    broad_suggestion_df$cluster,
    " -> ",
    broad_suggestion_df$broad_type_suggestion,
    " | margin=",
    round(
      broad_suggestion_df$score_margin,
      3
    ),
    " | ",
    broad_suggestion_df$suggestion_confidence
  ),
  "",
  paste0(
    "Probable author annotation available: ",
    author_annotation_available
  ),
  paste0(
    "Detected author annotation column: ",
    ifelse(
      is.na(author_annotation_col),
      "[none]",
      author_annotation_col
    )
  ),
  "",
  paste0(
    "Strict microglia candidates: ",
    n_strict,
    " | ",
    ifelse(
      n_strict > 0,
      paste(strict_clusters, collapse = ", "),
      "[none]"
    )
  ),
  paste0(
    "Extended myeloid clusters for review: ",
    n_extended,
    " | ",
    ifelse(
      n_extended > 0,
      paste(extended_clusters, collapse = ", "),
      "[none]"
    )
  ),
  "",
  "IMPORTANT:",
  "broad_type_suggestion is algorithm-assisted and is NOT the final manual annotation.",
  "Microglia must be separated from infiltrating macrophage/monocyte using multiple markers, cluster markers, time distribution, and original annotation if available.",
  "Activated microglia may downregulate P2ry12/Tmem119, so extended myeloid clusters must also be reviewed.",
  "",
  "STEP 1F COMPLETE = TRUE",
  "",
  "NEXT AFTER REVIEW:",
  "Finalize broad annotation and define clean resident-microglia candidate clusters.",
  "Then proceed to microglia-only re-clustering.",
  "============================================================"
)

writeLines(
  report,
  file.path(
    out_dir,
    "12_broad_annotation_report.txt"
  )
)

log_msg("Step 1F complete.")

cat("\n")
cat(paste(report, collapse = "\n"))
cat("\n\n")
cat(
  "✅ Step 1F 完成。请把关键表格和图片发给我判断，不要自动提取 microglia。\n"
)
