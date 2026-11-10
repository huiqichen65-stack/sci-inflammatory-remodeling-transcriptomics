# Injury-associated inflammatory-remodeling program after spinal cord injury

This repository contains the analysis code and public cell-level annotations used for the revised manuscript **“Integrated single-cell and spatial transcriptomics identify an injury-associated inflammatory-remodeling program detectable near spinal cord lesions at chronic stages.”**

The study reanalyzes public mouse spinal cord injury datasets to define an acute microglia-enriched inflammatory-remodeling transcriptional program in GSE162610 and assess the frozen program in independent GSE234774 single-nucleus and spatial data. The workflow includes single-cell atlas construction, myeloid-lineage adjudication, resident-state analysis, Monocle 3 trajectory analysis, branch-aware differential expression, SCENIC, CellChat, independent snRNA-seq validation, RCTD-based spatial analysis, sensitivity analyses, and final figure generation.

## Data availability

All source datasets analyzed in the study are publicly available. Accession numbers and their roles are summarized in [`public_datasets.csv`](public_datasets.csv).

The primary datasets are:

- **GSE162610** — discovery scRNA-seq dataset used for the acute atlas, myeloid refinement, state/trajectory analyses, SCENIC, and CellChat.
- **GSE234774** — independent snRNA-seq and spatial transcriptomic dataset used for external program assessment and spatial analyses.

Raw sequencing matrices and large intermediate analysis objects are not duplicated in this repository. The raw files should be downloaded from GEO and placed using the layout described in [`01_raw_data/README.md`](01_raw_data/README.md).

Public cell-level annotations requested for reproducibility are provided in [`metadata`](metadata):

- `GSE162610_cell_annotations.csv.gz` — 62,842 atlas cells with barcode, broad annotation, myeloid lineage, resident state, transition status, trajectory information, pseudotime, and formal branch labels where applicable.
- `GSE234774_myeloid_cell_annotations.csv.gz` — 25,205 nuclei from the project myeloid comparison table with published Tabulae Paralytica annotation layers and project myeloid/resident-lineage labels.

## Directory structure

Analysis scripts are organized by analytical component. Original numerical prefixes are retained where possible so that the released code remains traceable to the frozen analysis history.

1. [`scripts/atlas`](scripts/atlas) — GSE162610 input/metadata checks, quality control, doublet removal, normalization, Harmony integration, clustering, and broad atlas annotation. These analyses provide the basis for **Figure 1A–D** and Supplementary Figure S1.
2. [`scripts/myeloid`](scripts/myeloid) — myeloid refinement, lineage adjudication, resident-microglia isolation, state scoring, activated/DAM-like transition mapping, temporal differential expression, and candidate-program definition. These analyses provide the basis for **Figure 1E–H**, Figure 2A–B, and related supplementary analyses.
3. [`scripts/trajectory`](scripts/trajectory) — balanced root construction, Monocle 3 trajectory reconstruction, branch definition, progressive-divergence analysis, branch-specific differential expression, dynamic modules, enrichment, and pseudobulk ssGSEA. These analyses provide the basis for **Figure 2C–F** and trajectory-related supplementary figures/tables.
4. [`scripts/regulatory`](scripts/regulatory) — pySCENIC input preparation, GRNBoost2, cisTarget pruning, AUCell, RSS/CSI, regulon timing, and target-program regulatory closure. These analyses provide the basis for **Figure 3** and regulatory components of Figure 4.
5. [`scripts/cellchat`](scripts/cellchat) — CellChat input construction, the frozen time-course visualization run, the 7-dpi nboot=100 confirmation analysis, recipient-program integration, and mechanism-closure analyses. These analyses provide the basis for **Figure 4A–E**.
6. [`scripts/validation`](scripts/validation) — GSE234774 snRNA-seq input/QC, doublet removal, normalization/clustering, published-guided lineage annotation, frozen myeloid labels, AUCell/ssGSEA validation, and comparison with published myeloid annotations. These analyses provide the basis for **Figure 4F** and Supplementary Figure S8.
7. [`scripts/spatial`](scripts/spatial) — GSE234774 spatial input audit, capture-area asset linkage, QC, RCTD reference construction and mapping, target/shared-program spatial scoring, lesion-zone analyses, recipient-neighborhood closure, and registered-coordinate anatomical anchoring. These analyses provide the basis for **Figure 5** and spatial supplementary analyses.
8. [`scripts/sensitivity`](scripts/sensitivity) — lineage specificity, dissociation-signature exclusion, trajectory robustness, full-database CellChat sensitivity, receptor-depth/Tnf-source sensitivity, matched-null validation, and spatial cell-source sensitivity.
9. [`scripts/figures`](scripts/figures) — scripts used to generate the final main and supplementary figures from frozen upstream outputs.
10. [`scripts/helpers`](scripts/helpers) — public-metadata export, repository validation, pySCENIC Python helpers, and environment/compatibility support code.

The original-to-repository mapping is recorded in [`script_manifest.csv`](script_manifest.csv). Gaps in the original numbering correspond primarily to debugging, checkpoint-recovery, obsolete, or superseded scripts that are intentionally not part of the public analysis workflow.

## Running the analyses

Run scripts from the repository root. The project root can be set explicitly:

```r
Sys.setenv(GLIA_PROJECT_ROOT = normalizePath(".", winslash = "/"))
```

or in a shell before calling `Rscript`.

The GSE162610 workflow begins with:

```text
scripts/atlas/01_step0A_GSE162610_audit_manifest.R
```

and then proceeds approximately by the original numerical prefixes through the atlas, myeloid, trajectory, regulatory, and CellChat modules.

The independent GSE234774 workflow begins with:

```text
scripts/validation/51_step9A_GSE234774_rnaseq_input_metadata_audit.R
```

The spatial branch begins with:

```text
scripts/spatial/63_step12A_GSE234774_spatial_input_section_metadata_audit.R
```

The major workflow is:

```text
GSE162610 atlas
      |
      v
myeloid lineage + resident states
      |
      v
transition + trajectory + branch program
      |-------------------|
      v                   v
SCENIC                  CellChat
      \                   /
       \                 /
        ---- target program
                 |
                 v
       GSE234774 snRNA-seq
                 |
                 v
       GSE234774 spatial analysis
```

Downstream modules are complementary branches rather than one strictly serial script chain. Each script contains its own expected upstream files and frozen-count checks.

## Path configuration

Machine-specific project paths have been removed from the public copies. Released R scripts resolve the project root from `GLIA_PROJECT_ROOT`, with the repository root (`.`) as the default.

[`config.R`](config.R) documents the project, raw-data, result, figure, reference, and Python path conventions. For maximum compatibility with the frozen scripts, use the historical project subdirectories:

```text
01_raw_data/
04_results/
09_final_figures/
```

Large generated objects under `04_results/` and rendered figures under `09_final_figures/` are intentionally excluded from Git.

## Reproducibility information

Key analytical parameters are summarized in [`analysis_parameters.csv`](analysis_parameters.csv), observed software/package versions in [`software_versions.csv`](software_versions.csv), and explicit random-seed statements in [`random_seeds.csv`](random_seeds.csv).

The regulatory-network workflow used an isolated **Python 3.10 / pySCENIC 0.12.1 / arboreto 0.1.6 / ctxcore 0.2.0** environment. Environment specifications and the compatibility helpers used during the frozen run are retained under [`environment`](environment) and [`scripts/helpers/environment`](scripts/helpers/environment).

A repository-level validation script is provided:

```text
Rscript scripts/helpers/verify_release.R
```

It checks R-script parsing, required release files, absence of known private workstation paths, annotation row counts/uniqueness, and the absence of large analysis-object formats. It does **not** rerun the biological analyses.

## Notes on frozen analysis versions

The repository is intended to reproduce the analysis underlying the revised manuscript, not every exploratory or debugging run.

- The released Step 7B script is the final GRNBoost2 run corresponding to the frozen Step 7B outputs.
- The released Step 8A script is the final annotation-corrected CellChat preflight version corresponding to the frozen Step 8A outputs.
- The nboot=20 CellChat time-course script is retained because frozen Figure 4B uses the corresponding time-course CellChat objects; the prespecified 7-dpi confirmation analysis used nboot=100.
- Both the original Step 12F lesion-zone analysis and the revision-specific Step 12F within-section update are retained because they support different frozen manuscript outputs.
- Debug-only scripts, checkpoint-recovery scripts, false-stop repair scripts, abandoned histology-manual-annotation work, and superseded figure versions are not included in the main release.

Trajectory, regulatory-network, ligand-receptor, and spatial analyses are associative frameworks. The released code preserves the interpretation boundaries stated in the revised manuscript and does not imply developmental fate, cellular persistence, or causal signaling.
