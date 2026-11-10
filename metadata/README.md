# Public cell-level annotations

This directory contains cell-level annotations released with the analysis code.

## GSE162610 discovery dataset

`GSE162610_cell_annotations.csv.gz` contains one row for each of the 62,842 atlas cells retained after doublet removal. It links the cell barcode to the broad atlas annotation and, where applicable, the final myeloid lineage, resident-microglial state, transition status, trajectory pseudotime, and formal branch label.

Supporting provenance files document field definitions, mapping rules, source-file hashes, and export audits. The frozen formal branch labels are:

- `Shared_prebranch`
- `Alternative`
- `Target`

The raw branch-arm labels are also retained so that the public labels remain traceable to the original trajectory analysis.

## GSE234774 independent cohort

`GSE234774_myeloid_cell_annotations.csv.gz` contains the 25,205 nuclei in the project myeloid comparison table used to relate the frozen discovery program to published Tabulae Paralytica annotations and the project myeloid labels. Published layer-2 through layer-5 annotations are retained alongside the project `myeloid_final_label` and `resident_lineage_status`.

This file is annotation-focused: count depth, feature depth, and program-score columns from the internal comparison table are intentionally omitted because they are not needed to satisfy the public cell-annotation release.

## Scope

These files contain annotations and identifiers, not expression matrices. Raw sequencing matrices and large intermediate RDS/loom/BPCells objects are not included in this repository.
