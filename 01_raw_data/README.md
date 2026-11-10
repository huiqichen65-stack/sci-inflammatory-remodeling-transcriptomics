# External data placement

Raw sequencing matrices and large third-party reference resources are not stored in this repository.

For maximum compatibility with the frozen analysis scripts, place downloaded files under this directory using the project layout below:

```text
01_raw_data/
├── GSE162610/
├── GSE234774/
├── pySCENIC_resources/
└── reference_gene_sets/
```

The discovery scripts expect the public GSE162610 matrix/metadata files from GEO. The independent-validation and spatial scripts expect the GSE234774 snRNA-seq and processed spatial files from GEO.

The regulatory-network workflow additionally requires the mouse cisTarget resources referenced in the Step 7 scripts, including the mouse TF list, mm10 ranking database, and motif annotation table. Sensitivity analyses that use the Marsh et al. dissociation signature and offline KEGG mappings require the corresponding external reference files.

Accession numbers and study roles are listed in `../public_datasets.csv`.
