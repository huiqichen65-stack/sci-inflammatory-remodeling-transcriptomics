# Environment notes

The main R analyses were run with R 4.5.2. Package versions recovered from frozen analysis sessions are summarized in `../software_versions.csv`.

The SCENIC workflow used an isolated Python 3.10 environment with pySCENIC 0.12.1, arboreto 0.1.6, and ctxcore 0.2.0. The original environment specification files are retained here.

The `scripts/helpers/environment/` directory contains the environment/resource setup and compatibility helpers that were used during the frozen regulatory-network analysis. These helpers are separated from the biological analysis scripts because they are environment-specific support code.

External cisTarget ranking databases and motif resources are not redistributed in this repository.
