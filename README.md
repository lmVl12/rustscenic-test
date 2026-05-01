# Rustscenic multiomic pipeline

**RUSTSCENIC** (Single-Cell rEgulatory Network Inference and Clustering) is a modern, high-performance tool designed to significantly enhance the speed, scalability, and reproducibility of the single-cell RNA-sequencing data processing. For more detailed information about the tool, please refer to the [developer's repository](https://github.com/Ekin-Kahraman/rustscenic).

This repository contains reports and results from running rustscenic on various multiomic datasets.

The code for the **end-to-end multiome pipeline (RNA+ATACseq)** execution was adapted from this [script](https://github.com/Ekin-Kahraman/rustscenic/blob/main/validation/validate_multiome_e2e.py). Analysis was performed using publicly available datasets, specifically from [10x Genomics](https://www.10xgenomics.com).

## Performance
*Additional packages must be installed manually, as they are missing from the base setup (v0.3.3):

```
pip install igraph leidenalg
```
#### Pipeline runs

| Release  |      Dataset          |Cells |          Features (RNA/ATAC)          | TFs|Edges    | Regulons|GRN(wall,s)| AUCell(wall,s)| Topics (wall,s)|Total Time|RAM       |ARI   |Top Markers	|Full log    |
| ---------| --------------------- |------|---------------------------------------|----|---------|---------|-----------|---------------|----------------|----------|----------|------|------------------|------------|
| v0.3.3   | pbmc3k(smoke test e2e)|2,767 |RNA(2767, 21335) ATAC(451378, 81156)   |1467|2,367,067|   187   |   403.4   |0.4            |237.6           |641.4s    | 3.19 GB  |0.1190| SPIB, EBF1, PAX5 |[log file](https://github.com/lmVl12/rustscenic-test/blob/main/logs/pbmc3k_v0.3.3.log)|
| v0.3.3   | Human Brain Neclei 10k|10,198|RNA(10211, 32856) ATAC (644847, 130127)|1760|4,497,808|   250   |  2892.5   |3.3            |1510.5          |4406.3s   | 10.78 GB |0.5140|IKZF1, FLI1, CEBPD|[log file](https://github.com/lmVl12/rustscenic-test/blob/main/logs/pbmc3k_v0.3.3.log)|


## Runtime Environment (WSL2 Configuration)
- Virtual OS: Ubuntu 22.04.5 LTS LTS on WSL2
- Allocated RAM: 24.0 GB (out of 32.0 GB)
- Compute Threads: 16 Logical Processors (mapped to i7-14650HX)

## Structure
File Descriptions:
- `00_` prefixed files: Initial test runs used to verify the software installation and basic environment compatibility. Adapted from [example script](https://github.com/Ekin-Kahraman/rustscenic/blob/main/examples/pbmc3k_end_to_end.py).
- `prep_` scripts: Dedicated modules for raw data cleaning and normalization before the main integration. Adapted from [atac_fragments_to_matrix.py](https://github.com/Ekin-Kahraman/rustscenic/blob/main/examples/atac_fragments_to_matrix.py) and [prep_pbmc10k.py](https://github.com/Ekin-Kahraman/rustscenic/blob/main/validation/prep_pbmc10k.py).
- `multiome_` : The core pipeline.
- [reference_](https://github.com/Ekin-Kahraman/rustscenic/blob/main/validation/reference_pipeline_multiome.py) & [compare](https://github.com/Ekin-Kahraman/rustscenic/blob/main/validation/compare_pipelines_multiome.py)& [Dockerfile](https://github.com/Ekin-Kahraman/rustscenic/blob/main/validation/reference/Dockerfile): Included for future comparative testing. (also to be adapted from the developer's repo).
- `logs/`: Detailed execution summaries.
