# Rustscenic multiomic pipeline

**RUSTSCENIC** (Single-Cell rEgulatory Network Inference and Clustering) is a modern, high-performance tool designed to significantly enhance the speed, scalability, and reproducibility of the single-cell RNA-sequencing data processing. For more detailed information about the tool, please refer to the [developer's repository](https://github.com/Ekin-Kahraman/rustscenic).

This repository contains reports and results from running rustscenic on various multiomic datasets.

The code for the end-to-end multiome pipeline execution was adapted from this [script](https://github.com/Ekin-Kahraman/rustscenic/blob/main/validation/validate_multiome_e2e.py). Analysis was performed using publicly available datasets, specifically from [10x Genomics](https://www.10xgenomics.com).

## Performance

|      Dataset       | Cells |          Features (RNA/ATAC)        | TFs|Edges    | Regulons|GRN(wall,s)| AUCell(wall,s)| Topics (wall,s)|Total Time|ARI   |Top Markers	     |Full log    |
| ------------------ |-------|-------------------------------------|----|---------|---------|-----------|---------------|----------------|----------|------|-----------------|------------|
| PBMC 3k            |  2767 | RNA(2767, 21335) ATAC(451378, 81156)|1467|2,367,067|   187   |391.9      |0.4            |203.6           |595.9s    |0.1190| SPIB, EBF1, PAX5|[log file]()|
