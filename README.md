# Rustscenic multiomic pipeline

**RUSTSCENIC** (Single-Cell rEgulatory Network Inference and Clustering) is a modern, high-performance tool designed to significantly enhance the speed, scalability, and reproducibility of the single-cell RNA-sequencing data processing. For more detailed information about the tool, please refer to the [developer's repository](https://github.com/Ekin-Kahraman/rustscenic).

This repository contains reports and results from running rustscenic on various multiomic datasets.

The code for the **end-to-end multiome pipeline (RNA+ATACseq)** execution was adapted from this [script](https://github.com/Ekin-Kahraman/rustscenic/blob/main/validation/validate_multiome_e2e.py). Analysis was performed using publicly available datasets, specifically from [10x Genomics](https://www.10xgenomics.com).

## Performance

\*Additional packages must be installed manually, as they are missing from the base setup (v0.3.3) - fixed in v0.3.6:

```
pip install igraph leidenalg
```

#### Pipeline runs

|Release|Dataset|Cells, Features (RNA/ATAC)|TFs, Edges|Regulons/ eRegulons/ unique|Total time/ Per stage, s|RAM, GB| Markers|Metrics|Report, Notes&Findings|
|-|-|-|-|-|-|-|-|-|-|
|v0.3.3|pbmc3k (e2e)|Cells: 2,767 RNA (2767, 21335) ATAC (451378, 81156)|1,467/ 2,367,067|187|641.4 <br />GRN: 403.4<br />AUCell: 0.4<br />Topics: 237.6|3.19 |SPIB, EBF1, PAX5|ARI: 0.1190|[log file](https://github.com/lmVl12/rustscenic-test/blob/main/results/pbmc3k_e2e_v0.3.3/pbmc3k_v0.3.3.log)<br /> Initial test on known data|
|v0.3.3|Human Brain Nuclei 10k raw(e2e)|Cells: 10,198RNA (10211, 32856) ATAC (644847, 130127)|1,760/ 4,497,808|1,760/<br />250 unique|4406.3 <br />GRN: 2892.5<br />AUCell: 3.3<br />Topics: 1510.5|10.78 |IKZF1, FLI1, CEBPD|ARI: 0.5140|[log file](https://github.com/lmVl12/rustscenic-test/blob/main/results/human_brain10k_e2e_unfiltered_v0.3.3/human_brain10k_unfiltered_v0.3.3.log)<br />raw feature  matrix is used as an input (standard QC) - microglia markers recovered|
|v0.3.6|Human Brain Nuclei 10k raw (e2e)|Cells: 10,198 RNA (10211, 32856) ATAC (10198,127370)|1,760/ 4,497,808|1760/ 250 unique|4090.1 <br />GRN: 2749.2<br />AUCell: 3.93.1<br />Topics: 1337.8|10.72|IKZF1, FLI1, TFEC, HCLS1, NFATC2, RUNX1|ARI: 0.5139|[log file](https://github.com/lmVl12/rustscenic-test/blob/main/results/human_brain10k_e2e_v0.3.6/human_brain10k_unfiltered_v0.3.6.log) <br /> raw feature  matrix is used as an input (standard QC). Optimized time compared to v0.3.3, biology unchanged|
|v0.3.6|Human Brain Nuclei 10k filtered (e2e)|Cells: 9,665 RNA (9665, 32317) ATAC (9665, 123089)|1,748/ 4,467,176|1,748/ 120 unique|4777.7 <br />GRN: 3474.1<br />AUCell: 3.9<br />Topics: 1299.7|9.89|FLI1, IKZF1, RUNX1, HCLS1, SMAP2 |ARI: 05751|[log file](https://github.com/lmVl12/rustscenic-test/blob/main/results/human_brain10k_e2e_v0.3.6/human_brain10k_v0.3.6.log) <br /> used filtered feature matrix as an input (standard QC). The results of raw vs. filtered matrix are similar when standard QC is preformed. Decreased cell count did not speed up the process (longer GRN stage)|
|v0.4.1|Human Brain Nuclei 10k(e2e)|Cells: 8,215 RNA (8215, 32317) ATAC (8215,123089)|1748/ 4,293,902|1,748|GRN: 2144.098<br />AUCell: 1.916<br />Topics: 1056.69|9.08|SRRM3, CELF4, ZMAT4, ADARB1, BCL11A, ZBTB20|ARI: 0.45<br />Recovery: 0.1176 (\*from top-10 only)|[log file](https://github.com/lmVl12/rustscenic-test/blob/main/results/human_brain10k_e2e_v0.4.1/human_brain_10k_v0.4.1.log), [json file](https://github.com/lmVl12/rustscenic-test/blob/main/results/human_brain10k_e2e_v0.4.1/human_brain_10k_v0.4.1.json) <br /> used filtered feature matrix; additional subsetting to reduce microglia signal applied|
|v0.4.1|Lymphoma 14k (e2e)|Cells: 14,039 RNA (14039, 70132) ATAC (14039,70132)|1663/ 3,061,291|1,663|GRN: 2526.88<br />AUCell: 4.137<br />Topics: 954.58|8.4|GRHPR, POU2F2, PAX5, MEF2B, SPIB|ARI: 0.096|[log_file](https://github.com/lmVl12/rustscenic-test/blob/main/results/lymphoma14k_e2e_v0.4.1/lymphoma_14k_results_v0.4.1.log), [json_file](https://github.com/lmVl12/rustscenic-test/blob/main/results/lymphoma14k_e2e_v0.4.1/lymphoma_14k_results_v0.4.1.json) <br /> Natural high homigenicity of a sample might explain the low ARI|
|v0.4.6|Human Brain Nuclei 10k(full pipeline)|Cells: 10,014RNA (8215,32808) ATAC (8215,123089)|1693/ 4,314,539|1,693/ 1,693 eregulons|3295.42<br />GRN: 2071.98<br />AUCell: 2.488<br />Topics: 587.3496<br />Cistarget: 4.21<br />Enhancer: 188.38<br />Eregulons: 200.84|24,99|NEUROD2,           NEUROD6,           TBR1,           FOXG1,          BCL11A|Recovery: 0.9412 (16/17)|[json_file](https://github.com/lmVl12/rustscenic-test/blob/main/results/human_brain10k_full_pipeline_v0.4.6/human_brain_10k_v0.4.6.json)|



\*ARI is grn-based cell-type clustering (via top regulon activity) vs ATAC leiden

## Runtime Environment (WSL2 Configuration)

* Virtual OS: Ubuntu 22.04.5 LTS LTS on WSL2
* Allocated RAM: 24.0 GB out of 32.0 GB. (additional Swap memory was provisioned to support the full pipeline execution)
* Compute Threads: 16 Logical Processors (mapped to i7-14650HX)
* Python 3.12.13

## Structure

File Descriptions:
* `scripts\`:
  * `00_` prefixed files: Initial test runs used to verify the software installation and basic environment compatibility. Adapted from [example script](https://github.com/Ekin-Kahraman/rustscenic/blob/main/examples/pbmc3k_end_to_end.py).
  * `prep_` scripts: Dedicated modules for raw data cleaning and normalization before the main integration. Adapted from [atac\_fragments\_to\_matrix.py](https://github.com/Ekin-Kahraman/rustscenic/blob/main/examples/atac_fragments_to_matrix.py) and [prep\_pbmc10k.py](https://github.com/Ekin-Kahraman/rustscenic/blob/main/validation/prep_pbmc10k.py).
* `results\`: Detailed execution summaries, logs and output files.

