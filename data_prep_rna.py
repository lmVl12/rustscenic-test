"""Preprocess for GRN inference."""
from pathlib import Path
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
# Add a dictionary of datasets (to be changed every time)

CURRENT_RUN = "human_brain_10k"
LOG_DIR = Path("logs")
LOG_DIR.mkdir(exist_ok=True)
DATASETS = {
    "pbmc_3k": {
        "raw_h5": "data/pbmc_unsorted_3k_filtered_feature_bc_matrix.h5",
        "rna_out": "data/pbmc_3k_rna.h5ad",
        "qc": {
            "min_genes": 200,   
            "max_genes": 2500,  
            "max_mt": 20,       
            "min_cells_pct": 0.01
        }
    },
    "human_brain_10k": {
        "raw_h5": "data/10k_Human_Brain_MO_gemx_filtered_feature_bc_matrix.h5",
        "rna_out": "data/human_brain_10k_rna.h5ad",
        "qc": {
            "min_genes": 500,   
            "max_genes": 8000,  
            "max_mt": 5,       
            "min_cells_pct": 0.01
        }
    },
    "lymphoma_14k": {
        "raw_h5": "data/lymph_node_lymphoma_14k_filtered_feature_bc_matrix.h5",
        "rna_out": "data/lymphoma14k_rna.h5ad",
        "qc": {
            "min_genes": 200,   
            "max_genes": 6000,  
            "max_mt": 10,       
            "min_cells_pct": 0.01
        }
    }
}

src = Path(DATASETS[CURRENT_RUN]["raw_h5"])
# raw data loading
adata = sc.read_10x_h5(src)    
adata.var_names_make_unique()
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
n_cells_raw = adata.n_obs
n_genes_raw = adata.n_vars
print(f"raw: {adata.shape}")



fig, axes = plt.subplots(1, 3, figsize=(18, 6))
plt.subplots_adjust(wspace=0.4)
# Genes per cell 
sc.pl.violin(adata, "n_genes_by_counts", ax=axes[0], show=False, multi_panel=False)
axes[0].set_title("Genes per cell\n(min_genes / max_genes)", fontsize=12)
axes[0].set_ylabel("Gene counts")

# UMI counts 
sc.pl.violin(adata, "total_counts", ax=axes[1], show=False, multi_panel=False)
axes[1].set_title("UMI counts\n(Seq depth)", fontsize=12)
axes[1].set_ylabel("Total Counts")

# Mitochondrial percentage 
sc.pl.violin(adata, "pct_counts_mt", ax=axes[2], show=False, multi_panel=False)
axes[2].set_title("Mitochondrial RNA %\n(pct_mt)", fontsize=12)
axes[2].set_ylabel("% MT counts")
plt.savefig(LOG_DIR / f"{CURRENT_RUN}_qc_plots.png")
plt.close()

# standard QC + norm
sc.pp.filter_cells(adata, min_genes=DATASETS[CURRENT_RUN]["qc"]["min_genes"])
sc.pp.filter_cells(adata, max_genes=DATASETS[CURRENT_RUN]["qc"]["max_genes"])
sc.pp.filter_genes(adata, min_cells=3)
adata = adata[adata.obs["pct_counts_mt"] < DATASETS[CURRENT_RUN]["qc"]["max_mt"]].copy()

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

stats_file = LOG_DIR / f"{CURRENT_RUN}_summary.txt"
with open(stats_file, "w") as f:
    f.write(f"Dataset: {CURRENT_RUN}\n")
    f.write(f"Raw: {n_cells_raw} cells, {n_genes_raw} genes\n")
    f.write(f"Post-QC: {adata.n_obs} cells, {adata.n_vars} genes\n")
    f.write(f"Median genes/cell: {adata.obs['n_genes_by_counts'].median():.1f}\n")
    f.write(f"Output saved to: {DATASETS[CURRENT_RUN]['rna_out']}\n")


out = Path(DATASETS[CURRENT_RUN]["rna_out"])
adata.write_h5ad(out)
print(f"wrote {out}: {out.stat().st_size/1e6:.1f} MB")
print(f"[{CURRENT_RUN}] Done. Logs: {stats_file}, Plots: {LOG_DIR}/{CURRENT_RUN}_qc_plots.png")