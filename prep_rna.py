"""Preprocess PBMC-10k for GRN inference."""
from pathlib import Path
import scanpy as sc

# Add a dictionary of datasets (to be changed every time)

CURRENT_RUN = "pbmc_3k"

DATASETS = {
    "pbmc_3k": {
        "raw_h5": "data/pbmc_unsorted_3k_filtered_feature_bc_matrix.h5",
        "rna_out": "data/pbmc_3k_rna.h5ad"
    },
    "brain_10k": {
        "raw_h5": "data/10k_Human_Brain_MO_gemx_raw_feature_bc_matrix.h5",
        "rna_out": "data/human_brain_10k_rna.h5ad"
    }
}

src = Path(DATASETS[CURRENT_RUN]["raw_h5"])

adata = sc.read_10x_h5(src)
adata.var_names_make_unique()
print(f"raw: {adata.shape}")

# standard QC + norm
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
adata = adata[adata.obs["pct_counts_mt"] < 20].copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
print(f"post-QC: {adata.shape}")

out = Path(DATASETS[CURRENT_RUN]["rna_out"])
adata.write_h5ad(out)
print(f"wrote {out}: {out.stat().st_size/1e6:.1f} MB")
