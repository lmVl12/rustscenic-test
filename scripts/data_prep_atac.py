"""End-to-end example: fragments.tsv.gz + peaks.bed -> cells × peaks AnnData.
Workflow mirrors what you'd do on real 10x multiome output:
    1. Write tiny fragments.tsv.gz + peaks.bed to a temp dir
    2. Call rustscenic.preproc.fragments_to_matrix
    3. Verify shape, per-cell QC, and the expected counts
"""
from __future__ import annotations

import gzip
import tempfile
import time
from pathlib import Path
# Add a dictionary for easy rerun
CURRENT_RUN = "human_brain_10k"  # Change this key to switch datasets

DATASETS = {
    "pbmc_3k": {
        "fragments": "data/pbmc_unsorted_3k_atac_fragments.tsv.gz",
        "peaks": "data/pbmc_unsorted_3k_atac_peaks.bed",
        "atac_h5ad": "data/pbmc_3k_atac.h5ad",
        "rna_out": "data/pbmc_3k_rna.h5ad",
        "qc": {
            "min_cells_pct": 0.01,
            "n_top_peaks": 150000
        }
    },
    "human_brain_10k": {
        "fragments": "data/10k_Human_Brain_MO_gemx_atac_fragments.tsv.gz",
        "peaks": "data/10k_Human_Brain_MO_gemx_atac_peaks.bed",
        "atac_h5ad": "data/human_brain_10k_atac.h5ad",
        "rna_out": "data/human_brain_10k_rna.h5ad",
        "qc": {
            "min_cells_pct": 0.01,
            "n_top_peaks": 150000
        }
    },
    "lymphoma_14k": {
        "fragments": "data/lymph_node_lymphoma_14k_atac_fragments.tsv.gz",
        "peaks": "data/lymph_node_lymphoma_14k_atac_peaks.bed",
        "atac_h5ad": "data/lymphoma14k_atac.h5ad",
        "rna_out": "data/lymphoma14k_rna.h5ad",
        "qc": {
            "min_cells_pct": 0.01,
            "n_top_peaks": 150000
        }
    }
}

def main() -> int:
    import rustscenic.preproc
    import scanpy as sc

    fragments_path = Path(DATASETS[CURRENT_RUN]["fragments"])
    peaks_path     = Path(DATASETS[CURRENT_RUN]["peaks"])
    rna_path = Path(DATASETS[CURRENT_RUN]["rna_out"])
    t0 = time.perf_counter()
    adata = rustscenic.preproc.fragments_to_matrix(
        fragments_path, peaks_path
    )
    print(f"\nBuilt raw matrix in {(time.perf_counter()-t0)*1e3:.1f} ms")
    print(f"Initial shape: {adata.shape}")
# cells selected by RNA barcodes, then filter peaks by min_cells_pct, then select top variable peaks if >150k (just in case of large peak sets)
    print("Loading RNA barcodes for synchronization...")
    if rna_path.exists():
        rna_adata = sc.read_h5ad(rna_path)
        valid_barcodes = rna_adata.obs_names
        adata = adata[adata.obs_names.isin(valid_barcodes)].copy()
        print(f"Subsetted to {adata.n_obs} valid cells based on RNA matrix.")
    else:
        print("Warning: RNA file not found. Use manual thresholding.")

    min_cells = int(adata.n_obs * DATASETS[CURRENT_RUN]["qc"]["min_cells_pct"]) 
    sc.pp.filter_genes(adata, min_cells=min_cells)
    print(f"\nFinal shape:        {adata.shape}  (cells x peaks)")

    adata.write_h5ad(DATASETS[CURRENT_RUN]["atac_h5ad"])
    print(f"\nMatrix saved {DATASETS[CURRENT_RUN]['atac_h5ad']}")

    elapsed = time.perf_counter() - t0
    print(f"\nbuilt matrix in {elapsed*1e3:.1f} ms")
    print(f"cells:        {list(adata.obs_names)}")
    print(f"peaks:        {list(adata.var_names)}")
    print(f"\nper-cell QC (.obs):")
    print(adata.obs.to_string())
    print(adata.X[:5, :5].toarray())
    return 0 
    
if __name__ == "__main__":
    raise SystemExit(main())

import scanpy as sc
adata = sc.read_h5ad("data/human_brain_10k_atac.h5ad")

# 1. Скільки піків залишилось?
print(f"Кількість піків після фільтрації: {adata.n_vars}")

# 2. Чи є там взагалі щось, окрім 0?
print("Середня кількість піків на клітину:", adata.X.sum(axis=1).mean())