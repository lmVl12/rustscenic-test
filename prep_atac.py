"""End-to-end example: fragments.tsv.gz + peaks.bed -> cells × peaks AnnData.

Demonstrates rustscenic.preproc on synthetic data you can verify by eye.
No external downloads, no paths to configure — runs standalone.

Workflow mirrors what you'd do on real 10x multiome output:

    1. Write tiny fragments.tsv.gz + peaks.bed to a temp dir
    2. Call rustscenic.preproc.fragments_to_matrix
    3. Verify shape, per-cell QC, and the expected counts

Runtime: <1 second.
"""
from __future__ import annotations

import gzip
import tempfile
import time
from pathlib import Path

CURRENT_RUN = "human_brain_10k"

DATASETS = {
    "pbmc_3k": {
        "fragments": "data/pbmc_unsorted_3k_atac_fragments.tsv.gz",
        "peaks": "data/pbmc_unsorted_3k_atac_peaks.bed",
        "atac_h5ad": "data/pbmc_3k_atac.h5ad"
    },
    "human_brain_10k": {
        "fragments": "data/10k_Human_Brain_MO_gemx_atac_fragments.tsv.gz",
        "peaks": "data/10k_Human_Brain_MO_gemx_atac_peaks.bed",
        "atac_h5ad": "data/human_brain_10k_atac.h5ad"

    }
}

def main() -> int:
    import rustscenic.preproc

    fragments_path = Path(DATASETS[CURRENT_RUN]["fragments"])
    peaks_path     = Path(DATASETS[CURRENT_RUN]["peaks"])

    t0 = time.perf_counter()
    adata = rustscenic.preproc.fragments_to_matrix(
        fragments_path, peaks_path
    )
    elapsed = time.perf_counter() - t0
    print(f"\nbuilt matrix in {elapsed*1e3:.1f} ms")

    print(f"\nshape:        {adata.shape}  (cells x peaks)")
    print(f"cells:        {list(adata.obs_names)}")
    print(f"peaks:        {list(adata.var_names)}")
    print(f"\nper-cell QC (.obs):")
    print(adata.obs.to_string())
    print(f"\ncount matrix (dense view):")
    print(adata.obs.head().to_string())
    adata.write_h5ad(DATASETS[CURRENT_RUN]["atac_h5ad"])
    print("\nMatrix saved")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
