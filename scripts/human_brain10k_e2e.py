"""End-to-end rustscenic on paired 10x Multiome: grn → aucell (RNA) +
topics (ATAC). Does regulon activity from grn correlate with ATAC topic
assignment (should — same cells, same underlying cell-type structure)?
"""
import time
import json
import platform
from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
from sklearn.metrics import adjusted_rand_score
import resource

import rustscenic, rustscenic.grn, rustscenic.aucell, rustscenic.topics
# Add a dictionary for easy rerun 
CURRENT_RUN = "human_brain_10k" 

DATASETS = {
    "pbmc_3k": {
        "rna": "data/pbmc_3k_rna.h5ad",
        "atac": "data/pbmc_3k_atac.h5ad",
        "tfs": "data/tfs_hg38.txt",
        "out": "data/pbmc_3k_results.json"
    },
    "human_brain_10k": {
        "rna": "data/human_brain_10k_rna.h5ad",
        "atac": "data/human_brain_10k_atac.h5ad",
        "tfs": "data/tfs_hg38.txt",
        "out": "out/human_brain_10k_results.json"
    }
}

# Works independently once the CURRENT RUN is predefined
RNA  = Path(DATASETS[CURRENT_RUN]["rna"])
ATAC = Path(DATASETS[CURRENT_RUN]["atac"])
TFS  = Path(DATASETS[CURRENT_RUN]["tfs"])
OUT  = Path(DATASETS[CURRENT_RUN]["out"])


rna = ad.read_h5ad(RNA)
atac = ad.read_h5ad(ATAC)
print(f"RNA  {rna.shape}   ATAC  {atac.shape}")

# Intersect cells present in both modalities (should be ~all of them)
shared = sorted(set(rna.obs_names) & set(atac.obs_names))
rna = rna[shared].copy()
atac = atac[shared].copy()
print(f"shared cells: {len(shared)}")

# Prior knowledge - the microglia noise should be removed
microglia_markers = ['CX3CR1', 'TMEM119', 'P2RY12', 'AIF1','C1QB', 'S100A8', 'CSF1R', 'S100A9', 'BIN2']
# sort out the matrix with the microglia markers
present_markers = [g for g in microglia_markers if g in rna.var_names]

if present_markers:
    # count the mean value across the cells
    rna.obs['microglia_score'] = rna[:, present_markers].X.mean(axis=1).A1
    
    # 85% of most relevant 
    threshold = np.percentile(rna.obs['microglia_score'], 85) 
    
    # Keep neuronal markers 
    cells_to_keep = rna.obs['microglia_score'] < threshold
    
    rna = rna[cells_to_keep].copy()
    atac = atac[cells_to_keep].copy()
    print(f"Filtered out immune cells. Neurons remaining: {rna.n_obs}")



# --- stage 1: grn on RNA ---
elapsed_per_stage = {"preproc": 0.0}
tfs = [t for t in rustscenic.grn.load_tfs(TFS) if t in set(rna.var_names)]
print(f"\n--- grn (RNA) ---  cells={rna.n_obs}  genes={rna.n_vars}  tfs={len(tfs)}")
t0 = time.monotonic()
grn_df = rustscenic.grn.infer(rna, tfs, seed=777, n_estimators=300, early_stop_window=25)
t_grn = time.monotonic() - t0
elapsed_per_stage["grn"] = t_grn
print(f"  wall: {t_grn:.1f}s  edges: {len(grn_df)}")

# --- stage 2: aucell on RNA ---
# Build regulons from grn (top-50 targets per TF), filter to ≥10 genes
regs = []
for tf, grp in grn_df.groupby("TF"):
    top_targets = grp.nlargest(50, "importance")["target"].tolist()
    if len(top_targets) >= 10:
        regs.append((f"{tf}_regulon", top_targets))
print(f"\n--- aucell ---  regulons={len(regs)}")
t0 = time.monotonic()
auc = rustscenic.aucell.score(rna, regs, top_frac=0.05)
t_aucell = time.monotonic() - t0
elapsed_per_stage["aucell"] = t_aucell
print(f"  wall: {t_aucell:.1f}s  shape {auc.shape}")

# --- stage 3: topics on ATAC ---
K = 10  # set to 20, as on 10 it only identifies 6 clusters
print(f"\n--- topics (ATAC) ---  K={K}")
t0 = time.monotonic()
tres = rustscenic.topics.fit(
    atac, n_topics=K, n_passes=20, batch_size=256, seed=777,
    alpha=1.0/K, eta=1.0/K,
)
t_topics = time.monotonic() - t0
elapsed_per_stage["topics"] = t_topics
topic_assign = np.asarray([int(s.replace("Topic_", "")) for s in tres.cell_assignment().values])
print(f"  wall: {t_topics:.1f}s  unique top-1 topic: {len(set(topic_assign))}")

# --- cluster cells from ATAC for ground-truth cell type proxy ---
atac_norm = atac.copy()
atac_norm.X = atac_norm.X.astype(np.float32)
sc.pp.normalize_total(atac_norm); sc.pp.log1p(atac_norm)
sc.pp.highly_variable_genes(atac_norm, n_top_genes=5000)
sc.tl.pca(atac_norm, n_comps=30, mask_var="highly_variable")
sc.pp.neighbors(atac_norm, n_neighbors=15)
sc.tl.leiden(atac_norm, resolution=0.6, flavor="igraph", n_iterations=2, directed=False) # set res to 0.5 to get more clusters
cluster = atac_norm.obs["leiden"].astype(str).values
print(f"  atac leiden clusters: {len(set(cluster))}")

# --- cross-modal check: do regulon activities discriminate the same cell clusters
# that ATAC topics do?
auc_vals = auc.values  # cells x regulons
# For each cluster, mean regulon activity; find regulons differentiating clusters
from scipy.stats import f_oneway
top_discriminative_reg = []
for r_idx, reg_name in enumerate(auc.columns):
    groups = [auc_vals[cluster == c, r_idx] for c in np.unique(cluster)]
    try:
        F, p = f_oneway(*[g for g in groups if len(g) > 1])
        top_discriminative_reg.append((reg_name, F, p))
    except Exception:
        pass
top_discriminative_reg.sort(key=lambda x: -x[1])
print(f"\ntop-10 cluster-discriminative regulons (by ANOVA F):")
for name, F, p in top_discriminative_reg[:10]:
    print(f"  {name:25s}  F={F:>7.1f}  p={p:.2e}")

# --- total pipeline timing ---
t_total = t_grn + t_aucell + t_topics
print(f"\n=== TOTAL pipeline wall-clock: {t_total:.1f}s ({t_total/60:.1f} min) ===")
print(f"  grn:    {t_grn:>6.1f}s  ({100*t_grn/t_total:.0f}%)")
print(f"  aucell: {t_aucell:>6.1f}s  ({100*t_aucell/t_total:.0f}%)")
print(f"  topics: {t_topics:>6.1f}s  ({100*t_topics/t_total:.0f}%)")
print(f"\nARI of grn-based cell-type clustering (via top regulon activity) vs ATAC leiden:")
# assign cells by top regulon
cell_by_topreg = auc.idxmax(axis=1).values
mapped = np.unique(cell_by_topreg, return_inverse=True)[1]
print(f"  {adjusted_rand_score(cluster, mapped):.4f}  ({len(set(cell_by_topreg))} unique top regulons)")

peak_mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
report = {
    "test": "https://www.10xgenomics.com/datasets/multiome-gemx-10k-human-brain",
    "rustscenic_version": rustscenic.__version__,
    "api_call": "rustscenic.grn.infer(n_estimators=300) + rustscenic.aucell.score + rustscenic.topics.fit(K=10)",
    "dataset": {
        "name": CURRENT_RUN
    },
    "shapes": {
        "rna_post_qc": list(rna.shape),
        "atac_subset_to_rna_cells": list(atac.shape)
    },
    "metrics": {
        "ari": float(adjusted_rand_score(cluster, mapped)),
        "peak_rss_gb": float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024 / 1024),
        "n_regulons": len(regs)
    },
    "elapsed_per_stage": {
        "preproc": 0.0,
        "topics": float(t_topics),
        "grn": float(t_grn),
        "aucell": float(t_aucell),
        "cistarget": None,
        "enhancer": None,
        "eregulons": None
    },
    "env": {
        "python": platform.python_version(),
        "scanpy": sc.__version__,
        "anndata": ad.__version__,
        "os": platform.platform(),
        "cpu": platform.processor()
    },
    "headline_results": [
        {"name": name, "F": float(F), "p": float(p)} 
        for name, F, p in top_discriminative_reg[:10]
    ]
}

with open(OUT, "w") as f:
    json.dump(report, f, indent=4)

print(f"\n[DONE] JSON Report saved to: {OUT}")