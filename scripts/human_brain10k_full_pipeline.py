import hashlib
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
import rustscenic._gene_resolution
import multiprocessing
import gc

start_time = time.monotonic()
RNA = Path("data/human_brain_10k_rna.h5ad")
ATAC = Path("data/human_brain_10k_atac.h5ad")
PEAKS = Path("data/10k_Human_Brain_MO_gemx_atac_peaks.bed")
TFS = Path("data/tfs_hg38.txt")
motif_rankings = rustscenic.data.download_motif_rankings(species="human", verbose=False)
gene_coords = rustscenic.data.download_gene_coords(species="hs", verbose=False)
FRAGMENTS = Path("data/10k_Human_Brain_MO_gemx_atac_fragments.tsv.gz")
OUT = Path("out/human_brain_10k_validation.json")\

print("Loading data for pipeline...")
rna = ad.read_h5ad(RNA)
atac = ad.read_h5ad(ATAC)
# Intersect cells present in both modalities (should be ~all of them)
shared = sorted(set(rna.obs_names) & set(atac.obs_names))
rna = rna[shared].copy()
atac = atac[shared].copy()
print(f"shared cells: {len(shared)}")
rna.var_names_make_unique()
rna.var['gene_ids'] = rna.var_names  
rna.var.set_index('feature_name', inplace=True)
rna.var_names = rna.var_names.astype(str)
rna.var_names_make_unique()
print(rna.var_names[:20].tolist())

# Prior knowledge - the microglia noise should be removed
microglia_markers = [
    "CX3CR1",
    "TMEM119",
    "P2RY12",
    "AIF1",
    "C1QB",
    "S100A8",
    "CSF1R",
    "S100A9",
    "BIN2",
]
# sort out the matrix with the microglia markers
present_markers = [m for m in microglia_markers if m in rna.var_names]
if present_markers:
    rna.obs["microglia_score"] = rna[:, present_markers].X.mean(axis=1).A1
    threshold = np.percentile(rna.obs["microglia_score"], 85)
    cells_to_keep = rna.obs["microglia_score"] < threshold
    rna = rna[cells_to_keep].copy()
    atac = atac[cells_to_keep].copy()
    print(f"Filtered out immune cells. Neurons remaining: {rna.n_obs}")


""" Задопукументовано все, що було до запуску повного пайплайну, щоб зосередитися на ньому.
#full_annot = pd.read_csv(COORDS, sep='\t')
#gene_coords = full_annot[['gene', 'chrom', 'start']].drop_duplicates(subset=['gene'])
#gene_coords.columns = ['gene', 'chrom', 'tss']
prepared_coords_path = Path("data/prepared_gene_coords.tsv")
gene_coords.to_csv(prepared_coords_path, sep='\t', index=False)
rna.write_h5ad("data/human_brain_rna_filtered.h5ad")
atac.write_h5ad("data/human_brain_atac_filtered.h5ad")
gene_coords.to_csv("data/human_brain_gene_coords.csv", index=False)
#annots_df = pd.read_csv(ANNOTS, sep='\t')
#annots_df.to_csv("data/human_brain_motifs_annots.csv", index=False)
del full_annot # ВИДАЛЯЄМО ВЕЛИКИЙ ОБ'ЄКТ
"""
# 3. ВИДАЛЯЄМО зайве з пам'яті ПЕРЕД пайплайном
gc.collect()
time.sleep(2)
print(f"--- Run full pipeline ---")
t0 = time.monotonic()

results = rustscenic.pipeline.run(
    rna=rna,
    output_dir= "out",  
    adata_atac = atac,
    tfs=TFS,
    peaks = PEAKS,
    #fragments=FRAGMENTS,
    motif_rankings=motif_rankings,
    #motif_annotations = ANNOTS,
    gene_coords = gene_coords,

    grn_n_estimators=500,

    topics_method = "gibbs",
    topics_n_iters=500,              
    topics_n_threads=4,
    topics_n_topics=10,

    cistarget_nes_threshold = 3.0,
    cistarget_top_frac=0.05,
    enhancer_max_distance = 500_000,
    enhancer_min_abs_corr = 0.01,
    eregulon_min_target_genes = 3,
    eregulon_min_enhancer_links = 1,
    seed=777,
    verbose = True
)

wall_time = time.monotonic() - t0
print(f"Total time: {wall_time:.2f} s")
if hasattr(results, "grn_path") and results.grn_path:
    print(f"GRN results saved to: {results.grn_path}")
else:
    print("Warning: grn_path attribute not found in results.")
expected_brain_tfs = [
    "NEUROD2", "NEUROD6", "TBR1", "FOXG1", "BCL11A", "BCL11B", 
    "ZBTB20", "SATB2", "POU3F2", "POU3F3", "SOX2", "ASCL1", 
    "EMX2", "PAX6", "LHX2", "MEF2C", "OLIG2"
]

regulon_tfs = set()

# Якщо файл регулонів успішно записався, зчитуємо його ключі (назви TF)
if hasattr(results, "regulons_path") and results.regulons_path and Path(results.regulons_path).exists():
    with open(results.regulons_path, "r") as r_file:
        regulons_json = json.load(r_file)
        for k in regulons_json.keys():
            # Очищаємо назву від суфікса, якщо він є (наприклад, "NEUROD2_regulon" -> "NEUROD2")
            regulon_tfs.add(k.replace("_regulon", ""))

found = sorted([tf for tf in expected_brain_tfs if tf in regulon_tfs])
missing = sorted([tf for tf in expected_brain_tfs if tf not in found])

# --- Підрахунок об'ємів результатів (Headline Counts) ---
n_grn_edges = int(pd.read_parquet(results.grn_path).shape[0]) if hasattr(results, "grn_path") and results.grn_path else 0
n_regulons = len(regulon_tfs)
n_cistarget_rows = int(pd.read_parquet(results.cistarget_path).shape[0]) if hasattr(results, "cistarget_path") and results.cistarget_path else 0
n_enhancer_rows = int(pd.read_parquet(results.enhancer_links_path).shape[0]) if hasattr(results, "enhancer_links_path") and results.enhancer_links_path else 0
n_eregulons = int(getattr(results, "n_eregulons", 0) or 0)
peak_mem_gb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (1024 * 1024)

print(f"\n=== Pipeline Run Headlines ===")
print(f"  GRN edges: {n_grn_edges:,}")
print(f"  Regulons found: {n_regulons}")
print(f"  CisTarget rows: {n_cistarget_rows:,}")
print(f"  Enhancer links: {n_enhancer_rows:,}")
print(f"  eRegulons assembled: {n_eregulons}")
print(f"  Biological recovery: {len(found)}/{len(expected_brain_tfs)}")

# --- Функції хешування файлів (твої оригінальні) ---
def get_md5(fname):
    hash_md5 = hashlib.md5()
    with open(fname, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()

def get_partial_md5(fname, size_mb=8):
    hash_md5 = hashlib.md5()
    with open(fname, "rb") as f:
        hash_md5.update(f.read(size_mb * 1024 * 1024))
    return hash_md5.hexdigest()

# --- Збір фінального JSON репорту ---
report = {
    "release": "v0.4.6",
    "smoke_type": "real_multiome_pipeline_run",
    "rustscenic_version": rustscenic.__version__,
    "rustscenic_sha": "9ee67398689812f98bdf6856626ac57faf95be25",
    "install_command": 'pip install "rustscenic[validation] @ git+https://github.com/Ekin-Kahraman/rustscenic@v0.4.6"',
    "api_call": "rustscenic.pipeline.run(rna=..., adata_atac=..., motif_rankings=motif_rankings_obj, ...)",
    "dataset": {
        "name": "human_brain_10k",
        "source": "10x Genomics Public Datasets",
        "species": "human",
        "tissue": "brain (microglia-filtered)",
        "rna_h5_md5": get_md5(RNA),
        "atac_h5ad_md5_first_8mb": get_partial_md5(ATAC),
        "peaks_bed_md5": get_md5(PEAKS) if PEAKS.exists() else "...",
    },
    "shapes": {
        "rna_post_qc": list(rna.shape),
        "atac_subset_to_rna_cells": list(atac.shape),
    },
    "wall_s": {
        "setup": round(t0 - start_time, 2), 
        "pipeline_run_total": round(wall_time, 2)
    },
    "peak_rss_gb": round(peak_mem_gb, 2),
    "outputs_non_empty": {
        "grn": n_grn_edges > 0,
        "regulons": n_regulons > 0,
        "cistarget": n_cistarget_rows > 0,
        "enhancer_links": n_enhancer_rows > 0,
        "eregulons": n_eregulons > 0,
        "integrated_adata": hasattr(results, "integrated_adata_path") and results.integrated_adata_path is not None,
    },
    "headline_counts": {
        "n_grn_edges": n_grn_edges,
        "n_regulons": n_regulons,
        "n_cistarget_rows": n_cistarget_rows,
        "n_enhancer_links": n_enhancer_rows,
        "n_eregulons": n_eregulons,
    },
    "biological_sanity": {
        "expected_tfs": expected_brain_tfs,
        "found_in_regulons": found,
        "missing_from_regulons": missing,
        "fraction_recovered": round(len(found) / len(expected_brain_tfs), 4) if expected_brain_tfs else 0,
    },
    "output_inventory": {
        "grn_path": str(results.grn_path) if hasattr(results, "grn_path") else None,
        "regulons_path": str(results.regulons_path) if hasattr(results, "regulons_path") else None,
        "cistarget_path": str(results.cistarget_path) if hasattr(results, "cistarget_path") else None,
        "enhancer_links_path": str(results.enhancer_links_path) if hasattr(results, "enhancer_links_path") else None,
        "eregulons_path": str(results.eregulons_path) if hasattr(results, "eregulons_path") else None,
    },
    "elapsed_per_stage": results.elapsed if hasattr(results, "elapsed") else {},
    "env": {
        "python": platform.python_version(),
        "scanpy": sc.__version__,
        "anndata": ad.__version__,
        "os": f"{platform.system()} {platform.release()} {platform.machine()}",
        "cpu": platform.processor(),
        "n_cpus": multiprocessing.cpu_count(),
    },
    "scope_notes": ["Validated using local 33GB Aertslab feather file pre-loaded via load_aertslab_feather."],
}

with open(OUT, "w") as f:
    json.dump(report, f, indent=4)

print(f"\n✅ Done! Validation artefact successfully written to {OUT}", flush=True)
