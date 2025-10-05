# EC subcluster analysis full pipeline
# For main fig.8
# Author: Tianyue Ming, @CIBR, @PKU


# ----------------------------------------------
# 0. Import packages
# ----------------------------------------------
import scanpy as sc
import scvelo as scv
import pandas as pd
import anndata as ad
import numpy as np
import matplotlib.pyplot as plt
import scanpy.external as sce
import meld
import gseapy as gp
import scprep
from pathlib import Path
import os

# ----------------------------------------------
# 1. Global settings
# ----------------------------------------------
# File paths
project_dir = Path.cwd()
data_dir = project_dir / "data"
results_dir = project_dir / "results"
figures_dir = project_dir / "figures"

# Create folders if not exist
results_dir.mkdir(parents=True, exist_ok=True)
figures_dir.mkdir(parents=True, exist_ok=True)

# Input files
annot_file = data_dir / "GRCm38_annotion_M23.csv"
avm1_path = data_dir / "AVM1_filtered_feature_bc_matrix"
avm2_path = data_dir / "AVM2_filtered_feature_bc_matrix"
ctrl_path = data_dir / "Ctrl_filtered_feature_bc_matrix"

# Output intermediate files
raw_merged_h5ad = results_dir / "240328_minibend_AVM_WT_merge_rawCounts.h5ad"
filtered_h5ad = results_dir / "240328_minibend_AVM_WT_merge_filterMT_log1pCounts_doublecells_delete.h5ad"
processed_h5ad = results_dir / "240330_minibend_AVM_WT_merge_processed.h5ad"
final_h5ad = results_dir / "240402_Braf_AVM_EC_cluster_reanalysis.h5ad"

# Loom files for velocity
avm_loom_file = data_dir / "20240331_AVM_1_2_combined.loom"
ctrl_loom_file = data_dir / "XJABC170_2.loom"

# Markers for EC subclusters
ec_markers = ['Pecam1', 'Tie1', 'Myh11', 'Acta2', 'Tagln', 'Des', 'Pdgfrb', 'Col1a1', 'Col1a2', 'Srarp', 'Ackr1', 'Bmx']

# GSEA settings
gsea_params = dict(
    gene_sets='GO_Biological_Process_2021',
    permutation_type='phenotype',
    permutation_num=500,
    threads=10,
    method='signal_to_noise',
    ascending=False,
    seed=42,
    min_size=50,
    cutoff=0.05,
)

# Set scanpy figure output
sc.settings.figdir = str(figures_dir)

# ----------------------------------------------
# 2. Preprocess raw 10X data
# ----------------------------------------------
print("Reading 10X data...")
AVM_1 = sc.read_10x_mtx(avm1_path, var_names='gene_symbols')
AVM_2 = sc.read_10x_mtx(avm2_path, var_names='gene_symbols')
Ctrl = sc.read_10x_mtx(ctrl_path, var_names='gene_symbols')

# Add metadata
for adata, name in zip([AVM_1, AVM_2, Ctrl], ['AVM_1', 'AVM_2', 'Ctrl']):
    adata.obs['sample'] = name
    adata.obs['cell_id'] = [f"{i}_{name}" for i in adata.obs.index]
    adata.obs_names = adata.obs['cell_id']
    adata.obs['condition'] = 'AVM' if 'AVM' in name else 'Ctrl'

# Merge datasets
adata = ad.concat([AVM_1, AVM_2, Ctrl], join='inner', index_unique=None)

# Filter to protein-coding genes
annot = pd.read_csv(annot_file)
protein_genes = annot[annot['gene_type'] == 'protein_coding']
gene_map = dict(zip(protein_genes['gene_id'], protein_genes['gene_name']))

adata = adata[:, adata.var_names.isin(protein_genes['gene_id'])]
adata.var['gene_name'] = adata.var.index.map(gene_map)
adata.var_names = adata.var['gene_name']
adata.var_names_make_unique()

# QC filtering
print("QC filtering...")
sc.pp.filter_cells(adata, min_genes=300)
sc.pp.filter_genes(adata, min_cells=50)
adata.var['mt'] = adata.var_names.str.startswith('mt-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)

# Remove doublets
import scrublet as scr
print(f"Original cells: {adata.n_obs}")
sc.external.pp.scrublet(adata, random_state=112)
adata = adata[adata.obs['predicted_doublet'] == False].copy()
print(f"After scrublet: {adata.n_obs}")

# Further filter
adata = adata[adata.obs.pct_counts_mt < 15, :]
adata = adata[adata.obs.n_genes_by_counts < 4000, :]

# Normalize
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata

# Save processed raw file
adata.write(filtered_h5ad, compression='gzip')

# ----------------------------------------------
# 3. Cell type classification
# ----------------------------------------------
print("Clustering and annotating major cell types...")
adata = sc.read_h5ad(filtered_h5ad)

sc.pp.highly_variable_genes(adata, flavor='cell_ranger', n_top_genes=10000)
adata = adata[:, adata.var['highly_variable']]
sc.tl.pca(adata, svd_solver='arpack')
sce.pp.harmony_integrate(adata, key='sample')
sc.pp.neighbors(adata, n_neighbors=50, n_pcs=50, use_rep='X_pca_harmony')
sc.tl.umap(adata)
sc.tl.louvain(adata, resolution=1.0, key_added='louvain_1.0')

# Annotate major cell types
cluster_map_large = {
    "0": "MG", "1": "DC", "2": "OL", "3": "MG", "4": "MG", "5": "OL",
    "6": "EC", "7": "EC", "8": "SMC", "9": "EC", "10": "DC", "11": "T",
    "12": "NK", "13": "MP", "14": "EC", "15": "SMC", "16": "OL", "17": "AST",
    "18": "B", "19": "PC", "20": "B"
}
adata.obs['cell_type_large'] = adata.obs['louvain_1.0'].map(cluster_map_large).astype('category')

# Save processed
adata.write(processed_h5ad, compression='gzip')

# ----------------------------------------------
# 4. EC subcluster analysis
# ----------------------------------------------
print("Subclustering ECs...")
rawseq = sc.read_h5ad(filtered_h5ad)
rnaseq = sc.read_h5ad(processed_h5ad)
VS_data = rawseq[rnaseq.obs["cell_type_large"] == 'EC']

VS_data.uns['log1p']["base"] = None
sc.pp.highly_variable_genes(VS_data, n_top_genes=10000, flavor='seurat')
VS_data = VS_data[:, VS_data.var.highly_variable]
sc.tl.pca(VS_data)
sce.pp.harmony_integrate(VS_data, key='sample')
sc.pp.neighbors(VS_data, n_neighbors=50, n_pcs=50, use_rep='X_pca_harmony')
sc.tl.louvain(VS_data, resolution=1)
sc.tl.umap(VS_data)

# Rename subclusters
cluster_map = {
    "0": "cEC", "1": "aEC", "2": "vEC",
    "3": "cEC", "4": "cEC", "5": "cEC",
    "6": "aEC", "7": "aEC"
}
VS_data.obs['cell_type_vs'] = VS_data.obs['louvain'].map(cluster_map).astype('category')
VS_data.obs['condition'] = VS_data.obs['sample'].str.split('_').str[0]

# Save
VS_data.write(final_h5ad, compression='gzip')

# 5. Plot UMAP and markers
# ----------------------------------------------
sc.pl.umap(VS_data, color=['cell_type_vs'], save="_ECs_subcluster.pdf")
sc.pl.umap(VS_data, color=['sample'], save="_ECs_sample.pdf")
sc.pl.umap(VS_data, color=['condition'], save="_ECs_condition.pdf")
sc.pl.umap(VS_data, color=['Pecam1','Slc7a5','Ackr1','Bmx'], color_map='RdPu', save="_ECs_markers.pdf")

# ----------------------------------------------
# 6. Cell type proportion plot
# ----------------------------------------------
allcounts = pd.crosstab(VS_data.obs['cell_type_vs'], VS_data.obs['condition'])
allcounts = (((allcounts / allcounts.sum(axis=0)).T) / ((allcounts / allcounts.sum(axis=0)).sum(axis=1)).T).T * 100
allcounts = allcounts.sort_values(by=["AVM"])

plt.figure(dpi=300, figsize=(15,8))
ax = allcounts.plot.bar(stacked=True, color=['DarkRed', 'grey'], width=0.9, ylim=(0,100))
ax.legend(loc='upper left', bbox_to_anchor=(1,1))
plt.axhline(y=50, color='k', linestyle='--')
plt.ylabel('Cell Proportion (%)')
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.savefig(figures_dir / "cell_type_proportion.pdf")

# ----------------------------------------------
# 7. MELD analysis
# ----------------------------------------------
sample_density = meld.MELD().fit_transform(VS_data, VS_data.obs['condition'])
sample_likelihood = meld.utils.normalize_densities(sample_density)

for cond in VS_data.obs['condition'].unique():
    VS_data.obs[f"{cond}_MELD"] = list(sample_likelihood[cond])

metadata = VS_data.obs[['cell_type_vs', 'AVM_MELD', 'condition']]
means = metadata[metadata.condition == 'AVM'].groupby('cell_type_vs')['AVM_MELD'].mean()
ctrlmeans = metadata[metadata.condition == 'Ctrl'].groupby('cell_type_vs')['AVM_MELD'].mean()
maxdifflist = (means - ctrlmeans).sort_values().index.tolist()

changeorder = {maxdifflist[i]: sorted(maxdifflist)[i] for i in range(len(maxdifflist))}
metadata['cell_type_vs'] = [changeorder[i] for i in metadata['cell_type_vs']]

fig, ax = plt.subplots(figsize=(8,5), dpi=300)
scprep.plot.jitter(metadata['cell_type_vs'], metadata['AVM_MELD'], c=metadata['condition'], legend=True, plot_means=False, ylabel='Mean AVM Likelihood', ax=ax)
ax.scatter(means.index, means, color='red', edgecolor='k', s=80)
ax.scatter(ctrlmeans.index, ctrlmeans, color='grey', edgecolor='k', s=80)
ax.set_xticklabels(maxdifflist, rotation=90)
ax.set_ylim(0, 1)
plt.tight_layout()
plt.savefig(figures_dir / "VS_meld_analysis.svg")

# ----------------------------------------------
# 8. Velocity analysis
# ----------------------------------------------
def run_velocity(adata_file, loom_file, save_name):
    adata = scv.read(adata_file, cache=True)
    ldata = scv.read(loom_file, cache=True)
    adata.obs.index = [bc[:17] for bc in adata.obs.index]
    ldata.obs.index = [bc.split(':')[1][:17] for bc in ldata.obs.index]
    scv.utils.clean_obs_names(adata)
    scv.utils.clean_obs_names(ldata)
    adata = scv.utils.merge(adata, ldata)

    scv.pp.filter_genes(adata, min_shared_counts=20)
    scv.pp.normalize_per_cell(adata)
    scv.pp.filter_genes_dispersion(adata, n_top_genes=2000)
    scv.pp.moments(adata)
    scv.tl.velocity(adata)
    scv.tl.velocity_graph(adata)
    scv.pl.velocity_embedding_stream(adata, basis='umap', save=save_name)

run_velocity(rawseq_file, avm_loom_file, "_AVM_velocity_stream.pdf")
run_velocity(rawseq_file, ctrl_loom_file, "_Ctrl_velocity_stream.pdf")
# ----------------------------------------------
# 9. DEG analysis: AVM vs Ctrl in EC subclusters
# ----------------------------------------------
print("Performing differential expression analysis (DEG)...")

# Function to perform DEG for one subcluster
def deg_analysis(subset_name, adata_all, groupby='condition'):
    subset = adata_all[adata_all.obs['cell_type_vs'] == subset_name].copy()
    subset.uns['log1p']["base"] = None
    
    sc.tl.rank_genes_groups(
        subset,
        groupby=groupby,
        method='wilcoxon',
        key_added=f'deg_{subset_name}'
    )
    
    # Save DEG table
    result = subset.uns[f'deg_{subset_name}']
    groups = result['names'].dtype.names
    degs = pd.DataFrame({
        group + '_' + key: result[key][group]
        for group in groups for key in ['names', 'scores', 'pvals', 'pvals_adj', 'logfoldchanges']
    })
    degs.to_csv(results_dir / f"DEG_{subset_name}_AVM_vs_Ctrl.csv")

# ----------------------------------------------
# 10. Save h5ad
# ----------------------------------------------
VS_data.write(results_dir / "Braf_AVM_EC_cluster_reanalysis.h5ad", compression='gzip')---


