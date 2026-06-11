import scanpy as sc
sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()

import umap
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import pickle
import os
from multiprocessing import Pool

n_cores = int(os.environ["SLURM_JOB_CPUS_PER_NODE"])
print("Number of cores:", n_cores)

path = "/data/gent/vo/000/gvo00070/vsc43831/visium-hd/"

print("Reading anndata...")
adata = sc.read_10x_mtx(
    path + "data/pbmc3k_filtered_gene_bc_matrices/hg19/",  # the directory with the `.mtx` file
    var_names="gene_symbols",  # use gene symbols for the variable names (variables-axis index)
    cache=True,  # write a cache file for faster subsequent reading
)

print("Preprocessing data...")
adata.var_names_make_unique()  # this is unnecessary if using `var_names='gene_ids'` in `sc.read_10x_mtx`

sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

# annotate the group of mitochondrial genes as "mt"
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
)

adata = adata[adata.obs.n_genes_by_counts < 2500, :]
adata = adata[adata.obs.pct_counts_mt < 5, :].copy()

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

adata.raw = adata

adata = adata[:, adata.var.highly_variable]
sc.pp.regress_out(adata, ["total_counts", "pct_counts_mt"])
sc.pp.scale(adata, max_value=10)

sc.pp.pca(adata, n_comps = 50)
global pcs
pcs = np.copy(adata.obsm['X_pca'])

def compute_umap(params):
    n_pcs, n_neighbors, min_dist = params
    pcs_subset = pcs[:, :n_pcs]
    reducer = umap.UMAP(n_neighbors=n_neighbors, min_dist=min_dist, metric="euclidean",
                        random_state=0, n_jobs=1)
    embedding = reducer.fit_transform(pcs_subset)
    return (f"pcs_{n_pcs}_neighbors_{n_neighbors}_mindist_{min_dist}", embedding)

min_dist_list = [0, 0.1, 0.25, 0.5, 0.8, 0.99]
n_neighbors_list = [5, 10, 15, 20, 50, 100]
n_pcs_list = [15, 50]
embeddings = {}

for n_pcs in n_pcs_list:
    for n_neighbors in n_neighbors_list:
        for min_dist in min_dist_list:
            sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs = n_pcs)
            sc.tl.umap(adata, min_dist=min_dist)
            
            embeddings[f"pcs_{n_pcs}_neighbors_{n_neighbors}_mindist_{min_dist}"] = adata.obsm['X_umap']


# Save embeddings
with open(path + 'data/pbmc3k_embeddings.p', 'wb') as fp:
    pickle.dump(embeddings, fp, protocol=pickle.HIGHEST_PROTOCOL)

for n_pcs in n_pcs_list:

    # Plot the embeddings
    fig, axes = plt.subplots(len(n_neighbors_list), len(min_dist_list), figsize=(20, 16))

    for i in range(len(n_neighbors_list)):
        for j in range(len(min_dist_list)):
            
            embedding = embeddings[f"pcs_{n_pcs}_neighbors_{n_neighbors_list[i]}_mindist_{min_dist_list[j]}"]
            
            # Create pandas dataframe combining embedding and cell type annotation
            df = pd.DataFrame(embedding, columns=['X', 'Y'])

            sns.scatterplot(data=df, x='X', y='Y', ax=axes[i, j], s=0.5)
            
            # Turn off axis labels and text
            axes[i, j].get_xaxis().set_ticks([])
            axes[i, j].get_yaxis().set_ticks([])
            axes[i, j].set_xlabel('')
            axes[i, j].set_ylabel('')
            
            axes[i, j].set_title(f"n_neighbors={n_neighbors_list[i]}, min_dist={min_dist_list[j]}")
        
    plt.tight_layout()

    # Save plot
    fig.savefig(path + "data/pbmc3k_plots/umap_embeddings_{}pcs.png".format(n_pcs))

