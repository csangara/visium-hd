import os
import scanpy as sc
file_path = os.environ['VSC_DATA_VO_USER'] + '/spatial_datasets/Visium_HD_Mouse_Brain/square_008um/'

# Read anndata from 10x output
adata = sc.read_10x_h5(file_path + 'filtered_feature_bc_matrix.h5')

# Change var_names to gene ids
adata.var["gene_names"] = adata.var_names
adata.var_names = adata.var["gene_ids"]

adata.write(os.environ['VSC_DATA_VO_USER'] + '/rds/Visium_HD_Mouse_Brain_008um.h5ad')
