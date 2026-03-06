import spatialdata as sd
import anndata as ad
import spatialdata_io
import harpy as sp
import scanpy as sc
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("marker_gene_file", type=str, help="Path to the marker gene file (CSV format).")
parser.add_argument("-i", "--iter", help="Flag to use score_genes_iter instead of score_genes",
                    action="store_true")
parser.add_argument("-f", "--first_run", help="Flag to indicate that it's the first run",
                    action="store_true")                  
args = parser.parse_args()

sc.settings.figdir = "./plots/"

base_path = "/data/gent/gvo000/gvo00070/vsc43831/"
marker_gene_path = f"{base_path}visium-hd/liver_marker_genes/"

# Check that file exists
marker_gene_file = args.marker_gene_file
marker_gene_ext = marker_gene_file.split("_")[-1].split(".")[0]
output_layer_name = f"square_008um_atlas_{marker_gene_ext}_celltype_scores"
output_layer_name = f'{output_layer_name}{"_iter" if args.iter else ""}'
score_genes_func = sp.tb.score_genes_iter if args.iter else sp.tb.score_genes

print(output_layer_name)
print(score_genes_func)

import os
if not os.path.exists(marker_gene_path+marker_gene_file):
    raise FileNotFoundError(f"{marker_gene_file} does not exist.")

if args.first_run:
    sdata = spatialdata_io.visium_hd(base_path + 'spatial_datasets/Visium_HD_Liver_CAW009',
                                    fullres_image_file="microscope_image/CAW009_HE_hires.tif",
                                    annotate_table_by_labels=True)
    
    # We will change some AnnData properties to make it compatible with sparrow functions
    # See https://github.com/saeyslab/napari-sparrow/blob/main/src/sparrow/io/_visium_hd.py
    from harpy.utils._keys import _INSTANCE_KEY, _REGION_KEY
    from spatialdata.models import TableModel
    from spatialdata_io._constants._constants import VisiumHDKeys

    for tables_layer in [*sdata.tables]:
        adata = sdata[tables_layer]
        adata.var_names_make_unique()
        adata.X = adata.X.tocsc()

        adata.obs.rename(
            columns={VisiumHDKeys.REGION_KEY: _REGION_KEY, VisiumHDKeys.INSTANCE_KEY: _INSTANCE_KEY}, inplace=True
        )
        adata.uns.pop(TableModel.ATTRS_KEY)
        adata = TableModel.parse(
            adata,
            region_key=_REGION_KEY,
            region=adata.obs[_REGION_KEY].cat.categories.to_list(),
            instance_key=_INSTANCE_KEY,
        )

        del sdata[tables_layer]
        sdata[tables_layer] = adata

    for shapes_layer in [*sdata.shapes]:
        sdata[shapes_layer].index.name = _INSTANCE_KEY

    print(sdata)
    sdata.write(base_path + 'spatial_datasets/Visium_HD_Liver_CAW009/Visium_HD_Liver_CAW009_zarr')

    adata = sdata.tables['square_008um']
    print(adata)
    
    adata.var["mt"] = adata.var_names.str.startswith("Mt-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

    sc.pp.filter_cells(adata, min_counts=50)
    # sc.pp.filter_cells(adata, max_counts=4000)
    adata = adata[adata.obs["pct_counts_mt"] < 20].copy()
    print(f"#cells after MT filter: {adata.n_obs}")
    sc.pp.filter_genes(adata, min_cells=10)

    sc.pp.normalize_total(adata, inplace=True)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=2000)

    sc.pp.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    sc.tl.leiden(
        adata,
        key_added="clusters",
        flavor="igraph",
        directed=False,
        n_iterations=2
    )

    # Plot some stuff
    sc.pl.umap(adata, color=["total_counts", "n_genes_by_counts", "clusters"], wspace=0.4, save='_8um_metadata.png')
    sc.pl.umap(adata, color=['Glul', 'Hal'], wspace=0.4, save='_8um_zonation_markers.png')
    
    sdata.tables["square_008um_pp"] = adata
    sdata.write_element("square_008um_pp")
    
    # scale the data, but do not zero-center to keep sparsity
    sc.pp.scale(adata, zero_center=False)

    sdata.tables["square_008um_pp_scaled"] = adata
    sdata.write_element("square_008um_pp_scaled")

else:
    sdata = sd.read_zarr(base_path + 'spatial_datasets/Visium_HD_Liver_CAW009/Visium_HD_Liver_CAW009_zarr')

print(sdata)

print(f"Running {score_genes_func.__name__}...")
score_genes_func(sdata, labels_layer="CAW009_square_008um_labels",
                table_layer="square_008um_pp_scaled",
                output_layer=output_layer_name,
                path_marker_genes=f"{marker_gene_path}{marker_gene_file}",
                celltype_column='annotation')


adata = ad.read_zarr(f"{base_path}spatial_datasets/Visium_HD_Liver_CAW009/Visium_HD_Liver_CAW009_zarr/tables/{output_layer_name}")
adata.obs['annotation'].value_counts(normalize=True).to_csv(f'score_genes/{output_layer_name}_mean_props.csv')
adata.obs['annotation'].to_csv(f'score_genes/{output_layer_name}_props.csv')