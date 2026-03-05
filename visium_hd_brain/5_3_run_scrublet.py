import scanpy as sc
adata = sc.read_h5ad("/data/gent/vo/000/gvo00070/vsc43831/spatial_datasets/Visium_HD_MouseBrain/Visium_HD_MouseBrain_008um.h5ad")
sc.external.pp.scrublet(adata)
# Save obs dataframe to csv
adata.obs.to_csv("/data/gent/vo/000/gvo00070/vsc43831/visium-hd/visium_hd_brain/Visium_HD_MouseBrain_008um/scrublet_Visium_HD_MouseBrain_008um.csv")