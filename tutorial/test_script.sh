PATH_TO_SINGLECELL_SEURAT_OBJECT="data/scref_Lung_UBla/lung_combined_macs_DCs_seurat_annot_ident.rds"
PATH_TO_VISIUMHD="data/Visium_HD_Lung_ILA010/binned_outputs/square_032um"
ANNOTATION_COLUMN_IN_SINGLECELL_OBJECT="annot_ident"
NUM_CORES=1

Rscript runRCTD.R \
--sc_input $PATH_TO_SINGLECELL_SEURAT_OBJECT \
--sp_input $PATH_TO_VISIUMHD \
--annot $ANNOTATION_COLUMN_IN_SINGLECELL_OBJECT