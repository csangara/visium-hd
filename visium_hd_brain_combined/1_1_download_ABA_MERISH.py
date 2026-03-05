from pathlib import Path
from abc_atlas_access.abc_atlas_cache.abc_project_cache import AbcProjectCache
import os

# I just ran this on the command line of the HPC
# Had to update manifest though
#### DOWNLOAD MERFISH CELL TYPE METADATA TABLE ####
download_base = '/data/gent/vo/000/gvo00070/vsc43831/spatial_datasets/ABA_data'
abc_cache = AbcProjectCache.from_s3_cache(download_base)

cell = abc_cache.get_metadata_dataframe(
    directory='MERFISH-C57BL6J-638850',
    file_name='cell_metadata',
    dtype={"cell_label": str}
)
cell.set_index('cell_label', inplace=True)
print("Number of cells = ", len(cell))

# Get the cluster metadata table 
cluster_details = abc_cache.get_metadata_dataframe(
    directory='WMB-taxonomy',
    file_name='cluster_to_cluster_annotation_membership_pivoted',
    keep_default_na=False
)
cluster_details.set_index('cluster_alias', inplace=True)
cluster_details.head(5)

cluster_colors = abc_cache.get_metadata_dataframe(directory='WMB-taxonomy', file_name='cluster_to_cluster_annotation_membership_color')
cluster_colors.set_index('cluster_alias', inplace=True)
cluster_colors.head(5)

## Read in section, reconstructed and CCF coordinates for all cells ##
cell = abc_cache.get_metadata_dataframe(
    directory='MERFISH-C57BL6J-638850',
    file_name='cell_metadata_with_cluster_annotation')
cell.rename(columns={'x': 'x_section',
                     'y': 'y_section',
                     'z': 'z_section'},
            inplace=True)
cell.set_index('cell_label', inplace=True)
cell.head(5)

# Reconstructed coordinates
reconstructed_coords = abc_cache.get_metadata_dataframe(
    directory='MERFISH-C57BL6J-638850-CCF',
    file_name='reconstructed_coordinates',
    dtype={"cell_label": str}
)
reconstructed_coords.rename(columns={'x': 'x_reconstructed',
                                     'y': 'y_reconstructed',
                                     'z': 'z_reconstructed'},
                            inplace=True)
reconstructed_coords.set_index('cell_label', inplace=True)
reconstructed_coords.head(5)

# CCF coordinates
ccf_coords = abc_cache.get_metadata_dataframe(
    directory='MERFISH-C57BL6J-638850-CCF',
    file_name='ccf_coordinates',
    dtype={"cell_label": str}
)
ccf_coords.rename(columns={'x': 'x_ccf',
                           'y': 'y_ccf',
                           'z': 'z_ccf'},
                  inplace=True)
ccf_coords.drop(['parcellation_index'], axis=1, inplace=True)
ccf_coords.set_index('cell_label', inplace=True)
ccf_coords.head(5)

# Join cell metadata with reconstructed coordinates
cell_joined = cell.join(reconstructed_coords, how='inner')
cell_joined = cell_joined.join(ccf_coords, how='inner')
cell_joined.head(5)

# Read in the pivoted parcellation annotation information
parcellation_annotation = abc_cache.get_metadata_dataframe(directory='Allen-CCF-2020',
                                                           file_name='parcellation_to_parcellation_term_membership_acronym')
parcellation_annotation.set_index('parcellation_index', inplace=True)
parcellation_annotation.columns = ['parcellation_%s'% x for x in  parcellation_annotation.columns]
parcellation_annotation.head(5)

parcellation_color = abc_cache.get_metadata_dataframe(directory='Allen-CCF-2020',
                                                      file_name='parcellation_to_parcellation_term_membership_color')
parcellation_color.set_index('parcellation_index', inplace=True)
parcellation_color.columns = ['parcellation_%s'% x for x in  parcellation_color.columns]
parcellation_color.head(5)

cell_joined = cell_joined.join(parcellation_annotation, on='parcellation_index')
cell_joined = cell_joined.join(parcellation_color, on='parcellation_index')
cell_joined.head(5)

# Can actually skip all of this by downloading the data directly
final_cell_joined = abc_cache.get_metadata_dataframe(
    directory='MERFISH-C57BL6J-638850-CCF',
    file_name='cell_metadata_with_parcellation_annotation',
    dtype={"cell_label": str}
    )

# Zhuang dataset
datasets = ['Zhuang-ABCA-1', 'Zhuang-ABCA-2']
for d in datasets :

    cell = abc_cache.get_metadata_dataframe(
        directory=d,
        file_name='cell_metadata_with_cluster_annotation',
        dtype={"cell_label": str}
    )
    ccf_coordinates = abc_cache.get_metadata_dataframe(directory=f"{d}-CCF", file_name='ccf_coordinates')