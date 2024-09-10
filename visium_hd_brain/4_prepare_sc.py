import os
import anndata as ad
import pandas as pd

# Define the path to the subset folder
subset_folder = os.environ['VSC_DATA_VO_USER'] + '/ABA_data/spatial_datasets/subset/'
output_file = os.environ['VSC_DATA_VO_USER'] + '/ABA_data/spatial_datasets/WMB-10Xv3_subset.h5ad'

# List all .h5ad files in the subset folder
subset_files = sorted([f for f in os.listdir(subset_folder) if f.endswith('.h5ad')])

# Initialize an empty list to hold the AnnData objects
adatas = []

# Read each .h5ad file and append to the list
for subset_file in subset_files:
    file_path = os.path.join(subset_folder, subset_file)
    adata = ad.read_h5ad(file_path)
    adatas.append(adata)

# Concatenate all AnnData objects into a single AnnData object
combined_adata = ad.concat(adatas, join='outer', merge='same')

# Download metadata
downsampled_csv = os.environ['VSC_DATA_VO_USER'] + '/ABA_data/spatial_datasets/cell_metadata_with_cluster_annotation_downsampled.csv'
metadata = pd.read_csv(downsampled_csv)

# Ensure that the index of the metadata dataframe is the 'cell_label' column
metadata.set_index('cell_label', inplace=True)

# Identify overlapping columns
overlapping_columns = metadata.columns.intersection(combined_adata.obs.columns)
print(f"Overlapping columns: {overlapping_columns}")

# Rename overlapping columns in the metadata dataframe
metadata.rename(columns={col: f"{col}_metadata" for col in overlapping_columns}, inplace=True)

# Align the metadata indices with the obs_names of the AnnData object
metadata = metadata.loc[combined_adata.obs_names]

# Check if all indices are aligned
assert all(metadata.index == combined_adata.obs_names)

# Add metadata columns to the AnnData object
# Join using the indices
combined_adata.obs = combined_adata.obs.join(metadata)

# Save the combined AnnData object to a file
combined_adata.write(output_file)
print(f"Combined data saved to {output_file}")