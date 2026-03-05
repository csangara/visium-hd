## VisiumHD Brain FFPE analysis

In addition to the VisiumHD analysis, this repository also contains scripts to download the single-cell reference of the brain data, from the [Allen Brain Cell Atlas - Data Portal](https://alleninstitute.github.io/abc_atlas_access/intro.html). It would be useful to follow the [Getting started](https://alleninstitute.github.io/abc_atlas_access/notebooks/getting_started.html) tutorial to install the package used to download the data.

Scripts:
* 1_1: Download whole mouse brain (v3) single-cell data
* 1_2: Download VisiumHD data from 10x Genomics
* 2: Downsample single-cell data to about 86,000 cells
* 3_1: Create AnnData object of downsampled single-cell data
* 3_2: Create AnnData object of VisiumHD data
* 3_3: Create Seurat object of VisiumHD and single-cell data
* 4_1: Explore VisiumHD data
* 4_2: Use Stardist to segment nuclei
* 4_3: Load in the nuclei segmentation mask and visualize and calculate some metrics.
* 4_4: Annotate brain divisions in the VisiumHD data (using data from QUINT workflow)
* 5: Run RCTD, cell2location, and Scrublet. Note that cell2location does not work.
* 6: Visualize RCTD and scrublet results.
* 6_3: Compare doublet prediction of RCTD and scrublet with segmented nuclei
* ex: Compare predictions when using actual spot coordinates vs dummy coordinates, also compare predictions between doublet and full mode

Predicted proportions are compressed into `props.tar.gz` and plots are compressed into `plots.tar.gz`.
