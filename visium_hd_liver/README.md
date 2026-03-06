## VisiumHD Liver SCA002 analysis

This folder mostly contains scripts of me trying out things that didn't make it to the chapter. It also led to the deeper sequencing of the data (CAW009), as I found that I could not really get out that much information from this dataset.

Scripts:
* 1_1 - 1_4: Data exploration and data wrangling
* 1_5: Since this dataset was processed with SpaceRanger v3, there are no 32µm bins being provided, so I had to bin them myself.
* 1_6: Scripts for nuclei segmentation and binning, although this takes a lot of memory. 1_6_1 was done by Karen. 1_6_2 also creates a component of Figure 5.1b.
* 2_1: The various ways I tried to preprocess the 8µm bins. This starts from simply applying the single-cell analysis to the 8µm bins (2_1_1), but the resulting UMAP was just a blob that was highly influenced by the different zones. So, 2_1_2 - 2_1_7 is simply me trying different things, which did not end up working. First, I tried to subcluster the Leiden clusters, e.g., subclustering the central zone to see if we could get cell types from it. Then, I tried to use only a subset of genes (2_1_3 - 2_1_5). Finally, I used the 8µm bins that were aggregated based on the nuclei segmentation (2_1_6), and also tried using a subset of genes on that (2_1_7).
* 2_2: Preprocessing the 2µm bins. Since there were over 7 million bins, I did not go through with the UMAP creation. Instead, I also binned these onto the nuclei segmentation mask (2_2_2).
* 2_3: I also tried using bin2cell to segment and aggregate the cells, which gave an even weirder result.
* 3: Running RCTD and scrublet on the data.
* 4: Visualizing RCTD and scrublet results
* 5: An attempt to use a gene scoring approach rather than a deconvolution approach to annotate the bins. 5_1 creates the marker gene files, while 5_2 contains scripts to perform the gene scoring on 8µm bins.
* 6: Creates Figure 5.2.
* ex: Different ways to visualize the data using a scatterpie plot and a scattertriangle plot. I ended up not using them as I preferred the scatterbar plots.

Predicted proportions are compressed into `props.tar.gz` and the mean proportions of the gene scoring is compressed into `score_genes.tar.gz`.
