# Deconvolution of Visium HD data

We will use the tool RCTD -- now implemented as part of the spacexr package -- which is also what was done in the Visium HD announcement paper by Oliveria et al. This is mainly because it has been shown by several benchmarks to be among the top performers, and it is also relatively scalable. Another advantage is that it works well out of the box without any need for parameter tuning, and it is also (subjectively) an easy tool to run.

We highly recommend that users first follow the vignette of RCTD in the GitHub (https://github.com/dmcable/spacexr) to make sure that they understand the various functions and how the tool works. This tutorial will mostly focus on how to deploy the tool on Visium HD data.

Since Visium HD data has way more observations than normal Visium, **it is highly recommended (or even essential) that the tool is run on a computing cluster such as a HPC**. Running RCTD on 8µm bins can take up to a day!

## Requirements

### Packages

#### RCTD
Although RCTD is now available on Bioconductor, we will perform the installation from GitHub because that allows us to install the sped-up version that was implemented by the 10x Genomics team (dmcable/spacexr\#206)[https://github.com/dmcable/spacexr/pull/206]. Note that the speed up is only available for the `doublet_mode` which predicts at maximum 2 cell types, so this could be insufficient if your tissue is very densely packed together.

Install the sped-up version of RCTD on your computing cluster as follows: 

```
if(!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools") 
}
options(timeout = 600000000) # set this to avoid timeout error
devtools::install_github("jpromeror/spacexr", ref = "HD", build_vignettes = FALSE)
```

Note that for some clusters, you may not have installation rights and will have to ask the IT staff to install it for you.

#### Seurat (v5)

This is not strictly necessary since RCTD uses its own data constructors, but using Seurat will greatly facilitate the loading in of scRNA-seq and VisiumHD data. In theory SpatialExperiment or even loading in the h5 files directly is also possible by modifying step 1 and 2 of `runRCTD.R` accordingly.

### Example data

RCTD requires a single-cell reference and the spatial dataset to be deconvolved. The example mouse brain data can be downloaded at ZENODO_LINK. For the scRNA-seq reference, we will be using a subsample of the whole mouse brain dataset from the Allen Brain institute with 10000 genes and ~86,000 cells (the full dataset consists of over 2 million cells). The Visium HD data will be a subsample of the FFPE tissue section from [10x Genomics](https://www.10xgenomics.com/datasets/visium-hd-cytassist-gene-expression-libraries-of-mouse-brain-he-v4). We will be using the 16µm bins subsampled to a small ROI.

The expected results can also be downloaded in the same Zenodo link.

## Steps

It's important to know that the provided scripts are more of a guideline than a flexible pipeline. They have been written for  the UGent HPC cluster, which uses a Slurm scheduler as its backend and Torque as the user interface. Different computing clusters with different schedulers will likely require modification of the job submission script and command.

### 1. Provide input data

We have provided the `runRCTD.R` script which essentially performs the same steps from the [RCTD vignette](https://raw.githack.com/dmcable/spacexr/master/vignettes/spatial-transcriptomics.html). This will be used by the `submit_job.pbs` script which will submit the job to the cluster. This PBS file should be modified to include at least 3 arguments: the path to the scRNA-seq Seurat object (as an rds file), the cell type annotation column name, and the path to the Visium HD data. All possible arguments are:

```
Mandatory arguments:
--sc_input            - path to Seurat object as an rds file
--sp_input            - path to VisiumHD Seurat object as an rds file, OR the h5 file.
                        If the h5 file should be used, provide up to the directory containing the file, e.g.,                            `some/path/binned_outputs/square_008um`.
--annot               - name of metadata column in Seurat object that contains cell type annotation

Optional arguments:
--output_dir          - directory to save the results (default: current directory)
--output_prefix       - file prefix of the output (default: 'rctd_')
--num_cores           - number of cores to use (default: uses number of cores from the `SLURM_NTASKS` environment variable)
--gene_column         - which gene column to use for the VisiumHD data (default: 2, which uses gene symbols. 1 will return the ENSEMBL ids))
  
Example:
Rscript runRCTD.R \
  --sc_input path/to/seurat_object \
  --sp_input path/to/visiumhd \
  --annot celltype \
  --output_dir path/to/output/file \
  --output_prefix exp1_ \
  --num_cores 6
```

Don't forget that the gene names of the scRNA-seq and spatial data should match. In the example data, the ENSEMBL IDs are used. This required some extra data processing steps for the Visium HD data.

### 2. Submit job to the cluster

Once you have provided the necessary arguments, you can submit the job on the cluster. For PBS schedulers, this can be done with

```qsub submit_job.pbs```

We recommend that you test out the setup on the example data or a subset of your data first to see if everything works as expected, since the full 8µm dataset can take a long time. Don't forget to change the walltime and reduce memory limits accordingly to reduce the queue time.

The example data took two minutes with 4 cores on the [UGent HPC doduo cluster](https://docs.hpc.ugent.be/Linux/infrastructure/) which has processing architecture 2x 48-core AMD EPYC 7552 (Rome @ 2.2 GHz). The same cluster took XX hours on the 8µm FFPE brain data.

### 3. Read in and visualize results

You should have two files:
* `*_proportions.tsv` contains the bin-by-cell type proportion information. 
* `*_doublet_info.tsv` contains information on the quality of the predictions, i.e., "singlet", "doublet_certain", "doublet_uncertain", and "reject".

We provide the "read_results.Rmd" markdown file to show a way to read in the result files, as well as a helper function to create the following "scatterbarplot" visualization.