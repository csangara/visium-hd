# Deconvolution of Visium HD data

We will use the tool RCTD -- now implemented as part of the spacexr package -- which is also what was done in the Visium HD announcement paper by [Oliveira et al. (2025)](https://www.nature.com/articles/s41588-025-02193-3). This is because multiple benchmarking studies have found it to be among the best-perfoming deconvolution methods while remaining relatively scalable. Another advantage is that it performs well without parameter tuning and is relatively straightforward to run.

We highly recommend that users first follow the [RCTD vignette](https://raw.githack.com/dmcable/RCTD/dev/vignettes/spatial-transcriptomics.html) in the spacexr repository to make sure that they understand the various functions and how the tool works. This tutorial will mostly focus on how to deploy the tool on Visium HD data.

Since Visium HD datasets contain substantially more observations than those of standard Visium, **it is highly recommended that the tool is run on a computing cluster such as a HPC**.

## Requirements

### Packages

#### RCTD

Although RCTD is now available through Bioconductor, we recommend installing it from GitHub to access the optimized implementation developed by the 10x Genomics team ([dmcable/spacexr\#206](https://github.com/dmcable/spacexr/pull/206)). Note that this optimization is only available for `doublet_mode`, which limits predictions to a maximum of two cell types per bin.

Install the sped-up version of RCTD on your computing cluster as follows: 

```
if(!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools") 
}
options(timeout = 600000000) # set this to avoid timeout error
devtools::install_github("jpromeror/spacexr", ref = "HD", build_vignettes = FALSE)
```

#### Seurat (v5)

Seurat is not strictly required because RCTD uses its own data structures. However, we still recommend it because Seurat will simplify loading and preprocessing both scRNA-seq and Visium HD datasets. If other data formats must be used, you will have to modify Steps 1 and 2 in the script `runRCTD.R` accordingly.

### Example data

We provide a subset of the brain dataset we used in our paper in the following link: [https://zenodo.org/records/20716044]. You can download `example_data.zip` which contains the following:
* A scRNA-seq reference containing 1,933 cells, 10,000 genes, and 7 cell types (`class`). This is a subsample of the [whole mouse brain v3 dataset](https://alleninstitute.github.io/abc_atlas_access/descriptions/WMB-10Xv3.html) from the Allen Brain Institute.
* A VisiumHD object at 16µm resolution, subsampled to 289 bins and 6,951 genes. The original dataset is from the [10x Genomics mouse brain FFPE tissue](https://www.10xgenomics.com/datasets/visium-hd-cytassist-gene-expression-libraries-of-mouse-brain-he-v4).
* `.tsv` files containing the expected results

## Steps

It's important to know that the provided scripts are more of a guideline than a flexible pipeline. The provided scripts were written for the UGent HPC cluster. If your cluster uses a different scheduler or job submission system, you may need to modify the submission script accordingly.

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

Be mindful that the gene names used in the scRNA-seq and spatial dataset must match. In the example data, Ensembl IDs are used, which required some extra data loading steps for the Visium HD data.

### 2. Submit job to the cluster

Once you have provided the necessary arguments, you can submit the job on the cluster. In the UGent HPC, this is done with

```qsub submit_job.pbs```

We recommend  testing the workflow on the example data or a subset of your data first to see if everything works as expected (and request fewer resources and walltime to reduce the queue time).

The example data took two minutes with 4 cores on the [UGent HPC doduo cluster](https://docs.hpc.ugent.be/Linux/infrastructure/) which is equipped with 2x 48-core AMD EPYC 7552 processors (Rome @ 2.2 GHz). The same cluster took 9 hours on the full 8µm data.

If the job freezes when running on multiple cores but works on a single core, check if any of these issues are relevant:
* https://github.com/dmcable/spacexr/issues/82
* https://github.com/dmcable/spacexr/issues/215
* https://github.com/dmcable/spacexr/issues/227


### 3. Read in and visualize results

You should now have two files:
* `*_proportions.tsv` contains the bin-by-cell type proportions 
* `*_doublet_info.tsv` contains information about the confidence of the predictions, i.e., "singlet", "doublet_certain", "doublet_uncertain", and "reject".

We provide the "read_results.Rmd" markdown file to show a way to read in the result files, as well as a helper function to create the following "scatterbarplot" visualization.

## Full data

If you want to try out the pipeline on a full 8µm Visium HD data, we also provide this in the same [Zenodo link](https://zenodo.org/records/20716044) under `full_data.zip`. This contains:
* The scRNA-seq reference data we used in our paper, with 86,810 cells, 32,285 genes, and 34 cell types.
* The Visium HD object at 8µm resolution, consisting of 393,543 bins and 19,059 genes.

The same steps apply when using the full dataset. One difference is that the `runRCTD.R` script has to be modified at `create.RCTD` (line 87) to add the `CELL_MIN_INSTANCE=5` parameter, as follows:

```
RCTD_deconv <- create.RCTD(spatialRNA = spatialRNA_obj_visium,
                           reference = reference_obj,
                           max_cores = num_cores,
                           CELL_MIN_INSTANCE = 5)
```

This is because two of the cell types only contained five cells in the reference, whereas the default minimum number of cells required per cell type is 25. Alternatively, you can also remove these cell types from the reference dataset before running RCTD.