## Downsampling experiment on brain data

* 1_1: Generates 3 replicates of simulated spatial data from the scRNA-seq brain data (also used in Spotless). Each spot is randomly downsampled to have a UMI belonging to a Gaussian distribution of 1000 mean and 1000 sd. RCTD is run automatically on the generated datasets.
* 1_2: Runs NNLS on the data. To avoid this step, it is possible to also adapt `1_1_generate_and_run.config` to `methods = "rctd,nnls"` (instead of `methods = "rctd"`)
* 1_3: Runs Scrublet on the data
* 2_*: Evaluates results.
* ex: yaml files listing actual mean and standard deviation of actual Visium and VisiumHD data
   
