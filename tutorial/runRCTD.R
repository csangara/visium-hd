#!/usr/bin/env Rscript
Sys.setenv(RETICULATE_MINICONDA_ENABLED = "FALSE")

library(spacexr)
library(Matrix)
library(Seurat)
library(dplyr)

# Get user input
user_args <- R.utils::commandArgs(trailingOnly=TRUE, asValues=TRUE)

cat("User arguments:\n")
print(user_args)

#### 0. INITIALIZE ARGUMENTS ####
# Check file name
file_prefix <- ifelse(is.null(user_args$output_prefix), "rctd_", user_args$output_prefix)
filename_doublet_info <- paste0(file_prefix, "doublet_info.tsv")
filename_proportion <- paste0(file_prefix, "proportions.tsv")

# Output directory - check validity
if (!is.null(user_args$output_dir)){
  stopifnot("`output_dir` directory does not exist, please enter a valid directory." = dir.exists(user_args$output_dir))
  output_dir <- user_args$output_dir
} else {
  output_dir <- getwd()
  cat("No output directory provided, so file will be written at:", output_dir, "\n")
}
stopifnot("You do not have write access to the directory." = file.access(output_dir, 2) == 0)

# Get cores
if (!is.null(user_args$num_cores)){
  num_cores <- as.numeric(user_args$num_cores)
} else {
  num_cores <- parallel::detectCores()
}
cat("Using", num_cores, "cores...\n")

#### 1. CONSTRUCT SCRNA-SEQ RCTD OBJECT ####
cat("Reading input scRNA-seq reference from", user_args$sc_input, "\n")
seurat_obj_scRNA <- readRDS(user_args$sc_input)
ncelltypes <- length(unique(seurat_obj_scRNA[[user_args$annot, drop=TRUE]]))
cat("Found ", ncelltypes, "cell types in the reference.\n")

cat("Converting to Reference object...\n")
cell_types <- stringr::str_replace_all(seurat_obj_scRNA[[user_args$annot, drop=TRUE]],
                                       "[/ .]", "") # Replace prohibited characters
names(cell_types) <- colnames(seurat_obj_scRNA)
reference_obj <- Reference(counts = GetAssayData(seurat_obj_scRNA, layer="counts"),
                           cell_types = as.factor(cell_types))

#### 2. CONSTRUCT SPATIAL RCTD OBJECT ####
cat("Reading input spatial data from", user_args$sp_input, "\n")
sp_input <- user_args$sp_input

# If it is an rds file (assumed to be Seurat object)
if (grepl("*.rds$", basename(sp_input), ignore.case = TRUE)){
  spatial_data <- readRDS(sp_input)
} else {
  spatial_data <- Load10X_Spatial(user_args$sp_input)
}

cat("Converting spatial data to SpatialRNA object...\n")
# Check if there is images slot
use_fake_coords <- length(spatial_data@images) == 0
coords <- NULL
if (length(spatial_data@images)){
  coords <- GetTissueCoordinates(spatial_data) %>% 
    # Only get numeric columns
    dplyr::select(where(is.numeric))
}

spatialRNA_obj_visium <- SpatialRNA(coords = coords,
                                    counts = GetAssayData(spatial_data, layer="counts"),
                                    use_fake_coords = use_fake_coords)

#### 3. RUN RCTD ####
cat("Running RCTD with", num_cores, "cores...\n")
RCTD_deconv <- create.RCTD(spatialRNA = spatialRNA_obj_visium,
                           reference = reference_obj,
                           max_cores = num_cores)
RCTD_deconv <- run.RCTD(RCTD_deconv, doublet_mode = "doublet")

#### 4. EXPORT RESULTS ####
cat("Printing results...\n")

#### 4.1. OUTPUT DOUBLET LABELS AND LIKELIHOOD SCORE ####
write.table(RCTD_deconv@results$results_df %>% tibble::rownames_to_column("spot"),
            file=file.path(output_dir, filename_doublet_info),
            sep="\t", quote=FALSE, row.names=FALSE)


#### 4.2 OUTPUT PROPORTION MATRIX ####
# Get doublet proportions
weights_doublet <- RCTD_deconv@results$weights_doublet %>% 
  data.frame(., row.names = rownames(.)) %>% tibble::rownames_to_column("spot") %>%
  tidyr::pivot_longer(cols =-spot, names_to="type", values_to="proportion")

# Get doublet labels
labels_df <- RCTD_deconv@results$results_df %>%
  select(spot_class, first_type, scond_type) %>% 
  tibble::rownames_to_column("spot") %>% 
  tidyr::pivot_longer(-c(spot, spot_class), names_to="type", values_to="label")

deconv_matrix <- dplyr::left_join(weights_doublet, labels_df, by=c("spot", "type")) %>% 
  # Filter out second_type if the spot is classified as a singlet
  dplyr::filter(!(type == "scond_type" & spot_class == "singlet")) %>% 
  # Replace proportions with 1 if singlet
  dplyr::mutate(proportion = replace(proportion, spot_class == "singlet", 1)) %>% 
  tidyr::pivot_wider(id_cols = spot, names_from="label", values_from="proportion", values_fill = 0) %>% 
  tibble::column_to_rownames("spot")

# Print removed bins
if (nrow(deconv_matrix) != ncol(spatial_data)){
  message("The following rows were removed, due to low number of UMIs: ",
          paste0("'", colnames(spatial_data)[!colnames(spatial_data) %in% rownames(deconv_matrix)], "'", collapse=", "))
}

# Check for missing cell types
celltypes <- RCTD_deconv@cell_type_info$info[[2]]
if (!all(celltypes %in% colnames(deconv_matrix))){
  missing <- celltypes[which(!celltypes %in% colnames(deconv_matrix))]
  cat("Cell types with zero abundance in all spots:", paste(missing, collapse=", "), "\n")
  
  # Add columns of missing cell types
  deconv_matrix[, missing] <- 0
  deconv_matrix <- as.matrix(deconv_matrix)
}

write.table(deconv_matrix, file=file.path(output_dir, filename_proportion),
            sep="\t", quote=FALSE, row.names=TRUE)