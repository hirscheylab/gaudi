
library(tidyverse)

load("/Users/pol/Dropbox/aml_clusters/data/target_multiomics.RData")

omics <- list(expression_final[,-1], methylation_final[,-1], mirna_final[,-1])

data_target <- ubmi(omics)

# Function to calculate pairwise distances between two sets of embeddings
compute_inter_modal_distances <- function(embeddings1, embeddings2) {
  inter_modal_dist <- matrix(0, nrow(embeddings1), nrow(embeddings2))
  for(i in 1:nrow(embeddings1)) {
    for(j in 1:nrow(embeddings2)) {
      inter_modal_dist[i,j] <- sqrt(sum((embeddings1[i,] - embeddings2[j,])^2))
    }
  }
  return(inter_modal_dist)
}

# Compute inter-modal distances
# compute_inter_modal_distances <- function(embeddings1, embeddings2) {
#   as.vector(outer(1:nrow(embeddings1), 1:nrow(embeddings2), FUN = function(i, j) {
#     sqrt(sum((embeddings1[i, ] - embeddings2[j, ])^2))
#   }))
# }

# integrate_omics <- function(omics) {
#   n_datasets <- length(omics)
#   
#   # Step 1: Compute individual UMAP embeddings for each dataset
#   umap_list <- lapply(omics, function(data) uwot::umap(data))
#   
#   # Step 2: Calculate intra- and inter-modal distances
#   # Define the total size of the joint_dist_matrix
#   total_size <- sum(sapply(omics, nrow))
# 
#   # Initialize joint_dist_matrix
#   joint_dist_matrix <- matrix(0, nrow = total_size, ncol = total_size)
# 
#   # Fill joint_dist_matrix
#   n_datasets <- length(omics)
# 
#   # Fill joint_dist_matrix
#   start_row <- 1
# 
#   for(i in 1:(n_datasets - 1)) { # we use n_datasets-1 because the last dataset won't need to compare to any other datasets
# 
#     end_row <- start_row + nrow(omics[[i]]) - 1
# 
#     cat("\ni: ", i, "\nstart_row: ", start_row, "\nend_row: ", end_row, "\n")
# 
#     # Intra-modal distances
#     joint_dist_matrix[start_row:end_row, start_row:end_row] <- as.matrix(dist(omics[[i]]))
# 
#     for(j in (i + 1):n_datasets) {
# 
#       start_row_j <- end_row + 1
#       end_row_j <- start_row_j + nrow(omics[[j]]) - 1
# 
#       cat("\nj: ", j, "\nstart_row_j: ", start_row_j, "\nend_row_j: ", end_row_j, "\n")
# 
#       inter_modal_dist <- compute_inter_modal_distances(umap_list[[i]], umap_list[[j]])
# 
#       # Make sure the inter_modal_dist matrix dimensions match the target joint_dist_matrix dimensions
#       if (dim(inter_modal_dist)[1] != (end_row - start_row + 1) || dim(inter_modal_dist)[2] != (end_row_j - start_row_j + 1)) {
#         stop(paste0("Dimension mismatch at i=", i, ", j=", j))
#       }
# 
#       joint_dist_matrix[start_row:end_row, start_row_j:end_row_j] <- inter_modal_dist
#       joint_dist_matrix[start_row_j:end_row_j, start_row:end_row] <- t(inter_modal_dist)
# 
#     }
# 
#     start_row <- end_row + 1
#   }
#  
#   # Step 3: Apply UMAP on joint distance matrix
#   joint_umap <- uwot::umap(joint_dist_matrix)
#   
#   return(joint_umap)
# }

integrate_omics <- function(omics) {
  n_datasets <- length(omics)
  
  # Step 1: Compute individual UMAP embeddings for each dataset
  umap_list <- lapply(omics, function(data) uwot::umap(data))
  
  # Step 2: Calculate intra- and inter-modal distances
  # Define the total size of the joint_dist_matrix
  total_size <- sum(sapply(omics, nrow))
  
  # Initialize joint_dist_matrix
  joint_dist_matrix <- compute_inter_modal_distances(umap_list[[1]], umap_list[[2]])
  joint_dist_matrix <- compute_inter_modal_distances(umap_list[[1]], umap_list[[2]])
  
  # Step 3: Apply UMAP on joint distance matrix
  joint_umap <- uwot::umap(joint_dist_matrix)
  
  return(joint_umap)
}

# Usage
result <- integrate_omics(omics)

ggplot2::ggplot(as.data.frame(result), ggplot2::aes(V1, V2)) +
  ggplot2::geom_point(size = 3) +
  ggplot2::theme_bw()

factors_target <- data.frame(id = expression_final[,1], result)

####

library(MOFA2)

MOFAobject <- MOFA2::create_mofa(lapply(omics, function(data) t(data)))
DataOptions <- MOFA2::get_default_data_options(MOFAobject)
ModelOptions <- MOFA2::get_default_model_options(MOFAobject)
ModelOptions$num_factors <- 2
TrainOptions <- MOFA2::get_default_training_options(MOFAobject)

MOFAobject <- MOFA2::prepare_mofa(
  MOFAobject,
  data_options = DataOptions,
  model_options = ModelOptions,
  training_options = TrainOptions
)

MOFAobject <- MOFA2::run_mofa(MOFAobject) #, use_basilisk = TRUE)
factors_mofa <- MOFA2::get_factors(MOFAobject)

factorizations_mofa <- factors_mofa$group1

ggplot2::ggplot(as.data.frame(factorizations_mofa), ggplot2::aes(Factor1, Factor2)) +
  ggplot2::geom_point(size = 3) +
  ggplot2::theme_bw() 

