
#' Align Omics Datasets
#'
#' This function aligns multiple omics datasets by ensuring they have the same set of samples and removing features with zero standard deviation.
#' It also transposes each dataset so that samples become rows and features become columns, which is a common format for omics data analysis.
#'
#' @param omics A list of data frames or matrices, where each element represents a different omics dataset.
#'              Each dataset should have the same set of samples (columns) but can have different features (rows).
#'
#' @return A list of data frames or matrices, where each dataset is transposed and aligned in terms of samples.
#'         Features with zero standard deviation across samples are removed.
#'
#' @export
align_omics <- function(omics) {
  samples <- colnames(omics[[1]])
  for (j in 1:length(omics)) {
    samples <- intersect(samples, colnames(omics[[j]]))
  }
  for (j in 1:length(omics)) {
    omics[[j]] <- omics[[j]][, samples]
    omics[[j]] <- omics[[j]][, which(apply(omics[[j]], 2, sd) > 0)]
  }
  
  omics <- lapply(omics, FUN = function(x){t(x)})
  
  return(omics)
}

#' Clean Feature Names
#'
#' This function standardizes feature names by replacing certain characters with underscores.
#'
#' @param omics A list of data frames or matrices representing omics datasets.
#'              Each element in the list should have column names that potentially need cleaning.
#'
#' @return A list of data frames or matrices, similar to the input but with cleaned column names.
#'         All periods ('.') and hyphens ('-') in the column names are replaced with underscores ('_').
#'
#' @export
clean_feature_names <- function(omics) {
  omics <- lapply(omics, FUN = function(x){
    colnames(x) <- gsub("\\.", "_", colnames(x))
    colnames(x) <- gsub("\\-", "_", colnames(x))
    x <- as.matrix(x)
  })
  
  return(omics)
}

#' UMAP Factorization
#'
#' This function performs Uniform Manifold Approximation and Projection (UMAP) 
#' for dimensionality reduction, using specified parameters.
#'
#' @param umap_params A list of parameters to pass to the UMAP function from the `uwot` package.
#'                    If not specified, default parameters are used.
#'
#' @return A data frame containing the UMAP factorization results. 
#'         Each row corresponds to a data point, and the columns are the UMAP dimensions.
#'         
#' @export
umap_factorization <- function(umap_params = list()) {
  umap_factors <- do.call(uwot::umap, umap_params)
  umap_factors <- data.frame(umap_factors)
  
  return(umap_factors)
}

#' XGBoost Model with SHAP Values
#'
#' This function builds an XGBoost model using the provided training data and parameters and 
#' calculates SHAP (SHapley Additive exPlanations) values to interpret the model.
#'
#' @param x A `matrix` or data frame of predictor variables (features).
#' @param y A `vector` of response variables (target).
#' @param xgboost_params A `list` of parameters for the XGBoost model. 
#'                       Defaults to an empty list, which uses the XGBoost default parameters.
#'
#' @return A matrix of SHAP contribution values for each feature in the dataset.
#'         This matrix provides insights into how each feature contributes to the model's predictions.
#'
#' @export
xgboost_model <- function(x, y, xgboost_params = list()) {
  mod <- xgboost::xgboost(data = x,
                          label = as.matrix(y),
                          params = xgboost_params, 
                          nrounds = 10,
                          verbose = FALSE, 
                          nthread = parallel::detectCores() - 2,
                          early_stopping_rounds = 8)
  
  shap_values <- SHAPforxgboost::shap.values(xgb_model = mod, X_train = x)
  shap_contrib <- as.matrix(shap_values$mean_shap_score)

  return(shap_contrib)
}

#' Drop Specified Clusters
#'
#' This function removes specified clusters from an `UBMIObject`.
#'
#' @param object An `UBMIObject`.
#' @param clusters A vector of cluster IDs (numeric or character) that should be removed from the object.
#'                 Default is 0, which means it will remove the cluster(s) labeled as '0'.
#'
#' @return An object similar to the input but with the specified clusters and their corresponding data removed.
#'         The function performs in-place modification but also returns the modified object for convenience.
#'         If the object is not valid after modification, the function will not return it.
#'
#' @export
drop_clusters <- function(object,
                          clusters = 0) {
  
  nonzero <- which(!(object@factors$clust %in% clusters))
  
  object@factors <- object@factors[nonzero ,]
  object@clusters <- object@clusters[nonzero]
  object@metagenes_factor1 <- object@metagenes_factor1[nonzero ,]
  object@metagenes_factor2 <- object@metagenes_factor2[nonzero ,] 
  
  if(validObject(object))
    return(object)
}

#' Bootstrap Omics Data
#'
#' This function generates synthetic omics data samples by applying random perturbations to existing data.
#' It aims to augment the dataset by creating additional, slightly varied samples, which can be useful in various analyses like machine learning.
#'
#' @param data A data frame.
#' @param n The proportion of the original dataset size to generate as synthetic samples.
#'          For example, n = 1 will create an equal number of synthetic samples as in the original dataset.
#' @param perturbation A numeric value specifying the magnitude of random perturbation to apply to each feature of the synthetic samples.
#'                     Default is set to 0.1. The perturbation is uniformly applied within the range [-perturbation, +perturbation].
#'
#' @return A data frame or matrix combining the original and synthetic samples. 
#'         The number of rows will be the sum of the original dataset and the number of synthetic samples created.
#'
#' @export
bootstrap_omics <- function(data, n, perturbation = 0.1) {
  # Number of synthetic samples to create
  n_samples_to_create <- round(nrow(data) * n)
  
  # Create synthetic samples
  synthetic_indices <- sample(1:nrow(data), n_samples_to_create, replace = TRUE)
  synthetic_samples <- data[synthetic_indices ,]
  
  # Apply perturbation to the synthetic samples
  # The perturbation is applied by adding a random value between -perturbation and +perturbation
  perturbation_matrix <- runif(n = n_samples_to_create * ncol(data), min = -perturbation, max = perturbation)
  
  # Reshape the perturbation matrix to match the dimensions of synthetic_samples
  perturbation_matrix <- matrix(perturbation_matrix, nrow = n_samples_to_create, ncol = ncol(data))
  
  # Apply the perturbation
  synthetic_samples <- synthetic_samples + perturbation_matrix
  
  combined_data <- rbind(data, synthetic_samples)
  
  return(combined_data)
}

#' Reassign Cluster Zero
#'
#' This function reassigns data points in cluster 0 to the nearest cluster based on either centroid distance or individual point distance.
#'
#' @param data A data frame.
#' @param nearest_centroid Logical, defaulting to FALSE. If TRUE, reassigns points in cluster 0 to the nearest centroid 
#'                         of other clusters. If FALSE, reassigns based on the nearest point in any other cluster.
#'
#' @return A data frame identical to the input but with updated cluster assignments for points initially in cluster 0.
#'
#' @export
reassign_cluster_zero <- function(data, nearest_centroid = FALSE) {
  # Calculate the centroids of non-zero clusters
  non_zero_clusters <- data %>% 
    dplyr::filter(clust != 0)
  
  centroids <- non_zero_clusters %>%
    dplyr::group_by(clust) %>% 
    dplyr::summarise(centroid_x = median(UMAP1, na.rm = TRUE),
                     centroid_y = median(UMAP2, na.rm = TRUE))
  
  # Get data points in cluster 0
  cluster_zero_points <- data %>% 
    dplyr::filter(clust == 0)
  
  if (nrow(cluster_zero_points) > 0) {
    if (nearest_centroid) {
      # For each point in cluster 0, calculate the distance to each centroid and reassign
      for(i in 1:nrow(cluster_zero_points)) {
        point <- cluster_zero_points[i, c("UMAP1", "UMAP2")]
        distances <- apply(centroids[, -1], 1, function(centroid) sqrt((point[1] - centroid[1])^2 + (point[2] - centroid[2])^2))
        nearest_cluster <- centroids[which.min(unlist(distances)), "clust"]
        cluster_zero_points[i, "clust"] <- nearest_cluster
      }
    } else {
      # For each point in cluster 0, calculate the distance to each point and reassign
      for (i in 1:nrow(cluster_zero_points)) {
        point <- cluster_zero_points[i, c("UMAP1", "UMAP2")]
        distances <- apply(non_zero_clusters[, -3], 1, function(data) sqrt((point[1] - data[1])^2 + (point[2] - data[2])^2))
        nearest_cluster <- non_zero_clusters[which.min(unlist(distances)), "clust"]
        cluster_zero_points[i, "clust"] <- nearest_cluster
      }
    }
    
    # Step 4: Combine the data with re-assigned clusters and non-zero clusters
    data <- dplyr::bind_rows(non_zero_clusters, cluster_zero_points)
  }
  
  return(data)
}

