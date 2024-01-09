#' Align omics
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

#' Clean names
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

#' UMAP factors
#'
#' @export
umap_factorization <- function(umap_params = list()) {
  umap_factors <- do.call(uwot::umap, umap_params)
  umap_factors <- data.frame(umap_factors)
  
  return(umap_factors)
}

#' XGBOOST
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

#' Drop clusts
#'
#' @export
drop_clusters <- function(object,
                          clusters = 0,
                          ...) {
  
  nonzero <- which(!(object@factors$clust %in% clusters))
  
  object@factors <- object@factors[nonzero ,]
  object@clusters <- object@clusters[nonzero]
  object@metagenes_factor1 <- object@metagenes_factor1[nonzero ,]
  object@metagenes_factor2 <- object@metagenes_factor2[nonzero ,] 
  
  if(validObject(object))
    return(object)
}

#' Bootstrap
#'
#' @export
bootstrap_omics <- function(data, n, perturbation = 0.1) {
  # Step 1: Calculate the number of synthetic samples to create
  n_samples_to_create <- round(nrow(data) * n)
  
  # Step 2: Sample rows with replacement to create synthetic samples
  synthetic_indices <- sample(1:nrow(data), n_samples_to_create, replace = TRUE)
  synthetic_samples <- data[synthetic_indices ,]
  
  # Step 3: Apply perturbation to the synthetic samples
  # The perturbation is applied by adding a random value between -perturbation and +perturbation
  perturbation_matrix <- runif(n = n_samples_to_create * ncol(data), min = -perturbation, max = perturbation)
  
  # Reshape the perturbation matrix to match the dimensions of synthetic_samples
  perturbation_matrix <- matrix(perturbation_matrix, nrow = n_samples_to_create, ncol = ncol(data))
  
  # Apply the perturbation
  synthetic_samples <- synthetic_samples + perturbation_matrix
  
  combined_data <- rbind(data, synthetic_samples)
  
  return(combined_data)
}

#' Reassign
#'
#' @export
reassign_cluster_zero <- function(data, nearest_centroid = FALSE) {
  # Step 1: Calculate the centroids of non-zero clusters
  non_zero_clusters <- data %>% 
    dplyr::filter(clust != 0)
  
  centroids <- non_zero_clusters %>%
    dplyr::group_by(clust) %>% 
    dplyr::summarise(centroid_x = median(UMAP1, na.rm = TRUE),
                     centroid_y = median(UMAP2, na.rm = TRUE))
  
  # Step 2: Get data points in cluster 0
  cluster_zero_points <- data %>% 
    dplyr::filter(clust == 0)
  
  if (nrow(cluster_zero_points) > 0) {
    if (nearest_centroid) {
      # Step 3: For each point in cluster 0, calculate the distance to each centroid and reassign
      for(i in 1:nrow(cluster_zero_points)) {
        point <- cluster_zero_points[i, c("UMAP1", "UMAP2")]
        distances <- apply(centroids[, -1], 1, function(centroid) sqrt((point[1] - centroid[1])^2 + (point[2] - centroid[2])^2))
        nearest_cluster <- centroids[which.min(unlist(distances)), "clust"]
        cluster_zero_points[i, "clust"] <- nearest_cluster
      }
    } else {
      # Step 3: For each point in cluster 0, calculate the distance to each point and reassign
      for(i in 1:nrow(cluster_zero_points)) {
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

