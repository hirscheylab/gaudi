
#' UMAP-based Multi-omics Integration
#'
#' @description
#' \code{ubmi} performs unsupervised integration of multi-omics data. 
#' It utilizes Uniform Manifold Approximation and Projection (UMAP) for dimensionality reduction and HDBSCAN for clustering.
#' The method concatenates individual UMAP embeddings of each omics data type, followed by a second UMAP for integrated analysis.
#' Clustering is performed on the integrated data, and optionally, cluster 0 reassignment and feature computation (metagenes) can be done.
#'
#' @param omics A `list` of omics data frames or matrices.
#' @param umap_params Parameters for UMAP factorization of individual omics datasets.
#' @param umap_params_conc Parameters for UMAP factorization of concatenated omics data.
#' @param min_pts Minimum number of points per cluster in HDBSCAN clustering.
#' @param xgboost_params Parameters for XGBoost model used in metagenes computation.
#' @param compute_features Logical indicating whether to compute metagenes.
#' @param combine_omics Logical indicating whether to combine combined omics metagenes.
#' @param clean_feature_names Logical indicating whether to clean feature names in omics datasets.
#' @param samples_in_rows Logical indicating whether samples are in rows (default) or columns.
#' @param reassign_cluster_zero Logical indicating whether to reassign data points in cluster 0.
#'
#' @return An object of class `UBMIObject` containing the results of the multi-omics integration and analysis.
#'         This includes factors, clusters, silhouette scores, individual factors, metagenes, and the \code{ubmi} package version.
#'
#' @examples
#' omics_list <- list(data.frame(matrix(rnorm(200), ncol = 20)),
#'                    data.frame(matrix(rnorm(200), ncol = 20)))
#' ubmi_result <- ubmi(omics_list, samples_in_rows = FALSE)
#'
#' @export
ubmi <- function(omics,
                 umap_params = list(n_neighbors = 15, n_components = 4, pca = min(nrow(omics[[1]]), ncol(omics[[1]]))),
                 umap_params_conc = list(n_neighbors = 15, n_components = 2),
                 min_pts = 5,
                 xgboost_params = list(lambda = 0, eta = 0.5, gamma = 50, max_depth = 10, subsample = 0.95),
                 compute_features = TRUE,
                 combine_omics = FALSE,
                 clean_feature_names = FALSE,
                 samples_in_rows = TRUE,
                 reassign_cluster_zero = FALSE) {
  
  if (samples_in_rows) {
    omics <- lapply(omics, t)
  }
  
  omics <- align_omics(omics)
  
  if (clean_feature_names) {
    omics <- clean_feature_names(omics)
  }
  
  # Factorization
  message("Computing individual factorizations...")
  umap_factors <- lapply(omics, function(x) umap_factorization(umap_params = c(list(X = x), umap_params)))
  concatenated_factors <- list(dplyr::bind_cols(umap_factors, .name_repair = "unique_quiet"))
  message("Done!")
  
  message("Computing multi-omics factorization...")
  umap_integrated <- lapply(concatenated_factors, function(x) umap_factorization(umap_params = c(list(X = x), umap_params_conc)))
  message("Done!")
  
  # Clustering
  message("Computing multi-omics clustering...")
  if (is.null(min_pts)) {
    min_pts <- floor(0.03 * nrow(omics[[1]]))
  }
  min_pts <- max(min_pts, 2)

  hdbscan_labels <- lapply(umap_integrated, function(x) dbscan::hdbscan(x, minPts = min_pts))[[1]]
  
  umap_clusters <- data.frame(umap_integrated[[1]], clust = hdbscan_labels$cluster)
  colnames(umap_clusters)[1:(ncol(umap_clusters) - 1)] <- paste0("UMAP", 1:(ncol(umap_clusters) - 1))
  message("Done!")
  
  if (reassign_cluster_zero) {
    umap_clusters <- reassign_cluster_zero(umap_clusters, nearest_centroid = FALSE)
  }
  
  # Compute silhouette score
  if (length(unique(umap_clusters$clust)) > 1) {
    sil_score <- cluster::silhouette(umap_clusters$clust, dist(umap_clusters[, 1:2]))
    mean_sil_score <- mean(sil_score[, 3])
  } else {
    mean_sil_score <- NULL
  }
  
  # Metagenes
  if (compute_features) {
    message("Computing metagenes...")
    if (combine_omics) {
      features <- as.matrix(dplyr::bind_cols(omics, .name_repair = "unique_quiet"))
      omics <- c(list(features), omics) 
    }
    
    xgboost_fixed_params <- list(objective = "reg:squarederror")
    metagenes <- list()
    for (i in 1:length(omics)) {
      for (j in 1:(ncol(umap_clusters) - 1)) {
        metagenes_tmp <- xgboost_model(x = omics[[i]], 
                                       y = umap_clusters[,j], 
                                       xgboost_params = c(xgboost_fixed_params, xgboost_params)) %>% 
          as.data.frame() %>% 
          tibble::rownames_to_column("feature")
        if (j > 1) {
          metagenes_omic <- metagenes_omic %>%
            dplyr::inner_join(metagenes_tmp, by = "feature")
        } else {
          metagenes_omic <- metagenes_tmp
        }
      }
      colnames(metagenes_omic) <- c("feature", paste0("contrib", 1:(ncol(umap_clusters) - 1)))
      metagenes_omic <- metagenes_omic %>% 
        tibble::column_to_rownames("feature") %>% 
        dplyr::arrange(dplyr::desc(abs(contrib1)))
      metagenes[[i]] <- metagenes_omic
    }
    message("Done!")
  }
  
  ubmi_res <- new("UBMIObject",
                  factors = umap_clusters,
                  clusters = umap_clusters$clust,
                  silhouette_score = mean_sil_score,
                  individual_factors = umap_factors,
                  metagenes = if(compute_features) metagenes else list(),
                  ubmiVersion = as.character(packageVersion("ubmi"))
  )
  
  if (validObject(ubmi_res))
    message("Multi-omics integration complete!")
    return(ubmi_res)
}

