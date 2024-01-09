
#' UMAP-Based Multi-Omics Data Integration
#'
#' @description
#' \code{ubmi} performs unsupervised integration of multi-omics data. 
#' It utilizes Uniform Manifold Approximation and Projection (UMAP) for dimensionality reduction and HDBSCAN for clustering.
#' The method concatenates individual UMAP embeddings of each omics data type, followed by a second UMAP for integrated analysis.
#' Clustering is performed on the integrated data, and optionally, cluster 0 reassignment and feature computation can be done.
#'
#' @param omics A list of omics data matrices, with each matrix representing a different omics data type.
#' @param umap_params A list of parameters for individual UMAP embeddings. Defaults to \code{list(n_neighbors = 15, n_components = 4)}.
#' @param umap_params_conc A list of parameters for concatenated UMAP embeddings. Defaults to \code{list(n_neighbors = 15, n_components = 2)}.
#' @param min_pts The minimum number of points for a cluster in HDBSCAN. Defaults to 10.
#' @param xgboost_params A list of parameters for XGBoost modeling. Defaults to \code{list(lambda = 1, eta = 0.3, gamma = 100, max_depth = 10, subsample = 0.95)}.
#' @param compute_features Logical. Indicates whether to compute metagenes using XGBoost. Defaults to TRUE.
#' @param samples_in_rows Logical. Indicates whether samples are in rows (TRUE) or columns (FALSE) of the omics matrices. Defaults to TRUE.
#'
#' @return A \code{UBMIObject} containing the results of the integration, including factors, clusters, silhouette scores, and metagenes.
#'
#' @examples
#' # Assuming expr, meth, miRNA, and metab are your omics datasets
#' omics_data <- list(expr, meth, miRNA, metab)
#' ubmi_result <- ubmi(omics = omics_data)
#'
#' @export
ubmi <- function(omics,
                 umap_params = list(n_neighbors = 15, n_components = 4, pca = min(nrow(omics[[1]]), ncol(omics[[1]]))),
                 umap_params_conc = list(n_neighbors = 15, n_components = 2),
                 min_pts = 5,
                 xgboost_params = list(lambda = 0, eta = 0.5, gamma = 50, max_depth = 10, subsample = 0.95),
                 compute_features = TRUE,
                 samples_in_rows = TRUE) {
  
  if (samples_in_rows) {
    omics <- lapply(omics, t)
  }
  
  omics <- align_omics(omics)
  omics <- clean_feature_names(omics)
  
  # Factorization
  umap_factors <- lapply(omics, function(x) umap_factorization(umap_params = c(list(X = x), umap_params)))
  concatenated_factors <- list(dplyr::bind_cols(umap_factors, .name_repair = "unique_quiet"))
  message("Computing individual factorizations... OK!")
  
  umap_integrated <- lapply(concatenated_factors, function(x) umap_factorization(umap_params = c(list(X = x), umap_params_conc)))
  message("Computing multi-omics factorization... OK!")
  
  # Clustering
  if (is.null(min_pts)) {
    min_pts <- floor(0.03 * nrow(omics[[1]])) # 5/170 ~ 0.03 (empirical factor)
  }
  min_pts <- max(min_pts, 2)

  hdbscan_labels <- lapply(umap_integrated, function(x) dbscan::hdbscan(x, minPts = min_pts))[[1]]
  
  umap_clusters <- data.frame(umap_integrated[[1]], clust = hdbscan_labels$cluster)
  colnames(umap_clusters)[1:(ncol(umap_clusters) - 1)] <- paste0("UMAP", 1:(ncol(umap_clusters) - 1))
  message("Computing multi-omics clustering... OK!")
  
  # Compute silhouette score
  sil_score <- cluster::silhouette(umap_clusters$clust, dist(umap_clusters[, 1:2]))
  mean_sil_score <- mean(sil_score[, 3])
  
  # Metagenes
  if (compute_features) {
    # features <- as.matrix(dplyr::bind_cols(omics, .name_repair = "unique_quiet"))
    # omics <- c(list(features), omics)
    
    xgboost_fixed_params <- list(objective = "reg:squarederror")
    metagenes <- list()
    for (i in 1:length(omics)) {
      for (j in 1:(ncol(umap_clusters) - 1)) {
        metagenes_tmp <- xgboost_model(x = omics[[i]], y = umap_clusters[,j], xgboost_params = c(xgboost_fixed_params, xgboost_params)) %>% 
          as.data.frame() %>% 
          tibble::rownames_to_column("feature")
        if (j > 1) {
          metagenes_omic <- metagenes_omic %>% dplyr::inner_join(metagenes_tmp, by = "feature")
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
    message("Computing metagenes for each individual omics dataset... OK!")
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
    message("Integration complete!")
    return(ubmi_res)
}

