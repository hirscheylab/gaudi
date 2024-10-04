
#' Consensus Group Aggregation via UMAP Data Integration
#'
#' @description
#' \code{c_gaudi} performs unsupervised integration of multi-omics data. 
#' It utilizes Uniform Manifold Approximation and Projection (UMAP) for dimensionality reduction and HDBSCAN for clustering.
#' The method concatenates individual UMAP embeddings of each omics data type, followed by a second UMAP for integrated analysis.
#' Clustering is performed on the integrated data, and optionally, cluster 0 reassignment and feature computation (metagenes) can be done.
#'
#' @param omics A `list` of omics data frames or matrices.
#' @param n_max Numeric value indicating the total number of latent spaces to be computed (the final latent space will be based on all these spaces).
#' @param umap_params Parameters for UMAP factorization of individual omics datasets.
#' @param umap_params_conc Parameters for UMAP factorization of concatenated omics data.
#' @param min_pts Minimum number of points per cluster in HDBSCAN clustering.
#' @param xgboost_params Parameters for XGBoost model used in metagenes computation.
#' @param compute_features Logical indicating whether to compute metagenes.
#' @param combine_omics Logical indicating whether to combine combined omics metagenes.
#' @param clean_feature_names Logical indicating whether to clean feature names in omics datasets.
#' @param samples_in_rows Logical indicating whether samples are in rows (default) or columns.
#' @param reassign_cluster_zero Logical indicating whether to reassign data points in cluster 0.
#' @param method Character indicating the regression model used to extract features. Options are "xgboost" and "rf".
#' 
#' @return An object of class `GAUDIObject` containing the results of the multi-omics integration and analysis.
#'         This includes factors, clusters, silhouette scores, individual factors, metagenes, and the \code{gaudi} package version.
#'
#' @examples
#' omics_list <- list(data.frame(matrix(rnorm(200), ncol = 20)),
#'                    data.frame(matrix(rnorm(200), ncol = 20)))
#' c_gaudi_result <- c_gaudi(omics_list, n_max = 3, samples_in_rows = FALSE)
#'
#' @export
c_gaudi <- function(omics,
                    n_max = 10,
                    umap_params = list(n_neighbors = 15, n_components = 4),
                    umap_params_conc = list(n_neighbors = 15, n_components = 2),
                    min_pts = NULL,
                    xgboost_params = list(lambda = 0, eta = 1, gamma = 100, max_depth = 10, subsample = 0.95),
                    compute_features = TRUE,
                    combine_omics = FALSE,
                    clean_feature_names = FALSE,
                    samples_in_rows = TRUE,
                    reassign_cluster_zero = FALSE,
                    method = "xgboost") {
  
  if (samples_in_rows) {
    omics <- lapply(omics, t)
  }
  
  omics <- align_omics(omics)
  
  if (clean_feature_names) {
    omics <- clean_feature_names(omics)
  }
  
  message("Computing factorizations...")
  cgaudi_loop <- list()
  for (i in 1:n_max) {
    suppressMessages({
      cgaudi_loop[[i]] <- gaudi::gaudi(omics,
                                       umap_params = umap_params,
                                       umap_params_conc = umap_params_conc,
                                       min_pts = min_pts,
                                       compute_features = FALSE,
                                       combine_omics = FALSE,
                                       clean_feature_names = FALSE,
                                       samples_in_rows = TRUE,
                                       reassign_cluster_zero = FALSE)
    })
    message(paste0("Iteration ", i, "/", n_max, " completed."))
  }
  message("Done!")
  
  idx_max <- which.max(unlist(lapply(cgaudi_loop, function(x){x@silhouette_score})))
  
  # Reference
  reference <- cgaudi_loop[[idx_max]]@factors[,-3]
  
  message("Computing maximum similarity embedding across other configurations...")
  umap_list <- list()
  pvals <- list()
  for (i in 1:length(cgaudi_loop)) {
    umap_list[[i]] <- vegan::procrustes(reference, cgaudi_loop[[i]]@factors[,-3], symmetric = TRUE)$Yrot
    pvals[[i]] <- vegan::protest(reference, cgaudi_loop[[i]]@factors[,-3], symmetric = TRUE)$signif
  }
  
  non_random_pval <- median(unlist(pvals))
  
  comp1 <- data.frame(id = rownames(reference))
  comp2 <- data.frame(id = rownames(reference))
  
  for (i in 1:length(umap_list)) {
    comp1 <- cbind(comp1, umap_list[[i]][,1])
    comp2 <- cbind(comp2, umap_list[[i]][,2])
  }
  
  comp1 <- comp1[,-1] %>% 
    as.matrix() %>% 
    apply(1, median, na.rm = TRUE)
  
  comp2 <- comp2[,-1] %>% 
    as.matrix() %>% 
    apply(1, median, na.rm = TRUE)
  
  message("Done!")
  
  message("Computing multi-omics clustering...")
  if (is.null(min_pts)) {
    min_pts <- floor(0.03 * nrow(omics[[1]]))
  }
  min_pts <- max(min_pts, 2)
  
  consensus_umap <- data.frame(id = rownames(reference), UMAP1 = comp1, UMAP2 = comp2) %>%
    tibble::remove_rownames() %>% 
    tibble::column_to_rownames("id")
  
  hdbscan_labels <- dbscan::hdbscan(consensus_umap, minPts = min_pts)$cluster
  umap_clusters <- consensus_umap %>% 
    dplyr::mutate(clust = hdbscan_labels)
  
  if (reassign_cluster_zero & length(unique(umap_clusters$clust)) > 1) {
    umap_clusters <- reassign_cluster_zero(umap_clusters, nearest_centroid = FALSE)
  }
  message("Done!")
  
  # Compute silhouette score
  if (length(unique(umap_clusters$clust)) > 1) {
    sil_score <- cluster::silhouette(umap_clusters$clust, dist(umap_clusters[, 1:2]))
    mean_sil_score <- mean(sil_score[, 3])
  } else {
    mean_sil_score <- 0
  }
  
  # Metagenes
  if (compute_features) {
    message("Computing metagenes...")
    if (combine_omics) {
      suppressMessages({
        # Omic type effect correction (batch)
        batch <- generate_batch(unlist(lapply(omics, ncol)))
        features <- as.matrix(dplyr::bind_cols(omics, .name_repair = "unique_quiet"))
        features <- sva::ComBat(dat = features, batch = batch, mod = NULL)
        omics <- c(list(features), omics)
      })
    }
    
    xgboost_fixed_params <- list(objective = "reg:squarederror")
    metagenes <- list()
    for (i in 1:length(omics)) {
      for (j in 1:(ncol(umap_clusters) - 1)) {
        if (method == "xgboost") {
          metagenes_tmp <- xgboost_model(x = omics[[i]], 
                                         y = umap_clusters[,j], 
                                         xgboost_params = c(xgboost_fixed_params, xgboost_params))
        }  else if (method == "rf") {
          metagenes_tmp <- randomForest::importance(
            randomForest::randomForest(x = omics[[i]],
                                       y = umap_clusters[,j])
          )
        }
        metagenes_tmp <- metagenes_tmp %>% 
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
  
  gaudi_res <- new("GAUDIObject",
                   factors = umap_clusters,
                   clusters = umap_clusters$clust,
                   silhouette_score = mean_sil_score,
                   individual_factors = list(),
                   metagenes = if(compute_features) metagenes else list(),
                   gaudiVersion = as.character(packageVersion("gaudi"))
  )
  
  if (validObject(gaudi_res))
    message("Multi-omics integration complete!")
  return(gaudi_res)
}

