
ubmi <- function(omics,
                 umap_params = list(n_neighbors = 15, n_components = 4, pca = 50),
                 umap_params_conc = list(n_neighbors = 15, n_components = 2),
                 min_pts = 10,
                 xgboost_params = list(lambda = 1, eta = 0.3, gamma = 100, max_depth = 10, subsample = 0.95),
                 compute_features = TRUE,
                 ...) {
  
  omics <- align_omics(omics)
  omics <- clean_feature_names(omics)
  
  # Factorization
  umap_factors <- lapply(omics, function(x) umap_factorization(umap_params = c(list(X = x), umap_params)))
  concatenated_factors <- list(dplyr::bind_cols(umap_factors, .name_repair = "unique_quiet"))
  message("Computing individual factorizations... OK!")
  
  umap_integrated <- lapply(concatenated_factors, function(x) umap_factorization(umap_params = c(list(X = x), umap_params_conc)))
  message("Computing multi-omics factorization... OK!")
  
  # Clustering
  hdbscan_labels <- lapply(umap_integrated, function(x) dbscan::hdbscan(x, minPts = min_pts))[[1]]
  
  umap_clusters <- data.frame(umap_integrated[[1]], clust = hdbscan_labels$cluster)
  colnames(umap_clusters)[1:(ncol(umap_clusters) - 1)] <- paste0("UMAP", 1:(ncol(umap_clusters) - 1))
  message("Computing multi-omics clustering... OK!")
  
  if (compute_features) {
    # Metagenes
    message("Computing metagenes...")
    features <- as.matrix(dplyr::bind_cols(omics, .name_repair = "unique_quiet"))
    omics <- c(list(features), omics)
    
    xgboost_fixed_params <- list(objective = "reg:squarederror")
    
    xgboost_metagenesUMAP1 <- lapply(omics, function(x) xgboost_model(x, y = umap_clusters$UMAP1, xgboost_params = c(xgboost_fixed_params, xgboost_params)))
    message("Metagenes associated with the 1st dimension of the manifold... OK!")
    xgboost_metagenesUMAP2 <- lapply(omics, function(x) xgboost_model(x, y = umap_clusters$UMAP2, xgboost_params = c(xgboost_fixed_params, xgboost_params)))
    message("Metagenes associated with the 2nd dimension of the manifold... OK!")
    
    shaps1 <- data.frame(xgboost_metagenesUMAP1[[1]][[1]])
    if (ncol(shaps1) == 1) colnames(shaps1) <- xgboost_metagenesUMAP1[[1]][[2]][1]
    shaps2 <- data.frame(xgboost_metagenesUMAP2[[1]][[1]])
    if (ncol(shaps2) == 1) colnames(shaps2) <- xgboost_metagenesUMAP2[[1]][[2]][1]
  }
  
  ubmi_res <- new("UBMIObject",
                  factors = umap_clusters,
                  clusters = umap_clusters$clust,
                  single_factors = umap_factors,
                  metagenes_factor1 = if(compute_features) shaps1 else data.frame(),
                  metagenes_factor1_rank = if(compute_features) xgboost_metagenesUMAP1[[1]][[2]] else character(),
                  metagenes_factor2 = if(compute_features) shaps2 else data.frame(),      
                  metagenes_factor2_rank = if(compute_features) xgboost_metagenesUMAP2[[1]][[2]] else character(),        
                  single_metagenes_factor1 = if(compute_features) xgboost_metagenesUMAP1[-1] else list(),       
                  single_metagenes_factor2 = if(compute_features) xgboost_metagenesUMAP2[-1] else list(),
                  ubmiVersion = as.character(packageVersion("ubmi"))
  )
  
  if(validObject(ubmi_res))
    message("Integration complete!")
    return(ubmi_res)
}

