
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

umap_factorization <- function(omics, umap_params = list(X = omics)) {
  umap_neighbors <- floor(sqrt(nrow(omics)))
  umap_factors <- do.call(uwot::umap, umap_params)
  umap_factors <- data.frame(umap_factors)
  return(umap_factors)
}

ubmi <- function(omics,
                 umap_params = list(),
                 min_pts = 5,
                 ...) {
                   
  omics <- align_omics(omics) 
                   
  umap_factors <- lapply(omics, function(x) umap_factorization(x, umap_params = list(X = x)))
  concatenated_factors <- list(dplyr::bind_cols(umap_factors, .name_repair = "unique_quiet"))
  
  umap_integrated <- lapply(concatenated_factors, function(x) umap_factorization(x, umap_params = list(X = x)))

  hdbscan_labels <- lapply(umap_integrated, function(x) dbscan::hdbscan(x, minPts = min_pts))[[1]]

  umap_clusters <- data.frame(umap_integrated[[1]], clust = hdbscan_labels$cluster)
  colnames(umap_clusters)[1:2] <- c("UMAP1", "UMAP2")
  
  return(umap_clusters)
  
}

aaa <- ubmi(omics)





