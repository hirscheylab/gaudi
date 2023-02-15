
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

##

umap_factorization <- function(umap_params = list()) {
  umap_factors <- do.call(uwot::umap, umap_params)
  umap_factors <- data.frame(umap_factors)
  
  return(umap_factors)
}

##

xgboost_model <- function(x, y, xgboost_params = list()) {
  mod <- xgboost::xgboost(data = x,
                          label = as.matrix(y),
                          params = xgboost_params, 
                          nrounds = 10,
                          verbose = FALSE, 
                          nthread = parallel::detectCores() - 2,
                          early_stopping_rounds = 8)
  
  shap_values <- SHAPforxgboost::shap.values(xgb_model = mod, X_train = x)

  shap_contrib <- shap_values$shap_score
  ranked_col <- names(colMeans(abs(shap_contrib))[order(colMeans(abs(shap_contrib)), decreasing = TRUE)])
  
  shap_contrib_df <- data.frame(shap_contrib)
  shap_contrib_nonzero <- shap_contrib_df[, apply(shap_contrib_df, 2, function(x) !all(x == 0, na.rm = TRUE))]
  
  return(list(shap_contrib_nonzero, ranked_col))
}

##

ubmi <- function(omics,
                 umap_params = list(),
                 umap_params_conc = list(),
                 min_pts = 5,
                 xgboost_params = list(),
                 ...) {
  
  omics <- align_omics(omics)
  
  message("Computing single omic factorizations...")
  umap_factors <- lapply(omics, function(x) umap_factorization(umap_params = c(list(X = x), umap_params)))
  concatenated_factors <- list(dplyr::bind_cols(umap_factors, .name_repair = "unique_quiet"))
  message("Single omic factorizations done!")
  
  message("Computing multi-omics factorization...")
  umap_integrated <- lapply(concatenated_factors, function(x) umap_factorization(umap_params = c(list(X = x), umap_params_conc)))
  message("Multi-omics factorization done!")
  
  message("Computing clusters...")
  hdbscan_labels <- lapply(umap_integrated, function(x) dbscan::hdbscan(x, minPts = min_pts))[[1]]
  
  umap_clusters <- data.frame(umap_integrated[[1]], clust = hdbscan_labels$cluster)
  colnames(umap_clusters)[1:(ncol(umap_clusters) - 1)] <- paste0("UMAP", 1:(ncol(umap_clusters) - 1))
  message("Clusters done!")
  
  message("Computing metagenes...")
  features <- as.matrix(dplyr::bind_cols(omics, .name_repair = "unique_quiet"))
  omics <- c(list(features), omics)
    
  xgboost_fixed_params <- list(objective = "reg:squarederror", eta = 0.02, max_depth = 10, gamma = 0.01, subsample = 0.95)
  
  xgboost_metagenesUMAP1 <- lapply(omics, function(x) xgboost_model(x, y = umap_clusters$UMAP1, xgboost_params = c(xgboost_fixed_params, xgboost_params)))
  message("Factor 1 metagenes done!")
  xgboost_metagenesUMAP2 <- lapply(omics, function(x) xgboost_model(x, y = umap_clusters$UMAP2, xgboost_params = c(xgboost_fixed_params, xgboost_params)))
  message("Factor 2 metagenes done!")
  
  return(list(factorizations = list(umap_clusters, umap_factors),
              metagenes = list(xgboost_metagenesUMAP1, xgboost_metagenesUMAP2)))
}

test_ubmi <- ubmi(omics, 
                  umap_params = list(n_neighbors = 15, n_components = 4, pca = 50),
                  umap_params_conc = list(n_neighbors = 15, n_components = 2),
                  min_pts = 5,
                  xgboost_params = list())

test_ubmi$factorizations[[1]] #factors
test_ubmi$metagenes[[1]][[1]][[1]] #metagens
test_ubmi$metagenes[[1]][[1]][[2]] #ranks







