
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

umap_factorization <- function(umap_params = list()) {
  umap_factors <- do.call(uwot::umap, umap_params)
  umap_factors <- data.frame(umap_factors)
  
  return(umap_factors)
}

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

