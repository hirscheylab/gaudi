
align_omics <- function(omics) {
  samples <- colnames(omics[[1]])
  for (j in 1:length(omics)) {
    samples <- intersect(samples, colnames(omics[[j]]))
  }
  for (j in 1:length(omics)) {
    omics[[j]] <- omics[[j]][, samples]
    omics[[j]] <- omics[[j]][, which(apply(omics[[j]], 2, sd) > 0)]
  }
  
  # omics <- lapply(omics, FUN = function(x){t(x)})
  return(omics)
}

clean_feature_names <- function(omics) {
  omics <- lapply(omics, FUN = function(x){
    colnames(x) <- gsub("\\.", "_", colnames(x))
    colnames(x) <- gsub("\\-", "_", colnames(x))
    x <- as.matrix(x)
  })
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

