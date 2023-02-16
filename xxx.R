
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
  
  # Factorization
  umap_factors <- lapply(omics, function(x) umap_factorization(umap_params = c(list(X = x), umap_params)))
  concatenated_factors <- list(dplyr::bind_cols(umap_factors, .name_repair = "unique_quiet"))
  message("Computing single omics factorizations... OK!")
  
  umap_integrated <- lapply(concatenated_factors, function(x) umap_factorization(umap_params = c(list(X = x), umap_params_conc)))
  message("Computing multi-omics factorization... OK!")
  
  # Clustering
  hdbscan_labels <- lapply(umap_integrated, function(x) dbscan::hdbscan(x, minPts = min_pts))[[1]]
  
  umap_clusters <- data.frame(umap_integrated[[1]], clust = hdbscan_labels$cluster)
  colnames(umap_clusters)[1:(ncol(umap_clusters) - 1)] <- paste0("UMAP", 1:(ncol(umap_clusters) - 1))
  message("Computing multi-omics clusters... OK!")
  
  # Metagenes
  message("Computing metagenes...")
  features <- as.matrix(dplyr::bind_cols(omics, .name_repair = "unique_quiet"))
  omics <- c(list(features), omics)
    
  xgboost_fixed_params <- list(objective = "reg:squarederror", lambda = 1, eta = 0.3, gamma = 100, max_depth = 10, subsample = 0.95)
  
  xgboost_metagenesUMAP1 <- lapply(omics, function(x) xgboost_model(x, y = umap_clusters$UMAP1, xgboost_params = c(xgboost_fixed_params, xgboost_params)))
  message("Metagenes associated with factor 1... OK!")
  xgboost_metagenesUMAP2 <- lapply(omics, function(x) xgboost_model(x, y = umap_clusters$UMAP2, xgboost_params = c(xgboost_fixed_params, xgboost_params)))
  message("Metagenes associated with factor 2... OK!")
  
  return(list(factorizations = list(umap_clusters, umap_factors),
              metagenes = list(xgboost_metagenesUMAP1, xgboost_metagenesUMAP2)))
}

test_ubmi <- ubmi(omics, 
                  umap_params = list(n_neighbors = 15, n_components = 4, pca = 50),
                  umap_params_conc = list(n_neighbors = 15, n_components = 2),
                  min_pts = 10,
                  xgboost_params = list())

##

library(tidyverse)

plot_metagenes <- function(object,
                           component = 1,
                           top = 10,
                           clusters = FALSE,
                           ...) {
  
  factors <- object$factorizations[[1]]
  
  if(component == 1) {
    metagenes <- object$metagenes[[1]][[1]][[1]]
    ranks <- object$metagenes[[1]][[1]][[2]]
  } else {
    metagenes <- object$metagenes[[2]][[1]][[1]]
    ranks <- object$metagenes[[2]][[1]][[2]] 
  }
  
  shap_values_nonzero_long <- metagenes %>% 
    mutate(id = rownames(factors), clust = as.factor(factors$clust)) %>% 
    pivot_longer(cols = -c(id, clust)) %>%
    mutate(feature = case_when(name %in% ranks[1:top] ~ name,
                               !(name %in% ranks[1:top]) ~ "other"))
  
  # colors_raw <- ggsci::pal_npg()(length(unique(shap_values_nonzero_long$feature)) - 1)
  # names(colors_raw) <- unique(shap_values_nonzero_long$feature)[unique(shap_values_nonzero_long$feature) != "other"]
  # other_color <- "gray80"
  # names(other_color) <- "other"
  # manual_colors <- c(colors_raw, other_color)
    
  ggplot(shap_values_nonzero_long) +
    {if(!clusters) geom_col(aes(reorder(id, as.numeric(clust)), value, fill = feature))} +
    {if(clusters) geom_col(aes(reorder(id, as.numeric(clust)), value, fill = clust))} +
    geom_hline(yintercept = 0) +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major = element_blank()) +
    # {if(!clusters) scale_fill_manual(values = manual_colors)} +
    {if(clusters) scale_fill_manual(values = ggsci::pal_npg()(length(table(shap_values_nonzero_long$clust))))} +
    NULL
}

library(patchwork)

aa <- plot_metagenes(test_ubmi)
bb <- plot_metagenes(test_ubmi, component = 2)

aa / bb

(cc <- plot_metagenes(test_ubmi, clusters = TRUE))

##

plot_factors <- function(object,
                         ...) {
  
  factors <- object$factorizations[[1]]

  ggplot(factors, aes(UMAP1, UMAP2, fill = as.factor(clust))) +
    geom_point(pch = 21, size = 3, alpha = 0.8, color = "black") +
    theme_bw() +
    scale_fill_manual(values = ggsci::pal_npg()(length(table(factors$clust))))
}

(dd <- plot_factors(test_ubmi))

##

(dd/cc) | (aa/bb)

##

poma_object <- function(object,
                        omics,
                        ...) {
  
  omics <- align_omics(omics)
  
  object <- POMA::PomaSummarizedExperiment(target = data.frame(id = rownames(object$factorizations[[1]]),
                                                               cluster = as.factor(object$factorizations[[1]]$clust)),
                                           features = as.matrix(dplyr::bind_cols(omics, .name_repair = "unique_quiet"))) %>% 
    POMA::PomaNorm(method = "log_scaling")
  return(object)
}

my_poma <- poma_object(test_ubmi, omics)

plot_expressions <- function(object,
                             features = NULL,
                             theme_params = list(),
                             ...) { 
  
  object %>% 
    POMA::PomaBoxplots(group = "features", feature_name = features, theme_params = theme_params) 
}

c1_exp <- plot_expressions(my_poma, 
                           features = test_ubmi$metagenes[[1]][[1]][[2]][1:10],
                           theme_params = list(axis_x_rotate = TRUE))

c2_exp <- plot_expressions(my_poma, 
                           features = test_ubmi$metagenes[[2]][[1]][[2]][1:10],
                           theme_params = list(axis_x_rotate = TRUE))

c1_exp/c2_exp







