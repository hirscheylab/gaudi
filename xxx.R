
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
  
  ubmi_res <- new("UBMIObject",
                  
                  factors = umap_clusters,
                  clusters = umap_clusters$clust,
                  
                  single_factors = umap_factors,
                  
                  metagenes_factor1 = xgboost_metagenesUMAP1[[1]][[1]],
                  metagenes_factor1_rank = xgboost_metagenesUMAP1[[1]][[2]],
                  
                  metagenes_factor2 = xgboost_metagenesUMAP2[[1]][[1]],
                  metagenes_factor2_rank = xgboost_metagenesUMAP2[[1]][[2]],
                  
                  single_metagenes_factor1 = xgboost_metagenesUMAP1[-1],
                  single_metagenes_factor2 = xgboost_metagenesUMAP2[-1]
                  )
  
  if(validObject(ubmi_res))
    return(ubmi_res)
}

test_ubmi <- ubmi(omics, 
                  umap_params = list(n_neighbors = 15, n_components = 4, pca = 50),
                  umap_params_conc = list(n_neighbors = 15, n_components = 2),
                  min_pts = 10,
                  xgboost_params = list())

# saveRDS(test_ubmi, file = "test_ubmi.Rds")

##

library(tidyverse)

plot_metagenes <- function(object,
                           component = 1,
                           top = 10,
                           clusters = FALSE,
                           ...) {
  
  factors <- object@factors
  
  if(component == 1) {
    metagenes <- object@metagenes_factor1
    ranks <- object@metagenes_factor1_rank
  } else {
    metagenes <- object@metagenes_factor2
    ranks <- object@metagenes_factor2_rank
  }
  
  shap_values_nonzero_long <- metagenes %>% 
    mutate(id = rownames(factors), Cluster = as.factor(paste0("Cluster ", factors$clust))) %>% 
    pivot_longer(cols = -c(id, Cluster)) %>%
    mutate(Feature = case_when(name %in% ranks[1:top] ~ name,
                               !(name %in% ranks[1:top]) ~ "Others"))
  
  # colors_raw <- ggsci::pal_npg()(length(unique(shap_values_nonzero_long$Feature)) - 1)
  # names(colors_raw) <- unique(shap_values_nonzero_long$Feature)[unique(shap_values_nonzero_long$Feature) != "Others"]
  # other_color <- "gray80"
  # names(other_color) <- "Others"
  # manual_colors <- c(colors_raw, other_color)
    
  ggplot(shap_values_nonzero_long) +
    {if(!clusters) geom_col(aes(reorder(id, as.numeric(Cluster)), value, fill = Feature))} +
    {if(clusters) geom_col(aes(reorder(id, as.numeric(Cluster)), value, fill = Cluster))} +
    geom_hline(yintercept = 0) +
    labs(x = "Samples (Ranked by Cluster)",
         y = "SHAP Value") +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major = element_blank()) +
    {if(!clusters & component == 1) scale_fill_viridis_d()} +
    {if(!clusters & component != 1) scale_fill_viridis_d(option = "plasma")} +
    {if(clusters) scale_fill_manual(values = ggsci::pal_npg()(length(table(shap_values_nonzero_long$Cluster))))} +
    NULL
}

library(patchwork)

aa <- plot_metagenes(test_ubmi)
bb <- plot_metagenes(test_ubmi, component = 2)

(cc <- plot_metagenes(test_ubmi, clusters = TRUE))
(dd <- plot_metagenes(test_ubmi, clusters = TRUE, component = 2))

##

plot_factors <- function(object,
                         label_size = 0,
                         ...) {
  
  factors <- object@factors %>% 
    mutate(Cluster = as.factor(paste0("Cluster ", clust))) %>% 
    group_by(Cluster) %>% 
    mutate(cord1 = median(UMAP1),
           cord2 = median(UMAP2)) %>% 
    ungroup()

  ggplot(factors, aes(UMAP1, UMAP2)) +
    geom_point(aes(fill = Cluster), pch = 21, size = 3, alpha = 0.8, color = "black") +
    theme_bw() +
    {if(label_size != 0) geom_label(data = factors[!duplicated(factors$Cluster),], 
                                    aes(cord1, cord2, fill = Cluster, label = Cluster),
                                    show.legend = FALSE, size = label_size)} +
    scale_fill_manual(values = ggsci::pal_npg()(length(table(factors$Cluster))))
}

(ee <- plot_factors(test_ubmi, label_size = 4))

##
subtitle0 <- "2-dimensional manifold for 170 samples and 26236 features"
subtitle1 <- "SHAP values of the top features associated to first coordinate (x-axis)"
subtitle2 <- "SHAP values of the top features associated to second coordinate (y-axis)"
  
((ee + labs(subtitle = subtitle0) + theme(legend.position = "bottom")) /
    ((cc + labs(subtitle = subtitle1) + theme(legend.position = "none")) /
       (dd + labs(subtitle = subtitle2) + theme(legend.position = "none")))) | 
  ((aa + labs(subtitle = subtitle1)) / (bb + labs(subtitle = subtitle2)))

# ggsave(filename = "PoC_grid_AML.png", width = 15, height = 9)

##

poma_object <- function(object,
                        omics,
                        ...) {
  
  omics <- align_omics(omics)
  
  object <- POMA::PomaSummarizedExperiment(target = data.frame(id = rownames(object@factors),
                                                               cluster = as.factor(object@clusters)),
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
                           features = test_ubmi@metagenes_factor1_rank[1:10],
                           theme_params = list(axis_x_rotate = TRUE))

c2_exp <- plot_expressions(my_poma, 
                           features = test_ubmi@metagenes_factor2_rank[1:10],
                           theme_params = list(axis_x_rotate = TRUE))

c1_exp/c2_exp

