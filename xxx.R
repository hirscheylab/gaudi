
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
                 umap_params = list(n_neighbors = 15, n_components = 4, pca = 50),
                 umap_params_conc = list(n_neighbors = 15, n_components = 2),
                 min_pts = 10,
                 xgboost_params = list(lambda = 1, eta = 0.3, gamma = 100, max_depth = 10, subsample = 0.95),
                 ...) {
  
  omics <- align_omics(omics)
  
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
  
  # Metagenes
  message("Computing metagenes...")
  features <- as.matrix(dplyr::bind_cols(omics, .name_repair = "unique_quiet"))
  omics <- c(list(features), omics)
    
  xgboost_fixed_params <- list(objective = "reg:squarederror")
  
  xgboost_metagenesUMAP1 <- lapply(omics, function(x) xgboost_model(x, y = umap_clusters$UMAP1, xgboost_params = c(xgboost_fixed_params, xgboost_params)))
  message("Metagenes associated with the 1st dimension of the mainfold... OK!")
  xgboost_metagenesUMAP2 <- lapply(omics, function(x) xgboost_model(x, y = umap_clusters$UMAP2, xgboost_params = c(xgboost_fixed_params, xgboost_params)))
  message("Metagenes associated with the 2nd dimension of the mainfold... OK!")
  
  ubmi_res <- new("UBMIObject",
                  factors = umap_clusters,
                  clusters = umap_clusters$clust,
                  single_factors = umap_factors,
                  metagenes_factor1 = xgboost_metagenesUMAP1[[1]][[1]],
                  metagenes_factor1_rank = xgboost_metagenesUMAP1[[1]][[2]],
                  metagenes_factor2 = xgboost_metagenesUMAP2[[1]][[1]],
                  metagenes_factor2_rank = xgboost_metagenesUMAP2[[1]][[2]],
                  single_metagenes_factor1 = xgboost_metagenesUMAP1[-1],
                  single_metagenes_factor2 = xgboost_metagenesUMAP2[-1]#,
                  # ubmiVersion = packageVersion("ubmi")
                  )
  
  if(validObject(ubmi_res))
    return(ubmi_res)
}

test_ubmi <- ubmi(omics)

# saveRDS(test_ubmi, file = "test_ubmi.Rds")
test_ubmi <- readRDS("test_ubmi.Rds")

##

library(tidyverse)

plot_metagenes <- function(object,
                           component = 1,
                           top = 10,
                           clusters = FALSE,
                           ...) {
  
  factors <- object@factors
  
  if (component == 1) {
    metagenes <- object@metagenes_factor1
    ranks <- object@metagenes_factor1_rank[object@metagenes_factor1_rank %in% colnames(metagenes)]
  } else {
    metagenes <- object@metagenes_factor2
    ranks <- object@metagenes_factor2_rank[object@metagenes_factor2_rank %in% colnames(metagenes)]
  }
  
  shap_values_nonzero_long <- metagenes %>% 
    dplyr::mutate(id = rownames(factors), Cluster = as.factor(paste0("Cluster ", factors$clust))) %>% 
    tidyr::pivot_longer(cols = -c(id, Cluster)) %>%
    dplyr::mutate(Feature = dplyr::case_when(name %in% ranks[1:top] ~ name,
                                             !(name %in% ranks[1:top]) ~ "Other features"))
    
  ggplot2::ggplot(shap_values_nonzero_long) +
    {if(!clusters) ggplot2::geom_col(ggplot2::aes(reorder(id, as.numeric(Cluster)), value, fill = Feature))} +
    # {if(clusters) ggplot2::geom_col(ggplot2::aes(reorder(id, as.numeric(Cluster)), value, fill = Cluster))} +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::labs(x = "Samples (Ranked by Cluster)",
                  y = "SHAP Value") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = element_blank(),
                   axis.ticks.x = element_blank(),
                   panel.grid.major = element_blank()) +
    
    # {if(clusters) ggplot2::scale_fill_manual(values = ggsci::pal_npg()(length(table(shap_values_nonzero_long$Cluster))))} +
    
    {if(!clusters & length(ranks) >= top & top <= 10 & component == 1) 
      scale_fill_manual(breaks = c(ranks[1:top], "Other features"),
                        values = c(ggsci::pal_npg()(top), "gray90"))} +
    {if(!clusters & length(ranks) >= top & top <= 10 & component == 2) 
      scale_fill_manual(breaks = c(ranks[1:top], "Other features"),
                        values = c(ggsci::pal_simpsons()(top), "gray90"))} +
    
    {if(!clusters & length(ranks) >= top & top > 10 & component == 1) 
      scale_fill_manual(breaks = c(ranks[1:top], "Other features"),
                        values = c(viridis::viridis(top), "gray90"))} +
    {if(!clusters & length(ranks) >= top & top > 10 & component == 2) 
      scale_fill_manual(breaks = c(ranks[1:top], "Other features"),
                        values = c(viridis::viridis(top, option = "plasma"), "gray90"))} +

    {if(!clusters & length(ranks) < top & top > 10 & component == 1) 
      scale_fill_manual(breaks = c(ranks[1:top]),
                        values = c(viridis::viridis(length(unique(shap_values_nonzero_long$Feature)))))} +
    {if(!clusters & length(ranks) < top & top > 10 & component == 2) 
      scale_fill_manual(breaks = c(ranks[1:top]),
                        values = c(viridis::viridis(length(unique(shap_values_nonzero_long$Feature)), option = "plasma")))}
}

library(patchwork)

(aa <- plot_metagenes(test_ubmi, component = 1, top = 1000))
(bb <- plot_metagenes(test_ubmi, component = 2, top = 2000))

(cc <- plot_metagenes(test_ubmi, clusters = TRUE))
(dd <- plot_metagenes(test_ubmi, clusters = TRUE, component = 2))

##

plot_factors <- function(object,
                         label_size = 0,
                         ...) {
  
  factors <- object@factors %>% 
    dplyr::mutate(Cluster = as.factor(paste0("Cluster ", clust))) %>% 
    dplyr::group_by(Cluster) %>% 
    dplyr::mutate(cord1 = median(UMAP1),
                  cord2 = median(UMAP2)) %>% 
    dplyr::ungroup()

  ggplot2::ggplot(factors, ggplot2::aes(UMAP1, UMAP2)) +
    ggplot2::geom_point(ggplot2::aes(fill = Cluster), pch = 21, size = 3, alpha = 0.8, color = "black") +
    ggplot2::theme_bw() +
    {if(label_size != 0) ggplot2::geom_label(data = factors[!duplicated(factors$Cluster),], 
                                             ggplot2::aes(cord1, cord2, fill = Cluster, label = Cluster),
                                             show.legend = FALSE, size = label_size)} +
    ggplot2::scale_fill_manual(values = ggsci::pal_npg()(length(table(factors$Cluster))))
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

