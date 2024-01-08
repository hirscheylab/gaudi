
#' Plot Metagenes
#'
#' @export
plot_metagenes <- function(object,
                           component = 1,
                           top = 10,
                           cluster_line = FALSE,
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
                                             !(name %in% ranks[1:top]) ~ "Other features")) %>% 
    dplyr::mutate(Cluster = factor(Cluster, levels = c(paste0("Cluster ", min(factors$clust):max(factors$clust)))))
  
  if (top > length(ranks)) top <- length(ranks)
  
  ggplot2::ggplot(shap_values_nonzero_long) +
    ggplot2::geom_col(ggplot2::aes(reorder(id, as.numeric(Cluster)), value, fill = Feature)) +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::labs(x = "Samples (Ranked by Cluster)",
                  y = "SHAP Value") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_blank()) +
    
    {if(length(ranks) >= top & top <= 10 & component == 1) 
      scale_fill_manual(breaks = c(ranks[1:top], "Other features"),
                        values = c(ggsci::pal_npg()(top), "gray90"))} +
    {if(length(ranks) >= top & top <= 10 & component == 2) 
      scale_fill_manual(breaks = c(ranks[1:top], "Other features"),
                        values = c(ggsci::pal_simpsons()(top), "gray90"))} +
    
    {if(length(ranks) >= top & top > 10 & component == 1) 
      scale_fill_manual(breaks = c(ranks[1:top], "Other features"),
                        values = c(viridis::viridis(top), "gray90"))} +
    {if(length(ranks) >= top & top > 10 & component == 2) 
      scale_fill_manual(breaks = c(ranks[1:top], "Other features"),
                        values = c(viridis::viridis(top, option = "plasma"), "gray90"))} +
    
    {if(length(ranks) < top & top > 10 & component == 1) 
      scale_fill_manual(breaks = c(ranks[1:top]),
                        values = c(viridis::viridis(length(unique(shap_values_nonzero_long$Feature)))))} +
    {if(length(ranks) < top & top > 10 & component == 2) 
      scale_fill_manual(breaks = c(ranks[1:top]),
                        values = c(viridis::viridis(length(unique(shap_values_nonzero_long$Feature)), option = "plasma")))} +
    
    {if(cluster_line) ggplot2::geom_vline(xintercept = cumsum(as.numeric(table(factors$clust)))[-length(table(factors$clust))],
                                          linetype = "dashed", color = "gray50", linewidth = 0.5)}
}

#' Plot Metagenes Clusters
#'
#' @export
plot_metagenes_clusters <- function(object,
                                    component = 1,
                                    ...) {
  
  factors <- object@factors
  clust_num <- length(table(factors$clust))
  
  if (component == 1) {
    metagenes <- object@metagenes_factor1
  } else {
    metagenes <- object@metagenes_factor2
  }
  
  shap_values_nonzero_long <- metagenes %>% 
    dplyr::mutate(id = rownames(factors), Cluster = as.factor(paste0("Cluster ", factors$clust))) %>% 
    tidyr::pivot_longer(cols = -c(id, Cluster)) %>% 
    dplyr::group_by(Cluster) %>% 
    dplyr::summarise(median_shap = median(value)) %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(Cluster = factor(Cluster, levels = c(paste0("Cluster ", min(factors$clust):max(factors$clust)))))

  ggplot2::ggplot(shap_values_nonzero_long, ggplot2::aes(as.factor(as.numeric(gsub("Cluster ", "", Cluster))),
                                                         median_shap, fill = Cluster)) +
    ggplot2::geom_col(alpha = 0.8) +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::labs(x = NULL,
                  y = "Median SHAP Value") +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank()) +
    {if(clust_num <= 10) ggplot2::scale_fill_manual(values = ggsci::pal_jco()(clust_num))} +
    {if(clust_num > 10) ggplot2::scale_fill_viridis_d(option = "inferno")} +
    NULL
}
         
#' Plot Factors
#'
#' @export
plot_factors <- function(object,
                         label_size = 0,
                         draw_lines = FALSE,
                         ad_hoc_label = NULL,
                         ...) {
  
  factors <- object@factors
  clust_num <- length(table(factors$clust))
  ad_hoc_label_num <- length(table(ad_hoc_label))
  
  if (is.null(ad_hoc_label)) {
    plot_data <- factors %>% 
      dplyr::mutate(Label = as.factor(paste0("Cluster ", clust))) %>% 
      dplyr::group_by(Label) %>% 
      dplyr::mutate(cord1 = median(UMAP1), cord2 = median(UMAP2)) %>% 
      dplyr::ungroup() %>% 
      dplyr::mutate(Label = factor(Label, levels = c(paste0("Cluster ", min(factors$clust):max(factors$clust)))))
    
  } else {
    plot_data <- factors %>% 
      dplyr::mutate(Label = as.factor(ad_hoc_label)) %>% 
      dplyr::group_by(Label) %>% 
      dplyr::mutate(cord1 = median(UMAP1), cord2 = median(UMAP2)) %>% 
      dplyr::ungroup()
  }
  
  plot_complete <- ggplot2::ggplot(plot_data, ggplot2::aes(UMAP1, UMAP2)) +
    {if(draw_lines) geom_vline(xintercept = median(plot_data$UMAP1), linetype = "dashed", color = "black")} +
    {if(draw_lines) geom_hline(yintercept = median(plot_data$UMAP2), linetype = "dashed", color = "black")} +
    ggplot2::geom_point(ggplot2::aes(fill = Label), pch = 21, size = 3, alpha = 0.8, color = "black") +
    ggplot2::theme_bw() +
    {if(label_size != 0) ggplot2::geom_label(data = plot_data[!duplicated(plot_data$Label),], 
                                             ggplot2::aes(cord1, cord2, fill = Label, label = Label),
                                             color = "white", show.legend = FALSE, size = label_size)}
  
  if (is.null(ad_hoc_label)) {
    plot_complete +
      {if(clust_num <= 10) ggplot2::scale_fill_manual(values = ggsci::pal_jco()(clust_num))} +
      {if(clust_num > 10) ggplot2::scale_fill_viridis_d(option = "inferno")} +
      NULL
  } else {
    plot_complete +
      {if(ad_hoc_label_num <= 10) ggplot2::scale_fill_manual(values = ggsci::pal_jco()(ad_hoc_label_num))} +
      {if(ad_hoc_label_num > 10) ggplot2::scale_fill_viridis_d(option = "inferno")} +
      NULL
  }
  
}

#' Plot grid
#'
#' @export
plot_ubmi_grid <- function(object,
                           # top_features = 10,
                           cluster_label_size = 4,
                           ad_hoc_label = NULL,
                           draw_lines = TRUE,
                           ...) {
  
  top_features <- 10
  
  # ranks1 <- object@metagenes_factor1_rank[object@metagenes_factor1_rank %in% colnames(object@metagenes_factor1)]
  # ranks2 <- object@metagenes_factor2_rank[object@metagenes_factor2_rank %in% colnames(object@metagenes_factor2)]
  # 
  # ranks_max <- max(length(ranks1), length(ranks2))
  # 
  # if (top_features > ranks_max) top_features <- ranks_max
  
  subtitle_factors <- paste0("2-dimensional manifold (", nrow(object@factors), " samples and " , length(object@metagenes_factor1_rank), " features)")
  
  factors <- plot_factors(object, label_size = cluster_label_size, ad_hoc_label = ad_hoc_label, draw_lines = draw_lines) +
    ggplot2::theme(legend.position = "bottom") + 
    ggplot2::labs(subtitle = subtitle_factors)
  
  metaclusters1 <- plot_metagenes_clusters(object) + 
    ggplot2::theme(legend.position = "none") + 
    ggplot2::labs(subtitle = "SHAP values by cluster (1st dimension)",
                  x = "Cluster") 
  
  metaclusters2 <- plot_metagenes_clusters(object, component = 2) + 
    ggplot2::theme(legend.position = "none") + 
    ggplot2::labs(subtitle = "SHAP values by cluster (2nd dimension)",
                  x = "Cluster") 
  
  metagenes1 <- plot_metagenes(object, component = 1, top = top_features, cluster_line = TRUE) + 
    ggplot2::labs(subtitle = paste0("SHAP values of the top ", top_features, " features associated with the 1st dimension of the manifold")) 
  
  metagenes2 <- plot_metagenes(object, component = 2, top = top_features, cluster_line = TRUE) + 
    ggplot2::labs(subtitle = paste0("SHAP values of the top ", top_features, " features associated with the 2nd dimension of the manifold")) 
  
  (factors / (metaclusters1 | metaclusters2)) | 
    (metagenes1 / metagenes2)
  
}

