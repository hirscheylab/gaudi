
#' Plot GAUDI Factors
#'
#' This function creates a plot of factors derived from GAUDI.
#' It provides options to adjust label sizes, draw lines, and use ad-hoc labels for clustering.
#'
#' @param object A `GAUDIObject`.
#' @param label_size Numeric value specifying the size of labels on the plot. 
#'                   A value of 0 means labels are not drawn.
#' @param draw_lines Logical value indicating whether to draw dashed lines at the median of UMAP coordinates.
#' @param ad_hoc_label Optional vector of labels to be used instead of the default cluster labels.
#' @param palette Character string specifying the color palette for the plot.
#'                Supported palettes include 'inferno', 'plasma', and others from the viridis package.
#'                
#' @return A ggplot object representing the factor plot.
#'
#' @export
plot_factors <- function(object,
                         label_size = 0,
                         draw_lines = FALSE,
                         ad_hoc_label = NULL,
                         palette = "magma") {
  
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
    {if(draw_lines) ggplot2::geom_vline(xintercept = median(plot_data$UMAP1), linetype = "dashed", color = "black")} +
    {if(draw_lines) ggplot2::geom_hline(yintercept = median(plot_data$UMAP2), linetype = "dashed", color = "black")} +
    ggplot2::geom_point(ggplot2::aes(fill = Label), pch = 21, size = 3, alpha = 0.8, color = "black") +
    ggplot2::theme_bw() +
    {if(label_size != 0) ggplot2::geom_label(data = plot_data[!duplicated(plot_data$Label),], 
                                             ggplot2::aes(cord1, cord2, fill = Label, label = Label),
                                             color = "white", show.legend = FALSE, size = label_size)} +
    ggplot2::scale_fill_viridis_d(option = palette, begin = 0.2, end = 0.9) +
    NULL
  
  return(plot_complete)
}

#' Plot GAUDI Metagenes
#'
#' This function creates a plot of the top metagenes based on SHAP values derived from GAUDI's feature importance analysis.
#' It provides options to select which omics data and factor to use, the number of top features to display, and the color palette.
#'
#' @param object A `GAUDIObject`.
#' @param on_omics Numeric index indicating which omics dataset in the `object` to use for plotting.
#' @param on_factor Numeric index indicating which factor from the multi-omics integration to use for plotting.
#' @param top Integer specifying the number of top features to display in the plot.
#' @param palette Character string specifying the color palette for the plot.
#'                Supported palettes include 'inferno', 'plasma', and others from the viridis package.
#'
#' @return A ggplot object representing the metagenes plot with SHAP values.
#'
#' @export
plot_metagenes <- function(object,
                           on_omics = 1,
                           on_factor = 1,
                           top = 10,
                           palette = "plasma") {
  
  factors <- object@factors
  metagenes <- data.frame(shap = object@metagenes[[on_omics]][, on_factor])
  
  shap_values_nonzero_long <- metagenes %>% 
    dplyr::mutate(Feature = rownames(object@metagenes[[on_omics]])) %>% 
    dplyr::arrange(dplyr::desc(abs(shap)))
  
  if (top > length(shap_values_nonzero_long$shap[shap_values_nonzero_long$shap > 0])) {
    top <- length(shap_values_nonzero_long$shap[shap_values_nonzero_long$shap > 0])
  }
  
  shap_values_nonzero_long <- shap_values_nonzero_long %>% 
    dplyr::mutate(Feature = c(Feature[1:top], rep("Other features", dplyr::n() - top)),
                  Feature = factor(Feature, levels = Feature[1:top])) %>% 
    dplyr::slice(1:top)
  
  plot_complete <- ggplot2::ggplot(shap_values_nonzero_long) +
    ggplot2::geom_col(ggplot2::aes(reorder(Feature, -abs(as.numeric(shap))), shap, fill = Feature)) +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::labs(x = "Features (Ranked by SHAP)",
                  y = "SHAP Value") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_blank()) +
    ggplot2::scale_fill_viridis_d(option = palette, begin = 0.2, end = 0.9) +
    NULL
  
  return(plot_complete)
}

#' Plot GAUDI Grid with Factors and Metagenes
#'
#' This function creates a comprehensive grid plot for the GAUDI results. 
#' It combines plots of factors and top metagenes associated with different dimensions of the manifold.
#'
#' @param object A `GAUDIObject`.
#' @param top_features Integer specifying the number of top features to display in the metagenes plots.
#' @param on_omics Numeric index indicating which omics dataset in the `object` to use for plotting.
#' @param label_size Numeric value specifying the size of labels on the factors plot.
#' @param draw_lines Logical value indicating whether to draw dashed lines at the median of UMAP coordinates.
#' @param ad_hoc_label Optional vector of labels to be used instead of the default cluster labels.
#' @param annotations Logical indicating whether to add annotations to the plot.
#'
#' @return A patchwork grid plot combining the factors plot and two metagenes plots.
#'
#' @export
plot_gaudi_grid <- function(object,
                            top_features = 10,
                            on_omics = 1,
                            label_size = 4,
                            draw_lines = TRUE,
                            ad_hoc_label = NULL,
                            annotations = TRUE) {
  
  factors_plot <- plot_factors(object, label_size = max(1, label_size), draw_lines = draw_lines, 
                               ad_hoc_label = ad_hoc_label, palette = "magma") +
    ggplot2::theme(legend.position = "none",
                   legend.title = ggplot2::element_blank()) + 
    ggplot2::labs(subtitle = paste0("2-dimensional manifold (", nrow(object@factors), " samples and " ,
                                    sum(unlist(lapply(object@metagenes, nrow))), " features)"))
  
  metagenes_plot1 <- plot_metagenes(object, top = top_features, on_omics = on_omics, 
                                    on_factor = 1, palette = "plasma") + 
    ggplot2::labs(subtitle = paste0("SHAP values of the top ", top_features, 
                                    " features associated with the 1st dimension of the manifold"))
  
  metagenes_plot2 <- plot_metagenes(object, top = top_features, on_omics = on_omics, 
                                    on_factor = 2, palette = "mako") + 
    ggplot2::labs(subtitle = paste0("SHAP values of the top ", top_features, 
                                    " features associated with the 2nd dimension of the manifold"))
  
  plot_complete <- factors_plot | (metagenes_plot1 / metagenes_plot2)
  
  if (annotations) {
    plot_complete <- plot_complete +
      patchwork::plot_annotation(tag_levels = "A") 
  }
  
  return(plot_complete)
}

