
poma_object <- function(object,
                        omics,
                        ...) {
  
  omics <- align_omics(omics)
  omics <- clean_feature_names(omics)
  
  object <- POMA::PomaSummarizedExperiment(target = data.frame(id = rownames(object@factors),
                                                               cluster = as.factor(object@clusters)),
                                           features = as.matrix(dplyr::bind_cols(omics, .name_repair = "unique_quiet"))) %>% 
    POMA::PomaNorm(method = "log_scaling")
  
  return(object)
}

my_poma <- poma_object(target_ubmi, omics)

plot_expressions <- function(object,
                             features = NULL,
                             theme_params = list(),
                             ...) { 
  
  object %>% 
    POMA::PomaBoxplots(group = "features", feature_name = features, theme_params = theme_params) 
  
}

c1_exp <- plot_expressions(my_poma, 
                           features = target_ubmi@metagenes_factor1_rank[1:2],
                           theme_params = list(axis_x_rotate = TRUE))

c2_exp <- plot_expressions(my_poma, 
                           features = target_ubmi@metagenes_factor2_rank[1:10],
                           theme_params = list(axis_x_rotate = TRUE))

c1_exp/c2_exp

##

e <- t(SummarizedExperiment::assay(my_poma))
target <- SummarizedExperiment::colData(my_poma) %>% as.data.frame() %>% 
  tibble::rownames_to_column("ID") %>% dplyr::rename(Group = 2) %>% 
  dplyr::select(ID, Group)
data <- cbind(target, e)

  plot_data <- data %>% dplyr::select(-ID) %>% tidyr::pivot_longer(cols = -Group)

    plot_data <- plot_data %>% dplyr::filter(name %in% "CPNE8")

plot_data %>% 
  ggplot2::ggplot(ggplot2::aes(name, value, color = Group, fill = Group)) +
  ggplot2::geom_boxplot(alpha = 0.5) + 
  ggplot2::labs(x = NULL, y = "Value") + 
  theme_bw()






