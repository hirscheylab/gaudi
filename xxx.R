
test_ubmi <- readRDS("test_ubmi.Rds")

plot_ubmi_grid(test_ubmi)

  


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

