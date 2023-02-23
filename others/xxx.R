
dichotomize_clusters <- function(object,
                                 component = 1,
                                 ...) {
  
  factors <- object@factors
  
  if (component == 1) {
    metagenes <- object@metagenes_factor1
  } else {
    metagenes <- object@metagenes_factor2
  }
  
  groups <- metagenes %>% 
    dplyr::mutate(id = rownames(factors), Cluster = as.factor(paste0("Cluster ", factors$clust))) %>% 
    tidyr::pivot_longer(cols = -c(id, Cluster)) %>% 
    dplyr::group_by(Cluster) %>% 
    dplyr::summarise(median_shap = median(value)) %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(group = dplyr::case_when(median_shap > 0 ~ "high",
                                           median_shap <= 0 ~ "low"))
  
  group <- factors %>% 
    dplyr::mutate(Cluster = as.factor(paste0("Cluster ", clust))) %>% 
    dplyr::left_join(groups, by = "Cluster") %>% 
    dplyr::pull(group)
  
  return(group)
  
}

# dichotomize_clusters(ubmi_object, component = 2)

ubmi_gsea <- function(data,
                      groups = NULL,
                      species = "human", 
                      category = "H",
                      subcategory = NULL,
                      ...) {
  
  msigdb_paths <- msigdbr::msigdbr(species = species, category = category, subcategory = subcategory)
  hallmark_pathways <- split(x = msigdb_paths$gene_symbol, f = msigdb_paths$gs_name)
  
  design <- stats::model.matrix(~ 0 + as.factor(groups))
  colnames(design) <- levels(as.factor(groups))
  
  hallmark_indices_exp <- limma::ids2indices(hallmark_pathways, rownames(data))
  gsea_exp <- limma::camera(data, hallmark_indices_exp, design = design) 
  
  gsea_res <- gsea_exp %>% 
    dplyr::filter(FDR < 0.05) %>% 
    tibble::rownames_to_column()
  
  return(gsea_res)
}

# ubmi_gsea(omics[[1]], groups = dichotomize_clusters(ubmi_object, component = 1))
# ubmi_gsea(omics[[1]], groups = dichotomize_clusters(ubmi_object, component = 2))

####

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

my_poma <- poma_object(ubmi_object, omics)

plot_expressions <- function(object,
                             features = NULL,
                             theme_params = list(),
                             ...) { 
  
  object %>% 
    POMA::PomaBoxplots(group = "features", feature_name = features, theme_params = theme_params)
  
}

c1_exp <- plot_expressions(my_poma, 
                           features = ubmi_object@metagenes_factor1_rank[1],
                           theme_params = list(axis_x_rotate = TRUE))

c2_exp <- plot_expressions(my_poma, 
                           features = ubmi_object@metagenes_factor2_rank[1:2],
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






