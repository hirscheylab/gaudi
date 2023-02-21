
library(tidyverse)
library(patchwork)

# DepMap
load("/Users/pol/Dropbox/gmv_project/data/22Q2/multiomics_data_processed_all.RData")
colnames(methylation_clean)[2:ncol(methylation_clean)] <- paste0("tss_", colnames(methylation_clean)[2:ncol(methylation_clean)])
colnames(expression_clean)[2:ncol(expression_clean)] <- gsub("\\..*", "", colnames(expression_clean)[2:ncol(expression_clean)])

# expression <- t(expression_clean %>% column_to_rownames("id"))
# methylation <- t(methylation_clean %>% column_to_rownames("id"))
# mirna <- t(mirna_clean %>% column_to_rownames("id"))
# metabolomics <- t(metabolomics_clean %>% column_to_rownames("id"))

expression <- t(expression_clean %>% filter(id %in% yvar_all$id) %>% column_to_rownames("id"))
methylation <- t(methylation_clean %>% filter(id %in% yvar_all$id) %>% column_to_rownames("id"))
mirna <- t(mirna_clean %>% filter(id %in% yvar_all$id) %>% column_to_rownames("id"))
metabolomics <- t(metabolomics_clean %>% filter(id %in% yvar_all$id) %>% column_to_rownames("id"))

omics <- list(expression, methylation, mirna, metabolomics)

depmap_ubmi <- ubmi(omics, 
                    min_pts = 8,
                    xgboost_params = list(lambda = 1, eta = 0.3, gamma = 50, max_depth = 10, subsample = 0.95))
# saveRDS(depmap_ubmi, file = "depmap_ubmi_no_blood_skin.Rds")
ubmi_object <- readRDS("depmap_ubmi_no_blood_skin.Rds")

clean_object <- drop_clusters(ubmi_object, clusters = c(0)) # c(0, 7:21)
plot_ubmi_grid(clean_object, cluster_label_size = 3, ad_hoc_label = yvar_all$group[yvar_all$id %in% rownames(clean_object@factors)])
plot_ubmi_grid(clean_object, cluster_label_size = 3, ad_hoc_label = yvar_all$lineage[yvar_all$id %in% rownames(clean_object@factors)])
plot_ubmi_grid(clean_object, cluster_label_size = 0)

# TARGET
load("/Users/pol/Dropbox/aml_clusters/data/target_multiomics.RData")

expression <- t(expression_final %>% column_to_rownames("id"))
methylation <- t(methylation_final %>% column_to_rownames("id"))
mirna <- t(mirna_final %>% column_to_rownames("id"))

omics <- list(expression, methylation, mirna)

target_ubmi <- ubmi(omics, 
                    umap_params = list(n_neighbors = 11, 
                                       n_components = 10, pca = 50, min_dist = 0.01),
                    umap_params_conc = list(n_neighbors = 11, 
                                            n_components = 2, spread = 1, min_dist = 0.1),
                    min_pts = 5)
# saveRDS(target_ubmi, file = "target_ubmi.Rds")
# target_ubmi <- readRDS("target_ubmi.Rds")

clean_object <- drop_clusters(target_ubmi, clusters = c(0))
plot_ubmi_grid(clean_object, cluster_label_size = 3)






# ggsave(filename = "PoC_grid_AML.png", width = 15, height = 9)

