
# DepMap
load("/Users/pol/Dropbox/gmv_project/data/22Q2/multiomics_data_processed_all.RData")
colnames(methylation_clean)[2:ncol(methylation_clean)] <- paste0("tss_", colnames(methylation_clean)[2:ncol(methylation_clean)])

expression <- t(expression_clean %>% column_to_rownames("id"))
methylation <- t(methylation_clean %>% column_to_rownames("id"))
mirna <- t(mirna_clean %>% column_to_rownames("id"))
metabolomics <- t(metabolomics_clean %>% column_to_rownames("id"))

omics <- list(expression, methylation, mirna, metabolomics)

depmap_ubmi <- ubmi(omics)
# saveRDS(depmap_ubmi, file = "depmap_ubmi.Rds")
# ubmi_object <- readRDS("depmap_ubmi.Rds")

clean_object <- drop_clusters(ubmi_object, clusters = c(0)) # c(0, 7:21)
plot_ubmi_grid(clean_object, cluster_label_size = 3)

# TARGET
load("/Users/pol/Dropbox/aml_clusters/data/target_multiomics.RData")

expression <- t(expression_final %>% column_to_rownames("id"))
methylation <- t(methylation_final %>% column_to_rownames("id"))
mirna <- t(mirna_final %>% column_to_rownames("id"))

omics <- list(expression, methylation, mirna)

target_ubmi <- ubmi(omics)
# saveRDS(target_ubmi, file = "target_ubmi.Rds")
# target_ubmi <- readRDS("target_ubmi.Rds")

clean_object <- drop_clusters(target_ubmi, clusters = c(0))
plot_ubmi_grid(clean_object, cluster_label_size = 3)





# num_factors <- 10
# num_fact_conc <- 2
# pca_comps <- 50
# min_dist_conc <- 0.1
# min_dist_indv <- NULL
# spread_conc <- NULL
# min_pts <- 5

# ggsave(filename = "PoC_grid_AML.png", width = 15, height = 9)

