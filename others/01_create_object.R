
library(tidyverse)

load("/Users/pol/Dropbox/aml_clusters/data/target_multiomics.RData")

omics <- list(expression_final %>% tibble::column_to_rownames("id"), 
              methylation_final %>% tibble::column_to_rownames("id"), 
              mirna_final %>% tibble::column_to_rownames("id"))

# benchmarking_df <- data.frame(minpts = c(5, 5, 5, 6, 6, 6, 7, 7, 7, 8, 8, 8, 9, 9, 9, 10, 10, 10),
#                               pca = c(100, 75, 50, 100, 75, 50, 100, 75, 50, 100, 75, 50, 100, 75, 50, 100, 75, 50),
#                               mean_sil_score = NA)
# 
# for (i in 1:nrow(benchmarking_df)) {
#   ubmi_object <- ubmi(omics, 
#                       samples_in_rows = TRUE,
#                       compute_features = FALSE,
#                       umap_params = list(n_neighbors = 15, n_components = 10, pca = benchmarking_df$pca[i]),
#                       min_pts = benchmarking_df$minpts[i])
#   
#   sil_score <- cluster::silhouette(ubmi_object@factors$clust, dist(ubmi_object@factors[, 1:2]))
#   mean_sil_score <- mean(sil_score[, 3])
#   benchmarking_df$mean_sil_score[i] <- mean_sil_score
#   
#   print(i)
# }

# saveRDS(benchmarking_df, file = "others/benchmarking_df.Rds")

# Score > 0.7: A strong structure has been found.
# Score between 0.5 and 0.7: A reasonable structure has been found.
# Score < 0.25: The structure is weak and could be artificial. 

# ANOVA between clusters for each omic

sil_value <- 0
for (i in 1:6) {
  ubmi_object <- ubmi(omics, 
                      umap_params = list(n_neighbors = 8, n_components = 4, pca = 75, spread = 0.5), 
                      min_pts = 7, # 131*0.05
                      compute_features = TRUE
                      )
  
  sil_value_tmp <- ubmi_object@silhouette_score
  
  if (sil_value_tmp > sil_value) {
    saveRDS(ubmi_object, file = "others/optimized_ubmi.Rds")
    sil_value <- sil_value_tmp
  }
  
  print(paste0(i, ": ", sil_value_tmp))
}
  
ubmi_object <- readRDS(file = "others/optimized_ubmi.Rds")
plot_factors(ubmi_object)
# plot_metagenes(ubmi_object)

  