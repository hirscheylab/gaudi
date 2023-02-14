
library(tidyverse)
library(xgboost)
library(SHAPforxgboost)

out <- readRDS(file = "/Users/pol/Dropbox/compare_clusters/clone/results_target/target_ubmi_factorization.Rds")
factors <- out$factorizations[[1]][[1]]
clusters <- paste0("Cluster ", out$UBMI.clusters)

# load("/Users/pol/Dropbox/compare_clusters/clone/results20221114222238/results_clusters/survival_clusters.RData")
# out <- readRDS(file = "/Users/pol/Dropbox/compare_clusters/clone/results20221114222238/amlresults_out.Rds")
# factors <- out$factorizations[[7]][[1]]
# clusters <- paste0("Cluster ", out_clust_surv$aml$UBMI)

# UBMI Scatterplot
ubmi <- data.frame(factors, clust = clusters) %>%
  # filter(clust != "Cluster 0") %>%
  # filter(clust %in% c("Cluster 6", "Cluster 7")) %>%
  filter(!clust %in% c("Cluster 0", "Cluster 1")) %>%
  rownames_to_column("id")

dataX <- cbind(t(omics[[1]]), t(omics[[2]]), t(omics[[3]]))
# dataX <- cbind(t(omics[[1]]))

param_list <- list(objective = "reg:squarederror",  # For regression
                   eta = 0.02,
                   max_depth = 10,
                   gamma = 0.01,
                   subsample = 0.95)

mod <- xgboost::xgboost(data = dataX[rownames(dataX) %in% ubmi$id, ],
                        label = as.matrix(ubmi[, 2]), # Embedding 1 
                        params = param_list, nrounds = 10,
                        verbose = FALSE, nthread = parallel::detectCores() - 2,
                        early_stopping_rounds = 8)

# To return the SHAP values and ranked features by mean|SHAP|
shap_values <- shap.values(xgb_model = mod, X_train = dataX[rownames(dataX) %in% ubmi$id, ])
top_features <- 10

# The ranked features by mean |SHAP|
shap_contrib <- shap_values$shap_score
ranked_col <- names(colMeans(abs(shap_contrib))[order(colMeans(abs(shap_contrib)), decreasing = TRUE)])

shap_values_df <- data.frame(shap_values$shap_score)
shap_values_nonzero <- data.frame(id = ubmi$id, clust = ubmi$clust,
                                  shap_values_df[, apply(shap_values_df, 2, function(x) !all(x == 0, na.rm = TRUE))])

shap_values_nonzero_long <- shap_values_nonzero %>% 
  pivot_longer(cols = -c(id, clust)) %>%
  mutate(feature = case_when(name %in% ranked_col[1:top_features] ~ name,
                             !(name %in% ranked_col[1:top_features]) ~ "other"))

colors_raw <- ggsci::pal_npg()(length(unique(shap_values_nonzero_long$feature)) - 1)
names(colors_raw) <- unique(shap_values_nonzero_long$feature)[unique(shap_values_nonzero_long$feature) != "other"]
other_color <- "gray80"
names(other_color) <- "other"
manual_colors <- c(colors_raw, other_color)

ggplot(shap_values_nonzero_long, aes(reorder(id, as.numeric(gsub("Cluster ", "", clust))), value, fill = feature)) +
  geom_col() +
  geom_hline(yintercept = 0) +
  # facet_wrap(~ clust, scales = "free") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank()) +
  scale_fill_manual(values = manual_colors)

ggplot(shap_values_nonzero_long, aes(reorder(id, as.numeric(gsub("Cluster ", "", clust))), value, fill = clust)) +
  geom_col() +
  geom_hline(yintercept = 0) +
  # facet_wrap(~ clust, scales = "free") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank()) +
  scale_fill_manual(values = ggsci::pal_npg()(length(unique(shap_values_nonzero_long$clust))))

                     