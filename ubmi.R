
library(uwot)

folder <- "/Users/pol/Dropbox/compare_clusters/data/cancer/aml/"
file.names = c("log_exp", "methy", "log_mirna")
filtering = "none"
neighbors_conc <- NULL
num_fact_conc <- NULL
min_dist_conc <- NULL
spread_conc <- NULL

omics <- list()
for (i in 1:length(file.names)){
  omics[[i]] <- as.matrix(read.table(paste(folder, file.names[i], sep = "/"), sep = " ", row.names = 1, header = T))
}
##restricting to common samples and filtering
samples <- colnames(omics[[1]])
for (j in 1:length(omics)) {
  samples <- intersect(samples,colnames(omics[[j]]))
}
for (j in 1:length(omics)) {
  omics[[j]] <- omics[[j]][,samples]
  if (filtering != "none") {
    x <- apply(omics[[j]], 1, sd)
    x <- as.matrix(sort(x, decreasing = T))
    w <- which(x > 0)
    if (filtering == "stringent") {
      selected <- rownames(x)[1:min(w[length(w)],5000)]
    } else {
      selected <- rownames(x)[1:min(w[length(w)],6000)]
    }
    m <- match(rownames(omics[[j]]), selected)
    w <- which(!is.na(m))
    omics[[j]] <- omics[[j]][w,]
  } else {
    omics[[j]] <- omics[[j]][, which(apply(omics[[j]], 2, sd) > 0)]
  }
}

####

min_dist_indv <- 0.01
num.factors <- 2
pca_ubmi <- 50

umap_factor <- function(x) {
  umap_neighbors <- floor(sqrt(nrow(t(x))))
  factorization <- uwot::umap(t(x), 
                              n_neighbors = umap_neighbors, 
                              n_components = num.factors, 
                              pca = pca_ubmi,
                              min_dist = min_dist_indv)
  factorization <- as.data.frame(factorization)
  return(factorization)
}

factorizations_umap <- lapply(omics, function(x) umap_factor(x))
concatenated_embeddings <- dplyr::bind_cols(factorizations_umap)

if (is.null(num_fact_conc)) {
  num_fact_conc <- num.factors
}
if (is.null(neighbors_conc)) {
  neighbors_conc <- floor(sqrt(nrow(concatenated_embeddings)))
}
if (is.null(min_dist_conc)) {
  min_dist_conc <- 0.01 # uwot default
}
if (is.null(spread_conc)) {
  spread_conc <- 1 # uwot default
}

concatenated_umap <- uwot::umap(concatenated_embeddings, 
                                n_neighbors = neighbors_conc, 
                                n_components = num_fact_conc,
                                min_dist = min_dist_conc,
                                spread = spread_conc)

####
library(xgboost)
library(SHAPforxgboost)

param_list <- list(objective = "reg:squarederror",  # For regression
                   eta = 0.02,
                   max_depth = 10,
                   gamma = 0.01,
                   subsample = 0.95)

dataX <- t(omics[[1]])

mod <- xgboost::xgboost(data = dataX, 
                        label = as.matrix(concatenated_umap[,1]), 
                        params = param_list, nrounds = 10,
                        verbose = FALSE, nthread = parallel::detectCores() - 2,
                        early_stopping_rounds = 8)

# To return the SHAP values and ranked features by mean|SHAP|
shap_values <- shap.values(xgb_model = mod, X_train = dataX)

# The ranked features by mean |SHAP|
shap_values$mean_shap_score

# To prepare the long-format data:
shap_long <- shap.prep(xgb_model = mod, X_train = dataX)
# is the same as: using given shap_contrib
shap_long <- shap.prep(shap_contrib = shap_values$shap_score, X_train = dataX)
# **SHAP summary plot**
shap.plot.summary(shap_long)

plot_data <- shap.prep.stack.data(shap_contrib = shap_values$shap_score, top_n = 4, n_groups = 6)
# you may choose to zoom in at a location, and set y-axis limit using `y_parent_limit`  
shap.plot.force_plot(plot_data#, 
                     # zoom_in_location = 10,
                     # y_parent_limit = c(-0.1,0.1)
                     )


# if(min_pts == "auto") {
#   min_pts <- floor(0.03 * nrow(t(omics[[1]]))) # 5/170 ~ 0.03 (empirical factor)
# }
# 
# hdbscan_labels <- dbscan::hdbscan(concatenated_umap, minPts = min_pts)
# umap_factors <- data.frame(concatenated_umap)
# 
# umap_metagenes <- function(x) {
#   umap_neighbors <- floor(sqrt(nrow(x)))
#   metagenes <- uwot::umap(x, n_neighbors = umap_neighbors, n_components = num.factors)
#   metagenes <- as.matrix(metagenes)
#   return(metagenes)
# }

# metagenes_umap <- lapply(omics, function(x) umap_metagenes(x))





