
library(tidyverse)

load("/Users/pol/Dropbox/modp/data/multiomics_train_prism.RData")
load("/Users/pol/Dropbox/modp/data/multiomics_test_prism.RData")

prism <- prism_train %>% 
  dplyr::bind_rows(prism_test) %>% 
  tibble::column_to_rownames("id") %>%
  as.matrix()

mirna <- mirna_train %>% 
  dplyr::bind_rows(mirna_test) %>% 
  tibble::column_to_rownames("id") %>%
  as.matrix()

metabolomics <- metabolomics_train %>%
  dplyr::bind_rows(metabolomics_test) %>% 
  tibble::column_to_rownames("id") %>%
  as.matrix()

methylation <- methylation_train %>% 
  dplyr::bind_rows(methylation_test) %>% 
  tibble::column_to_rownames("id") %>%
  as.matrix()

expression <- expression_train %>% 
  dplyr::bind_rows(expression_test) %>% 
  tibble::column_to_rownames("id") %>%
  as.matrix()

achilles <- achilles_train %>% 
  dplyr::bind_rows(achilles_test) %>% 
  tibble::column_to_rownames("id") %>%
  as.matrix()

colnames(achilles) <- paste0("dep_", colnames(achilles))

##

omics <- list(t(prism), t(achilles), t(expression), t(methylation), t(mirna), t(metabolomics))
ubmi_object <- ubmi(omics, compute_features = FALSE, min_pts = 10,
                    umap_params = list(metric = "cosine"),
                    umap_params_conc = list(metric = "cosine", n_components = 2))

plot_factors(ubmi_object)

ggplot(ubmi_object@factors, aes(UMAP1, UMAP2, color = lineage)) +
  geom_point() +
  theme_bw() +
  gghighlight::gghighlight(lineage == "skin")

find_segments <- function(gene, min_percent = 5) {
  
  spline_num <- 0
  while (spline_num < 3) {
    
    var <- rlang::sym(colnames(gene)[2])
    gene <- gene %>% 
      dplyr::arrange(!!var)
    
    y <- gene[,2]
    x <- 1:length(y)
    
    # cubic spline
    y_spline <- spline(x, y, n = 100, method = "natural")
    
    second_deriv <- diff(diff(y_spline$y)) / ((x[2]-x[1])^2)
    
    # Find the inflection points
    inflection_points <- c()
    for (i in 2:(length(second_deriv) - 1)) {
      # Check if the second derivative changes sign
      if (sign(second_deriv[i]) != sign(second_deriv[i-1]) & sign(second_deriv[i]) != sign(second_deriv[i-1])) {
        inflection_points <- c(inflection_points, x[i])
      }
    }
    
    # Sort the inflection points
    inflection_points <- sort(inflection_points)
    
    min_inf_point <- inflection_points[inflection_points >= min_percent][1]
    max_inf_point <- inflection_points[inflection_points >= (100 - min_percent)][1]
    
    dd <- data.frame(x = 1:100, y = y_spline$y) %>% 
      filter(x == min_inf_point | x == max_inf_point)
    
    dd2 <- data.frame(x, y) %>% 
      dplyr::mutate(tail = case_when(y <= min(dd$y) ~ "lower",
                                     y >= max(dd$y) ~ "upper",
                                     TRUE ~ "middle")) %>% 
      dplyr::mutate(id = gene[,1])
    
    spline_num <- length(table(dd2$tail))
    min_percent <- min_percent + 1
  }
  
  return(dd2)
  
}

molecular_features <- function(drug, 
                               dependencies = achilles,
                               expression = expression,
                               methylation = methylation,
                               miRNA = mirna,
                               metabolomics = metabolomics,
                               drugs = prism,
                               min_percent_cells = 5,
                               fdr_cutoff = 0.05) {
  
  # segments <- find_segments(dependencies[, colnames(dependencies) %in% c("id", gene)], min_percent = min_percent_cells) %>% 
  #   dplyr::filter(tail != "middle") %>% 
  #   dplyr::arrange(id)
  
  segments <- find_segments(drugs[, colnames(drugs) %in% c("id", drug)], min_percent = min_percent_cells) %>% 
    dplyr::filter(tail != "middle") %>% 
    dplyr::arrange(id)
  
  # drugs <- drugs[drugs$id %in% segments$id, -1]
  # methylation <- methylation[methylation$id %in% segments$id, -1]
  # mirna <- miRNA[miRNA$id %in% segments$id, -1]
  # metabolomics <- metabolomics_clean[metabolomics_clean$id %in% segments$id, -1]
  
  # xvar <- cbind(expression, methylation, mirna, metabolomics)
  xvar <- drugs[drugs$id %in% segments$id, -1]
  
  ## Limma
  design <- stats::model.matrix(~ 0 + as.factor(segments$tail))
  colnames(design) <- levels(as.factor(segments$tail))
  cont.matrix <- limma::makeContrasts(contrasts = "upper-lower", levels = design)
  model <- limma::lmFit(t(xvar), design)
  model <- limma::contrasts.fit(model, cont.matrix)
  modelstats <- limma::eBayes(model)
  
  res <- limma::topTable(modelstats, number = ncol(xvar), coef = "upper-lower", sort.by = "p", adjust.method = "fdr") %>%
    tibble::rownames_to_column("feature") %>% 
    filter(adj.P.Val < fdr_cutoff) %>% 
    dplyr::as_tibble()
  
  return(res)
  
}

molecular_features(gene = "ACTG1")

segments <- find_segments(dependencies[, colnames(dependencies) %in% c("id", gene)], min_percent = min_percent_cells) %>% 
  dplyr::filter(tail != "middle") %>% 
  dplyr::arrange(id) %>% 
  left_join(prism, by = "id")

ggplot(segments, aes(tail, asaraldehyde)) +
  geom_boxplot()


for (i in colnames(achilles)[-1]) {
  print(i)
  print(molecular_features(gene = i, fdr_cutoff = 0.05))
}



####
x <- prism[,-1]
y <- achilles$SDHB

mod <- xgboost::xgboost(data = as.matrix(x),
                        label = y,
                        lambda = 0.5, eta = 0.3, gamma = 10, max_depth = 10, subsample = 0.95,
                        nrounds = 10,
                        verbose = FALSE,
                        nthread = parallel::detectCores() - 2,
                        early_stopping_rounds = 8)

shap_values <- SHAPforxgboost::shap.values(xgb_model = mod, X_train = as.matrix(x))

shap_contrib <- shap_values$shap_score
ranked_col <- names(colMeans(abs(shap_contrib))[order(colMeans(abs(shap_contrib)), decreasing = TRUE)])

shap_contrib_df <- data.frame(shap_contrib)
shap_contrib_nonzero <- shap_contrib_df[, apply(shap_contrib_df, 2, function(x) !all(x == 0, na.rm = TRUE))]





