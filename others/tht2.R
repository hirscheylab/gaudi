
library(tidyverse)

achilles <- vroom::vroom("/Users/pol/Dropbox/gmv_project/data/22Q2/CRISPR_gene_effect.csv", delim = ",")
expression <- vroom::vroom("/Users/pol/Dropbox/gmv_project/data/22Q2/Expression_22Q2_Public.csv", delim = ",")
methylation <- vroom::vroom("/Users/pol/Dropbox/gmv_project/data/22Q2/Methylation_(1kb_upstream_TSS).csv", delim = ",")
mirna <- vroom::vroom("/Users/pol/Dropbox/gmv_project/data/22Q2/miRNA_Expression.csv", delim = ",")
metabolomics <- vroom::vroom("/Users/pol/Dropbox/modp/data/Metabolomics.csv", delim = ",") %>% janitor::clean_names()

cv_compute <- function(data) {
  cvs <- data.frame(feature = colnames(data), cv = apply(data, 2, function(x){sd(x)/mean(x)}))
  cvs <- cvs[order(-cvs$cv),]
  return(cvs)
} 

top_cv_features <- function(data) {
  idx <- PCEmisc::optimal_elbow(data)
  features <- data[1:idx,] 
  return(features)
}

select_features <- function(data) {
  data_proc <- PCEmisc::cleaning(data, knn = TRUE, norm = FALSE, removeZeros = TRUE)
  cvs <- cv_compute(data_proc[,-1])
  top_cvs <- top_cv_features(cvs)
  return(top_cvs$feature)
}

achilles_features <- select_features(achilles)
achilles_features <- sub(" .*", "", achilles_features)
expression_features <- select_features(expression)
methylation_features <- select_features(methylation)
methylation_features <- paste0("tss_", sub("_.*", "", methylation_features))
mirna_features <- select_features(mirna)
metabolomics_features <- select_features(metabolomics)

####
load("/Users/pol/Dropbox/gmv_project/data/22Q2/multiomics_data_processed_all.RData")
colnames(methylation_clean)[2:ncol(methylation_clean)] <- paste0("tss_", colnames(methylation_clean)[2:ncol(methylation_clean)])
colnames(expression_clean)[2:ncol(expression_clean)] <- gsub("\\..*", "", colnames(expression_clean)[2:ncol(expression_clean)])
colnames(achilles_clean)[2:ncol(achilles_clean)] <- sub(" .*", "", colnames(achilles_clean)[2:ncol(achilles_clean)])

partial_pca <- function(data, features) {
  sample_names <- data[,1]
  data <- data[,-1]
  data_top <- data[, colnames(data) %in% features]
  data_rest <- data[, !colnames(data) %in% features]
  
  data_rest_pcs <- data_rest %>% 
    prcomp() %>% 
    purrr::pluck("x") %>% 
    as.data.frame() %>% 
    dplyr::select(PC1:PC50)
  
  data <- data.frame(id = sample_names, data_top, data_rest_pcs)

  return(data)
  
}

# achilles_clean <- partial_pca(achilles_clean, achilles_features)
expression_clean <- partial_pca(expression_clean, expression_features)
methylation_clean <- partial_pca(methylation_clean, methylation_features)
mirna_clean <- partial_pca(mirna_clean, mirna_features)
metabolomics_clean <- partial_pca(metabolomics_clean, metabolomics_features)

####
sample_info <- sample_info %>%
  dplyr::select(DepMap_ID, cell_name, age, lineage, lineage_subtype) %>% 
  dplyr::filter(!(lineage %in% c("blood", "lymphocyte", "plasma_cell", "skin")))
# dplyr::filter(lineage %in% c("thyroid"))

expression <- t(expression_clean %>% filter(id %in% sample_info$DepMap_ID) %>% column_to_rownames("id"))
methylation <- t(methylation_clean %>% filter(id %in% sample_info$DepMap_ID) %>% column_to_rownames("id"))
mirna <- t(mirna_clean %>% filter(id %in% sample_info$DepMap_ID) %>% column_to_rownames("id"))
metabolomics <- t(metabolomics_clean %>% filter(id %in% sample_info$DepMap_ID) %>% column_to_rownames("id"))
# achilles <- t(achilles_clean %>% filter(id %in% sample_info$DepMap_ID) %>% column_to_rownames("id"))

omics <- list(expression, methylation, mirna, metabolomics) #, achilles)

ubmi_object <- ubmi(omics, compute_features = FALSE, min_pts = 10,
                    umap_params = list(metric = "cosine"),
                    umap_params_conc = list(metric = "cosine", n_components = 2))

# ubmi3d <- cbind(ubmi_object@factors, lineage = sample_info$lineage)
# options(warn = -1)
# plot_ly(data = ubmi3d, x = ~UMAP1, y = ~UMAP2, z = ~UMAP3, color = ~lineage) %>% 
#   add_markers()

plot_factors(ubmi_object)
plot_factors(ubmi_object, ad_hoc_label = sample_info$lineage, label_size = 3)

####
ubmi_object_thyroid <- data.frame(ubmi_object@factors, 
                                  lineage = sample_info$lineage, 
                                  lineage_subtype = sample_info$lineage_subtype,
                                  id = rownames(ubmi_object@factors))

ggplot(ubmi_object_thyroid, aes(UMAP1, UMAP2, color = lineage_subtype)) +
  geom_point(size = 3) +
  theme_bw() +
  ggrepel::geom_text_repel(data = ubmi_object_thyroid[ubmi_object_thyroid$lineage == "thyroid",], aes(label = id)) +
  gghighlight::gghighlight(lineage_subtype %in% c("thyroid_carcinoma", "thyroid_squamous"), use_direct_label = FALSE) +
  theme(legend.title = element_blank(),
        legend.position = c(0.2, 0.9))

####
pdb <- ubmi_object_thyroid %>% 
  filter(lineage_subtype == "thyroid_carcinoma") %>% 
  mutate(clust = paste0("clust", clust)) %>% 
  dplyr::select(id, clust) %>% 
  mutate(clust = ifelse(id %in% c("ACH-000174", "ACH-000716"), "clust1", "clust2")) # REMOVE

omicsdata <- bind_cols(expression_clean, methylation_clean[,-1], mirna_clean[,-1], metabolomics_clean[,-1]) %>% 
  filter(id %in% pdb$id) %>%
  left_join(pdb, by = "id") %>%
  column_to_rownames("id") %>% 
  dplyr::select(clust, everything()) %>% 
  as.data.frame()

omicsdata <- smotefamily::SMOTE(X = omicsdata[,-1], target = as.factor(omicsdata[,1]), K = 1)$data

poma_obj <- POMA::PomaSummarizedExperiment(target = data.frame(ID = paste0("sample_", 1:nrow(omicsdata)),
                                                               cluster = omicsdata$class), 
                                           features = omicsdata[,-ncol(omicsdata)])

limma_res <- POMA::PomaLimma(poma_obj, contrast = "clust1-clust2")

POMA::PomaBoxplots(poma_obj, group = "features", feature_name = limma_res$feature[1:20], theme_params = list(axis_x_rotate = TRUE))

####
pathways <- readRDS(file = "/Users/pol/Dropbox/ddh_multiquery/msigs_ddh_list_large.Rds")

MOSEA <- function(data = NULL, # No ID column, just genes
                  groups = NULL,
                  pathways = NULL, # DDH pathways in the long format (63824 entries)
                  adjPval_cutoff = 0.05) {
  
  ddh_pathways <- split(x = pathways$gene_symbol, f = pathways$gs_name)
  
  design <- stats::model.matrix( ~ 0 + groups, data = data)
  ddh_indices <- limma::ids2indices(ddh_pathways, rownames(t(data)))
  gsea_res_dep <- limma::camera(t(data), ddh_indices, design = design) 
  
  gsea_res <- gsea_res_dep %>% 
    dplyr::filter(FDR < adjPval_cutoff) %>% 
    tibble::rownames_to_column("GeneSet")
  
  return(gsea_res)
  
}

aaa <- MOSEA(data = omicsdata[, colnames(omicsdata) %in% colnames(expression_clean[,-1])], 
             groups = omicsdata$class, pathways = pathways, adjPval_cutoff = 1)

bbb <- MOSEA(data = omicsdata[, paste0("tss_", colnames(omicsdata)) %in% colnames(methylation_clean[,-1])], 
             groups = omicsdata$class, pathways = pathways, adjPval_cutoff = 1)


####
load("/Users/pol/Dropbox/gmv_project/data/22Q2/multiomics_data_processed_all.RData")
colnames(methylation_clean)[2:ncol(methylation_clean)] <- paste0("tss_", colnames(methylation_clean)[2:ncol(methylation_clean)])
colnames(expression_clean)[2:ncol(expression_clean)] <- gsub("\\..*", "", colnames(expression_clean)[2:ncol(expression_clean)])
colnames(achilles_clean)[2:ncol(achilles_clean)] <- sub(" .*", "", colnames(achilles_clean)[2:ncol(achilles_clean)])

# sample_info <- sample_info %>%
#   dplyr::select(DepMap_ID, cell_name, age, lineage, lineage_subtype) %>% 
#   dplyr::filter(!(lineage %in% c("blood", "lymphocyte", "plasma_cell", "skin")))
# 
# expression <- t(expression_clean %>% filter(id %in% sample_info$DepMap_ID) %>% column_to_rownames("id"))
# methylation <- t(methylation_clean %>% filter(id %in% sample_info$DepMap_ID) %>% column_to_rownames("id"))
# mirna <- t(mirna_clean %>% filter(id %in% sample_info$DepMap_ID) %>% column_to_rownames("id"))
# metabolomics <- t(metabolomics_clean %>% filter(id %in% sample_info$DepMap_ID) %>% column_to_rownames("id"))
# achilles <- t(achilles_clean %>% filter(id %in% sample_info$DepMap_ID) %>% column_to_rownames("id"))

mf_object <- data.frame(id = rownames(ubmi_object@factors), ubmi_object@factors)

omicsdata <- bind_cols(expression_clean, methylation_clean[,-1], mirna_clean[,-1], metabolomics_clean[,-1]) %>% 
  filter(id %in% mf_object$id)

# find_segments <- function(gene, min_percent = 5) {
#   
#   spline_num <- 0
#   while (spline_num < 3) {
#     
#     var <- rlang::sym(colnames(gene)[2])
#     gene <- gene %>% 
#       dplyr::arrange(!!var)
#     
#     y <- gene[,2]
#     x <- 1:length(y)
#     
#     # cubic spline
#     y_spline <- spline(x, y, n = 100, method = "natural")
#     
#     second_deriv <- diff(diff(y_spline$y)) / ((x[2]-x[1])^2)
#     
#     # Find the inflection points
#     inflection_points <- c()
#     for (i in 2:(length(second_deriv) - 1)) {
#       # Check if the second derivative changes sign
#       if (sign(second_deriv[i]) != sign(second_deriv[i-1]) & sign(second_deriv[i]) != sign(second_deriv[i-1])) {
#         inflection_points <- c(inflection_points, x[i])
#       }
#     }
#     
#     # Sort the inflection points
#     inflection_points <- sort(inflection_points)
#     
#     min_inf_point <- inflection_points[inflection_points >= min_percent][1]
#     max_inf_point <- inflection_points[inflection_points >= (100 - min_percent)][1]
#     
#     dd <- data.frame(x = 1:100, y = y_spline$y) %>% 
#       filter(x == min_inf_point | x == max_inf_point)
#     
#     dd2 <- data.frame(x, y) %>% 
#       dplyr::mutate(tail = case_when(y <= min(dd$y) ~ "lower",
#                                      y >= max(dd$y) ~ "upper",
#                                      TRUE ~ "middle")) %>% 
#       dplyr::mutate(id = gene[,1])
#     
#     spline_num <- length(table(dd2$tail))
#     min_percent <- min_percent + 1
#   }
#   
#   return(dd2)
#   
# }

molecular_features <- function(gene, 
                               cluster = 1,
                               dependencies = achilles_clean,
                               omics = omicsdata,
                               min_percent_cells = 5,
                               fdr_cutoff = 0.05) {
  
  # segments <- find_segments(dependencies[, colnames(dependencies) %in% c("id", gene)], min_percent = min_percent_cells) %>%
  #   dplyr::filter(tail != "middle") %>%
  #   dplyr::arrange(id)
  
  segments <- dependencies[, colnames(dependencies) %in% c("id", gene)] %>%
    dplyr::filter(id %in% mf_object$id[mf_object$clust == cluster]) %>%
    dplyr::arrange(id) %>% 
    rename(gene = 2)
    
  xvar <- omics[omics$id %in% segments$id, -1]
  
  ## Limma
  design <- stats::model.matrix(~ 0 + segments$gene)
  model <- limma::lmFit(t(xvar), design)
  modelstats <- limma::eBayes(model)
  
  res <- limma::topTable(modelstats, number = ncol(xvar), sort.by = "p", adjust.method = "fdr") %>%
    tibble::rownames_to_column("feature") %>% 
    filter(adj.P.Val < fdr_cutoff) %>% 
    dplyr::as_tibble()
  
  return(res)
  
}

molecular_features(gene = "ACTG1")

# 
# ggplot(segments, aes(tail, asaraldehyde)) +
#   geom_boxplot()


for (i in colnames(achilles_clean)[-1]) {
  print(i)
  print(molecular_features(gene = i, fdr_cutoff = 0.05, cluster = 32))
}

test <- mf_object %>% left_join(sample_info, by = c("id" = "DepMap_ID")) %>% filter(clust == 32)
omicsdata_test <- omicsdata[omicsdata$id %in% test$id,]
omicsdata_test <- omicsdata_test %>% 
  left_join(achilles_clean %>% 
              dplyr::select(id, dep_PSMG2 = PSMG2, dep_PPP4R1 = PPP4R1, dep_MYL12B = MYL12B), 
            by = "id")

plot(omicsdata_test$ACTG1, omicsdata_test$dep_PSMG2)












