
load("/Users/pol/Dropbox/gmv_project/data/22Q2/multiomics_data_processed_all.RData")
colnames(methylation_clean)[2:ncol(methylation_clean)] <- paste0("tss_", colnames(methylation_clean)[2:ncol(methylation_clean)])
colnames(expression_clean)[2:ncol(expression_clean)] <- gsub("\\..*", "", colnames(expression_clean)[2:ncol(expression_clean)])

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

ubmi_object <- ubmi(omics, compute_features = FALSE, min_pts = 4)

plot_factors(ubmi_object)

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

##
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

##
source("/Users/pol/Dropbox/ddh_multiquery/02_algorithm.R")
load("/Users/pol/Dropbox/ddh_multiquery/data/master_final.RData")
pathways <- readRDS(file = "/Users/pol/Dropbox/ddh_multiquery/msigs_ddh_list_large.Rds")

pathways <- pathways %>% 
  dplyr::select(pathway = gs_name, gene = gene_symbol) %>% 
  group_by(pathway) %>%
  dplyr::mutate(count = n()) %>%
  nest(data = gene) %>%
  ungroup()

my_pathways <- tibble(pathway = "custom_pol", 
                      gene = limma_res$feature[1:50]) %>%
  group_by(pathway) %>%
  dplyr::mutate(count = n()) %>%
  nest(data = gene) %>%
  ungroup()

pathways <- bind_rows(pathways, my_pathways)

dep_res_cca <- functional_cca(target = "custom_pol", pathways = pathways, dataset = dependencies, min_genes = 3)

##
text <- sub(".+?_", "", dep_res_cca$Pathway)
text <- gsub("_", " ", text)
cat(text[1:10], sep = ", ")








