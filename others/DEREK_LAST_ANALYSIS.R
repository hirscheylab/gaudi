
library(tidyverse)

load("/Users/pol/Dropbox/gmv_project/data/22Q2/exp_mth_mir_norm_TRUE.RData")

sample_info_aml <- sample_info %>%
  dplyr::select(DepMap_ID, lineage_subtype) %>%
  dplyr::filter(lineage_subtype == "AML") %>%
  dplyr::pull(DepMap_ID)

omics <- list(expression_clean %>% dplyr::filter(id %in% sample_info_aml) %>% tibble::column_to_rownames("id"), 
              methylation_clean %>% dplyr::filter(id %in% sample_info_aml) %>% tibble::column_to_rownames("id"), 
              mirna_clean %>% dplyr::filter(id %in% sample_info_aml) %>% tibble::column_to_rownames("id"))

colnames(omics[[2]]) <- paste0("TSS_", colnames(omics[[2]]))

ubmi_aml <- ubmi(omics,
                 umap_params = list(n_neighbors = 15, n_components = 4, pca = 30),
                 umap_params_conc = list(n_neighbors = 5, n_components = 2, min_dist = 0.01),
                 min_pts = 5,
                 xgboost_params = list(lambda = 0.1, eta = 0.3, gamma = 10, max_depth = 10, subsample = 0.95),
                 compute_features = TRUE,
                 samples_in_rows = TRUE)

plot_ubmi_grid(ubmi_aml)

# save(ubmi_aml, omics, file = "/Users/pol/Desktop/MANDEL/dz_aml_depmap_multiomics.RData")

load("/Users/pol/Desktop/MANDEL/dz_aml_depmap_multiomics.RData")

features <- dplyr::bind_cols(omics)
metadata <- data.frame(sample = rownames(ubmi_aml@factors), cluster = factor(paste0("cluster_", ubmi_aml@factors$clust)))

poma_object <- POMA::PomaCreateObject(metadata = metadata, features = features)

limma_res <- poma_object %>% 
  POMA::PomaLimma(contrast = "cluster_1-cluster_2")

# poma_object %>% 
#   POMA::PomaVolcano(method = "ttest", 
#                     pval = "adjusted",
#                     pval_cutoff = 0.05,
#                     adjust = "fdr",
#                     log2fc_cutoff = 4,
#                     labels = TRUE)

pathways <- ddh::get_content("universal_pathways", dataset = TRUE)
ddh_pathways <- split(x = pathways$gene_symbol, f = pathways$gs_name)

# GSEA
ordered_genes <- limma_res %>% dplyr::pull(logFC) # NOT AUTO - DRAFT
names(ordered_genes) <- limma_res %>% dplyr::pull(feature) %>% toupper() # NOT AUTO - DRAFT
ordered_genes <- sort(ordered_genes, decreasing = TRUE)

gsea_res_dep <- fgsea::fgsea(pathways = ddh_pathways, stats = ordered_genes, eps = 1e-50,
                             minSize = 2, maxSize = length(ordered_genes) - 1, scoreType = "std")

gsea_res <- gsea_res_dep %>%
  dplyr::rowwise() %>%
  dplyr::mutate(leadingEdge = paste0(leadingEdge, collapse = ", ")) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(Direction = ifelse(NES > 0, "Up", "Down"),
                `Gene Set` = gsub("_.*", "", pathway)) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(pathway = gsub(paste0(`Gene Set`, "_"), "", pathway)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(pathway = stringr::str_replace_all(pathway, "_", " ")) %>%
  dplyr::mutate(pathway = stringr::str_to_title(pathway)) %>%
  dplyr::select(Pathway = pathway, NES, Pval = pval, adjPval = padj,
                Direction, `Pathway Size` = size, Genes = leadingEdge, `Gene Set`) %>%
  dplyr::arrange(Pval)

# gsea_res %>%
#   dplyr::select(Pathway, NES, Pval, adjPval, Direction) %>%
#   dplyr::slice(1:20) %>%
#   dplyr::mutate_if(is.numeric, ~ signif(., digits = 3)) %>%
#   kableExtra::kbl(booktabs = TRUE) %>%
#   kableExtra::kable_styling()

plot_data <- gsea_res %>%
  dplyr::filter(`Gene Set` %in% c("GOBP", "GOMF")) %>% 
  dplyr::group_by(`Gene Set`) %>% 
  dplyr::slice(1:20) %>%
  dplyr::ungroup() %>% 
  dplyr::mutate(Pathway = stringr::str_sub(Pathway, end = 40L))

ggplot(plot_data) +
  # geom_point(aes(-log10(adjPval), reorder(Pathway, -log10(adjPval)), fill = `Gene Set`), pch = 21, alpha = 0.8) +
  # geom_point(aes(-log10(adjPval), reorder(Pathway, -log10(adjPval)), size = NES, fill = `Gene Set`), pch = 21, alpha = 0.8, show.legend = FALSE) +
  geom_col(aes(-log10(adjPval), reorder(Pathway, -log10(adjPval)), fill = `Gene Set`), alpha = 0.9) +
  labs(x = "-log10(FDR)",
       y = NULL) +
  POMA::theme_poma() +
  theme(text = element_text(size = 15),
        legend.position = "bottom") +
  facet_grid(~ Direction) +
  scale_size_continuous(range = c(3, 8)) +
  POMA::scale_fill_poma_d()

aml_factors <- ubmi_aml@factors %>% rownames_to_column("id") %>% left_join(sample_info, by = c("id" = "DepMap_ID")) %>% 
  dplyr::select(id, UMAP1, UMAP2, clust, cell_name, age)

# save(aml_factors, limma_res, gsea_res, file = "/Users/pol/Desktop/MANDEL/dz_aml_depmap_multiomics_results.RData")

load(file = "/Users/pol/Desktop/MANDEL/dz_aml_depmap_multiomics_results.RData")
load("/Users/pol/Dropbox/gmv_project/data/22Q2/multiomics_data_processed_all.RData")

library(POMA)

poma_obj <- PomaCreateObject(metadata = aml_factors %>% 
                               dplyr::filter(id %in% achilles_clean[,1]) %>% 
                               dplyr::mutate(clust = paste0("cluster_", clust)) %>% 
                               dplyr::select(id, clust),
                             features = achilles_clean %>% 
                               dplyr::filter(id %in% aml_factors[,1]) %>% 
                               dplyr::select(-id)
                             )

limma_res2 <- poma_obj %>% 
  PomaLimma(contrast = "cluster_1-cluster_2") %>% 
  dplyr::mutate(feature = gsub("_.*", "", feature))

ordered_genes <- limma_res2 %>% dplyr::pull(logFC) # NOT AUTO - DRAFT
names(ordered_genes) <- limma_res2 %>% dplyr::pull(feature) %>% toupper() # NOT AUTO - DRAFT
ordered_genes <- sort(ordered_genes, decreasing = TRUE)

gsea_res_dep <- fgsea::fgsea(pathways = ddh_pathways, stats = ordered_genes, eps = 1e-50,
                             minSize = 2, maxSize = length(ordered_genes) - 1, scoreType = "std")

gsea_res <- gsea_res_dep %>%
  dplyr::rowwise() %>%
  dplyr::mutate(leadingEdge = paste0(leadingEdge, collapse = ", ")) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(Direction = ifelse(NES > 0, "Up", "Down"),
                `Gene Set` = gsub("_.*", "", pathway)) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(pathway = gsub(paste0(`Gene Set`, "_"), "", pathway)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(pathway = stringr::str_replace_all(pathway, "_", " ")) %>%
  dplyr::mutate(pathway = stringr::str_to_title(pathway)) %>%
  dplyr::select(Pathway = pathway, NES, Pval = pval, adjPval = padj,
                Direction, `Pathway Size` = size, Genes = leadingEdge, `Gene Set`) %>%
  dplyr::arrange(Pval)

# save(limma_res2, gsea_res, file = "/Users/pol/Desktop/yabadabaaaa.RData")

