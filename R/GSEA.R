# GSEA
library(clusterProfiler)
library(enrichplot)
# need: named vector in decreasing order, gene | fold change
# extracts 
gsea <- function(object, 
                 index, 
                 UMAP_index = 1, 
                 organism = "org.Hs.eg.db", 
                 gseGO_params = list(ont = "ALL",
                                     keyType = "SYMBOL",
                                     minGSSize = 3,
                                     maxGSSize = 800,
                                     pvalueCutoff = 0.05,
                                     verbose = TRUE,
                                     OrgDb = organism,
                                     pAdjustMethod = "none"), 
                 shorten_gene_names = FALSE, 
                 norm = TRUE
) {
  library(organism, character.only = TRUE)
  
  gene_list <- object@metagenes[[index]][, UMAP_index]
  names(gene_list) <- rownames(object@metagenes[[index]])
  
  gene_list <- gene_list[gene_list != 0]
  
  if (shorten_gene_names) {
    # Convert the gene_list to a data frame for easier manipulation
    gene_df <- data.frame(
      gene = gsub("\\..*", "", names(gene_list)),
      value = gene_list,
      stringsAsFactors = FALSE
    )
    
    # Combine duplicates and calculate the median
    combined_gene_df <- gene_df %>%
      group_by(gene) %>%
      summarise(value = median(value))
    
    # Convert the data frame back to a named vector
    gene_list <- setNames(combined_gene_df$value, combined_gene_df$gene)
    
    gene_list <- sort(gene_list, decreasing = TRUE)
  }
  
  gene_list <<- gene_list
  
  if (norm) {
    scaled <- as.vector(scale(gene_list))
    names(scaled) <- names(gene_list)
    gene_list <- scaled
  }
  
  # Set default parameters for gseGO
  default_params <- list(
    geneList = gene_list,
    ont = "ALL",
    keyType = "SYMBOL",
    minGSSize = 3,
    maxGSSize = 800,
    pvalueCutoff = 0.05,
    verbose = TRUE,
    OrgDb = organism,
    pAdjustMethod = "none"
  )
  
  # Merge default parameters with user-provided parameters
  gseGO_params <- modifyList(default_params, gseGO_params)
  
  # Perform GSEA
  gse <- do.call(gseGO, gseGO_params)
  
  gse
}

## Using the function

load("/Users/home/TCGA/Cancer Plots/CUBMI 2024-07-17 15-20-41/ubmi_results")
res <- gsea(sorted_ubmi_results[["gbm"]], 2, 1, shorten_gene_names = TRUE)

library(tidyverse)
dotplot(res, showCategory=5, split=".sign") + 
  facet_grid(.~.sign) + 
  theme(
    axis.text.x = element_text(size = 8),  # Adjust font size for x-axis text
    axis.text.y = element_text(size = 8),  # Adjust font size for y-axis text
    text = element_text(size = 10)  # Adjust font size for other text elements
  )
enrichplot::gseaplot2(res, geneSetID = 1, pvalue_table = TRUE)
heatplot(res, foldChange=gene_list, showCategory=10)
