#' Enrichment Analysis on UBMI Object
#'
#' This function performs enrichment analysis on UBMI results using the 
#' 'clusterProfiler' package's Gene Set Enrichment Analysis (GSEA). It provides 
#' options to select which omics data and factor to use, to shorten gene names
#' into 'SYMBOL' format, and to modify the GSEA parameters. 
#'
#' @param object An UBMIObject. 
#' @param on_omics Numeric index indicating which omics dataset in the `object` to use for plotting.
#' @param on_factor Numeric index indicating which factor from the multi-omics integration to use for plotting.
#' @param organism Character string specifying the organism database. Default is "org.Hs.eg.db".
#' @param gseGO_params List of parameters to be passed to the `gseGO` function. Default is a list with parameters:
#' \describe{
#'   \item{ont}{Ontology, one of "BP", "MF", and "CC" subontologies, or "ALL" for all three}
#'   \item{keyType}{Key type of gene}
#'   \item{minGSSize}{Minimum gene set size}
#'   \item{maxGSSize}{Maximum gene set size}
#'   \item{pvalueCutoff}{P-value cutoff}
#'   \item{verbose}{Whether to print messages}
#'   \item{OrgDb}{Organism database, defaults to the value of `organism`}
#'   \item{pAdjustMethod}{P-value adjustment method}
#' }
#' @param shorten_gene_names Logical indicating whether to shorten gene names by removing everything after the dot. 'BRCA2.675' becomes 'BRCA2' for example. 
#'
#' @note This function assigns a value to the global environment using the `<<-` operator.
#'
#' @return A 'gseaResult' object containing the enrichment results. 
#' 
#' @export
ubmi_enrichment <- function(object, 
                            on_omics, 
                            on_factor = 1, 
                            organism = "org.Hs.eg.db", 
                            gseGO_params = list(ont = "ALL",
                                                keyType = "SYMBOL",
                                                minGSSize = 2,
                                                maxGSSize = 200,
                                                pvalueCutoff = 0.05,
                                                verbose = TRUE,
                                                OrgDb = organism,
                                                pAdjustMethod = "none"), 
                            shorten_gene_names = FALSE) {
  
  # Load the organism database dynamically
  library(organism, character.only = TRUE)
  
  # Extract the gene list based on the specified omics and factor indices
  gene_list <- object@metagenes[[on_omics]][, on_factor]
  names(gene_list) <- rownames(object@metagenes[[on_omics]])
  
  # Remove genes with zero values
  gene_list <- gene_list[gene_list != 0]
  
  if (shorten_gene_names) {
    # Shorten gene names by removing everything after the dot
    gene_df <- data.frame(
      gene = gsub("\\..*", "", names(gene_list)),
      value = gene_list,
      stringsAsFactors = FALSE
    )
    
    # Combine duplicates and calculate the median value
    combined_gene_df <- gene_df %>%
      dplyr::group_by(gene) %>%
      dplyr::summarise(value = median(value))
    
    # Convert the data frame back to a named vector
    gene_list <- setNames(combined_gene_df$value, combined_gene_df$gene)
    
    # Sort the gene list in decreasing order
    gene_list <- sort(gene_list, decreasing = TRUE)
  }
  
  # Set default parameters for gseGO
  default_params <- list(
    geneList = gene_list,
    ont = "ALL",
    keyType = "SYMBOL",
    minGSSize = 2,
    maxGSSize = 200,
    pvalueCutoff = 0.05,
    verbose = TRUE,
    OrgDb = organism,
    pAdjustMethod = "none"
  )
  
  # Override default parameters with user-provided parameters
  gseGO_params <- modifyList(default_params, gseGO_params)
  
  # Perform GSEA
  gse <- do.call(clusterProfiler::gseGO, gseGO_params)
  
  return(gse)
}