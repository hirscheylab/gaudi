#' Enrichment Analysis on UBMI Object
#'
#' This function performs enrichment analysis on UBMI results using the 
#' 'clusterProfiler' package's Gene Set Enrichment Analysis (GSEA). 
#'
#' @param object An UBMIObject. 
#' @param on_omics Numeric index indicating which omics dataset in the `object` to use for plotting. Names of metagenes here should be \strong{gene symbols}. 
#' @param on_factor Numeric index indicating which factor from the multi-omics integration to use for plotting.
#' @param organism Character string specifying the organism database. Default is "org.Hs.eg.db".
#' @param ontology One of "BP", "MF", and "CC" subontologies, or "ALL" for all three
#' @param minGSSize Minimum gene set size
#' @param maxGSSize Maximum gene set size
#'
#' @return A data frame containing the enrichment results. 
#' 
#' @export
ubmi_enrichment <- function(object, 
                            on_omics, 
                            on_factor = 1, 
                            organism = org.Hs.eg.db, 
                            ontology = "ALL", 
                            minGSSize = 2, 
                            maxGSSize = 200) {
  
  # Check if the organism library is installed
  organism_name <- deparse(substitute(organism))
  if (!requireNamespace(organism_name, quietly = TRUE)) {
    stop("Annotation data package '", organism_name, "' is not installed.")
  }
  
  # Extract the gene list based on the specified omics and factor indices
  gene_list <- object@metagenes[[on_omics]][, on_factor]
  names(gene_list) <- rownames(object@metagenes[[on_omics]])
  
  # Remove genes with zero values
  gene_list <- gene_list[gene_list != 0]
  
  # Perform GSEA Analysis
  gse <- clusterProfiler::gseGO(geneList = gene_list, 
                                ont = ontology, 
                                OrgDb = organism, 
                                keyType = "SYMBOL", 
                                minGSSize = minGSSize, 
                                maxGSSize = maxGSSize, 
                                pvalueCutoff = 1, 
                                pAdjustMethod = "none" 
                                ) 
  
  return(gse@result)
}
