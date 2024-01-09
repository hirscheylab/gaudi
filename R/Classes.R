#' UBMIObject Class
#'
#' @description
#' A class representing the output of the UMAP-based Multi-omics Integration (UBMI) process.
#' The `UBMIObject` class encapsulates various results and metadata from the UBMI analysis, 
#' including integrated factors, clustering results, silhouette scores, and metagene associations.
#'
#' @slot factors A `data.frame` containing the UMAP factors for the integrated multi-omics data.
#' @slot clusters A numeric vector representing the cluster assignments for each sample.
#' @slot silhouette_score A numeric value representing the average silhouette score of the clustering.
#' @slot individual_factors A list of UMAP factors for individual omics datasets.
#' @slot metagenes A list containing metagenes for individual omics datasets.
#' @slot ubmiVersion A character string indicating the version of the UBMI package used for the analysis.
#'
#' @details
#' The `UBMIObject` class is designed to store comprehensive results from the UBMI analysis, providing an integrated view of multi-omics data.
#' It allows for easy access to both combined and individual omics data results, facilitating further analysis and visualization.
#'
#' @export
setClass("UBMIObject",
         representation(
           factors = "data.frame",
           clusters = "numeric",
           silhouette_score = "numeric",
           individual_factors = "list",
           metagenes = "list",
           ubmiVersion = "character"),
         prototype(
           factors = data.frame(),
           clusters = numeric(),
           silhouette_score = numeric(),
           individual_factors = list(),
           metagenes = list(),
           ubmiVersion = character()))
