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
#' @slot single_factors A list of UMAP factors for individual omics datasets.
#' @slot metagenes_factor1 A `data.frame` containing the metagenes associated with the first dimension of the UMAP manifold.
#' @slot metagenes_factor1_rank A character vector indicating the ranking of metagenes associated with the first dimension.
#' @slot metagenes_factor2 A `data.frame` containing the metagenes associated with the second dimension of the UMAP manifold.
#' @slot metagenes_factor2_rank A character vector indicating the ranking of metagenes associated with the second dimension.
#' @slot single_metagenes_factor1 A list containing metagenes for individual omics datasets associated with the first dimension.
#' @slot single_metagenes_factor2 A list containing metagenes for individual omics datasets associated with the second dimension.
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
           single_factors = "list",
           metagenes_factor1 = "data.frame",
           metagenes_factor1_rank = "character",
           metagenes_factor2 = "data.frame",
           metagenes_factor2_rank = "character",
           single_metagenes_factor1 = "list",
           single_metagenes_factor2 = "list",
           ubmiVersion = "character"),
         prototype(
           factors = data.frame(),
           clusters = numeric(),
           silhouette_score = numeric(),
           single_factors = list(),
           metagenes_factor1 = data.frame(),
           metagenes_factor1_rank = character(),
           metagenes_factor2 = data.frame(),
           metagenes_factor2_rank = character(),
           single_metagenes_factor1 = list(),
           single_metagenes_factor2 = list(),
           ubmiVersion = character()))
