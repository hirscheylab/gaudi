#' GAUDIObject Class
#'
#' @description
#' A class representing the output of the Group Aggregation via UMAP Data Integration (GAUDI) method.
#' The `GAUDIObject` class encapsulates various results and metadata from the GAUDI analysis, 
#' including integrated factors, clustering results, silhouette scores, and metagene associations.
#'
#' @slot factors A `data.frame` containing the UMAP factors for the integrated multi-omics data.
#' @slot clusters A numeric vector representing the cluster assignments for each sample.
#' @slot silhouette_score A numeric value representing the average silhouette score of the clustering.
#' @slot individual_factors A list of UMAP factors for individual omics datasets.
#' @slot metagenes A list containing metagenes for individual omics datasets.
#' @slot gaudiVersion A character string indicating the version of the GAUDI package used for the analysis.
#'
#' @details
#' The `GAUDIObject` class is designed to store comprehensive results from the GAUDI analysis, providing an integrated view of multi-omics data.
#' It allows for easy access to both combined and individual omics data results, facilitating further analysis and visualization.
#'
#' @export
setClass("GAUDIObject",
         representation(
           factors = "data.frame",
           clusters = "numeric",
           silhouette_score = "numeric",
           individual_factors = "list",
           metagenes = "list",
           gaudiVersion = "character"),
         prototype(
           factors = data.frame(),
           clusters = numeric(),
           silhouette_score = numeric(),
           individual_factors = list(),
           metagenes = list(),
           gaudiVersion = character()))
