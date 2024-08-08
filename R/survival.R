# Helper function to validate and convert vectors to data frames
convert_to_df <- function(vec, name) {
  if (!is.null(names(vec))) {
    return(data.frame(
      names = names(vec),
      !!name := as.numeric(vec),
      stringsAsFactors = FALSE
    ))
  }
  stop(paste(name, "must be a named vector or data frame with a column named", name))
}

# Function to check if inputs are correct and if so merge into one data frame
validate_and_merge_inputs <- function(object, time, censored) {
  
  # UBMI Object validation
  if (!inherits(object, "UBMIObject")) {
    stop("The 'object' parameter must be of class 'UBMIObject'.")
  }
  
  # Extract clusters from the object
  clusters <- object@factors[, "clust", drop = FALSE]
  
  # Input type validation and conversion
  time <- if (is.vector(time)) convert_to_df(time, "time") else {
    if ("time" %in% colnames(time)) time else stop("time data frame should have a column named \"time\"")
  }
  
  censored <- if (is.vector(censored)) convert_to_df(censored, "censored") else {
    if ("censored" %in% colnames(censored)) censored else stop("censored data frame should have a column named \"censored\"")
  }
  
  # Convert row names to columns
  clusters <- clusters %>% tibble::rownames_to_column(var = "names")
  time <- time %>% tibble::rownames_to_column(var = "names")
  censored <- censored %>% tibble::rownames_to_column(var = "names")
  
  time$time <- as.numeric(time$time)
  censored$censored <- as.numeric(censored$censored)
  
  # Merge data frames
  merged_df <- clusters %>%
    dplyr::left_join(time %>% dplyr::select(names, time), by = "names") %>%
    dplyr::left_join(censored %>% dplyr::select(names, censored), by = "names") %>% 
    tibble::column_to_rownames(var = "names")
  
  return(merged_df)
}

# Function to compute centroids of clusters
compute_centroids <- function(df) {
  num_cols <- ncol(df) - 1  # Number of dimensions (excluding 'clust' column)
  centroid_cols <- names(df)[1:num_cols]
  
  centroids <- df %>%
    dplyr::group_by(clust) %>%
    dplyr::summarise(across(all_of(centroid_cols), median, .names = "centroid_{col}"), .groups = 'drop')
  
  return(centroids)
}

# Function to compute the Euclidean distance between two centroids
compute_distance <- function(centroid1, centroid2) {
  distance <- sqrt(sum((centroid1 - centroid2)^2))
  return(distance)
}

#' Kaplan-Meier Survival Plots
#'
#' This function creates a Kaplan-Meier survival plot and allows the user to
#' select which clusters to plot. 
#'
#' @param object A UBMIObject containing cluster data.
#' @param time A named vector or data frame with a column named "time". 
#' @param censored A named vector or data frame with a column named "censored". 
#' @param show_clusters A numeric list of which clusters to plot, by default 
#' plots all the clusters (except noise, cluster 0). 
#'
#' @return A Kaplan-Meier plot, of ggplot2 type. 
#' 
#' @export
plot_survival <- function(object, time, censored, show_clusters = NULL) {
  
  # Join UBMI object with time and censored into data frame
  merged_df <- validate_and_merge_inputs(object, time, censored)
  
  # Filter out clusters with noise and remove NA values
  cleaned_df <- merged_df %>%
    dplyr::filter(clust != 0) %>% 
    tidyr::drop_na()
  
  # Filter for only selected clusters
  if (!is.null(show_clusters)) {
    cleaned_df <- dplyr::filter(cleaned_df, clust %in% show_clusters) 
  }
  
  # Convert cluster to factor and sort levels
  cleaned_df$clust <- factor(cleaned_df$clust, levels = sort(unique(cleaned_df$clust)))
  
  # Perform survival analysis
  model_event <- survival::survfit(survival::Surv(time, censored) ~ clust, data = cleaned_df)
  pval <- signif(survminer::surv_pvalue(model_event, data = cleaned_df)$pval, digits = 3)
  
  # Generate Kaplan-Meier plot
  kaplan_meier_plot <- survminer::ggsurvplot(
    model_event,
    palette = viridis::viridis(length(levels(cleaned_df$clust)), end = 0.8),
    surv.median.line = "hv",
    data = cleaned_df,
    legend.title = "Cluster",
    legend.labs = levels(cleaned_df$clust)
  )
  
  # Add title and subtitle to the plot
  kaplan_meier_plot$plot <- kaplan_meier_plot$plot +
    ggplot2::labs(
      title = paste("Cluster Distribution (", length(levels(cleaned_df$clust)), " clusters)"),
      subtitle = paste("p-value =", pval)
    )
  
  return(kaplan_meier_plot)
}

#' Kaplan-Meier Survival Summary Table
#'
#' This function provides information on pairwise comparisons between clusters
#' of the UBMI object. 
#'
#' @param object A UBMIObject containing cluster data.
#' @param time A named vector or data frame with a column named "time". 
#' @param censored A named vector or data frame with a column named "censored". 
#' 
#' Centroids are calculated using median, quality is formed by a 80-20 split of 
#' log p-value and scaled distance
#'
#' @return A data frame, with p-value of the survival comparison, distance 
#' between clusters in the UBMI manifold, and quality, a non-interpretable 
#' metric to assist in prioritizing clusters to target. 
#' 
#' @export
#' 
#' @examples
#' table <- get_pairwise_survival_data(ubmi_obj, time, censored)
#' plot_survival(ubmi_obj, time, censored, table[["clusters"]][[1]])
#' 
get_pairwise_survival_data <- function(object, time, censored) {
  
  # Join UBMI object with time and censored into data frame
  merged_df <- validate_and_merge_inputs(object, time, censored)
  
  # Initialize variables
  results <- list()
  
  all_clusters <- unique(object@factors$clust)
  all_clusters <- all_clusters[all_clusters != 0]
  
  # Generate all possible pairs of clusters
  all_pairs <- combn(all_clusters, 2, simplify = FALSE) %>% 
    lapply(sort)
  
  # Calculate centroids for each cluster
  centroids <- compute_centroids(object@factors)
  
  # Loop through all pairs and calculate survival analysis and distances
  for (selected_clusters in all_pairs) {
    # Filter out clusters with noise and remove NA values
    cleaned_df <- merged_df %>%
      dplyr::filter(clust != 0) %>%
      dplyr::filter(clust %in% selected_clusters) %>% 
      tidyr::drop_na()
    
    # Convert cluster to factor and sort levels
    cleaned_df$clust <- factor(cleaned_df$clust, levels = sort(unique(cleaned_df$clust)))
    
    # Perform survival analysis
    model_event <- survival::survfit(survival::Surv(time, censored) ~ clust, data = cleaned_df)
    pval <- signif(survminer::surv_pvalue(model_event, data = cleaned_df)$pval, digits = 3)
    
    # Compute distance between the centroids of the selected clusters
    centroid1 <- centroids %>% 
      dplyr::filter(clust == selected_clusters[1]) %>% 
      dplyr::select(starts_with("centroid_"))
    centroid2 <- centroids %>% 
      dplyr::filter(clust == selected_clusters[2]) %>% 
      dplyr::select(starts_with("centroid_"))
    distance <- compute_distance(as.numeric(centroid1), as.numeric(centroid2))
    
    # Store the result in the list
    results[[paste(selected_clusters, collapse = "_")]] <- list(
      clusters = I(selected_clusters), 
      pval = pval,
      distance = distance
    )
  }
  
  # Convert results to a data frame for easier manipulation
  results_df <- lapply(results, function(x) tibble::tibble(clusters = list(x$clusters), pval = x$pval, distance = x$distance)) %>% 
    dplyr::bind_rows(.id = "rowname") %>% 
    tibble::column_to_rownames("rowname")
  
  # Calculate quality metric
  results_df <- results_df %>%
    mutate(
      scaled_distance = (distance - min(distance)) / (max(distance) - min(distance)),
      n_log_p = log10(pval) / log10(min(pval)), # -log p_val divided by maximum -log p
      quality = 0.8 * n_log_p + 0.2 * scaled_distance
    )
  
  results_df <- results_df %>% 
    mutate(
      scaled_distance = NULL, 
      n_log_p = NULL
    )
  
  # Sort results by the quality metric
  sorted_results <- results_df %>%
    dplyr::arrange(desc(quality))
  
  return (sorted_results)
}