#' Kaplan-Meier Survival Plots (Interactive)
#'
#' This function creates a Kaplan-Meier survival plot and allows the user to
#' select which clusters to plot. 
#'
#' @param object A UBMIObject containing cluster data.
#' @param time A named vector or data frame with a column named "time". 
#' @param censored A named vector or data frame with a column named "censored". 
#' @param select_clusters if TRUE, will open a Shiny window to allow the user to
#' select which clusters they would like to plot
#'
#' @return A Kaplan-Meier plot, of ggplot2 type. 
#' 
#' @export
plot_survival <- function(object, time, censored, select_clusters = FALSE) {
  
  # UBMI Object validation
  if (!inherits(object, "UBMIObject")) {
    stop("The 'object' parameter must be of class 'UBMIObject'.")
  }
  
  # Extract clusters from the object
  clusters <- object@factors[, "clust", drop = FALSE]
  
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
  
  if (!select_clusters) {
    # Filter out clusters with noise and remove NA values
    cleaned_df <- merged_df %>%
      dplyr::filter(clust != 0) %>%
      tidyr::drop_na()
    
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
  
  # Determine the range of y-axis based on all clusters
  full_data <- merged_df %>% dplyr::filter(clust %in% unique(clusters$clust))
  full_model_event <- survival::survfit(survival::Surv(time, censored) ~ clust, data = full_data)
  
  xlim <- range(full_model_event$time)  # Fix x-axis range
  
  interactive_survival <- function(merged_df, selected_clusters, xlim) {
    
    # Filter out clusters with noise and remove NA values
    cleaned_df <- merged_df %>%
      dplyr::filter(clust %in% selected_clusters) %>%
      tidyr::drop_na()
    
    # Convert cluster to factor and sort levels
    cleaned_df$clust <- factor(cleaned_df$clust, levels = sort(unique(cleaned_df$clust)))
    
    # Perform survival analysis
    model_event <- survival::survfit(survival::Surv(time, censored) ~ clust, data = cleaned_df)
    pval <- signif(survminer::surv_pvalue(model_event, data = cleaned_df)$pval, digits = 3)
    
    # Generate Kaplan-Meier plot with fixed plot size
    kaplan_meier_plot <- survminer::ggsurvplot(
      model_event,
      palette = viridis::viridis(length(levels(cleaned_df$clust)), end = 0.8),
      surv.median.line = "hv",
      data = cleaned_df,
      legend.title = "Cluster",
      legend.labs = levels(cleaned_df$clust),
      ylim = c(0, 1),  # Fix y-axis range
      xlim = xlim
    )
    
    # Adjust the legend and plot margins
    kaplan_meier_plot$plot <- kaplan_meier_plot$plot +
      ggplot2::theme(
        legend.position = "right",
        legend.box.margin = margin(0, 20, 0, 0),  # Adjust margin around legend box
        legend.key.size = unit(1, "cm"),           # Adjust size of legend keys
        legend.text = element_text(size = 10),     # Adjust size of legend text
        plot.margin = margin(10, 50, 10, 10)       # Adjust plot margins
      ) +
      ggplot2::labs(
        title = paste("Cluster Distribution (", length(levels(cleaned_df$clust)), " clusters)"),
        subtitle = paste("p-value =", pval)
      )
    
    return(kaplan_meier_plot)
  }
  
  
  # Define UI for the application
  ui <- fluidPage(
    titlePanel("Interactive Kaplan-Meier Plot"),
    
    sidebarLayout(
      sidebarPanel(
        checkboxGroupInput("clusters", "Select Clusters:",
                           choices = NULL, # Choices will be updated in server
                           selected = NULL), 
        actionButton(inputId = "close", label = "Save and Close")
      ),
      
      mainPanel(
        plotOutput("kmPlot")
      )
    )
  )
  
  # Define server logic
  server <- function(input, output, session) {
    
    # Update the cluster choices based on data
    observe({
      clusters <- unique(object@factors$clust)
      clusters <- clusters[clusters != 0]  # Exclude cluster 0
      clusters <- sort(clusters)
      updateCheckboxGroupInput(session, "clusters", choices = clusters, selected = clusters)
    })
    
    output$kmPlot <- renderPlot({
      req(input$clusters)  # Ensure at least one cluster is selected
      
      # Generate the plot based on selected clusters
      plot <- interactive_survival(merged_df, input$clusters, xlim)
      print(plot$plot)
    })
    
    observeEvent(input$close, {
      stopApp(input$clusters)
    })
  }
  
  # Run the application
  app <- shinyApp(ui = ui, server = server)
  result <- runApp(app)
  
  plot <- interactive_survival(merged_df, result, xlim)
  
  return(plot$plot)
}