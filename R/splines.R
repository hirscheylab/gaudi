
find_segments <- function(gene, min_percent = 10) {
  
  spline_num <- 0
  while (spline_num < 3) {
    y <- sort(gene)
    x <- 1:length(y)
    
    # cubic spline
    y_spline <- spline(x, y, n = 100, method = "natural")
    
    second_deriv <- diff(diff(y_spline$y)) / ((x[2]-x[1])^2)
    
    # Find the inflection points
    inflection_points <- c()
    for (i in 2:(length(second_deriv) - 1)) {
      # Check if the second derivative changes sign
      if (sign(second_deriv[i]) != sign(second_deriv[i-1]) & sign(second_deriv[i]) != sign(second_deriv[i-1])) {
        inflection_points <- c(inflection_points, x[i])
      }
    }
    
    # Sort the inflection points
    inflection_points <- sort(inflection_points)
    
    min_inf_point <- inflection_points[inflection_points >= min_percent][1]
    max_inf_point <- inflection_points[inflection_points >= (100 - min_percent)][1]
    
    dd <- data.frame(x = 1:100, y = y_spline$y) %>% 
      filter(x == min_inf_point | x == max_inf_point)
    
    dd2 <- data.frame(x, y) %>% 
      dplyr::mutate(tail = case_when(y <= min(dd$y) ~ "lower",
                                     y >= max(dd$y) ~ "upper",
                                     TRUE ~ "middle")) 
    
    spline_num <- length(table(dd2$tail))
    min_percent <- min_percent + 1
  }
  
  ggplot(dd2, aes(x, y, color = tail)) +
    geom_point() +
    labs(x = paste0("Cell Rank (", table(dd2$tail)[1], " sensitive cells vs ", table(dd2$tail)[3], " resistant cells)"),
         y = "Dependency Score") +
    theme_classic() +
    theme(legend.position = "none") +
    scale_color_manual(values = c("lower" = "red", "middle" = "gray", "upper" = "red"))
  
}

a <- find_segments(gene = achilles_clean$`TP53 (7157)`, min_percent = 5) + labs(title = "TP53")
b <- find_segments(gene = achilles_clean$`SDHB (6390)`, min_percent = 5) + labs(title = "SDHB")
c <- find_segments(gene = achilles_clean$`AASS (10157)`, min_percent = 5) + labs(title = "AASS")
d <- find_segments(gene = achilles_clean$`SIRT4 (23409)`, min_percent = 5) + labs(title = "SIRT4")

library(patchwork)

(a + b) / (c + d)

