
min_cells <- 10

y <- sort(achilles_clean$`SDHB (6390)`)
x <- 1:length(y)

# cubic spline
y_spline <- spline(x, y, n = 100, method = "natural")$y

####

second_deriv <- diff(diff(y_spline)) / ((x[2]-x[1])^2)

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

min_inf_point <- inflection_points[inflection_points >= min_cells][1]
max_inf_point <- inflection_points[inflection_points >= (100 - min_cells)][1]

# # Divide the curve into 3 segments based on the inflection points
# segments <- list()
# if (length(inflection_points) == 0) {
#   # No inflection points, divide the curve into 3 equal segments
#   segment_length <- length(x) / 3
#   segments[[1]] <- list(x = x[1:segment_length], y = y[1:segment_length])
#   segments[[2]] <- list(x = x[(segment_length + 1):(2*segment_length)], y = y[(segment_length + 1):(2*segment_length)])
#   segments[[3]] <- list(x = x[(2*segment_length + 1):length(x)], y = y[(2*segment_length + 1):length(x)])
# } else if (length(inflection_points) == 1) {
#   # One inflection point, divide the curve into two segments around the inflection point
#   segment_length <- which(x == inflection_points)
#   segments[[1]] <- list(x = x[1:segment_length], y = y[1:segment_length])
#   segments[[2]] <- list(x = x[(segment_length + 1):length(x)], y = y[(segment_length + 1):length(x)])
# } else {
#   # Two or more inflection points, divide the curve into three segments around the inflection points
#   segments[[1]] <- list(x = x[1:which(x == inflection_points[1])], y = y[1:which(x == inflection_points[1])])
#   segments[[2]] <- list(x = x[(which(x == inflection_points[1]) + 1):which(x == inflection_points[2])], y = y[(which(x == inflection_points[1]) + 1):which(x == inflection_points[2])])
#   segments[[3]] <- list(x = x[(which(x == inflection_points[2]) + 1):length(x)], y = y[(which(x == inflection_points[2]) + 1):length(x)])
# }
# 
# # Plot the segments
# colors <- c("red", "green", "purple")
# for (i in 1:length(segments)) {
#   lines(segments[[i]]$x, segments[[i]]$y, col = colors[i])
# }      


dd <- data.frame(x = 1:100, y_spline) %>% 
  filter(x == min_inf_point | x == max_inf_point)

dd2 <- data.frame(x, y) %>% 
  dplyr::mutate(tail = case_when(y <= min(dd$y_spline) ~ "lower",
                                 y >= max(dd$y_spline) ~ "upper",
                                 TRUE ~ "middle"
                                 ))

ggplot(dd2, aes(x, y, color = tail)) +
  geom_point()

table(dd2$tail)




