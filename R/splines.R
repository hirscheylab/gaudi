
min_cells <- 10 # 10%

y <- sort(achilles_clean$`SDHB (6390)`)
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

min_inf_point <- inflection_points[inflection_points >= min_cells][1]
max_inf_point <- inflection_points[inflection_points >= (100 - min_cells)][1]

dd <- data.frame(x = 1:100, y = y_spline$y) %>% 
  filter(x == min_inf_point | x == max_inf_point)

dd2 <- data.frame(x, y) %>% 
  dplyr::mutate(tail = case_when(y <= min(dd$y) ~ "lower",
                                 y >= max(dd$y) ~ "upper",
                                 TRUE ~ "middle"))

ggplot(dd2, aes(x, y, color = tail)) +
  geom_point()

table(dd2$tail)




