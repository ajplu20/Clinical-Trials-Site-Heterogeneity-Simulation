predict_mortality_from_data <- function(df, min_mortality = 0.20, max_mortality = 0.6, individual_sd = 0.01) {
  sites <- unique(df$site)
  n_sites <- length(sites)
  
  # Determine site-specific mean mortality values
  if (n_sites == 1) {
    mortality_values <- mean(c(min_mortality, max_mortality))
  } else {
    mortality_values <- seq(min_mortality, max_mortality, length.out = n_sites)
  }
  
  # Create site -> mean mortality mapping
  mortality_map <- setNames(mortality_values, sites)
  
  # Generate random mortality probability per row
  means <- mortality_map[as.character(df$site)]
  predictions <- rnorm(n = nrow(df), mean = means, sd = individual_sd)
  
  # Clip values to [0, 1]
  predictions <- pmin(pmax(predictions, 0), 1)
  
  return(predictions)
}

# Example dataframe
df <- data.frame(site =c("A","B"))

# Apply function
df$mortality_probability <- predict_mortality_from_data(df)
sd(df$mortality_probability)




#  # Create a small dummy dataframe
#  df_test <- data.frame(
#    cal_agey = c(67, 54, 81),
#    sex = c("F", "M", "F"),
#    vap = c(0, 1, 1),
#    hpd_admreason = c("GIT", "GIT", "OTH"),
#    comorbidities_CCI = c(1, 2, 0),
#    site = c("BD 001", "BD 001", "HK 000")
#  )

# Run prediction
#  preds <- predict_mortality_from_data(df_test)

# Show result
#  print(preds)


