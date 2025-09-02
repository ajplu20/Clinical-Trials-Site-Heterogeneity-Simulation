
#import functions
source("data generation bootstrapping.R")
source("data generation distribution.R")


library(ggplot2)
library(dplyr)



plot_country_continuous_comparison <- function(synthetic_df, observed_df, site, variable = "cal_agey") {
  country <- substr(site, 1, 2)
  
  observed_subset <- observed_df %>%
    filter(country == !!country)
  
  synthetic_subset <- synthetic_df %>%
    filter(site == !!site)
  
  observed_subset$source <- "Observed (Country)"
  synthetic_subset$source <- "Synthetic (Site)"
  
  combined <- bind_rows(
    observed_subset[, c(variable, "source")],
    synthetic_subset[, c(variable, "source")]
  )
  
  ggplot(combined, aes_string(x = variable, fill = "source")) +
    geom_density(alpha = 0.4, adjust = 1.2) +
    labs(
      title = paste("Distribution of", variable, "for site", site, "vs country", country),
      x = variable,
      y = "Density"
    ) +
    theme_minimal()
}


plot_country_categorical_comparison <- function(synthetic_df, observed_df, site, variable = "sex") {
  country <- substr(site, 1, 2)
  
  observed_subset <- observed_df %>%
    filter(country == !!country)
  
  synthetic_subset <- synthetic_df %>%
    filter(site == !!site)
  
  observed_subset$source <- "Observed (Country)"
  synthetic_subset$source <- "Synthetic (Site)"
  
  combined <- bind_rows(
    observed_subset[, c(variable, "source")],
    synthetic_subset[, c(variable, "source")]
  )
  
  combined %>%
    group_by(source, .data[[variable]]) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(source) %>%
    mutate(percent = count / sum(count)) %>%
    ggplot(aes_string(x = variable, y = "percent", fill = "source")) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(
      title = paste("Proportion of", variable, "for site", site, "vs country", country),
      y = "Proportion",
      x = variable
    ) +
    theme_minimal()
}



sites <- c("CN 001", "VN 002")
site_c <- c(2500,1500)
synthetic_data <- generate_random_dataset_bootstrap(sites, site_c)
synthetic_data <- generate_random_dataset_distribution(sites, site_c)
table(synthetic_data$site)
plot_country_continuous_comparison(synthetic_data, df, site = "CN 001", variable = "cal_agey")

plot_country_categorical_comparison(synthetic_data, df, site = "CN 001", variable = "hpd_admreason")


