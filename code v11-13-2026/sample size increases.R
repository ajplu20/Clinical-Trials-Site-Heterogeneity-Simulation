#load in files
source("simulate trial.R")

#####################run simulations
run_main_simulation <- function(
    total_sample_size   = 964,
    max_sample_multiply = 10,
    num_sites           = 5,
    iterations          = 10,
    generation_method   = "distribution",  # or "bootstrap"
    treatment_proportion = 0.5
) {
  # Pre-allocate list for speed
  results_list <- vector("list", max_sample_multiply)
  
  for (mult in seq_len(max_sample_multiply)) {
    current_sample_size <- total_sample_size * mult
    cat("Running simulation with sample size", current_sample_size, "...\n")
    
    trial_results <- simulate_multiple_trials(
      num_sites             = num_sites,
      specification         = "random_selection",
      iterations            = iterations,
      total_sample_size     = current_sample_size,
      site_counts_per_site  = "same",
      treatment_proportions = treatment_proportion,
      generation_method     = generation_method
    )
    
    # Store results: one row per multiplier
    results_list[[mult]] <- data.frame(
      sample_size     = current_sample_size,
      stringsAsFactors = FALSE
    )
    results_list[[mult]]$simulation_data <- list(trial_results)
  }
  
  results_df <- do.call(rbind, results_list)
  rownames(results_df) <- NULL
  return(results_df)
}
###################running and saving the simulation########################
results_df <- run_main_simulation(total_sample_size   = 840,
                                  max_sample_multiply = 10,
                                  num_sites           = 2,
                                  iterations          = 1000,
                                  generation_method   = "distribution",  # or "bootstrap"
                                  treatment_proportion = 0.5)

colnames(results_df)
colnames(results_df$simulation_data[[1]]) #output of site number
results_df$simulation_data[[1]]$z_test[1] #output of the first set number's first repeat of simulation
saveRDS(results_df, file = "results_df_pw__25_75_unif_-0.7_0.9_x10.rds")
results_df <- readRDS("results_df_pw_naive_50_Normal_0.1_0.05.rds")



##############plotting#####################
library(ggplot2)
library(dplyr)
library(tidyr)

plot_test_power_by_sample_size <- function(results_df) {
  # Define test columns to summarize
  test_names <- c(
    "naive_t_test",
    "fix_effect_model_test",
    "mixed_effect_model_test",
    "mantel_haenszel_test"
    # Add more tests if needed
  )
  
  # Build summary dataframe: one row per (sample_size, test)
  summary_df <- lapply(seq_len(nrow(results_df)), function(i) {
    sample_size <- results_df$sample_size[i]
    sim_data <- results_df$simulation_data[[i]]
    
    means <- sapply(test_names, function(test) {
      mean(sim_data[[test]], na.rm = TRUE)
    })
    
    data.frame(sample_size = sample_size, t(means))
  }) %>% bind_rows()
  
  # Convert to long format for ggplot
  plot_data <- summary_df %>%
    pivot_longer(
      cols = all_of(test_names),
      names_to = "test_type",
      values_to = "mean_result"
    )
  
  # Plot
  ggplot(plot_data, aes(x = sample_size, y = mean_result, color = test_type)) +
    geom_line(size = 1.2) +
    geom_point(size = 2) +
    labs(
      title = "Mean Power of Each Test by Total Sample Size (2 sites)",
      x = "Total Sample Size",
      y = "Proportion of Simulations Passed (Power)",
      color = "Test Type"
    ) +
    ylim(0.25, 0.75) +
    theme_minimal() +
    theme(text = element_text(size = 12))
}

# Usage:
plot_test_power_by_sample_size(results_df)
