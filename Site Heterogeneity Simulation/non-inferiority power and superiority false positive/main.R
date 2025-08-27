#load in files
source("simulate trial.R")

#####################run simulations
run_main_simulation <- function(
    total_sample_size = 964,
    max_sites = 5,
    iterations = 10,
    generation_method = "distribution",  # or "bootstrap"
    treatment_proportion = 0.5,
    individual_sd = 0.01,
    treatment_sd = 0.03,
    recruitment_sd = 0.04
) {
  results_list <- list()
  
  for (num_sites in 2:max_sites) {
    cat("Running simulation with", num_sites, "sites...\n")
    
    # Run full set of simulations for this site count
    trial_results <- simulate_multiple_trials(
      num_sites = num_sites,
      specification = "random_selection",
      iterations = iterations,
      total_sample_size = total_sample_size,
      treatment_proportions = treatment_proportion,
      generation_method = generation_method,
      individual_sd = individual_sd,
      treatment_sd = treatment_sd,
      recruitment_sd = recruitment_sd
    )
    
    # Store as a row with two columns: num_sites, simulation_data
    results_list[[num_sites]] <- data.frame(
      num_sites = num_sites,
      stringsAsFactors = FALSE
    )
    results_list[[num_sites]]$simulation_data <- list(trial_results)
  }
  
  # Combine all rows into a single dataframe
  results_df <- do.call(rbind, results_list)
  rownames(results_df) <- NULL
  
  return(results_df)
}


###################running and saving the simulation, none inferiority fp ########################
# results_df <- run_main_simulation(total_sample_size = 824,
#                                   max_sites = 25,
#                                   iterations = 2500,
#                                   generation_method = "distribution",  # or "bootstrap"
#                                   treatment_proportion = 0.5,
#                                   individual_sd = 0.1,
#                                   treatment_sd = 0.0,
#                                   recruitment_sd = 0.0
# )
# saveRDS(results_df, file = "non-inferiority_individualsd=0.1_treatmentsd = 0.0_recruitmentsd=0.0.rds")

# results_df <- run_main_simulation(total_sample_size = 824,
#                                   max_sites = 25,
#                                   iterations = 2500,
#                                   generation_method = "distribution",  # or "bootstrap"
#                                   treatment_proportion = 0.5,
#                                   individual_sd = 0.05,
#                                   treatment_sd = 0.1,
#                                   recruitment_sd = 0.0
# )
# saveRDS(results_df, file = "non-inferiority_individualsd=0.05_treatmentsd = 0.1_recruitmentsd=0.0.rds.rds")

# results_df <- run_main_simulation(total_sample_size = 824,
#                                   max_sites = 25,
#                                   iterations = 2500,
#                                   generation_method = "distribution",  # or "bootstrap"
#                                   treatment_proportion = 0.5,
#                                   individual_sd = 0.05,
#                                   treatment_sd = 0.0,
#                                   recruitment_sd = 0.1
# )
# saveRDS(results_df, file = "non-inferiority_individualsd=0.05_treatmentsd = 0.0_recruitmentsd=0.1.rds.rds")

# results_df <- run_main_simulation(total_sample_size = 824,
#                                   max_sites = 25,
#                                   iterations = 2500,
#                                   generation_method = "distribution",  # or "bootstrap"
#                                   treatment_proportion = 0.5,
#                                   individual_sd = 0.05,
#                                   treatment_sd = 0.1,
#                                   recruitment_sd = 0.1
# )
# saveRDS(results_df, file = "non-inferiority_individualsd=0.05_treatmentsd = 0.1_recruitmentsd=0.1.rds")
###################running and saving the simulation, superiority fp########################
source("simulate trial.R")
results_df <- run_main_simulation(total_sample_size = 840,
                                  max_sites = 35,
                                  iterations = 2500,
                                  generation_method = "distribution",  # or "bootstrap"
                                  treatment_proportion = 0.5,
                                  individual_sd = 0.1,
                                  treatment_sd = 0.0,
                                  recruitment_sd = 0.0
)
saveRDS(results_df, file = "superiority_individualsd=0.1_treatmentsd = 0.0_recruitmentsd=0.0.rds")

results_df <- run_main_simulation(total_sample_size = 840,
                                  max_sites = 35,
                                  iterations = 2500,
                                  generation_method = "distribution",  # or "bootstrap"
                                  treatment_proportion = 0.5,
                                  individual_sd = 0.05,
                                  treatment_sd = 0.1,
                                  recruitment_sd = 0.0
)
saveRDS(results_df, file = "superiority_individualsd=0.05_treatmentsd = 0.1_recruitmentsd=0.0.rds.rds")

results_df <- run_main_simulation(total_sample_size = 840,
                                  max_sites = 35,
                                  iterations = 2500,
                                  generation_method = "distribution",  # or "bootstrap"
                                  treatment_proportion = 0.5,
                                  individual_sd = 0.05,
                                  treatment_sd = 0.0,
                                  recruitment_sd = 0.1
)
saveRDS(results_df, file = "superiority_individualsd=0.05_treatmentsd = 0.0_recruitmentsd=0.1.rds.rds")

results_df <- run_main_simulation(total_sample_size = 840,
                                  max_sites = 35,
                                  iterations = 2500,
                                  generation_method = "distribution",  # or "bootstrap"
                                  treatment_proportion = 0.5,
                                  individual_sd = 0.05,
                                  treatment_sd = 0.1,
                                  recruitment_sd = 0.1
)
saveRDS(results_df, file = "superiority_individualsd=0.05_treatmentsd = 0.1_recruitmentsd=0.1.rds")

###################running and saving the simulation, none inferiority power ########################
source("simulate trial non-inferiority.R")
results_df <- run_main_simulation(total_sample_size = 824,
                                  max_sites = 35,
                                  iterations = 2500,
                                  generation_method = "distribution",  # or "bootstrap"
                                  treatment_proportion = 0.5,
                                  individual_sd = 0.1,
                                  treatment_sd = 0.0,
                                  recruitment_sd = 0.0
)
saveRDS(results_df, file = "non-inferiority_individualsd=0.1_treatmentsd = 0.0_recruitmentsd=0.0.rds")

results_df <- run_main_simulation(total_sample_size = 824,
                                  max_sites = 35,
                                  iterations = 2500,
                                  generation_method = "distribution",  # or "bootstrap"
                                  treatment_proportion = 0.5,
                                  individual_sd = 0.05,
                                  treatment_sd = 0.1,
                                  recruitment_sd = 0.0
)
saveRDS(results_df, file = "non-inferiority_individualsd=0.05_treatmentsd = 0.1_recruitmentsd=0.0.rds.rds")

results_df <- run_main_simulation(total_sample_size = 824,
                                  max_sites = 35,
                                  iterations = 2500,
                                  generation_method = "distribution",  # or "bootstrap"
                                  treatment_proportion = 0.5,
                                  individual_sd = 0.05,
                                  treatment_sd = 0.0,
                                  recruitment_sd = 0.1
)
saveRDS(results_df, file = "non-inferiority_individualsd=0.05_treatmentsd = 0.0_recruitmentsd=0.1.rds.rds")

results_df <- run_main_simulation(total_sample_size = 824,
                                  max_sites = 35,
                                  iterations = 2500,
                                  generation_method = "distribution",  # or "bootstrap"
                                  treatment_proportion = 0.5,
                                  individual_sd = 0.05,
                                  treatment_sd = 0.1,
                                  recruitment_sd = 0.1
)
saveRDS(results_df, file = "non-inferiority_individualsd=0.05_treatmentsd = 0.1_recruitmentsd=0.1.rds")
