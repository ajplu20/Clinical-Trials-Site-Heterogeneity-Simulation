source("data generation naive.R")   #Generating a none informative dataset, for speed only.
#print("gen naive done")

#source("select sites.R")
source("select sites naive.R") #generate fake sites instead of using the acornhai dataset.
#print("help done")

source("generate_site_counts_per_site.R")

source("predict mortality naive 40%.R") #used if we want custom defined mortality instead of predicted from dataset


source("generate trial effect non-inferiority.R") #
#print("naive done")

source("perform_naive_t_test.R") 
#print("t done")

source("perform_fix_effect_model_test.R")
#print("fixed done")

source("perform_mixed_effect_model_test.R")
#print("mixed done")

source("perform_mantel_haenszel_test.R")
#print("mh done")

#####################################simulation function for 1 trial ###########################
library(stats)
#model <- readRDS("glmer_model2.rds")
#print("in")
simulate_trial <- function(
    sites,
    site_counts,
    treatment_proportions,
    generation_method = "distribution",
    trial_effect,                # numeric vector, one entry per site
    individual_sd = 0.01
) {
  ## ---------- 1. Input checks ----------
  if (length(sites) != length(site_counts) ||
      length(sites) != length(treatment_proportions) ||
      length(sites) != length(trial_effect)) {
    stop("sites, site_counts, treatment_proportions, and trial_effect must all have the same length.")
  }
  
  ## Give every helper vector the site names so we can look them up later
  names(site_counts)          <- sites
  names(treatment_proportions) <- sites
  names(trial_effect)          <- sites 
  
  ## ---------- 2. Generate synthetic patients ----------
  patient_data <- switch(tolower(generation_method),
                         bootstrap    = generate_random_dataset_bootstrap(sites, site_counts),
                         distribution = generate_random_dataset_distribution(sites, site_counts),
                         stop("generation_method must be 'bootstrap' or 'distribution'"))
  
  ## ---------- 3. Assign experiment vs control ----------
  patient_data <- patient_data |>
    dplyr::group_by(site)   |>
    dplyr::mutate(row_id = dplyr::row_number()) |>
    dplyr::ungroup()
  
  experiment_vector <- rep(NA_integer_, nrow(patient_data))
  
  for (site in sites) {
    site_rows <- which(patient_data$site == site)
    n_site    <- length(site_rows)
    prop      <- treatment_proportions[site]
    
    n_experiment <- ceiling(prop       * n_site)
    n_control    <- ceiling((1 - prop) * n_site)
    total_needed <- n_experiment + n_control
    
    ## Top‑up with resampled rows if rounding made us short
    if (total_needed > n_site) {
      extra_data <- dplyr::slice_sample(
        patient_data[site_rows, ], n = total_needed - n_site, replace = TRUE)
      extra_data$row_id <- max(patient_data$row_id, na.rm = TRUE) + seq_len(nrow(extra_data))
      patient_data <- dplyr::bind_rows(patient_data, extra_data)
      site_rows    <- which(patient_data$site == site)   # refresh indices
    }
    
    chosen <- sample(site_rows, n_experiment)
    experiment_vector[chosen]               <- 1
    experiment_vector[setdiff(site_rows, chosen)] <- 0
  }
  patient_data$experiment <- experiment_vector
  
  ## ---------- 4. Baseline mortality probabilities ----------
  predictors <- c("site","cal_agey","sex","vap","hpd_admreason","comorbidities_CCI")
  patient_data <- dplyr::filter(patient_data, stats::complete.cases(dplyr::across(dplyr::all_of(predictors))))
  #print(head(patient_data))
  #print(patient_data$site)
  patient_data$mortality_probability <- predict_mortality_from_data(patient_data, individual_sd = individual_sd)
  #print(head(patient_data$mortality_probability))
  
  ## ---------- 5. Apply site‑specific trial effect ----------
  patient_data$mortality_probability <- mapply(
    function(prob, exp, site) {
      if (exp == 1L) {
        # reduce probability, floor at 0
        min(max(prob - trial_effect[site], 0),1)
      } else prob
    },
    prob = patient_data$mortality_probability,
    exp  = patient_data$experiment,
    site = patient_data$site
  )
  
  ## ---------- 6. Simulate death outcome ----------
  
  #print("mortality is")
  #print(patient_data$mortality_probability)
  hold <- stats::rbinom(nrow(patient_data), size = 1, prob = patient_data$mortality_probability)
  #print("minmax is")
  #print(min(1,max(0,patient_data$mortality_probability)))
  #print("mortality_probabilit is")
  #print(patient_data$mortality_probability)
  #if(any(is.na(hold))){
  #  print(min(1,max(0,patient_data$mortality_probability)))
  #  print(patient_data$mortality_probability)
  #  print(max(0,patient_data$mortality_probability))
  #  print("hold_NA")
  #}
  patient_data$death <- hold
  ## ---------- 7. Clean ups ----------
  patient_data <- dplyr::select(patient_data, -row_id)
  return(patient_data)
}



#########example usage
## 1.Inputs -------------------------------------------------
sites                 <- c("CN 001", "JP 001", "PH 001")
site_counts           <- c(1000, 1000, 1000)
treatment_proportions <- c(0.5, 0.5, 0.5)   # 50%, 60%, 40% experiment
trial_effect          <- c(-0.1, -0.1, -0.1)  # site‑specific risk reduction
#trial_effect          <- c(-0, -0, -0)  # site‑specific risk reduction


## 3.Run one simulated trial -------------------------------
set.seed(123)
sim_data <- simulate_trial(
  sites                 = sites,
  site_counts           = site_counts,
  treatment_proportions = treatment_proportions,
  generation_method     = "distribution",
  trial_effect          = trial_effect,
  individual_sd = 0.05
)
df <- sim_data
#sd(sim_data[sim_data$site == "CN 001" & sim_data$experiment == 0,]$mortality_probability)

## 4.Inspect ------------------------------------------------
head(sim_data)
table(sim_data$site, sim_data$experiment)
aggregate(death ~ site + experiment, data = sim_data, mean)
mean(sim_data$mortality_probability[sim_data$site =="JP 001" & sim_data$experiment ==0])


###########################run multiple trials, with dataframe output as well###########
#site_count_list <- c()

simulate_multiple_trials <- function(
    num_sites,
    specification = "random_selection",
    iterations = 1000,
    total_sample_size = 858,
    treatment_proportions = 0.5,
    generation_method = "distribution",
    individual_sd = 0.01,
    treatment_sd = 0.01,
    recruitment_sd = 0.01
) {
  
  # Preallocate result dataframe with an extra column for trial data
  results <- data.frame(
    difference = integer(iterations),
    naive_t_test = integer(iterations) ,
    fix_effect_model_test = integer(iterations),
    mixed_effect_model_test = integer(iterations),
    mantel_haenszel_test = integer(iterations),
    reject = integer(iterations)
  )
  
  # Add an empty list-column to store full trial data
  results$trial_data <- vector("list", iterations)
  
  # Handle flexible site_counts
  
  
  
  # Handle flexible treatment proportions
  if (length(treatment_proportions) == 1) {
    treatment_props <- rep(treatment_proportions, num_sites)
  } else if (length(treatment_proportions) == num_sites) {
    treatment_props <- treatment_proportions
  } else {
    stop("treatment_proportions must be either a single number or a vector of length num_sites.")
  }
  
  #trial_effect_list <- c()
  site_count_list <- c()
  for (i in 1:iterations) {
    selected_sites <- select_sites(num_sites, specification)
    
    site_counts <- generate_site_counts_per_site(selected_sites, total_sample_size, recruitment_sd)
    
    effects <- generate_trial_effect(num_sites, treatment_sd)
    #trial_effect_list <- append(trial_effect_list, effects)
    site_count_list <- append(site_count_list, site_counts)
    #print(site_counts)
    trial_data <- simulate_trial(
      sites = selected_sites,
      site_counts = site_counts,
      treatment_proportions = treatment_props,
      generation_method = generation_method,
      trial_effect = effects
    )
    
    # find the proportion of death in control and experiement group
    p_death_control <- mean(trial_data$death[trial_data$experiment == 0])
    #print(p_death_control)
    p_death_treatment <- mean(trial_data$death[trial_data$experiment == 1])
    
    results$difference[i] <-  p_death_treatment - p_death_control
    results$naive_t_test[i] <- perform_naive_t_test(trial_data, alpha = 0.05, trial_type = "non-inferiority")
    results$fix_effect_model_test[i]  <- perform_fix_effect_model_test(trial_data, alpha = 0.05, trial_type = "non-inferiority",basket = FALSE)
    results$mixed_effect_model_test[i]<- perform_mixed_effect_model_test(trial_data, alpha = 0.05, trial_type = "non-inferiority",basket = FALSE)
    results$mantel_haenszel_test[i]   <- perform_mantel_haenszel_test(trial_data, alpha = 0.05, trial_type = "non-inferiority-naive")
    results$trial_data[[i]] <- trial_data
  }
  difference_sd <- sd(results$difference)
  for (i in 1:iterations) {
    results$reject[i] <- (results$difference[i] + difference_sd * qnorm(0.975) > 0.1)
  }
  
  
  #hist(site_count_list)
  #print(sd(site_count_list))
  return(results)
}


##########testing simulate multiple trials###########################
# Simulate 500 trials using 6 sites, each with 120 patients, with random site selection
# result_df <- simulate_multiple_trials(
#   num_sites = 50,
#   specification = "random_selection",
#   iterations = 250,
#   total_sample_size = 824,
#   generation_method = "distribution",
#   individual_sd = 0.05,
#   treatment_sd = 0.1,
#   recruitment_sd = 0.2
# )



# View estimated power (proportion of trials rejecting H0)
#mean(result_df$naive_t_test)
#mean(result_df$fix_effect_model_test)
#mean(result_df$mixed_effect_model_test)
#mean(result_df$mantel_haenszel_test)

# mean(result_df$z_test)


#num_sites = 2
#specification = "random_selection"
#iterations = 500
#total_sample_size = 824
#generation_method = "distribution"
#individual_sd = 0.05
#treatment_sd = 0.0
#recruitment_sd = 0.00
#treatment_proportions = 0.5

