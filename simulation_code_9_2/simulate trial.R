# ---- Common sources (data gen, site selection, counts, tests) ----
source("data generation naive.R")      # generate_random_dataset_*()
source("select sites naive.R")         # select_sites()
source("generate_site_counts_per_site.R")

source("perform_naive_t_test.R")
source("perform_fix_effect_model_test.R")
source("perform_mixed_effect_model_test.R")
source("perform_mantel_haenszel_test.R")
source("generate trial effect.R")
source("predict mortality naive.R")
library(stats)

# ---- Single-trial simulator (shared) ----
simulate_trial <- function(
    sites,
    site_counts,
    treatment_proportions,
    generation_method = "distribution",
    trial_effect,                     # numeric vector, one per site (positive lowers risk)
    individual_sd = 0.01,
    min_mortality = 0.20,
    max_mortality = 0.60
) {
  if (length(sites) != length(site_counts) ||
      length(sites) != length(treatment_proportions) ||
      length(sites) != length(trial_effect)) {
    stop("sites, site_counts, treatment_proportions, and trial_effect must all have the same length.")
  }
  
  names(site_counts)           <- sites
  names(treatment_proportions) <- sites
  names(trial_effect)          <- sites
  
  patient_data <- switch(tolower(generation_method),
                         bootstrap    = generate_random_dataset_bootstrap(sites, site_counts),
                         distribution = generate_random_dataset_distribution(sites, site_counts),
                         stop("generation_method must be 'bootstrap' or 'distribution'"))
  
  patient_data <- patient_data |>
    dplyr::group_by(site) |>
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
    
    if (total_needed > n_site) {
      extra_data <- dplyr::slice_sample(
        patient_data[site_rows, ], n = total_needed - n_site, replace = TRUE)
      extra_data$row_id <- max(patient_data$row_id, na.rm = TRUE) + seq_len(nrow(extra_data))
      patient_data <- dplyr::bind_rows(patient_data, extra_data)
      site_rows    <- which(patient_data$site == site)
    }
    
    chosen <- sample(site_rows, n_experiment)
    experiment_vector[chosen] <- 1L
    experiment_vector[setdiff(site_rows, chosen)] <- 0L
  }
  patient_data$experiment <- experiment_vector
  
  # Baseline mortality (site gradient + noise), clipped
  predictors <- c("site","cal_agey","sex","vap","hpd_admreason","comorbidities_CCI")
  patient_data <- dplyr::filter(patient_data, stats::complete.cases(dplyr::across(dplyr::all_of(predictors))))
  
  patient_data$mortality_probability <- predict_mortality_from_data(
    patient_data,
    min_mortality = min_mortality,
    max_mortality = max_mortality,
    individual_sd = individual_sd
  )
  
  # Apply site-specific effect to treatment arm: subtract (positive = benefit)
  patient_data$mortality_probability <- mapply(
    function(prob, exp, site) {
      if (exp == 1L) {
        min(max(prob - trial_effect[site], 0), 1)
      } else prob
    },
    prob = patient_data$mortality_probability,
    exp  = patient_data$experiment,
    site = patient_data$site
  )
  
  patient_data$death <- stats::rbinom(nrow(patient_data), size = 1,
                                      prob = patient_data$mortality_probability)
  dplyr::select(patient_data, -row_id)
}

# ---- Multi-trial simulator (unified for superiority & non-inferiority) ----
simulate_multiple_trials <- function(
    num_sites,
    specification = "random_selection",
    iterations = 1000,
    total_sample_size = 900,
    treatment_proportions = 0.5,
    generation_method = "distribution",
    # Mortality profile
    min_mortality = 0.20,
    max_mortality = 0.60,
    individual_sd = 0.01,
    # Effect generator
    mu = 0.10,
    treatment_sd = 0.00,  #how much treatment effect varies from site to site
    # Recruitment variability
    recruitment_sd = 0.01, #how much recruitment imbalance exists in the number of patients per site? 0 = balanced recruitment from all sites
    # Design intent
    design_type = c("superiority","non-inferiority"),
    alpha = 0.05,
    ni_margin = 0.10,            # risk-difference NI margin (treat - control)
    # Output toggles
    save_trial_data = TRUE
) {
  design_type <- match.arg(design_type)
  
  # Flexible treatment proportions
  if (length(treatment_proportions) == 1) {
    treatment_props <- rep(treatment_proportions, num_sites)
  } else if (length(treatment_proportions) == num_sites) {
    treatment_props <- treatment_proportions
  } else {
    stop("treatment_proportions must be either a single number or a vector of length num_sites.")
  }
  
  # Build result frame
  base_cols <- list(
    naive_t_test = integer(iterations),
    fix_effect_model_test = integer(iterations),
    mixed_effect_model_test = integer(iterations),
    mantel_haenszel_test = integer(iterations)
  )
  if (design_type == "non-inferiority") {
    base_cols$difference <- numeric(iterations)  # p_treat - p_control
    base_cols$reject     <- integer(iterations)  # NI success via Wald CI
  }
  
  results <- as.data.frame(base_cols, stringsAsFactors = FALSE)
  if (save_trial_data) results$trial_data <- vector("list", iterations)
  
  for (i in seq_len(iterations)) {
    results$state_before[[i]] <- .Random.seed
    selected_sites <- select_sites(num_sites, specification)
    site_counts <- generate_site_counts_per_site(selected_sites, total_sample_size, recruitment_sd)
    effects <- generate_trial_effect(num_sites, mu = mu, treatment_sd = treatment_sd)
    
    trial_data <- simulate_trial(
      sites = selected_sites,
      site_counts = site_counts,
      treatment_proportions = treatment_props,
      generation_method = generation_method,
      trial_effect = effects,
      individual_sd = individual_sd,
      min_mortality = min_mortality,
      max_mortality = max_mortality
    )
    results$state_after[[i]] <- .Random.seed
    # Per-trial summaries for NI (risk difference + Wald CI decision)
    if (design_type == "non-inferiority") {
      n_t <- sum(trial_data$experiment == 1L)
      n_c <- sum(trial_data$experiment == 0L)
      p_t <- mean(trial_data$death[trial_data$experiment == 1L])
      p_c <- mean(trial_data$death[trial_data$experiment == 0L])
      diff <- p_t - p_c
      se   <- sqrt(p_t*(1-p_t)/n_t + p_c*(1-p_c)/n_c)
      z    <- qnorm(1 - alpha/2)
      ci_hi <- diff + z*se
      results$difference[i] <- diff
      # Reject inferiority (i.e., declare NI) if upper CI bound > NI margin
      results$reject[i] <- as.integer(ci_hi > ni_margin)
    }
    
    # Common tests (pass intent + alpha)
    tt <- if (design_type == "non-inferiority") "non-inferiority" else "superiority"
    mh_tt <- if (design_type == "non-inferiority") "non-inferiority-naive" else "superiority"
    
    results$naive_t_test[i]           <- perform_naive_t_test(trial_data, alpha = alpha, trial_type = tt)
    results$fix_effect_model_test[i]  <- perform_fix_effect_model_test(trial_data, alpha = alpha, trial_type = tt, basket = FALSE)
    results$mixed_effect_model_test[i]<- perform_mixed_effect_model_test(trial_data, alpha = alpha, trial_type = tt, basket = FALSE)
    results$mantel_haenszel_test[i]   <- perform_mantel_haenszel_test(trial_data, alpha = alpha, trial_type = mh_tt)
    
    if (save_trial_data) results$trial_data[[i]] <- trial_data
  }
  
  results
}


# -------------------------- Example usage --------------------------

# # Superiority (e.g., 25–75% )
# set.seed(123)
# res_sup <- simulate_multiple_trials(
#   num_sites = 3,
#   specification = "random_selection",
#   iterations = 2,
#   total_sample_size = 1000,
#   generation_method = "distribution",
#   min_mortality = 0.25, max_mortality = 0.75,   
#   individual_sd = 0.05,
#   mu = 0.10, treatment_sd = 0.00,  #for false positive, mu = 0
#   recruitment_sd = 0.00,
#   design_type = "superiority",
#   alpha = 0.05
# )
# # Take the first simulated trial's data
# trial1 <- res_sup$trial_data[[1]]
#
# # Table of mean *mortality_probability* by site and arm
# mortality_prob_table <- aggregate(mortality_probability ~ site + experiment, 
#                                   data = trial1, 
#                                   mean)
#
# # Table of observed death rate by site and arm
# death_rate_table <- aggregate(death ~ site + experiment, 
#                               data = trial1, 
#                               mean)
#
# cat("Mean baseline mortality probability:\n")
# print(mortality_prob_table)
#
# cat("\nObserved death rates:\n")
# print(death_rate_table)
#
#
# # Non-inferiority (≈40% baseline profile)
# set.seed(123)
# res_ni <- simulate_multiple_trials(
#   num_sites = 3,
#   specification = "random_selection",
#   iterations = 2,
#   total_sample_size = 900,
#   generation_method = "distribution",
#   min_mortality = 0.2, max_mortality = 0.6,   
#   individual_sd = 0.05,
#   mu = 0.00, treatment_sd = 0.02,               # “true” diff ~ 0 #for false positive, mu = -0.1
#   recruitment_sd = 0.10,
#   design_type = "non-inferiority",
#   alpha = 0.05, ni_margin = 0.10
# )
# # Take the first simulated trial's data (non-inferiority)
# trial1_ni <- res_ni$trial_data[[1]]
#
# # Table of mean *mortality_probability* by site and arm
# mortality_prob_table_ni <- aggregate(mortality_probability ~ site + experiment, 
#                                      data = trial1_ni, 
#                                      mean)
#
# # Table of observed death rate by site and arm
# death_rate_table_ni <- aggregate(death ~ site + experiment, 
#                                  data = trial1_ni, 
#                                  mean)
#
# cat("Mean baseline mortality probability (non-inferiority):\n")
# print(mortality_prob_table_ni)
#
# cat("\nObserved death rates (non-inferiority):\n")
# print(death_rate_table_ni)
#
# mean(res_ni$naive_t_test)
# mean(res_ni$fix_effect_model_test)
# mean(res_ni$mixed_effect_model_test)
# mean(res_ni$mantel_haenszel_test)
