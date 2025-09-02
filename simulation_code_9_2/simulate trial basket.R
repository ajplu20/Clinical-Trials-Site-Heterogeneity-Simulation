# ---- Common sources (data gen, site selection, counts, tests) ----
source("data generation naive.R")      # generate_random_dataset_*()
source("select sites naive.R")         # select_sites()
source("generate_site_counts_per_site.R")

source("perform_naive_t_test.R")
source("perform_fix_effect_model_test.R")
source("perform_mixed_effect_model_test.R")
source("perform_mantel_haenszel_test.R")

# unified effect + mortality helpers
source("generate trial effect.R")      # generate_trial_effect(num_sites, mu, treatment_sd)
source("predict mortality naive.R")    # predict_mortality_from_data(df, min_mortality, max_mortality, individual_sd)

# basket helpers
source("generate_syndrome_effect_df.R")
source("generate_syndrome_prop_df.R")
source("split_into_groups.R")

library(stats)

# ---- Single-trial simulator (BASKET-capable, shared) ----
simulate_trial <- function(
    sites,
    site_counts,
    treatment_proportions,
    generation_method = "distribution",
    trial_effect,                     # numeric vector, one per site (positive lowers risk)
    individual_sd = 0.01,
    min_mortality = 0.20,
    max_mortality = 0.60,
    # basket pieces
    syndrome_effect_df,
    syndrome_proportion_df,
    stratify_by_syndrome = FALSE
) {
  # basic checks
  if (length(sites) != length(site_counts) ||
      length(sites) != length(treatment_proportions) ||
      length(sites) != length(trial_effect)) {
    stop("sites, site_counts, treatment_proportions, and trial_effect must all have the same length.")
  }
  
  # name helper vectors
  names(site_counts)            <- sites
  names(treatment_proportions)  <- sites
  names(trial_effect)           <- sites
  
  # 2) Generate patients
  patient_data <- switch(tolower(generation_method),
                         bootstrap    = generate_random_dataset_bootstrap(sites, site_counts),
                         distribution = generate_random_dataset_distribution(sites, site_counts),
                         stop("generation_method must be 'bootstrap' or 'distribution'")
  )
  
  # 3) Assign arm (+ optional syndrome-stratified allocation)
  patient_data <- patient_data |>
    dplyr::group_by(site) |>
    dplyr::mutate(row_id = dplyr::row_number()) |>
    dplyr::ungroup()
  
  experiment_vector <- rep(NA_integer_, nrow(patient_data))
  syndrome_vector   <- rep(NA,            nrow(patient_data))  # stores syndrome labels
  
  for (site in sites) {
    site_rows <- which(patient_data$site == site)
    n_site    <- length(site_rows)
    prop      <- treatment_proportions[site]
    
    if (stratify_by_syndrome) {
      # arm sizes by syndrome (round up)
      ctrl_counts <- ceiling(syndrome_proportion_df[site,] * (1 - prop) * n_site)
      exp_counts  <- ceiling(syndrome_proportion_df[site,] * (    prop) * n_site)
      n_control   <- sum(ctrl_counts)
      n_experiment<- sum(exp_counts)
    } else {
      # total per-syndrome -> then split into arms
      total_counts <- ceiling(syndrome_proportion_df[site,] * n_site)
      n_experiment <- ceiling(prop       * sum(total_counts))
      n_control    <- ceiling((1 - prop) * sum(total_counts))
    }
    total_needed <- n_control + n_experiment
    
    # top up if rounding exceeded available rows
    if (total_needed > n_site) {
      extra <- dplyr::slice_sample(patient_data[site_rows,], n = total_needed - n_site, replace = TRUE)
      extra$row_id <- max(patient_data$row_id, na.rm = TRUE) + seq_len(nrow(extra))
      patient_data <- dplyr::bind_rows(patient_data, extra)
      site_rows    <- which(patient_data$site == site) # refresh indices
    }
    
    # assign treatment
    chosen <- sample(site_rows, n_experiment)
    experiment_vector[chosen] <- 1L
    experiment_vector[setdiff(site_rows, chosen)] <- 0L
    
    # assign syndrome labels
    if (stratify_by_syndrome) {
      # split EXP rows into syndrome buckets
      exp_groups  <- split_into_groups(chosen, as.numeric(exp_counts))
      # split CTRL rows into syndrome buckets
      ctrl_groups <- split_into_groups(setdiff(site_rows, chosen), as.numeric(ctrl_counts))
      k <- 1
      for (synd in colnames(syndrome_proportion_df)) {
        if (length(exp_groups)  >= k) syndrome_vector[exp_groups[[k]]]  <- synd
        if (length(ctrl_groups) >= k) syndrome_vector[ctrl_groups[[k]]] <- synd
        k <- k + 1
      }
    } else {
      # no stratified allocation: assign syndrome for all rows first
      total_counts <- ceiling(syndrome_proportion_df[site,] * length(site_rows))
      groups <- split_into_groups(site_rows, as.numeric(total_counts))
      k <- 1
      for (synd in colnames(syndrome_proportion_df)) {
        if (length(groups) >= k) syndrome_vector[groups[[k]]] <- synd
        k <- k + 1
      }
    }
  }
  
  patient_data$experiment <- experiment_vector
  patient_data$syndrome   <- syndrome_vector
  
  # 4) Baseline mortality (site gradient + noise), clipped
  predictors <- c("site","cal_agey","sex","vap","hpd_admreason","comorbidities_CCI")
  patient_data <- dplyr::filter(patient_data, stats::complete.cases(dplyr::across(dplyr::all_of(predictors))))
  
  patient_data$mortality_probability <- predict_mortality_from_data(
    patient_data,
    min_mortality = min_mortality,
    max_mortality = max_mortality,
    individual_sd = individual_sd
  )
  
  # 5) Apply site + syndrome effects on treatment arm
  patient_data$mortality_probability <- mapply(
    function(prob, exp, site, synd) {
      if (exp == 1L) {
        # subtract (positive values = benefit), floor/ceil to [0,1]
        min(max(prob - trial_effect[site] - syndrome_effect_df[site, synd], 0), 1)
      } else prob
    },
    prob = patient_data$mortality_probability,
    exp  = patient_data$experiment,
    site = patient_data$site,
    synd = patient_data$syndrome
  )
  
  # 6) Simulate outcome
  patient_data$death <- stats::rbinom(nrow(patient_data), size = 1, prob = patient_data$mortality_probability)
  
  # 7) Clean
  dplyr::select(patient_data, -row_id)
}

# ---- Multi-trial simulator (unified: superiority & NI; basket-aware) ----
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
    treatment_sd = 0.00,
    # Recruitment variability
    recruitment_sd = 0.01,
    # Design intent
    design_type = c("superiority","non-inferiority"),
    alpha = 0.05,
    ni_margin = 0.10,                 # risk-difference NI margin (treat - control)
    # Basket knobs
    syndrome_num = 2,
    syndrome_prop_var = 0,            # 0 = even allocation across syndromes, defines how imbalanced syndrom is in a site
    syndrome_effect = c(0.1, -0.1),   # impact syndrom on treatment effect,  0.1 is further reduce mortality by 0.1.
    syndrome_effect_variation = c(0, 0),  #how much syndrom impact on treatment varies across sites. 0 is no variation.
    stratify_by_syndrome = FALSE,     # do we stratify by syndrom in addition to site when allocating to treatment vs experiment?
    # Output toggles
    save_trial_data = TRUE,
    save_syndrome_frames = TRUE
) {
  design_type <- match.arg(design_type)
  
  # treatment props vector
  if (length(treatment_proportions) == 1) {
    treatment_props <- rep(treatment_proportions, num_sites)
  } else if (length(treatment_proportions) == num_sites) {
    treatment_props <- treatment_proportions
  } else {
    stop("treatment_proportions must be either a single number or a vector of length num_sites.")
  }
  
  # results frame
  base_cols <- list(
    naive_t_test = integer(iterations),
    fix_effect_model_test = integer(iterations),
    mixed_effect_model_test = integer(iterations),
    mantel_haenszel_test = integer(iterations)
  )
  if (design_type == "non-inferiority") {
    base_cols$difference <- numeric(iterations)  # p_t - p_c
    base_cols$reject     <- integer(iterations)  # NI via Wald CI
  }
  results <- as.data.frame(base_cols, stringsAsFactors = FALSE)
  if (save_trial_data)        results$trial_data              <- vector("list", iterations)
  if (save_syndrome_frames) { results$syndrome_effect_df      <- vector("list", iterations)
  results$syndrome_proportion_df  <- vector("list", iterations) }
  
  for (i in seq_len(iterations)) {
    results$state_before[[i]] <- .Random.seed
    # pick sites & counts
    selected_sites <- select_sites(num_sites, specification)
    
    # make per-site syndrome effects & mix
    # (use provided helpers; length of syndrome_effect controls number of syndromes if not using syndrome_num)
    synd_eff_df  <- generate_syndrome_effect_df(selected_sites, syndrome_effect, syndrome_effect_variation)
    synd_prop_df <- generate_syndrome_prop_df(selected_sites, syndrome_effect, syndrome_prop_var)
    if (save_syndrome_frames) {
      results$syndrome_effect_df[[i]]     <- synd_eff_df
      results$syndrome_proportion_df[[i]] <- synd_prop_df
    }
    
    site_counts <- generate_site_counts_per_site(selected_sites, total_sample_size, recruitment_sd)
    effects     <- generate_trial_effect(num_sites, mu = mu, treatment_sd = treatment_sd)
    
    # run one trial
    trial_data <- simulate_trial(
      sites = selected_sites,
      site_counts = site_counts,
      treatment_proportions = treatment_props,
      generation_method = generation_method,
      trial_effect = effects,
      individual_sd = individual_sd,
      min_mortality = min_mortality,
      max_mortality = max_mortality,
      syndrome_effect_df = synd_eff_df,
      syndrome_proportion_df = synd_prop_df,
      stratify_by_syndrome = stratify_by_syndrome
    )
    results$state_after[[i]] <- .Random.seed
    # NI summaries (risk-diff Wald CI)
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
      results$reject[i]     <- as.integer(ci_hi > ni_margin)
    }
    
    # pass design to tests; basket=TRUE for FE/Mixed
    tt   <- if (design_type == "non-inferiority") "non-inferiority" else "superiority"
    mh_t <- if (design_type == "non-inferiority") "non-inferiority-naive" else "superiority"
    
    results$naive_t_test[i]            <- perform_naive_t_test(trial_data, alpha = alpha, trial_type = tt)
    results$fix_effect_model_test[i]   <- perform_fix_effect_model_test(trial_data, alpha = alpha, trial_type = tt, basket = TRUE)
    results$mixed_effect_model_test[i] <- perform_mixed_effect_model_test(trial_data, alpha = alpha, trial_type = tt, basket = TRUE)
    results$mantel_haenszel_test[i]    <- perform_mantel_haenszel_test(trial_data, alpha = alpha, trial_type = mh_t)
    
    if (save_trial_data) results$trial_data[[i]] <- trial_data
  }
  
  results
}

# # -------------------------- Example usage --------------------------
# # # Superiority (e.g., 25–75%)
# res_sup_basket <- simulate_multiple_trials(
#   num_sites = 3,
#   specification = "random_selection",
#   iterations = 1,
#   total_sample_size = 900,
#   generation_method = "distribution",
#   min_mortality = 0.25, max_mortality = 0.75,
#   individual_sd = 0.05,
#   mu = 0.1, treatment_sd = 0.00,
#   recruitment_sd = 0.00,
#   design_type = "superiority",
#   alpha = 0.05,
#   syndrome_num = 2,
#   syndrome_prop_var = 0,
#   syndrome_effect = c(0.10, -0.10),
#   syndrome_effect_variation = c(0, 0),
#   stratify_by_syndrome = FALSE,
#   save_trial_data = TRUE,
#   save_syndrome_frames = TRUE
# )
# 
# trial1 <- res_sup_basket$trial_data[[1]]
# 
# site_arm_summary <- aggregate(
#   cbind(mortality_probability, death) ~ site + experiment,
#   data = trial1,
#   FUN = mean
# )
# 
# names(site_arm_summary)[names(site_arm_summary) == "mortality_probability"] <- "mean_mortality_prob"
# names(site_arm_summary)[names(site_arm_summary) == "death"]                 <- "mean_death_rate"
# 
# cat("=== SITE × ARM (overall) ===\n")
# print(site_arm_summary[order(site_arm_summary$site, site_arm_summary$experiment), ])
# 
# syndrome_arm_summary <- aggregate(
#   cbind(mortality_probability, death) ~ syndrome + experiment,
#   data = trial1,
#   FUN = mean
# )
# 
# names(syndrome_arm_summary)[names(syndrome_arm_summary) == "mortality_probability"] <- "mean_mortality_prob"
# names(syndrome_arm_summary)[names(syndrome_arm_summary) == "death"]                 <- "mean_death_rate"
# 
# cat("\n=== SYNDROME × ARM (overall) ===\n")
# print(syndrome_arm_summary[order(syndrome_arm_summary$syndrome, syndrome_arm_summary$experiment), ])
# 
# # # Non-inferiority (≈40% baseline profile)
# res_ni_basket <- simulate_multiple_trials(
#   num_sites = 3,
#   specification = "random_selection",
#   iterations = 1,
#   total_sample_size = 900,
#   generation_method = "distribution",
#   min_mortality = 0.20, max_mortality = 0.60,
#   individual_sd = 0.05,
#   mu = 0.00, treatment_sd = 0.00,
#   recruitment_sd = 0.00,
#   design_type = "non-inferiority",
#   alpha = 0.05, ni_margin = 0.10,
#   syndrome_num = 2,
#   syndrome_prop_var = 0,
#   syndrome_effect = c(0.10, -0.10),
#   syndrome_effect_variation = c(0, 0),
#   stratify_by_syndrome = FALSE,
#   save_trial_data = TRUE,
#   save_syndrome_frames = TRUE
# )
# 
# trial1_ni <- res_ni_basket$trial_data[[1]]
# 
# site_arm_summary_ni <- aggregate(
#   cbind(mortality_probability, death) ~ site + experiment,
#   data = trial1_ni,
#   FUN = mean
# )
# names(site_arm_summary_ni)[names(site_arm_summary_ni) == "mortality_probability"] <- "mean_mortality_prob"
# names(site_arm_summary_ni)[names(site_arm_summary_ni) == "death"]                 <- "mean_death_rate"
# 
# cat("=== (NI) SITE × ARM (overall) ===\n")
# print(site_arm_summary_ni[order(site_arm_summary_ni$site, site_arm_summary_ni$experiment), ])
# 
# syndrome_arm_summary_ni <- aggregate(
#   cbind(mortality_probability, death) ~ syndrome + experiment,
#   data = trial1_ni,
#   FUN = mean
# )
# names(syndrome_arm_summary_ni)[names(syndrome_arm_summary_ni) == "mortality_probability"] <- "mean_mortality_prob"
# names(syndrome_arm_summary_ni)[names(syndrome_arm_summary_ni) == "death"]                 <- "mean_death_rate"
# 
# cat("\n=== (NI) SYNDROME × ARM (overall) ===\n")
# print(syndrome_arm_summary_ni[order(syndrome_arm_summary_ni$syndrome, syndrome_arm_summary_ni$experiment), ])