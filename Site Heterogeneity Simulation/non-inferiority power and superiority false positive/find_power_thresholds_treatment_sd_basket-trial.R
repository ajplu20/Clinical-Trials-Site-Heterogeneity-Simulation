# basket_parallel.R
# Parallel outer-loop across treatment_sd values for basket simulations
options(parallelly.availableCores.methods = c("env", "pbs", "system", "fallback"))
options(parallelly.cgroups.mounts.parse = FALSE)
Sys.setenv(R_PARALLELLY_AVAILABLECORES_METHODS = "env,pbs,system,fallback")

cat("=== R starting ===\n")

library(future)
library(future.apply)


# ----- One task = run the full inner site loop for a single treatment_sd -----
run_one_sd_basket <- function(
    sd_val,
    total_sample_size,
    individual_sd,
    recruitment_sd,
    stratify_by_syndrome,
    goal,                 # "fp" or "power"
    sim_file              # path to the appropriate simulate_*_basket_trial.R
) {
  # Load the right simulator INSIDE each worker (once per task)
  # Instead of source(sim_file), do an if/else:
  if (sim_file == "simulate superiority basket trial.R") {
    source("simulate superiority basket trial.R")
    message("Loaded superiority simulator")
  } else if (sim_file == "simulate non-inferiority basket trial.R") {
    source("simulate non-inferiority basket trial.R")
    message("Loaded non-inferiority simulator")
  } else {
    stop("sim_file must be 'superiority' or 'non-inferiority'")
  }
  
  threshold_met <- list(
    naive_t_test = NA_integer_,
    fix_effect_model_test = NA_integer_,
    mixed_effect_model_test = NA_integer_,
    mantel_haenszel_test = NA_integer_
  )
  
  for (site_num in 2:100) {
    result_df <- simulate_multiple_trials(
      num_sites = site_num,
      specification = "random_selection",
      iterations = 2500,
      total_sample_size = total_sample_size,
      generation_method = "distribution",
      individual_sd = individual_sd,
      treatment_sd  = sd_val,
      recruitment_sd = recruitment_sd,
      syndrome_num = 2,
      syndrome_prop_var = 0.05,                  # even allocation
      syndrome_effect = c(0.1, -0.1),            # syndrome effect
      syndrome_effect_variation = c(0.05, 0.05), # fixed per site
      stratify_by_syndrome = stratify_by_syndrome
    )
    
    powers <- colMeans(
      result_df[, c("naive_t_test","fix_effect_model_test","mixed_effect_model_test","mantel_haenszel_test")],
      na.rm = TRUE
    )
    
    for (test_name in names(threshold_met)) {
      if (is.na(threshold_met[[test_name]]) && !is.na(powers[test_name])) {
        if (goal == "fp"    && powers[test_name] < 0.05) threshold_met[[test_name]] <- site_num
        if (goal == "power" && powers[test_name] > 0.90) threshold_met[[test_name]] <- site_num
      }
    }
    if (all(!is.na(unlist(threshold_met)))) break
  }
  
  data.frame(
    treatment_sd = sd_val,
    naive_t_test = threshold_met$naive_t_test,
    fix_effect_model_test = threshold_met$fix_effect_model_test,
    mixed_effect_model_test = threshold_met$mixed_effect_model_test,
    mantel_haenszel_test = threshold_met$mantel_haenszel_test
  )
}

# ----- Utility to run a whole scenario in parallel over an sd list -----
run_scenario <- function(sd_list, total_n, indiv_sd, recruit_sd, stratify, goal, sim_file, out_file) {
  message(sprintf(
    "Running %s: total_n=%d, indiv_sd=%.3f, recruit_sd=%.3f, stratify=%s, |sd|=%d",
    goal, total_n, indiv_sd, recruit_sd, stratify, length(sd_list)
  ))
  
  # Use exactly one worker per sd (won't exceed PBS_NCPUS)
  ncores_avail <- as.integer(Sys.getenv("PBS_NCPUS", unset = length(sd_list)))
  workers <- min(ncores_avail, length(sd_list))
  plan(multisession, workers = workers)
  
  res_list <- future_lapply(
    sd_list,
    function(sd) run_one_sd_basket(
      sd_val = sd,
      total_sample_size = total_n,
      individual_sd = indiv_sd,
      recruitment_sd = recruit_sd,
      stratify_by_syndrome = stratify,
      goal = goal,
      sim_file = sim_file
    ),
    future.seed = TRUE
  )
  
  out <- do.call(rbind, res_list)
  print(out)
  saveRDS(out, file = out_file)
  invisible(out)
}

# ======================
# sd lists
# ======================
sd_list_NI <- c(0.005, 0.01, 0.015, 0.02, 0.025, 0.035, 0.05, 0.075, 0.1)
sd_list_SUP <- c(0.005, 0.01, 0.015, 0.02, 0.025, 0.035, 0.05, 0.075, 0.1)

# ======================
# Non-inferiority (you passed goal = "power" in your snippet)
# ======================
run_scenario(
  sd_list = sd_list_NI,
  total_n = 824, indiv_sd = 0.05, recruit_sd = 0.0,
  stratify = TRUE, goal = "power",
  sim_file = "simulate non-inferiority basket trial.R",
  out_file = "non-inferiority_basket_stratify_treatmentsd_individualsd=0.05_recruitmentsd=0.0.rds"
)

run_scenario(
  sd_list = sd_list_NI,
  total_n = 824, indiv_sd = 0.10, recruit_sd = 0.10,
  stratify = TRUE, goal = "power",
  sim_file = "simulate non-inferiority basket trial.R",
  out_file = "non-inferiority_basket_stratify_treatmentsd_individualsd=0.1_recruitmentsd=0.1.rds"
)

# ======================
# Superiority (you passed goal = "fp" in your snippet)
# ======================
run_scenario(
  sd_list = sd_list_SUP,
  total_n = 840, indiv_sd = 0.05, recruit_sd = 0.0,
  stratify = TRUE, goal = "fp",
  sim_file = "simulate superiority basket trial.R",
  out_file = "superiority_basket_stratify_treatmentsd_individualsd=0.05_recruitmentsd=0.0.rds"
)

run_scenario(
  sd_list = sd_list_SUP,
  total_n = 840, indiv_sd = 0.10, recruit_sd = 0.10,
  stratify = TRUE, goal = "fp",
  sim_file = "simulate superiority basket trial.R",
  out_file = "superiority_basket_stratify_treatmentsd_individualsd=0.1_recruitmentsd=0.1.rds"
)
