# individualsd_twoways_parallel.R
# Parallel outer-loop over individual_sd for Superiority (power) and Non-inferiority (FP)
options(parallelly.availableCores.methods = c("env", "pbs", "system", "fallback"))
options(parallelly.cgroups.mounts.parse = FALSE)
Sys.setenv(R_PARALLELLY_AVAILABLECORES_METHODS = "env,pbs,system,fallback")

cat("=== R starting ===\n")

library(future)
library(future.apply)

# ---- Run ONE individual_sd end-to-end (keeps your early-stop) ----
run_one_individual <- function(
    indiv_sd,            # swept value
    total_n,             # 840 (superiority) or 824 (non-inferiority)
    criterion,           # "power>0.90" or "fp<0.05"
    sim_file             # "simulate trial.R" or "simulate trial non-inferiority.R"
) {
  # Load the right simulator in the worker so simulate_multiple_trials() is correct

  # Use label instead of filename
  if (sim_file == "simulate trial.R") {
    source("simulate trial.R")
    message("Loaded superiority simulator")
  } else if (sim_file == "simulate trial non-inferiority.R") {
    source("simulate trial non-inferiority.R")
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
      total_sample_size = total_n,
      generation_method = "distribution",
      individual_sd = indiv_sd,  # swept
      treatment_sd = 0,
      recruitment_sd = 0
    )
    
    powers <- colMeans(
      result_df[, c("naive_t_test","fix_effect_model_test","mixed_effect_model_test","mantel_haenszel_test")],
      na.rm = TRUE
    )
    
    for (nm in names(threshold_met)) {
      if (is.na(threshold_met[[nm]]) && !is.na(powers[nm])) {
        if (criterion == "power>0.90" && powers[nm] > 0.90) threshold_met[[nm]] <- site_num
        if (criterion == "fp<0.05"    && powers[nm] < 0.05) threshold_met[[nm]] <- site_num
      }
    }
    
    if (all(!is.na(unlist(threshold_met)))) break
  }
  
  data.frame(
    individual_sd = indiv_sd,
    naive_t_test = threshold_met$naive_t_test,
    fix_effect_model_test = threshold_met$fix_effect_model_test,
    mixed_effect_model_test = threshold_met$mixed_effect_model_test,
    mantel_haenszel_test = threshold_met$mantel_haenszel_test
  )
}

# ---- Run a whole scenario in parallel over a vector of individual_sd ----
run_scenario_individual <- function(individual_sd_list, total_n, criterion, sim_file, out_file) {
  # One worker per sd (capped by PBS_NCPUS)
  n_avail <- as.integer(Sys.getenv("PBS_NCPUS", unset = length(individual_sd_list)))
  workers <- min(n_avail, length(individual_sd_list))
  plan(multisession, workers = workers)
  
  message(sprintf(
    "Scenario: sim=%s criterion=%s total_n=%d workers=%d |sd|=%d",
    sim_file, criterion, total_n, workers, length(individual_sd_list)
  ))
  
  res_list <- future_lapply(
    individual_sd_list,
    function(sd) run_one_individual(sd, total_n, criterion, sim_file),
    future.seed = TRUE  # independent, reproducible RNG per task
  )
  
  out <- do.call(rbind, res_list)
  print(out)
  saveRDS(out, file = out_file)
  invisible(out)
}

# ======================
# Grid
# ======================
individual_sd_list <- c(0.005, 0.01, 0.015, 0.02, 0.025, 0.035, 0.05, 0.075,
                        0.1, 0.15, 0.2, 0.3)

# ======================
# Superiority (POWER > 0.90), simulator = "simulate trial.R"
# ======================
run_scenario_individual(
  individual_sd_list,
  total_n = 840,
  criterion = "power>0.90",
  sim_file = "simulate trial.R",
  out_file = "individualsd_treatmentsd=0_recruitmentsd=0.rds"
)

# ======================
# Non-inferiority (FP < 0.05), simulator = "simulate trial non-inferiority.R"
# ======================
run_scenario_individual(
  individual_sd_list,
  total_n = 824,
  criterion = "fp<0.05",
  sim_file = "simulate trial non-inferiority.R",
  out_file = "non-inferiority_individualsd_treatmentsd=0_recruitmentsd=0.rds"
)
