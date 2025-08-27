# individual_parallel.R
# Parallel outer-loop over individual_sd values
options(parallelly.availableCores.methods = c("env", "pbs", "system", "fallback"))
options(parallelly.cgroups.mounts.parse = FALSE)
Sys.setenv(R_PARALLELLY_AVAILABLECORES_METHODS = "env,pbs,system,fallback")

cat("=== R starting ===\n")

library(future)
library(future.apply)


# ---- One task = run ONE individual_sd value with early-stop ----
run_one_sd <- function(sd_val, total_n, criterion, sim_file) {
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
      individual_sd = sd_val,
      treatment_sd = 0,
      recruitment_sd = 0
    )
    
    powers <- colMeans(
      result_df[, c("naive_t_test","fix_effect_model_test","mixed_effect_model_test","mantel_haenszel_test")],
      na.rm = TRUE
    )
    
    for (nm in names(threshold_met)) {
      if (is.na(threshold_met[[nm]]) && !is.na(powers[nm])) {
        if (criterion == "fp<0.05"    && powers[nm] < 0.05) threshold_met[[nm]] <- site_num
        if (criterion == "power>0.90" && powers[nm] > 0.90) threshold_met[[nm]] <- site_num
      }
    }
    
    if (all(!is.na(unlist(threshold_met)))) break
  }
  
  data.frame(
    individual_sd = sd_val,
    naive_t_test = threshold_met$naive_t_test,
    fix_effect_model_test = threshold_met$fix_effect_model_test,
    mixed_effect_model_test = threshold_met$mixed_effect_model_test,
    mantel_haenszel_test = threshold_met$mantel_haenszel_test
  )
}

# ---- Run a whole scenario in parallel ----
run_scenario <- function(indiv_sd_list, total_n, criterion, sim_file, out_file) {
  n_avail <- as.integer(Sys.getenv("PBS_NCPUS", unset = length(indiv_sd_list)))
  workers <- min(n_avail, length(indiv_sd_list))
  plan(multisession, workers = workers)
  
  message(sprintf(
    "Scenario: total_n=%d, sim=%s, criterion=%s, workers=%d, |sd|=%d",
    total_n, sim_file, criterion, workers, length(indiv_sd_list)
  ))
  
  res_list <- future_lapply(
    indiv_sd_list,
    function(sd) run_one_sd(sd, total_n, criterion, sim_file),
    future.seed = TRUE
  )
  
  out <- do.call(rbind, res_list)
  print(out)
  saveRDS(out, out_file)
  invisible(out)
}

# ======================
# Define the grid
# ======================
indiv_sd_list <- c(0.005, 0.01, 0.015, 0.02, 0.025, 0.035, 0.05, 0.075,
                        0.1, 0.15, 0.2, 0.3)

# ======================
# Superiority (criterion: false positive < 0.05)
# ======================
run_scenario(
  indiv_sd_list,
  total_n = 840,
  criterion = "fp<0.05",
  sim_file = "simulate trial.R",
  out_file = "individualsd_treatmentsd=0_recruitmentsd=0.rds"
)

# ======================
# Non-inferiority (criterion: power > 0.90)
# ======================
run_scenario(
  indiv_sd_list,
  total_n = 824,
  criterion = "power>0.90",
  sim_file = "simulate trial non-inferiority.R",
  out_file = "non-inferiority_individualsd_treatmentsd=0_recruitmentsd=0.rds"
)
