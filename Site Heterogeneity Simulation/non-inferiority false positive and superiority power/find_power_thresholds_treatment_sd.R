# trial_treatmentsd_parallel.R
# Parallel outer-loop over treatment_sd values (non-basket trials)

## -------- Robust header: avoid cgroups bug & log environment --------
options(parallelly.availableCores.methods = c("env", "pbs", "system", "fallback"))
options(parallelly.cgroups.mounts.parse = FALSE)
Sys.setenv(R_PARALLELLY_AVAILABLECORES_METHODS = "env,pbs,system,fallback")

cat("=== R starting ===\n")

library(future)
library(future.apply)

# ---- Run ONE treatment_sd end-to-end (preserves early-stop) ----
run_one_treatmentsd <- function(
    sd_val, total_n, individual_sd, recruitment_sd,
    sim_file,              # "simulate trial non-inferiority.R" or "simulate trial.R"
    criterion              # "fp<0.05" (NI) or "power>0.90" (Super)
) {
  # Load the right simulator inside the worker so simulate_multiple_trials() is correct
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
      individual_sd = individual_sd,
      treatment_sd  = sd_val,
      recruitment_sd = recruitment_sd
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
    treatment_sd = sd_val,
    naive_t_test = threshold_met$naive_t_test,
    fix_effect_model_test = threshold_met$fix_effect_model_test,
    mixed_effect_model_test = threshold_met$mixed_effect_model_test,
    mantel_haenszel_test = threshold_met$mantel_haenszel_test
  )
}

# ---- Run a whole scenario in parallel over a vector of treatment_sd ----
run_scenario_treatmentsd <- function(treatment_sd_list, total_n, individual_sd, recruitment_sd,
                                     sim_file, criterion, out_file) {
  # One worker per treatment_sd (bounded by PBS_NCPUS)
  n_avail <- as.integer(Sys.getenv("PBS_NCPUS", unset = length(treatment_sd_list)))
  workers <- min(n_avail, length(treatment_sd_list))
  plan(multisession, workers = workers)
  
  message(sprintf(
    "Scenario: sim=%s criterion=%s total_n=%d indiv_sd=%.3f recruit_sd=%.3f workers=%d |sd|=%d",
    sim_file, criterion, total_n, individual_sd, recruitment_sd, workers, length(treatment_sd_list)
  ))
  
  res_list <- future_lapply(
    treatment_sd_list,
    function(sd)
      run_one_treatmentsd(
        sd_val = sd, total_n = total_n,
        individual_sd = individual_sd, recruitment_sd = recruitment_sd,
        sim_file = sim_file, criterion = criterion
      ),
    future.seed = TRUE    # independent, reproducible RNG per task
  )
  
  out <- do.call(rbind, res_list)
  print(out)
  saveRDS(out, file = out_file)
  invisible(out)
}

# ======================
# Your grid
# ======================
treatment_sd_list <- c(0.005, 0.01, 0.015, 0.02, 0.025, 0.035, 0.05, 0.075, 0.1)

# ======================
# Non-inferiority (FP < 0.05), simulator = NI
# ======================
run_scenario_treatmentsd(
  treatment_sd_list,
  total_n = 824, individual_sd = 0.05, recruitment_sd = 0.0,
  sim_file = "simulate trial non-inferiority.R",
  criterion = "fp<0.05",
  out_file = "non-inferiority_treatmentsd_individualsd=0.05_recruitmentsd=0.0.rds"
)

run_scenario_treatmentsd(
  treatment_sd_list,
  total_n = 824, individual_sd = 0.10, recruitment_sd = 0.10,
  sim_file = "simulate trial non-inferiority.R",
  criterion = "fp<0.05",
  out_file = "non-inferiority_treatmentsd_individualsd=0.1_recruitmentsd=0.1.rds"
)

# ======================
# Superiority (Power > 0.90), simulator = Super
# ======================
run_scenario_treatmentsd(
  treatment_sd_list,
  total_n = 840, individual_sd = 0.05, recruitment_sd = 0.0,
  sim_file = "simulate trial.R",
  criterion = "power>0.90",
  out_file = "superiority_treatmentsd_individualsd=0.05_recruitmentsd=0.0.rds"
)

run_scenario_treatmentsd(
  treatment_sd_list,
  total_n = 840, individual_sd = 0.10, recruitment_sd = 0.10,
  sim_file = "simulate trial.R",
  criterion = "power>0.90",
  out_file = "superiority_treatmentsd_individualsd=0.1_recruitmentsd=0.1.rds"
)
