# recruitmentsd_parallel.R
# Parallel outer-loop over recruitment_sd for NON-BASKET trials (NI + Super)
options(parallelly.availableCores.methods = c("env", "pbs", "system", "fallback"))
options(parallelly.cgroups.mounts.parse = FALSE)
Sys.setenv(R_PARALLELLY_AVAILABLECORES_METHODS = "env,pbs,system,fallback")

cat("=== R starting ===\n")

library(future)
library(future.apply)

# ---- Run ONE recruitment_sd end-to-end (keeps early-stop) ----
run_one_recruitmentsd <- function(
    recruit_sd,          # swept value
    total_n,             # 824 (NI) or 840 (Super)
    individual_sd,       # fixed (0.05 in your examples)
    treatment_sd,        # fixed (0 in your examples)
    criterion,           # "fp<0.05" (NI) or "power>0.90" (Super)
    sim_file             # "simulate trial non-inferiority.R" or "simulate trial.R"
) {
  
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
      individual_sd = individual_sd,  # fixed
      treatment_sd  = treatment_sd,   # fixed
      recruitment_sd = recruit_sd     # swept
    )
    
    powers <- colMeans(
      result_df[, c("naive_t_test","fix_effect_model_test","mixed_effect_model_test","mantel_haenszel_test")],
      na.rm = TRUE
    )
    
    for (nm in names(threshold_met)) {
      if (is.na(threshold_met[[nm]]) && !is.na(powers[nm])) {
        if (criterion == "fp<0.05"    && powers[nm] < 0.05)  threshold_met[[nm]] <- site_num
        if (criterion == "power>0.90" && powers[nm] > 0.90)  threshold_met[[nm]] <- site_num
      }
    }
    
    if (all(!is.na(unlist(threshold_met)))) break
  }
  
  data.frame(
    recruitment_sd = recruit_sd,
    naive_t_test = threshold_met$naive_t_test,
    fix_effect_model_test = threshold_met$fix_effect_model_test,
    mixed_effect_model_test = threshold_met$mixed_effect_model_test,
    mantel_haenszel_test = threshold_met$mantel_haenszel_test
  )
}

# ---- Run a whole scenario in parallel over a vector of recruitment_sd ----
run_scenario_recruitmentsd <- function(
    recruitment_sd_list, total_n, individual_sd, treatment_sd, criterion, sim_file, out_file
) {
  # One worker per sd (bounded by PBS_NCPUS if set)
  n_avail <- as.integer(Sys.getenv("PBS_NCPUS", unset = length(recruitment_sd_list)))
  workers <- min(n_avail, length(recruitment_sd_list))
  plan(multisession, workers = workers)
  
  message(sprintf(
    "Scenario: sim=%s criterion=%s total_n=%d indiv_sd=%.3f treat_sd=%.3f workers=%d |sd|=%d",
    sim_file, criterion, total_n, individual_sd, treatment_sd, workers, length(recruitment_sd_list)
  ))
  
  res_list <- future_lapply(
    recruitment_sd_list,
    function(rs) run_one_recruitmentsd(
      recruit_sd = rs,
      total_n = total_n,
      individual_sd = individual_sd,
      treatment_sd  = treatment_sd,
      criterion = criterion,
      sim_file = sim_file
    ),
    future.seed = TRUE  # independent, reproducible RNG substreams per task
  )
  
  out <- do.call(rbind, res_list)
  print(out)
  saveRDS(out, file = out_file)
  invisible(out)
}

# ======================
# Grid to sweep
# ======================
recruitment_sd_list <- c(0.005, 0.01, 0.015, 0.02, 0.025, 0.035, 0.05, 0.075,
                         0.1, 0.15, 0.2, 0.3)

# ======================
# Non-inferiority (FP < 0.05), simulator = NI file
# ======================
run_scenario_recruitmentsd(
  recruitment_sd_list,
  total_n = 824,
  individual_sd = 0.05,
  treatment_sd = 0.0,
  criterion = "fp<0.05",
  sim_file = "simulate trial non-inferiority.R",
  out_file = "non-inferiority_recruitmentsd_individualsd=0.05_treatmentsd=0.0.rds"
)

# ======================
# Superiority (Power > 0.90), simulator = Super file
# ======================
run_scenario_recruitmentsd(
  recruitment_sd_list,
  total_n = 840,
  individual_sd = 0.05,
  treatment_sd = 0.0,
  criterion = "power>0.90",
  sim_file = "simulate trial.R",
  out_file = "recruitmentsd_individualsd=0.05_treatmentsd=0.0.rds"
)
