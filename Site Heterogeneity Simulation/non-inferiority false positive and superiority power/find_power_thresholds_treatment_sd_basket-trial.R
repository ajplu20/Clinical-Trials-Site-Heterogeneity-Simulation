# basket_treatmentsd_parallel.R
# Parallel outer-loop over treatment_sd values (basket trials)
options(parallelly.availableCores.methods = c("env", "pbs", "system", "fallback"))
options(parallelly.cgroups.mounts.parse = FALSE)
Sys.setenv(R_PARALLELLY_AVAILABLECORES_METHODS = "env,pbs,system,fallback")

cat("=== R starting ===\n")

library(future)
library(future.apply)

# ---- Run ONE treatment_sd end-to-end (keeps early-stop) ----
run_one_treatmentsd_basket <- function(
    sd_val,
    total_sample_size,
    individual_sd,
    recruitment_sd,
    stratify_by_syndrome,
    goal,        # "fp" (<0.05) or "power" (>0.90)
    sim_file     # "simulate superiority basket trial.R" or "simulate non-inferiority basket trial.R"
) {
  # Load the correct simulator inside the worker so simulate_multiple_trials() is the right one
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
      syndrome_prop_var = 0.05,           # even allocation
      syndrome_effect = c(0.1, -0.1),     # syndrome effect
      syndrome_effect_variation = c(0.05, 0.05),
      stratify_by_syndrome = stratify_by_syndrome
    )
    
    powers <- colMeans(
      result_df[, c("naive_t_test","fix_effect_model_test","mixed_effect_model_test","mantel_haenszel_test")],
      na.rm = TRUE
    )
    
    for (nm in names(threshold_met)) {
      if (is.na(threshold_met[[nm]]) && !is.na(powers[nm])) {
        if (goal == "fp"    && powers[nm] < 0.05) threshold_met[[nm]] <- site_num
        if (goal == "power" && powers[nm] > 0.90) threshold_met[[nm]] <- site_num
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

# ---- Run a full scenario in parallel over a vector of treatment_sd ----
run_scenario_treatmentsd <- function(treatment_sd_list, total_n, individual_sd, recruitment_sd,
                                     stratify, goal, sim_file, out_file) {
  # One worker per sd, bounded by PBS_NCPUS if set
  n_avail <- as.integer(Sys.getenv("PBS_NCPUS", unset = length(treatment_sd_list)))
  workers <- min(n_avail, length(treatment_sd_list))
  plan(multisession, workers = workers)
  
  message(sprintf(
    "Scenario: sim=%s goal=%s total_n=%d indiv_sd=%.3f recruit_sd=%.3f stratify=%s workers=%d |sd|=%d",
    sim_file, goal, total_n, individual_sd, recruitment_sd, stratify, workers, length(treatment_sd_list)
  ))
  
  res_list <- future_lapply(
    treatment_sd_list,
    function(sd)
      run_one_treatmentsd_basket(
        sd_val = sd,
        total_sample_size = total_n,
        individual_sd = individual_sd,
        recruitment_sd = recruitment_sd,
        stratify_by_syndrome = stratify,
        goal = goal,
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
# Your grids
# ======================
treatment_sd_list <- c(0.005, 0.01, 0.015, 0.02, 0.025, 0.035, 0.05, 0.075, 0.1)

# ======================
# Non-inferiority runs (goal = "fp"), simulator = non-inferiority basket
# ======================
run_scenario_treatmentsd(
  treatment_sd_list,
  total_n = 824,
  individual_sd = 0.05, recruitment_sd = 0.0,
  stratify = TRUE, goal = "fp",
  sim_file = "simulate non-inferiority basket trial.R",
  out_file = "non-inferiority_basket_stratify_treatmentsd_individualsd=0.05_recruitmentsd=0.00.rds"
)

run_scenario_treatmentsd(
  treatment_sd_list,
  total_n = 824,
  individual_sd = 0.10, recruitment_sd = 0.10,
  stratify = TRUE, goal = "fp",
  sim_file = "simulate non-inferiority basket trial.R",
  out_file = "non-inferiority_basket_stratify_treatmentsd_individualsd=0.1_recruitmentsd=0.1.rds"
)

# ======================
# Superiority runs (goal = "power"), simulator = superiority basket
# ======================
run_scenario_treatmentsd(
  treatment_sd_list,
  total_n = 840,
  individual_sd = 0.05, recruitment_sd = 0.0,
  stratify = TRUE, goal = "power",
  sim_file = "simulate superiority basket trial.R",
  out_file = "superiority_basket_stratify_treatmentsd_individualsd=0.05_recruitmentsd=0.0.rds"
)

run_scenario_treatmentsd(
  treatment_sd_list,
  total_n = 840,
  individual_sd = 0.10, recruitment_sd = 0.10,
  stratify = TRUE, goal = "power",
  sim_file = "simulate superiority basket trial.R",
  out_file = "superiority_basket_stratify_treatmentsd_individualsd=0.1_recruitmentsd=0.1.rds"
)
