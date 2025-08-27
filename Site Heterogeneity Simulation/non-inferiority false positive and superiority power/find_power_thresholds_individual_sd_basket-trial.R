# basket_individualsd_parallel.R
# Parallel outer-loop over individual_sd values for basket trials
options(parallelly.availableCores.methods = c("env", "pbs", "system", "fallback"))
options(parallelly.cgroups.mounts.parse = FALSE)
Sys.setenv(R_PARALLELLY_AVAILABLECORES_METHODS = "env,pbs,system,fallback")

cat("=== R starting ===\n")

library(future)
library(future.apply)


# ---- Run ONE individual_sd end-to-end (preserves early-stop) ----
run_one_individualsd_basket <- function(
    indiv_sd,                 # swept value
    total_sample_size,        # 840 (superiority) or 824 (non-inferiority)
    treatment_sd,             # fixed (e.g., 0.00)
    recruitment_sd,           # fixed
    stratify_by_syndrome,     # TRUE/FALSE
    goal,                     # "fp" (<0.05) or "power" (>0.90)
    sim_file                  # "simulate superiority basket trial.R" or "simulate non-inferiority basket trial.R"
) {
  # Load the correct simulator inside the worker
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
      individual_sd = indiv_sd,         # <-- swept here
      treatment_sd  = treatment_sd,     # <-- fixed
      recruitment_sd = recruitment_sd,
      syndrome_num = 2,
      syndrome_prop_var = 0.05,         # even allocation
      syndrome_effect = c(0.1, -0.1),
      syndrome_effect_variation = c(0.05, 0.05),
      stratify_by_syndrome = stratify_by_syndrome
    )
    
    powers <- colMeans(
      result_df[, c("naive_t_test","fix_effect_model_test","mixed_effect_model_test","mantel_haenszel_test")],
      na.rm = TRUE
    )
    
    for (nm in names(threshold_met)) {
      if (is.na(threshold_met[[nm]]) && !is.na(powers[nm])) {
        if (goal == "fp"    && powers[nm] < 0.05)  threshold_met[[nm]] <- site_num
        if (goal == "power" && powers[nm] > 0.90)  threshold_met[[nm]] <- site_num
      }
    }
    
    if (all(!is.na(unlist(threshold_met)))) break
  }
  
  data.frame(
    individual_sd = indiv_sd,  # <-- correct label
    naive_t_test = threshold_met$naive_t_test,
    fix_effect_model_test = threshold_met$fix_effect_model_test,
    mixed_effect_model_test = threshold_met$mixed_effect_model_test,
    mantel_haenszel_test = threshold_met$mantel_haenszel_test
  )
}

# ---- Run a whole scenario in parallel over a vector of individual_sd ----
run_scenario_individualsd <- function(individual_sd_list, total_n, treatment_sd, recruitment_sd,
                                      stratify, goal, sim_file, out_file) {
  # One worker per sd (bounded by PBS_NCPUS if set)
  n_avail <- as.integer(Sys.getenv("PBS_NCPUS", unset = length(individual_sd_list)))
  workers <- min(n_avail, length(individual_sd_list))
  plan(multisession, workers = workers)
  
  message(sprintf(
    "Scenario: sim=%s goal=%s total_n=%d treat_sd=%.3f recruit_sd=%.3f stratify=%s workers=%d |sd|=%d",
    sim_file, goal, total_n, treatment_sd, recruitment_sd, stratify, workers, length(individual_sd_list)
  ))
  
  res_list <- future_lapply(
    individual_sd_list,
    function(sd)
      run_one_individualsd_basket(
        indiv_sd = sd,
        total_sample_size = total_n,
        treatment_sd = treatment_sd,
        recruitment_sd = recruitment_sd,
        stratify_by_syndrome = stratify,
        goal = goal,
        sim_file = sim_file
      ),
    future.seed = TRUE   # independent, reproducible RNG per task
  )
  
  out <- do.call(rbind, res_list)
  print(out)
  saveRDS(out, file = out_file)
  invisible(out)
}

# ======================
# Your grid (was named treatment_sd_list in your snippet, but it’s individual_sd values)
# ======================
individual_sd_list <- c(0.005, 0.01, 0.015, 0.02, 0.025, 0.035, 0.05, 0.075,
                        0.1, 0.15, 0.2, 0.3)

# ======================
# Superiority (you set goal = "power" in your call)
# ======================
run_scenario_individualsd(
  individual_sd_list,
  total_n = 840,
  treatment_sd = 0.00,
  recruitment_sd = 0.0,
  stratify = TRUE,
  goal = "power",
  sim_file = "simulate superiority basket trial.R",
  out_file = "superiority_basket_stratify_individualsd_treatmentsd=0.0_recruitmentsd=0.0.rds"
)

# ======================
# Non-inferiority (you set goal = "fp" in your call)
# ======================
run_scenario_individualsd(
  individual_sd_list,
  total_n = 824,
  treatment_sd = 0.00,
  recruitment_sd = 0.0,
  stratify = TRUE,
  goal = "fp",
  sim_file = "simulate non-inferiority basket trial.R",
  out_file = "non-inferiority_basket_stratify_individualsd_treatmentsd=0.0_recruitmentsd=0.0.rds"
)
