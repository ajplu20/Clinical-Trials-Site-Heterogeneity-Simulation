# basket_recruitment_parallel.R
# Parallel outer-loop over recruitment_sd values (basket trials)
options(parallelly.availableCores.methods = c("env", "pbs", "system", "fallback"))
options(parallelly.cgroups.mounts.parse = FALSE)
Sys.setenv(R_PARALLELLY_AVAILABLECORES_METHODS = "env,pbs,system,fallback")


cat("=== R starting ===\n")
library(future)
library(future.apply)


# ---- One task = run ONE recruitment_sd with early-stop preserved ----
run_one_recruitment_basket <- function(
    recruit_sd,               # swept
    total_sample_size,        # 840 (superiority) or 824 (non-inferiority)
    individual_sd,            # fixed (e.g., 0.05)
    treatment_sd,             # fixed (e.g., 0.00)
    stratify_by_syndrome,
    goal,                     # "fp" (<0.05) or "power" (>0.90)
    sim_file                  # "simulate superiority basket trial.R" or "simulate non-inferiority basket trial.R"
) {
  # Load the right simulator inside the worker so simulate_multiple_trials() is the correct one
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
  
  # Your original code uses 2:100 for superiority, 2:100 or 2:30 for NI; we can use 2:100 safely
  for (site_num in 2:100) {
    result_df <- simulate_multiple_trials(
      num_sites = site_num,
      specification = "random_selection",
      iterations = 2500,
      total_sample_size = total_sample_size,
      generation_method = "distribution",
      individual_sd = individual_sd,     # fixed
      treatment_sd  = treatment_sd,      # fixed
      recruitment_sd = recruit_sd,       # swept
      syndrome_num = 2,
      syndrome_prop_var = 0.05,
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
        if (goal == "fp"    && powers[nm] < 0.05) threshold_met[[nm]] <- site_num
        if (goal == "power" && powers[nm] > 0.90) threshold_met[[nm]] <- site_num
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
run_scenario_recruitment <- function(recruitment_sd_list, total_n, individual_sd, treatment_sd,
                                     stratify, goal, sim_file, out_file) {
  # One worker per sd (capped by PBS_NCPUS)
  n_avail <- as.integer(Sys.getenv("PBS_NCPUS", unset = length(recruitment_sd_list)))
  workers <- min(n_avail, length(recruitment_sd_list))
  plan(multisession, workers = workers)
  
  message(sprintf(
    "Scenario: sim=%s goal=%s total_n=%d indiv_sd=%.3f treat_sd=%.3f stratify=%s workers=%d |sd|=%d",
    sim_file, goal, total_n, individual_sd, treatment_sd, stratify, workers, length(recruitment_sd_list)
  ))
  
  res_list <- future_lapply(
    recruitment_sd_list,
    function(rs) run_one_recruitment_basket(
      recruit_sd = rs,
      total_sample_size = total_n,
      individual_sd = individual_sd,
      treatment_sd  = treatment_sd,
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
# Grid
# ======================
recruitment_sd_list <- c(0.005, 0.01, 0.015, 0.02, 0.025, 0.035, 0.05, 0.075, 0.1, 0.15, 0.2, 0.3)

# ======================
# Superiority (find FP): goal="fp", sim = superiority basket
# ======================
run_scenario_recruitment(
  recruitment_sd_list,
  total_n = 840,
  individual_sd = 0.05,
  treatment_sd = 0.00,
  stratify = TRUE,
  goal = "fp",
  sim_file = "simulate superiority basket trial.R",
  out_file = "superiority_basket_stratify_recruitmentsd_individualsd=0.05_treatmentsd=0.0.rds"
)

# ======================
# Non-inferiority (find POWER): goal="power", sim = non-inferiority basket
# ======================
run_scenario_recruitment(
  recruitment_sd_list,
  total_n = 824,
  individual_sd = 0.05,
  treatment_sd = 0.00,
  stratify = TRUE,
  goal = "power",
  sim_file = "simulate non-inferiority basket trial.R",
  out_file = "non-inferiority_basket_stratify_recruitmentsd_individualsd=0.05_treatmentsd=0.0.rds"
)
