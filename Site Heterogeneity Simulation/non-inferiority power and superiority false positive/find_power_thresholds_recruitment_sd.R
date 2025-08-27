# recruitment_parallel.R
# Parallel outer-loop over recruitment_sd values (non-inferiority & superiority)
options(parallelly.availableCores.methods = c("env", "pbs", "system", "fallback"))
options(parallelly.cgroups.mounts.parse = FALSE)
Sys.setenv(R_PARALLELLY_AVAILABLECORES_METHODS = "env,pbs,system,fallback")

cat("=== R starting ===\n")

library(future)
library(future.apply)


# ---- One task = run ONE recruitment_sd value end-to-end (keeps early-stop) ----
run_one_recruitment <- function(
    recruit_sd,            # swept value
    total_n,               # 824 (NI) or 840 (Super)
    individual_sd,         # fixed (e.g., 0.05)
    treatment_sd,          # fixed (e.g., 0)
    max_site,              # 30 (NI in your code) or 100 (Super)
    goal,                  # "power" (>0.90) or "fp" (<0.05)
    sim_file               # simulator to source: NI or Super
) {
  # Load the correct simulator inside each worker
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
  
  for (site_num in 2:max_site) {
    result_df <- simulate_multiple_trials(
      num_sites = site_num,
      specification = "random_selection",
      iterations = 2500,
      total_sample_size = total_n,
      generation_method = "distribution",
      individual_sd = individual_sd,   # fixed
      treatment_sd  = treatment_sd,    # fixed
      recruitment_sd = recruit_sd      # swept
    )
    
    powers <- colMeans(
      result_df[, c("naive_t_test","fix_effect_model_test","mixed_effect_model_test","mantel_haenszel_test")],
      na.rm = TRUE
    )
    
    for (nm in names(threshold_met)) {
      if (is.na(threshold_met[[nm]]) && !is.na(powers[nm])) {
        if (goal == "power" && powers[nm] > 0.90) threshold_met[[nm]] <- site_num
        if (goal == "fp"    && powers[nm] < 0.05) threshold_met[[nm]] <- site_num
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

# ---- Run a whole scenario in parallel ----
run_scenario_recruitment <- function(recruitment_sd_list, total_n, individual_sd, treatment_sd,
                                     max_site, goal, sim_file, out_file) {
  # One worker per recruitment_sd (bounded by PBS_NCPUS if present)
  n_avail <- as.integer(Sys.getenv("PBS_NCPUS", unset = length(recruitment_sd_list)))
  workers <- min(n_avail, length(recruitment_sd_list))
  plan(multisession, workers = workers)
  
  message(sprintf(
    "Scenario: sim=%s goal=%s total_n=%d indiv_sd=%.3f treat_sd=%.3f max_site=%d workers=%d |sd|=%d",
    sim_file, goal, total_n, individual_sd, treatment_sd, max_site, workers, length(recruitment_sd_list)
  ))
  
  res_list <- future_lapply(
    recruitment_sd_list,
    function(rs) run_one_recruitment(
      recruit_sd = rs,
      total_n = total_n,
      individual_sd = individual_sd,
      treatment_sd  = treatment_sd,
      max_site = max_site,
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
# Non-inferiority (find POWER): goal="power", sim_file = NI
# Your original NI code used max site_num = 30
# ======================
run_scenario_recruitment(
  recruitment_sd_list,
  total_n = 824,
  individual_sd = 0.05,
  treatment_sd = 0.0,
  max_site = 30,
  goal = "power",
  sim_file = "simulate trial non-inferiority.R",
  out_file = "non-inferiority_recruitmentsd_individualsd=0.05_treatmentsd=0.0.rds"
)

# ======================
# Superiority (find FP): goal="fp", sim_file = Super
# Your original Super code used max site_num = 100
# ======================
run_scenario_recruitment(
  recruitment_sd_list,
  total_n = 840,
  individual_sd = 0.05,
  treatment_sd = 0.0,
  max_site = 100,
  goal = "fp",
  sim_file = "simulate trial.R",
  out_file = "recruitmentsd_individualsd=0.05_treatmentsd=0.0.rds"
)
