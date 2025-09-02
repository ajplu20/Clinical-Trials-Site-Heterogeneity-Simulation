# ================================
# Parallel Sweep Runner (Non-Basket) — HPC-safe & portable
# ================================

# ---- Deterministic core detection (env -> PBS -> system -> fallback) ----
options(parallelly.availableCores.methods = c("env", "pbs", "system", "fallback"))
options(parallelly.cgroups.mounts.parse = FALSE)
Sys.setenv(R_PARALLELLY_AVAILABLECORES_METHODS = "env,pbs,system,fallback")

# ---- Keep math libs single-threaded per worker to avoid oversubscription ----
Sys.setenv(
  OMP_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1"
)

library(future)
library(future.apply)
library(parallelly)

source("simulate trial.R")  # provides simulate_multiple_trials()

# ---- Run ONE value for the chosen sweep parameter end-to-end ----
run_one_value <- function(
    sweep_which = c("treatment_sd","individual_sd","recruitment_sd"),
    sweep_value,
    total_n, iterations = 1000,
    individual_sd = 0.05, treatment_sd = 0.00, recruitment_sd = 0.00,
    design_type = c("superiority","non-inferiority"),
    criterion   = c("power>0.90", "fp<0.05"),
    alpha = 0.05, ni_margin = 0.10,
    mu, min_mortality, max_mortality,
    generation_method = "distribution",
    treatment_proportion = 0.5
) {
  sweep_which <- match.arg(sweep_which)
  design_type <- match.arg(design_type)
  criterion   <- match.arg(criterion)
  
  if (sweep_which == "treatment_sd")   treatment_sd  <- sweep_value
  if (sweep_which == "individual_sd")  individual_sd <- sweep_value
  if (sweep_which == "recruitment_sd") recruitment_sd <- sweep_value
  
  threshold_met <- list(
    naive_t_test = NA_integer_,
    fix_effect_model_test = NA_integer_,
    mixed_effect_model_test = NA_integer_,
    mantel_haenszel_test = NA_integer_
  )
  
  for (site_num in 2:100) {
    result_df <- simulate_multiple_trials(
      num_sites             = site_num,
      specification         = "random_selection",
      iterations            = iterations,
      total_sample_size     = total_n,
      treatment_proportions = treatment_proportion,
      generation_method     = generation_method,
      min_mortality         = min_mortality,
      max_mortality         = max_mortality,
      individual_sd         = individual_sd,
      mu                    = mu,
      treatment_sd          = treatment_sd,
      recruitment_sd        = recruitment_sd,
      design_type           = design_type,
      alpha                 = alpha,
      ni_margin             = ni_margin,
      save_trial_data       = FALSE
    )
    
    powers <- colMeans(
      result_df[, c("naive_t_test","fix_effect_model_test",
                    "mixed_effect_model_test","mantel_haenszel_test")],
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
  
  out <- data.frame(
    sweep_param = sweep_which,
    sweep_value = sweep_value,
    naive_t_test = threshold_met$naive_t_test,
    fix_effect_model_test = threshold_met$fix_effect_model_test,
    mixed_effect_model_test = threshold_met$mixed_effect_model_test,
    mantel_haenszel_test = threshold_met$mantel_haenszel_test,
    stringsAsFactors = FALSE
  )
  names(out)[names(out) == "sweep_value"] <- sweep_which
  out
}

# ---- Run a sweep over vector of values ----
run_scenario_sweep <- function(
    sweep_which = c("treatment_sd","individual_sd","recruitment_sd"),
    sweep_values,
    total_n, iterations = 1000,
    individual_sd = 0.05, treatment_sd = 0.00, recruitment_sd = 0.00,
    design_type = c("superiority","non-inferiority"),
    criterion   = c("power>0.90", "fp<0.05"),
    alpha = 0.05, ni_margin = 0.10,
    mu, min_mortality, max_mortality,
    generation_method = "distribution",
    treatment_proportion = 0.5
) {
  sweep_which <- match.arg(sweep_which)
  design_type <- match.arg(design_type)
  criterion   <- match.arg(criterion)
  
  rows <- lapply(
    sweep_values,
    function(val) run_one_value(
      sweep_which = sweep_which, sweep_value = val,
      total_n = total_n, iterations = iterations,
      individual_sd = individual_sd, treatment_sd = treatment_sd, recruitment_sd = recruitment_sd,
      design_type = design_type, criterion = criterion,
      alpha = alpha, ni_margin = ni_margin,
      mu = mu, min_mortality = min_mortality, max_mortality = max_mortality,
      generation_method = generation_method,
      treatment_proportion = treatment_proportion
    )
  )
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

# =========================
# Hyperparameters
# =========================
ITERATIONS   <- 1000
SWEEP_VALUES <- c(0.005, 0.01, 0.015, 0.02, 0.025, 0.035, 0.05, 0.075, 0.10)

SUP_MIN_MORT <- 0.25; SUP_MAX_MORT <- 0.75
NI_MIN_MORT  <- 0.20; NI_MAX_MORT  <- 0.60

SUP_MU_POWER <-  0.10
SUP_MU_FP    <-  0.00
NI_MU_POWER  <-  0.00
NI_MU_FP     <- -0.10

SUP_SS_POWER <- 840
SUP_SS_FP    <- 840
NI_SS_POWER  <- 824
NI_SS_FP     <- 824

# =========================
# Define all 16 scenarios
# =========================
scenarios <- list(
  # NI Power
  list(file="non-inferiority_power/individualsd_caseA.rds", args=list("individual_sd", SWEEP_VALUES, NI_SS_POWER, ITERATIONS, 0,0,0, "non-inferiority","power>0.90", mu=NI_MU_POWER, min_mortality=NI_MIN_MORT, max_mortality=NI_MAX_MORT)),
  list(file="non-inferiority_power/treatmentsd_caseA.rds", args=list("treatment_sd", SWEEP_VALUES, NI_SS_POWER, ITERATIONS, 0.05,0,0, "non-inferiority","power>0.90", mu=NI_MU_POWER, min_mortality=NI_MIN_MORT, max_mortality=NI_MAX_MORT)),
  list(file="non-inferiority_power/treatmentsd_caseB.rds", args=list("treatment_sd", SWEEP_VALUES, NI_SS_POWER, ITERATIONS, 0.05,0,0.10, "non-inferiority","power>0.90", mu=NI_MU_POWER, min_mortality=NI_MIN_MORT, max_mortality=NI_MAX_MORT)),
  list(file="non-inferiority_power/recruitsd_caseA.rds",   args=list("recruitment_sd", SWEEP_VALUES, NI_SS_POWER, ITERATIONS, 0.05,0,0, "non-inferiority","power>0.90", mu=NI_MU_POWER, min_mortality=NI_MIN_MORT, max_mortality=NI_MAX_MORT)),
  
  # NI FP
  list(file="non-inferiority_fp/individualsd_caseA.rds", args=list("individual_sd", SWEEP_VALUES, NI_SS_FP, ITERATIONS, 0,0,0, "non-inferiority","fp<0.05", mu=NI_MU_FP, min_mortality=NI_MIN_MORT, max_mortality=NI_MAX_MORT)),
  list(file="non-inferiority_fp/treatmentsd_caseA.rds", args=list("treatment_sd", SWEEP_VALUES, NI_SS_FP, ITERATIONS, 0.05,0,0, "non-inferiority","fp<0.05", mu=NI_MU_FP, min_mortality=NI_MIN_MORT, max_mortality=NI_MAX_MORT)),
  list(file="non-inferiority_fp/treatmentsd_caseB.rds", args=list("treatment_sd", SWEEP_VALUES, NI_SS_FP, ITERATIONS, 0.05,0,0.10, "non-inferiority","fp<0.05", mu=NI_MU_FP, min_mortality=NI_MIN_MORT, max_mortality=NI_MAX_MORT)),
  list(file="non-inferiority_fp/recruitsd_caseA.rds",   args=list("recruitment_sd", SWEEP_VALUES, NI_SS_FP, ITERATIONS, 0.05,0,0, "non-inferiority","fp<0.05", mu=NI_MU_FP, min_mortality=NI_MIN_MORT, max_mortality=NI_MAX_MORT)),
  
  # SUP Power
  list(file="superiority_power/individualsd_caseA.rds", args=list("individual_sd", SWEEP_VALUES, SUP_SS_POWER, ITERATIONS, 0,0,0, "superiority","power>0.90", mu=SUP_MU_POWER, min_mortality=SUP_MIN_MORT, max_mortality=SUP_MAX_MORT)),
  list(file="superiority_power/treatmentsd_caseA.rds", args=list("treatment_sd", SWEEP_VALUES, SUP_SS_POWER, ITERATIONS, 0.05,0,0, "superiority","power>0.90", mu=SUP_MU_POWER, min_mortality=SUP_MIN_MORT, max_mortality=SUP_MAX_MORT)),
  list(file="superiority_power/treatmentsd_caseB.rds", args=list("treatment_sd", SWEEP_VALUES, SUP_SS_POWER, ITERATIONS, 0.05,0,0.10, "superiority","power>0.90", mu=SUP_MU_POWER, min_mortality=SUP_MIN_MORT, max_mortality=SUP_MAX_MORT)),
  list(file="superiority_power/recruitsd_caseA.rds",   args=list("recruitment_sd", SWEEP_VALUES, SUP_SS_POWER, ITERATIONS, 0.05,0,0, "superiority","power>0.90", mu=SUP_MU_POWER, min_mortality=SUP_MIN_MORT, max_mortality=SUP_MAX_MORT)),
  
  # SUP FP
  list(file="superiority_fp/individualsd_caseA.rds", args=list("individual_sd", SWEEP_VALUES, SUP_SS_FP, ITERATIONS, 0,0,0, "superiority","fp<0.05", mu=SUP_MU_FP, min_mortality=SUP_MIN_MORT, max_mortality=SUP_MAX_MORT)),
  list(file="superiority_fp/treatmentsd_caseA.rds", args=list("treatment_sd", SWEEP_VALUES, SUP_SS_FP, ITERATIONS, 0.05,0,0, "superiority","fp<0.05", mu=SUP_MU_FP, min_mortality=SUP_MIN_MORT, max_mortality=SUP_MAX_MORT)),
  list(file="superiority_fp/treatmentsd_caseB.rds", args=list("treatment_sd", SWEEP_VALUES, SUP_SS_FP, ITERATIONS, 0.05,0,0.10, "superiority","fp<0.05", mu=SUP_MU_FP, min_mortality=SUP_MIN_MORT, max_mortality=SUP_MAX_MORT)),
  list(file="superiority_fp/recruitsd_caseA.rds",   args=list("recruitment_sd", SWEEP_VALUES, SUP_SS_FP, ITERATIONS, 0.05,0,0, "superiority","fp<0.05", mu=SUP_MU_FP, min_mortality=SUP_MIN_MORT, max_mortality=SUP_MAX_MORT))
)

# =========================
# Run all scenarios in parallel (PBS-aware)
# =========================
n_avail <- parallelly::availableCores()  # auto-reads PBS_NCPUS on Vanda
workers <- min(n_avail, length(scenarios))
cat("Workers available:", n_avail, " | Using:", workers, "\n")

plan(multisession, workers = workers)
on.exit({ plan(sequential) }, add = TRUE)

results_list <- future_lapply(
  scenarios,
  function(sc) {
    message("Worker PID: ", Sys.getpid(),
            " | Host: ", Sys.info()[["nodename"]],
            " | File: ", sc$file)
    
    out <- do.call(run_scenario_sweep, sc$args)
    
    # Robust save (temp file then atomic rename)
    tmp <- paste0(sc$file, ".tmp_", Sys.getpid())
    saveRDS(out, tmp)
    file.rename(tmp, sc$file)
    
    data.frame(file = sc$file, rows = nrow(out), pid = Sys.getpid())
  },
  future.seed = 123   # reproducible, independent streams per iteration
)

results <- do.call(rbind, results_list)
print(results)
