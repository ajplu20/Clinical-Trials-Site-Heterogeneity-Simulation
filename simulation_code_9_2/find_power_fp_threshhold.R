# ================================
# Simple (non-basket) sweep runner
# ================================

source("simulate trial.R")  # uses simulate_multiple_trials(...)

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

# ---- Run a whole scenario over a vector of values ----
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
ITERATIONS   <- 1
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
# NON-INFERIORITY — POWER
# =========================

ni_power_ind   <- run_scenario_sweep("individual_sd", SWEEP_VALUES, NI_SS_POWER, ITERATIONS, 0,0,0, "non-inferiority","power>0.90", mu=NI_MU_POWER, min_mortality=NI_MIN_MORT, max_mortality=NI_MAX_MORT)
saveRDS(ni_power_ind, "non-inferiority_power/individualsd_caseA.rds")

ni_power_trt_a <- run_scenario_sweep("treatment_sd", SWEEP_VALUES, NI_SS_POWER, ITERATIONS, 0.05,0,0, "non-inferiority","power>0.90", mu=NI_MU_POWER, min_mortality=NI_MIN_MORT, max_mortality=NI_MAX_MORT)
saveRDS(ni_power_trt_a, "non-inferiority_power/treatmentsd_caseA.rds")

ni_power_trt_b <- run_scenario_sweep("treatment_sd", SWEEP_VALUES, NI_SS_POWER, ITERATIONS, 0.05,0,0.10, "non-inferiority","power>0.90", mu=NI_MU_POWER, min_mortality=NI_MIN_MORT, max_mortality=NI_MAX_MORT)
saveRDS(ni_power_trt_b, "non-inferiority_power/treatmentsd_caseB.rds")

ni_power_rec   <- run_scenario_sweep("recruitment_sd", SWEEP_VALUES, NI_SS_POWER, ITERATIONS, 0.05,0,0, "non-inferiority","power>0.90", mu=NI_MU_POWER, min_mortality=NI_MIN_MORT, max_mortality=NI_MAX_MORT)
saveRDS(ni_power_rec, "non-inferiority_power/recruitsd_caseA.rds")

# =========================
# NON-INFERIORITY — FP
# =========================

ni_fp_ind   <- run_scenario_sweep("individual_sd", SWEEP_VALUES, NI_SS_FP, ITERATIONS, 0,0,0, "non-inferiority","fp<0.05", mu=NI_MU_FP, min_mortality=NI_MIN_MORT, max_mortality=NI_MAX_MORT)
saveRDS(ni_fp_ind, "non-inferiority_fp/individualsd_caseA.rds")

ni_fp_trt_a <- run_scenario_sweep("treatment_sd", SWEEP_VALUES, NI_SS_FP, ITERATIONS, 0.05,0,0, "non-inferiority","fp<0.05", mu=NI_MU_FP, min_mortality=NI_MIN_MORT, max_mortality=NI_MAX_MORT)
saveRDS(ni_fp_trt_a, "non-inferiority_fp/treatmentsd_caseA.rds")

ni_fp_trt_b <- run_scenario_sweep("treatment_sd", SWEEP_VALUES, NI_SS_FP, ITERATIONS, 0.05,0,0.10, "non-inferiority","fp<0.05", mu=NI_MU_FP, min_mortality=NI_MIN_MORT, max_mortality=NI_MAX_MORT)
saveRDS(ni_fp_trt_b, "non-inferiority_fp/treatmentsd_caseB.rds")

ni_fp_rec   <- run_scenario_sweep("recruitment_sd", SWEEP_VALUES, NI_SS_FP, ITERATIONS, 0.05,0,0, "non-inferiority","fp<0.05", mu=NI_MU_FP, min_mortality=NI_MIN_MORT, max_mortality=NI_MAX_MORT)
saveRDS(ni_fp_rec, "non-inferiority_fp/recruitsd_caseA.rds")

# =========================
# SUPERIORITY — POWER
# =========================

sup_power_ind   <- run_scenario_sweep("individual_sd", SWEEP_VALUES, SUP_SS_POWER, ITERATIONS, 0,0,0, "superiority","power>0.90", mu=SUP_MU_POWER, min_mortality=SUP_MIN_MORT, max_mortality=SUP_MAX_MORT)
saveRDS(sup_power_ind, "superiority_power/individualsd_caseA.rds")

sup_power_trt_a <- run_scenario_sweep("treatment_sd", SWEEP_VALUES, SUP_SS_POWER, ITERATIONS, 0.05,0,0, "superiority","power>0.90", mu=SUP_MU_POWER, min_mortality=SUP_MIN_MORT, max_mortality=SUP_MAX_MORT)
saveRDS(sup_power_trt_a, "superiority_power/treatmentsd_caseA.rds")

sup_power_trt_b <- run_scenario_sweep("treatment_sd", SWEEP_VALUES, SUP_SS_POWER, ITERATIONS, 0.05,0,0.10, "superiority","power>0.90", mu=SUP_MU_POWER, min_mortality=SUP_MIN_MORT, max_mortality=SUP_MAX_MORT)
saveRDS(sup_power_trt_b, "superiority_power/treatmentsd_caseB.rds")

sup_power_rec   <- run_scenario_sweep("recruitment_sd", SWEEP_VALUES, SUP_SS_POWER, ITERATIONS, 0.05,0,0, "superiority","power>0.90", mu=SUP_MU_POWER, min_mortality=SUP_MIN_MORT, max_mortality=SUP_MAX_MORT)
saveRDS(sup_power_rec, "superiority_power/recruitsd_caseA.rds")

# =========================
# SUPERIORITY — FP
# =========================

sup_fp_ind   <- run_scenario_sweep("individual_sd", SWEEP_VALUES, SUP_SS_FP, ITERATIONS, 0,0,0, "superiority","fp<0.05", mu=SUP_MU_FP, min_mortality=SUP_MIN_MORT, max_mortality=SUP_MAX_MORT)
saveRDS(sup_fp_ind, "superiority_fp/individualsd_caseA.rds")

sup_fp_trt_a <- run_scenario_sweep("treatment_sd", SWEEP_VALUES, SUP_SS_FP, ITERATIONS, 0.05,0,0, "superiority","fp<0.05", mu=SUP_MU_FP, min_mortality=SUP_MIN_MORT, max_mortality=SUP_MAX_MORT)
saveRDS(sup_fp_trt_a, "superiority_fp/treatmentsd_caseA.rds")

sup_fp_trt_b <- run_scenario_sweep("treatment_sd", SWEEP_VALUES, SUP_SS_FP, ITERATIONS, 0.05,0,0.10, "superiority","fp<0.05", mu=SUP_MU_FP, min_mortality=SUP_MIN_MORT, max_mortality=SUP_MAX_MORT)
saveRDS(sup_fp_trt_b, "superiority_fp/treatmentsd_caseB.rds")

sup_fp_rec   <- run_scenario_sweep("recruitment_sd", SWEEP_VALUES, SUP_SS_FP, ITERATIONS, 0.05,0,0, "superiority","fp<0.05", mu=SUP_MU_FP, min_mortality=SUP_MIN_MORT, max_mortality=SUP_MAX_MORT)
saveRDS(sup_fp_rec, "superiority_fp/recruitsd_caseA.rds")
