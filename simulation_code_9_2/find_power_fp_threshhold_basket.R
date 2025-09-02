# ================================
# Basket sweep runner
# ================================

source("simulate trial basket.R")  # expects simulate_multiple_trials(...)

# ---- Run ONE value for the chosen sweep parameter end-to-end ----
run_one_value_basket <- function(
    sweep_which = c("treatment_sd","individual_sd","recruitment_sd"),
    sweep_value,
    total_n, iterations = 1000,
    individual_sd = 0.05, treatment_sd = 0.00, recruitment_sd = 0.00,
    design_type = c("superiority","non-inferiority"),
    criterion   = c("power>0.90", "fp<0.05"),
    alpha = 0.05, ni_margin = 0.10,
    mu, min_mortality, max_mortality,
    syndrome_num = 2,
    syndrome_prop_var = 0.05,
    syndrome_effect = c(0.1, -0.1),
    syndrome_effect_variation = c(0.05, 0.05),
    stratify_by_syndrome = FALSE,
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
      syndrome_num          = syndrome_num,
      syndrome_prop_var     = syndrome_prop_var,
      syndrome_effect       = syndrome_effect,
      syndrome_effect_variation = syndrome_effect_variation,
      stratify_by_syndrome  = stratify_by_syndrome,
      save_trial_data       = FALSE,
      save_syndrome_frames  = FALSE
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

# ---- Run a whole scenario over a vector ----
run_scenario_sweep_basket <- function(
    sweep_which, sweep_values,
    total_n, iterations = 1000,
    individual_sd = 0.05, treatment_sd = 0.00, recruitment_sd = 0.00,
    design_type = c("superiority","non-inferiority"),
    criterion   = c("power>0.90", "fp<0.05"),
    alpha = 0.05, ni_margin = 0.10,
    mu, min_mortality, max_mortality,
    syndrome_num = 2,
    syndrome_prop_var = 0.05,
    syndrome_effect = c(0.1, -0.1),
    syndrome_effect_variation = c(0.05, 0.05),
    stratify_by_syndrome = FALSE
) {
  rows <- lapply(
    sweep_values,
    function(val) run_one_value_basket(
      sweep_which = sweep_which, sweep_value = val,
      total_n = total_n, iterations = iterations,
      individual_sd = individual_sd, treatment_sd = treatment_sd, recruitment_sd = recruitment_sd,
      design_type = design_type, criterion = criterion,
      alpha = alpha, ni_margin = ni_margin,
      mu = mu, min_mortality = min_mortality, max_mortality = max_mortality,
      syndrome_num = syndrome_num,
      syndrome_prop_var = syndrome_prop_var,
      syndrome_effect = syndrome_effect,
      syndrome_effect_variation = syndrome_effect_variation,
      stratify_by_syndrome = stratify_by_syndrome
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

# Basket-specific knobs
BASKET_SYNDROME_NUM              <- 2
BASKET_SYNDROME_PROP_VAR         <- 0.05
BASKET_SYNDROME_EFFECT           <- c(0.1, -0.1)
BASKET_SYNDROME_EFFECT_VARIATION <- c(0.05, 0.05)
BASKET_STRATIFY_BY_SYNDROME      <- TRUE

# ========================= 
# NON-INFERIORITY — POWER
# =========================

ni_power_ind <- run_scenario_sweep_basket(
  sweep_which   = "individual_sd", sweep_values  = SWEEP_VALUES,
  total_n       = NI_SS_POWER, iterations = ITERATIONS,
  individual_sd = 0, treatment_sd = 0, recruitment_sd = 0,
  design_type   = "non-inferiority", criterion = "power>0.90",
  mu            = NI_MU_POWER, min_mortality = NI_MIN_MORT, max_mortality = NI_MAX_MORT,
  syndrome_num  = BASKET_SYNDROME_NUM, syndrome_prop_var = BASKET_SYNDROME_PROP_VAR,
  syndrome_effect = BASKET_SYNDROME_EFFECT, syndrome_effect_variation = BASKET_SYNDROME_EFFECT_VARIATION,
  stratify_by_syndrome = BASKET_STRATIFY_BY_SYNDROME
)
saveRDS(ni_power_ind, "non-inferiority_power/bas_individualsd_caseA.rds")

ni_power_trt_a <- run_scenario_sweep_basket(
  sweep_which   = "treatment_sd", sweep_values  = SWEEP_VALUES,
  total_n       = NI_SS_POWER, iterations = ITERATIONS,
  individual_sd = 0.05, treatment_sd = 0, recruitment_sd = 0,
  design_type   = "non-inferiority", criterion = "power>0.90",
  mu            = NI_MU_POWER, min_mortality = NI_MIN_MORT, max_mortality = NI_MAX_MORT,
  syndrome_num  = BASKET_SYNDROME_NUM, syndrome_prop_var = BASKET_SYNDROME_PROP_VAR,
  syndrome_effect = BASKET_SYNDROME_EFFECT, syndrome_effect_variation = BASKET_SYNDROME_EFFECT_VARIATION,
  stratify_by_syndrome = BASKET_STRATIFY_BY_SYNDROME
)
saveRDS(ni_power_trt_a, "non-inferiority_power/bas_treatmentsd_caseA.rds")

ni_power_trt_b <- run_scenario_sweep_basket(
  sweep_which   = "treatment_sd", sweep_values  = SWEEP_VALUES,
  total_n       = NI_SS_POWER, iterations = ITERATIONS,
  individual_sd = 0.05, treatment_sd = 0, recruitment_sd = 0.10,
  design_type   = "non-inferiority", criterion = "power>0.90",
  mu            = NI_MU_POWER, min_mortality = NI_MIN_MORT, max_mortality = NI_MAX_MORT,
  syndrome_num  = BASKET_SYNDROME_NUM, syndrome_prop_var = BASKET_SYNDROME_PROP_VAR,
  syndrome_effect = BASKET_SYNDROME_EFFECT, syndrome_effect_variation = BASKET_SYNDROME_EFFECT_VARIATION,
  stratify_by_syndrome = BASKET_STRATIFY_BY_SYNDROME
)
saveRDS(ni_power_trt_b, "non-inferiority_power/bas_treatmentsd_caseB.rds")

ni_power_rec <- run_scenario_sweep_basket(
  sweep_which   = "recruitment_sd", sweep_values  = SWEEP_VALUES,
  total_n       = NI_SS_POWER, iterations = ITERATIONS,
  individual_sd = 0.05, treatment_sd = 0, recruitment_sd = 0,
  design_type   = "non-inferiority", criterion = "power>0.90",
  mu            = NI_MU_POWER, min_mortality = NI_MIN_MORT, max_mortality = NI_MAX_MORT,
  syndrome_num  = BASKET_SYNDROME_NUM, syndrome_prop_var = BASKET_SYNDROME_PROP_VAR,
  syndrome_effect = BASKET_SYNDROME_EFFECT, syndrome_effect_variation = BASKET_SYNDROME_EFFECT_VARIATION,
  stratify_by_syndrome = BASKET_STRATIFY_BY_SYNDROME
)
saveRDS(ni_power_rec, "non-inferiority_power/bas_recruitsd_caseA.rds")

# =========================
# NON-INFERIORITY — FP
# =========================

ni_fp_ind <- run_scenario_sweep_basket(
  sweep_which   = "individual_sd", sweep_values  = SWEEP_VALUES,
  total_n       = NI_SS_FP, iterations = ITERATIONS,
  individual_sd = 0, treatment_sd = 0, recruitment_sd = 0,
  design_type   = "non-inferiority", criterion = "fp<0.05",
  mu            = NI_MU_FP, min_mortality = NI_MIN_MORT, max_mortality = NI_MAX_MORT,
  syndrome_num  = BASKET_SYNDROME_NUM, syndrome_prop_var = BASKET_SYNDROME_PROP_VAR,
  syndrome_effect = BASKET_SYNDROME_EFFECT, syndrome_effect_variation = BASKET_SYNDROME_EFFECT_VARIATION,
  stratify_by_syndrome = BASKET_STRATIFY_BY_SYNDROME
)
saveRDS(ni_fp_ind, "non-inferiority_fp/bas_individualsd_caseA.rds")

ni_fp_trt_a <- run_scenario_sweep_basket(
  sweep_which   = "treatment_sd", sweep_values  = SWEEP_VALUES,
  total_n       = NI_SS_FP, iterations = ITERATIONS,
  individual_sd = 0.05, treatment_sd = 0, recruitment_sd = 0,
  design_type   = "non-inferiority", criterion = "fp<0.05",
  mu            = NI_MU_FP, min_mortality = NI_MIN_MORT, max_mortality = NI_MAX_MORT,
  syndrome_num  = BASKET_SYNDROME_NUM, syndrome_prop_var = BASKET_SYNDROME_PROP_VAR,
  syndrome_effect = BASKET_SYNDROME_EFFECT, syndrome_effect_variation = BASKET_SYNDROME_EFFECT_VARIATION,
  stratify_by_syndrome = BASKET_STRATIFY_BY_SYNDROME
)
saveRDS(ni_fp_trt_a, "non-inferiority_fp/bas_treatmentsd_caseA.rds")

ni_fp_trt_b <- run_scenario_sweep_basket(
  sweep_which   = "treatment_sd", sweep_values  = SWEEP_VALUES,
  total_n       = NI_SS_FP, iterations = ITERATIONS,
  individual_sd = 0.05, treatment_sd = 0, recruitment_sd = 0.10,
  design_type   = "non-inferiority", criterion = "fp<0.05",
  mu            = NI_MU_FP, min_mortality = NI_MIN_MORT, max_mortality = NI_MAX_MORT,
  syndrome_num  = BASKET_SYNDROME_NUM, syndrome_prop_var = BASKET_SYNDROME_PROP_VAR,
  syndrome_effect = BASKET_SYNDROME_EFFECT, syndrome_effect_variation = BASKET_SYNDROME_EFFECT_VARIATION,
  stratify_by_syndrome = BASKET_STRATIFY_BY_SYNDROME
)
saveRDS(ni_fp_trt_b, "non-inferiority_fp/bas_treatmentsd_caseB.rds")

ni_fp_rec <- run_scenario_sweep_basket(
  sweep_which   = "recruitment_sd", sweep_values  = SWEEP_VALUES,
  total_n       = NI_SS_FP, iterations = ITERATIONS,
  individual_sd = 0.05, treatment_sd = 0, recruitment_sd = 0,
  design_type   = "non-inferiority", criterion = "fp<0.05",
  mu            = NI_MU_FP, min_mortality = NI_MIN_MORT, max_mortality = NI_MAX_MORT,
  syndrome_num  = BASKET_SYNDROME_NUM, syndrome_prop_var = BASKET_SYNDROME_PROP_VAR,
  syndrome_effect = BASKET_SYNDROME_EFFECT, syndrome_effect_variation = BASKET_SYNDROME_EFFECT_VARIATION,
  stratify_by_syndrome = BASKET_STRATIFY_BY_SYNDROME
)
saveRDS(ni_fp_rec, "non-inferiority_fp/bas_recruitsd_caseA.rds")

# =========================
# SUPERIORITY — POWER
# =========================

sup_power_ind <- run_scenario_sweep_basket(
  sweep_which   = "individual_sd", sweep_values  = SWEEP_VALUES,
  total_n       = SUP_SS_POWER, iterations = ITERATIONS,
  individual_sd = 0, treatment_sd = 0, recruitment_sd = 0,
  design_type   = "superiority", criterion = "power>0.90",
  mu            = SUP_MU_POWER, min_mortality = SUP_MIN_MORT, max_mortality = SUP_MAX_MORT,
  syndrome_num  = BASKET_SYNDROME_NUM, syndrome_prop_var = BASKET_SYNDROME_PROP_VAR,
  syndrome_effect = BASKET_SYNDROME_EFFECT, syndrome_effect_variation = BASKET_SYNDROME_EFFECT_VARIATION,
  stratify_by_syndrome = BASKET_STRATIFY_BY_SYNDROME
)
saveRDS(sup_power_ind, "superiority_power/bas_individualsd_caseA.rds")

sup_power_trt_a <- run_scenario_sweep_basket(
  sweep_which   = "treatment_sd", sweep_values  = SWEEP_VALUES,
  total_n       = SUP_SS_POWER, iterations = ITERATIONS,
  individual_sd = 0.05, treatment_sd = 0, recruitment_sd = 0,
  design_type   = "superiority", criterion = "power>0.90",
  mu            = SUP_MU_POWER, min_mortality = SUP_MIN_MORT, max_mortality = SUP_MAX_MORT,
  syndrome_num  = BASKET_SYNDROME_NUM, syndrome_prop_var = BASKET_SYNDROME_PROP_VAR,
  syndrome_effect = BASKET_SYNDROME_EFFECT, syndrome_effect_variation = BASKET_SYNDROME_EFFECT_VARIATION,
  stratify_by_syndrome = BASKET_STRATIFY_BY_SYNDROME
)
saveRDS(sup_power_trt_a, "superiority_power/bas_treatmentsd_caseA.rds")

sup_power_trt_b <- run_scenario_sweep_basket(
  sweep_which   = "treatment_sd", sweep_values  = SWEEP_VALUES,
  total_n       = SUP_SS_POWER, iterations = ITERATIONS,
  individual_sd = 0.05, treatment_sd = 0, recruitment_sd = 0.10,
  design_type   = "superiority", criterion = "power>0.90",
  mu            = SUP_MU_POWER, min_mortality = SUP_MIN_MORT, max_mortality = SUP_MAX_MORT,
  syndrome_num  = BASKET_SYNDROME_NUM, syndrome_prop_var = BASKET_SYNDROME_PROP_VAR,
  syndrome_effect = BASKET_SYNDROME_EFFECT, syndrome_effect_variation = BASKET_SYNDROME_EFFECT_VARIATION,
  stratify_by_syndrome = BASKET_STRATIFY_BY_SYNDROME
)
saveRDS(sup_power_trt_b, "superiority_power/bas_treatmentsd_caseB.rds")

sup_power_rec <- run_scenario_sweep_basket(
  sweep_which   = "recruitment_sd", sweep_values  = SWEEP_VALUES,
  total_n       = SUP_SS_POWER, iterations = ITERATIONS,
  individual_sd = 0.05, treatment_sd = 0, recruitment_sd = 0,
  design_type   = "superiority", criterion = "power>0.90",
  mu            = SUP_MU_POWER, min_mortality = SUP_MIN_MORT, max_mortality = SUP_MAX_MORT,
  syndrome_num  = BASKET_SYNDROME_NUM, syndrome_prop_var = BASKET_SYNDROME_PROP_VAR,
  syndrome_effect = BASKET_SYNDROME_EFFECT, syndrome_effect_variation = BASKET_SYNDROME_EFFECT_VARIATION,
  stratify_by_syndrome = BASKET_STRATIFY_BY_SYNDROME
)
saveRDS(sup_power_rec, "superiority_power/bas_recruitsd_caseA.rds")

# =========================
# SUPERIORITY — FP
# =========================

sup_fp_ind <- run_scenario_sweep_basket(
  sweep_which   = "individual_sd", sweep_values  = SWEEP_VALUES,
  total_n       = SUP_SS_FP, iterations = ITERATIONS,
  individual_sd = 0, treatment_sd = 0, recruitment_sd = 0,
  design_type   = "superiority", criterion = "fp<0.05",
  mu            = SUP_MU_FP, min_mortality = SUP_MIN_MORT, max_mortality = SUP_MAX_MORT,
  syndrome_num  = BASKET_SYNDROME_NUM, syndrome_prop_var = BASKET_SYNDROME_PROP_VAR,
  syndrome_effect = BASKET_SYNDROME_EFFECT, syndrome_effect_variation = BASKET_SYNDROME_EFFECT_VARIATION,
  stratify_by_syndrome = BASKET_STRATIFY_BY_SYNDROME
)
saveRDS(sup_fp_ind, "superiority_fp/bas_individualsd_caseA.rds")

sup_fp_trt_a <- run_scenario_sweep_basket(
  sweep_which   = "treatment_sd", sweep_values  = SWEEP_VALUES,
  total_n       = SUP_SS_FP, iterations = ITERATIONS,
  individual_sd = 0.05, treatment_sd = 0, recruitment_sd = 0,
  design_type   = "superiority", criterion = "fp<0.05",
  mu            = SUP_MU_FP, min_mortality = SUP_MIN_MORT, max_mortality = SUP_MAX_MORT,
  syndrome_num  = BASKET_SYNDROME_NUM, syndrome_prop_var = BASKET_SYNDROME_PROP_VAR,
  syndrome_effect = BASKET_SYNDROME_EFFECT, syndrome_effect_variation = BASKET_SYNDROME_EFFECT_VARIATION,
  stratify_by_syndrome = BASKET_STRATIFY_BY_SYNDROME
)
saveRDS(sup_fp_trt_a, "superiority_fp/bas_treatmentsd_caseA.rds")

sup_fp_trt_b <- run_scenario_sweep_basket(
  sweep_which   = "treatment_sd", sweep_values  = SWEEP_VALUES,
  total_n       = SUP_SS_FP, iterations = ITERATIONS,
  individual_sd = 0.05, treatment_sd = 0, recruitment_sd = 0.10,
  design_type   = "superiority", criterion = "fp<0.05",
  mu            = SUP_MU_FP, min_mortality = SUP_MIN_MORT, max_mortality = SUP_MAX_MORT,
  syndrome_num  = BASKET_SYNDROME_NUM, syndrome_prop_var = BASKET_SYNDROME_PROP_VAR,
  syndrome_effect = BASKET_SYNDROME_EFFECT, syndrome_effect_variation = BASKET_SYNDROME_EFFECT_VARIATION,
  stratify_by_syndrome = BASKET_STRATIFY_BY_SYNDROME
)
saveRDS(sup_fp_trt_b, "superiority_fp/bas_treatmentsd_caseB.rds")

sup_fp_rec <- run_scenario_sweep_basket(
  sweep_which   = "recruitment_sd", sweep_values  = SWEEP_VALUES,
  total_n       = SUP_SS_FP, iterations = ITERATIONS,
  individual_sd = 0.05, treatment_sd = 0, recruitment_sd = 0,
  design_type   = "superiority", criterion = "fp<0.05",
  mu            = SUP_MU_FP, min_mortality = SUP_MIN_MORT, max_mortality = SUP_MAX_MORT,
  syndrome_num  = BASKET_SYNDROME_NUM, syndrome_prop_var = BASKET_SYNDROME_PROP_VAR,
  syndrome_effect = BASKET_SYNDROME_EFFECT, syndrome_effect_variation = BASKET_SYNDROME_EFFECT_VARIATION,
  stratify_by_syndrome = BASKET_STRATIFY_BY_SYNDROME
)
saveRDS(sup_fp_rec, "superiority_fp/bas_recruitsd_caseA.rds")
