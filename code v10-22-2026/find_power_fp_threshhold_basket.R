# =========================================
# Basket Trial — LOESS + CI Crossing Logic
# =========================================

source("simulate trial basket.R")  # provides simulate_multiple_trials()

# -------------------------
# Hyperparameters (tunable)
# -------------------------
ITERATIONS   <- 1000
SWEEP_VALUES <- c(0.02, 0.04, 0.06, 0.08, 0.10)

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

# LOESS defaults (change here as needed)
LOESS_SHOW_SE   <- FALSE       # for plotting parity only; CI below uses predict(..., se=TRUE)
LOESS_SPAN      <- 0.6
LOESS_DEGREE    <- 2
LOESS_FAMILY    <- "gaussian"
LOESS_LINEWIDTH <- 1

# Basket knobs (tune as needed)
BASKET_SYNDROME_NUM              <- 2
BASKET_SYNDROME_PROP_VAR         <- 0.10
BASKET_SYNDROME_EFFECT           <- c(0.05, -0.05)
BASKET_SYNDROME_EFFECT_VARIATION <- c(0.025, 0.025)
BASKET_STRATIFY_BY_SYNDROME      <- FALSE

# --------------------------------
# Utilities (same as non-basket)
# --------------------------------

# LOESS fit + 95% CI for mean curve
.fit_loess_with_ci <- function(x, y,
                               span = LOESS_SPAN,
                               degree = LOESS_DEGREE,
                               family = LOESS_FAMILY) {
  df <- data.frame(x = x, y = y)
  df <- df[order(df$x), , drop = FALSE]
  fit <- stats::loess(y ~ x, data = df, span = span, degree = degree,
                      family = family, control = loess.control(surface = "interpolate"))
  pred <- predict(fit, newdata = data.frame(x = df$x), se = TRUE)
  out <- data.frame(
    x      = df$x,
    fit    = as.numeric(pred$fit),
    sefit  = as.numeric(pred$se.fit)
  )
  out$lower95 <- out$fit - 1.96 * out$sefit
  out$upper95 <- out$fit + 1.96 * out$sefit
  out
}

# First x where CI meets criterion
.first_crossing <- function(loess_df, criterion = c("power>0.90", "fp<0.05")) {
  criterion <- match.arg(criterion)
  idx <- if (criterion == "power>0.90") {
    which(loess_df$upper95 > 0.90)
  } else {
    which(loess_df$lower95 < 0.05)
  }
  if (length(idx) == 0) return(NA_integer_)
  as.integer(loess_df$x[min(idx)])
}

# Mapping: for non-treatment sweeps, choose sigma s.t. qnorm(0.975, 0, sigma) = 0.95 * sweep_value
.mapped_sigma_from_95pct <- function(val) {
  if (is.na(val) || val <= 0) return(NA_real_)
  z975 <- qnorm(0.975)
  (0.95 * val) / z975
}

.map_sd <- function(sweep_which, sweep_value) {
  if (sweep_which == "treatment_sd") {
    list(treatment_sd = sweep_value,
         individual_sd = NA_real_,
         recruitment_sd = NA_real_)
  } else if (sweep_which == "individual_sd") {
    sig <- .mapped_sigma_from_95pct(sweep_value)
    list(treatment_sd = NA_real_,
         individual_sd = sig,
         recruitment_sd = NA_real_)
  } else { # recruitment_sd
    sig <- .mapped_sigma_from_95pct(sweep_value)
    list(treatment_sd = NA_real_,
         individual_sd = NA_real_,
         recruitment_sd = sig)
  }
}

# --------------------------------------------------------
# Run ONE basket sweep value — LOESS + CI crossing logic
# --------------------------------------------------------
run_one_value_basket <- function(
    sweep_which = c("treatment_sd","individual_sd","recruitment_sd"),
    sweep_value,
    total_n, iterations = 1000,
    individual_sd = 0.05, treatment_sd = 0.00, recruitment_sd = 0.00,
    design_type = c("superiority","non-inferiority"),
    criterion   = c("power>0.90", "fp<0.05"),
    alpha = 0.05, ni_margin = 0.10,
    mu, min_mortality, max_mortality,
    syndrome_num = BASKET_SYNDROME_NUM,
    syndrome_prop_var = BASKET_SYNDROME_PROP_VAR,
    syndrome_effect = BASKET_SYNDROME_EFFECT,
    syndrome_effect_variation = BASKET_SYNDROME_EFFECT_VARIATION,
    stratify_by_syndrome = BASKET_STRATIFY_BY_SYNDROME,
    generation_method = "distribution",
    treatment_proportion = 0.5,
    site_grid = seq(2, 100, by = 2)  # match your original basket spacing
) {
  sweep_which <- match.arg(sweep_which)
  design_type <- match.arg(design_type)
  criterion   <- match.arg(criterion)
  
  # Apply mapping rule
  mapped <- .map_sd(sweep_which, sweep_value)
  if (!is.na(mapped$treatment_sd))  treatment_sd  <- mapped$treatment_sd
  if (!is.na(mapped$individual_sd)) individual_sd <- mapped$individual_sd
  if (!is.na(mapped$recruitment_sd)) recruitment_sd <- mapped$recruitment_sd
  
  methods <- c("naive_t_test","fix_effect_model_test","mixed_effect_model_test","mantel_haenszel_test")
  per_site <- setNames(vector("list", length(methods)), methods)
  for (m in methods) per_site[[m]] <- data.frame(site = integer(0), mean = numeric(0))
  
  # Simulate each site count in site_grid and collect column means
  for (site_num in site_grid) {
    result_df <- simulate_multiple_trials(
      num_sites                 = site_num,
      specification             = "random_selection",
      iterations                = iterations,
      total_sample_size         = total_n,
      treatment_proportions     = treatment_proportion,
      generation_method         = generation_method,
      min_mortality             = min_mortality,
      max_mortality             = max_mortality,
      individual_sd             = individual_sd,
      mu                        = mu,
      treatment_sd              = treatment_sd,
      recruitment_sd            = recruitment_sd,
      design_type               = design_type,
      alpha                     = alpha,
      ni_margin                 = ni_margin,
      syndrome_num              = syndrome_num,
      syndrome_prop_var         = syndrome_prop_var,
      syndrome_effect           = syndrome_effect,
      syndrome_effect_variation = syndrome_effect_variation,
      stratify_by_syndrome      = stratify_by_syndrome,
      save_trial_data           = FALSE,
      save_syndrome_frames      = FALSE
    )
    
    cm <- colMeans(result_df[, methods, drop = FALSE], na.rm = TRUE)
    for (m in methods) {
      per_site[[m]] <- rbind(per_site[[m]], data.frame(site = site_num, mean = as.numeric(cm[[m]])))
    }
  }
  
  # LOESS + CI crossing
  thresholds <- list()
  for (m in methods) {
    dfm <- per_site[[m]]
    if (length(unique(dfm$site)) < 4) { thresholds[[m]] <- NA_integer_; next }
    lo <- .fit_loess_with_ci(x = dfm$site, y = dfm$mean,
                             span = LOESS_SPAN, degree = LOESS_DEGREE, family = LOESS_FAMILY)
    thresholds[[m]] <- .first_crossing(lo, criterion = criterion)
  }
  
  # Output format identical to your original
  out <- data.frame(
    sweep_param = sweep_which,
    sweep_value = sweep_value,
    naive_t_test = thresholds[["naive_t_test"]],
    fix_effect_model_test = thresholds[["fix_effect_model_test"]],
    mixed_effect_model_test = thresholds[["mixed_effect_model_test"]],
    mantel_haenszel_test = thresholds[["mantel_haenszel_test"]],
    stringsAsFactors = FALSE
  )
  names(out)[names(out) == "sweep_value"] <- sweep_which
  out
}

# --------------------------------------------------------
# Run a basket sweep over vector of values (LOESS logic)
# --------------------------------------------------------
run_scenario_sweep_basket <- function(
    sweep_which, sweep_values,
    total_n, iterations = 1000,
    individual_sd = 0.05, treatment_sd = 0.00, recruitment_sd = 0.00,
    design_type = c("superiority","non-inferiority"),
    criterion   = c("power>0.90", "fp<0.05"),
    alpha = 0.05, ni_margin = 0.10,
    mu, min_mortality, max_mortality,
    syndrome_num = BASKET_SYNDROME_NUM,
    syndrome_prop_var = BASKET_SYNDROME_PROP_VAR,
    syndrome_effect = BASKET_SYNDROME_EFFECT,
    syndrome_effect_variation = BASKET_SYNDROME_EFFECT_VARIATION,
    stratify_by_syndrome = BASKET_STRATIFY_BY_SYNDROME,
    generation_method = "distribution",
    treatment_proportion = 0.5,
    site_grid = c(2, seq(5, 50, by = 5))
) {
  design_type <- match.arg(design_type)
  criterion   <- match.arg(criterion)
  
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
      stratify_by_syndrome = stratify_by_syndrome,
      generation_method = generation_method,
      treatment_proportion = treatment_proportion,
      site_grid = site_grid
    )
  )
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

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
