# for each value in the vector, the goal is the output a dataframe of the number of sites required to hit 90% power or 5% fp depending on user input
#the way to interpret the vector is that, rather than just plug it in for sweep value, if it is treatment sd, plug in for sweep value
#if it is recruitment or individual, we find the 95% interval for treatment sd, then find the coresponding recruitmentsd/individualsd with the same interval width
#then we run for that
# the number in first number column that originally showed sweep values is still the sweep value we enter, disregard transformation
#logic for each sweep value: for each site, simulate_multiple_trials, save the colmeans of that and the site number, run up to 100 sites,
#after that, fit a loess curve for each analysis method, treating colmean as y and site number as x (ie naive t, fixed, mixed, mantel haenzel), with stipulated parameters shown
#we create 95% confidence interval around the curve for each method, not prediction interval, meanign interval. for each method, the number we record down for the final dataframe is the first site for which
#the confidence interval's upper bound > 0.9 for the power case, or if the user input is false positive case, lower bound < 0.05. 

source("simulate trial.R")  # uses simulate_multiple_trials(...)

# =========================
# Hyperparameters
# =========================
ITERATIONS   <- 1000
SWEEP_VALUES <- c(0.02, 0.04, 0.06, 0.08, 0.1)

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
LOESS_SHOW_SE  <- FALSE        # used only for plotting UI parity; not used in CI computation below
LOESS_SPAN     <- 0.6
LOESS_DEGREE   <- 2
LOESS_FAMILY   <- "gaussian"   # keep "gaussian" to enable SE computation
LOESS_LINEWIDTH<- 1



# =========================
# Utilities
# =========================

# Fit loess on (x, y) and return a data.frame with fit and 95% CI for the mean curve
.fit_loess_with_ci <- function(x, y, span = LOESS_SPAN, degree = LOESS_DEGREE, family = LOESS_FAMILY) {
  df <- data.frame(x = x, y = y)
  # Order by x to avoid interpolation weirdness
  df <- df[order(df$x), , drop = FALSE]
  
  fit <- stats::loess(y ~ x, data = df, span = span, degree = degree, family = family, control = loess.control(surface = "interpolate"))
  pred <- predict(fit, newdata = data.frame(x = df$x), se = TRUE)  # se=TRUE gives SE for mean curve
  
  out <- data.frame(
    x     = df$x,
    fit   = as.numeric(pred$fit),
    sefit = as.numeric(pred$se.fit)
  )
  out$lower95 <- out$fit - 1.96 * out$sefit
  out$upper95 <- out$fit + 1.96 * out$sefit
  out
}

# Given a loess-with-CI data.frame and a threshold rule, find the first x meeting it
.first_crossing <- function(loess_df, criterion = c("power>0.90", "fp<0.05")) {
  criterion <- match.arg(criterion)
  if (criterion == "power>0.90") {
    idx <- which(loess_df$upper95 > 0.90)
  } else {
    idx <- which(loess_df$lower95 < 0.05)
  }
  if (length(idx) == 0) return(NA_integer_)
  as.integer(loess_df$x[min(idx)])
}

# Map SD for sweeps:
# - treatment_sd: use as-is
# - individual_sd / recruitment_sd: choose sigma so that
#   qnorm(0.975) * sigma = 0.95 * sweep_value
.mapped_sigma_from_95pct <- function(val) {
  if (is.na(val) || val <= 0) return(NA_real_)
  z975 <- qnorm(0.975)  # ~1.959964
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
  } else { # "recruitment_sd"
    sig <- .mapped_sigma_from_95pct(sweep_value)
    list(treatment_sd = NA_real_,
         individual_sd = NA_real_,
         recruitment_sd = sig)
  }
}


# =========================
# Run ONE sweep value end-to-end (new logic)
# =========================
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
    treatment_proportion = 0.5,
    max_sites = 50
) {
  sweep_which <- match.arg(sweep_which)
  design_type <- match.arg(design_type)
  criterion   <- match.arg(criterion)
  
  # Interpret the vector as requested:
  # - If sweeping treatment_sd: plug in directly.
  # - If sweeping individual/recruitment: choose SD so its 95% width matches (scalers adjustable).
  mapped <- .map_sd(sweep_which, sweep_value)
  # Overwrite only if specified by mapping; otherwise use incoming defaults
  if (!is.na(mapped$treatment_sd))  treatment_sd  <- mapped$treatment_sd
  if (!is.na(mapped$individual_sd)) individual_sd <- mapped$individual_sd
  if (!is.na(mapped$recruitment_sd)) recruitment_sd <- mapped$recruitment_sd
  
  # Collect per-site means for each analysis method
  methods <- c("naive_t_test","fix_effect_model_test","mixed_effect_model_test","mantel_haenszel_test")
  per_site <- setNames(vector("list", length(methods)), methods)
  for (m in methods) per_site[[m]] <- data.frame(site = integer(0), mean = numeric(0))
  
  # Simulate across sites 2..max_sites, store column means
  for (site_num in seq(2, max_sites, by = 2)) {
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
    
    col_means <- colMeans(result_df[, methods, drop = FALSE], na.rm = TRUE)
    for (m in methods) {
      per_site[[m]] <- rbind(per_site[[m]], data.frame(site = site_num, mean = as.numeric(col_means[[m]])))
    }
  }
  
  # Fit LOESS + CI for each method and compute first crossing per the criterion
  thresholds <- list()
  for (m in methods) {
    dfm <- per_site[[m]]
    # Handle degenerate case: not enough unique sites
    if (length(unique(dfm$site)) < 4) {
      thresholds[[m]] <- NA_integer_
      next
    }
    lo <- .fit_loess_with_ci(x = dfm$site, y = dfm$mean,
                             span = LOESS_SPAN, degree = LOESS_DEGREE, family = LOESS_FAMILY)
    thresholds[[m]] <- .first_crossing(lo, criterion = criterion)
  }
  
  # Same output format as before
  out <- data.frame(
    sweep_param = sweep_which,
    sweep_value = sweep_value,
    naive_t_test = thresholds[["naive_t_test"]],
    fix_effect_model_test = thresholds[["fix_effect_model_test"]],
    mixed_effect_model_test = thresholds[["mixed_effect_model_test"]],
    mantel_haenszel_test = thresholds[["mantel_haenszel_test"]],
    stringsAsFactors = FALSE
  )
  names(out)[names(out) == "sweep_value"] <- sweep_which  # keep the first numeric column name aligned with the sweep param
  out
}

# =========================
# Run a whole scenario over a vector of values (new logic)
# =========================
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
    treatment_proportion = 0.5,
    max_sites = 100
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
      treatment_proportion = treatment_proportion,
      max_sites = max_sites
    )
  )
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

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
