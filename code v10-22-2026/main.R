# load in files
source("simulate trial.R")
set.seed(123)

##################### run simulations (simple, non-basket)
run_main_simulation <- function(
    total_sample_size    = 964,
    max_sites            = 5,
    iterations           = 10,
    generation_method    = "distribution",  # or "bootstrap"
    treatment_proportion = 0.5,
    individual_sd        = 0.01,
    treatment_sd         = 0.03,
    recruitment_sd       = 0.04,
    # --- NEW knobs to reflect the combined simple simulator ---
    design_type          = c("superiority","non-inferiority"),
    mu                   = 0.10,            # set per scenario later
    min_mortality        = 0.25,            # e.g., 0.25–0.75 (naive) or 0.20–0.60 (40%)
    max_mortality        = 0.75,
    alpha                = 0.05,
    ni_margin            = 0.10,
    save_trial_data      = TRUE
) {
  design_type <- match.arg(design_type)
  
  results_list <- list()
  
  for (num_sites in 2:max_sites) {
    cat("Running simulation with", num_sites, "sites...\n")
    
    trial_results <- simulate_multiple_trials(
      num_sites              = num_sites,
      specification          = "random_selection",
      iterations             = iterations,
      total_sample_size      = total_sample_size,
      treatment_proportions  = treatment_proportion,
      generation_method      = generation_method,
      min_mortality          = min_mortality,
      max_mortality          = max_mortality,
      individual_sd          = individual_sd,
      mu                     = mu,
      treatment_sd           = treatment_sd,
      recruitment_sd         = recruitment_sd,
      design_type            = design_type,
      alpha                  = alpha,
      ni_margin              = ni_margin,
      save_trial_data        = save_trial_data
    )
    
    # store one row per site count, same structure as before
    row <- data.frame(
      num_sites = num_sites,
      stringsAsFactors = FALSE
    )
    row$simulation_data <- list(trial_results)
    results_list[[length(results_list) + 1]] <- row
  }
  
  results_df <- do.call(rbind, results_list)
  rownames(results_df) <- NULL
  results_df
}


# =========================
# Hyperparameters (easy to tweak)
# =========================
MAX_SITES   <- 3
ITERATIONS  <- 5
GEN_METHOD  <- "distribution"
TRT_PROP    <- 0.5
ALPHA       <- 0.05
NI_MARGIN   <- 0.10

# Mortality bands
SUP_MIN_MORT <- 0.25; SUP_MAX_MORT <- 0.75
NI_MIN_MORT  <- 0.20; NI_MAX_MORT  <- 0.60

# Mean treatment effects (mu)
SUP_MU_POWER <-  0.10  # superiority power
SUP_MU_FP    <-  0.00  # superiority false positive
NI_MU_POWER  <-  0.00  # non-inferiority power
NI_MU_FP     <- -0.10  # non-inferiority false positive

# Sample sizes
SUP_SS_POWER <- 840
SUP_SS_FP    <- 840
NI_SS_POWER  <- 824
NI_SS_FP     <- 824

# SD cases in the order: (individual_sd, treatment_sd, recruitment_sd)
CASE1_ind_sd <- 0.10; CASE1_trt_sd <- 0.00; CASE1_rec_sd <- 0.00
CASE2_ind_sd <- 0.05; CASE2_trt_sd <- 0.10; CASE2_rec_sd <- 0.00
CASE3_ind_sd <- 0.05; CASE3_trt_sd <- 0.00; CASE3_rec_sd <- 0.10
CASE4_ind_sd <- 0.05; CASE4_trt_sd <- 0.10; CASE4_rec_sd <- 0.10

# =========================
# Output folders
# =========================
dir.create("superiority_power",       showWarnings = FALSE, recursive = TRUE)
dir.create("superiority_fp",          showWarnings = FALSE, recursive = TRUE)
dir.create("non-inferiority_power",   showWarnings = FALSE, recursive = TRUE)
dir.create("non-inferiority_fp",      showWarnings = FALSE, recursive = TRUE)

set.seed(123)

# =========================
# SUPERIORITY — POWER (mu = +0.10), mortality 0.25–0.75
# =========================


results_df <- run_main_simulation(
  total_sample_size    = SUP_SS_POWER,
  max_sites            = MAX_SITES,
  iterations           = ITERATIONS,
  generation_method    = GEN_METHOD,
  treatment_proportion = TRT_PROP,
  individual_sd        = CASE1_ind_sd, treatment_sd = CASE1_trt_sd, recruitment_sd = CASE1_rec_sd,
  design_type          = "superiority",
  mu                   = SUP_MU_POWER,
  min_mortality        = SUP_MIN_MORT, max_mortality = SUP_MAX_MORT,
  alpha                = ALPHA
)
saveRDS(results_df, file = "superiority_power/sup_power_case1.rds")

results_df <- run_main_simulation(
  total_sample_size    = SUP_SS_POWER,
  max_sites            = MAX_SITES,
  iterations           = ITERATIONS,
  generation_method    = GEN_METHOD,
  treatment_proportion = TRT_PROP,
  individual_sd        = CASE2_ind_sd, treatment_sd = CASE2_trt_sd, recruitment_sd = CASE2_rec_sd,
  design_type          = "superiority",
  mu                   = SUP_MU_POWER,
  min_mortality        = SUP_MIN_MORT, max_mortality = SUP_MAX_MORT,
  alpha                = ALPHA
)
saveRDS(results_df, file = "superiority_power/sup_power_case2.rds")

results_df <- run_main_simulation(
  total_sample_size    = SUP_SS_POWER,
  max_sites            = MAX_SITES,
  iterations           = ITERATIONS,
  generation_method    = GEN_METHOD,
  treatment_proportion = TRT_PROP,
  individual_sd        = CASE3_ind_sd, treatment_sd = CASE3_trt_sd, recruitment_sd = CASE3_rec_sd,
  design_type          = "superiority",
  mu                   = SUP_MU_POWER,
  min_mortality        = SUP_MIN_MORT, max_mortality = SUP_MAX_MORT,
  alpha                = ALPHA
)
saveRDS(results_df, file = "superiority_power/sup_power_case3.rds")

results_df <- run_main_simulation(
  total_sample_size    = SUP_SS_POWER,
  max_sites            = MAX_SITES,
  iterations           = ITERATIONS,
  generation_method    = GEN_METHOD,
  treatment_proportion = TRT_PROP,
  individual_sd        = CASE4_ind_sd, treatment_sd = CASE4_trt_sd, recruitment_sd = CASE4_rec_sd,
  design_type          = "superiority",
  mu                   = SUP_MU_POWER,
  min_mortality        = SUP_MIN_MORT, max_mortality = SUP_MAX_MORT,
  alpha                = ALPHA
)
saveRDS(results_df, file = "superiority_power/sup_power_case4.rds")

# =========================
# SUPERIORITY — FALSE POSITIVE (mu = 0.00), mortality 0.25–0.75
# =========================

results_df <- run_main_simulation(
  total_sample_size    = SUP_SS_FP,
  max_sites            = MAX_SITES,
  iterations           = ITERATIONS,
  generation_method    = GEN_METHOD,
  treatment_proportion = TRT_PROP,
  individual_sd        = CASE1_ind_sd, treatment_sd = CASE1_trt_sd, recruitment_sd = CASE1_rec_sd,
  design_type          = "superiority",
  mu                   = SUP_MU_FP,
  min_mortality        = SUP_MIN_MORT, max_mortality = SUP_MAX_MORT,
  alpha                = ALPHA
)
saveRDS(results_df, file = "superiority_fp/sup_fp_case1.rds")

results_df <- run_main_simulation(
  total_sample_size    = SUP_SS_FP,
  max_sites            = MAX_SITES,
  iterations           = ITERATIONS,
  generation_method    = GEN_METHOD,
  treatment_proportion = TRT_PROP,
  individual_sd        = CASE2_ind_sd, treatment_sd = CASE2_trt_sd, recruitment_sd = CASE2_rec_sd,
  design_type          = "superiority",
  mu                   = SUP_MU_FP,
  min_mortality        = SUP_MIN_MORT, max_mortality = SUP_MAX_MORT,
  alpha                = ALPHA
)
saveRDS(results_df, file = "superiority_fp/sup_fp_case2.rds")

results_df <- run_main_simulation(
  total_sample_size    = SUP_SS_FP,
  max_sites            = MAX_SITES,
  iterations           = ITERATIONS,
  generation_method    = GEN_METHOD,
  treatment_proportion = TRT_PROP,
  individual_sd        = CASE3_ind_sd, treatment_sd = CASE3_trt_sd, recruitment_sd = CASE3_rec_sd,
  design_type          = "superiority",
  mu                   = SUP_MU_FP,
  min_mortality        = SUP_MIN_MORT, max_mortality = SUP_MAX_MORT,
  alpha                = ALPHA
)
saveRDS(results_df, file = "superiority_fp/sup_fp_case3.rds")

results_df <- run_main_simulation(
  total_sample_size    = SUP_SS_FP,
  max_sites            = MAX_SITES,
  iterations           = ITERATIONS,
  generation_method    = GEN_METHOD,
  treatment_proportion = TRT_PROP,
  individual_sd        = CASE4_ind_sd, treatment_sd = CASE4_trt_sd, recruitment_sd = CASE4_rec_sd,
  design_type          = "superiority",
  mu                   = SUP_MU_FP,
  min_mortality        = SUP_MIN_MORT, max_mortality = SUP_MAX_MORT,
  alpha                = ALPHA
)
saveRDS(results_df, file = "superiority_fp/sup_fp_case4.rds")

# =========================
# NON-INFERIORITY — POWER (mu = 0.00), mortality 0.20–0.60, N = 840
# =========================

results_df <- run_main_simulation(
  total_sample_size    = NI_SS_POWER,
  max_sites            = MAX_SITES,
  iterations           = ITERATIONS,
  generation_method    = GEN_METHOD,
  treatment_proportion = TRT_PROP,
  individual_sd        = CASE1_ind_sd, treatment_sd = CASE1_trt_sd, recruitment_sd = CASE1_rec_sd,
  design_type          = "non-inferiority",
  mu                   = NI_MU_POWER,
  min_mortality        = NI_MIN_MORT, max_mortality = NI_MAX_MORT,
  alpha                = ALPHA, ni_margin = NI_MARGIN
)
saveRDS(results_df, file = "non-inferiority_power/ni_power_case1.rds")

results_df <- run_main_simulation(
  total_sample_size    = NI_SS_POWER,
  max_sites            = MAX_SITES,
  iterations           = ITERATIONS,
  generation_method    = GEN_METHOD,
  treatment_proportion = TRT_PROP,
  individual_sd        = CASE2_ind_sd, treatment_sd = CASE2_trt_sd, recruitment_sd = CASE2_rec_sd,
  design_type          = "non-inferiority",
  mu                   = NI_MU_POWER,
  min_mortality        = NI_MIN_MORT, max_mortality = NI_MAX_MORT,
  alpha                = ALPHA, ni_margin = NI_MARGIN
)
saveRDS(results_df, file = "non-inferiority_power/ni_power_case2.rds")

results_df <- run_main_simulation(
  total_sample_size    = NI_SS_POWER,
  max_sites            = MAX_SITES,
  iterations           = ITERATIONS,
  generation_method    = GEN_METHOD,
  treatment_proportion = TRT_PROP,
  individual_sd        = CASE3_ind_sd, treatment_sd = CASE3_trt_sd, recruitment_sd = CASE3_rec_sd,
  design_type          = "non-inferiority",
  mu                   = NI_MU_POWER,
  min_mortality        = NI_MIN_MORT, max_mortality = NI_MAX_MORT,
  alpha                = ALPHA, ni_margin = NI_MARGIN
)
saveRDS(results_df, file = "non-inferiority_power/ni_power_case3.rds")

results_df <- run_main_simulation(
  total_sample_size    = NI_SS_POWER,
  max_sites            = MAX_SITES,
  iterations           = ITERATIONS,
  generation_method    = GEN_METHOD,
  treatment_proportion = TRT_PROP,
  individual_sd        = CASE4_ind_sd, treatment_sd = CASE4_trt_sd, recruitment_sd = CASE4_rec_sd,
  design_type          = "non-inferiority",
  mu                   = NI_MU_POWER,
  min_mortality        = NI_MIN_MORT, max_mortality = NI_MAX_MORT,
  alpha                = ALPHA, ni_margin = NI_MARGIN
)
saveRDS(results_df, file = "non-inferiority_power/ni_power_case4.rds")

# =========================
# NON-INFERIORITY — FALSE POSITIVE (mu = -0.10), mortality 0.20–0.60, N = 824
# =========================

results_df <- run_main_simulation(
  total_sample_size    = NI_SS_FP,
  max_sites            = MAX_SITES,
  iterations           = ITERATIONS,
  generation_method    = GEN_METHOD,
  treatment_proportion = TRT_PROP,
  individual_sd        = CASE1_ind_sd, treatment_sd = CASE1_trt_sd, recruitment_sd = CASE1_rec_sd,
  design_type          = "non-inferiority",
  mu                   = NI_MU_FP,
  min_mortality        = NI_MIN_MORT, max_mortality = NI_MAX_MORT,
  alpha                = ALPHA, ni_margin = NI_MARGIN
)
saveRDS(results_df, file = "non-inferiority_fp/ni_fp_case1.rds")

results_df <- run_main_simulation(
  total_sample_size    = NI_SS_FP,
  max_sites            = MAX_SITES,
  iterations           = ITERATIONS,
  generation_method    = GEN_METHOD,
  treatment_proportion = TRT_PROP,
  individual_sd        = CASE2_ind_sd, treatment_sd = CASE2_trt_sd, recruitment_sd = CASE2_rec_sd,
  design_type          = "non-inferiority",
  mu                   = NI_MU_FP,
  min_mortality        = NI_MIN_MORT, max_mortality = NI_MAX_MORT,
  alpha                = ALPHA, ni_margin = NI_MARGIN
)
saveRDS(results_df, file = "non-inferiority_fp/ni_fp_case2.rds")

results_df <- run_main_simulation(
  total_sample_size    = NI_SS_FP,
  max_sites            = MAX_SITES,
  iterations           = ITERATIONS,
  generation_method    = GEN_METHOD,
  treatment_proportion = TRT_PROP,
  individual_sd        = CASE3_ind_sd, treatment_sd = CASE3_trt_sd, recruitment_sd = CASE3_rec_sd,
  design_type          = "non-inferiority",
  mu                   = NI_MU_FP,
  min_mortality        = NI_MIN_MORT, max_mortality = NI_MAX_MORT,
  alpha                = ALPHA, ni_margin = NI_MARGIN
)
saveRDS(results_df, file = "non-inferiority_fp/ni_fp_case3.rds")

results_df <- run_main_simulation(
  total_sample_size    = NI_SS_FP,
  max_sites            = MAX_SITES,
  iterations           = ITERATIONS,
  generation_method    = GEN_METHOD,
  treatment_proportion = TRT_PROP,
  individual_sd        = CASE4_ind_sd, treatment_sd = CASE4_trt_sd, recruitment_sd = CASE4_rec_sd,
  design_type          = "non-inferiority",
  mu                   = NI_MU_FP,
  min_mortality        = NI_MIN_MORT, max_mortality = NI_MAX_MORT,
  alpha                = ALPHA, ni_margin = NI_MARGIN
)
saveRDS(results_df, file = "non-inferiority_fp/ni_fp_case4.rds")
