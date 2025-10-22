# ==========================================
# BASKET trial runner (PARALLEL per 4 cases)
# ==========================================

# load basket-capable simulator
source("simulate trial basket.R")

# ---- HPC-friendly parallel setup (same style as your working code) ----
options(parallelly.availableCores.methods = c("env", "pbs", "system", "fallback"))
options(parallelly.cgroups.mounts.parse = FALSE)
Sys.setenv(R_PARALLELLY_AVAILABLECORES_METHODS = "env,pbs,system,fallback")

suppressPackageStartupMessages({
  library(future)
  library(future.apply)
})

set.seed(123)

##################### run simulations (basket)
run_main_simulation_basket <- function(
    total_sample_size    = 964,
    max_sites            = 5,
    iterations           = 10,
    generation_method    = "distribution",  # or "bootstrap"
    treatment_proportion = 0.5,
    individual_sd        = 0.01,
    treatment_sd         = 0.03,
    recruitment_sd       = 0.04,
    # design knobs
    design_type          = c("superiority","non-inferiority"),
    mu                   = 0.10,
    min_mortality        = 0.25,
    max_mortality        = 0.75,
    alpha                = 0.05,
    ni_margin            = 0.10,
    save_trial_data      = TRUE,
    # ---- basket knobs ----
    syndrome_num = 2,
    syndrome_prop_var = 0.05,
    syndrome_effect = c(0.1, -0.1),
    syndrome_effect_variation = c(0.05, 0.05),
    stratify_by_syndrome = FALSE,
    save_syndrome_frames = TRUE
) {
  design_type <- match.arg(design_type)
  results_list <- list()
  
  for (num_sites in 2:max_sites) {
    cat("Running BASKET simulation with", num_sites, "sites...\n")
    
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
      # basket
      syndrome_num                 = syndrome_num,
      syndrome_prop_var            = syndrome_prop_var,
      syndrome_effect              = syndrome_effect,
      syndrome_effect_variation    = syndrome_effect_variation,
      stratify_by_syndrome         = stratify_by_syndrome,
      save_trial_data              = save_trial_data,
      save_syndrome_frames         = save_syndrome_frames
    )
    
    row <- data.frame(num_sites = num_sites, stringsAsFactors = FALSE)
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
MAX_SITES   <- 35
ITERATIONS  <- 2500
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

# Basket knobs (per your request)
BASKET_SYNDROME_NUM              <- 2
BASKET_SYNDROME_PROP_VAR         <- 0.05
BASKET_SYNDROME_EFFECT           <- c(0.1, -0.1)
BASKET_SYNDROME_EFFECT_VARIATION <- c(0.05, 0.05)
BASKET_STRATIFY_BY_SYNDROME      <- FALSE
BASKET_SAVE_SYNDROME_FRAMES      <- TRUE

# =========================
# Output folders (same four; filenames get a `basket_` prefix)
# =========================
#dir.create("superiority_power",       showWarnings = FALSE, recursive = TRUE)
#dir.create("superiority_fp",          showWarnings = FALSE, recursive = TRUE)
#dir.create("non-inferiority_power",   showWarnings = FALSE, recursive = TRUE)
#dir.create("non-inferiority_fp",      showWarnings = FALSE, recursive = TRUE)

# =========================
# Helper: run 4 SD cases IN PARALLEL for one basket scenario
# =========================
run_four_cases_parallel_basket <- function(
    folder, file_prefix,
    total_sample_size, design_type, mu,
    min_mortality, max_mortality,
    alpha = ALPHA, ni_margin = NI_MARGIN,
    max_sites = MAX_SITES, iterations = ITERATIONS,
    gen_method = GEN_METHOD, trt_prop = TRT_PROP,
    # basket knobs
    syndrome_num = BASKET_SYNDROME_NUM,
    syndrome_prop_var = BASKET_SYNDROME_PROP_VAR,
    syndrome_effect = BASKET_SYNDROME_EFFECT,
    syndrome_effect_variation = BASKET_SYNDROME_EFFECT_VARIATION,
    stratify_by_syndrome = BASKET_STRATIFY_BY_SYNDROME,
    save_syndrome_frames = BASKET_SAVE_SYNDROME_FRAMES,
    # parallel controls
    workers = 4,                            # 4 workers, one per case
    seed_base = 123
) {
  cases <- list(
    list(case_id=1, ind=CASE1_ind_sd, trt=CASE1_trt_sd, rec=CASE1_rec_sd),
    list(case_id=2, ind=CASE2_ind_sd, trt=CASE2_trt_sd, rec=CASE2_rec_sd),
    list(case_id=3, ind=CASE3_ind_sd, trt=CASE3_trt_sd, rec=CASE3_rec_sd),
    list(case_id=4, ind=CASE4_ind_sd, trt=CASE4_trt_sd, rec=CASE4_rec_sd)
  )
  
  avail <- parallelly::availableCores()
  plan(multisession, workers = min(workers, avail, length(cases)))
  
  invisible(future_lapply(
    cases,
    function(cfg) {
      # re-source in worker so simulate_multiple_trials() is available
      source("simulate trial basket.R")
      
      res <- run_main_simulation_basket(
        total_sample_size    = total_sample_size,
        max_sites            = max_sites,
        iterations           = iterations,
        generation_method    = gen_method,
        treatment_proportion = trt_prop,
        individual_sd        = cfg$ind,
        treatment_sd         = cfg$trt,
        recruitment_sd       = cfg$rec,
        design_type          = design_type,
        mu                   = mu,
        min_mortality        = min_mortality,
        max_mortality        = max_mortality,
        alpha                = alpha,
        ni_margin            = ni_margin,
        # basket knobs
        syndrome_num                 = syndrome_num,
        syndrome_prop_var            = syndrome_prop_var,
        syndrome_effect              = syndrome_effect,
        syndrome_effect_variation    = syndrome_effect_variation,
        stratify_by_syndrome         = stratify_by_syndrome,
        save_syndrome_frames         = save_syndrome_frames
      )
      
      out_path <- file.path(folder, sprintf("basket_%s_case%d.rds", file_prefix, cfg$case_id))
      saveRDS(res, file = out_path)
      message(sprintf("Saved: %s", out_path))
      TRUE
    },
    future.seed = seed_base
  ))
  
  plan(sequential)
}

# =========================
# Run scenarios; within each, the 4 cases run in parallel
# =========================

# SUPERIORITY — POWER (mu = +0.10), mortality 0.25–0.75
run_four_cases_parallel_basket(
  folder = "superiority_power", file_prefix = "sup_power",
  total_sample_size = SUP_SS_POWER, design_type = "superiority", mu = SUP_MU_POWER,
  min_mortality = SUP_MIN_MORT, max_mortality = SUP_MAX_MORT
)

# SUPERIORITY — FALSE POSITIVE (mu = 0.00), mortality 0.25–0.75
run_four_cases_parallel_basket(
  folder = "superiority_fp", file_prefix = "sup_fp",
  total_sample_size = SUP_SS_FP, design_type = "superiority", mu = SUP_MU_FP,
  min_mortality = SUP_MIN_MORT, max_mortality = SUP_MAX_MORT
)

# NON-INFERIORITY — POWER (mu = 0.00), mortality 0.20–0.60, N = 840
run_four_cases_parallel_basket(
  folder = "non-inferiority_power", file_prefix = "ni_power",
  total_sample_size = NI_SS_POWER, design_type = "non-inferiority", mu = NI_MU_POWER,
  min_mortality = NI_MIN_MORT, max_mortality = NI_MAX_MORT
)

# NON-INFERIORITY — FALSE POSITIVE (mu = -0.10), mortality 0.20–0.60, N = 824
run_four_cases_parallel_basket(
  folder = "non-inferiority_fp", file_prefix = "ni_fp",
  total_sample_size = NI_SS_FP, design_type = "non-inferiority", mu = NI_MU_FP,
  min_mortality = NI_MIN_MORT, max_mortality = NI_MAX_MORT
)
