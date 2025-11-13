# Parallel debugging script for simulate_power_fp
# Requires: future, future.apply
library(future)
library(future.apply)

debug_replicates_parallel <- function(reps = 1000, n_trials = 1000, seed = NULL, workers = 10) {
  # Preserve current plan and restore on exit
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  
  # Windows -> multisession; Unix -> multicore
  if (.Platform$OS.type == "windows") {
    future::plan(multisession, workers = workers)
  } else {
    future::plan(multicore, workers = workers)
  }
  
  if (!is.null(seed)) set.seed(seed)
  
  idx <- seq_len(reps)
  
  res_list <- future.apply::future_lapply(
    X = idx,
    FUN = function(i) {
      # IMPORTANT: each worker needs the function definition
      source("simulate_trial_rshiny.R")
      
      r <- simulate_power_fp(
        sample_size    = 824,
        n_sites        = 3,
        mortality_mode = "per_site",
        mortality_list = c(0.2, 0.4, 0.6),
        mortality_mean = NA,
        mortality_sd   = NA,
        individual_sd  = 0.05,
        recruit_mode   = "balanced",
        recruit_sd     = NA,
        recruit_props  = numeric(0),
        te_mode        = "uniform",
        te_const       = NA,
        te_mean        = 0.1,
        te_spread      = 0.2,
        te_sd          = NA,
        is_ni          = TRUE,
        ni_margin      = 0.1,
        n_trials       = n_trials,
        alpha          = 0.05,
        
      )
      
      c(power = as.numeric(r$power), fpr = as.numeric(r$fpr))  # keep as percentages
    },
    future.seed = TRUE  # parallel-safe RNG (L'Ecuyer-CMRG)
  )
  
  mat <- do.call(rbind, res_list)
  pow <- mat[, "power"]
  fpr <- mat[, "fpr"]
  
  summarize <- function(x) {
    m  <- mean(x, na.rm = TRUE)
    sd <- stats::sd(x, na.rm = TRUE)
    se <- sd / sqrt(length(x))
    ci <- m + c(-1, 1) * 1.96 * se
    data.frame(
      mean_pct     = m,
      sd_pct       = sd,
      se_pct       = se,
      ci95_low_pct = ci[1],
      ci95_high_pct= ci[2],
      reps         = length(x)
    )
  }
  
  list(
    reps_power_pct = pow,
    reps_fpr_pct   = fpr,
    summary_power  = summarize(pow),
    summary_fpr    = summarize(fpr)
  )
}

# ---- Run and save ----
out <- debug_replicates_parallel(reps = 1000, n_trials = 1000, seed = 123, workers = 10)
saveRDS(out, file = "debug_replicates_out_ni.rds")
out <- readRDS("debug_replicates_out_ni.rds")

print(out$summary_power)
print(out$summary_fpr)
# out$reps_power_pct  # vector of the 1000 power estimates (percent)
