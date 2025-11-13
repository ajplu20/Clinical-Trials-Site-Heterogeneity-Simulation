# Debugging script for simulate_power_fp

source("simulate_trial_rshiny.R")

debug_replicates <- function(reps = 100, n_trials = 1000, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  pow <- numeric(reps)
  fpr <- numeric(reps)
  
  for (i in seq_len(reps)) {
    print(i)
    res <- simulate_power_fp(
      sample_size    = 840,
      n_sites        = 3,
      mortality_mode = "per_site",
      mortality_list = c(0.25, 0.5, 0.75),
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
      is_ni          = FALSE,
      ni_margin      = NA,
      n_trials       = n_trials,
      alpha          = 0.05
    )
    
    pow[i] <- as.numeric(res$power)
    fpr[i] <- as.numeric(res$fpr)
  }
  
  summarize <- function(x) {
    m  <- mean(x, na.rm = TRUE)
    sd <- stats::sd(x, na.rm = TRUE)
    se <- sd / sqrt(length(x))
    ci <- m + c(-1, 1) * 1.96 * se
    data.frame(
      mean_pct = m,
      sd_pct = sd,
      se_pct = se,
      ci95_low_pct = ci[1],
      ci95_high_pct = ci[2],
      reps = length(x)
    )
  }
  
  list(
    reps_power_pct = pow,
    reps_fpr_pct   = fpr,
    summary_power  = summarize(pow),
    summary_fpr    = summarize(fpr)
  )
}

# Run the debug experiment
out <- debug_replicates(reps = 1000, n_trials = 1000, seed = 123)
saveRDS(out, file = "debug_replicates_out.rds")

print(out$summary_power)
print(out$summary_fpr)
out$reps_power_pct

#acceptable se in power estimation variation in power
sd(out$reps_fpr_pct)#is the current interval

#acceptable se in power estimate. 

