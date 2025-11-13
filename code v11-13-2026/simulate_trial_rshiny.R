source("perform_naive_t_test.R")
# ============================================================
# Helpers
# ============================================================

# Truncated normal via simple rejection sampling
rtruncnorm01 <- function(n, mean, sd, a = 0, b = 1, max_iter = 1e6) {
  out <- numeric(n)
  k <- 0L
  i <- 1L
  while (i <= n) {
    k <- k + 1L
    if (k > max_iter) stop("rtruncnorm01: exceeded max_iter; try smaller sd or adjust bounds.")
    x <- rnorm(1L, mean, sd)
    if (x >= a && x <= b) {
      out[i] <- x
      i <- i + 1L
    }
  }
  out
}

# Truncated to [-1, 1] for treatment effects
rtruncnorm11 <- function(n, mean, sd, a = -1, b = 1, max_iter = 1e6) {
  out <- numeric(n)
  k <- 0L
  i <- 1L
  while (i <= n) {
    k <- k + 1L
    if (k > max_iter) stop("rtruncnorm11: exceeded max_iter; try smaller sd or adjust bounds.")
    x <- rnorm(1L, mean, sd)
    if (x >= a && x <= b) {
      out[i] <- x
      i <- i + 1L
    }
  }
  out
}

# Randomized rounding to non-negative integers that sum to N
randomized_round_to_n <- function(props, N) {
  props <- pmax(props, 0)
  if (sum(props) == 0) {
    # uniform fallback
    props <- rep(1/length(props), length(props))
  } else {
    props <- props / sum(props)
  }
  raw <- props * N
  base <- floor(raw)
  short <- N - sum(base)
  # distribute the remaining "short" by highest fractional parts
  frac <- raw - base
  if (short > 0) {
    idx <- order(frac, decreasing = TRUE)[seq_len(short)]
    base[idx] <- base[idx] + 1L
  } else if (short < 0) {
    idx <- order(frac, decreasing = FALSE)[seq_len(abs(short))]
    base[idx] <- pmax(0L, base[idx] - 1L)
  }
  base
}

# Randomized 50/50 split of total n_obs: with 50% chance, ceil for treatment; else floor
randomize_50_50 <- function(n_obs) {
  n_treat <- if (runif(1) < 0.5) ceiling(n_obs / 2) else floor(n_obs / 2)
  n_ctrl  <- n_obs - n_treat
  c(n_treat = n_treat, n_ctrl = n_ctrl)
}

# ============================================================
# Core simulation used inside each trial replicate
# ============================================================

simulate_one_trial_df <- function(
    sample_size,
    n_sites,
    mortality_mode = c("per_site", "mean_sd"),
    mortality_list = numeric(0),
    mortality_mean = 0.1,
    mortality_sd = 0.03,
    individual_sd = 0.05,
    recruit_mode = c("balanced", "imbalance_sd", "custom_props"),
    recruit_sd = 0.05,
    recruit_props = numeric(0),
    te_mode = c("constant", "uniform", "normal"),
    te_const = -0.02,
    te_mean = 0.02,
    te_spread = 0.10,
    te_sd = 0.01,
    override_te_all_zero = FALSE
) {
  mortality_mode <- match.arg(mortality_mode)
  recruit_mode   <- match.arg(recruit_mode)
  te_mode        <- match.arg(te_mode)
  
  # 1) Site mortality vector (length n_sites)
  if (mortality_mode == "per_site") {
    if (length(mortality_list) != n_sites) stop("mortality_list length must equal n_sites.")
    site_mort <- pmin(1, pmax(0, as.numeric(mortality_list)))
  } else {
    # mean/sd -> truncated normal in [0,1]
    site_mort <- rtruncnorm01(n_sites, mean = mortality_mean, sd = mortality_sd, a = 0, b = 1)
  }
  
  # 2) Patients per site (length n_sites)
  if (recruit_mode == "balanced") {
    # start with equal split; randomized rounding of any remainder
    props <- rep(1 / n_sites, n_sites)
    n_per_site <- randomized_round_to_n(props, sample_size)
  } else if (recruit_mode == "imbalance_sd") {
    # Perturb around equal share with Normal(mean=1/n_sites, sd=recruit_sd), clamp to >=0, renormalize
    base <- rnorm(n_sites, mean = 1 / n_sites, sd = recruit_sd)
    base <- pmax(base, 0)
    if (sum(base) == 0) base <- rep(1 / n_sites, n_sites)
    base <- base / sum(base)
    n_per_site <- randomized_round_to_n(base, sample_size)
  } else {
    # custom proportions provided by user (assumed validated)
    if (length(recruit_props) != n_sites) stop("recruit_props length must equal n_sites.")
    base <- as.numeric(recruit_props)
    base <- pmax(base, 0)
    if (sum(base) == 0) base <- rep(1 / n_sites, n_sites)
    base <- base / sum(base)
    n_per_site <- randomized_round_to_n(base, sample_size)
  }
  
  # 3) Treatment effect per site (length n_sites), then apply overrides for FPR/Power scenarios
  if (override_te_all_zero) {
    te_const <- 0
    te_mean <- 0
  } 
  if (te_mode == "constant") {
    site_te <- rep(te_const, n_sites)
  } else if (te_mode == "uniform") {
    lo <- te_mean - te_spread / 2
    hi <- te_mean + te_spread / 2
    # allow TE in [-1, 1]
    lo <- max(-1, lo)
    hi <- min( 1, hi)
    if (lo > hi) lo <- hi
    site_te <- runif(n_sites, min = lo, max = hi)
  } else {
    # normal truncated to [-1, 1]
    site_te <- rtruncnorm11(n_sites, mean = te_mean, sd = te_sd, a = -1, b = 1)
  }
  
  
  # ---- Build patient-level dataframe ----
  site_id <- rep(seq_len(n_sites), times = n_per_site)
  n_obs   <- length(site_id)
  
  # Randomize treatment 1/2 vs 1/2 (global), with randomized ceil/floor split
  split <- randomize_50_50(n_obs)
  treat_idx <- sample.int(n_obs, size = split["n_treat"], replace = FALSE)
  experiment <- integer(n_obs)
  experiment[treat_idx] <- 1L
  
  # Base patient mortality ~ N(site_mort[site], individual_sd), clamped to [0,1]
  base_mu <- site_mort[site_id]
  mu <- pmin(1, pmax(0, rnorm(n_obs, mean = base_mu, sd = individual_sd)))
  
  # Apply treatment effect (subtract if assigned to treatment), then clamp to [0,1]
  te_by_patient <- site_te[site_id]
  mu_post <- mu
  mu_post[experiment == 1L] <- mu_post[experiment == 1L] - te_by_patient[experiment == 1L]
  mu_post <- pmin(1, pmax(0, mu_post))
  
  # Realize deaths (Bernoulli)
  death <- rbinom(n_obs, size = 1L, prob = mu_post)
  
  data.frame(
    site       = factor(site_id),
    mortality  = mu_post,
    experiment = experiment,
    death      = death
  )
}

# ============================================================
# Main entry: compute Power and FPR
# ============================================================

simulate_power_fp <- function(
    sample_size,
    n_sites,
    mortality_mode = c("per_site", "mean_sd"),
    mortality_list = numeric(0),
    mortality_mean = 0.10,
    mortality_sd = 0.03,
    individual_sd = 0.05,
    recruit_mode = c("balanced", "imbalance_sd", "custom_props"),
    recruit_sd = 0.05,
    recruit_props = numeric(0),
    te_mode = c("constant", "uniform", "normal"),
    te_const = -0.02,
    te_mean = 0.02,
    te_spread = 0.10,
    te_sd = 0.01,
    is_ni = FALSE,
    ni_margin = 0.10,
    n_trials = 1000,
    alpha = 0.05,
    seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)
  
  mortality_mode <- match.arg(mortality_mode)
  recruit_mode   <- match.arg(recruit_mode)
  te_mode        <- match.arg(te_mode)
  
  # Wrapper to run one replicate given override for TE=0 or not
  run_repl <- function(override_zero) {
    df <- simulate_one_trial_df(
      sample_size    = sample_size,
      n_sites        = n_sites,
      mortality_mode = mortality_mode,
      mortality_list = mortality_list,
      mortality_mean = mortality_mean,
      mortality_sd   = mortality_sd,
      individual_sd  = individual_sd,
      recruit_mode   = recruit_mode,
      recruit_sd     = recruit_sd,
      recruit_props  = recruit_props,
      te_mode        = te_mode,
      te_const       = te_const,
      te_mean        = te_mean,
      te_spread      = te_spread,
      te_sd          = te_sd,
      override_te_all_zero = override_zero
    )
    # Use the user's t-test function (assumed available in env)
    if (is_ni) {
      perform_naive_t_test(df, alpha = alpha, trial_type = "non-inferiority", ni_margin = ni_margin)
    } else {
      perform_naive_t_test(df, alpha = alpha, trial_type = "superiority")
    }
  }
  
  # Logic per your spec:
  # Superiority: Power = user TE; FPR = TE=0
  # Non-inferiority: Power = TE=0;  FPR = user TE (i.e., inferior scenario)
  if (!is_ni) {
    # Superiority
    rejections_power <- replicate(n_trials, run_repl(override_zero = FALSE))
    rejections_fpr   <- replicate(n_trials, run_repl(override_zero = TRUE))
  } else {
    # Non-inferiority
    rejections_power <- replicate(n_trials, run_repl(override_zero = TRUE))   # TE = 0
    rejections_fpr   <- replicate(n_trials, run_repl(override_zero = FALSE))  # user TE (possibly < 0 -> inferior)
  }
  
  power <- mean(rejections_power) * 100
  fpr   <- mean(rejections_fpr)   * 100
  
  list(power = power, fpr = fpr)
}
