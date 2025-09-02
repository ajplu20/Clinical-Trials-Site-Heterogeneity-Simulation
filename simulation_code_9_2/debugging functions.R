
library(dplyr)
library(readr)

df <- read_csv("final_cleaned_data.csv")
valid_sites <- c(
  "BD 001", "BN 001", "CN 001", "CN 002", "HK 000", "ID 004", "ID 005", "ID 006",
  "IN 004", "IN 007", "JP 001", "KH 002", "MG 001", "MG 000", "MY 001", "MY 002",
  "MY 003", "MY 004", "MY 005", "NP 001", "PH 001", "PH 002", "PK 001", "PK 002",
  "SG 002", "PK 003", "SG 001", "SL 001", "TH 001", "TH 002", "TH 004", "TH 005",
  "TL 001", "TW 001", "VN 001", "VN 004", "VN 005"
)


#######test generated distribution is correct (ie, matchs acornhai distribution)

#find mean mortality in the simulate data for each site
compute_site_mortality_means <- function(results_df) {
  # Initialize an empty list for each site
  site_mortality_lists <- setNames(vector("list", length(valid_sites)), valid_sites)
  counter <- 1
  # Loop through each simulation data frame
  for (sim in results_df$simulation_data) {
    print(counter)
    counter <- counter + 1
    # Loop through each trial data frame inside the simulation
    for (trial_df in sim$trial_data) {
      # Subset to control group (experiment == 0)
      control_df <- subset(trial_df, experiment == 0)
      
      # Loop through each row of the control group
      for (i in seq_len(nrow(control_df))) {
        site <- control_df$site[i]
        prob <- control_df$mortality_probability[i]
        
        # Append probability to corresponding site's list
        if (site %in% valid_sites) {
          site_mortality_lists[[site]] <- c(site_mortality_lists[[site]], prob)
        }
      }
    }
  }
  
  # Compute mean mortality probability for each site
  site_means <- sapply(site_mortality_lists, function(probs) {
    if (length(probs) > 0) mean(probs) else NA
  })
  
  return(site_means)
}

# Example usage:
#valid_sites <- paste("MY", formatC(1:37, width = 3, flag = "0"))
site_mortality_means <- compute_site_mortality_means(results_df)
print(site_mortality_means)
saveRDS(site_mortality_means, file = "site_mortality_means_fp.rds")

site_mortality_means <- readRDS("site_mortality_means_pw.rds")
#find the sites with mean moretality difference larger than a certain number
compare_site_mortality <- function(df, site_mortality_means, max_diff = 0.05) {
  # Step 1: Compute observed mean mortality from df
  observed_means <- df %>%
    group_by(site) %>%
    summarise(observed_mean = mean(mortality, na.rm = TRUE)) %>%
    ungroup()
  
  # Step 2: Convert site_mortality_means into a dataframe
  sim_means <- data.frame(site = names(site_mortality_means),
                          simulated_mean = as.numeric(site_mortality_means),
                          stringsAsFactors = FALSE)
  
  # Step 3: Join observed and simulated means
  comparison <- left_join(observed_means, sim_means, by = "site") %>%
    mutate(diff = abs(observed_mean - simulated_mean)) %>%
    filter(!is.na(diff) & diff > max_diff)
  
  return(comparison)
}

# Example usage:
outlier_sites_df <- compare_site_mortality(df, site_mortality_means, max_diff = 0.01)



#plot acornhai mortality per site vs observed mortality per site:
library(dplyr)
library(ggplot2)

graph_site_mortality <- function(df, site_mortality_means) {
  # Step 1: Compute observed mean mortality from df
  observed_means <- df %>%
    group_by(site) %>%
    summarise(observed_mean = mean(mortality, na.rm = TRUE)) %>%
    ungroup()
  
  # Step 2: Convert site_mortality_means into a dataframe
  sim_means <- data.frame(site = names(site_mortality_means),
                          simulated_mean = as.numeric(site_mortality_means),
                          stringsAsFactors = FALSE)
  
  # Step 3: Merge both
  comparison <- left_join(observed_means, sim_means, by = "site") %>%
    filter(!is.na(simulated_mean))
  
  # Step 4: Plot observed vs simulated
  ggplot(comparison, aes(x = simulated_mean, y = observed_mean)) +
    geom_point() +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    labs(
      title = "Observed vs Simulated Site-Level Mortality",
      x = "Simulated Mortality (results_df)",
      y = "Observed Mortality (From ACORNHAI dataset)"
    ) +
    theme_minimal()
}


graph_site_mortality(df, site_mortality_means)








#########plotting differnece in mortality between experiments and control group for number of sites:
library(dplyr)
library(ggplot2)

plot_mortality_diff_by_sites <- function(results_df) {
  # Initialize list to store per-trial mortality differences
  diff_list <- list()
  
  # Loop through each simulation
  for (i in seq_along(results_df$simulation_data)) {
    sim <- results_df$simulation_data[[i]]
    num_sites <- results_df$num_sites[i]
    
    # Loop through each trial_data dataframe
    for (trial_df in sim$trial_data) {
      # Compute mean deaths in control and treatment arms
      mean_control <- mean(trial_df$death[trial_df$experiment == 0], na.rm = TRUE)
      mean_treatment <- mean(trial_df$death[trial_df$experiment == 1], na.rm = TRUE)
      
      # Append result
      diff_list[[length(diff_list) + 1]] <- data.frame(
        num_sites = num_sites,
        mortality_diff = mean_treatment - mean_control
      )
    }
  }
  
  # Combine into one dataframe
  plot_data <- bind_rows(diff_list)
  plot_data %>% group_by(num_sites) %>% summarise(mean_diff = mean(mortality_diff, na.rm = TRUE))
  
  # Generate the boxplot
  ggplot(plot_data, aes(x = factor(num_sites), y = mortality_diff)) +
    geom_boxplot() +
    labs(
      title = "Mortality Difference (Treatment - Control) by Number of Sites",
      x = "Number of Sites",
      y = "Mortality Rate Difference"
    )
}



plot_mortality_diff_by_sites(results_df)





###########make confidence interval of difference per run per site number
library(dplyr)
library(lme4)


# -------------- helper: Welch-t CI width (difference in risk) -----------------
t_width <- function(x_treat, x_ctrl, ci_level = 0.95) {
  if (length(x_treat) == 0 || length(x_ctrl) == 0) return(NA_real_)
  tt <- t.test(x_treat, x_ctrl,
               var.equal   = FALSE,              # Welch
               alternative = "two.sided")       # two-sided for finite width
  diff(range(tt$conf.int))
}

# --------------------------- main aggregator ----------------------------------
compute_avg_ci_width_by_num_sites <- function(results_df, ci_level = 0.95) {
  alpha   <- 1 - ci_level
  z_crit  <- qnorm(1 - alpha/2)                 # 1.96 for 95%
  
  run_ci_widths <- vector("list", length(results_df$simulation_data))
  
  for (i in seq_along(results_df$simulation_data)) {
    sim       <- results_df$simulation_data[[i]]
    num_sites <- results_df$num_sites[i]
    print(num_sites)
    t_w     <- c()
    fix_w   <- c()
    mixed_w <- c()
    mh_w    <- c()                              # NEW: Mantel-Haenszel
    
    for (trial_df in sim$trial_data) {
      ## ---------------- t-test (risk difference) ----------------
      death_ctrl  <- trial_df$death[trial_df$experiment == 0]
      death_treat <- trial_df$death[trial_df$experiment == 1]
      t_w         <- c(t_w, t_width(death_treat, death_ctrl, ci_level))
      
      ## ------------- fixed-effect logistic ----------------------
      trial_df$site <- factor(trial_df$site)
      
      # ---- Fixed effects model (logistic regression with site dummies) ----
      fmod <- glm(death ~ experiment + site, data = trial_df, family = binomial())
      
      if (!is.na(coef(fmod)["experiment"])) {
        coef_fix <- coef(fmod)["experiment"]
        se_fix   <- summary(fmod)$coeff["experiment", "Std. Error"]
        
        ci_low_fix <- coef_fix - z_crit * se_fix
        ci_high_fix <- coef_fix + z_crit * se_fix
        
        fix_w <- c(fix_w, exp(ci_high_fix) - exp(ci_low_fix))  # width on OR scale
      } else {
        fix_w <- c(fix_w, NA_real_)
      }
      
      # ---- Mixed-effects model (GLMM with random intercepts by site) ----
      mmod <- tryCatch(
        glmer(death ~ experiment + (1 | site),
              data = trial_df, family = binomial(),
              control = glmerControl(optimizer = "bobyqa", calc.derivs = FALSE)),
        error = function(e) NULL
      )
      
      if (!is.null(mmod) && !is.na(fixef(mmod)["experiment"])) {
        coef_mix <- fixef(mmod)["experiment"]
        se_mix   <- summary(mmod)$coeff["experiment", "Std. Error"]
        
        ci_low_mix <- coef_mix - z_crit * se_mix
        ci_high_mix <- coef_mix + z_crit * se_mix
        
        mixed_w <- c(mixed_w, exp(ci_high_mix) - exp(ci_low_mix))  # width on OR scale
      } else {
        mixed_w <- c(mixed_w, NA_real_)
      }
      
      ## ------------- Mantel-Haenszel odds ratio -----------------
      mh_try <- tryCatch({
        tab3d   <- xtabs(~ experiment + death + site, data = trial_df)
        mh      <- mantelhaen.test(tab3d,
                                   alternative = "two.sided",
                                   correct = FALSE)
        diff(range(mh$conf.int))                  # width on OR scale
      }, error = function(e) NA_real_)
      
      mh_w <- c(mh_w, mh_try)
    } # end trial loop
    
    run_ci_widths[[i]] <- data.frame(
      num_sites       = num_sites,
      t_test_width    = mean(t_w,     na.rm = TRUE),
      fixed_width     = mean(fix_w,   na.rm = TRUE),
      mixed_width     = mean(mixed_w, na.rm = TRUE),
      mh_width        = mean(mh_w,    na.rm = TRUE)   # NEW
    )
  } # end sim loop
  
  bind_rows(run_ci_widths) %>%
    group_by(num_sites) %>%
    summarise(across(
      c(t_test_width, fixed_width, mixed_width, mh_width),
      ~ mean(.x, na.rm = TRUE),
      .names = "avg_{.col}"
    ), .groups = "drop")
}


conf_int <- compute_avg_ci_width_by_num_sites(results_df)
library(tidyr)
library(ggplot2)

# Assuming your dataframe is named conf_int

# Pivot all y-variables into long format
conf_int_long <- conf_int %>%
  pivot_longer(
    cols = -num_sites,            # all columns except num_sites
    names_to = "method",          # name for the new column holding old column names
    values_to = "ci_width"        # name for the values
  )

# Plot all lines
ggplot(conf_int_long, aes(x = num_sites, y = ci_width, color = method)) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  ylim(0, 0.5) +  
  labs(
    title = "95% Confidence Interval Widths by Method and Number of Sites",
    x = "Number of Sites",
    y = "Average CI Width",
    color = "Method"
  ) +
  theme_minimal()


##check site is indeed randomly assigned and different for each one
table(results_df$simulation_data[[5]]$trial_data[[5]]$site)
table(results_df$simulation_data[[5]]$trial_data[[50]]$site)
table(results_df$simulation_data[[3]]$trial_data[[5]]$site)
table(results_df$simulation_data[[3]]$trial_data[[50]]$site)


##check equal division between sites
table(results_df$simulation_data[[5]]$trial_data[[6]]$experiment[
  results_df$simulation_data[[5]]$trial_data[[6]]$site == "MY 004"
])



#show mean mortalty per site of acornhai data
mortality_frame <- 
  df %>%
  group_by(site) %>%
  summarise(mean_mortality = mean(mortality, na.rm = TRUE)) %>%
  arrange(site)

mortality_frame






#######test equal assignment to sites
results_df <- readRDS("results_df_pw.rds")

i <- 24
j <- 500
table(results_df$simulation_data[[i]]$trial_data[[j]]$experiment,results_df$simulation_data[[i]]$trial_data[[j]]$site)
