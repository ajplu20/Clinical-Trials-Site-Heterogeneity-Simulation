library(tern)
perform_mantel_haenszel_test <- function(df, alpha = 0.05, trial_type = "superiority", ni_margin = 0.1) {
  # Ensure proper types
  df$experiment <- as.factor(df$experiment)
  df$death <- as.factor(df$death)
  df$site <- as.factor(df$site)

  # Run Mantel-Haenszel test
  
  if (trial_type == "superiority") {
    # Drop sites with <= 1 patient
    site_sizes <- table(df$site)
    bad_sites  <- names(site_sizes[site_sizes <= 1])
    if (length(bad_sites) > 0) {
      df <- df[!(df$site %in% bad_sites), , drop = FALSE]
      df <- droplevels(df)  # <<< crucial: remove empty site levels
    }
    
    # If only one site remains, fall back to naive t-test
    if (nlevels(df$site) <= 1) {
      # t-test expects numeric 0/1, not factors
      df$experiment <- as.numeric(as.character(df$experiment))
      df$death      <- as.numeric(as.character(df$death))
      return(perform_naive_t_test(df, alpha = alpha, trial_type = "superiority"))
    }
    
    # Build contingency table and run MH
    table_3d <- xtabs(~ experiment + death + site, data = df)
    storage.mode(table_3d) <- "double"
    mh_result <- mantelhaen.test(table_3d, alternative = "less", correct = FALSE)
    return(as.numeric(mh_result$p.value < alpha))
  }
  else{
    #print("non-inferiority MH")
    #print(df$death)
    rsp <- df$death == 1
    #print(rsp)
    grp <- df$experiment == 1
    #print(grp)
    result <- prop_diff_cmh(rsp = rsp, grp = grp, strata = df$site, conf_level = 1-(alpha*2))
    return(as.numeric(result$diff_ci[2] < ni_margin))
  } 
}

#output <- perform_mantel_haenszel_test(sim_data)
# --- Minimal stress test dataset with 3 sites ---
#df_test <- data.frame(
#  site       = c(rep("Site1", 4), rep("Site2", 1), rep("Site3", 1)),
#  experiment = c(0,0,1,1,   0,  0),   # per site: control vs experiment
#  death      = c(0,1,0,1,   1,   0)    # arbitrary outcomes
#)

#print(df_test)

# Run Mantel–Haenszel superiority test
#res_sup <- perform_mantel_haenszel_test(df_test, alpha = 0.05, trial_type = "superiority")
#cat("Superiority test result:", res_sup, "\n")

# Run Mantel–Haenszel non-inferiority test
#res_ni <- perform_mantel_haenszel_test(df_test, alpha = 0.05, trial_type = "non-inferiority")
#cat("Non-inferiority test result:", res_ni, "\n")
