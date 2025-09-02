library(tern)
perform_mantel_haenszel_test <- function(df, alpha = 0.05, trial_type = "superiority") {
  # Ensure proper types
  df$experiment <- as.factor(df$experiment)
  df$death <- as.factor(df$death)
  df$site <- as.factor(df$site)
  
  # Create a 3D contingency table: treatment × outcome × site
  table_3d <- xtabs(~ experiment + death + site, data = df)
  storage.mode(table_3d) <- "double"
  
  #print(table_3d)
  # Run Mantel-Haenszel test
  
  if (trial_type == "superiority") {
    #print("superiority MH")
    mh_result <- mantelhaen.test(table_3d, alternative = "less", correct = FALSE) #no correction cause our sample size is pretty large usually
    
    # Check if we reject null hypothesis
    reject_null <- mh_result$p.value < alpha
    
    #return(list(
    #  passed = reject_null,
    #  p_value = mh_result$p.value,
    #  common_odds_ratio = mh_result$estimate,
    #  test_statistic = mh_result$statistic,
    #  alpha = alpha,
    #  details = mh_result
    #))
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
    return(result$diff_ci[2] < 0.1)
  } 
}

#output <- perform_mantel_haenszel_test(sim_data)
