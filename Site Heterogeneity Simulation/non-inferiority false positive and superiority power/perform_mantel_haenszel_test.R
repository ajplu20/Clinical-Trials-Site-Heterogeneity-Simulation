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
  else if (trial_type == "non-inferiority-naive"){
    mh_result <- mantelhaen.test(table_3d, alternative = "less", correct = FALSE, conf.level = 1-alpha) #no correction cause our sample size is pretty large usually
    return(mh_result$conf.int[2] < 1.5)
  }
  else{
    # Initialize vectors
    risk_diffs <- c()
    weights <- c()
    
    # Loop over strata (sites)
    for (s in dimnames(table_3d)$site) {
      tab <- table_3d[, , s]
      
      # Ensure 2x2 matrix format even if levels are missing
      a <- ifelse("1" %in% rownames(tab) && "1" %in% colnames(tab), tab["1", "1"], 0)
      #print(a)
      b <- ifelse("1" %in% rownames(tab) && "0" %in% colnames(tab), tab["1", "0"], 0)
      #print(b)
      c <- ifelse("0" %in% rownames(tab) && "1" %in% colnames(tab), tab["0", "1"], 0)
      #print(c)
      d <- ifelse("0" %in% rownames(tab) && "0" %in% colnames(tab), tab["0", "0"], 0)
      #print(d)
      
      n1 <- a + b  # treatment group size
      n0 <- c + d  # control group size
      
      if (n1 > 0 && n0 > 0) {
        p1 <- a / n1
        p0 <- c / n0
        rd <- p1 - p0
        var <- (p1 * (1 - p1) / n1) + (p0 * (1 - p0) / n0)
        w <- 1 / var
        risk_diffs <- c(risk_diffs, rd)
        weights <- c(weights, w)
      }
    }
    
    # Combine stratum-wise risk differences using inverse-variance weighting
    w_sum <- sum(weights)
    #print(w_sum)
    rd_pooled <- sum(weights * risk_diffs) / w_sum
    se_pooled <- sqrt(1 / w_sum) #cause in inverse variance weighting, var(RD) = 1/w_sum
    
    # Confidence interval
    ci_upper <- rd_pooled + qnorm(alpha,0,1, lower.tail = FALSE) * se_pooled
    #print("risk diff is")
    #print(rd_pooled)
    #print("se diff is")
    #print(se_pooled)
    
    # Output
    #cat("Mantel-Haenszel adjusted risk difference:", round(rd_pooled, 4), "\n")
    #print(ci_upper)    
    return(ci_upper < 0.1)
  } 
}

#output <- perform_mantel_haenszel_test(sim_data)
