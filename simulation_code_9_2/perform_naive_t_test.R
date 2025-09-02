perform_naive_t_test <- function(df, alpha = 0.05, trial_type = "superiority") {
  # Make sure experiment and death are numeric
  df$experiment <- as.numeric(df$experiment)
  df$death <- as.numeric(df$death)
  
  # Subset deaths by group
  death_control <- df$death[df$experiment == 0]
  death_treatment <- df$death[df$experiment == 1]
  
  if (trial_type == "superiority") {
    #print("superiority T test")
    # Perform Welch's t-test (unequal variances), one-sided: treatment < control
    t_result <- t.test(death_treatment, death_control,
                       alternative = "less",  # treatment group < control
                       var.equal = FALSE)     # Welch's t-test
    
    return(as.numeric(t_result$p.value < alpha))    
  }
  
  else {
    #print("non-inferiority T test")
    # Perform Welch's t-test (unequal variances), one-sided: treatment < control
    t_result <- t.test(death_treatment,death_control,
                       alternative = "less",  # treatment group < control
                       var.equal = FALSE,
                       conf.level = 1-alpha)     # Welch's t-test
    
    return(as.numeric(t_result$conf.int[2] < 0.1))    
      
  }

}

