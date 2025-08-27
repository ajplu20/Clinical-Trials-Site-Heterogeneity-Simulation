library(marginaleffects)

perform_fix_effect_model_test <- function(df, alpha = 0.05, trial_type = "superiority", basket = FALSE) {
  # Ensure variables are of correct type
  df$death <- as.numeric(df$death)
  df$experiment <- as.numeric(df$experiment)
  df$site <- as.factor(df$site)
  
  if (trial_type == "superiority" && !basket) {
    #print("fixed effect simple superiority")
    model <- glm(death ~ experiment + site, data = df, family = binomial())
    
    # Extract treatment coefficient and p-value
    coef_summary <- summary(model)$coefficients
    treat_coef <- coef_summary["experiment", "Estimate"]
    treat_pval <- coef_summary["experiment", "Pr(>|z|)"]
    
    # Check if effect is negative and p-value < alpha
    result <- (treat_coef < 0) & ((treat_pval/2) < alpha)
    #print(result)
    #return(list(
    #  passed = result,
    #  coef = treat_coef,
    #  p_value = treat_pval,
    #  model_summary = coef_summary
    #))
    return (as.numeric(result))
  }
  else if (trial_type == "superiority" && basket){
    #print("fixed effect basket superiority")
    model <- glm(death ~ experiment + site + syndrome, data = df, family = binomial())
    
    # Extract treatment coefficient and p-value
    coef_summary <- summary(model)$coefficients
    treat_coef <- coef_summary["experiment", "Estimate"]
    treat_pval <- coef_summary["experiment", "Pr(>|z|)"]
    
    # Check if effect is negative and p-value < alpha
    result <- (treat_coef < 0) & ((treat_pval/2) < alpha)
    #print(result)
    #return(list(
    #  passed = result,
    #  coef = treat_coef,
    #  p_value = treat_pval,
    #  model_summary = coef_summary
    #))
    return (as.numeric(result))
  }
  else if (trial_type != "superiority" && !basket){
    #print("fixed effect simple non-inferiority")
    model <- glm(death ~ experiment + site, data = df, family = binomial())
    coef <- avg_comparisons(model, variables = "experiment", type = "response")
    treat_coef <- coef$estimate
    treat_SE  <- coef$std.error
    
    z <- qnorm(1 - alpha)
    upper <- treat_coef + z * treat_SE
    return(as.numeric(upper < 0.1)) # non-inferiority test
    
  }
  else {
    #print("fixed effect basket non-inferiority")
    model <- glm(death ~ experiment + site + syndrome, data = df, family = binomial())
    coef <- avg_comparisons(model, variables = "experiment", type = "response")
    treat_coef <- coef$estimate
    treat_SE  <- coef$std.error
    z <- qnorm(1 - alpha)
    upper <- treat_coef + z * treat_SE
    return(as.numeric(upper < 0.1)) # non-inferiority test
  } 
}
    
    
    
    


