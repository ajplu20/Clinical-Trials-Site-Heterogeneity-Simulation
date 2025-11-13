library(lme4)
library(marginaleffects)
options(marginaleffects_safe = FALSE)

perform_mixed_effect_model_test <- function(df, alpha = 0.05, trial_type = "superiority", basket = FALSE, ni_margin = 0.1) {
  # Ensure variables are of correct type
  df$death <- as.numeric(df$death)
  df$experiment <- as.numeric(df$experiment)
  df$site <- as.factor(df$site)
  
  if (trial_type == "superiority" && !basket) {
    #print("mixed effect simple superiority")
    
    # Fit fixed-effects logistic regression
    model <- glmer(death ~ experiment + (1 | site), data = df, family = binomial())
    
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
  else if (trial_type == "superiority" && basket) {
    #print("mixed effect basket superiority")
    
    # Fit fixed-effects logistic regression
    if (length(unique(df$syndrome)) > 1) {
      model <- glmer(death ~ experiment + (1 | site) + syndrome,
                     data = df, family = binomial())
    } else {
      model <- glmer(death ~ experiment + (1 | site),
                     data = df, family = binomial())
    }
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
      #print("mixed effect simple non-inferiority")
    
      model <- glmer(death ~ experiment + (1 | site), data = df, family = binomial())
      coef <- avg_comparisons(model, variables = "experiment", type = "response")
      treat_coef <- coef$estimate
      treat_SE  <- coef$std.error
      
      z <- qnorm(1 - alpha)
      upper <- treat_coef + z * treat_SE
      return(as.numeric(upper < ni_margin)) # non-inferiority test
  }
  else {
    #print("mixed effect basket non-inferiority")
    
    if (length(unique(df$syndrome)) > 1) {
      model <- glmer(death ~ experiment + (1 | site) + syndrome,
                     data = df, family = binomial())
    } else {
      model <- glmer(death ~ experiment + (1 | site),
                     data = df, family = binomial())
    }
    coef <- avg_comparisons(model, variables = "experiment", type = "response")
    treat_coef <- coef$estimate
    treat_SE  <- coef$std.error
    
    z <- qnorm(1 - alpha)
    upper <- treat_coef + z * treat_SE
    return(as.numeric(upper < ni_margin)) # non-inferiority test
  }
}

#outcome <- perform_mixed_effect_model_test(sim_data)
