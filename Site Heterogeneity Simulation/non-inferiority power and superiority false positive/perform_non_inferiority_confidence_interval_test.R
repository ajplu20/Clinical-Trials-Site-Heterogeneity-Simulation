perform_non_inferiority_confidence_interval_test <- function(df, alpha = 0.05) {
  # find the proportion of death in control and experiement group
  n_control <- length(df$death[df$experiment == 0])
  p_death_control <- mean(df$death[df$experiment == 0])
  
  n_treatment <- length(df$death[df$experiment == 1])
  p_death_treatment <- mean(df$death[df$experiment == 1])
  
  p_diff <- p_death_treatment - p_death_control
  
  SE <-  sqrt((p_death_control * (1-p_death_control) / n_control) 
  + (p_death_treatment * (1-p_death_treatment) / n_treatment))
  #print(SE)
  return((p_diff + SE * qnorm(alpha,0,1, lower.tail = FALSE)) > 0.1)
}
