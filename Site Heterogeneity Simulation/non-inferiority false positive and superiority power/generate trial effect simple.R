#library(truncnorm)

#change this function to test different distributions of mean
generate_trial_effect <- function(num_sites, treatment_sd = 0) {

  mu <- 0.1
  return (runif(num_sites, min = mu - treatment_sd, max = mu + treatment_sd))
  
}
