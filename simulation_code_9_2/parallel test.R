library(future)
library(future.apply)

# Plan for parallel execution
plan(multisession, workers = 3)

# Function to test RNG streams
test_rng <- function(seed = 123) {
  res <- future_lapply(1:3, function(i) {
    rnorm(5, mean = 0, sd = 1)   # draw 5 N(0,1)
  }, future.seed = seed)
  
  return(res)
}

# Run the test
out <- test_rng(123)
print(out)
