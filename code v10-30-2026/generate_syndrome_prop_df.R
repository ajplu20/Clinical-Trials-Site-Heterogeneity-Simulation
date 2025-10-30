generate_syndrome_prop_df <- function(selected_sites, syndrome_effect, syndrome_prop_var) {
  n_rows <- length(selected_sites)
  n_cols <- length(syndrome_effect)
  mu     <- 1 / n_cols
  
  out <- matrix(NA_real_, n_rows, n_cols)
  
  for (i in seq_len(n_rows)) {
    # 1) Draw, 2) clip negatives to 0
    x <- rnorm(n_cols, mean = mu, sd = syndrome_prop_var)
    x <- pmax(x, 0)
    
    s <- sum(x)
    if (s == 0) {
      # Extreme imbalance: make a random 1-hot vector
      j <- sample.int(n_cols, size = 1)
      x <- rep(0, n_cols)
      x[j] <- 1
    } else {
      # Normalize to proportions
      x <- x / s
    }
    out[i, ] <- x
  }
  
  rownames(out) <- selected_sites
  colnames(out) <- paste0("syndrome_type_", seq_len(n_cols))
  as.data.frame(out)
}


#test
#selected_sites        <- c("A","B","C","D")
#syndrome_effect        <- list(0,0,0)   # length only matters
#syndrome_prop_var     <- 0.5

#df <- generate_syndrome_prop_df(selected_sites, syndrome_effect, syndrome_prop_var)
#df
