generate_syndrome_prop_df <- function(selected_sites, syndrome_effect, syndrome_prop_var){
  df <- {                       # --- build the data-frame inline
    n_rows  <- length(selected_sites)
    n_cols  <- length(syndrome_effect)          # same #cols as before
    mu      <- 1 / n_cols                      # target mean per draw
    
    out <- matrix(NA_real_, n_rows, n_cols)    # storage
    
    for (i in seq_len(n_rows)) {               # ----- row loop
      remain <- 1                              # mass still to allocate
      
      for (j in seq_len(n_cols)) {             # --- column loop
        
        if (j == n_cols) {                     # last column → close row
          out[i, j] <- remain
          break
        }
        
        x <- rnorm(1, mu, syndrome_prop_var)   # candidate draw
        x <- max(0, x)                         # no negatives
        
        if (x > remain) {                      # would overshoot → split
          fill <- remain / (n_cols - j + 1)
          out[i, j:n_cols] <- fill
          break
        }
        
        out[i, j] <- x                         # accept draw
        remain    <- remain - x
      }
    }
    
    rownames(out) <- selected_sites
    setNames(as.data.frame(out),
             paste0("syndrome_type_", seq_len(n_cols)))
  }
  
  return(df)
}

#test
#selected_sites        <- c("A","B","C","D")
#syndrome_effect        <- list(0,0,0)   # length only matters
#syndrome_prop_var     <- 0.5

#df <- generate_syndrome_prop_df(selected_sites, syndrome_effect, syndrome_prop_var)
#df
