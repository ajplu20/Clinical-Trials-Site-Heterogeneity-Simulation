generate_syndrome_effect_df <- function(selected_sites, syndrome_effect, syndrome_effect_variation){
  df <- {out <- sapply(seq_along(syndrome_effect),
                       \(j) rnorm(length(selected_sites),
                                  mean = syndrome_effect[[j]],
                                  sd   = syndrome_effect_variation[[j]]));
  rownames(out) <- selected_sites;
  setNames(as.data.frame(out),
           paste0("syndrome_type_", seq_along(syndrome_effect)))}
  
  return(df)
}

#test
#selected_sites <- c("site_A", "site_B", "site_C")
#syndrome_effect <- c(0.2, -0.2)
#syndrome_effect_variation <- c(0.0, 0.0)

#generate_syndrome_effect_df(selected_sites, syndrome_effect, syndrome_effect_variation)
