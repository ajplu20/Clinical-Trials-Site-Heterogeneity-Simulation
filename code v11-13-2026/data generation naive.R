generate_random_dataset_distribution <- function(sites, site_counts) {
  if (length(sites) != length(site_counts)) {
    stop("Length of sites and site_counts must be equal.")
  }
  
  all_rows <- vector("list", length(sites))
  
  for (i in seq_along(sites)) {
    n <- site_counts[i]
    #print("n is")
    #print(n)
    site <- sites[i]
    
    site_df <- data.frame(
      cal_agey = rep(0, n),
      sex = rep("M", n),
      vap = rep(0, n),
      hpd_admreason = rep("GENERIC", n),
      comorbidities_CCI = rep(0, n),
      site = rep(site, n),
      stringsAsFactors = FALSE
    )
    
    all_rows[[i]] <- site_df
  }
  
  return(do.call(rbind, all_rows))
}

# Example: simulate 5 patients from "BD 001", 3 from "JP 001"
sites <- c("BD 001", "JP 001")
site_counts <- c(5, 3)

synthetic_df <- generate_random_dataset_distribution(sites, site_counts)
print(synthetic_df)
