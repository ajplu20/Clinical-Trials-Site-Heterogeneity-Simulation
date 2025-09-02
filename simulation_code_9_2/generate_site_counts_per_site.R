generate_site_counts_per_site <- function(site_list, total_sample, recruitment_sd = 0) {
  n_sites <- length(site_list)
  expected_prop <- 1 / n_sites
  
  # Draw noisy proportions
  raw_props <- rnorm(n_sites, mean = expected_prop, sd = recruitment_sd)
  raw_props <- pmax(raw_props, 0.0001)  # Avoid negative proportions
  
  # Normalize so they sum to 1
  normalized_props <- raw_props / sum(raw_props)
  
  # Allocate initial samples (rounding up to avoid 0)
  site_allocations <- ceiling(normalized_props * total_sample)
  
  # Cumulative sum check: go through and fix if overflow
  csum <- 0
  for (i in 1:(n_sites - 1)) {
    csum <- csum + site_allocations[i]
    if (csum >= total_sample) {
      # Set remaining sites to 1, including the last one
      site_allocations[(i+1):(n_sites - 1)] <- 1
      break
    }
  }
  
  # Final site: assign remaining sample, or at least 1
  csum <- sum(site_allocations[1:(n_sites - 1)])
  site_allocations[n_sites] <- max(total_sample - csum, 1)
  
  return(site_allocations)
}




#sites_to_allocate <- c("VN 001", "SG 001","SG 003")
#total_sample_size <- 2000
#allocated_counts <- c()
#allocated_counts <- generate_site_counts_per_site(sites_to_allocate, total_sample_size, 0.5)

#for (i in 1:10000) {
#  allocated_counts <- append(allocated_counts,generate_site_counts_per_site(sites_to_allocate, total_sample_size, 0.04))
#}
  
#print(sd(allocated_counts))
