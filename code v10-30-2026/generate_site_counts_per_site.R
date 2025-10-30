generate_site_counts_per_site <- function(site_list, total_sample, recruitment_sd = 0) {
  n_sites <- length(site_list)
  expected_prop <- 1 / n_sites
  
  # Draw noisy proportions
  raw_props <- rnorm(n_sites, mean = expected_prop, sd = recruitment_sd)
  raw_props <- pmax(raw_props, 0.0001)  # Avoid negative proportions
  
  # Normalize so they sum to 1
  normalized_props <- raw_props / sum(raw_props)
  
  # Suppose normalized_props * total_sample gives fractional allocations
  raw_allocations <- normalized_props * total_sample
  
  # Random Bernoulli(0.5) for each element
  rand_flip <- rbinom(length(raw_allocations), size = 1, prob = 0.5)
  
  # Apply ceiling if 1, floor if 0
  site_allocations <- ifelse(rand_flip == 1,
                             ceiling(raw_allocations),
                             floor(raw_allocations))
  
  
  return(site_allocations)
}




#sites_to_allocate <- c("VN 001", "SG 001","SG 003")
#total_sample_size <- 2000
#allocated_counts <- c()
#allocated_counts <- generate_site_counts_per_site(sites_to_allocate, total_sample_size, 0.5)
#allocated_counts
#sum(allocated_counts)
#for (i in 1:10000) {
#  allocated_counts <- append(allocated_counts,generate_site_counts_per_site(sites_to_allocate, total_sample_size, 0.04))
##}
  
#print(sd(allocated_counts))
