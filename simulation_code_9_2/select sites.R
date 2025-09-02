valid_sites <- c(
  "BD 001", "BN 001", "CN 001", "CN 002", "HK 000", "ID 004", "ID 005", "ID 006",
  "IN 004", "IN 007", "JP 001", "KH 002", "MG 001", "MG 000", "MY 001", "MY 002",
  "MY 003", "MY 004", "MY 005", "NP 001", "PH 001", "PH 002", "PK 001", "PK 002",
  "SG 002", "PK 003", "SG 001", "SL 001", "TH 001", "TH 002", "TH 004", "TH 005",
  "TL 001", "TW 001", "VN 001", "VN 004", "VN 005"
)



select_sites <- function(num_sites, specification = "random_selection") {
  
  if (specification == "random_selection") {
    #print("inside")
    if (num_sites > length(valid_sites)) {
      #print("error")
      stop("Error: Number of sites requested exceeds number of available sites.")
    }
    return(sample(valid_sites, num_sites, replace = FALSE))
  }
  
  if (specification == "different_country") {
    # Extract country codes (first 2 characters)
    country_codes <- substr(valid_sites, 1, 2)
    country_site_map <- split(valid_sites, country_codes)
    #print(length(country_site_map))
    if (num_sites > length(country_site_map)) {
      stop("Error: Number of sites requested exceeds number of distinct countries.")
    }
    
    selected_countries <- sample(names(country_site_map), num_sites)
    selected_sites <- sapply(selected_countries, function(cty) {
      sample(country_site_map[[cty]], 1)
    })
    
    return(unname(selected_sites))
  }
  
  stop("Error: Unknown specification. Use 'random_selection' or 'different_country'.")
}

