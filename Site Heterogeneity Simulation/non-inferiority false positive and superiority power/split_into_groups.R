split_into_groups <- function(vec, group_sizes) {
  shuffled <- sample(vec)  # Randomly shuffle input
  total <- length(vec)
  
  # Ensure group_sizes isn't too large
  if (sum(group_sizes) >= total) {
    group_sizes[length(group_sizes)] <- total - sum(group_sizes[-length(group_sizes)])
  }
  
  # Compute group boundaries
  split_points <- cumsum(group_sizes)
  start_points <- c(1, head(split_points, -1) + 1)
  
  groups <- Map(function(start, end) shuffled[start:end],
                start_points,
                pmin(split_points, total))
  
  # Add remaining elements to the last group
  if (split_points[length(split_points)] < total) {
    groups[[length(groups)]] <- c(groups[[length(groups)]],
                                  shuffled[(split_points[length(split_points)] + 1):total])
  }
  
  return(groups)
}

#vec <- 1:20
#group_sizes <- c(5, 7, 4)  # Only 16 total, but vec has 20
#group_sizes <- c(5, 7, 15)  # Only 27 total, but vec has 20

#split_into_groups(vec, group_sizes)
