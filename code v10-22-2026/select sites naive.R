
select_sites <- function(num_sites, specification = "random_selection") {
  return(paste0("site", seq_len(num_sites)))
}
