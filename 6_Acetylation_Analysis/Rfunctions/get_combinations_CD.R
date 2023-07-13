get_combinations <- function(my_conditions){
  combinations <- c()
  n_comb <- length(my_conditions)
  for (a in c(1:(n_comb-1))){
    for (b in c((a+1):n_comb)){
      combinations = c(combinations, paste(my_conditions[a], "-", my_conditions[b]))
    }
  }
  return(combinations)
}
