
rho_test_perm = function(x, y, R = 1000, verbose = FALSE, return_dist = TRUE, use_rank = TRUE){
  
  if(use_rank) rho_f = rho_rank
  if(!use_rank) rho_f = rho
  
  rho_minus_ones = NA
  
  if(return_dist){
    rho_minus_ones = rho_f(x, y, return_dist = TRUE)
  }
  
  rho_minus_one_mean = rho_f(x, y)
  
  rho_minus_one_perm = rep(0, R)
  
  for(i in 1:R){
    
    rho_minus_one_perm[i] = rho_f(sample(x), sample(y))
    
    if(verbose) cat("\r", round(100*i/R, 2), "%")
    
  }
  
  p_value = mean(rho_minus_one_perm > rho_minus_one_mean)
  
  return(list(rho_minus_one_mean = rho_minus_one_mean,
              rho_minus_one_sample = rho_minus_ones,
              p_value = p_value)
         )
  
}