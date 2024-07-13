
rho = function(x, y, return_dist = F, m = 10^4, ad_hoc_check = T, scale_monotone = T){
  
  n = length(x)
  
  d_cond = matrix(0, n, 4)
  
  d_marg = matrix(0, n, 4)
  
  for(i in 1:n){
    
    d_cond[i, 1] = mean(y[x <= x[i]] >= y[i])
    d_cond[i, 2] = mean(y[x >= x[i]] >= y[i])
    d_cond[i, 3] = mean(y[x <= x[i]] <= y[i])
    d_cond[i, 4] = mean(y[x >= x[i]] <= y[i])
    
  }
  
  #d_marg[, 1] = colMeans(outer(y, y, `>=`))
  #d_marg[, 2] = d_marg[, 1]
  #d_marg[, 3] = colMeans(outer(y, y, `<=`))
  #d_marg[, 4] = d_marg[, 3]
  
  y_rank_rev = rank(-y, ties.method = "max")
  y_rank = rank(y, ties.method = "max")
  
  d_marg[, c(1, 2)] = y_rank_rev/n
  d_marg[, c(3, 4)] = y_rank/n
  
  a = d_cond/d_marg # accept ratio
  
  a = as.vector(a)
  
  i = 2:(n-1)
  
  d = 1/2 + 1/n + (1/2)*sum(1/(i*(n + 1 - i))) # when n -> infty, d -> 1/2
  
  #d = 0.5
  
  d = scale_monotone*d
  
  if(return_dist){
    
    p = rep(0, m)
    
    for(i in 1:m){
      
      P = (1 - mean(a > runif(length(a)), na.rm = T))/(1 - d) # scaled proportion 
      
      p[i] = P - ad_hoc_check*((P > 1)*(P - 1)) 
      
    }
    
  }
  
  if(!return_dist){
    
    b = a < 1
    
    b = a*b
    
    b[b == 0] = NA
    
    p = (1 - (sum(a >= 1) + sum(b, na.rm = T))/(4*n))/(1 - d) # scaled proportion 
    
    p = p - ad_hoc_check*((p > 1)*(p - 1))
    
  }
  
  return(p)
  
}
