
clusters = function(n, s){
  
  ex = s*rnorm(n)
  ey = s*rnorm(n)
  cx = sample(c(-1, 1), size = n, replace = T)
  cy = sample(c(-1, 1), size = n, replace = T)
  
  ind = sample(c(1, 0), size = n, replace = T)
  
  cx = cx*ind
  cy = cy*ind
  
  u = cx + ex
  v = cy + ey
  
  return(list(x = u, y = v))
  
}