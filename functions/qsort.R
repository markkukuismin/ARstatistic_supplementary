qsort = function(y){
  N = length(y)
  if(N <= 1){
    y
  }else{
    m = sample(1:N, 1)
    k = 1
    b = 0
    v1 = numeric()
    v2 = numeric()
    vm = numeric()
    while(k <= N){
      if(y[k] < y[m]){
        v1 = c(v1, y[k])
        b = b + 1
      }
      if(y[k] == y[m]){
        vm = c(vm, y[k])
        z[y[k]] <<- z[y[k]] + b
        b = b + 1
      }
      if(y[k] > y[m]){
        v2 = c(v2, y[k])
        z[y[k]] <<- z[y[k]] + b
      }
      k = k + 1
    }
    v1 = qsort(v1)
    v2 = qsort(v2)
    y = c(v1, vm, v2)
    y
  }
}
