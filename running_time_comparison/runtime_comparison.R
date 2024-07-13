
# Examine running times of (i) distance correlation,
# (ii) MIC, (iii) xi_cor, (iv) naive a-r method, and 
# (v) a-r method using quicksort.

# There are no ties in data

library(dplyr)
library(XICOR)
library(minerva)
library(energy)
library(ggplot2)
library(dplyr)

source("functions/rho.R")
#source("functions/rho_rank.R")
source("functions/rho_rank_qs.R")

M = 10

N = c(1000, 2000, 5000, 10^4)

set.seed(1)

y = rnorm(max(N))
x = rnorm(max(N))

comp_time_dcor = comp_time_ar = comp_time_ar_rank = comp_time_xicor =
  comp_time_mic = matrix(0, M, length(N))

j = 0

for(n in N){
  
  i = 1
  
  j = j + 1
  
  while(i <= M){
    
    comp_time_dcor[i, j] = unname(system.time(
      energy::dcor(x[1:n], y[1:n])
    )[3])
    
    comp_time_ar[i, j] = unname(system.time(
      rho(x[1:n], y[1:n])
    )[3])
    
    comp_time_ar_rank[i, j] = unname(system.time(
      rho_rank_qs(x[1:n], y[1:n])
    )[3])
    
    comp_time_xicor[i, j] = unname(system.time(
      XICOR::calculateXI(x[1:n], y[1:n])
    )[3])
    
    comp_time_mic[i, j] = unname(system.time(
      minerva::mine(x[1:n], y[1:n])
    )[3])
    
    i = i + 1
    
  }
  
  cat("\r", n)
  
}

method = rep(c("ar", "ar (rank)", "dcor", "MIC", "xicor"),
             each = length(N)*M)

comp_times = data.frame(times = c(comp_time_ar,
                                  comp_time_ar_rank,
                                  comp_time_dcor,
                                  comp_time_mic,
                                  comp_time_xicor),
                        Method = method,
                        n = rep(rep(N, each = M), 5)
                        )

# write.table(comp_times, 
#             "running_time_comparison/runtimes.txt",
#             quote = T, row.names = F)

comp_times$Method = as.factor(comp_times$Method)

data = comp_times %>% 
  dplyr::group_by(Method, n) %>%
  dplyr::summarise(time = mean(times))

p = ggplot(data, aes(n, time, 
                     color = Method, 
                     group = Method,
                     shape = Method)) +
  geom_line(aes(linetype = Method), linewidth = 1.5) + 
  geom_point(size = 4) +
  ylab("time (seconds)") +
  theme(text = element_text(size = 30))  

p

ggsave(filename="running_time_comparison/runtimes.pdf", 
       plot = p, 
       device = cairo_pdf, 
       width = 400, 
       height = 297, 
       units = "mm")

ggsave(filename="running_time_comparison/runtimes.png", 
       plot = p, 
       device = png, 
       width = 400, 
       height = 297, 
       units = "mm",
       dpi = 700)
