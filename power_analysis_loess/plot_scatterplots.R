
# 24.04.2023

# Illustrate association between x and y when y is
# generated using a simple random walk.

library(ggplot2)
library(dplyr)
library(gridExtra)
library(grid)
library(ggh4x)

n = 160

n_train = 20

n_noise = 10

n_loess = 1:100 # How many "smoothed" random walks are examined

num_noise = seq(0, 1, length.out = n_noise)

set.seed(1)

Y = list()

for(i in n_loess){
  
  # A simple random walk,
  
  Y[[i]] = sample(c(-1, 1), n_train, TRUE)
  
}

# Plot couple random walks with different noise levels,

D = c(3, 31, 93, 98)

L = c(1, 2, 5, 10)

x_train = seq(-1, 1, length.out = n_train)

y_data = c()
x_data = c()

for(d in D){
  
  for(l in L){
    
    y = cumsum(Y[[d]])
    
    train_data = data.frame(y = y, x = x_train)
    
    lo = loess(y ~ x, span = 0.5, data = train_data)
    
    x = runif(n, -1, 1)
    
    newdata = data.frame(x = x)
    
    y = predict(lo, newdata = newdata)
    
    ep = rnorm(n)
    
    y = y + 10*ep*num_noise[l]

    #plot(x, y)
    
    y_data = c(y_data, y)
    
    x_data = c(x_data , x)
    
    
  }
  
}

Noise = rep(rep(L, each = n), length(D))

data = tibble(y = y_data, 
              x = x_data,
              `Random walk no.` = rep(D, each = n*length(L)),
              `Noise level =` = Noise
              )

label_both_v2 = function(labels, multi_line = TRUE, sep = " "){
  label_both(labels, multi_line, sep)
}

p = ggplot(data, aes(x = x, y = y)) +
  geom_point() + 
  ggh4x::facet_grid2(`Random walk no.` ~ `Noise level =`, 
                     labeller = label_both_v2,
                     scales = "free_y", 
                     independent = "y") +
  theme(text = element_text(size = 20),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
        ) +
  labs(x = NULL, y = NULL)

ggsave(filename="power_analysis_loess/rw_simulations.pdf", 
       plot = p, 
       device = cairo_pdf, 
       width = 400, 
       height = 250, 
       units = "mm")

ggsave(filename="power_analysis_loess/rw_simulations.png", 
       plot = p, 
       device = png, 
       width = 400, 
       height = 250, 
       units = "mm",
       dpi = 700)
