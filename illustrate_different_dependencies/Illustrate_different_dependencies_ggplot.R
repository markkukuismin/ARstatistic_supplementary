
# 28.9.2022

# Illustrate the value of AR probability with different types of dependency 
# between x and y.

# The original script is from distance correlation wikipedia page

#Description	
# English: Recreation of File:Correlation_examples2.svg using the distance 
# correlation method (dcor from the R energy library)
# Date: 2 March 2012
# Source: Created in R using a modification of DenisBoigelot's script at File:Correlation_examples2.svg
# Author: Naught101

#Title: An example of the correlation of x and y for various distributions of (x,y) pairs
#Tags: Mathematics; Statistics; Correlation
#Author: Denis Boigelot
#Packets needed : mvtnorm (rmvnorm), RSVGTipsDevice (devSVGTips)
#How to use: output()
#
#This is an translated version in R of an Matematica 6 code by Imagecreator.

# First, the distance correlation

library(mvtnorm)
library(energy)
library(tidyverse)
library(gridExtra)

source("functions/rho.R")


MvNormal = function(n = 1000, cor_p = 0.8) {
  X = tibble(
    cor = character(),
    x = numeric(),
    y = numeric()
  )
  for (i in cor_p) {
    sd = matrix(c(1, i, i, 1), ncol = 2)
    xx = rmvnorm(n, c(0, 0), sd)
    x_res = tibble::tibble(cor = paste0(i, "MvNormal"), x = xx[, 1], y = xx[, 2])
    X = bind_rows(X, x_res)
  }
  return(X)
}

rotation = function(t, X) return(X %*% matrix(c(cos(t), sin(t), -sin(t), cos(t)), ncol = 2))

RotNormal = function(n = 1000, t = pi/2) {
  X = tibble(
    cor = character(),
    x = numeric(),
    y = numeric()
  )
  sd = matrix(c(1, 1, 1, 1), ncol = 2)
  x = rmvnorm(n, c(0, 0), sd)
  for (i in t){
    xx = rotation(i, x)
    x_res = tibble::tibble(cor = paste0("RotNormal", i), x = xx[, 1], y = xx[, 2])
    X = bind_rows(X, x_res)
  }
  return(X)
}

Others = function(n = 1000, spread = 1) {
  
  X = tibble(
    cor = character(),
    x = numeric(),
    y = numeric()
  )
  
  # x = seq(-1, 1, length.out = n) #generates a time series (serial dependence)
  x = runif(n, -1, 1) #same image, not a time series (no serial dependence)
  y = 4 * (x^2 - 1/2)^2 + runif(n, -1, 1) /3 * spread
  x_res = tibble::tibble(cor = paste0("Others", 1, spread), x = x, y = y)
  X = bind_rows(X, x_res)
  
  y = runif(n, -1, 1)
  xy = rotation(-pi/8, cbind(x,y))
  x_res = tibble::tibble(cor = paste0("Others", 2, spread), x = xy[,1], y = xy[,2]*spread)
  X = bind_rows(X, x_res)
  
  xy = rotation(-pi/8, xy)
  x_res = tibble::tibble(cor = paste0("Others", 3, spread), x = xy[,1], y = xy[,2]*spread)
  X = bind_rows(X, x_res)
  
  y = 2*x^2 + runif(n, -1, 1) * spread
  x_res = tibble::tibble(cor = paste0("Others", 4, spread), x = x, y = y)
  X = bind_rows(X, x_res)
  
  y = (x^2 + runif(n, 0, 1/2) * spread) * sample(seq(-1, 1, 2), n, replace = TRUE)
  x_res = tibble::tibble(cor = paste0("Others", 5, spread), x = x, y = y)
  X = bind_rows(X, x_res)
  
  y = cos(x*pi) + rnorm(n, 0, 1/8) * spread
  x = sin(x*pi) + rnorm(n, 0, 1/8) * spread
  x_res = tibble::tibble(cor = paste0("Others", 6, spread), x = x, y = y)
  X = bind_rows(X, x_res)
  
  xy1 = rmvnorm(n/4, c( 3,  3), diag(2)*spread)
  xy2 = rmvnorm(n/4, c(-3,  3), diag(2)*spread)
  xy3 = rmvnorm(n/4, c(-3, -3), diag(2)*spread)
  xy4 = rmvnorm(n/4, c( 3, -3), diag(2)*spread)
  x = c(xy1[, 1], xy2[, 1], xy3[, 1], xy4[, 1])
  y = c(xy1[, 2], xy2[, 2], xy3[, 2], xy4[, 2])
  x_res = tibble::tibble(cor = paste0("Others", 7, spread), x = x, y = y)
  X = bind_rows(X, x_res)
  
  return(X)
  
}

set.seed(1)

d1 = MvNormal(800, c(1.0, 0.8, 0.4, 0.0, -0.4, -0.8, -1.0))
d2 = RotNormal(200, c(0, pi/12, pi/6, pi/4, pi/2-pi/6, pi/2-pi/12, pi/2))
d3 = Others(800)
d4 = Others(800, 0.3)

Data = bind_rows(d1, d2, d3, d4)

# Running one realization of the A-R procedure is faster,

c_ar = Data %>% 
  group_by(cor) %>% 
  summarise(p = rho(x,y)) %>%
  add_column(method = "A-R")

c_dcor = Data %>% 
  group_by(cor) %>% 
  #summarise(p = round(dcor(x,y), 1)) %>% 
  summarise(p = dcor(x,y)) %>% 
  add_column(method = "dcor")

Data_ar = dplyr::left_join(Data, c_ar)

Data_dcor = dplyr::left_join(Data, c_dcor)

Data = bind_rows(Data_ar, Data_dcor)

# The next modification is a silly one: when the normal data is rotated 45 degrees
# clockwise; RotNormal(200, t = pi/4). This corresponds to the horizontal line of
# the original figure found also on the Wikipage. But, when the simulated data is
# extracted and the correlation is calculated, the corresponding correlation is not actually
# zero but something else,

xy = RotNormal(200, t = pi/4)

xy = xy[, 2:3]

cor(xy)

# See also the third line in the function in "Myplot" on the wikipage,
#MyPlot <- function(xy, xlim = c(-4, 4), ylim = c(-4, 4), eps = 1e-15) {
#  title = round(cor(xy[,1], xy[,2]), 1)
#  if (sd(xy[,2]) < eps) title = "" # corr. coeff. is undefined

# The correlation is just treated as zero, because the other random variable is 
# treated as constant - although it strictly speaking isn't,

apply(xy, 2, function(x)round(sd(x), 2))

# The A-R accept ratio is also zero, if we treat the other random variable as a
# constant (if you can call it in this case a random variable...),

x_test = rnorm(200)

y_test = rep(0, 200)

rho(x_test, y_test)
rho(y_test, x_test)
cor(x_test, y_test)


# Thus, the results are modified to make them be in line with the basic idea
# of the original wiki figure,

Data[str_detect(Data$cor, regex("RotNormal0.785", ignore_case = TRUE)), "p"] = 0

# The rest of the code is for plotting the image

plot_list = list()

diff_dep = unique(Data$cor)

for(i in 1:length(diff_dep)){
  
  AR_cor = filter(Data, cor == diff_dep[i], method == "A-R") %>% dplyr::select(p)
  AR_cor = AR_cor$p[1]
  AR_cor = format(round(AR_cor, digits = 2), nsmall = 2)
  
  dcor_cor = filter(Data, cor == diff_dep[i], method == "dcor") %>% dplyr::select(p)
  dcor_cor = dcor_cor$p[1]
  dcor_cor = format(round(dcor_cor, digits = 2), nsmall = 2)
  
  RN = str_detect(diff_dep[i], regex("RotNormal", ignore_case = TRUE))
  
  RN2 = diff_dep[i] %in% c("Others20.3", "Others30.3")
  
  plot_list[[i]] = filter(Data, cor == diff_dep[i]) %>% 
    ggplot(data = ., aes(x = x, y = y)) +
    geom_point(size = 1, alpha = 1/10) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          plot.title = element_text(hjust = 0.5)
    ) +
    {if(RN) xlim(-4, 4)} +
    {if(RN) ylim (-4, 4)} +
    {if(RN2) xlim(-1.4, 1.4)} +
    {if(RN2) ylim (-1, 1)} +
    labs(title = bquote(~bold(.(AR_cor))~"("*.(dcor_cor)*")"))
  
}

main_plot = do.call("grid.arrange", c(plot_list, nrow = 4))

ggsave(filename="illustrate_different_dependencies/AR_and_dcor.png",
       plot = main_plot, 
       device = png,
       width = 250, 
       height = 150, 
       units = "mm")

ggsave(filename="illustrate_different_dependencies/AR_and_dcor.pdf",
       plot = main_plot, 
       device = cairo_pdf,
       width = 250, 
       height = 150, 
       units = "mm")
