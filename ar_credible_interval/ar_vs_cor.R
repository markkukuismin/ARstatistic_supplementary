
# 12.4.2024

library(foreach)
library(doParallel)
library(ggplot2)
library(tidyverse)

cl = makePSOCKcluster(5)
registerDoParallel(cl)

source("functions/rho_rank.R")

set.seed(1)

r = c(0, 0.01, 0.15, 0.3, 0.5, 0.7, 0.9, 1)

n = 150

M = 100

LU = c()

for(j in 1:length(r)){
  
  Sigma = matrix(r[j], 2, 2)
  
  diag(Sigma) = 1
  
  LU_temp <- foreach(i = 1:M, .combine = rbind) %dopar% {
    
    X = MASS::mvrnorm(n, mu = rep(0, 2), Sigma = Sigma)
    
    x = X[, 1]
    y = X[, 2]
    
    lu = rho_rank(x, y, return_dist = T)
    
    lu = HDInterval::hdi(lu, credMass = 0.95)
    
    c(lu, rho_rank(x, y)) 
    
  }
  
  LU = cbind(LU, LU_temp)
  
}

#stopCluster(cl)

CIs = colMeans(LU)

u_ind = stringr::str_detect(names(CIs), "[u]")
l_ind = stringr::str_detect(names(CIs), "[l]")

df_ar = data.frame(low = CIs[l_ind], 
                   up = CIs[u_ind], 
                   est = CIs[!u_ind & !l_ind], 
                   Method = "ar",
                   cor = r)

##

LU = c()

for(j in 1:length(r)){
  
  Sigma = matrix(r[j], 2, 2)
  
  diag(Sigma) = 1
  
  LU_temp <- foreach(i = 1:M, .combine = rbind) %dopar% {
    
    X = MASS::mvrnorm(n, mu = rep(0, 2), Sigma = Sigma)
    
    x = X[, 1]
    y = X[, 2]
    
    lu = c(cor.test(x, y)$conf.int)
    
    c(lu, cor(x, y)) 
    
  }
  
  LU = cbind(LU, LU_temp)
  
}

stopCluster(cl)

colnames(LU) = rep(c("low", "up", "est"), length(r))

CIs = colMeans(LU)

u_ind = stringr::str_detect(names(CIs), "[u]")
l_ind = stringr::str_detect(names(CIs), "[l]")

df_cor = data.frame(low = CIs[l_ind], 
                    up = CIs[u_ind], 
                    est = CIs[!u_ind & !l_ind], 
                    Method = "cor",
                    cor = r)

df = rbind(df_ar, df_cor)

df = df %>% 
  group_by(Method, cor) %>%
  mutate(width = up - low) %>%
  mutate(width = format(round(width, digits = 3), nsmall = 3) )

df$cor = as.factor(df$cor)

p = ggplot(df, aes(cor, est)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = low, ymax = up), width=.1) +
  facet_grid(cols = vars(Method)) +
  geom_text(aes(label = width, y = up + 0.05), size = 5) +
  xlab("Population correlation value") +
  ylab("Estimate") +
  theme(text = element_text(size = 20))

p

pdf_file = paste0("ar_credible_interval/ar_vs_cor_", n, ".pdf")

ggsave(filename = pdf_file,
       plot = p, 
       device = cairo_pdf, 
       width = 400, 
       height = 297, 
       units = "mm")

png_file = paste0("ar_credible_interval/ar_vs_cor_", n, ".png")

ggsave(filename = png_file, 
       plot = p, 
       device = png, 
       width = 400, 
       height = 297, 
       units = "mm",
       dpi = 700)
