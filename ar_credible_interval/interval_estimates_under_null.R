
# 2.4.2023

# Illustrate how the sample size affect to the
# ar association statistic values under null

library(tidyverse)
library(energy)
library(SIT)
library(XICOR)
library(minerva)

source("functions/rho_rank.R")
source("functions/rdc.R")

ar_hdi = function(x, alpha = 0.95){
  
  c(HDInterval::hdi(x, credMass = alpha))
  
}

d_cor = function(data, indices){
  
  return(energy::dcor(data[indices, 1], data[indices, 2]))
  
}

sit_cor = function(data, indices){
  
  return(SIT::sitcor(data[indices, 1], data[indices, 2], c = 16))
  
}

xi_cor = function(data, indices){
  
  return(XICOR::xicor(data[indices, 1], data[indices, 2]))
  
}

mic_cor = function(data, indices){
  
  return(minerva::mine(data[indices, 1], data[indices, 2])$MIC)
  
}

rdc_cor = function(data, indices){
  
  return(rdc(data[indices, 1], data[indices, 2]))
  
}

cor_ci = function(bs_cor, alpha = 0.95){
  
  c(boot::boot.ci(bs_cor, conf = alpha, type = "perc")$percent[1, 4],
    boot::boot.ci(bs_cor, conf = alpha, type = "perc")$percent[1, 5])
  
}

set.seed(1)

test = rho_rank(runif(100), runif(100), return_dist = TRUE)

ar_hdi(test)

n = 1000

x = rnorm(n)

y = rnorm(n)

N = c(20, 50, 100, 200, 500, 1000)

l_u = rep(c("lower", "upper"), length(N))

ar_temp = sapply(N, function(m) rho_rank(x[1:m], y[1:m], return_dist = TRUE))
    
ar_mean = sapply(N, function(m) rho_rank(x[1:m], y[1:m]))
cor_m = sapply(N, function(m) cor(x[1:m], y[1:m]))
dcor_m = sapply(N, function(m) energy::dcor(x[1:m], y[1:m]))
sitcor_m = sapply(N, function(m) SIT::sitcor(x[1:m], y[1:m], c = 16))
mic_m = sapply(N, function(m) minerva::mine(x[1:m], y[1:m])$MIC)
xicor_m = sapply(N, function(m) xicor(x[1:m], y[1:m]))
rdc_m = sapply(N, function(m) rdc(x[1:m], y[1:m]))
    
ar_res = apply(ar_temp, 2, ar_hdi)

dcor_res = cor_res = sitcor_res = mic_res = xicor_res = rdc_res = 
  matrix(0, nrow = 2, ncol = length(N))

rownames(dcor_res) = rownames(cor_res) = 
  rownames(sitcor_res) = rownames(mic_res) = 
  rownames(xicor_res) = rownames(rdc_res) = c("lower", "upper")

v = 1

for(kk in N){
  
  dcor_temp = boot::boot(data = cbind(x[1:kk], y[1:kk]), statistic = d_cor, R = 1000)
  cor_temp = cor.test(x[1:kk], y[1:kk])
  sitcor_temp = boot::boot(data = cbind(x[1:kk], y[1:kk]), statistic = sit_cor, R = 1000)
  xicor_temp = boot::boot(data = cbind(x[1:kk], y[1:kk]), statistic = xi_cor, R = 1000)
  rdc_temp = boot::boot(data = cbind(x[1:kk], y[1:kk]), statistic = rdc_cor, R = 1000)
  mic_temp = boot::boot(data = cbind(x[1:kk], y[1:kk]), statistic = mic_cor, R = 1000)
  
  dcor_res[, v] = cor_ci(dcor_temp) 
  cor_res[, v] = c(cor_temp$conf.int)
  sitcor_res[, v] = cor_ci(sitcor_temp)
  xicor_res[, v] = cor_ci(xicor_temp)
  rdc_res[, v] = cor_ci(rdc_temp)
  mic_res[, v] = cor_ci(mic_temp)
  
  v = v + 1
  
}
  
ar_df = data.frame(value = ar_mean,
                   lower = ar_res["lower", ],
                   upper = ar_res["upper", ],
                   n = N,
                   Method = "ar")

cor_df = data.frame(value = cor_m,
                    lower = cor_res["lower", ],
                    upper = cor_res["upper", ],
                    n = N,
                    Method = "cor")

dcor_df = data.frame(value = dcor_m,
                      lower = dcor_res["lower", ],
                      upper = dcor_res["upper", ],
                      n = N,
                      Method = "dcor")

sitcor_df = data.frame(value = sitcor_m,
                       lower = sitcor_res["lower", ],
                       upper = sitcor_res["upper", ],
                       n = N,
                       Method = "SIT")

xicor_df = data.frame(value = xicor_m,
                       lower = xicor_res["lower", ],
                       upper = xicor_res["upper", ],
                       n = N,
                       Method = "xicor")

rdc_df = data.frame(value = rdc_m,
                      lower = rdc_res["lower", ],
                      upper = rdc_res["upper", ],
                      n = N,
                      Method = "rdc")

mic_df = data.frame(value = mic_m,
                    lower = mic_res["lower", ],
                    upper = mic_res["upper", ],
                    n = N,
                    Method = "MIC")

df = dplyr::bind_rows(ar_df, cor_df, dcor_df, sitcor_df, xicor_df, 
                      rdc_df, mic_df)

##

df$n = as.factor(df$n)
df$Method = as.factor(df$Method)

p = ggplot(df, 
           aes(x = n,
               y = value, 
               ymin = lower,
               ymax = upper,
               color = Method)
) +
  geom_point(size = 1, position = position_dodge(0.25)) +
  geom_linerange(position = position_dodge(0.25))

p + xlab("Sample size") +
  ylab("Value under null") +
  ylim(-0.45, 1)

##

p = ggplot(df, aes(x = n, y = value, color = Method)) +
  geom_point(size = 2, position = position_dodge(0.25)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2,
                position = position_dodge(0.25),
                linewidth = 1)

p = p +
  xlab("Sample size") +
  ylab("Value under null") +
  ylim(-0.45, 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(text = element_text(size = 20))

p

ggsave("ar_credible_interval/interval_estimates_null.pdf",
       plot = p, 
       device = cairo_pdf, 
       width = 400, 
       height = 297, 
       units = "mm")

ggsave(filename="ar_credible_interval/interval_estimates_null.png", 
       plot = p, 
       device = png, 
       width = 400, 
       height = 297, 
       units = "mm",
       dpi = 700)
