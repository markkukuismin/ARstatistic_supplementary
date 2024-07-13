
# 9.4.2024

# Careless use of bootstrap might distort the relationship between variables

library(minerva)
library(energy)
library(XICOR) 
library(ggplot2)
library(dplyr)
library(SIT)

source("functions/rho.R")
source("functions/rho_rank.R")
source("functions/rho_test_perm.R")
source("functions/rdc.R")

x = 1:10

y = rep(x, each = length(x))

x = rep(x, times = length(x))

n = length(x)

grid_plot = ggplot(data = data.frame(x = x, y = y), aes(x = x, y = y)) +
  geom_point() +
  theme(axis.title=element_text(size=14, face="bold"))

grid_plot

ggsave(filename="bootstrap_example/grid.pdf", 
       plot = grid_plot, 
       device = cairo_pdf, 
       width = 400, 
       height = 297, 
       units = "mm")

ggsave(filename="bootstrap_example/grid.png", 
       plot = grid_plot, 
       width = 250, 
       height = 150, 
       units = "mm")

# Everything is OK when we test the independence of the variables,

set.seed(1)

dc_hat = energy::dcor.test(x, y, R = 500)

dc_hat

dc_pvalue = dc_hat$p.value

dc_hat = dc_hat$estimates["dCor"]

MIC_hat = minerva::mine(x, y)$MIC

rdc_hat = rdc(x, y)

xi_hat = XICOR::xicor(x, y, pvalue = T)

xi_hat

xi_pvalue = xi_hat$pval

xi_hat = xi_hat$xi

sit_hat = sitcor(x, y, pvalue = T)

sit_hat$sitcor

sit_hat$pval

d = rho_test_perm(x, y, R = 1000, use_rank = FALSE)
d$rho_minus_one_mean
d$p_value

# What about if we want to estimate confidence intervals?

dist_cor = function(data, indices){
  
  return(energy::dcor(data[indices, 1], data[indices, 2]))
  
}

xi_cor = function(data, indices){
  
  return(XICOR::xicor(data[indices, 1], data[indices, 2]))
  
}

MIC = function(data, indices){
  
  return(minerva::mine(data[indices, 1], data[indices, 2])$MIC)
  
}

sit_cor = function(data, indices){
  
  return(SIT::sitcor(data[indices, 1], data[indices, 2]))
  
}

rdc_b = function(data, indices){
  k=20
  s=1/6
  x = data[indices, 1]
  y = data[indices, 2]
  x <- cbind(apply(as.matrix(x),2,function(u)rank(u)/length(u)),1)
  y <- cbind(apply(as.matrix(y),2,function(u)rank(u)/length(u)),1)
  x <- s/ncol(x)*x%*%matrix(rnorm(ncol(x)*k),ncol(x))
  y <- s/ncol(y)*y%*%matrix(rnorm(ncol(y)*k),ncol(y))
  cancor(cbind(sin(x),1),cbind(sin(y),1))$cor[1]
}

bs_dcor = boot::boot(data = cbind(x, y), statistic = dist_cor, R=1000)
bs_xicor = boot::boot(data = cbind(x, y), statistic = xi_cor, R=1000)
bs_MIC = boot::boot(data = cbind(x, y), statistic = MIC, R=1000)
bs_rdc = boot::boot(data = cbind(x, y), statistic = rdc_b, R=1000)
bs_sit = boot::boot(data = cbind(x, y), statistic = sit_cor, R=1000)

CI = data.frame(method = c("dcor", "xi_cor", "MIC", "rdc", "SIT", "AR"), 
                point.estimates = c(dc_hat, xi_hat, MIC_hat, rdc_hat, sit_hat$sitcor, d$rho_minus_one_mean),
                lower.bound = c(
                  boot::boot.ci(bs_dcor, conf = 0.95, type = "perc")$percent[1, 4],
                  boot::boot.ci(bs_xicor, conf = 0.95, type = "perc")$percent[1, 4],
                  boot::boot.ci(bs_MIC, conf = 0.95, type = "perc")$percent[1, 4],
                  boot::boot.ci(bs_rdc, conf = 0.95, type = "perc")$percent[1, 4],
                  boot::boot.ci(bs_sit, conf = 0.95, type = "perc")$percent[1, 4],
                  HDInterval::hdi(d$rho_minus_one_sample, credMass = 0.95)[1]
                ),
                upper.bound = c(
                  boot::boot.ci(bs_dcor, conf = 0.95, type = "perc")$percent[1, 5],
                  boot::boot.ci(bs_xicor, conf = 0.95, type = "perc")$percent[1, 5],
                  boot::boot.ci(bs_MIC, conf = 0.95, type = "perc")$percent[1, 5],
                  boot::boot.ci(bs_rdc, conf = 0.95, type = "perc")$percent[1, 5],
                  boot::boot.ci(bs_sit, conf = 0.95, type = "perc")$percent[1, 5],
                  HDInterval::hdi(d$rho_minus_one_sample, credMass = 0.95)[2]
                ),
                p.values = c(dc_pvalue, xi_pvalue, NA, NA, sit_hat$pval, d$p_value)
                )

CI %>% dplyr::mutate_if(is.numeric, ~round(., 3)) # Not in line with the data...

cor.test(x, y)

# Add some noise and re-run the analysis,

x = x + rnorm(n, sd = 0.5)
y = y + rnorm(n, sd = 0.5)

grid_plot_noise = ggplot(data = data.frame(x = x, y = y), aes(x = x, y = y)) +
  geom_point() +
  theme(axis.title=element_text(size=14, face="bold"))

grid_plot_noise

ggsave(filename="bootstrap_example/grid_noise.pdf", 
       plot = grid_plot_noise, 
       device = cairo_pdf, 
       width = 400, 
       height = 297, 
       units = "mm")

ggsave(filename="bootstrap_example/grid_noise.png", 
       plot = grid_plot_noise, 
       width = 250, 
       height = 150, 
       units = "mm")

dc_hat = energy::dcor.test(x, y, R = 500)

dc_hat

dc_pvalue = dc_hat$p.value

dc_hat = dc_hat$estimates["dCor"]

MIC_hat = minerva::mine(x, y)$MIC

rdc_hat = rdc(x, y)

xi_hat = XICOR::xicor(x, y, pvalue = T)

xi_hat

xi_pvalue = xi_hat$pval

xi_hat = xi_hat$xi

sit_hat = sitcor(x, y, pvalue = TRUE)

d = rho_test_perm(x, y, R = 1000, use_rank = FALSE)
d$rho_minus_one_mean
d$p_value

bs_dcor = boot::boot(data = cbind(x, y), statistic = dist_cor, R=1000)
bs_xicor = boot::boot(data = cbind(x, y), statistic = xi_cor, R=1000)
bs_MIC = boot::boot(data = cbind(x, y), statistic = MIC, R=1000)
bs_rdc = boot::boot(data = cbind(x, y), statistic = rdc_b, R=1000)
bs_sit = boot::boot(data = cbind(x, y), statistic = sit_cor, R=1000)

CI_noise = data.frame(method = c("dcor", "xi_cor", "MIC", "rdc", "SIT", "AR"), 
                      point.estimates = c(dc_hat, xi_hat, MIC_hat, rdc_hat, sit_hat$sitcor, d$rho_minus_one_mean),
                      lower.bound = c(
                        boot::boot.ci(bs_dcor, conf = 0.95, type = "perc")$percent[1, 4],
                        boot::boot.ci(bs_xicor, conf = 0.95, type = "perc")$percent[1, 4],
                        boot::boot.ci(bs_MIC, conf = 0.95, type = "perc")$percent[1, 4],
                        boot::boot.ci(bs_rdc, conf = 0.95, type = "perc")$percent[1, 4],
                        boot::boot.ci(bs_sit, conf = 0.95, type = "perc")$percent[1, 4],
                        HDInterval::hdi(d$rho_minus_one_sample, credMass = 0.95)[1]
                      ),
                      upper.bound = c(
                        boot::boot.ci(bs_dcor, conf = 0.95, type = "perc")$percent[1, 5],
                        boot::boot.ci(bs_xicor, conf = 0.95, type = "perc")$percent[1, 5],
                        boot::boot.ci(bs_MIC, conf = 0.95, type = "perc")$percent[1, 5],
                        boot::boot.ci(bs_rdc, conf = 0.95, type = "perc")$percent[1, 5],
                        boot::boot.ci(bs_sit, conf = 0.95, type = "perc")$percent[1, 5],
                        HDInterval::hdi(d$rho_minus_one_sample, credMass = 0.95)[2]
                      ),
                      p.values = c(dc_pvalue, xi_pvalue, NA, NA, sit_hat$pval, d$p_value)
                      )

CI_noise %>% dplyr::mutate_if(is.numeric, ~round(., 3)) # Not in line with the data

# Bivariate normal mixture (see Chatterjee 2021),

n = 200

X1 = MASS::mvrnorm(n, mu = rep(0, 2), diag(1, 2))

X2 = MASS::mvrnorm(n, mu = rep(5, 2), diag(1, 2))

p = rbinom(n, 1, prob = 0.5)

X = p*X1 + (1 - p)*X2

x = X[, 1]
y = X[, 2]

ggplot(data = data.frame(x = x, y = y), aes(x = x, y = y)) +
  geom_point() +
  theme(axis.title=element_text(size=14, face="bold"))

dc_hat = energy::dcor.test(x, y, R = 500)

dc_hat

dc_pvalue = dc_hat$p.value

dc_hat = dc_hat$estimates["dCor"]

MIC_hat = minerva::mine(x, y)$MIC

rdc_hat = rdc(x, y)

xi_hat = XICOR::xicor(x, y, pvalue = T)

xi_hat

xi_pvalue = xi_hat$pval

xi_hat = xi_hat$xi

sit_hat = sitcor(x, y, pvalue = TRUE)

sit_hat$sitcor

sit_hat$pval

d = rho_test_perm(x, y, R = 1000, use_rank = TRUE)
d$rho_minus_one_mean
d$p_value

bs_dcor = boot::boot(data = cbind(x, y), statistic = dist_cor, R=1000)
bs_xicor = boot::boot(data = cbind(x, y), statistic = xi_cor, R=1000)
bs_MIC = boot::boot(data = cbind(x, y), statistic = MIC, R=500)
bs_rdc = boot::boot(data = cbind(x, y), statistic = rdc_b, R=1000)
bs_sit = boot::boot(data = cbind(x, y), statistic = sit_cor, R=1000)

CI_bi_norm = data.frame(method = c("dcor", "xi_cor", "MIC", "rdc", "SIT", "AR"),
                        point.estimates = c(dc_hat, xi_hat, MIC_hat, rdc_hat, sit_hat$sitcor, d$rho_minus_one_mean),
                        lower.bound = c(
                          boot::boot.ci(bs_dcor, conf = 0.95, type = "perc")$percent[1, 4],
                          boot::boot.ci(bs_xicor, conf = 0.95, type = "perc")$percent[1, 4],
                          boot::boot.ci(bs_MIC, conf = 0.95, type = "perc")$percent[1, 4],
                          boot::boot.ci(bs_rdc, conf = 0.95, type = "perc")$percent[1, 4],
                          boot::boot.ci(bs_sit, conf = 0.95, type = "perc")$percent[1, 4],
                          HDInterval::hdi(d$rho_minus_one_sample, credMass = 0.95)[1]
                        ),
                        upper.bound = c(
                          boot::boot.ci(bs_dcor, conf = 0.95, type = "perc")$percent[1, 5],
                          boot::boot.ci(bs_xicor, conf = 0.95, type = "perc")$percent[1, 5],
                          boot::boot.ci(bs_MIC, conf = 0.95, type = "perc")$percent[1, 5],
                          boot::boot.ci(bs_rdc, conf = 0.95, type = "perc")$percent[1, 5],
                          boot::boot.ci(bs_sit, conf = 0.95, type = "perc")$percent[1, 5],
                          HDInterval::hdi(d$rho_minus_one_sample, credMass = 0.95)[2]
                        ),
                        p.values = c(dc_pvalue, xi_pvalue, NA, NA, sit_hat$pval, d$p_value)
                        )



CI_bi_norm %>% dplyr::mutate_if(is.numeric, ~round(., 3)) # Not in line with the data

