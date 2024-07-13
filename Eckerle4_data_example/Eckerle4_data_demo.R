
# 5.04.2024

library(NISTnls)
library(energy)
library(minerva)
library(XICOR)
library(HDInterval)

source("functions/rho.R")
source("functions/rho_rank.R")
source("functions/rho_test_perm.R")

data("Eckerle4")

plot(Eckerle4$x, Eckerle4$y)

set.seed(1)

rho_minus_one_mean = rho(Eckerle4$x, Eckerle4$y)

M = 999

rho_minus_one_perm = rep(0, M)

for(i in 1:M){
  
  rho_minus_one_perm[i] = rho(sample(Eckerle4$x), sample(Eckerle4$y))
  
  cat("\r", i)
  
}

round(max(rho_minus_one_perm), 3)

test = rho_test_perm(Eckerle4$x, Eckerle4$y, verbose = T)

hist(rho_minus_one_perm, probability = T, xlim = c(0, max(c(rho_minus_one_perm, rho_minus_one_mean))),
     ylim = c(0, 15))

rho_minus_one_0 = mean(rho_minus_one_perm)

abline(v = rho_minus_one_0, lty = 2, lwd = 2, col = "red")

abline(v = rho_minus_one_mean, lwd = 2, col = "blue")

# mean,

round(rho_minus_one_mean, 3)

# Highest density interval,

round(HDInterval::hdi(rho(Eckerle4$x, Eckerle4$y,
                          return_dist = T)), 3)

round(HDInterval::hdi(test$rho_minus_one_sample), 3)

# p-value,

mean(rho_minus_one_perm > rho_minus_one_mean)

test$p_value

# Approximate p-value,

n = nrow(Eckerle4)

i = 2:(n-1)

d = 1/2 + 1/n + (1/2)*sum(1/(i*(n + 1 - i)))

#d = 0.5

rho_0_mean = 1 - (1 - d)*rho_minus_one_0

m_0 = (1 - rho_0_mean)/(1 - d)

s_0 = sqrt(rho_0_mean*(1 - rho_0_mean)/((1 - d)^2*4*n))

z = seq(0, 1, length.out = 200)

lines(z, dnorm(z, mean = m_0, sd = s_0), lwd = 2)

pnorm(rho_minus_one_mean, mean = m_0, sd = s_0, lower.tail = F)

# dcor,

energy::dcor.test(Eckerle4$x, Eckerle4$y, R = M)

# xi_cor,

XICOR::xicor(Eckerle4$x, Eckerle4$y, pvalue = T)

# MIC

minerva::mine(Eckerle4$x, Eckerle4$y)$MIC

# residuals

b1 = 1.55438
b2 = 4.08883
b3 = 451.541

res = (b1/b2)*exp(-(Eckerle4$x - b3)^2/(2*b2^2)) - Eckerle4$y

plot(Eckerle4$x, res, type = "b")
abline(h = 0)

##

rho_minus_ones = rho(Eckerle4$y, res, return_dist = T)

rho_minus_one_mean = rho(Eckerle4$y, res)

rho_minus_one_perm = rep(0, M)

for(i in 1:M){
  
  rho_minus_one_perm[i] = rho(sample(Eckerle4$y), sample(res))
  
  cat("\r", i)
  
}

hist(rho_minus_one_perm, probability = T, xlim = c(0, max(c(rho_minus_one_perm, rho_minus_one_mean))),
     ylim = c(0, 15))

rho_minus_one_0 = mean(rho_minus_one_perm)

abline(v = rho_minus_one_0, lty = 2, lwd = 2, col = "red")

abline(v = rho_minus_one_mean, lwd = 2, col = "blue")

# mean,

round(rho_minus_one_mean, 3)

# Highest density interval,

round(HDInterval::hdi(rho_minus_ones), 3)

# p-value,

mean(rho_minus_one_perm > rho_minus_one_mean)

# Approximate p-value,

rho_0_mean = 1 - (1 - d)*rho_minus_one_0

m_0 = (1 - rho_0_mean)/(1 - d)

s_0 = sqrt(rho_0_mean*(1 - rho_0_mean)/((1 - d)^2*4*n))

lines(z, dnorm(z, mean = m_0, sd = s_0), lwd = 2)

pnorm(rho_minus_one_mean, mean = m_0, sd = s_0, lower.tail = F)

# dcor,

energy::dcor.test(Eckerle4$y, res, R = M)

# xi_cor

XICOR::xicor(Eckerle4$y, res, pvalue = T)

# MIC

minerva::mine(Eckerle4$y, res)$MIC

# correlation

cor.test(Eckerle4$y, res)

