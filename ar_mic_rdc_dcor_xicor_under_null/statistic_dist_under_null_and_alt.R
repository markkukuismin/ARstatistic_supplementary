
# 2.8.2022

# rdc, ar, dcor and MIC distributions under null,

library(minerva) # Reshef et al (2011) DOI: 10.1126/science.120543
library(energy)
library(ggplot2)
library(ggridges)
library(dplyr)
library(XICOR) # Chatterjee (2021) DOI: 10.1080/01621459.2020.1758115
library(pwr)
library(MASS)
library(SIT)

source("functions/rdc.R")
source("functions/rho.R")

M = 1000

# Assume that the effect size is small. How different are the average values 
# under these values under the null hypothesis?

r = 0.1

res_temp = pwr::pwr.r.test(r = r, power = 0.8)

n = ceiling(res_temp$n)

ar_null = cor_null = rdc_null = dcor_null = 
  MIC_null = xicor_null = sit_null = rep(0, M)

ar_alt = cor_alt = rdc_alt = dcor_alt = 
  MIC_alt = xicor_alt = sit_alt = rep(0, M)

Sigma = diag(1, 2)

Sigma[1, 2] = Sigma[2, 1] = r

set.seed(1)

for(i in 1:M){
  
  # Under null,
  
  x = rnorm(n)
  
  y = rnorm(n)
  
  ar_null[i] = rho(x, y)
  
  cor_null[i] = cor(x, y, method = "pearson")
  
  rdc_null[i] = rdc(x, y)
  
  dcor_null[i] = dcor(x, y)
  
  MIC_null[i] = mine(x, y)$MIC
  
  xicor_null[i] = XICOR::calculateXI(x, y)
  
  sit_null[i] = SIT::sitcor(x, y, c = 16)
  
  # When effect size is weak, r = 0.1,
  
  X = mvrnorm(n = n, mu = rep(0, 2), Sigma = Sigma)
  
  x = X[, 1]
  
  y = X[, 2]
  
  ar_alt[i] = rho(x, y)
  
  cor_alt[i] = cor(x, y, method = "pearson")
  
  rdc_alt[i] = rdc(x, y)
  
  dcor_alt[i] = dcor(x, y)
  
  MIC_alt[i] = mine(x, y)$MIC

  xicor_alt[i] = max(XICOR::calculateXI(x, y), XICOR::calculateXI(y, x))
  
  sit_alt[i] = max(SIT::sitcor(x, y, c = 16), SIT::sitcor(y, x, c = 16))
  
  cat("\r", round(100*i/M, 2), "%")
  
}

Methods = rep(c("ar", "cor", "rdc", "dcor", "MIC", "xi_cor", "SIT"), each = M)

Methods = c(Methods, Methods)

value = c(ar_null, cor_null, rdc_null, dcor_null, 
          MIC_null, xicor_null, sit_null,
          ar_alt, cor_alt, rdc_alt, dcor_alt, 
          MIC_alt, xicor_alt, sit_alt)

null_or_alt = rep(c("null", "alt"), 
                  each = length(unique(Methods))*M)

Data = tibble(value = value,
              Method = Methods,
              null_or_alt = null_or_alt)

write.table(Data, 
            "ar_mic_rdc_dcor_xicor_under_null/under_null_and_alternative_results.txt", 
            row.names = F)

# Order densities so that the statistic with mean closest to zero under null is on bottom,

D = Data %>% group_by(Method, null_or_alt) %>% summarise(s_mean = mean(value))

D = D %>% filter(null_or_alt == "null")

s_order = D$Method[order(abs(D$s_mean))]

Data$Method = factor(Data$Method, levels = s_order)

main_plot = ggplot(data = Data, aes(x = value, y = Method, fill = null_or_alt)) + 
  ggridges::geom_density_ridges(alpha = 0.5, size = 0.6) + 
  geom_vline(xintercept = 0, linewidth = 1.5, linetype = "dashed") +
  scale_fill_manual(values = c("#F8766D", "white"), name = "Hypothesis", labels = c("Alternative", "Null")) + 
  xlab("Value") +
  theme(legend.position="bottom")

main_plot

main_plot = main_plot + 
  ggplot2::scale_fill_manual(values = c("gray40", "white"), name = "Hypothesis", labels = c("Alternative", "Null"))

main_plot

ggsave(filename="ar_mic_rdc_dcor_xicor_under_null/statistic_density_under_null_and_alt.pdf", 
       plot = main_plot, 
       device = cairo_pdf, 
       width = 250, 
       height = 150, 
       units = "mm")

ggsave(filename="ar_mic_rdc_dcor_xicor_under_null/statistic_density_under_null_and_alt.png", 
       plot = main_plot, 
       width = 250, 
       height = 150, 
       units = "mm")
