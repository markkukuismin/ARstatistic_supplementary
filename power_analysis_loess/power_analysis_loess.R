
# 05.03.2024

# Examine the power compared to other measures, such as (i) Pearson correlation,
# (ii) distance correlation (iii) MIC, (iv) Randomized Dependence Coefficient
# (RDC), (v) xi_cor, (vi) Hoeffding, (vii) Refined Hoeffding,
# (viii) Bergsma-Dassios-Yanagimoto coefficient (BDY), (ix) S_ADP test, (x) dHSIC,
# (xi) Sliced Independence Test (SIT)

library(minerva)
library(energy)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(grid)
library(XICOR) # Chatterjee (2021) DOI: 10.1080/01621459.2020.1758115
library(independence) # Removed from CRAN (last visit 3/5/2024), package available from the archive
library(TauStar) # Removed from CRAN (last visit 3/5/2024), package available from the archive
library(SIT)
library(HHG)
library(dHSIC)

source("functions/rdc.R")
source("functions/rho_rank.R")

# We write the following function because MIC glitches a small percentage of 
# the time, and we do not wish to average over those trials

notNA.greater <- function(a,b){
  ind <- which(!is.na(a))
  pow <- sum(a[ind] > b)/length(ind)
  return(pow)
}

n = 160

n_train = 20

n_noise = 10

num_noise = seq(0, 1, length.out = n_noise)

M = 500

val_null_ar = val_null_cor = val_null_dcor = 
  val_null_mine = val_null_rdc = val_null_xicor = 
  val_null_Hoeff = val_null_RHoeff = val_null_BDY =
  val_null_dHSIC = val_null_SIT = rep(NA, M)  # Vectors holding the null "correlations"

val_ar = val_cor = val_dcor = 
  val_mine = val_rdc = val_xicor = 
  val_Hoeff = val_RHoeff = val_BDY =
  val_Sadp = val_dHSIC = val_SIT = rep(NA, M) # Vectors holding the alternative "correlations"

n_loess = 1:100 # How many "smoothed" random walks are examined

power_ar = power_cor = power_dcor = power_mine = power_rdc = power_xicor = 
  power_Hoeff = power_RHoeff = power_BDY = power_Sadp = power_dHSIC = 
  power_SIT = array(NA, c(length(n_loess), n_noise)) 

set.seed(1)

# Null table for the S_ADP test. Create null table for aggregation by summation,

null_table_mxl = HHG::Fast.independence.test.nulltable(n, variant = 'ADP-EQP-ML', 
                                                       nr.atoms = 40, nr.perm = 1000)

Y = list()

for(i in n_loess){
  
  # A simple random walk,
  
  Y[[i]] = sample(c(-1, 1), n_train, TRUE)
  
}

k = 1

for(l in 1:length(num_noise)){
  
  j = 1
  
  for(d in n_loess){
    
    for(i in 1:M){
      
      x_train = seq(-1, 1, length.out = n_train)
      
      y = cumsum(Y[[d]])
      
      train_data = data.frame(y = y, x = x_train)
      
      lo = loess(y ~ x, span = 0.5, data = train_data)
      
      x = runif(n, -1, 1)
      
      newdata = data.frame(x = x)
      
      y = predict(lo, newdata = newdata)
      
      ep = rnorm(n)
      
      y = y + 10*ep*num_noise[l]
      
      # Under null,
      
      x = runif(n, -1, 1)
      
      val_null_ar[i] = rho_rank(x, y)
      val_null_cor[i] = (cor(x,y))^2            # Calculate the correlation
      val_null_dcor[i] = dcor(x,y)              # Calculate dcor
      val_null_mine[i] = minerva::mine(x, y)$MIC # Calculate mic
      val_null_rdc[i] = rdc(x, y)
      #val_null_xicor[i] = xicor(x, y) # running xicor messes the seed number...
      val_null_xicor[i] = XICOR::calculateXI(x, y)
      Dn = independence::hoeffding.D.test(x, y, precision = 1)$Dn
      val_null_Hoeff[i] = Dn
      Rn = independence::hoeffding.refined.test(x, y, precision = 1)$Rn
      val_null_RHoeff[i] = Rn
      val_null_BDY[i] = 12*(Dn + 2*Rn) #independence::tau.star.test(x, y, precision = 1)$Tn
      val_null_dHSIC[i] = dHSIC::dhsic(x, y)$dHSIC
      val_null_SIT[i] = SIT::sitcor(x, y, c = 16)
      
    }
    
    # we remove the mic trials which glitch
    
    val_null_mine = val_null_mine[which(!is.na(val_null_mine))]
    
    ## Next we calculate rejection cutoffs
    
    cut_ar = quantile(val_null_ar, .95)
    cut_cor = quantile(val_null_cor, .95)
    cut_dcor = quantile(val_null_dcor, .95)
    cut_mine = quantile(val_null_mine, .95)
    cut_rdc = quantile(val_null_rdc, .95)
    cut_xicor = quantile(val_null_xicor, .95)
    cut_Hoeff = quantile(val_null_Hoeff, .95)
    cut_RHoeff = quantile(val_null_RHoeff, .95)
    cut_BDY = quantile(val_null_BDY, .95)
    cut_dHSIC = quantile(val_null_dHSIC, .95)
    cut_SIT = quantile(val_null_SIT, .95)
    
    # Next we simulate the data again, this time under the alternative
    
    for(i in 1:M){
      
      y = cumsum(Y[[d]])
      
      lo = loess(y ~ x, span = 0.5, data = train_data)
      
      x = runif(n, -1, 1)
      
      newdata = data.frame(x = x)
      
      y = predict(lo, newdata = newdata)
      
      ep = rnorm(n)
      
      y = y + 10*ep*num_noise[l]
      
      val_ar[i] = rho_rank(x, y)
      val_cor[i] = (cor(x,y))^2            # Calculate the correlation
      val_dcor[i] = dcor(x,y)              # Calculate dcor
      val_mine[i] = minerva::mine(x, y)$MIC # Calculate mic
      val_rdc[i] = rdc(x, y)
      #val_xicor[i] = xicor(x, y)
      val_xicor[i] = XICOR::calculateXI(x, y)
      Dn = independence::hoeffding.D.test(x, y, precision = 1)$Dn
      val_Hoeff[i] = Dn
      Rn = independence::hoeffding.refined.test(x, y, precision = 1)$Rn
      val_RHoeff[i] = Rn
      val_BDY[i] = 12*(Dn + 2*Rn) #independence::tau.star.test(x, y, precision = 1)$Tn
      val_Sadp[i] = HHG::Fast.independence.test(x, y, null_table_mxl, combining.type = "MinP")$MinP.pvalue
      val_dHSIC[i] = dHSIC::dhsic(x, y)$dHSIC
      val_SIT[i] = SIT::sitcor(x, y, c = 16)
      
    }
    
    ## Now we estimate the power as the number of alternative statistics exceeding our estimated cutoffs
    
    power_ar[j, l] = sum(val_ar > cut_ar)/M
    power_cor[j, l] = sum(val_cor > cut_cor)/M
    power_dcor[j, l] = sum(val_dcor > cut_dcor)/M
    power_mine[j, l] = notNA.greater(val_mine, cut_mine)
    power_rdc[j, l] = sum(val_rdc > cut_rdc)/M
    power_xicor[j, l] = sum(val_xicor > cut_xicor)/M
    power_Hoeff[j, l]  = sum(val_Hoeff > cut_Hoeff)/M
    power_RHoeff[j, l] = sum(val_RHoeff > cut_RHoeff)/M
    power_BDY[j, l] = sum(val_BDY > cut_BDY)/M
    power_Sadp[j, l] = sum(val_Sadp < 0.05)/M
    power_dHSIC[j, l] = sum(val_dHSIC > cut_dHSIC)/M
    power_SIT[j, l] = sum(val_SIT > cut_SIT)/M
    
    j = j + 1
    
  }
  
  cat("\r", l)
  
}

## The rest of the code is for saving results

ar_d = data.frame(power = c(t(power_ar)), noise = rep(1:n_noise, times = length(n_loess)), 
                  n_loess = rep(n_loess, each = n_noise), method = "ar")

cor_d = data.frame(power = c(t(power_cor)), noise = rep(1:n_noise, times = length(n_loess)), 
                   n_loess = rep(n_loess, each = n_noise), method = "cor")

dcor_d = data.frame(power = c(t(power_dcor)), noise = rep(1:n_noise, times = length(n_loess)), 
                    n_loess = rep(n_loess, each = n_noise), method = "dcor")

MIC_d = data.frame(power = c(t(power_mine)), noise = rep(1:n_noise, times = length(n_loess)), 
                   n_loess = rep(n_loess, each = n_noise), method = "MIC")

rdc_d = data.frame(power = c(t(power_rdc)), noise = rep(1:n_noise, times = length(n_loess)), 
                   n_loess = rep(n_loess, each = n_noise), method = "rdc")

xicor_d = data.frame(power = c(t(power_xicor)), noise = rep(1:n_noise, times = length(n_loess)), 
                     n_loess = rep(n_loess, each = n_noise), method = "xi_cor")

Hoeff_d = data.frame(power = c(t(power_Hoeff)), noise = rep(1:n_noise, times = length(n_loess)), 
                     n_loess = rep(n_loess, each = n_noise), method = "Hoeff")

RHoeff_d = data.frame(power = c(t(power_RHoeff)), noise = rep(1:n_noise, times = length(n_loess)), 
                      n_loess = rep(n_loess, each = n_noise), method = "RHoeff")

BDY_d = data.frame(power = c(t(power_BDY)), noise = rep(1:n_noise, times = length(n_loess)), 
                   n_loess = rep(n_loess, each = n_noise), method = "BDY")

Sadp_d = data.frame(power = c(t(power_Sadp)), noise = rep(1:n_noise, times = length(n_loess)), 
                   n_loess = rep(n_loess, each = n_noise), method = "S_ADP")

dHSIC_d = data.frame(power = c(t(power_dHSIC)), noise = rep(1:n_noise, times = length(n_loess)), 
                    n_loess = rep(n_loess, each = n_noise), method = "dHSIC")

SIT_d = data.frame(power = c(t(power_SIT)), noise = rep(1:n_noise, times = length(n_loess)), 
                   n_loess = rep(n_loess, each = n_noise), method = "SIT")


Results = ar_d %>% 
  dplyr::full_join(cor_d) %>% 
  dplyr::full_join(dcor_d) %>% 
  dplyr::full_join(MIC_d) %>%
  dplyr::full_join(rdc_d) %>%
  dplyr::full_join(xicor_d) %>%
  dplyr::full_join(Hoeff_d) %>%
  dplyr::full_join(RHoeff_d) %>%
  dplyr::full_join(BDY_d) %>%
  dplyr::full_join(Sadp_d) %>%
  dplyr::full_join(dHSIC_d) %>%
  dplyr::full_join(SIT_d)

write.table(Results, "power_analysis_loess/power_results_loess.txt", quote = T, row.names = F)

#Results = read.table("power_analysis_loess/power_results_loess.txt", header = T)

