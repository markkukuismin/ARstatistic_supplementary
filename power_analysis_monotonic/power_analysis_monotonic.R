
# 07.03.2024

# Examine the power compared to other measures, such as (i) Pearson correlation,
# (ii) distance correlation (iii) MIC, (iv) Randomized Dependence Coefficient
# (RDC), (v) Spearman, (vi) xi_cor, (vii) Hoeffding, (viii) Refined Hoeffding,
# (ix) Bergsma-Dassios-Yanagimoto coefficient (BDY), (x) S_ADP test, (xi) dHSIC,
# (xii) Sliced Independence Test (SIT)

# Update 18.7.2022: the output probability of the A-R method is scaled by the 
# expected value of the accept ratio which is derived when a strictly monotonic
# dependency between x and y is examined. This is to
# make the interpretation of the AR-dependency statistic easier; if the accept
# proportion/probability is 1, then x and y are strictly monotonically dependent.

# Update 20.9.2022: added xi_cor

# Update 20.10.2022: added rho_dist

# Update: 16.3.2023: rho_dist replaced with expected value

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

## We write the following function because MIC glitches a small percentage of the time, and we do not wish to average over those trials

notNA.greater <- function(a,b){
  ind <- which(!is.na(a))
  pow <- sum(a[ind] > b)/length(ind)
  return(pow)
}

source("functions/rdc.R")
source("functions/rho_rank.R")

n = 320

#noise = 3
noise = 9

num_noise = 30

M = 500

val_null_ar = val_null_cor = val_null_dcor = 
  val_null_mine = val_null_rdc = val_null_xicor = 
  val_null_Hoeff = val_null_RHoeff = val_null_BDY =
  val_null_dHSIC = val_null_SIT = val_null_spear =
  rep(NA, M)  # Vectors holding the null "correlations"

val_ar = val_cor = val_dcor = 
  val_mine = val_rdc = val_xicor = 
  val_Hoeff = val_RHoeff = val_BDY =
  val_Sadp = val_dHSIC = val_SIT = val_spear = 
  rep(NA, M) # Vectors holding the alternative "correlations"
 
D_type = c("exp", "log")

power_ar = power_cor = power_dcor = power_mine = power_rdc = power_xicor = 
  power_Hoeff = power_RHoeff = power_BDY = power_Sadp = power_dHSIC = 
  power_SIT = power_spear = array(NA, c(length(D_type), num_noise)) 

set.seed(1)

# Null table for the S_ADP test. Create null table for aggregation by summation,

null_table_mxl = HHG::Fast.independence.test.nulltable(n, variant = 'ADP-EQP-ML', 
                                                       nr.atoms = 40, nr.perm = 1000)

for(l in 1:num_noise){
  
  j = 1
  
  for(d_type in D_type){
    
    for(i in 1:M){
      
      x = runif(n)
      
      if(d_type == "exp") y = exp(x) + noise*(l/num_noise)*rnorm(n)
      
      if(d_type == "log") y = log(x) + noise*(l/num_noise)*rnorm(n)
      
      # Under null,
      
      x = runif(n)
      
      val_null_ar[i] = rho_rank(x, y)
      val_null_cor[i] = (cor(x,y))^2            # Calculate the correlation
      val_null_dcor[i] = dcor(x,y)              # Calculate dcor
      val_null_mine[i] = minerva::mine(x, y)$MIC # Calculate mic
      val_null_rdc[i] = rdc(x, y)
      val_null_xicor[i] = XICOR::calculateXI(x, y)
      Dn = independence::hoeffding.D.test(x, y, precision = 1)$Dn
      val_null_Hoeff[i] = Dn
      Rn = independence::hoeffding.refined.test(x, y, precision = 1)$Rn
      val_null_RHoeff[i] = Rn
      val_null_BDY[i] = 12*(Dn + 2*Rn) #independence::tau.star.test(x, y, precision = 1)$Tn
      val_null_dHSIC[i] = dHSIC::dhsic(x, y)$dHSIC
      val_null_SIT[i] = SIT::sitcor(x, y, c = 16)
      val_null_spear[i] = (cor(x, y, method = "spearman"))^2
      
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
    cut_spear = quantile(val_null_spear, .95)
    
    # Next we simulate the data again, this time under the alternative
    
    for(i in 1:M){
      
      x = runif(n)
      
      if(d_type == "exp") y = exp(x) + noise*(l/num_noise)*rnorm(n)
      
      if(d_type == "log") y = log(x) + noise*(l/num_noise)*rnorm(n)
      
      val_ar[i] = rho_rank(x, y)
      val_cor[i] = (cor(x,y))^2            # Calculate the correlation
      val_dcor[i] = dcor(x,y)              # Calculate dcor
      val_mine[i] = minerva::mine(x, y)$MIC # Calculate mic
      val_rdc[i] = rdc(x, y)
      val_xicor[i] = XICOR::calculateXI(x, y)
      Dn = independence::hoeffding.D.test(x, y, precision = 1)$Dn
      val_Hoeff[i] = Dn
      Rn = independence::hoeffding.refined.test(x, y, precision = 1)$Rn
      val_RHoeff[i] = Rn
      val_BDY[i] = 12*(Dn + 2*Rn) #independence::tau.star.test(x, y, precision = 1)$Tn
      val_Sadp[i] = HHG::Fast.independence.test(x, y, null_table_mxl, combining.type = "MinP")$MinP.pvalue
      val_dHSIC[i] = dHSIC::dhsic(x, y)$dHSIC
      val_SIT[i] = SIT::sitcor(x, y, c = 16)
      val_spear[i] = (cor(x, y, method = "spearman"))^2
      
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
    power_spear[j, l] = sum(val_spear > cut_spear)/M
    
    j = j + 1
    
  }
  
  cat("\r", l)
  
}

## The rest of the code is for plotting the image


ar_d = data.frame(power = c(t(power_ar)), noise = rep(1:num_noise, times = length(D_type)), 
                  d_type = rep(D_type, each = num_noise), method = "ar")

cor_d = data.frame(power = c(t(power_cor)), noise = rep(1:num_noise, times = length(D_type)), 
                  d_type = rep(D_type, each = num_noise), method = "cor")

dcor_d = data.frame(power = c(t(power_dcor)), noise = rep(1:num_noise, times = length(D_type)), 
                  d_type = rep(D_type, each = num_noise), method = "dcor")

MIC_d = data.frame(power = c(t(power_mine)), noise = rep(1:num_noise, times = length(D_type)), 
                  d_type = rep(D_type, each = num_noise), method = "MIC")

rdc_d = data.frame(power = c(t(power_rdc)), noise = rep(1:num_noise, times = length(D_type)), 
                   d_type = rep(D_type, each = num_noise), method = "rdc")

spear_d = data.frame(power = c(t(power_spear)), noise = rep(1:num_noise, times = length(D_type)), 
                     d_type = rep(D_type, each = num_noise), method = "spearman")

xicor_d = data.frame(power = c(t(power_xicor)), noise = rep(1:num_noise, times = length(D_type)), 
                     d_type = rep(D_type, each = num_noise), method = "xi_cor")

Hoeff_d = data.frame(power = c(t(power_Hoeff)), noise = rep(1:num_noise, times = length(D_type)), 
                     d_type = rep(D_type, each = num_noise), method = "Hoeff")

RHoeff_d = data.frame(power = c(t(power_RHoeff)), noise = rep(1:num_noise, times = length(D_type)), 
                      d_type = rep(D_type, each = num_noise), method = "RHoeff")

BDY_d = data.frame(power = c(t(power_BDY)), noise = rep(1:num_noise, times = length(D_type)), 
                   d_type = rep(D_type, each = num_noise), method = "BDY")

Sadp_d = data.frame(power = c(t(power_Sadp)), noise = rep(1:num_noise, times = length(D_type)), 
                   d_type = rep(D_type, each = num_noise), method = "S_ADP")

dHSIC_d = data.frame(power = c(t(power_dHSIC)), noise = rep(1:num_noise, times = length(D_type)), 
                    d_type = rep(D_type, each = num_noise), method = "dHSIC")

SIT_d = data.frame(power = c(t(power_SIT)), noise = rep(1:num_noise, times = length(D_type)), 
                   d_type = rep(D_type, each = num_noise), method = "SIT")



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
  dplyr::full_join(SIT_d) %>%
  dplyr::full_join(spear_d)

write.table(Results, "power_analysis_monotonic/monotonic_power_results.txt", quote = T, row.names = F)

Results = Results %>% dplyr::rename(Method = method, Power = power, Noise = noise)

p1 = ggplot(data = subset(Results, d_type %in% "exp"), 
            aes(x = Noise, y = Power, group = Method, color = Method)) +
  geom_line(aes(linetype = Method)) +
  geom_point(aes(shape = Method)) +
  scale_linetype_manual(values = c(1:6, 1:6, 1)) + 
  #scale_shape_manual(values = c(16, 17, 15, 3, 7, 8, 10)) +
  scale_shape_manual(values = c(15:25, 3, 7)) +
  labs(title = "exp(x)")

p2 = ggplot(data = subset(Results, d_type %in% "log"), 
            aes(x = Noise, y = Power, group = Method, color = Method)) +
  geom_line(aes(linetype = Method)) +
  geom_point(aes(shape = Method)) +
  scale_linetype_manual(values = c(1:6, 1:6, 1)) +
  #scale_shape_manual(values = c(16, 17, 15, 3, 7, 8, 10)) +
  scale_shape_manual(values = c(15:25, 3, 7)) +
  labs(title = "log(x)")

power_plot = grid.arrange(p1, p2, ncol = 2, nrow = 1)

# Make a nicer plot,

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

p1 = p1 + theme(legend.position = "bottom") + ggplot2::guides(linetype = guide_legend(nrow = 2))

p1

legend = get_legend(p1)

p1 = p1 + theme(legend.position="none") + ylab(element_blank())
p2 = p2 + theme(legend.position="none") + ylab(element_blank())

power_plot = grid.arrange(arrangeGrob(p1, p2, legend, ncol = 2, nrow = 2, 
                                      layout_matrix = matrix(c(1, 2, 
                                                               3, 3), ncol = 2, byrow = T),
                                      widths = c(2.7, 2.7), heights = c(2.5, 0.2),
                                      left = textGrob("Power", rot = 90, vjust = 1))
                          )

ggsave(filename="power_analysis_monotonic/monotonic_ggplot_power.pdf", 
       plot = power_plot, 
       device = cairo_pdf, 
       width = 400, 
       height = 297, 
       units = "mm")

ggsave(filename="power_analysis_monotonic/monotonic_ggplot_power.png", 
       plot = power_plot, 
       width = 270, 
       height = 150, 
       units = "mm")
