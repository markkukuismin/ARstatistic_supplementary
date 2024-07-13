
# 06.03.2024

# Examine the power compared to other measures, such as (i) Pearson correlation,
# (ii) distance correlation (iii) MIC, (iv) Randomized Dependence Coefficient
# (RDC), (v) xi_cor, (vi) Hoeffding, (vii) Refined Hoeffding,
# (viii) Bergsma-Dassios-Yanagimoto coefficient (BDY), (ix) S_ADP test, (x) dHSIC,
# (xi) Sliced Independence Test (SIT)

# Update 18.7.2022: the output probability of the A-R method is scaled by a 
# probability which is derived when a strictly  monotonic dependency between x 
# and y is examined. This is to make the interpretation of the AR-dependency 
# probability easier; if the prob. is 1, then x and y are strictly monotonically
# dependent.

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
#source("functions/rho.R")
source("functions/rho_rank.R")
source("functions/clusters.R")

# We write the following function because MIC glitches a small percentage of 
# the time, and we do not wish to average over those trials

notNA.greater <- function(a,b){
  ind <- which(!is.na(a))
  pow <- sum(a[ind] > b)/length(ind)
  return(pow)
}

n = 25

noise = 3

num_noise = 30

M = 500

val_null_ar = val_null_cor = val_null_dcor = 
  val_null_mine = val_null_rdc = val_null_xicor = 
  val_null_Hoeff = val_null_RHoeff = val_null_BDY = 
  val_null_dHSIC = val_null_SIT = rep(NA, M)  # Vectors holding the null "correlations"

val_ar = val_cor = val_dcor = 
  val_mine = val_rdc = val_xicor = 
  val_Hoeff = val_RHoeff = val_BDY =
  val_Sadp = val_dHSIC = val_SIT = rep(NA, M) # Vectors holding the alternative "correlations"

D_type = c("linear", "quadratic", "cubic", "sine period 1/2", "sine period 1/8",
           "x^(1/4)", "circle", "step", "cross", "clusters", "two parabola", "spiral")

power_ar = power_cor = power_dcor = power_mine = power_rdc = power_xicor = 
  power_Hoeff = power_RHoeff = power_BDY = power_Sadp = power_dHSIC = 
  power_SIT = array(NA, c(length(D_type), num_noise)) 

set.seed(1)

# Null table for the S_ADP test. Create null table for aggregation by summation,

null_table_mxl = HHG::Fast.independence.test.nulltable(n, variant = 'ADP-EQP-ML', 
                                                       nr.atoms = n, nr.perm = 1000)

for(l in 1:num_noise){
  
  j = 1
  
  for(d_type in D_type){
    
    for(i in 1:M){
      
      x = runif(n)
      
      if(d_type == "sine period 1/2") y = sin(4*pi*x) + 2*noise * (l/num_noise)*rnorm(n)
      
      if(d_type == "sine period 1/8") y = sin(16*pi*x) + noise*(l/num_noise)*rnorm(n)
      
      if(d_type == "linear") y = x + noise*(l/num_noise)*rnorm(n) 
      
      if(d_type == "quadratic") y = 4*(x - 0.5)^2 + noise*(l/num_noise)*rnorm(n)
      
      if(d_type == "cubic") y = 128*(x - 1/3)^3 -48*(x - 1/3)^3 - 12*(x - 1/3) + 10*noise*(l/num_noise)*rnorm(n)
      
      if(d_type == "circle") y = (2*rbinom(n, 1, 0.5) - 1)*(sqrt(1 - (2*x - 1)^2)) + noise/4*l/num_noise*rnorm(n)
      
      if(d_type == "step") y = (x > 0.5) + noise*5*l/num_noise*rnorm(n)
      
      if(d_type == "x^(1/4)") y = x^(1/4) + noise*(l/num_noise)*rnorm(n)
      
      if(d_type == "cross") y = (2*rbinom(n, 1, 0.5) - 1)*(x - 0.5) + noise/4*(l/num_noise)*rnorm(n)
      
      if(d_type == "clusters"){
        yx = clusters(n, s = noise/3*(l/num_noise)) 
        x = yx$x
        y = yx$y
      }
      
      if(d_type == "two parabola") y = (4*(x - 0.5)^2 + noise/2*(l/num_noise)*rnorm(n))*(2*rbinom(n, 1, 0.5) - 1)
      
      if(d_type == "spiral"){
        theta = runif(n)
        x = 0.5*5*pi*theta*cos(5*pi*theta)
        y = 0.5*5*pi*theta*sin(5*pi*theta) + 2*noise*(l/num_noise)*rnorm(n)
      }
      
      # Under null,
      
      x = runif(n)
      
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
      
      x = runif(n)
      
      if(d_type == "sine period 1/2") y = sin(4*pi*x) + 2*noise * (l/num_noise)*rnorm(n)
      
      if(d_type == "sine period 1/8") y = sin(16*pi*x) + noise*(l/num_noise)*rnorm(n)
      
      if(d_type == "linear") y = x + noise*(l/num_noise)*rnorm(n) 
      
      if(d_type == "quadratic") y = 4*(x - 0.5)^2 + noise*(l/num_noise)*rnorm(n)
      
      if(d_type == "cubic") y = 128*(x - 1/3)^3 -48*(x - 1/3)^3 - 12*(x - 1/3) + 10*noise*(l/num_noise)*rnorm(n)
      
      if(d_type == "circle") y = (2*rbinom(n, 1, 0.5) - 1)*(sqrt(1 - (2*x - 1)^2)) + noise/4*l/num_noise*rnorm(n)
      
      if(d_type == "step") y = (x > 0.5) + noise*5*l/num_noise*rnorm(n)
      
      if(d_type == "x^(1/4)") y = x^(1/4) + noise*(l/num_noise)*rnorm(n)
      
      if(d_type == "cross") y = (2*rbinom(n, 1, 0.5) - 1)*(x - 0.5) + noise/4*(l/num_noise)*rnorm(n)
      
      if(d_type == "clusters"){
        yx = clusters(n, s = noise/3*(l/num_noise)) 
        x = yx$x
        y = yx$y
      }
      
      if(d_type == "two parabola") y = (4*(x - 0.5)^2 + noise/2*(l/num_noise)*rnorm(n))*(2*rbinom(n, 1, 0.5) - 1)
      
      if(d_type == "spiral"){
        theta = runif(n)
        x = 0.5*5*pi*theta*cos(5*pi*theta)
        y = 0.5*5*pi*theta*sin(5*pi*theta) + 2*noise*(l/num_noise)*rnorm(n)
      }
      
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
  dplyr::full_join(SIT_d)

write.table(Results, "power_analysis_small_sample/power_results.txt", quote = T, row.names = F)

#Results = read.table("power_analysis_small_sample/power_results.txt", header = T)

Results = Results %>% dplyr::rename(Method = method, Power = power, Noise = noise)

Results$Method = as.factor(Results$Method)

p1 = ggplot(data = subset(Results, d_type %in% "linear"), 
            aes(x = Noise, y = Power, color = Method)) +
  geom_line(aes(linetype = Method)) +
  geom_point(aes(shape = Method)) +
  scale_shape_manual(values = c(15:25, 3)) +
  labs(title = "Linear") +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks = element_blank())

p2 = ggplot(data = subset(Results, d_type %in% "quadratic"), 
            aes(x = Noise, y = Power, color = Method)) +
  geom_line(aes(linetype = Method)) +
  geom_point(aes(shape = Method)) +
  scale_shape_manual(values = c(15:25, 3)) +
  labs(title = "Quadratic")  +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks = element_blank())

p3 = ggplot(data = subset(Results, d_type %in% "cubic"), 
            aes(x = Noise, y = Power, color = Method)) +
  geom_line(aes(linetype = Method)) +
  geom_point(aes(shape = Method)) +
  scale_shape_manual(values = c(15:25, 3)) +
  labs(title = "Cubic")  +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks = element_blank())

p4 = ggplot(data = subset(Results, d_type %in% "sine period 1/8"), 
            aes(x = Noise, y = Power, color = Method)) +
  geom_line(aes(linetype = Method)) +
  geom_point(aes(shape = Method)) +
  scale_shape_manual(values = c(15:25, 3)) +
  labs(title = "Sine period 1/8")  +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks = element_blank())

p5 = ggplot(data = subset(Results, d_type %in% "sine period 1/2"), 
            aes(x = Noise, y = Power, color = Method)) +
  geom_line(aes(linetype = Method)) +
  geom_point(aes(shape = Method)) +
  scale_shape_manual(values = c(15:25, 3)) +
  labs(title = "Sine period 1/2")  +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks = element_blank())

p6 = ggplot(data = subset(Results, d_type %in% "x^(1/4)"), 
            aes(x = Noise, y = Power, color = Method)) +
  geom_line(aes(linetype = Method)) +
  geom_point(aes(shape = Method)) +
  scale_shape_manual(values = c(15:25, 3)) +
  #labs(title = "x^(1/4)")  +
  labs(title = expression(x^{1/4})) +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks = element_blank())

p7 = ggplot(data = subset(Results, d_type %in% "circle"), 
            aes(x = Noise, y = Power, color = Method)) +
  geom_line(aes(linetype = Method)) +
  geom_point(aes(shape = Method)) +
  scale_shape_manual(values = c(15:25, 3)) +
  labs(title = "Circle")

p8 = ggplot(data = subset(Results, d_type %in% "step"), 
            aes(x = Noise, y = Power, color = Method)) +
  geom_line(aes(linetype = Method)) +
  geom_point(aes(shape = Method)) +
  scale_shape_manual(values = c(15:25, 3)) +
  labs(title = "Step")

p9 = ggplot(data = subset(Results, d_type %in% "cross"), 
            aes(x = Noise, y = Power, color = Method)) +
  geom_line(aes(linetype = Method)) +
  geom_point(aes(shape = Method)) +
  scale_shape_manual(values = c(15:25, 3)) +
  labs(title = "Cross")

p10 = ggplot(data = subset(Results, d_type %in% "clusters"), 
             aes(x = Noise, y = Power, color = Method)) +
  geom_line(aes(linetype = Method)) +
  geom_point(aes(shape = Method)) +
  scale_shape_manual(values = c(15:25, 3)) +
  labs(title = "5 clusters")

p11 = ggplot(data = subset(Results, d_type %in% "two parabola"), 
             aes(x = Noise, y = Power, color = Method)) +
  geom_line(aes(linetype = Method)) +
  geom_point(aes(shape = Method)) +
  scale_shape_manual(values = c(15:25, 3)) +
  labs(title = "Two parabola")

p12 = ggplot(data = subset(Results, d_type %in% "spiral"), 
             aes(x = Noise, y = Power, color = Method)) +
  geom_line(aes(linetype = Method)) +
  geom_point(aes(shape = Method)) +
  scale_shape_manual(values = c(15:25, 3)) +
  labs(title = "Spiral")

power_plot = grid.arrange(p1, p2, p3, 
                          p4, p5, p6, 
                          p7, p8, p9, 
                          p10, p11, p12,
                          ncol = 3, nrow = 4)
# Make a nicer plot,

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

p1 = p1 + theme(legend.position = "bottom") + ggplot2::guides(shape = guide_legend(nrow = 1))

p1

legend = get_legend(p1)

p1 = p1 + theme(legend.position="none") + ylab(element_blank())
p2 = p2 + theme(legend.position="none") + ylab(element_blank())
p3 = p3 + theme(legend.position="none") + ylab(element_blank())
p4 = p4 + theme(legend.position="none") + ylab(element_blank())
p5 = p5 + theme(legend.position="none") + ylab(element_blank())
p6 = p6 + theme(legend.position="none") + ylab(element_blank())
p7 = p7 + theme(legend.position="none") + ylab(element_blank())
p8 = p8 + theme(legend.position="none") + ylab(element_blank())
p9 = p9 + theme(legend.position="none") + ylab(element_blank())
p10 = p10 + theme(legend.position="none") + ylab(element_blank())
p11 = p11 + theme(legend.position="none") + ylab(element_blank())
p12 = p12 + theme(legend.position="none") + ylab(element_blank())

power_plot = grid.arrange(arrangeGrob(p1, p2, p3, 
                                      p4, p5, p6, 
                                      p7, p8, p9, 
                                      p10, p11, p12,
                                      legend, 
                                      ncol = 3, 
                                      nrow = 5, 
                                      layout_matrix = matrix(c(1, 2, 3,
                                                               4, 5, 6, 
                                                               7, 8, 9, 
                                                               10, 11, 12,
                                                               13, 13, 13), ncol = 3, byrow = T),
                                      widths = c(2.7, 2.7, 2.7), heights = c(2.5, 2.5, 2.5, 2.5, 0.2),
                                      left = textGrob("Power", rot = 90, vjust = 1))
)

ggsave(filename="power_analysis_small_sample/ggplot_power.pdf", 
       plot = power_plot, 
       device = cairo_pdf, 
       width = 400, 
       height = 297, 
       units = "mm")

ggsave(filename="power_analysis_small_sample/ggplot_power.png", 
       plot = power_plot, 
       device = png, 
       width = 400, 
       height = 297, 
       units = "mm",
       dpi = 700)
