
# 21.4.2023

library(dplyr)
library(pracma)
library(ggplot2)

Results = read.table("power_analysis_loess/power_results_loess.txt", header = T)

n_loess = unique(Results$n_loess)

method = unique(Results$method)

x = unique(Results$noise)

D = Results %>% group_by(n_loess, method) %>% summarise(AUC = pracma::trapz(x, power))

D = D %>% arrange(method)

D$AUC = D$AUC/10 # 10 is the theoretical max AUC. Thus AUC in D is a proportion of AUC of the max AUC

D_cor = D %>% filter(method == "cor")

D = D %>% group_by(method) %>% mutate(relative_AUC = AUC/D_cor$AUC) 

D = D %>% group_by(method) %>% mutate(sum_AUC = sum(AUC), mean_AUC = mean(AUC), sd_AUC = sd(AUC))

D2 = D %>% group_by(method) %>% summarise(sum_AUC = unique(sum_AUC), mean_AUC = unique(mean_AUC), sd_AUC = unique(sd_AUC))

D2 %>% arrange(desc(sum_AUC))

# The percentage of proportional AUC with respect to the method returning the max prop. AUC 

D2 %>% arrange(desc(sum_AUC)) %>% mutate(relative_AUC = sum_AUC/max(sum_AUC), CV = sd_AUC/mean_AUC)

AUC_boxplot = ggplot(data = D, 
                     mapping = aes(x = reorder(method, AUC, FUN = median, decreasing = TRUE), 
                                   y = AUC, group = method)) +
  geom_boxplot() +
  labs(x = "Method",
       y = "AUC proportional to the max AUC") +
  theme(text = element_text(size = 25))  

AUC_boxplot

ggsave(filename="power_analysis_loess/AUC_boxplot_loess.pdf", 
       plot = AUC_boxplot, 
       device = cairo_pdf, 
       width = 420, 
       height = 297, 
       units = "mm")

ggsave(filename="power_analysis_loess/AUC_boxplot_loess.png", 
       plot = AUC_boxplot, 
       device = png, 
       width = 420, 
       height = 297, 
       units = "mm",
       dpi = 700)
