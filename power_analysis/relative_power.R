
# 7.3.2024

library(dplyr)
library(pracma)
library(ggplot2)

Results = read.table("power_analysis/power_results.txt", header = T)

Results_monotonic = read.table("power_analysis_monotonic/monotonic_power_results.txt", header = T)

Results = dplyr::bind_rows(Results, Results_monotonic)

Results = Results %>% dplyr::filter(method != "spearman")

rm(Results_monotonic)

Results = Results %>% 
  mutate(
    d_type = case_match(d_type, 
                        "clusters" ~ "5 clusters",
                        .default = d_type)
    )

d_type = unique(Results$d_type)

method = unique(Results$method)

x = unique(Results$noise)

D = Results %>% group_by(d_type, method) %>% summarise(AUC = pracma::trapz(x, power))

D = D %>% arrange(method)

D$AUC = D$AUC/30 # 30 is the theoretical max AUC. Thus AUC in D is a proportion of AUC of the max AUC

D_cor = D %>% filter(method == "cor")

D = D %>% group_by(method) %>% mutate(relative_AUC = AUC/D_cor$AUC) 

ggplot(data = D, mapping = aes(x = 1, y = relative_AUC, group = method)) +
  geom_point(aes(color = method, shape = method), size = 3) +
  scale_shape_manual(values = c(15:25, 3)) +
  ggforce::facet_row(vars(d_type), scales = "free", space = "free") + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  geom_hline(yintercept = 1) +
  ylab("AUC relative to cor")



D = D %>% group_by(method) %>% mutate(sum_AUC = sum(AUC), mean_AUC = mean(AUC), sd_AUC = sd(AUC))

D2 = D %>% group_by(method) %>% summarise(sum_AUC = unique(sum_AUC), mean_AUC = unique(mean_AUC), sd_AUC = unique(sd_AUC))

D2 %>% arrange(desc(sum_AUC))

# The percentage of proportional AUC with respect to the method returning the max prop. AUC 

D2 %>% arrange(desc(sum_AUC)) %>% mutate(relative_AUC = sum_AUC/max(sum_AUC), CV = sd_AUC/mean_AUC)

AUC_boxplot = ggplot(data = D, mapping = aes(x = reorder(method, AUC, FUN = median, decreasing = TRUE), 
                                             y = AUC, group = method)) +
  geom_boxplot() +
  labs(x = "Method",
       y = "AUC proportional to the max AUC") +
  theme(text = element_text(size = 25))  

AUC_boxplot

ggsave(filename="power_analysis/AUC_boxplot.pdf", 
       plot = AUC_boxplot, 
       device = cairo_pdf, 
       width = 400, 
       height = 297, 
       units = "mm")

ggsave(filename="power_analysis/AUC_boxplot.png", 
       plot = AUC_boxplot, 
       device = png, 
       width = 420, 
       height = 297, 
       units = "mm",
       dpi = 700)

# How much lower is the power of the A-R method in the linear case compared to 
# Pearson and distance correlation?

100*(1 - D$AUC[D$method == "ar" & D$d_type == "linear"]/D$AUC[D$method == "cor" & D$d_type == "linear"])
100*(1 - D$AUC[D$method == "dcor" & D$d_type == "linear"]/D$AUC[D$method == "cor" & D$d_type == "linear"])
100*(1 - D$AUC[D$method == "S_ADP" & D$d_type == "linear"]/D$AUC[D$method == "cor" & D$d_type == "linear"])


100*(1 - D$AUC[D$method == "ar" & D$d_type == "linear"]/D$AUC[D$method == "dcor" & D$d_type == "linear"])

# How much lower is the power of the A-R method in the quadratic case compared to 
# Pearson and distance correlation?

D$AUC[D$method == "ar" & D$d_type == "quadratic"]/D$AUC[D$method == "cor" & D$d_type == "quadratic"]
D$AUC[D$method == "dcor" & D$d_type == "quadratic"]/D$AUC[D$method == "cor" & D$d_type == "quadratic"]
D$AUC[D$method == "S_ADP" & D$d_type == "quadratic"]/D$AUC[D$method == "cor" & D$d_type == "quadratic"]

D$AUC[D$method == "ar" & D$d_type == "quadratic"]/D$AUC[D$method == "dcor" & D$d_type == "quadratic"]

#

ar_vs_cor = D$AUC[D$method == "ar"]/D$AUC[D$method == "cor"]
names(ar_vs_cor) = D$d_type[D$method == "ar"]
ar_vs_cor
         
SADP_vs_cor = D$AUC[D$method == "S_ADP"]/D$AUC[D$method == "cor"]
names(SADP_vs_cor) = D$d_type[D$method == "S_ADP"]
SADP_vs_cor

dcor_vs_cor = D$AUC[D$method == "dcor"]/D$AUC[D$method == "cor"]
names(dcor_vs_cor) = D$d_type[D$method == "dcor"]
dcor_vs_cor

ar_vs_dcor = D$AUC[D$method == "ar"]/D$AUC[D$method == "dcor"]
names(ar_vs_dcor) = D$d_type[D$method == "ar"]
ar_vs_dcor   

ar_vs_SADP = D$AUC[D$method == "ar"]/D$AUC[D$method == "S_ADP"]
names(ar_vs_SADP) = D$d_type[D$method == "ar"]
ar_vs_SADP   

##

Results = read.table("power_analysis/power_results.txt", header = T)

Results_monotonic = read.table("power_analysis_monotonic/monotonic_power_results.txt", header = T)

Results = dplyr::bind_rows(Results, Results_monotonic)

rm(Results_monotonic)

D = Results %>% group_by(d_type, method) %>% summarise(AUC = pracma::trapz(x, power))

D$AUC = D$AUC/30

D %>% 
  filter(d_type == "linear") %>%
  select(method, AUC) %>%
  arrange(desc(AUC))

D %>% 
  filter(d_type == "exp") %>%
  select(method, AUC) %>%
  arrange(desc(AUC))

D %>% 
  filter(d_type == "log") %>%
  select(method, AUC) %>%
  arrange(desc(AUC))


# Do all methods always find the association when noise is small?

D_always_find = Results %>% group_by(d_type, method) %>% dplyr::filter(power == 1 & noise <= 1)

D_always_find %>% 
  group_by(method, noise) %>% 
  summarise(n()) %>% 
  dplyr::arrange(dplyr::desc(`n()`)) %>%
  print(n = 30)
