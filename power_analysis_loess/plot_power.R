
# 26.4.2023

library(dplyr)
library(pracma)
library(ggplot2)

Results = read.table("power_analysis_loess/power_results_loess.txt", header = T)

D = Results %>% group_by(noise, method) %>% summarise(mean_power = mean(power))

colnames(D) = c("Noise", "Method", "mean_power")

mean_power_lineplot = ggplot(data = D, 
                             aes(x = Noise, y = mean_power, color = Method)) +
  geom_line(aes(linetype = Method), linewidth = 1.5) +
  geom_point(aes(shape = Method), size = 4.5) +
  #scale_shape_manual(values = c(16, 17, 15, 3, 7, 8)) +
  scale_shape_manual(values = c(15:25, 3)) +
  #scale_color_manual(values = viridis::mako(12)) +
  ylab("Mean power")

mean_power_lineplot = mean_power_lineplot + 
  theme(legend.position = "bottom") + 
  ggplot2::guides(shape = guide_legend(nrow = 2)) +
  scale_x_continuous(breaks = 1:10,
                     labels = as.character(1:10)) +
  theme(text = element_text(size = 20)) 

mean_power_lineplot

ggsave(filename="power_analysis_loess/mean_power_lineplot.pdf", 
       plot = mean_power_lineplot, 
       device = cairo_pdf, 
       width = 400, 
       height = 297, 
       units = "mm")

ggsave(filename="power_analysis_loess/mean_power_lineplot.png", 
       plot = mean_power_lineplot, 
       device = png, 
       width = 250, 
       height = 120, 
       units = "mm",
       dpi = 700)
