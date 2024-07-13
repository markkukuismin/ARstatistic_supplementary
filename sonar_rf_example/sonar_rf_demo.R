
# 6.4.2024

library(mlbench)
library(caret)
library(energy)
library(XICOR)
library(minerva)
library(SIT)
library(tidyverse)
library(doParallel)

source("functions/rho_rank.R")
source("functions/rdc.R")
source("functions/findCorrelation_length.R")

set.seed(1)

data("Sonar")

X = Sonar[, -ncol(Sonar)]
y = Sonar[, ncol(Sonar)]

# First, use all covariates

cl = makePSOCKcluster(5)
registerDoParallel(cl)

control = trainControl(method = "repeatedcv", 
                       number = 10, 
                       repeats = 50)

rf_all = caret::train(X, 
                      y, 
                      method = "rf",
                      metric = "Accuracy",
                      ntree = 1000,
                      trControl = control)

rf_all

##

# Feature selection: does it increase the accuracy?

# Pearson correlation

C_cor = cor(X)

# Other association statistics

p = ncol(X)

C_ar = C_dcor = C_xicor = C_rdc = C_sit = C_mic = matrix(0, p, p)

C_ar[upper.tri(C_ar)] = 1

ind = which(C_ar == 1, arr.ind = T)

for(k in 1:nrow(ind)){
  
  C_ar[ind[k, 1], ind[k, 2]] = rho_rank(X[, ind[k, 1]],
                                        X[, ind[k, 2]])
  
  C_dcor[ind[k, 1], ind[k, 2]] = dcor(X[, ind[k, 1]],
                                      X[, ind[k, 2]])
  
  C_xicor[ind[k, 1], ind[k, 2]] = max(calculateXI(X[, ind[k, 1]],
                                                  X[, ind[k, 2]]),
                                      calculateXI(X[, ind[k, 2]],
                                                  X[, ind[k, 1]])
                                      )
  
  C_rdc[ind[k, 1], ind[k, 2]] = rdc(X[, ind[k, 1]],
                                    X[, ind[k, 2]])
  
  C_sit[ind[k, 1], ind[k, 2]] = max(sitcor(X[, ind[k, 1]],
                                           X[, ind[k, 2]], c = 16),
                                    sitcor(X[, ind[k, 2]],
                                           X[, ind[k, 1]], c = 16)
                                    )
  
  C_mic[ind[k, 1], ind[k, 2]] = mine(X[, ind[k, 1]],
                                     X[, ind[k, 2]])$MIC
  
  cat("\r", k)
  
}

C_ar = C_ar + t(C_ar)
diag(C_ar) = 1

C_dcor = C_dcor + t(C_dcor)
diag(C_dcor) = 1

C_xicor = C_xicor + t(C_xicor)
diag(C_xicor) = 1

C_rdc = C_rdc + t(C_rdc)
diag(C_rdc) = 1

C_sit = C_sit + t(C_sit)
diag(C_sit) = 1

C_mic = C_mic + t(C_mic)
diag(C_mic) = 1

##

# Find cutoff values which return the same number of features,

f_cor = findCorrelation_length(C_cor)
f_ar = findCorrelation_length(C_ar)
f_dcor = findCorrelation_length(C_dcor)
f_xicor = findCorrelation_length(C_xicor)
f_rdc = findCorrelation_length(C_rdc)
f_sit = findCorrelation_length(C_sit)
f_mic = findCorrelation_length(C_mic)

# Select the number of features and the corresponding cutoffs,

ab = Reduce(intersect, list(f_cor$feature_nmb, f_ar$feature_nmb,
                            f_dcor$feature_nmb, f_xicor$feature_nmb,
                            f_rdc$feature_nmb, f_sit$feature_nmb,
                            f_mic$feature_nmb))

ab

nmb_ind = c(52, 40, 30, 20, 10, 5)

cor_thr = f_cor$feature_nmb %in% nmb_ind
cor_thr = f_cor$thr[cor_thr]

ar_thr = f_ar$feature_nmb %in% nmb_ind
ar_thr = f_ar$thr[ar_thr]

dcor_thr = f_dcor$feature_nmb %in% nmb_ind
dcor_thr = f_dcor$thr[dcor_thr]

xicor_thr = f_xicor$feature_nmb %in% nmb_ind
xicor_thr = f_xicor$thr[xicor_thr]

rdc_thr = f_rdc$feature_nmb %in% nmb_ind
rdc_thr = f_rdc$thr[rdc_thr]

sit_thr = f_sit$feature_nmb %in% nmb_ind
sit_thr = f_sit$thr[sit_thr]

mic_thr = f_mic$feature_nmb %in% nmb_ind
mic_thr = f_mic$thr[mic_thr]

thrs = data.frame(cor = cor_thr, ar = ar_thr, dcor = dcor_thr, 
                  xicor = xicor_thr, rdc = rdc_thr, SIT = sit_thr,
                  MIC = mic_thr)

C = list(cor = C_cor, ar = C_ar, dcor = C_dcor, 
         xicor = C_xicor, rdc = C_rdc, SIT = C_sit, 
         MIC = C_mic)

acc = Method = c()

for(j in 1:length(C)){
  
  for(i in 1:length(nmb_ind)){
    
    thr_temp = thrs[, names(C)[j]]
    
    temp_features = findCorrelation(C[[j]], cutoff = thr_temp[i])
    
    rf_temp = caret::train(X[, temp_features, drop = FALSE], 
                           y, 
                           method = "rf",
                           metric = "Accuracy",
                           ntree = 1000,
                           trControl = control)
    
    acc = c(acc, max(rf_temp$results["Accuracy"]))
    
    Method = c(Method, names(C)[j])
    
  }
  
  cat("\r", j)
  
}

stopCluster(cl)

df = data.frame(acc = acc,
                Method = as.factor(Method),
                nmb = as.factor(nmb_ind))

write.table(df, "sonar_rf_example/sonar_rf_results.txt", 
            row.names = FALSE)

p_acc = ggplot(df, aes(nmb, acc, colour = Method, shape = Method)) + 
  geom_point(position = position_dodge(0.5), size = 5) +
  xlab("Number of covariates") +
  ylab("Training accuracy") +
  geom_hline(yintercept = max(rf_all$results["Accuracy"]),
             linetype = "dashed", linewidth = 1.5) +
  theme(text = element_text(size = 20)) +
  scale_shape_manual(values = c(15:21))

p_acc

ggsave(filename="sonar_rf_example/train_accuracy.pdf", 
       plot = p_acc, 
       device = cairo_pdf, 
       width = 400, 
       height = 297, 
       units = "mm")

ggsave(filename="sonar_rf_example/train_accuracy.png", 
       plot = p_acc, 
       device = png, 
       width = 400, 
       height = 297, 
       units = "mm",
       dpi = 700)

# Almost all statistics increase the training accuracy of RF

df %>% group_by(Method) %>% 
  summarise(m = max(acc)) %>%
  arrange(desc(m))

rf_all

# Illustrate the association between different patterns,

#idcor = which(C_dcor == max(C_dcor[upper.tri(C_dcor)]), 
#              arr.ind = TRUE)

#iar = which(C_ar == max(C_ar[upper.tri(C_ar)]), 
#            arr.ind = TRUE)

idcor = matrix(c(18, 17, 17, 18), 2, 2)
iar = matrix(c(16, 15, 15, 16), 2, 2)

x = c(X[, idcor[1, 1]], X[, idcor[1, 2]])
y = c(X[, iar[1, 1]], X[, iar[1, 2]])

n = nrow(X)

Method = rep(c("dcor", "ar"), each = n)

df = data.frame(x = x, y = y, Method = Method)

p1 = ggplot(subset(df, Method == "dcor"), aes(x, y)) +
  geom_point() +
  xlab(paste0("Spectral ", idcor[1, 1])) +
  ylab(paste0("Spectral ", idcor[1, 2])) +
  geom_smooth() +
  geom_text(x = 0.2, y = 1, 
            label = round(max(C_dcor[upper.tri(C_dcor)]), 3),
            size = 9) +
  ggtitle("(A)") +
  theme(plot.title = element_text(size = 20, face = "bold"),
        text = element_text(size = 20))

p2 = ggplot(subset(df, Method == "ar"), aes(x, y)) +
  geom_point() +
  xlab(paste0("Spectral ", iar[1, 1])) +
  ylab(paste0("Spectral ", iar[1, 2])) +
  geom_smooth() +
  geom_text(x = 0.2, y = 1, 
            label = round(max(C_ar[upper.tri(C_ar)]), 3),
            size = 9) +
  ggtitle("(B)") +
  theme(plot.title = element_text(size = 20, face = "bold"),
        text = element_text(size = 20))

p_final = gridExtra::grid.arrange(p1, p2, ncol = 2)

ggsave(filename="sonar_rf_example/spectral_scatterplot.pdf", 
       plot = p_final, 
       device = cairo_pdf, 
       width = 400, 
       height = 297, 
       units = "mm")

ggsave(filename="sonar_rf_example/spectral_scatterplot.png", 
       plot = p_final, 
       device = png, 
       width = 400, 
       height = 297, 
       units = "mm",
       dpi = 700)
