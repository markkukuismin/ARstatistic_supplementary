
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

control = trainControl(method = "repeatedcv", 
                       number = 10, 
                       repeats = 50)

##

p = ncol(X)

C_ar = C_dcor = C_xicor = C_rdc = C_sit = C_mic = matrix(0, p, p)

C_ar[upper.tri(C_ar)] = 1

ind = which(C_ar == 1, arr.ind = T)

train_i = createDataPartition(y, 
                              p = 0.7, 
                              list = FALSE)

x_train_temp = X[train_i, ]

for(k in 1:nrow(ind)){
  
  C_ar[ind[k, 1], ind[k, 2]] = rho_rank(x_train_temp[, ind[k, 1]],
                                        x_train_temp[, ind[k, 2]])
  
  C_dcor[ind[k, 1], ind[k, 2]] = dcor(x_train_temp[, ind[k, 1]],
                                      x_train_temp[, ind[k, 2]])
  
  C_xicor[ind[k, 1], ind[k, 2]] = max(calculateXI(x_train_temp[, ind[k, 1]],
                                                  x_train_temp[, ind[k, 2]]),
                                      calculateXI(x_train_temp[, ind[k, 2]],
                                                  x_train_temp[, ind[k, 1]])
                                      )
  
  C_rdc[ind[k, 1], ind[k, 2]] = rdc(x_train_temp[, ind[k, 1]],
                                    x_train_temp[, ind[k, 2]])
  
  C_sit[ind[k, 1], ind[k, 2]] = max(sitcor(x_train_temp[, ind[k, 1]],
                                           x_train_temp[, ind[k, 2]], c = 16),
                                    sitcor(x_train_temp[, ind[k, 2]],
                                           x_train_temp[, ind[k, 1]], c = 16)
                                    )
  
  C_mic[ind[k, 1], ind[k, 2]] = mine(x_train_temp[, ind[k, 1]],
                                     x_train_temp[, ind[k, 2]])$MIC
  
  cat("\r", k)
  
}

C_cor = cor(x_train_temp)

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

# Select the cutoff which returns approximately 20 predictors,

Reduce(intersect, list(f_cor$feature_nmb, f_ar$feature_nmb,
                       f_dcor$feature_nmb, f_xicor$feature_nmb,
                       f_rdc$feature_nmb, f_sit$feature_nmb,
                       f_mic$feature_nmb))

nmb_ind = 20

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

M = 50

C = list(cor = C_cor, ar = C_ar, dcor = C_dcor,
         xicor = C_xicor, rdc = C_rdc, 
         SIT = C_sit, MIC = C_mic)

thrs = list(cor = cor_thr, ar = ar_thr, dcor = dcor_thr,
            xicor = xicor_thr, rdc = rdc_thr, SIT = sit_thr, 
            MIC = mic_thr)

train_test = createDataPartition(y, times = M, p = 0.7, list = FALSE)

Method = c()

acc_all_train = c()
acc_all_test = acc_all_test_low = acc_all_test_up = c()

acc_train = c()
acc_test = acc_test_low = acc_test_up = c()

cl = makePSOCKcluster(5)
registerDoParallel(cl)

for(j in 1:length(C)){
  
  for(i in 1:M){
    
    x_train = X[train_test[, i], ]
    y_train = y[train_test[, i]]
    
    x_test = X[-train_test[, i], ]
    y_test = y[-train_test[, i]]
    
    # All predictors,
    
    if(j == 1){
      
      rf_temp = caret::train(x_train, 
                             y_train, 
                             method = "rf",
                             trControl = control,
                             ntree = 1000)
      
      acc_all_train = c(acc_all_train, 
                        max(rf_temp$results["Accuracy"]))
      
      temp_predictions = predict(rf_temp, x_test)
      
      d_temp = confusionMatrix(temp_predictions, y_test)
      
      acc_all_test = c(acc_all_test, d_temp$overall["Accuracy"])
      acc_all_test_low = c(acc_all_test_low, d_temp$overall["AccuracyLower"])
      acc_all_test_up = c(acc_all_test_up, d_temp$overall["AccuracyUpper"])
      
    }
    
    # Subset of features,
    
    thr_temp = thrs[[names(C)[[j]]]]
    
    temp_features = findCorrelation(C[[j]], cutoff = thr_temp)
    
    rf_temp = caret::train(x_train[, temp_features, drop = FALSE], 
                           y_train, 
                           method = "rf",
                           trControl = control,
                           ntree = 1000)
    
    acc_train = c(acc_train, max(rf_temp$results["Accuracy"]))
    
    temp_predictions = predict(rf_temp, x_test[, temp_features, drop = FALSE])
    
    d_temp = confusionMatrix(temp_predictions, y_test)
    
    acc_test = c(acc_test, d_temp$overall["Accuracy"])
    acc_test_low = c(acc_test_low, d_temp$overall["AccuracyLower"])
    acc_test_up = c(acc_test_up, d_temp$overall["AccuracyUpper"])
    
    Method = c(Method, names(C)[j])
    
  }
  
  cat("\r", j)
  
}

stopCluster(cl)

df = data.frame(acc_train = acc_train, 
                acc_test = acc_test,
                test_lower = acc_test_low,
                test_upper = acc_test_up,
                Method = as.factor(Method)
)

df_all = data.frame(acc_train = acc_all_train, 
                    acc_test = acc_all_test,
                    test_lower = acc_all_test_low,
                    test_upper = acc_all_test_up,
                    Method = "Full-model")

write.table(df, "sonar_rf_example/sonar_rf_predict_results.txt", 
            row.names = FALSE)

write.table(df_all, "sonar_rf_example/sonar_rf_predict_results_full_model.txt", 
            row.names = FALSE)

df_comb = rbind(df, df_all)

df_comb$Method = as.factor(df_comb$Method)

df_comb %>%
  group_by(Method) %>%
  select(acc_train) %>%
  mutate(avg_train_acc = mean(acc_train)) %>%
  arrange(desc(avg_train_acc)) %>%
  distinct(avg_train_acc)

p_acc_train = ggplot(df_comb, 
                     aes(reorder(Method, acc_train, FUN = median, decreasing = TRUE), acc_train)) + 
  geom_boxplot() +
  xlab("Method") +
  ylab("Train accuracy") +
  ggtitle("(A)") +
  theme(plot.title = element_text(size = 20, face = "bold"),
        text = element_text(size = 20)) +
  ylim(0.65, 0.95)

p_acc_train

df_temp = df_comb %>%
  group_by(Method) %>%
  mutate(avg_pred_acc = mean(acc_test),
         lower = mean(test_lower),
         upper = mean(test_upper)) %>%
  select(avg_pred_acc, lower, upper) %>%
  arrange(desc(avg_pred_acc)) %>%
  distinct(avg_pred_acc, .keep_all = TRUE)

df_temp

p_acc_test = ggplot(df_temp, 
                    aes(reorder(Method, avg_pred_acc, FUN = mean, decreasing = TRUE), avg_pred_acc)) +
  geom_pointrange(aes(ymin = lower, ymax = upper)) +
  xlab("Method") +
  ylab("Test accuracy") +
  ggtitle("(B)") +
  theme(plot.title = element_text(size = 20, face = "bold"),
        text = element_text(size = 20)) +
  ylim(0.65, 0.95)

p_acc_test

final_plot = gridExtra::grid.arrange(p_acc_train, p_acc_test, ncol = 2)

ggsave(filename="sonar_rf_example/test_accuracy.pdf", 
       plot = final_plot, 
       device = cairo_pdf, 
       width = 400, 
       height = 297, 
       units = "mm")

ggsave(filename="sonar_rf_example/test_accuracy.png", 
       plot = final_plot, 
       device = png, 
       width = 400, 
       height = 297, 
       units = "mm",
       dpi = 700)
