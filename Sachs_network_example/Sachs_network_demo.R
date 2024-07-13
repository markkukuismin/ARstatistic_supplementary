
# 8.11.2022

library(GGMnonreg)
library(energy)
library(minerva)
library(XICOR)
library(HDInterval)
library(igraph)
library(MASS)

source("functions/rho.R")

data("Sachs")

head(Sachs)

dim(Sachs)

set.seed(1)

#Sachs_npn = huge.npn(Sachs)

#Sachs_npn = Sachs[sample(1:nrow(Sachs), 1000), ]

p = ncol(Sachs)

#colnames(Sachs_npn) = colnames(Sachs)

C_rho_minus_one = matrix(0, p, p)

C_rho_minus_one[upper.tri(C_rho_minus_one)] = 1

ind = which(C_rho_minus_one == 1, arr.ind = T)

C_edge_list = data.frame(row_col = ind, ar = 0, 
                         hdi_95_lower = 0, 
                         hdi_95_upper = 0)

C_edge_list$hdi_95_lower[5] = 1

for(i in 1:nrow(ind)){
  
  temp = rho(Sachs[, ind[i, 1]], Sachs[, ind[i, 2]],
             return_dist = T)
  
  c_hdi = HDInterval::hdi(temp)

  rho_mean = rho(Sachs[, ind[i, 1]], Sachs[, ind[i, 2]])
    
  C_rho_minus_one[ind[i, 1], ind[i, 2]] = rho_mean
  
  C_edge_list$ar[i] = rho_mean
  
  C_edge_list$hdi_95_lower[i] = as.numeric(c_hdi[1])
  
  C_edge_list$hdi_95_upper[i] = as.numeric(c_hdi[2])
  
  cat("\r", i)
  
}

write.table(C_edge_list, "Sachs_network_example/C_ar_edge_list.txt")

C_rho_minus_one = C_rho_minus_one + t(C_rho_minus_one)

colnames(C_rho_minus_one) = rownames(C_rho_minus_one) = colnames(Sachs)

write.table(C_rho_minus_one, "Sachs_network_example/C_ar.txt")

# dcor,

C_dcor = matrix(0, ncol = ncol(Sachs), nrow = ncol(Sachs))

for(i in 1:nrow(ind)){
  
  C_dcor[ind[i, 1], ind[i, 2]] = dcor(Sachs[, ind[i, 1]], Sachs[, ind[i, 2]])

  cat("\r", i)
    
}

C_dcor = C_dcor + t(C_dcor)

colnames(C_dcor) = rownames(C_dcor) = colnames(Sachs)

write.table(C_dcor, "Sachs_network_example/C_dcor.txt")

##

# correlation,

C = cor(Sachs)

######

lambda = C_rho_minus_one[upper.tri(C_rho_minus_one)]

lambda = sort(lambda)

sparsity_ar = c()

for(i in 1:length(lambda)){
  
  sparsity_ar[i] = 1 - sum(C_rho_minus_one[upper.tri(C_rho_minus_one)] <= lambda[i])/(p*(p-1)/2)
  
}

##

lambda_dcor = C_dcor[upper.tri(C_dcor)]

lambda_dcor = sort(lambda_dcor)

sparsity_dcor = c()

for(i in 1:length(lambda_dcor)){
  
  sparsity_dcor[i] = 1 - sum(C_dcor[upper.tri(C_dcor)] <= lambda_dcor[i])/(p*(p-1)/2)
  
}

##

lambda_cor = C[upper.tri(C)]

lambda_cor = sort(abs(lambda_cor))

sparsity_cor = c()

for(i in 1:length(lambda_cor)){
  
  sparsity_cor[i] = 1 - sum(abs(C[upper.tri(C)]) <= lambda_cor[i])/(p*(p-1)/2)
  
}

all.equal(sparsity_ar, sparsity_dcor)

all.equal(sparsity_ar, sparsity_cor)

all.equal(sparsity_dcor, sparsity_cor)

##

sparsity = sparsity_ar

sparsity = sparsity[which(sparsity_dcor == sparsity_ar)]

names(sparsity) = which(sparsity_dcor == sparsity_ar)

# Determine the threshold based on the sparsity level,

which(sparsity_ar*(11*10/2) == 20) 

thr = sparsity[35]

thr = as.numeric(names(thr))

C_rho_sparse = C_rho_minus_one
C_dcor_sparse = C_dcor
C_sparse = C

C_rho_sparse[C_rho_minus_one <= lambda[thr]] = 0
C_dcor_sparse[C_dcor <= lambda_dcor[thr]] = 0
C_sparse[abs(C) <= lambda_cor[thr]] = 0

##

# a-r

A_rho = C_rho_sparse

A_rho[A_rho != 0] = 1

G_ar = graph_from_adjacency_matrix(A_rho, mode = "undirected", diag = F)

gsize(G_ar)

coords = igraph::layout_in_circle(G_ar)

##

# dcor

A_dcor = C_dcor_sparse

A_dcor[A_dcor != 0] = 1

G_dcor = graph_from_adjacency_matrix(A_dcor, mode = "undirected", diag = F)

gsize(G_dcor)

coords_dcor = igraph::layout_in_circle(G_dcor)

# Plot graphs

png("Sachs_network_example/ar_vs_dcor_graph.png",
    res = 300,
    width = 4000, height = 3000)

par(mfrow = c(1, 2), mar = c(0, 0.1, 0.1, 0.1))

plot(G_ar, layout = coords, vertex.size = 20, 
     vertex.color = "white",
     edge.color = "black",
     vertex.label.color = "black")

title("(A)", cex.main = 2, adj = 0, line = -5) # AR

##

plot(G_dcor, layout = coords_dcor, vertex.size = 20, 
     vertex.color = "white",
     edge.color = "black",
     vertex.label.color = "black")

title("(B)", cex.main = 2, adj = 0, line = -5) # dcor

dev.off()

pdf("Sachs_network_example/ar_vs_dcor_graph.pdf",
    width = 8, height = 7)

par(mfrow = c(1, 2), mar = c(0, 0.1, 0.1, 0.1))

plot(G_ar, layout = coords, vertex.size = 20, 
     vertex.color = "white",
     edge.color = "black",
     vertex.label.color = "black")

title("(A)", cex.main = 2, adj = 0, line = -5) # AR

##

plot(G_dcor, layout = coords_dcor, vertex.size = 20, 
     vertex.color = "white",
     edge.color = "black",
     vertex.label.color = "black")

title("(B)", cex.main = 2, adj = 0, line = -5) # dcor

dev.off()

##

# Correlation graph

A_cor = C_sparse

A_cor[A_cor != 0] = 1

G_cor = graph_from_adjacency_matrix(A_cor, mode = "undirected", diag = F)

coords = igraph::layout_in_circle(G_cor)

plot(G_cor, layout = coords, main = "cor")

##

degree(G_ar)

degree(G_dcor)

degree(G_cor)

gsize(G_ar)

gsize(G_dcor)

gsize(G_cor)

##

G_dif = igraph::difference(G_ar, G_dcor)

coords = igraph::layout_in_circle(G_dif)

plot(G_dif, layout = coords, main = "A-R vs. dcor graph: edges only in A-R graph")

#

G_dif_cor_vs_dcor = igraph::difference(G_cor, G_dcor)

coords_c = igraph::layout_in_circle(G_dif_cor_vs_dcor)

plot(G_dif_cor_vs_dcor, layout = coords, main = "cor vs. dcor graph: edges only in cor graph")

#

A_dif = igraph::get.adjacency(G_dif)

A_dif = as.matrix(A_dif)

C_rho_minus_one*(A_dif == 1)

#plot(Sachs[, "PIP2"], Sachs[, "PIP3"], cex = 0.5)

#z = kde2d(Sachs[, "PIP2"], Sachs[, "PIP3"], n = 200)

#contour(z, add = T, col = "red", lwd = 2)

j = which(degree(G_dif) != 0)

j = names(j)

A_temp = A_dif

only_in_AR = A_temp*C_rho_minus_one

panel_ar = function(x, y, ...) {
  i = names(which(colSums(round(Sachs - x)) == 0))
  ii = names(which(colSums(round(Sachs - y)) == 0))
  horizontal = (par("usr")[1] + par("usr")[2]) / 2
  vertical = (par("usr")[3] + par("usr")[4]) / 2
  if(only_in_AR[i, ii] == 0){
    text(horizontal, vertical, format(C_rho_minus_one[i, ii], digits = 2, nsmall = 2), cex = 2) 
  }
  if(only_in_AR[i, ii] != 0){
    text(horizontal, vertical, format(C_rho_minus_one[i, ii], digits = 2, nsmall = 2), cex = 4) 
  }
  text(horizontal, 0.3*vertical, paste0("(", format(C_dcor[i, ii], digits = 2, nsmall = 2), ")"), cex = 2)
}

# The corresponding distance correlation value is shown in parentheses

png("Sachs_network_example/ar_Sachs_scatterplot.png",
    res = 300,
    width = 5000, height = 4000)

plot(Sachs[, j], pch = ".", cex = 2, lower.panel = panel_ar, xaxt = "n", yaxt = "n")

dev.off()

pdf("Sachs_network_example/ar_Sachs_scatterplot.pdf",
    width = 9, height = 8)

plot(Sachs[, j], pch = ".", cex = 2, lower.panel = panel_ar, xaxt = "n", yaxt = "n")

dev.off()

lambda[thr]

##

G_dif_dcor = igraph::difference(G_dcor, G_ar)

coords_dcor = igraph::layout_in_circle(G_dif_dcor)

plot(G_dif_dcor, layout = coords_dcor, main = "Difference graph: edges only in dcor graph")

A_dif_dcor = igraph::get.adjacency(G_dif_dcor)

A_dif_dcor = as.matrix(A_dif_dcor)

C_dcor*(A_dif_dcor == 1)

#plot(Sachs[, "PKA"], Sachs[, "Plcg"], cex = 0.5)

#z = kde2d(Sachs[, "PKA"], Sachs[, "Plcg"], n = 200)

#contour(z, add = T, col = "red", lwd = 2)

j = which(degree(G_dif_dcor) != 0)

j = names(j)

#panel_dcor = function(x, y, ...) {
#  i = names(which(colSums(round(Sachs - x)) == 0))
#  j = names(which(colSums(round(Sachs - y)) == 0))
#  horizontal = (par("usr")[1] + par("usr")[2]) / 2
#  vertical = (par("usr")[3] + par("usr")[4]) / 2
#  text(horizontal, vertical, format(C_dcor[i, j], digits=2))
#}

A_temp = A_dif_dcor

only_in_dcor = A_temp*C_dcor

panel_dcor = function(x, y, ...) {
  i = names(which(colSums(round(Sachs - x)) == 0))
  ii = names(which(colSums(round(Sachs - y)) == 0))
  horizontal = (par("usr")[1] + par("usr")[2]) / 2
  vertical = (par("usr")[3] + par("usr")[4]) / 2
  if(only_in_dcor[i, ii] == 0){
    text(horizontal, vertical, format(C_dcor[i, ii], digits = 2, nsmall = 2)) 
  }
  if(only_in_dcor[i, ii] != 0){
    text(horizontal, vertical, format(C_dcor[i, ii], digits = 2, nsmall = 2), cex = 2) 
  }
  text(horizontal, 0.3*vertical, paste0("(", format(C_rho_minus_one[i, ii], digits = 2, nsmall = 2), ")"))
}

# The corresponding acceptance ratio statistic value is shown in parentheses

plot(Sachs[, j], pch = ".", cex = 2, lower.panel = panel_dcor, xaxt = "n", yaxt = "n")

##

which(C_rho_minus_one == max(C_rho_minus_one), arr.ind = T)

plot(Sachs[, "Mek"], Sachs[, "Raf"])
max(C_rho_minus_one)

which(C_rho_minus_one == min(C_rho_minus_one[upper.tri(C_rho_minus_one)]), arr.ind = T)

plot(Sachs[, "PIP3"], Sachs[, "Raf"])
min(C_rho_minus_one[upper.tri(C_rho_minus_one)])

##

which(C_dcor == max(C_dcor), arr.ind = T)

max(C_dcor)

which(C_dcor == min(C_dcor[upper.tri(C_dcor)]), arr.ind = T)

min(C_dcor[upper.tri(C_dcor)])

##

abs(C_edge_list$hdi_95_upper[c(5, 21, 38, 49, 55)] - C_edge_list$hdi_95_lower[c(5, 21, 38, 49, 55)]) < 0.013
