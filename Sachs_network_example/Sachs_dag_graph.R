
# Create an undirected version of the directed acyclic graph (dag) reported in
# Sachs et al (2005),

library(GGMnonreg)

data("Sachs")

A_dag = matrix(0, ncol(Sachs), ncol(Sachs))

colnames(A_dag) = rownames(A_dag) = colnames(Sachs)

A_dag[c("Raf", "Mek", "Jnk", "P38", "PKA"), "PKC"] = 1

A_dag[c("Raf", "Mek", "Erk", "Jnk", "P38", "Akt"), "PKA"] = 1

A_dag[c("Mek"), "Raf"] = 1

A_dag[c("Erk"), "Mek"] = 1

A_dag[c("PIP2"), "Plcg"] = 1

A_dag[c("PIP3"), "Plcg"] = 1

A_dag[c("PIP2"), "PIP3"] = 1

A_dag[c("Akt"), "Erk"] = 1

A_dag[c("Akt"), "Erk"] = 1

A_dag[c("PKC"), "PIP2"] = 1

A_dag[c("PKC"), "Plcg"] = 1

A_dag[c("Akt"), "PIP3"] = 1

##

# dag is used only as a shorthand for the graph of Sachs et al (2005),

G_dag = graph_from_adjacency_matrix(A_dag, mode = "undirected", diag = F)

gsize(G_dag)

coords = igraph::layout_in_circle(G_dag)

plot(G_dag, layout = coords, main = "Sachs graph")

##

# Note! Run first the script "Sachs_network_demo.R"

G_dif = igraph::difference(G_ar, G_dag)

coords = igraph::layout_in_circle(G_dif)

plot(G_dif, layout = coords, main = "A-R vs. DAG graph: edges only in A-R graph")

##

G_dif = igraph::difference(G_dcor, G_ar)

coords = igraph::layout_in_circle(G_dif)

plot(G_dif, layout = coords, main = "dcor vs. DAG graph: edges only in dcor graph")
