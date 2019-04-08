#############################################################################
# M1 MABS Math - Projet Math_Bio 2014-2015                                  #
# Professeur : Roland Barriot                                               #
# Etudiant : BENABDERRAHMANE MOHAMMED                                       # 
#############################################################################


# #Intro ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
# Le projet consiste à implémenter en R l’algorithme MCL de partitionnement de graphe (clustering) : 
# à partir d’un graphe et d’un paramètre nommé inflate factor (noté IF par la suite), 
# l’algorithme découpe le graphe en clusters. 
# 
# Un graphe est un objet mathématique composé de sommets (vertex) reliés par des arètes (edge). 
# Une des représentations informatiques consiste en une matrice carrée nommée matrice d’adjacence. 
# Une cellule de cette matrice contient le poids de l’arète reliant deux sommets ou zéro en absence de lien.
# 
# On distingue les graphes orientés (directed graph) des graphes non orientés (undirected graph). 
# Dans un graphe orienté, les arètes sont appelées arcs et ont un sens (de A vers B). 
# La matrice d’adjacence n’est donc plus obligatoirement symétrique et un choix arbitraire est fait concernant son interprétation:
# par exemple, les sommets sources correspondent aux lignes et les sommets destinations correspondent aux colonnes.  
# 
# 
# Principes et algorithme de MCl
# Une bonne méthode de clustering maximise la similarité des objets appartenant à un même cluster (cohésion)
# la dissemblance entre objets appartenant à des clusters différents (séparation).
# 
# En ce qui concerne le partitionnement de graphes, une bonne méthode vise à ce que les clusters obtenus 
# 
# aient un maximum d’arètes ou d’arcs reliant les sommets d’un même cluster
# aient un minimum d’arètes ou d’arcs reliant les sommets de clusters différents
# 
# Pour identifier les clusters, MCL se base sur la remarque suivante : 
# si on place un marcheur aléatoire sur un sommet d’un cluster, 
# il a plus de chances (en se promenant aléatoirement d’un sommet à un autre) de rester dans le cluster (plus grande densité de liens) que d’en sortir.
# 
# Le procédé s’appuie donc sur le calcul des probabilités de marches aléatoires (probabilité de passer par un sommet et probabilité d’emprunter un arc) dans un graphe donné.
# Le calcul s’effectue sur des matrices de Markov qui représentent les probabilités de transition d’un sommet à un autre ;
# il s’agit donc d’une matrice d’adjacence dont la somme des lignes vaut 1 (la somme des probabilités de sortie d’un sommet vaut 1). 
# 
# 
# L’algorithme MCL simule des marches aléatoires dans un graphe en alternant deux oprations appelées expansion et inflation. 
# L’expansion correspond au calcul de marches aléatoires de grande longueur et coincide à élever une matrice stochastique à une certaine puissance. 
# Concrètement, l’expansion consiste à multiplier la matrice avec elle-même avec le produit matriciel.
# L’inflation consiste à exagérer les probabilités de marches au sein d’un même cluster et à atténuer les probabilités de marches inter-clusters
# En pratique, l’inflation consiste à élever chaque cellule de la matrice à une certaine puissance (inflate factor) puis à normaliser les valeurs obtenues afin que la matrice soit stochastique.
# 
# 
# Au final, à force de répéter les opération d’expansion et d’inflation sur la matrice de transition (et donc sur le graphe), 
# celui-ci est décomposé en différentes composantes (composantes connexes) déconnectées les unes des autres 
# et qui correspondent aux clusters. En d’autres termes, l’algorithme converge et la matrice n’évolue plus.
# 
# 
# L’algorithme est donc le suivant :
#   
#   M : matrice d’adjacence
# 
# ajouter les boucles à M
# M_1 la matrice stochastique obtenue à partir de M
# tant que on observe un changement dans la matrice
# faire :
#   M_2 ←← expansion(M_1)
# M_2 ←← inflation(M_2, inflate_factor)
# déterminer s’il y a eu un changement
# M_1 ←← M_2




# Acquisition des données 
install.packages('igraph',dependencies = TRUE)
library(igraph)
g= read.graph('http://silico.biotoul.fr/enseignement/m1/math/projet/toy.tgr',directed=FALSE)
plot(g)


# Matrice d'adjacence du graph 
M = get.adjacency(g)
M= as.matrix(M)

add_loops <- function(M){
  diag(M) = 1
  return(M)
}

Make_stochastic <- function(M) {
  sum_row=apply(M,1,sum)
  M = M / sum_row
  return ( M )
}

Expansion = function (M ) {
  M = M %*% M
  return (M)
}


inflation = function(M,inflate) {
  M = M^inflate
  M = M / apply(M,1,sum)
  return (M)
}



# Fonction Chaos 
Chaos = function(M){
  max=apply(M,1,max)
  sum_sq=apply( M^2, 1, sum )
  row_chaos= max-sum_sq
  return(row_chaos)
}
########

########
MCL= function(M,inflate) {
  M2=add_loops(M)
  M2=Make_stochastic(M2)
  change = 1
  while (change > 0.001) {
    M2=Expansion(M2)
    M2=inflation(M2,inflate)
    change=Chaos(M2)
  }
  return(round(M2,3))
}
########


# Ignore warnings ... 
options(warn=-1)

??plot()
#~~~~~~~~~~~~~
# Tests : Graphes avec différentes valeurs inflate : 2 , 4 , 6 , 8 , 10, 15, 20 , 100  
#~~~~~~~~~~~~~

#1. inflate = 2
m = MCL(as.matrix(get.adjacency(g)), inflate=2)
g2=graph.adjacency(m, mode='undirected', weighted=T)
clusters(g2)
plot(g2, vertex.color=clusters(g2)$membership, edge.width=E(g)$weight)

#2. inflate = 4
m = MCL(as.matrix(get.adjacency(g)), inflate=4)
g2=graph.adjacency(m, mode='undirected', weighted=T)
clusters(g2)
plot(g2, vertex.color=clusters(g2)$membership, edge.width=E(g)$weight)

#3. inflate = 6
m = MCL(as.matrix(get.adjacency(g)), inflate=6)
g2=graph.adjacency(m, mode='undirected', weighted=T)
clusters(g2)
plot(g2, vertex.color=clusters(g2)$membership, edge.width=E(g)$weight)

#4. inflate = 8
m = MCL(as.matrix(get.adjacency(g)), inflate=8)
g2=graph.adjacency(m, mode='undirected', weighted=T)
clusters(g2)
plot(g2, vertex.color=clusters(g2)$membership, edge.width=E(g)$weight)

#5. inflate = 10
m = MCL(as.matrix(get.adjacency(g)), inflate=10)
g2=graph.adjacency(m, mode='undirected', weighted=T)
clusters(g2)
plot(g2, vertex.color=clusters(g2)$membership, edge.width=E(g)$weight)

#6. inflate = 15
m = MCL(as.matrix(get.adjacency(g)), inflate=15)
g2=graph.adjacency(m, mode='undirected', weighted=T)
clusters(g2)
plot(g2, vertex.color=clusters(g2)$membership, edge.width=E(g)$weight)

#7. inflate = 20
m = MCL(as.matrix(get.adjacency(g)), inflate=20)
g2=graph.adjacency(m, mode='undirected', weighted=T)
clusters(g2)
plot(g2, vertex.color=clusters(g2)$membership, edge.width=E(g)$weight)

#8. inflate = 100
m = MCL(as.matrix(get.adjacency(g)), inflate=100)
g2=graph.adjacency(m, mode='undirected', weighted=T)
clusters(g2)
plot(g2, vertex.color=clusters(g2)$membership, edge.width=E(g)$weight)

########################################################











