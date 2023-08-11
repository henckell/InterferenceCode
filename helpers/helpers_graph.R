## Graph function
A_2_d_lattice <- function(n, usedirect=TRUE){
  two_d_lattice_graph <- make_lattice(c(round(sqrt(n)),round(sqrt(n))), directed=usedirect)
  A <- get.adjacency(two_d_lattice_graph, type="both",attr=NULL, names=F, sparse=T)
  A <- t(A)
  A[sqrt(nrow(A))+2,2] <- 0
  A[2,sqrt(nrow(A))+2] <- 1
  A[2,1] <- 0
  A[1,2] <- 1
  
  return(A)
}


# 3-d-lattice with directed edges
A_3_d_lattice <- function(n, usedirect=TRUE){
  net <- make_lattice(c(round(n^(1/3)),round(n^(1/3)),round(n^(1/3))), directed=usedirect)
  A <- get.adjacency(net, type="both",attr=NULL, names=F, sparse=T)
  
  return(A)
}

## Erdoes-Renyi G(n,p) random graph
A_random_graph <- function(n, prob){
  # construct Adjacency matrix A with uniformly with bernoulli distr.
  A_entries_nonsym <- rbinom(n*n,size=1, prob = prob)
  A_matrix_nonsym <- matrix(A_entries_nonsym, ncol =n)
  sym_structure <- lower.tri(A_matrix_nonsym, diag = FALSE)
  sym_structure <- matrix(as.numeric(sym_structure),ncol= n,nrow=n)
  ind <- which(sym_structure ==0)
  A_matrix_nonsym[ind]<- 0
  A <- A_matrix_nonsym + t(A_matrix_nonsym) # this is now a symmetric adj. matrix
}

## Watts-Strogatz small-world graph
## Has size^dim many nodes
## eg. n =1000 nodes is received by dim=3, size=10
A_Watts_Strogatz <- function(size, p, mean.degree){
  net <- sample_smallworld(dim=1,size=size, nei=mean.degree, p=p, loops = FALSE, multiple = FALSE)
  A <- get.adjacency(net, type="both",attr=NULL, names=F, sparse=T)
  return(A)
}

A_family_graph <- function(n) {
  # split units into disjoint sets of random size 1 to 6
  # within each set, have a full graph
  group_assignment <- unlist(sapply(seq_len(n / 2), function(i) {
    rep(i, sample(1:6, 1))
  }))[seq_len(n)]
  if (group_assignment[n - 1] != group_assignment[n]) {
    group_assignment[n] <- group_assignment[n - 1]
  }
  A <- matrix(0, n, n)
  for (j in unique(group_assignment)) {
    #print(paste0("j= ",j))
    
    inds <- which(group_assignment == j)
    if(length(inds)>1){
      combs <- combn(inds, 2)
      for (k in seq_len(ncol(combs))) {
        A[combs[1, k], combs[2, k]] <- A[combs[2, k], combs[1, k]] <- 1
        #print(paste0(combs[1,k],",",combs[2,k]))
      }
    }
    
  }
  #plot(graph_from_adjacency_matrix(adjmatrix = A, mode = "undirected"))
  return(A)
}