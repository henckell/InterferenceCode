### feature functions

# feature X1: fraction of treated neighbours
feat_X1 <- function(A, treat){
  A_tilde <- as.matrix(normalize_A(A))
  feat <- c(A_tilde%*%treat)
  feat[which(c(feat)=="NaN")] <- 0
  return(feat)
}

# feature X2: indicator if at least 50% of neighbours are treated
feat_X2 <- function(A, treat){
  num_neigs <- rowSums(A)
  feat <- ifelse(A%*%treat >=0.5*num_neigs,1,0) 
  return(feat)
}

# feature X3: fraction of treated neighbours of neighbours
feat_X3 <- function(A, W,typeofgraph="3dlatt"){
  A_blub <- A%*%A
  #directed 3d-latt
  if(typeofgraph=="3dlatt"){
    A_blub <- t(Matrix(apply(A_blub,1, function(r) ifelse(r>=1, 1,0)), sparse=T))
    
  }else{
    A_blub <- Matrix(apply(A_blub,1, function(r) ifelse(r>=1, 1,0)), sparse=T)
  }
  
  diag(A_blub) <-0
  feat <- feat_X1(A_blub, W)
  return(feat)
}
