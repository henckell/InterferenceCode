# Helpers: Simulation Study for the Paper
# 27.06.2022

# simulate.tau <- function(setting, nrep_graph, typeofgraph,error.type,
#                          type.causal.effect,B,prob,prob.rewiring,
#                          const, nval,
#                          nrep.omega,
#                          n.cores.nrep){
#   
#   foldername <- paste0(typeofgraph,"_", type.causal.effect, "_setting", setting)
#   
#   if (file.exists(foldername)) {
#     setwd(foldername)
#     
#   } else {
#     dir.create(foldername)
#     setwd(foldername)
#     
#   }
#   
#   set.seed(11)
#   seeds <- sample(c(0:10000),nrep_graph,replace = FALSE)
#   params <- get(paste0("setting", setting))
#   beta_0 <- params$beta_0
#   beta_1 <- params$beta_1
#   
#   Results <- mclapply(1:nrep_graph, function(s){
#     
#     set.seed(seeds[s])
#     print(paste0("nrep_graph ",s))
#     
#     
#     if (typeofgraph=="chain"){
#       A <- A_chain_graph(nval)
#       A <- Matrix(A, sparse = T)
#     } else if (typeofgraph=="ring"){
#       A <- A_ring_graph(nval)
#       A <- Matrix(A, sparse = T)
#     } else if (typeofgraph=="2dlatt"){
#       A <- A_2_d_lattice(nval)
#       #A <- as.matrix(A)
#       #A <- Matrix(A, sparse = T)
#     } else if (typeofgraph=="3dlatt"){
#       A <- A_3_d_lattice(nval)
#       #A <- as.matrix(A)
#       #A <- Matrix(A, sparse = T)
#     } else if (typeofgraph=="rand_pfix"){
#       A <- A_random_graph(nval, prob)
#       A <- Matrix(A, sparse = T)
#     } else if (typeofgraph=="rand_npfix"){
#       A <- A_random_graph(nval, const/nval)
#       A <- Matrix(A, sparse = T)
#     } else if (typeofgraph=="rand_npfix_growing"){
#       A <- A_random_graph(nval, const*nval^(-4/5))
#       A <- Matrix(A, sparse = T)
#     } else if (typeofgraph=="WS"){
#       mean.degree <- const
#       A <- A_Watts_Strogatz(nval, prob.rewiring, mean.degree)
#       #A <- as.matrix(A)
#       #A <- Matrix(A, sparse = T)
#     } else if (typeofgraph=="WS_growing"){
#       mean.degree <- const*nval^(-4/5)
#       A <- A_Watts_Strogatz(nval, prob.rewiring, mean.degree)
#       #A <- as.matrix(A)
#       #A <- Matrix(A, sparse = T)
#     } else if (typeofgraph=="family"){
#       A <- A_family_graph(nval)
#     } else {print("graph type not known")}
#     
#     n <- dim(A)[1] #effective nval (not always possible to construct certain graph of size nval)
# 
#     feat_names <- c("featX2", "featX3")
#     # simulate omega
#     omegas <- simulate.omega(feat_names, type.causal.effect,n,A, nrep.omega)
#     
#     omega_1 <- unname(c(1, omegas$omega_1))
#     omega_0 <- unname(c(1, omegas$omega_0))
#     
#     tau <- as.numeric(omega_1%*%beta_1-omega_0%*%beta_0)
#     
#     return(list(tau=tau, omega_1=omega_1, omega_0=omega_0))
#     
#   }, mc.cores=n.cores.nrep)#11
#   
#   print(paste0("typeofgraph", typeofgraph, "_nval",nval))
#   
#   filename <- paste0("taus_",typeofgraph,"_", type.causal.effect, "_setting", setting, "_nval", nval)
#   
#   filename <- gsub("[.]","",filename)
#   assign(filename, Results)
#   save(Results, file=paste(filename, ".Rda",sep=""))
#   
#   setwd("..")
#   
# }

IV_estimate <- function(K,X, Y){
  alpha_hat <- solve(t(K)%*%X)%*%t(K)%*%Y
  colnames(alpha_hat) <- NULL
  rownames(alpha_hat) <- NULL
  #betas <- reparam_coefs(alpha_hat)
  #beta_hat_0 <- betas$beta_hat_0
  #beta_hat_1 <- betas$beta_hat_1
  #return(list(beta_hat_0=beta_hat_0, beta_hat_1=beta_hat_1))
  return(list(alpha_hat=alpha_hat))
}

reparam_coefs <- function(alpha_hat){
  length_betas <- length(alpha_hat)/2
  beta_hat_0 <- c(alpha_hat[1:length_betas])
  beta_hat_1 <- c(alpha_hat[(length_betas+1):length(alpha_hat)])+c(alpha_hat[1:length_betas])
  return(list(beta_hat_0=beta_hat_0, beta_hat_1=beta_hat_1))
}

freedman_bootstrap_var <- function(X,IVs,feat,
                                   Y,B,omega_1,omega_0,
                                   type){
  
  K <-  as.matrix(cbind(feat, IVs))
  
  alpha_hat <- solve(t(K)%*%X)%*%t(K)%*%Y
  
  Y_hat <- X%*%alpha_hat
  resids <- Y - Y_hat
  
  meat <- lapply(1:nrow(X), function(i){
    return(cbind(IVs[i,])%*%IVs[i,])
    
  })
  meat <-  Reduce("+", meat) / length(meat)
  
  vec.meat <-  lapply(1:nrow(X), function(i){
    return(cbind(IVs[i,])%*%resids[i])
    
  })
  vec.meat <-  Reduce("+", vec.meat) / length(vec.meat)
  
  resids.transf <- sapply(1:nrow(X), function(i){
    resids[i] - IVs[i,]%*%meat%*% vec.meat
    
    
  })
  
  GATE.est.b <- c()
  

  b <- 1
  while(b < B){
    
    inds <- sample(1:nrow(X), replace = TRUE)
    X.b <- X[inds,]
    feats.b <- feat[inds,]
    IVs.b <- IVs[inds,]
    resids.transf.b <- resids.transf[inds]
    K.b <- as.matrix(cbind(feats.b, IVs.b))
    
    Y.b <- X.b%*%alpha_hat +resids.transf.b
    
    
    alpha_hat.b <- tryCatch({ests <- IV_estimate(K = K.b, X = X.b,Y=Y.b) 
    ests},
    error = function(e) NA)
    
    if(any(is.na(alpha_hat.b))){
      print("bootstrap IV not invertible.")
    }else{
      
      if(type=="IV2"){
        betas <- reparam_coefs(alpha_hat.b$alpha_hat)
        beta_1_hat <- betas$beta_hat_1
        beta_0_hat <- betas$beta_hat_0  
        
      }else if (type=="IV3"){
        betas <- reparam_coefs(alpha_hat.b$alpha_hat[2:(length(alpha_hat.b$alpha_hat))])
        beta_1_hat <- betas$beta_hat_1
        beta_0_hat <- betas$beta_hat_0  
      }
      
      GATE.est.b[b] <- c(omega_1[1:length( beta_0_hat )]%*%beta_1_hat - 
                           omega_0[1:length( beta_0_hat)]%*%beta_0_hat)
      b <- b+1
    }

  }
  
  return(GATE.est.b)
}

case_bootstrap_var <- function(X_OLS,Y,B,omega_1,omega_0){
  
  alpha_hat <- solve(t(X_OLS)%*%X_OLS)%*%t(X_OLS)%*%Y
  resids <- Y - X_OLS%*%alpha_hat
  
  GATE.est.b <- c()
  b <- 1
  while(b < B){
    inds <- sample(1:nrow(X_OLS), replace = TRUE)
    X_OLS.b <- X_OLS[inds,]
    resids.b <- resids[inds]
    Y.b <- X_OLS.b%*%alpha_hat + resids.b
    
    
    alpha_hat.b <- tryCatch({ests <- solve(t(X_OLS.b)%*%X_OLS.b)%*%t(X_OLS.b)%*%Y.b 
    ests},
    error = function(e) NA)
    
    if(any(is.na(alpha_hat.b))){
      print("bootstrap OLS not invertible.")
    }else{
    
    colnames(alpha_hat.b) <- NULL
    rownames(alpha_hat.b) <- NULL
    
    beta_1_hat <- reparam_coefs(alpha_hat.b[1:(length(alpha_hat.b)-2)])$beta_hat_1
    beta_0_hat <- reparam_coefs(alpha_hat.b[1:(length(alpha_hat.b)-2)])$beta_hat_0
    
    GATE.est.b[b] <- c(omega_1[1:length( beta_0_hat )]%*%beta_1_hat - 
                         omega_0[1:length( beta_0_hat)]%*%beta_0_hat)
    b <- b+1
  }
  }
  return(GATE.est.b)
}




sigm <- function(x){
  return(1/(1+exp(-x)))
}


normalize_A <- function(A){
  weights <- 1/rowSums(A)
  weights[which(weights=="Inf")] <- 0
  weight_matrix <- matrix(rep(weights, dim(A)[1]), byrow=F, ncol=dim(A)[1])
  A_norm <- A*weight_matrix
  return(A_norm)
}



# feature X1: indicator if at least 50% of neighbours are treated
feat_X1 <- function(A, treat){
  num_neigs <- rowSums(A)
  feat <- ifelse(A%*%treat >=0.5*num_neigs,1,0) 
  return(feat)
}

# feature X2: fraction of treated neighbours
feat_X2 <- function(A, treat){
  A_tilde <- as.matrix(normalize_A(A))
  feat <- c(A_tilde%*%treat)
  feat[which(c(feat)=="NaN")] <- 0
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
  feat <- feat_X2(A_blub, W)
  return(feat)
}

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

build_X_K_multiple_C <- function(Z, feat, treat, A, C, model){
  
  # SUTVA and no confounding
  X_OLS1 <- as.matrix(cbind(rep(1, length(treat)), treat))
  X_OLS1 <- unname(X_OLS1)
  
  # STUVA and observed confounding
  
  if(model==1){
    X_OLS2 <- as.matrix(cbind(rep(1, length(treat)), treat, C[,1],C[,2]))
    X_OLS2 <- unname(X_OLS2)
  }else if (model==2){
    X_OLS2 <- as.matrix(cbind(rep(1, length(treat)), treat, C[,1]))
    X_OLS2 <- unname(X_OLS2)
  }else if (model==3){
    X_OLS2 <- as.matrix(cbind(rep(1, length(treat)), treat, C[,1],C[,4]))
    X_OLS2 <- unname(X_OLS2)
  }
 
  
  
  # Interference with true features, no confounding
  X_OLS3 <- as.matrix(cbind(feat, feat* treat))
  X_OLS3 <- unname(X_OLS3)
  
  
  # Interference with true features, controlling for confounding with C
  if(model==1){
    X_OLS4 <- as.matrix(cbind(feat, feat* treat, C[,1],C[,2]))
    X_OLS4 <- unname(X_OLS4)
    
  }else if (model==2){
    X_OLS4 <- as.matrix(cbind(feat, feat* treat, C[,1]))
    X_OLS4 <- unname(X_OLS4)
  }else if (model==3){
    X_OLS4 <- as.matrix(cbind(feat, feat* treat, C[,1], C[,4]))
    X_OLS4 <- unname(X_OLS4)
  }
 
  # X_OLS4.simple <- as.matrix(cbind(feat, feat* treat, C1,C2))
  # X_OLS4.simple <- unname(X_OLS4.simple)
  
  
  # Using Z and K as instruments for W and XW, 
  
 
  
  
  return(list(X_OLS1=X_OLS1,
              X_OLS2=X_OLS2, 
              X_OLS3=X_OLS3, 
              X_OLS4=X_OLS4
              ))
}

build_X_K_multiple_C_EATE <- function(Z, feat, treat, A, C1,C2, model){
  
  # SUTVA and no confounding
  X_OLS1 <- as.matrix(cbind(rep(1, length(treat)), treat))
  X_OLS1 <- unname(X_OLS1)
  
  # STUVA and observed confounding
  
  if(model==1){
    X_OLS2 <- as.matrix(cbind(rep(1, length(treat)), treat, C1,C2))
    X_OLS2 <- unname(X_OLS2)
    X_OLS2 <- as.matrix(cbind(rep(1, length(treat)), treat, C1))
    X_OLS2 <- unname(X_OLS2)
  }
  
  
  # Interference with true features, no confounding
  X_OLS3 <- as.matrix(cbind(rep(1, length(treat)), treat, feat))
  X_OLS3 <- unname(X_OLS3)
  
  
  # Interference with true features, controlling for confounding with C
  if(model==1){
    X_OLS4 <- as.matrix(cbind(rep(1, length(treat)), treat, feat, C1, C2))
    X_OLS4 <- unname(X_OLS4)
    
  }else if (model==2){
    X_OLS4 <-  as.matrix(cbind(rep(1, length(treat)), treat, feat, C1))
    X_OLS4 <- unname(X_OLS4)
  }
  
  
  return(list(X_OLS1=X_OLS1,X_OLS2=X_OLS2, X_OLS3=X_OLS3, 
              X_OLS4=X_OLS4
              #X_OLS4.simple=X_OLS4.simple, 
              #X_IV1=X_IV1, K_IV1=K_IV1,X_IV2= X_IV2, 
              #K_IV2=K_IV2
              #X_IV3= X_IV3, K_IV3= K_IV3
  ))
}





data_model1_multiple_C <- function(A,beta_0, beta_1, delta_C1,delta_C2, delta_C1C2,
                                    sigma_Y,pai, eta,
                                  typeofgraph, use.only.feat1, 
                                   nval.eff, 
                                   error.type,
                                   do.intervene=NULL){
  
  
  if(error.type=="rnorm"){
    C2 <- rnorm(nval.eff,0,sigma_C)
    C1 <- delta_C1C2*C2 + rnorm(nval.eff,0,sigma_C)
    Z <- rnorm(nval.eff,0,sigma_Z)
    
  }else if (error.type =="runif"){
    # make sure that the sd is sigma_Y, that is sqrt(1/12(b-a)^2)
    C2 <- -2 +runif(nval.eff, -sqrt(12*sigma_C^2)/2, sqrt(12*sigma_C^2)/2)
    C1 <- delta_C1C2*C2 + runif(nval.eff, -sqrt(12*sigma_C^2)/2, sqrt(12*sigma_C^2)/2)
    Z <- 0.5+ runif(nval.eff, -sqrt(12*sigma_Z^2)/2, sqrt(12*sigma_Z^2)/2)
    
  }else if (error.type=="rt"){
    # make sure that the sd is sigma_Y, that is sqrt(df/(df-2))
    C2 <- rt(nval.eff, df=(2*sigma_C^2/(sigma_C^2-1)))
    C1 <- delta_C1C2*C2 + rt(nval.eff, df=(2*sigma_C^2/(sigma_C^2-1)))
    Z <- rt(nval.eff, df=(2*sigma_Z^2/(sigma_Z^2-1)))
    
  }else if(error.type=="chisq"){
    # make sure to have var sigma_Y^2 and mean 0
    C2 <- rchisq(nval.eff, df= sigma_C^2/2)- sigma_C^2/2
    C1 <- delta_C1C2*C2 + rchisq(nval.eff, df= sigma_C^2/2)- sigma_C^2/2
    Z <- rchisq(nval.eff, df= sigma_Z^2/2)- sigma_Z^2/2
    
  }
  
  
  if(is.null(do.intervene)){
    sigm_eval <- sigm(eta_C1*C1 + eta_Z*Z)
    treat <- rbinom(n = nval.eff, size = 1, prob = sigm_eval)
  }else {
    treat <- rbinom(nval.eff, size=1, prob=do.intervene)
  }
  

  # U <- runif(nval.eff, 0,1)
  #treat <- ifelse(U<sigm_eval, 1,0)
  #treat <- ifelse(sigm_eval<0.5, 0,1)
  #table(treat)
  num0 <- unname(table(treat)[1])
  num1 <- unname(table(treat)[2])
  
  print(paste0("num0: ", num0, " and num1 ", num1))
  
  
  if(typeofgraph=="family" | use.only.feat1 ==TRUE){
    featX2 <- feat_X2(A, treat)
    featX2 <- as.numeric(featX2)
    feat <- cbind(rep(1, length(treat)),featX2)
    
  }else{
    featX2 <- feat_X2(A, treat)
    featX2 <- as.numeric(featX2)
    featX3 <- feat_X3(A, treat)
    featX3 <- as.numeric(featX3)
    feat <- cbind(rep(1, length(treat)),featX2, featX3)
  }
  
  n <- length(treat)
  p <- ncol(feat)
  
  if(error.type=="rnorm"){
    error.Y <- rnorm(n,0,sigma_Y)
  }else if (error.type =="runif"){
    # make sure that the sd is sigma_Y, that is sqrt(1/12(b-a)^2)
    error.Y <- runif(n, -sqrt(12*sigma_Y^2)/2, sqrt(12*sigma_Y^2)/2)
  }else if (error.type=="rt"){
    # make sure that the sd is sigma_Y, that is sqrt(df/(df-2))
    error.Y <- rt(n, df=(2*sigma_Y^2/(sigma_Y^2-1)))
  }else if(error.type=="chisq"){
    # make sure that the sd is sigma_Y, that is sqrt(df/(df-2))
    error.Y <- rchisq(n, df= sigma_Y^2/2)-sigma_Y^2/2
  }
  
  
  Y <- treat *feat%*%matrix(beta_1, ncol=1)[1:p] + 
    (1-treat)*feat%*%matrix(beta_0, ncol=1)[1:p] +
    delta_C1*C1 + delta_C2*C2+ error.Y 
  
 
  
  if(is.null(do.intervene)){
    # simulate omega
    feat_names <- colnames(feat)[-1]
    
    omegas <- simulate.omega(feat_names, pai, eta,A)
    
    omega_1 <- unname(omegas$omega_1)
    omega_0 <- unname(omegas$omega_0)
    
    tau <- as.numeric(omega_1%*%beta_1[1:p] + omega_0%*%beta_0[1:p])
    
    
    return(list(Y=as.numeric(Y), tau=tau,
                feat=feat, omega_1=omega_1, omega_0=omega_0, 
                C1=C1, C2=C2, Z=Z, treat=treat, num0=num0, num1=num1))
  }else{
    return(Y)
  }

}


data_model2_multiple_C <- function(A,
                                   beta_0, 
                                   beta_1, 
                                   delta_C1,
                                   delta_C2, 
                                   delta_C1C2,
                                   sigma_Y,
                                   pai, 
                                   eta,
                                  typeofgraph, 
                                  use.only.feat1, 
                                   nval.eff, 
                                   error.type, 
                                   do.intervene = NULL){
  
  
  if(error.type=="rnorm"){
    C2 <- rnorm(nval.eff,0,sigma_C)
    C1 <- delta_C1C2*C2 + rnorm(nval.eff,0,sigma_C)
    Z <- rnorm(nval.eff,0,sigma_Z)
    
  }else if (error.type =="runif"){
    # make sure that the sd is sigma_Y, that is sqrt(1/12(b-a)^2)
    C2 <- -2 +runif(nval.eff, -sqrt(12*sigma_C^2)/2, sqrt(12*sigma_C^2)/2)
    C1 <- delta_C1C2*C2 + runif(nval.eff, -sqrt(12*sigma_C^2)/2, sqrt(12*sigma_C^2)/2)
    Z <- 0.5+ runif(nval.eff, -sqrt(12*sigma_Z^2)/2, sqrt(12*sigma_Z^2)/2)
    
  }else if (error.type=="rt"){
    # make sure that the sd is sigma_Y, that is sqrt(df/(df-2))
    C2 <- rt(nval.eff, df=(2*sigma_C^2/(sigma_C^2-1)))
    C1 <- delta_C1C2*C2 + rt(nval.eff, df=(2*sigma_C^2/(sigma_C^2-1)))
    Z <- rt(nval.eff, df=(2*sigma_Z^2/(sigma_Z^2-1)))
    
  }else if(error.type=="chisq"){
    # make sure to have var sigma_Y^2 and mean 0
    C2 <- rchisq(nval.eff, df= sigma_C^2/2)- sigma_C^2/2
    C1 <- delta_C1C2*C2 + rchisq(nval.eff, df= sigma_C^2/2)- sigma_C^2/2
    Z <- rchisq(nval.eff, df= sigma_Z^2/2)- sigma_Z^2/2
    
  }
  

  if(is.null(do.intervene)){
    sigm_eval <- sigm(eta_C1*C1 + eta_Z*Z)
    treat <- rbinom(n = nval.eff, size = 1, prob = sigm_eval)
  }else {
    treat <- rbinom(nval.eff, size=1, prob=do.intervene)
  }
  
 
  # U <- runif(nval.eff, 0,1)
  #treat <- ifelse(U<sigm_eval, 1,0)
  #treat <- ifelse(sigm_eval<0.5, 0,1)
  #table(treat)
  num0 <- unname(table(treat)[1])
  num1 <- unname(table(treat)[2])
  #print(table(treat))
  
  #print(paste0("num0: ", num0, " and num1 ", num1))
  
  
  if(typeofgraph=="family" | use.only.feat1 ==TRUE){
    featX2 <- feat_X2(A, treat)
    featX2 <- as.numeric(featX2)
    feat <- cbind(rep(1, length(treat)),featX2)
    
  }else{
    featX2 <- feat_X2(A, treat)
    featX2 <- as.numeric(featX2)
    featX3 <- feat_X3(A, treat)
    featX3 <- as.numeric(featX3)
    feat <- cbind(rep(1, length(treat)),featX2, featX3)
  }
  
  n <- length(treat)
  p <- ncol(feat)
  
  if(error.type=="rnorm"){
    error.Y <- rnorm(n,0,sigma_Y)
  }else if (error.type =="runif"){
    # make sure that the sd is sigma_Y, that is sqrt(1/12(b-a)^2)
    error.Y <- runif(n, -sqrt(12*sigma_Y^2)/2, sqrt(12*sigma_Y^2)/2)
  }else if (error.type=="rt"){
    # make sure that the sd is sigma_Y, that is sqrt(df/(df-2))
    error.Y <- rt(n, df=(2*sigma_Y^2/(sigma_Y^2-1)))
  }else if(error.type=="chisq"){
    # make sure that the sd is sigma_Y, that is sqrt(df/(df-2))
    error.Y <- rchisq(n, df= sigma_Y^2/2)-sigma_Y^2/2
  }
  
  
  Y <- treat *feat%*%matrix(beta_1, ncol=1)[1:p] + 
    (1-treat)*feat%*%matrix(beta_0, ncol=1)[1:p] +
    delta_C2*C2 + error.Y 
  
  if(is.null(do.intervene)){
    feat_names <- colnames(feat)[-1]
    
    # simulate omega
    omegas <- simulate.omega(feat_names, pai=pai, eta=eta,A, typeofgraph)
    
    omega_1 <- unname(omegas$omega_1)
    omega_0 <- unname(omegas$omega_0)
    
    tau <- as.numeric(omega_1%*%beta_1[1:p] + omega_0%*%beta_0[1:p])
    
    
    return(list(Y=as.numeric(Y), tau=tau,
                feat=feat, omega_1=omega_1, omega_0=omega_0, 
                C1=C1, C2=C2,Z=Z, treat=treat, num0=num0, num1=num1))
  }else{
    return(Y)
  }
  

}


data_model3_multiple_C <- function(A,
                                   beta_0, 
                                   beta_1, 
                                   delta_C1,
                                   delta_C2, 
                                   delta_C1C2,
                                   sigma_Y,
                                   pai, 
                                   eta,
                                   typeofgraph, 
                                   use.only.feat1, 
                                   nval.eff, 
                                   error.type, 
                                   do.intervene = NULL){
  
  
  if(error.type=="rnorm"){
    C2 <- rnorm(nval.eff,0,sigma_C)
    C1 <- delta_C1C2*C2 + rnorm(nval.eff,0,sigma_C)
    Z <- rnorm(nval.eff,0,sigma_Z)
    
  }else if (error.type =="runif"){
    # make sure that the sd is sigma_Y, that is sqrt(1/12(b-a)^2)
    C2 <- -2 +runif(nval.eff, -sqrt(12*sigma_C^2)/2, sqrt(12*sigma_C^2)/2)
    C3 <- 1 +runif(nval.eff, -sqrt(12*sigma_C^2)/2, sqrt(12*sigma_C^2)/2)
    C4 <- 1.3*C3 +runif(nval.eff, -sqrt(12*sigma_C^2)/2, sqrt(12*sigma_C^2)/2)
    C1 <- 2*C2 -1.5*C4+ runif(nval.eff, -sqrt(12*sigma_C^2)/2, sqrt(12*sigma_C^2)/2)
    Z <- 0.5+ runif(nval.eff, -sqrt(12*sigma_Z^2)/2, sqrt(12*sigma_Z^2)/2)
    
  }else if (error.type=="rt"){
    # make sure that the sd is sigma_Y, that is sqrt(df/(df-2))
    C2 <- rt(nval.eff, df=(2*sigma_C^2/(sigma_C^2-1)))
    C1 <- delta_C1C2*C2 + rt(nval.eff, df=(2*sigma_C^2/(sigma_C^2-1)))
    Z <- rt(nval.eff, df=(2*sigma_Z^2/(sigma_Z^2-1)))
    
  }else if(error.type=="chisq"){
    # make sure to have var sigma_Y^2 and mean 0
    C2 <- rchisq(nval.eff, df= sigma_C^2/2)- sigma_C^2/2
    C1 <- delta_C1C2*C2 + rchisq(nval.eff, df= sigma_C^2/2)- sigma_C^2/2
    Z <- rchisq(nval.eff, df= sigma_Z^2/2)- sigma_Z^2/2
    
  }
  
  
  if(is.null(do.intervene)){
    sigm_eval <- sigm(eta_C1*C1 + eta_Z*C3)
    treat <- rbinom(n = nval.eff, size = 1, prob = sigm_eval)
  }else {
    treat <- rbinom(nval.eff, size=1, prob=do.intervene)
  }
  
  
  # U <- runif(nval.eff, 0,1)
  #treat <- ifelse(U<sigm_eval, 1,0)
  #treat <- ifelse(sigm_eval<0.5, 0,1)
  #table(treat)
  num0 <- unname(table(treat)[1])
  num1 <- unname(table(treat)[2])
  #print(table(treat))
  
  #print(paste0("num0: ", num0, " and num1 ", num1))
  
  
  if(typeofgraph=="family" | use.only.feat1 ==TRUE){
    featX2 <- feat_X2(A, treat)
    featX2 <- as.numeric(featX2)
    feat <- cbind(rep(1, length(treat)),featX2)
    
  }else{
    featX2 <- feat_X2(A, treat)
    featX2 <- as.numeric(featX2)
    featX3 <- feat_X3(A, treat)
    featX3 <- as.numeric(featX3)
    feat <- cbind(rep(1, length(treat)),featX2, featX3)
  }
  
  n <- length(treat)
  p <- ncol(feat)
  
  if(error.type=="rnorm"){
    error.Y <- rnorm(n,0,sigma_Y)
  }else if (error.type =="runif"){
    # make sure that the sd is sigma_Y, that is sqrt(1/12(b-a)^2)
    error.Y <- runif(n, -sqrt(12*sigma_Y^2)/2, sqrt(12*sigma_Y^2)/2)
  }else if (error.type=="rt"){
    # make sure that the sd is sigma_Y, that is sqrt(df/(df-2))
    error.Y <- rt(n, df=(2*sigma_Y^2/(sigma_Y^2-1)))
  }else if(error.type=="chisq"){
    # make sure that the sd is sigma_Y, that is sqrt(df/(df-2))
    error.Y <- rchisq(n, df= sigma_Y^2/2)-sigma_Y^2/2
  }
  
  
  Y <- treat *feat%*%matrix(beta_1, ncol=1)[1:p] + 
    (1-treat)*feat%*%matrix(beta_0, ncol=1)[1:p] +
    delta_C2*C2 + 3*C4+ error.Y 
  
  if(is.null(do.intervene)){
    feat_names <- colnames(feat)[-1]
    
    # simulate omega
    omegas <- simulate.omega(feat_names, pai=pai, eta=eta,A, typeofgraph)
    
    omega_1 <- unname(omegas$omega_1)
    omega_0 <- unname(omegas$omega_0)
    
    tau <- as.numeric(omega_1%*%beta_1[1:p] + omega_0%*%beta_0[1:p])
    
    
    return(list(Y=as.numeric(Y), tau=tau,
                feat=feat, omega_1=omega_1, omega_0=omega_0, 
                C1=C1, C2=C2,C3=C3, C4=C4,Z=Z, treat=treat, num0=num0, num1=num1))
  }else{
    return(Y)
  }
  
  
}

data_model2_multiple_C_EATE <- function(A,
                                   beta_0, 
                                   beta_1, 
                                   delta_C1,
                                   delta_C2, 
                                   delta_C1C2,
                                   sigma_Y,
                                   pai, 
                                   eta,
                                   typeofgraph, 
                                   use.only.feat1, 
                                   nval.eff, 
                                   error.type, 
                                   do.intervene = NULL){
  
  
  if(error.type=="rnorm"){
    C2 <- -2 + rnorm(nval.eff,0,sigma_C)
    C1 <- delta_C1C2*C2 + rnorm(nval.eff,0,sigma_C)
    Z <- 0.5 + rnorm(nval.eff,0,sigma_Z)
    
  }else if (error.type =="runif"){
    # make sure that the sd is sigma_Y, that is sqrt(1/12(b-a)^2)
    C2 <- -2 +runif(nval.eff, -sqrt(12*sigma_C^2)/2, sqrt(12*sigma_C^2)/2)
    C1 <- delta_C1C2*C2 + runif(nval.eff, -sqrt(12*sigma_C^2)/2, sqrt(12*sigma_C^2)/2)
    Z <- 0.5+ runif(nval.eff, -sqrt(12*sigma_Z^2)/2, sqrt(12*sigma_Z^2)/2)
    
  }else if (error.type=="rt"){
    # make sure that the sd is sigma_Y, that is sqrt(df/(df-2))
    C2 <- -2 + rt(nval.eff, df=(2*sigma_C^2/(sigma_C^2-1)))
    C1 <- delta_C1C2*C2 + rt(nval.eff, df=(2*sigma_C^2/(sigma_C^2-1)))
    Z <- 0.5 +rt(nval.eff, df=(2*sigma_Z^2/(sigma_Z^2-1)))
    
  }else if(error.type=="chisq"){
    # make sure to have var sigma_Y^2 and mean 0
    C2 <- -2 + rchisq(nval.eff, df= sigma_C^2/2)- sigma_C^2/2
    C1 <- delta_C1C2*C2 + rchisq(nval.eff, df= sigma_C^2/2)- sigma_C^2/2
    Z <- 0.5 + rchisq(nval.eff, df= sigma_Z^2/2)- sigma_Z^2/2
    
  }
  
  
  if(is.null(do.intervene)){
    sigm_eval <- sigm(eta_C1*C1 + eta_Z*Z)
    treat <- rbinom(n = nval.eff, size = 1, prob = sigm_eval)
  }else {
    treat <- rbinom(nval.eff, size=1, prob=do.intervene)
  }
  
  
  # U <- runif(nval.eff, 0,1)
  #treat <- ifelse(U<sigm_eval, 1,0)
  #treat <- ifelse(sigm_eval<0.5, 0,1)
  #table(treat)
  num0 <- unname(table(treat)[1])
  num1 <- unname(table(treat)[2])
  #print(table(treat))
  
  #print(paste0("num0: ", num0, " and num1 ", num1))
  
  
  if(typeofgraph=="family" | use.only.feat1 ==TRUE){
    featX2 <- feat_X2(A, treat)
    featX2 <- as.numeric(featX2)
    feat <- cbind(featX2)
    
  }else{
    featX2 <- feat_X2(A, treat)
    featX2 <- as.numeric(featX2)
    featX3 <- feat_X3(A, treat)
    featX3 <- as.numeric(featX3)
    feat <- cbind(featX2, featX3)
  }
  
  n <- length(treat)
 
  
  if(error.type=="rnorm"){
    error.Y <- rnorm(n,0,sigma_Y)
  }else if (error.type =="runif"){
    # make sure that the sd is sigma_Y, that is sqrt(1/12(b-a)^2)
    error.Y <- runif(n, -sqrt(12*sigma_Y^2)/2, sqrt(12*sigma_Y^2)/2)
  }else if (error.type=="rt"){
    # make sure that the sd is sigma_Y, that is sqrt(df/(df-2))
    error.Y <- rt(n, df=(2*sigma_Y^2/(sigma_Y^2-1)))
  }else if(error.type=="chisq"){
    # make sure that the sd is sigma_Y, that is sqrt(df/(df-2))
    error.Y <- rchisq(n, df= sigma_Y^2/2)-sigma_Y^2/2
  }
  
  
  design.mat <- cbind(rep(1,length(treat)),treat,feat)
  
  p <- ncol(design.mat)
  
  Y <- design.mat%*%matrix(beta_1, ncol=1)[1:p] + 
    delta_C2*C2 + error.Y 
  
  if(is.null(do.intervene)){
    
    
    tau <- beta_1[2]
    
    
    return(list(Y=as.numeric(Y), tau=tau,
                feat=feat, 
                C1=C1, C2=C2,Z=Z, treat=treat, num0=num0, num1=num1))
  }else{
    return(Y)
  }
  
  
}

# W <- C1 -> Y, W <- C2 -> C2, P -> Y

# gen_Y_conf_multiple_C <- function(A, feat, treat, beta_0, beta_1, 
#                                   delta_C1,delta_C2,delta_P, 
#                                   C1,C2,P, 
#                                   sigma_Y, error.type,type.causal.effect,
#                                   nrep.omega){
#   n <- length(treat)
#   p <- ncol(feat)
#   
#   if(error.type=="rnorm"){
#     error.Y <- rnorm(n,0,sigma_Y)
#   }else if (error.type =="runif"){
#     # make sure that the sd is sigma_Y, that is sqrt(1/12(b-a)^2)
#     error.Y <- runif(n, -sqrt(12*sigma_Y^2)/2, sqrt(12*sigma_Y^2)/2)
#   }else if (error.type=="rt"){
#     # make sure that the sd is sigma_Y, that is sqrt(df/(df-2))
#     error.Y <- rt(n, df=(2*sigma_Y^2/(sigma_Y^2-1)))
#   }else if(error.type=="chisq"){
#     # make sure that the sd is sigma_Y, that is sqrt(df/(df-2))
#     error.Y <- rchisq(n, df= sigma_Y^2/2)-sigma_Y^2/2
#   }
#   
#   
#   Y <- treat *feat%*%matrix(beta_1, ncol=1)[1:p] + 
#     (1-treat)*feat%*%matrix(beta_0, ncol=1)[1:p] +
#     delta_C1*C1 + delta_C2*C2+ error.Y 
#   
#   feat_names <- colnames(feat)[-1]
#   
#   # simulate omega
#   omegas <- simulate.omega(feat_names, type.causal.effect,n,A, nrep.omega)
#   
#   omega_1 <- unname(c(1, omegas$omega_1))
#   omega_0 <- unname(c(1, omegas$omega_0))
#   
#   tau <- as.numeric(omega_1%*%beta_1[1:p]+omega_0%*%beta_0[1:p])
#   
#   return(list(Y=Y, tau=tau,omega_1=omega_1, omega_0=omega_0))
#   
# }

simulate.omega <- function(feat_names, pai, eta,A, typeofgraph){
  
  n <- nrow(A)
  
  if(typeofgraph=="2dlatt" && feat_names==c("featX2", "featX3")){
    
    E_pai <- c(1, pai, pai)
    E_eta <- c(1,eta,eta)
    
    omega_1 <- pai*(E_pai) - eta*E_eta
    omega_0 <- (1-pai)*E_pai - (1-eta)*E_eta
    
  }else{
    
    # this only works if the features are linear in Wi!!!!!!
    E_pai <- sapply(feat_names, FUN = function(nam){
      mean(get(paste0(substring(nam,1,4), "_", substring(nam,5,nchar(nam))))(A,rep(pai,n)))
    })
    
    E_pai <- c(1, E_pai)
    
    E_eta <- sapply(feat_names, FUN = function(nam){
      mean(get(paste0(substring(nam,1,4), "_", substring(nam,5,nchar(nam))))(A,rep(eta,n)))
    })
    
    E_eta <- c(1, E_eta)
    
    
    omega_1 <- pai*(E_pai) - eta*E_eta
    omega_0 <- (1-pai)*E_pai - (1-eta)*E_eta
  }
 

  
  return(list(omega_0=omega_0, omega_1=omega_1)) 
}


sandwich.var.est <- function(Z, Y, alpha.hat, omega_0, omega_1, type){
  
  # Z is n x p
  
  eps.diag.mat <- diag(c(Y - Z%*%alpha.hat)^2)
  
  n <- nrow(Z)
  p <- ncol(Z)
  bread <- 1/n*t(Z)%*%Z
  meat <- 1/n*t(Z)%*%eps.diag.mat%*%Z
  alpha.hat.var <- solve(bread)%*%meat%*%solve(bread)
  
  
  if(type=="OLS1"){
    
    len.alpha_hat.half <- length(alpha.hat)/2
    
  }else if(type=="OLS2"){
    
    if(model==1){
      len.alpha_hat.half <- (length(alpha.hat)-2)/2
      
    }else if (model==2){
      len.alpha_hat.half <- (length(alpha.hat)-1)/2
      
    }else if (model==3){
      len.alpha_hat.half <- (length(alpha.hat)-2)/2
      
    }
    
  }else if(type=="OLS3"){
    len.alpha_hat.half <- length(alpha.hat)/2
    
  } else if(type=="OLS4"){
    
    if(model==1){
      len.alpha_hat.half <- (length(alpha.hat)-2)/2
      
    }else if (model==2){
      len.alpha_hat.half <- (length(alpha.hat)-1)/2
    }else if (model==3){
      
      len.alpha_hat.half <- (length(alpha.hat)-2)/2
      
    }
  }
  
  
  
  var.coef.mat <- c(omega_0[1:len.alpha_hat.half]+ 
                      omega_1[1:len.alpha_hat.half], 
                    omega_1[1:len.alpha_hat.half], rep(0, (p-2*len.alpha_hat.half)))
  
  var.est.tau <- as.numeric(t(var.coef.mat)%*%alpha.hat.var%*%var.coef.mat)
  
  return(var.est.tau)
  
  
}





simulation <- function(  beta_0, 
                         beta_1, 
                         delta_C1, 
                         delta_C2, 
                         delta_C1C2, 
                         nrep_graph,
                         nrep_data,
                         n.cores.data,
                         n.cores.graph,
                         typeofgraph,
                         error.type,
                         pai, 
                         eta,
                         prob,
                         prob.rewiring,
                         sigma_Y, 
                         const, 
                         eta_C1, 
                         eta_C2, 
                         sigma_C, 
                         sigma_Z,
                         eta_Z, 
                         nval,
                         growth.rate,
                         growth.rate.WS, 
                         use.only.feat1, 
                         model,
                         eff,
                         useseed){
  
  set.seed(useseed)
  seeds <- sample(c(0:10000),nrep_graph,replace = FALSE)
  
  # is the same for all data reps, also when directed
  if(typeofgraph =="2dlatt"){
    A <- A_2_d_lattice(nval)
  }
  
  #for(nval in nvals){
  Results <- mclapply(1:nrep_graph, FUN = function(s){
    
    set.seed(seeds[s])
    print(paste0("nrep_graph ",s))
  
    if(typeofgraph != "2dlatt"){
      if (typeofgraph=="chain"){
        A <- A_chain_graph(nval)
        A <- Matrix(A, sparse = T)
      } else if (typeofgraph=="ring"){
        A <- A_ring_graph(nval)
        A <- Matrix(A, sparse = T)
      } else if (typeofgraph=="2dlatt"){
        A <- A_2_d_lattice(nval)
        #A <- as.matrix(A)
        #A <- Matrix(A, sparse = T)
      } else if (typeofgraph=="3dlatt"){
        A <- A_3_d_lattice(nval, usedirect=TRUE)
        #A <- as.matrix(A)
        #A <- Matrix(A, sparse = T)
      } else if (typeofgraph=="rand_pfix"){
        A <- A_random_graph(nval, prob)
        A <- Matrix(A, sparse = T)
      } else if (typeofgraph=="rand_npfix"){
        A <- A_random_graph(nval, const/nval)
        A <- Matrix(A, sparse = T)
      } else if (typeofgraph=="rand_npfix_growing"){
        A <- A_random_graph(nval, const*nval^(-growth.rate))
        A <- Matrix(A, sparse = T)
      } else if (typeofgraph=="WS"){
        mean.degree <- const
        A <- A_Watts_Strogatz(nval, prob.rewiring, mean.degree)
        #A <- as.matrix(A)
        #A <- Matrix(A, sparse = T)
      } else if (typeofgraph=="WS_growing"){
        mean.degree <- const*nval^(growth.rate.WS)
        A <- A_Watts_Strogatz(nval, prob.rewiring, mean.degree)
        #A <- as.matrix(A)
        #A <- Matrix(A, sparse = T)
      } else if (typeofgraph=="family"){
        A <- A_family_graph(nval)
      } else {print("graph type not known")}
    }

    
    mclapply(1:nrep_data, function(s){
      
     
      
      ########################################################################################
      ############################### true tau computation ####################################
      ########################################################################################
    
      
      
      nval.eff <- dim(A)[1] #effective nval (not always possible to construct certain graph of size nval)
      
      
      if(eff == "EATE"){
        dat <- get(paste0("data_model2_multiple_C_EATE"))(A=A,
                                                              beta_0=beta_0, 
                                                              beta_1=beta_1, 
                                                              delta_C1=delta_C1,
                                                              delta_C2=delta_C2,  
                                                              delta_C1C2=delta_C1C2,
                                                              sigma_Y= sigma_Y,
                                                              pai=pai,
                                                              eta=eta,
                                                              typeofgraph= typeofgraph, 
                                                              use.only.feat1=use.only.feat1, 
                                                              nval.eff=nval.eff, 
                                                              error.type=error.type,
                                                              do.intervene=NULL)
        
        num0 <- dat$num0
        num1 <- dat$num1
        C1 <- dat$C1
        C2 <- dat$C2
        treat <- dat$treat
        Z <- dat$Z
        feat <- dat$feat
        p <- ncol(feat)
        tau <- dat$tau
        Y <- dat$Y         
        num.betas <- ncol(feat)+2

        
      }else{
        dat <- get(paste0("data_model", model,"_multiple_C"))(A=A,
                                                              beta_0=beta_0, 
                                                              beta_1=beta_1, 
                                                              delta_C1=delta_C1,
                                                              delta_C2=delta_C2,  
                                                              delta_C1C2=delta_C1C2,
                                                              sigma_Y= sigma_Y,
                                                              pai=pai,
                                                              eta=eta,
                                                              typeofgraph= typeofgraph, 
                                                              use.only.feat1=use.only.feat1, 
                                                              nval.eff=nval.eff, 
                                                              error.type=error.type,
                                                              do.intervene=NULL)
        
        num0 <- dat$num0
        num1 <- dat$num1
        C1 <- dat$C1
        C2 <- dat$C2
        if(model %in% c(1,2)){
          C <- cbind(C1, C2)
          
        }else if(model==3){
          C3 <- dat$C3
          C4 <- dat$C4
          C <- cbind(C1,C2,C3,C4)
        }
        
        
        treat <- dat$treat
        Z <- dat$Z
        feat <- dat$feat
        p <- ncol(feat)
        tau <- dat$tau
        Y <- dat$Y
        omega_1 <- unname(dat$omega_1)
        omega_0 <- unname(dat$omega_0)
        num.betas <- length(omega_0)
        
        
      }
   
     
    
      ######################################################################################
      ############## Estimation ############################################################
      ######################################################################################
      
      if(eff=="EATE"){
        dat_est <-  build_X_K_multiple_C_EATE(Z, feat, treat, A, C1,C2, model=model)
      }else{
        dat_est <- build_X_K_multiple_C(Z, feat, treat, A, C, model=model)
        
      }
      
      X_OLS1 <- dat_est$X_OLS1
      X_OLS2 <- dat_est$X_OLS2
      X_OLS3 <- dat_est$X_OLS3
      X_OLS4 <- dat_est$X_OLS4
 
      
      
      ############################# OLS1  ###########################################################
      
      alpha_OLS1 <- tryCatch({alpha_hat <- solve(t(X_OLS1)%*%X_OLS1)%*%t(X_OLS1)%*%Y
      alpha_hat},
      error = function(e) NA)

      if(eff=="EATE"){
        if(any(is.na(alpha_OLS1))){
          est_OLS1 <- NA
          print("OLS1 singular")
        }else{
          est_OLS1 <- c(alpha_OLS1)[2]
          
        }
        
      }else{
        if(any(is.na(alpha_OLS1))){
          beta_hat_1.OLS1 <- NA
          beta_hat_0.OLS1 <- NA
          est_OLS1 <- NA
          var_OLS1.est <- NA
          print("OLS1 singular")
        }else{
          colnames(alpha_OLS1) <- NULL
          rownames(alpha_OLS1) <- NULL
          
          beta_hat_1.OLS1 <- reparam_coefs(alpha_OLS1)$beta_hat_1
          beta_hat_0.OLS1 <- reparam_coefs(alpha_OLS1)$beta_hat_0
          
          est_OLS1<- c(omega_1[1:length( beta_hat_0.OLS1 )]%*%beta_hat_1.OLS1 +
                         omega_0[1:length( beta_hat_0.OLS1)]%*%beta_hat_0.OLS1)
          
          var_OLS1.est <- NA
          
          # var_OLS1.est <- sandwich.var.est(Z=X_OLS1, Y=Y, alpha.hat=alpha_OLS1, 
          #                                  omega_0=omega_0, 
          #                                  omega_1=omega_1, type="OLS1")
        }
        
      }
            
     
      ############################# OLS2 ###########################################################
      
      
      alpha_OLS2 <- tryCatch({alpha_hat <- solve(t(X_OLS2)%*%X_OLS2)%*%t(X_OLS2)%*%Y
      alpha_hat},
      error = function(e) NA)
      
      if(eff=="EATE"){
        if(any(is.na(alpha_OLS2))){
          est_OLS2 <- NA
          print("OLS2 singular")
        }else{
          est_OLS2 <- c(alpha_OLS2)[2]
          
        }
        
      }else{

      if(any(is.na(alpha_OLS2))){
        beta_hat_1.OLS2 <- NA
        beta_hat_0.OLS2 <- NA
        est_OLS2 <- NA
        var_OLS2.est <- NA
        print("OLS2 singular")
      }else{
        colnames(alpha_OLS2) <- NULL
        rownames(alpha_OLS2) <- NULL
        
        if(model %in% c(1,3)){
          beta_hat_1.OLS2 <- reparam_coefs(alpha_OLS2[1:(length(alpha_OLS2)-2)])$beta_hat_1
          beta_hat_0.OLS2 <- reparam_coefs(alpha_OLS2[1:(length(alpha_OLS2)-2)])$beta_hat_0
          
        }else if (model==2){
          beta_hat_1.OLS2 <- reparam_coefs(alpha_OLS2[1:(length(alpha_OLS2)-1)])$beta_hat_1
          beta_hat_0.OLS2 <- reparam_coefs(alpha_OLS2[1:(length(alpha_OLS2)-1)])$beta_hat_0
          
        }
        
        est_OLS2 <- c(omega_1[1:length( beta_hat_0.OLS2 )]%*%beta_hat_1.OLS2 +
                       omega_0[1:length( beta_hat_0.OLS2)]%*%beta_hat_0.OLS2)
        
        var_OLS2.est <- NA
        # var_OLS2.est <- sandwich.var.est(Z=X_OLS2, Y=Y, alpha.hat=alpha_OLS2, 
        #                                  omega_0=omega_0, 
        #                                  omega_1=omega_1, type="OLS2")
      }
      }
      
      ############################# OLS3 ###########################################################
      
      
      alpha_OLS3 <- tryCatch({alpha_hat <- solve(t(X_OLS3)%*%X_OLS3)%*%t(X_OLS3)%*%Y
      alpha_hat},
      error = function(e) NA)
      
      if(eff=="EATE"){
        if(any(is.na(alpha_OLS3))){
          est_OLS3 <- NA
          print("OLS3 singular")
        }else{
          est_OLS3 <- c(alpha_OLS3)[2]
          
        }
        
      }else{

      if(any(is.na(alpha_OLS3))){
        beta_hat_1.OLS3 <- NA
        beta_hat_0.OLS3 <- NA
        est_OLS3 <- NA
        var_OLS3.est <- NA
        print("OLS3 singular")
      }else{
        colnames(alpha_OLS3) <- NULL
        rownames(alpha_OLS3) <- NULL
        
        beta_hat_1.OLS3 <- reparam_coefs(alpha_OLS3)$beta_hat_1
        beta_hat_0.OLS3 <- reparam_coefs(alpha_OLS3)$beta_hat_0
        
        est_OLS3<- c(omega_1[1:length( beta_hat_0.OLS3 )]%*%beta_hat_1.OLS3 +
                       omega_0[1:length( beta_hat_0.OLS3)]%*%beta_hat_0.OLS3)
        
        # var_OLS3.est <- sandwich.var.est(Z=X_OLS3, Y=Y, alpha.hat=alpha_OLS3, 
        #                                  omega_0=omega_0, 
        #                                  omega_1=omega_1, type="OLS3")
        
        var_OLS3.est<- NA
      }
      }
      
      ############################# OLS4  ###########################################################
      
      
      alpha_OLS4 <- tryCatch({alpha_hat <- solve(t(X_OLS4)%*%X_OLS4)%*%t(X_OLS4)%*%Y
      alpha_hat},
      error = function(e) NA)
      
      if(eff=="EATE"){
        if(any(is.na(alpha_OLS4))){
          est_OLS4 <- NA
          print("OLS4 singular")
        }else{
          est_OLS4 <- c(alpha_OLS4)[2]
          
        }
        
      }else{

      if(any(is.na(alpha_OLS4))){
        beta_hat_1.OLS4 <- NA
        beta_hat_0.OLS4 <- NA
        est_OLS4 <- NA
        var_OLS4.est <- NA
        print("OLS4 singular")
      }else{
        colnames(alpha_OLS4) <- NULL
        rownames(alpha_OLS4) <- NULL
        
        if(model%in% c(1,3)){
          beta_hat_1.OLS4 <- reparam_coefs(alpha_OLS4[1:(length(alpha_OLS4)-2)])$beta_hat_1
          beta_hat_0.OLS4 <- reparam_coefs(alpha_OLS4[1:(length(alpha_OLS4)-2)])$beta_hat_0   
        }else if (model==2){
          beta_hat_1.OLS4 <- reparam_coefs(alpha_OLS4[1:(length(alpha_OLS4)-1)])$beta_hat_1
          beta_hat_0.OLS4 <- reparam_coefs(alpha_OLS4[1:(length(alpha_OLS4)-1)])$beta_hat_0
        }
       
        
        est_OLS4<- c(omega_1[1:length( beta_hat_0.OLS4 )]%*%beta_hat_1.OLS4 +
                       omega_0[1:length( beta_hat_0.OLS4)]%*%beta_hat_0.OLS4)
        
        var_OLS4.est <- sandwich.var.est(Z=X_OLS4, Y=Y, alpha.hat=alpha_OLS4, 
                                         omega_0=omega_0, 
                                         omega_1=omega_1, type="OLS4")
        
      }
      }
        
      if(eff=="EATE"){
        return(list(
          nval.eff = nval.eff,
          num0=num0,
          num1=num1,
          est_OLS1 = est_OLS1,
          est_OLS2 = est_OLS2,
          est_OLS3 = est_OLS3,
          est_OLS4 = est_OLS4,
         tau = tau
        )
        )
      }else{
        return(list(
          omega_1 = omega_1,
          omega_0 = omega_0,
          nval.eff = nval.eff,
          num0=num0,
          num1=num1,
          est_OLS1 = est_OLS1,
          est_OLS2 = est_OLS2,
          est_OLS3 = est_OLS3,
          est_OLS4 = est_OLS4,
          beta_hat.OLS1 = list(beta_hat_1=beta_hat_1.OLS1,beta_hat_0=beta_hat_0.OLS1),
          beta_hat.OLS2 = list(beta_hat_1=beta_hat_1.OLS2,beta_hat_0=beta_hat_0.OLS2),
          beta_hat.OLS3 = list(beta_hat_1=beta_hat_1.OLS3,beta_hat_0=beta_hat_0.OLS3),
          beta_hat.OLS4 = list(beta_hat_1=beta_hat_1.OLS4,beta_hat_0=beta_hat_0.OLS4),
          tau = tau,
          var_OLS1.est = var_OLS1.est,
          var_OLS2.est = var_OLS2.est,
          var_OLS3.est = var_OLS3.est,
          var_OLS4.est = var_OLS4.est
        )
        )
      }
      
    }, mc.cores=n.cores.data)#11
    }, mc.cores=n.cores.graph)
    print(paste0("typeofgraph", typeofgraph, "_nval",nval))
    
    filename <- paste0(useseed,"estimation_",typeofgraph,"_pi",pai ,"_eta",eta, "_model", model, "_nval", nval)
    
    filename <- gsub("[.]","",filename)
    assign(filename, Results)
    save(Results, file=paste(filename, ".Rda",sep=""))
    
    
  }#, mc.cores=n.cores.nval)#8

  
#}
