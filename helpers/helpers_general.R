## generic helpers

IV_estimate <- function(K,X, Y){
  alpha_hat <- solve(t(K)%*%X)%*%t(K)%*%Y
  colnames(alpha_hat) <- NULL
  rownames(alpha_hat) <- NULL
  return(list(alpha_hat=alpha_hat))
}

reparam_coefs <- function(alpha_hat){
  length_betas <- length(alpha_hat)/2
  beta_hat_0 <- c(alpha_hat[1:length_betas])
  beta_hat_1 <- c(alpha_hat[(length_betas+1):length(alpha_hat)])+c(alpha_hat[1:length_betas])
  return(list(beta_hat_0=beta_hat_0, beta_hat_1=beta_hat_1))
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

est_comp_prep <- function(feat, treat, A, C){
  
  # SUTVA and no confounding
  X_OLS1 <- as.matrix(cbind(rep(1, length(treat)), treat))
  X_OLS1 <- unname(X_OLS1)
  
  X_OLS2 <- as.matrix(cbind(rep(1, length(treat)), treat, C))
  X_OLS2 <- unname(X_OLS2)
  
  # Interference with true features, no confounding
  X_OLS3 <- as.matrix(cbind(feat, feat* treat))
  X_OLS3 <- unname(X_OLS3)
  
  X_OLS4 <- as.matrix(cbind(feat, feat* treat, C))
  X_OLS4 <- unname(X_OLS4)
  
  return(list(X_OLS1=X_OLS1,
              X_OLS2=X_OLS2, 
              X_OLS3=X_OLS3, 
              X_OLS4=X_OLS4
  ))
}
