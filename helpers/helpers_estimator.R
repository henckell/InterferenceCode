omega_comp <- function(pai, eta, feat_functions,A,B){
  
  n <- nrow(A)
  
  E_pai <- E_eta <- matrix(ncol=length(feat_functions))
  
  for(i in 1:B){
    E_pai <- rbind(E_pai,unlist(lapply(feat_functions, FUN = function(f){
      mean(f(A,rbinom(n,1,pai)))
    })))
    
    E_eta <- rbind(E_eta,unlist(lapply(feat_functions, FUN = function(f){
      mean(f(A,rbinom(n,1,eta)))
    })))
  }
  
  E_pai <- c(1, colMeans(as.matrix(E_pai[-1,])))
  
  E_eta <- c(1, colMeans(as.matrix(E_eta[-1,])))
  
  
  omega_1 <- pai*(E_pai) - eta*E_eta
  omega_0 <- (1-pai)*E_pai - (1-eta)*E_eta
  
  return(list(omega_1=omega_1,omega_0=omega_0))
}

sandwich.var.est <- function(Z, Y, alpha.hat, omega_0, omega_1, adj){
  
  eps.diag.mat <- diag(c(Y - Z%*%alpha.hat)^2)
  
  n <- nrow(Z)
  p <- ncol(Z)
  bread <- 1/n*t(Z)%*%Z
  meat <- 1/n*t(Z)%*%eps.diag.mat%*%Z
  alpha.hat.var <- solve(bread)%*%meat%*%solve(bread)
  
  var.coef.mat <- c(omega_0+omega_1, 
                    omega_1, 
                    rep(0, length(adj)))
  
  var.est.tau <- as.numeric(t(var.coef.mat)%*%alpha.hat.var%*%var.coef.mat)
  
  return(var.est.tau)}

reparam_coefs <- function(alpha_hat){
  length_betas <- length(alpha_hat)/2
  beta_hat_0 <- c(alpha_hat[1:length_betas])
  beta_hat_1 <- c(alpha_hat[(length_betas+1):length(alpha_hat)])+c(alpha_hat[1:length_betas])
  return(list(beta_hat_0=beta_hat_0, beta_hat_1=beta_hat_1))
}