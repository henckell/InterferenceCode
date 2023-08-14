### function implementing the estimator proposed in "A Graphical Approach to Treatment Effect Estimation with Observational Network Data"

### Parameters
## W: treatment vector
## Y: outcome vector
## C: matrix of covariates
## adj: vector of indicators showing which covariates to adjust for
## features:  list of feature functions
## A: adjacency matrix for the network graph
## Y: outcome vector
## B: number of repetitions used to approximate the omega weights 
## should be as high as possible but imposes large computational burden

### Returns
## hat_tau: estimate of target global treatment effect
## hat_var_tau: estimated variance of estimator for asy. valid CIs




estimator <- function(pai, eta, W,Y,C,adj,feat_functions,A,B=100){
  
  omega <- omega_comp(pai, eta, feat_functions,A,B)
  
  X <- do.call(cbind,lapply(feat_functions, FUN = function(f){
    f(A,W)}))
  
  design <- cbind(1,X,W,W*X,C[,adj])
  
  alpha_full <- lm(Y~design-1)$coef
  alpha <- alpha_full[1:(2+2*length(features))]
  
  beta <- reparam_coefs(alpha)
  
  estimate<- omega$omega_1%*%beta$beta_hat_1 +
    omega$omega_0%*%beta$beta_hat_0
  
  var_estimate <- sandwich.var.est(Z=design, Y=Y, alpha.hat=alpha_full, 
                                   omega_0=omega$omega_0, 
                                   omega_1=omega$omega_1, adj=adj)
  
  return(list(hat_tau=estimate,hat_var_tau=var_estimate))
}


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


