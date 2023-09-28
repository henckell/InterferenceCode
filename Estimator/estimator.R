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

setwd("~/GitHub/InvarianceCode")

source("helpers/helpers_estimator.R")


estimator <- function(pai, eta, W,Y,C,adj,feat_functions,A,B=100,factor=c()){
  
  omega <- omega_comp(pai, eta, feat_functions,A,B)
  
  X <- do.call(cbind,lapply(feat_functions, FUN = function(f){
    f(A,W)}))
  
  design <- cbind(1,X,W,W*X,as.matrix(C[,adj]))
  
  df <- as.data.frame(cbind(Y,design))
  
  if(length(factor)!=0){
  for(i in 1:length(factor)){
    df[,factor[i]] <- as.factor(df[,factor[i]])
  }}
  
  df <- dummy_cols(df)
  
  df <- df[,!names(df)%in% factor]
  
  alpha_full <- lm(Y~.-1.,data=df)$coef
  alpha <- alpha_full[1:(2+2*length(feat_functions))]
  
  beta <- reparam_coefs(alpha)
  
  estimate<- omega$omega_1%*%beta$beta_hat_1 +
    omega$omega_0%*%beta$beta_hat_0
  
  var_estimate <- sandwich.var.est(Z=df, Y=Y, alpha.hat=alpha_full, 
                                   omega_0=omega$omega_0, 
                                   omega_1=omega$omega_1, adj=adj)
  
  return(list(hat_tau=estimate,hat_var_tau=var_estimate))
}
s




