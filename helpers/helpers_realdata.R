feat_X1 <- function(A, treat){
  A_tilde <- as.matrix(normalize_A(A))
  feat <- c(A_tilde%*%treat)
  feat[which(c(feat)=="NaN")] <- 0
  return(feat)
}

normalize_A <- function(A){
  weights <- 1/rowSums(A)
  weights[which(weights=="Inf")] <- 0
  weight_matrix <- matrix(rep(weights, dim(A)[1]), byrow=F, ncol=dim(A)[1])
  A_norm <- A*weight_matrix
  return(A_norm)
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



create.A.init = function(){
  
  # Output: adjacency matrix of swiss cantons with the same ordering as our data
  
  # Names
  cantons = c("GE", "VD", "FR", "NE", "VS", "BE", "JU", "SO", "BL", "BS", "AG", 
              "LU", "OW", "NW", "UR", "TI", "ZG", "SZ", "GL", "SG", "ZH", "SH",
              "AR", "AI",  "GR", "TG")
  
  # Initialize
  A.init = matrix(NA, nrow = 26, ncol = 26)
  colnames(A.init) = cantons
  rownames(A.init) = cantons
  
  # GE
  A.init[1,2] = 1
  
  # VD
  A.init[2,c(1,4,3,5,6)] = 1
  
  # FR
  A.init[3, c(2,6, 4)] = 1
  
  # NE
  A.init[4, c(7,6,2,3)] = 1
  
  # VS
  A.init[5, c(2,6,15,16)] = 1
  
  # BE
  A.init[6,c(3,2,5,15,13,14,12,8,7,4)] = 1
  
  # JU
  A.init[7,c(9,10,8,6,4)] = 1
  
  # SO
  A.init[8,c(9,11,6,7)] = 1
  
  # BL
  A.init[9,c(8,11,10,7)] = 1
  
  # BS
  A.init[10,c(7,9)] = 1
  
  # AG
  A.init[11,c(9,8,12, 17,21)] = 1
  
  # LU
  A.init[12,c(11,17,18,14,13,6)] = 1
  
  # OW
  A.init[13,c(12,6,14)] = 1
  
  # NW
  A.init[14,c(15,13,12,18,6)] = 1
  
  # UR
  A.init[15, c(14,6,5,16,25,19,18)] = 1
  
  # TI
  A.init[16, c(5,15, 25)] = 1
  
  # ZG
  A.init[17, c(18,21,11,12)] = 1
  
  # SZ
  A.init[18,c(17,14,15,19,20,21,12)] = 1
  
  # GL
  A.init[19,c(18,15,25,20)] = 1
  
  # SG
  A.init[20, c(25,24,23,26,21,18,19)] = 1
  
  # ZH
  A.init[21,c(20,26,22,11,17,18)] = 1
  
  # SH
  A.init[22,c(21,26)] = 1
  
  # AR
  A.init[23, c(20,24)] = 1
  
  # AI
  A.init[24, c(20,23)] = 1
  
  # GR
  A.init[25,c(20,19,15,16)] = 1
  
  # TG
  A.init[26,c(22,21,20)] = 1
  
  # NA to 0
  A.init[is.na(A.init)] = 0
  
  # Return adjacency matrix
  return(A.init)
}

reparam_coefs <- function(alpha_hat){
  length_betas <- length(alpha_hat)/2
  beta_hat_0 <- c(alpha_hat[1:length_betas])
  beta_hat_1 <- c(alpha_hat[(length_betas+1):length(alpha_hat)])+c(alpha_hat[1:length_betas])
  return(list(beta_hat_0=beta_hat_0, beta_hat_1=beta_hat_1))
}