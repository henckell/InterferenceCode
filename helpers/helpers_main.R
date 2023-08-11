

data_gen <- function(A,beta, 
                      B_C,eta_C,gamma_C,delta_C,sigma_C,
                      delta_Y,sigma_Y,
                      typeofgraph, 
                      features=features, 
                      nval.eff, 
                      error.type, 
                      do.intervene = NULL){
  
  if(error.type_C=="rnorm"){
    C <- DataGenNorm(nval.eff,B_C,delta_C,sigma_C)
  }else if (error.type =="runif"){
    # make sure that the sd is sigma_Y, that is sqrt(1/12(b-a)^2)
    C <- DataGenUnif(nval.eff,B_C,delta_C,-sqrt(12*sigma_C^2)/2,sqrt(12*sigma_C^2)/2)
  }else if (error.type=="rt"){
    # make sure that the sd is sigma_Y, that is sqrt(df/(df-2))
    C <- DataGenStu(nval.eff,B_C,delta_C,sigma_C^2)
  }
  
  if(is.null(do.intervene)){
    sigm_eval <- sigm(C%*%eta_C)
    treat <- rbinom(n = nval.eff, size = 1, prob = sigm_eval)
  }else {
    treat <- rbinom(nval.eff, size=1, prob=do.intervene)
  }
  
  num0 <- unname(table(treat)[1])
  num1 <- unname(table(treat)[2])
  
  featX1 <- feat_X1(A, treat)
  featX1 <- as.numeric(featX1)
  featX2 <- feat_X2(A, treat)
  featX2 <- as.numeric(featX2)
  featX3 <- feat_X3(A, treat)
  featX3 <- as.numeric(featX3)
  feat <- cbind(rep(1,length(treat)),featX1,featX2,featX3)

  features <- 1+features
  features <- c(1,features)
  
  feat <- feat[,features]
  
  n <- length(treat)
  p <- ncol(feat)
  
  if(error.type_Y=="rnorm"){
    error.Y <- delta_Y + rnorm(n,mean_Y,sigma_Y)
  }else if (error.type_Y =="runif"){
    error.Y <- delta_Y + runif(n, -sqrt(12*sigma_Y^2)/2, sqrt(12*sigma_Y^2)/2)
  }else if (error.type_Y=="rt"){
    error.Y <- delta_Y + rt(n, df=5) * sqrt(sigma_Y^2 * (5-2)/5)
  }
  
  Y <- treat *feat%*%matrix(beta[,"beta_1"], ncol=1)[1:p] + 
    (1-treat)*feat%*%matrix(beta[,"beta_0"], ncol=1)[1:p] +
    C%*%gamma_C + error.Y 
  
  if(is.null(do.intervene)){
    feat_names <- colnames(feat)[-1]
    
    # simulate omega
    omegas <- simulate.omega(feat_names, pai=pai, eta=eta,A)
    
    omega_1 <- unname(omegas$omega_1)
    omega_0 <- unname(omegas$omega_0)
    
    tau <- as.numeric(omega_1%*%beta_1[1:p] + omega_0%*%beta_0[1:p])
    
    
    return(list(Y=as.numeric(Y), tau=tau,
                feat=feat, omega_1=omega_1, omega_0=omega_0, 
                C=C, treat=treat, num0=num0, num1=num1))
  }else{
    return(Y)
  }
}

simulate.omega <- function(feat_names, pai, eta,A){

  n <- nrow(A)
  
  if("featX2" %in% feat_names){
    
    E_pai <- numeric(length(feat_names))
    E_eta <- numeric(length(feat_names))
    
    for(i in 1:100){
    E_pai <- E_pai + sapply(feat_names, FUN = function(nam){
      mean(get(paste0(substring(nam,1,4), "_", substring(nam,5,nchar(nam))))(A,rbinom(n,1,pai)))
    })
    
    E_eta <- E_eta + sapply(feat_names, FUN = function(nam){
      mean(get(paste0(substring(nam,1,4), "_", substring(nam,5,nchar(nam))))(A,rbinom(n,1,eta)))
    })
    }
    
    E_pai <- c(1, E_pai/100)
    
    E_eta <- c(1, E_eta/100)
    
    
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

simulation <- function(beta, 
                       B_C,
                       eta_C,
                       gamma_C,
                       delta_C,
                       sigma_C,
                       delta_Y,
                       sigma_Y,
                       nrep_graph,
                       nrep_data,
                       nval=nval,
                       typeofgraph,
                       error.types,
                       pai=pai,
                       eta=eta,
                       prob =prob, ## double check what this and const do
                       const=const,
                       prob.rewiring=prob.rewiring,
                       growth.rate=growth.rate,
                       growth.rate.WS=growth.rate.WS,
                       features=features,
                       useseed=useseed,
                       n.cores.data=n.cores.data,
                       n.cores.graph=n.cores.graph){
  
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
      } else if (typeofgraph=="3dlatt"){
        A <- A_3_d_lattice(nval, usedirect=TRUE)
      } else if (typeofgraph=="rand_pfix"){
        A <- A_random_graph(nval, prob)
        A <- Matrix(A, sparse = T)
      } else if (typeofgraph=="rand_npfix"){
        A <- A_random_graph(nval, const/nval)
        A <- Matrix(A, sparse = T)
      } else if (typeofgraph=="rand_npfix_growing"){
        A <- A_random_graph(nval, nval^(-growth.rate))
        A <- Matrix(A, sparse = T)
      } else if (typeofgraph=="WS"){
        mean.degree <- const
        A <- A_Watts_Strogatz(nval, prob.rewiring, mean.degree)
      } else if (typeofgraph=="WS_growing"){
        mean.degree <- const*nval^(growth.rate.WS)
        A <- A_Watts_Strogatz(nval, prob.rewiring, mean.degree)
      } else if (typeofgraph=="family"){
        A <- A_family_graph(nval)
      } else {print("graph type not known")}
    }

    
    mclapply(1:nrep_data, function(s){
      
      ############################### true tau computation ####################################
      
      nval.eff <- dim(A)[1] #effective nval (not always possible to construct certain graph of size nval)
      
      dat <- data_gen(A,beta, 
                      B_C,eta_C,gamma_C,delta_C,sigma_C,
                      delta_Y,sigma_Y,
                      typeofgraph,
                      features, 
                      nval.eff, 
                      error.type, 
                      do.intervene)
      
      
        
        num0 <- dat$num0
        num1 <- dat$num1
        C <- dat$C
        Z <- dat$Z
        treat <- dat$treat
        feat <- dat$feat
        p <- ncol(feat)
        tau <- dat$tau
        Y <- dat$Y
        omega_1 <- unname(dat$omega_1)
        omega_0 <- unname(dat$omega_0)
        num.betas <- length(omega_0)
    
      ######################################################################################
      ############## Estimation ############################################################
      ######################################################################################
      

      dat_est <- est_comp_prep(feat, treat, A, C[,adj])
      
      X_OLS1 <- dat_est$X_OLS1
      X_OLS2 <- dat_est$X_OLS2
      X_OLS3 <- dat_est$X_OLS3
      X_OLS4 <- dat_est$X_OLS4
      
      ############################# OLS1  ###########################################################
      
      alpha_OLS1 <- tryCatch({alpha_hat <- solve(t(X_OLS1)%*%X_OLS1)%*%t(X_OLS1)%*%Y
      alpha_hat},
      error = function(e) NA)

      
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
        }
  
      ############################# OLS2 ###########################################################
      
      alpha_OLS2 <- tryCatch({alpha_hat <- solve(t(X_OLS2)%*%X_OLS2)%*%t(X_OLS2)%*%Y
      alpha_hat},
      error = function(e) NA)

      if(any(is.na(alpha_OLS2))){
        beta_hat_1.OLS2 <- NA
        beta_hat_0.OLS2 <- NA
        est_OLS2 <- NA
        var_OLS2.est <- NA
        print("OLS2 singular")
      }else{
        colnames(alpha_OLS2) <- NULL
        rownames(alpha_OLS2) <- NULL

        beta_hat_1.OLS2 <- reparam_coefs(alpha_OLS2[1:(length(alpha_OLS2)-length(adj))])$beta_hat_1
        beta_hat_0.OLS2 <- reparam_coefs(alpha_OLS2[1:(length(alpha_OLS2)-length(adj))])$beta_hat_0
        
        est_OLS2 <- c(omega_1[1:length( beta_hat_0.OLS2 )]%*%beta_hat_1.OLS2 +
                       omega_0[1:length( beta_hat_0.OLS2)]%*%beta_hat_0.OLS2)
        
        var_OLS2.est <- NA
      }
      
      ############################# OLS3 ###########################################################
      
      
      alpha_OLS3 <- tryCatch({alpha_hat <- solve(t(X_OLS3)%*%X_OLS3)%*%t(X_OLS3)%*%Y
      alpha_hat},
      error = function(e) NA)

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
        
        var_OLS3.est<- NA
      }
      
      ############################# OLS4  ###########################################################
    
      alpha_OLS4 <- tryCatch({alpha_hat <- solve(t(X_OLS4)%*%X_OLS4)%*%t(X_OLS4)%*%Y
      alpha_hat},
      error = function(e) NA)

      if(any(is.na(alpha_OLS4))){
        beta_hat_1.OLS4 <- NA
        beta_hat_0.OLS4 <- NA
        est_OLS4 <- NA
        var_OLS4.est <- NA
        print("OLS4 singular")
      }else{
        colnames(alpha_OLS4) <- NULL
        rownames(alpha_OLS4) <- NULL
        
        beta_hat_1.OLS4 <- reparam_coefs(alpha_OLS4[1:(length(alpha_OLS4)-length(adj))])$beta_hat_1
        beta_hat_0.OLS4 <- reparam_coefs(alpha_OLS4[1:(length(alpha_OLS4)-length(adj))])$beta_hat_0   
        
        est_OLS4<- c(omega_1[1:length( beta_hat_0.OLS4 )]%*%beta_hat_1.OLS4 +
                       omega_0[1:length( beta_hat_0.OLS4)]%*%beta_hat_0.OLS4)
        
        var_OLS4.est <- sandwich.var.est(Z=X_OLS4, Y=Y, alpha.hat=alpha_OLS4, 
                                         omega_0=omega_0, 
                                         omega_1=omega_1, adj=adj)
        
      }
        
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
      
    }, mc.cores=n.cores.data)#11
    }, mc.cores=n.cores.graph)
    print(paste0("typeofgraph", typeofgraph, "_nval",nval))
    
    filename <- paste0("estimation_",typeofgraph,"_nval", nval)
    
    filename <- gsub("[.]","",filename)
    assign(filename, Results)
    save(Results, file=paste(filename, ".Rda",sep=""))
}
