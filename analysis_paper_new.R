# Analysis: Simulation Study for the Paper
# 19.06.2023

################################### simulate tau
nvals <- c(300, 600, 1200,2400, 4800)#,9600)

setwd("~/GitHub/InvarianceCode/")


################################## estimation

library("latex2exp")
library(scales)

typeofgraph <- "2dlatt" #"rand_npfix" #"family","rand_npfix" #"WS_growing", "rand_npfix_growing", "family"
growth.rate <- 1/4#7/8#7/8#9/10#1/2#7/8 #3/4

eff<- "global"#"EATE"


nvals <- c(300, 600, 1200,2400, 4800)
dateused <- "2023-06-17"


use.only.feat1 <- 0
#confint.constr <-   "boot" #"normal" 
model <- 2
pai <- 0.5
eta <- 0.1
dateused <- "2023-06-17"
const <-10#


generate_plot <- function(typeofgraph,
                          nvals,
                          dateused, 
                          parts,
                          use.only.feat1, 
                          growth.rate.WS
                          ){
  

  if(typeofgraph=="WS_growing"){
    foldername <- paste0(dateused,"global_model",model,"_growth.rate.WS",round(growth.rate.WS,3),
                         "use.only.feat1",use.only.feat1,"WS_growing_pi",pai,"_eta",eta)
    
  }else{
    foldername <- paste0(dateused,"global_model",model,typeofgraph,"_growth.rate.WS",round(growth.rate,3),
                         "use.only.feat1",use.only.feat1,typeofgraph,"_pi",pai,"_eta",eta)
    
  }

 
    
  setwd(foldername)
  tau <- list()
  Res <- list()

  pval.normal_ols1 <- list()
  pval.normal_ols2 <- list()
  pval.normal_ols3 <- list()
  pval.normal_ols4 <- list()
  
  nvals.names <- nvals
  
  # effective nvals vary for the grid
  if(typeofgraph =="3dlatt"){
    nvals <- c(343, 512, 1331, 2197, 4913)
    nvals <- nvals[1:length(nvals.names)]
  }else if (typeofgraph =="2dlatt"){
    nvals <- c(289, 576, 1225, 2401, 4761)
    nvals <- nvals[1:length(nvals.names)]
    
  }
  
  
  
  
  for(i in 1:length(nvals)){
    nval <- nvals.names[i]
  
       
        if(typeofgraph=="family"){
          filename <- paste0(eff,"estimation_",typeofgraph,"_pi",ifelse(pai==1,1,10)*pai,"_eta",10*eta,"_model",model,"_nval",nval,".Rda")
          
        }else{
          filename <- paste0("1estimation_",typeofgraph,"_pi0",ifelse(pai==1,1,10)*pai,"_eta0",10*eta,"_model",model,"_nval",nval,".Rda")
          
        }
      load(filename)
      
      pval.normal_ols1[[i]] <- sapply(1:length(Results), function(k){

        dat <- sqrt(nvals[i])*(unlist(Results[[k]])[which(names(unlist(Results[[k]]))=="est_OLS1")])
        shapiro.test(dat )$p.value

      })



      pval.normal_ols2[[i]] <- sapply(1:length(Results), function(k){


        shapiro.test(sqrt(nvals[i])*(unlist(Results[[k]])[which(names(unlist(Results[[k]]))=="est_OLS2")]) )$p.value

      })

      pval.normal_ols3[[i]] <- sapply(1:length(Results), function(k){


        shapiro.test(sqrt(nvals[i])*(unlist(Results[[k]])[which(names(unlist(Results[[k]]))=="est_OLS3")]) )$p.value

      })

      pval.normal_ols4[[i]] <- sapply(1:length(Results), function(k){


        shapiro.test(sqrt(nvals[i])*(unlist(Results[[k]])[which(names(unlist(Results[[k]]))=="est_OLS4")]) )$p.value

      })


      Results <- unlist(Results, recursive=FALSE)
      
      
      
     
      
      taus <- unlist(do.call(rbind, Results)[,"tau"])
      tau[[i]] <- taus
      Res[[i]] <- unlist(Results)
      
      
      rm(Results)
      print(paste0("nval: ",nvals[i]))
    
    
    
  }
  
  pval.normal_ols1 <- do.call(cbind, pval.normal_ols1)
  pval.normal_ols1 <- data.frame(pval.normal_ols1)
  
  ecdf.ols1 <- sapply(1:ncol(pval.normal_ols1), function(i){
   
     ecdf(pval.normal_ols1[, i])
    
  })
    
   
  pval.normal_ols1 <- sapply(1:ncol(pval.normal_ols1), function(i){

    mean(pval.normal_ols1[,i]<=0.05)
  })
  names(pval.normal_ols1) <- paste0(nvals)
  
  
  
  pval.normal_ols2 <- do.call(cbind, pval.normal_ols2)
  pval.normal_ols2 <- data.frame(pval.normal_ols2)
  
  ecdf.ols2 <- sapply(1:ncol(pval.normal_ols2), function(i){
    
    ecdf(pval.normal_ols2[, i])
    
  })
  
  pval.normal_ols2 <- sapply(1:ncol(pval.normal_ols2), function(i){

    mean(pval.normal_ols2[,i]<=0.05)
  })
  names(pval.normal_ols2) <- paste0(nvals)
  
  
  pval.normal_ols3 <- do.call(cbind, pval.normal_ols3)
  pval.normal_ols3 <- data.frame(pval.normal_ols3)
  
  ecdf.ols3 <- sapply(1:ncol(pval.normal_ols3), function(i){
    
    ecdf(pval.normal_ols3[, i])
    
  })
  
  pval.normal_ols3 <- sapply(1:ncol(pval.normal_ols3), function(i){

    mean(pval.normal_ols3[,i]<=0.05)
  })
  names(pval.normal_ols3) <- paste0(nvals)
  
  
  pval.normal_ols4 <- do.call(cbind, pval.normal_ols4)
  pval.normal_ols4 <- data.frame(pval.normal_ols4)
  
  ecdf.ols4 <- sapply(1:ncol(pval.normal_ols4), function(i){
    
    ecdf(pval.normal_ols4[, i])
    
  })
  
  pval.normal_ols4 <- sapply(1:ncol(pval.normal_ols4), function(i){

    mean(pval.normal_ols4[,i]<=0.05)
  })
  names(pval.normal_ols4) <- paste0(nvals)

    
  boxplot(tau)

  nrep <- length(tau[[i]])
  
  # tau <- unlist(sapply(1:length(nvals), function(k){
  #   mean(tau[[k]])
  # }))
  
  ############# analysis
  
  #rmse_iv1 <- matrix(nrow=1, ncol=length(nvals))
  #rmse_iv2 <- matrix(nrow=1, ncol=length(nvals))
  # rmse_iv3 <- matrix(nrow=1, ncol=length(nvals))
  rmse_ols1 <- matrix(nrow=1, ncol=length(nvals))
  rmse_ols2 <- matrix(nrow=1, ncol=length(nvals))
  rmse_ols3 <- matrix(nrow=1, ncol=length(nvals))
  rmse_ols4 <- matrix(nrow=1, ncol=length(nvals))
  # rmse_ols4.simple <- matrix(nrow=1, ncol=length(nvals))
  
  
  #coverage_iv1 <- matrix(nrow=1, ncol=length(nvals))
  #coverage_iv2 <- matrix(nrow=1, ncol=length(nvals))
  # coverage_iv3 <- matrix(nrow=1, ncol=length(nvals))
  # coverage_ols1 <- matrix(nrow=1, ncol=length(nvals))
  # coverage_ols2 <- matrix(nrow=1, ncol=length(nvals))
  # coverage_ols3 <- matrix(nrow=1, ncol=length(nvals))
  # coverage_ols4 <- matrix(nrow=1, ncol=length(nvals))
  # coverage_ols4.simple <- matrix(nrow=1, ncol=length(nvals))
  
  
  #bias_iv1 <- matrix(nrow=1, ncol=length(nvals))
  #bias_iv2 <- matrix(nrow=1, ncol=length(nvals))
  # bias_iv3 <- matrix(nrow=1, ncol=length(nvals))
  bias_ols1 <- matrix(nrow=1, ncol=length(nvals))
  bias_ols2 <- matrix(nrow=1, ncol=length(nvals))
  bias_ols3 <- matrix(nrow=1, ncol=length(nvals))
  bias_ols4 <- matrix(nrow=1, ncol=length(nvals))
  # bias_ols4.simple <- matrix(nrow=1, ncol=length(nvals))
  
  
  #var_iv1 <- matrix(nrow=1, ncol=length(nvals))
  #var_iv2 <- matrix(nrow=1, ncol=length(nvals))
  # var_iv3 <- matrix(nrow=1, ncol=length(nvals))
  var_ols1 <- matrix(nrow=1, ncol=length(nvals))
  var_ols2 <- matrix(nrow=1, ncol=length(nvals))
  var_ols3 <- matrix(nrow=1, ncol=length(nvals))
  var_ols4 <- matrix(nrow=1, ncol=length(nvals))
  # var_ols4.simple <- matrix(nrow=1, ncol=length(nvals))
  
  var.est_ols1 <- matrix(nrow=1, ncol=length(nvals))
  var.est_ols2 <- matrix(nrow=1, ncol=length(nvals))
  var.est_ols3 <- matrix(nrow=1, ncol=length(nvals))
  var.est_ols4 <- matrix(nrow=1, ncol=length(nvals))
  
  var.var.est_ols1 <- matrix(nrow=1, ncol=length(nvals))
  var.var.est_ols2 <- matrix(nrow=1, ncol=length(nvals))
  var.var.est_ols3 <- matrix(nrow=1, ncol=length(nvals))
  var.var.est_ols4 <- matrix(nrow=1, ncol=length(nvals))
  

    #var_ols4.est <- matrix(nrow=1, ncol=length(nvals))
    #var_IV2.est <- matrix(nrow=1, ncol=length(nvals))
    # var_ols4.est.simple <- matrix(nrow=1, ncol=length(nvals))
    #var_iv2.est.freedman <- matrix(nrow=1, ncol=length(nvals))
    #var_iv3.est.freedman <- matrix(nrow=1, ncol=length(nvals))
  

  
  frac_0 <- matrix(nrow=1, ncol=length(nvals))


  num.res <- length(unlist(Res[[1]]))/nrep
  
  
  for (i in 1:length(nvals)){
    
    #rmse_iv1[i] <- sqrt(mean((unname(Res[[i]][which(names(unlist(Res[[i]]))=="est_IV1")])-tau[length(nvals)])^2, na.rm=T))
    #rmse_iv2[i] <- sqrt(mean((unname(Res[[i]][which(names(unlist(Res[[i]]))=="est_IV2")])-tau[length(nvals)])^2, na.rm=T))
    # rmse_iv3[i] <- sqrt(mean((unname(Res[[i]][which(names(unlist(Res[[i]]))=="est_IV3")])-tau[length(nvals)])^2, na.rm=T))
    rmse_ols1[i] <- sqrt(mean((unname(Res[[i]][which(names(unlist(Res[[i]]))=="est_OLS1")])-tau[[i]])^2, na.rm=T))
    rmse_ols2[i] <- sqrt(mean((unname(Res[[i]][which(names(unlist(Res[[i]]))=="est_OLS2")])-tau[[i]])^2, na.rm=T))
    rmse_ols3[i] <- sqrt(mean((unname(Res[[i]][which(names(unlist(Res[[i]]))=="est_OLS3")])-tau[[i]])^2, na.rm=T))
    rmse_ols4[i] <- sqrt(mean((unname(Res[[i]][which(names(unlist(Res[[i]]))=="est_OLS4")])-tau[[i]])^2, na.rm=T))
    # rmse_ols4.simple[i] <- sqrt(mean((unname(Res[[i]][which(names(unlist(Res[[i]]))=="est_OLS4.simple")])-tau[length(nvals)])^2, na.rm=T))
    
    frac_0[i] <- mean(as.numeric(unname(Res[[i]][which(names(unlist(Res[[i]]))=="num0")])))/nvals[i]
    
    # est <- as.numeric(unname(Res[[i]][which(names(unlist(Res[[i]]))=="est_IV1")]))
    # bias_iv1[i] <- mean(est, na.rm=TRUE)- tau[length(nvals)]
    # var_iv1[i] <- var(est, na.rm=TRUE)
    # rm(est)
    # 
    # est <- as.numeric(unname(Res[[i]][which(names(unlist(Res[[i]]))=="est_IV2")]))
    # bias_iv2[i] <- mean(est, na.rm=TRUE)- tau[length(nvals)]
    # var_iv2[i] <- var(est, na.rm=TRUE)
    # rm(est)
    
    # 
    # est <- as.numeric(unname(Res[[i]][which(names(unlist(Res[[i]]))=="est_IV3")]))
    # bias_iv3[i] <- mean(est, na.rm=TRUE)- tau[length(nvals)]
    # var_iv3[i] <- var(est, na.rm=TRUE)
    
    est <- as.numeric(unname(Res[[i]][which(names(unlist(Res[[i]]))=="est_OLS1")]))
    bias_ols1[i] <- mean(est- tau[[i]],na.rm=TRUE)
    var_ols1[i] <- var(est, na.rm=TRUE)
  
    var.est_ols1[i] <-  1/nvals[i]*mean(Res[[i]][which(names(unlist(Res[[i]]))=="var_OLS1.est")])
    var.var.est_ols1[i] <-  var(1/nvals[i]*Res[[i]][which(names(unlist(Res[[i]]))=="var_OLS1.est")])
    
    
    rm(est)
    
    est <- as.numeric(unname(Res[[i]][which(names(unlist(Res[[i]]))=="est_OLS2")]))
    bias_ols2[i] <- mean(est- tau[[i]],na.rm=TRUE)
    var_ols2[i] <- var(est, na.rm=TRUE)
    #pval.normal_ols2[i] <- shapiro.test(sqrt(nvals[i])*(est-tau[length(nvals)]) )$p.value
    
    var.est_ols2[i] <-  1/nvals[i]*mean(Res[[i]][which(names(unlist(Res[[i]]))=="var_OLS2.est")])
    var.var.est_ols2[i] <-  var(1/nvals[i]*Res[[i]][which(names(unlist(Res[[i]]))=="var_OLS2.est")])
    
    
    rm(est)
    
    est <- as.numeric(unname(Res[[i]][which(names(unlist(Res[[i]]))=="est_OLS3")]))
    bias_ols3[i] <- mean(est- tau[[i]],na.rm=TRUE)
    var_ols3[i] <- var(est, na.rm=TRUE)
    
    var.est_ols3[i] <- 1/nvals[i]*mean(Res[[i]][which(names(unlist(Res[[i]]))=="var_OLS3.est")])
    
    var.var.est_ols3[i] <-  var(1/nvals[i]*Res[[i]][which(names(unlist(Res[[i]]))=="var_OLS3.est")])
    
    
    #pval.normal_ols3[i] <- shapiro.test(sqrt(nvals[i])*(est-tau[length(nvals)]) )$p.value
    rm(est)
    
    est <- as.numeric(unname(Res[[i]][which(names(unlist(Res[[i]]))=="est_OLS4")]))
    bias_ols4[i] <- mean(est- tau[[i]],na.rm=TRUE)
    var_ols4[i] <- var(est, na.rm=TRUE)
    
    var.est_ols4[i] <-  1/nvals[i]*mean(Res[[i]][which(names(unlist(Res[[i]]))=="var_OLS4.est")])
    
    var.var.est_ols4[i] <-  var(1/nvals[i]*Res[[i]][which(names(unlist(Res[[i]]))=="var_OLS4.est")])
    
    rm(est)
  
  



    if(typeofgraph =="WS_growing"){
      
        if(growth.rate.WS == 1/2){
          ylim.rmse <- c(0,20)
          ylim.bias <- c(-1.5,1.5)
          ylim.var <- c(-8,7)
        }else if(growth.rate.WS == 1/4){
            ylim.rmse <- c(0,2.5)
            ylim.bias <- c(-1.5,1)
            ylim.var <- c(-8,1)
        }else if(growth.rate.WS %in% c(1/15, 1/17)){
          ylim.rmse <- c(0,1.6)
          ylim.bias <- c(-1.5,1.2)
          ylim.var <- c(-7,0)
        }
    
    }else if(typeofgraph == "family"){
      ylim.rmse <- c(0,1.8)
      ylim.bias <- c(-1.9,1.3)
      ylim.var <- c(-7,0)
    }else if(typeofgraph=="rand_npfix_growing"){
      ylim.rmse <- c(0,1.55)
      ylim.bias <- c(-1.5,1.3)
      ylim.var <- c(-7,0)
      
    }
    

  
  if(typeofgraph =="WS_growing"){
    filename.save <- paste0("growth.rate.WS", round(growth.rate.WS,3),typeofgraph,"_pi",pai, "_eta",eta )
  }else if (typeofgraph=="family"){
    filename.save <- paste0(typeofgraph,"_pi",pai, "_eta",eta )
  }else if(typeofgraph=="rand_npfix"){
    filename.save <- paste0(typeofgraph,"_const", const, "_pi",pai, "_eta",eta )
    
  }else if (typeofgraph =="3dlatt"){
    filename.save <- paste0(typeofgraph,"_pi",pai, "_eta",eta )
    
  }else if(typeofgraph =="2dlatt"){
    filename.save <- paste0(typeofgraph,"_pi",pai, "_eta",eta )
    
  }
  
  pdf(file = paste0(filename.save , ".pdf"), width = 11, height = 4.5)
  
  par(oma = c(4, 1, 1, 1))
  
  #par(mar=c(1.5,1.5,1,1))
  
  par(mfrow=c(1,3))
  
  #### RMSE
  
  # plot(rmse_iv1[1,]~nvals, ylim=ylim.rmse, type = "l", main =TeX("RMSE"), xlab=TeX("$N$"),
  #      col= "snow4",lwd=2, ylab="", cex.lab = 1.2, cex.main=1.2, frame.plot = FALSE, axes=F)

  
  plot(rmse_ols1[1,]~nvals, ylim=c(0,5), type = "l", main =paste0("RMSE"), xlab=TeX("$N$"),
       col= "tomato",lwd=2, ylab="", cex.lab = 1.2, cex.main=1.2, frame.plot = FALSE, axes=F)
  
    box(bty="l")
  # axis(1, at=nvals,las=2, labels=nvals, 
  #      las = 2, cex.lab = 1.2,cex.axis = 1.25)
  # 
  
  axis(1, at=seq(300, 5000, by=1000),las=2, labels=seq(300, 5000, by=1000),
       las = 2, cex.lab = 1.2,cex.axis = 1)
  #lines(rmse_iv2[1,]~nvals, col="red4",lwd=2, lty=2)
  #lines(rmse_ols1[1,]~nvals, col="tomato",lwd=2, lty=3)
  lines(rmse_ols2[1,]~nvals, col="skyblue3",lwd=2, lty=4)
  lines(rmse_ols3[1,]~nvals, col="palegreen4",lwd=2, lty=5)
  lines(rmse_ols4[1,]~nvals, col="hotpink3",lwd=2, lty=6)
  # lines(rmse_ols4.simple[1,]~nvals, col="darkblue",lwd=2, lty=7)
  # lines(rmse_iv3[1,]~nvals, col="yellow2",lwd=2, lty=8)
  
  
  abline(a=0,b=0,lty=3)
  axis(2, cex.axis=1)
  #axis(1, cex.axis=1.25)
  
 
  
  ### Bias
  
  # plot(bias_iv1[1,]~nvals,
  #      ylim=ylim.bias, type = "l", main =TeX("Bias"), xlab=TeX("$N$"),
  #      col= "snow4",lwd=2, ylab="",xaxt="n",cex.lab = 1.2, cex.main=1.2, frame.plot = FALSE, axes=F)
  
  plot(bias_ols1[1,]~nvals,
       ylim=c(0,5), type = "l", main =paste0("Bias"), xlab=TeX("$N$"),
       col= "tomato",lwd=2, ylab="",xaxt="n",cex.lab = 1.2, cex.main=1.2, frame.plot = FALSE, axes=F)
  
  box(bty="l")
  # axis(1, at=nvals,las=2, labels=nvals, 
  #      las = 2, cex.lab = 1.2,cex.axis = 1.25)
  
  axis(1, at=seq(300, 5000, by=1000),las=2, labels=seq(300, 5000, by=1000),
       las = 2, cex.lab = 1.2,cex.axis = 1)
  
  #lines(bias_iv2[1,]~nvals, col="red4",lwd=2, lty=2)
  #lines(bias_ols1[1,]~nvals, col="tomato",lwd=2, lty=3)
  lines(bias_ols2[1,]~nvals, col="skyblue3",lwd=2, lty=4)
  lines(bias_ols3[1,]~nvals, col="palegreen4",lwd=2, lty=5)
  lines(bias_ols4[1,]~nvals, col="hotpink3",lwd=2, lty=6)
  # lines(bias_ols4.simple[1,]~nvals, col="darkblue",lwd=2, lty=7)
  # lines(bias_iv3[1,]~nvals, col="yellow2",lwd=2, lty=8)
  
  
  abline(a=0,b=0,lty=3)
  axis(2, cex.axis=1)
  #axis(1, cex.axis=1.25)
  
  ### log(Var)
  # plot(log(var_iv1[1,])~log(nvals),
  #      ylim=ylim.var, type = "l", main =TeX("$\\log(Variance)$"), xlab=TeX("$\\log(N)$"),
  #      col= "snow4",lwd=2, ylab="", cex.lab = 1.2, cex.main=1.2, frame.plot = FALSE, axes=F)
  
  plot(log(var_ols1[1,])~log(nvals),
       ylim=c(0,5), type = "l", main =paste0("log(Variance)"), xlab=TeX("$log(N)$"),
       col= "tomato",lwd=2, ylab="", cex.lab = 1.2, cex.main=1.2, frame.plot = FALSE, axes=F)
  
  box(bty="l")
  # axis(1, at=log(nvals),las=2, labels=log(nvals), 
  #      las = 2, cex.lab = 1.2,cex.axis = 1.25)
  axis(1, at=c(6,6.5,7,7.5,8,8.5,9),las=2, labels=c(6, 6.5,7,7.5,8,8.5,9),
       las = 2, cex.lab = 1.2,cex.axis = 1)
  
  #lines(log(var_iv2[1,])~log(nvals), col="red4",lwd=2, lty=2)
  #lines(log(var_ols1[1,])~log(nvals), col="tomato",lwd=2, lty=3)
  lines(log(var_ols2[1,])~log(nvals), col="skyblue3",lwd=2, lty=4)
  lines(log(var_ols3[1,])~log(nvals), col="palegreen4",lwd=2, lty=5)
  lines(log(var_ols4[1,])~log(nvals), col="hotpink3",lwd=2, lty=6)
  
  
  lines(log(var.est_ols1[1,])~log(nvals), col="tomato")
  lines(log(var.est_ols2[1,])~log(nvals), col="skyblue3")
  lines(log(var.est_ols3[1,])~log(nvals), col="palegreen4")
  lines(log(var.est_ols4[1,])~log(nvals), col="hotpink3")
  
  
  # lines(log(var_ols4.simple[1,])~log(nvals), col="darkblue",lwd=2, lty=7)
  #lines(log(var_ols4.est[1,])~log(nvals), col="hotpink3",lwd=2, lty=3)
  #lines(log(var_IV2.est[1,])~log(nvals), col="black",lwd=2, lty=6)
  #lines(log(var_ols4.est.simple[1,])~log(nvals), col="darkblue",lwd=2, lty=3)
  #lines(log(var_iv2.est.freedman[1,])~log(nvals), col="red4",lwd=2, lty=3)
  # lines(log(var_iv3[1,])~log(nvals), col="yellow2",lwd=2, lty=8)
  # lines(log(var_iv3.est.freedman[1,])~log(nvals), col="yellow2",lwd=2, lty=3)
  
  #summary(lm(log(var_ols4[1,])~log(nvals)))
  # shapiro.test(est)
  
  abline(a=0,b=0,lty=3)
  axis(2, cex.axis=1)
  #axis(1, cex.axis=1.25)
  
  ######## asymptotic normality
  
# 
#   plot(pval.normal_ols1~nvals,
#        ylim=c(0,0.3), type = "l", main =paste0(TeX("Normality: "),"P[p-value <= 0.05]"), xlab=TeX("$N$"),
#        col= "tomato",lwd=2, ylab="", cex.lab = 1.2, cex.main=1.2, frame.plot = FALSE, axes=F)
#   
#   box(bty="l")
#  
# 
#   lines(pval.normal_ols2~nvals, col=alpha("skyblue3",0.4),lwd=2, lty=4)
#   lines(pval.normal_ols3~nvals, col="palegreen4",lwd=2, lty=5)
#   lines(pval.normal_ols4~nvals, col="hotpink3",lwd=2, lty=6)
#  
#   abline(h=0.05,lty=3)
#   axis(2, cex.axis=1)
#   axis(1, at=seq(300, 5000, by=1000),las=2, labels=seq(300, 5000, by=1000),
#        las = 2, cex.lab = 1.2,cex.axis = 1)

  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  
  legtext <- c(TeX("OLS1($W$)"),TeX("OLS2($W$, $C_2$)"),TeX("OLS3($W$, $X$, $O$)"),
               TeX("OLS4($W$, $X$, $O$,$C_2$)"))
               # TeX("IV1(C1,C2)"), 
               # TeX("IV2(C1,C2,X)"))
  
  xcoords <- c(1, 2, 3, 4)#,5,6)
  secondvector <- (1:length(legtext))-1
  textwidths <- xcoords/secondvector # this works for all but the first element
  textwidths[1] <- 0 
  
  legend("bottom", legtext, 
         #xpd = TRUE,
         horiz = TRUE, 
         #text.width=textwidths, 
         text.width=c(0.16, 0.18, 0.17,0.3),# 0.19),#, 0.272, 0.24),
         inset = c(0, 0), 
         x.intersp=1, 
         xjust=0,
         yjust =0, bty = "n", 
         lty = c(1, 4,5,6),
         col = c("tomato","skyblue3", "palegreen4","hotpink3"),#, "snow4", "red4"), 
  cex = 1.1, lwd=2)
  

  #  horiz = FALSE, inset = c(0.5,-0.01),"bottom"x.intersp=0.1,
  
  dev.off()
  
  
  setwd("..")
  
  # 
  # library(EnvStats)
  # library(scales)
  # col = c("tomato","skyblue3", "palegreen4","hotpink3")
  # #ylab="Fn(x)",ylab="x",
  # 
  # plot(unlist(ecdf.ols1[[1]]),col = "tomato",  lwd=0.01)
  # lines(ecdf.ols2[[1]],col="skyblue3", lwd=0.01)
  # lines(ecdf.ols3[[1]],col="palegreen4", lwd=0.01)
  # lines(ecdf.ols4[[1]],col="hotpink3", lwd=0.01)
  # 
  # plot(unlist(ecdf.ols1[[4]]),col = "tomato",  lwd=0.01)
  # lines(ecdf.ols2[[4]],col="skyblue3", lwd=0.01)
  # lines(ecdf.ols3[[4]],col="palegreen4", lwd=0.01)
  # lines(ecdf.ols4[[4]],col="hotpink3", lwd=0.01)
  # 
  # qqplot(pval.normal_ols4[,5])
  # 
  # 
  # par(oma = c(4, 1, 1, 1))
  # 
  # #par(mar=c(1.5,1.5,1,1))
  # 
  # par(mfrow=c(1,2))
  # 
  # qqPlot(pval.normal_ols4[,1], distribution = "unif", 
  #        param.list = list(min=0, max=1), lty=2
  #        )
  # 
  # 
  # qqPlot(pval.normal_ols4[,2], distribution = "unif", 
  #        param.list = list(min=0, max=1), lty=2
  # )
  # 
  # qqPlot(pval.normal_ols4[,3], distribution = "unif", 
  #        param.list = list(min=0, max=1), lty=2
  # )
  # 
  # qqPlot(pval.normal_ols4[,4], distribution = "unif", 
  #        param.list = list(min=0, max=1), lty=2
  # )
  # 
  # qqPlot(pval.normal_ols4[,1], distribution = "unif", 
  #        param.list = list(min=0, max=1), lty=2
  # )
  # 
  # sapply(1:1000, function(k){
  #   rand.samp <- runif(500, min=0, max=1)
  #   
  #   qqline(rand.samp, distribution = qunif, 
  #           col=alpha("blue", 0.1))
  #   
  #   
  # })
  # 
  # 
  # qqPlot(pval.normal_ols4[,5], distribution = "unif", 
  #        param.list = list(min=0, max=1), add=T
  # )
  # 
  # qq(pval.normal_ols4[,5],distribution = qunif, col="red")
  # 
  # 
  # 
  ###################################
  
  
  pdf(file = paste0("ecdf",typeofgraph , ".pdf"), width = 10, height = 4.5)
  
  par(mfrow=c(1,2))
  
  
  
  for( i in c(1,5)){
    plot(ecdf.ols4[[i]], col="tomato", main =paste0("N = ",nvals[[i]],
                                                    ": ecdf of p-values"), 
         cex.lab = 1, cex.main=1, frame.plot = FALSE, axes=F)
    
    
    box(bty="l")
    
    
    axis(1, at=seq(0, 1, by=0.2),las=2, labels=seq(0, 1, by=0.2),
         las = 2, cex.lab = 1,cex.axis = 1)
    
    axis(2, at=seq(0, 1, by=0.2),las=2, labels=seq(0, 1, by=0.2),
         cex.lab = 1,cex.axis = 1)
    
    for(k in 1:100){
      rand.samp <- runif(500, min=0, max=1)
      
      lines(ecdf(rand.samp), 
            col="lightgray")
      
    }
    
    lines(ecdf.ols4[[i]], col="tomato")
    
    

    
  }
  

  
  dev.off()
  
  
  # plot(var_iv1[1,]~nvals,
  #      ylim=ylim.var, type = "l", main =TeX("Standard Deviation"), xlab=TeX("$n$"),
  #      col= "snow4",lwd=2, ylab="", cex.lab = 1.5, cex.main=1.5, frame.plot = FALSE, axes=F)
  # box(bty="l")
  # lines(var_iv2[1,]~nvals, col="red4",lwd=2, lty=2)
  # lines(var_ols1[1,]~nvals, col="tomato",lwd=2, lty=3)
  # lines(var_ols2[1,]~nvals, col="skyblue3",lwd=2, lty=4)
  # lines(var_ols3[1,]~nvals, col="palegreen4",lwd=2, lty=5)
  # lines(var_ols4[1,]~nvals, col="hotpink3",lwd=2, lty=6)
  # lines(var_ols4.simple[1,]~nvals, col="darkblue",lwd=2, lty=6)
  # lines(sqrt(var_OLS4.est[1,])~nvals, col="black",lwd=2, lty=6)
  # lines(sqrt(var_IV2.est[1,])~nvals, col="black",lwd=2, lty=6)
  # lines(sqrt(var_OLS4.est.simple[1,])~nvals, col="darkblue",lwd=2, lty=6)
  # lines(sqrt(var_IV2.est.freedman[1,])~nvals, col="yellow",lwd=2, lty=6)
  
  # if(type %in% c("CiaufWj", "CiaufYj")){
  #   lines(var_iv4[1,]~nvals, col="darkblue",lwd=2, lty=7)
  # }
  #if(type=="multipleC"){
  #lines(var_ols4.simple[1,]~nvals, col="darkblue",lwd=2, lty=7)
  #}
  # 
  # 
  # abline(a=0,b=0,lty=3)
  # axis(2, cex.axis=1.25)
  # axis(1, cex.axis=1.25)
  
  
}


generate_plot(typeofgraph=typeofgraph,
              nvals=nvals,
              dateused=dateused, 
              parts=parts,
              use.only.feat1 = use.only.feat1,
              growth.rate.WS = growth.rate.WS 
              )



###########



