library(latex2exp)
library(ggplot2)
library(ggpubr)
library(xtable)


typeofgraphs <-  c("rand_pfix","rand_npfix","rand_npfix_growing","family","2dlatt")

nvals <- c(300,600,1200,2400,4800)

rmse_var.est <- matrix(numeric(25),ncol=5)

tau <- list()
Res <- list()

for(j in 1:length(typeofgraphs)){
  
  typeofgraph <- typeofgraphs[j]

for(i in 1:length(nvals)){
  
  nval <- nvals[i]
  
  filename <- paste0("SimulationStudy/results_data/model_", typeofgraph,"/estimation_",typeofgraph,"_nval",nval,".Rda")
  
  
  load(filename)
  
  
  Results <- unlist(Results, recursive=FALSE)
  
  
  taus <- unlist(do.call(rbind, Results)[,"tau"])
  tau[[i]] <- taus
  Res[[i]] <- unlist(Results)
  
  
  rm(Results)
  print(paste0("nval: ",nvals[i]))
}

if(typeofgraph=="2dlatt"){
  nvals <- c(289, 576, 1225, 2401, 4761)
}

rmse_ols1 <- matrix(nrow=1, ncol=length(nvals))
rmse_ols2 <- matrix(nrow=1, ncol=length(nvals))
rmse_ols3 <- matrix(nrow=1, ncol=length(nvals))
rmse_ols4 <- matrix(nrow=1, ncol=length(nvals))

bias_ols1 <- matrix(nrow=1, ncol=length(nvals))
bias_ols2 <- matrix(nrow=1, ncol=length(nvals))
bias_ols3 <- matrix(nrow=1, ncol=length(nvals))
bias_ols4 <- matrix(nrow=1, ncol=length(nvals))

var_ols1 <- matrix(nrow=1, ncol=length(nvals))
var_ols2 <- matrix(nrow=1, ncol=length(nvals))
var_ols3 <- matrix(nrow=1, ncol=length(nvals))
var_ols4 <- matrix(nrow=1, ncol=length(nvals))

var.est_ols1 <- matrix(nrow=1, ncol=length(nvals))
var.est_ols2 <- matrix(nrow=1, ncol=length(nvals))
var.est_ols3 <- matrix(nrow=1, ncol=length(nvals))
var.est_ols4 <- matrix(nrow=1, ncol=length(nvals))

var.var.est_ols1 <- matrix(nrow=1, ncol=length(nvals))
var.var.est_ols2 <- matrix(nrow=1, ncol=length(nvals))
var.var.est_ols3 <- matrix(nrow=1, ncol=length(nvals))
var.var.est_ols4 <- matrix(nrow=1, ncol=length(nvals))


for (i in 1:length(nvals)){
  rmse_ols1[i] <- sqrt(mean((unname(Res[[i]][which(names(unlist(Res[[i]]))=="est_OLS1")])-tau[[i]])^2, na.rm=T))
  rmse_ols2[i] <- sqrt(mean((unname(Res[[i]][which(names(unlist(Res[[i]]))=="est_OLS2")])-tau[[i]])^2, na.rm=T))
  rmse_ols3[i] <- sqrt(mean((unname(Res[[i]][which(names(unlist(Res[[i]]))=="est_OLS3")])-tau[[i]])^2, na.rm=T))
  rmse_ols4[i] <- sqrt(mean((unname(Res[[i]][which(names(unlist(Res[[i]]))=="est_OLS4")])-tau[[i]])^2, na.rm=T))
  
  # frac_0[i] <- mean(as.numeric(unname(Res[[i]][which(names(unlist(Res[[i]]))=="num0")])))/nvals[i]
  
  est <- as.numeric(unname(Res[[i]][which(names(unlist(Res[[i]]))=="est_OLS1")]))
  bias_ols1[i] <- mean(est- tau[[i]],na.rm=TRUE)
  var_ols1[i] <- var(est, na.rm=TRUE)
  var.est_ols1[i] <-  1/nvals[i]*mean(Res[[i]][which(names(unlist(Res[[i]]))=="var_OLS1.est")])
  var.var.est_ols1[i] <-  var(1/nvals[i]*Res[[i]][which(names(unlist(Res[[i]]))=="var_OLS1.est")])
  rm(est)
  
  est <- as.numeric(unname(Res[[i]][which(names(unlist(Res[[i]]))=="est_OLS2")]))
  bias_ols2[i] <- mean(est- tau[[i]],na.rm=TRUE)
  var_ols2[i] <- var(est, na.rm=TRUE)
  var.est_ols2[i] <-  1/nvals[i]*mean(Res[[i]][which(names(unlist(Res[[i]]))=="var_OLS2.est")])
  var.var.est_ols2[i] <-  var(1/nvals[i]*Res[[i]][which(names(unlist(Res[[i]]))=="var_OLS2.est")])
  rm(est)
  
  est <- as.numeric(unname(Res[[i]][which(names(unlist(Res[[i]]))=="est_OLS3")]))
  bias_ols3[i] <- mean(est- tau[[i]],na.rm=TRUE)
  var_ols3[i] <- var(est, na.rm=TRUE)
  var.est_ols3[i] <- 1/nvals[i]*mean(Res[[i]][which(names(unlist(Res[[i]]))=="var_OLS3.est")])
  var.var.est_ols3[i] <-  var(1/nvals[i]*Res[[i]][which(names(unlist(Res[[i]]))=="var_OLS3.est")])
  rm(est)
  
  est <- as.numeric(unname(Res[[i]][which(names(unlist(Res[[i]]))=="est_OLS4")]))
  bias_ols4[i] <- mean(est- tau[[i]],na.rm=TRUE)
  var_ols4[i] <- var(est, na.rm=TRUE)
  var.est_ols4[i] <-  1/nvals[i]*mean(Res[[i]][which(names(unlist(Res[[i]]))=="var_OLS4.est")])
  var.var.est_ols4[i] <-  var(1/nvals[i]*Res[[i]][which(names(unlist(Res[[i]]))=="var_OLS4.est")])
  rm(est)}



for(i in 1:length(nvals)){
  
  rmse_var.est[i,j] <- mean((Res[[i]][which(names(unlist(Res[[i]]))=="var_OLS4.est")]-nvals[i]*var_ols4[i])^2)
  
}
}

xtable(rmse_var.est[,c(2,4,5)])
