library(latex2exp)
library(ggplot2)
library(ggpubr)

setwd("~/GitHub/InvarianceCode")

typeofgraph <- "rand_npfix" #  "rand_pfix","rand_npfix","rand_npfix_growing","family","2dlatt"

nvals <- c(300, 600, 1200,2400,4800)


filename.save <- paste0("plot_",typeofgraph)

tau <- list()
Res <- list()



for(i in 1:length(nvals)){
  
  nval <- nvals[i]
  
  filename <- paste0("model_", typeofgraph,"/estimation_",typeofgraph,"_nval",nval,".Rda")
  
  
  load(filename)
  
  
  Results <- unlist(Results, recursive=FALSE)
  
  
  taus <- unlist(do.call(rbind, Results)[,"tau"])
  tau[[i]] <- taus
  Res[[i]] <- unlist(Results)
  
  
  rm(Results)
  print(paste0("nval: ",nvals[i]))
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

pdf(file = paste0(filename.save , ".pdf"), width = 11, height = 4.5)

estimator <- as.factor(c(rep(1,length(nvals)),rep(2,length(nvals)),rep(3,length(nvals)),rep(4,length(nvals))))

data_rmse <- data.frame(RMSE=c(rmse_ols1,rmse_ols2,rmse_ols3,rmse_ols4),
                   estimator=estimator,n=rep(nvals,4))
data_bias <- data.frame(Bias=c(bias_ols1,bias_ols2,bias_ols3,bias_ols4),
                        estimator=estimator,n=rep(nvals,4))
data_var <- data.frame(Variance=c(var_ols1,var_ols2,var_ols3,var_ols4,var.est_ols4),
                        estimator=as.factor(c(estimator,rep(5,length(nvals)))),n=rep(nvals,5))

g_rmse <- ggplot(data = data_rmse, aes(x = n, y = RMSE,group=estimator)) +  ylim(0, 1) + xlim(-500,5500) +
  geom_point(aes(shape = estimator,col=estimator)) + theme(legend.position = "none") + geom_line(aes(linetype=estimator,col=estimator))

g_bias <- ggplot(data = data_bias, aes(x = n, y = Bias)) +  ylim(-1, 1) + xlim(-500,5500) + 
  geom_point(aes(shape = estimator,col=estimator)) + theme(legend.position = "none") + geom_line(aes(linetype=estimator,col=estimator))

g_var <- ggplot(data = data_var, aes(x = log(n), y = log(Variance))) +  ylim(-8, -1) + xlim(5.2,8.8) +
  geom_point(aes(shape = estimator,col=estimator)) + theme(legend.position = "none") + geom_line(aes(linetype=estimator,col=estimator))



g <- ggarrange(g_rmse,g_bias,g_var,ncol=3,common.legend = TRUE, legend="bottom")

annotate_figure(g, top = text_grob("Erd\U0151s-R\U00E9nyi networks I(N,10/N)",
                                      color = "black", face = "bold", size = 14))
 g
  dev.off()





