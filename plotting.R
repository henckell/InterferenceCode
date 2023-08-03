library(latex2exp)

setwd("~/GitHub/InvarianceCode")

typeofgraph <- "rand_npfix_growing" #  "rand_pfix","rand_npfix","rand_npfix_growing","family","2dlatt"

nvals <- c(300, 600, 1200,2400,4800)


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

par(oma = c(4, 1, 1, 1))

#par(mar=c(1.5,1.5,1,1))

par(mfrow=c(1,3))

plot(rmse_ols1[1,]~nvals, ylim=c(0,1), type = "l", main =paste0("RMSE"), xlab=TeX("$N$"),
     col= "tomato",lwd=2, ylab="", cex.lab = 1.2, cex.main=1.2, frame.plot = FALSE, axes=F)

box(bty="l")
axis(1, at=seq(300, 5000, by=1000),las=2, labels=seq(300, 5000, by=1000),
     las = 2, cex.lab = 1.2,cex.axis = 1)
lines(rmse_ols2[1,]~nvals, col="skyblue3",lwd=2, lty=4)
lines(rmse_ols3[1,]~nvals, col="palegreen4",lwd=2, lty=5)
lines(rmse_ols4[1,]~nvals, col="hotpink3",lwd=2, lty=6)
abline(a=0,b=0,lty=3)
axis(2, cex.axis=1)



### Bias
plot(bias_ols1[1,]~nvals,
     ylim=c(-2,1), type = "l", main =paste0("Bias"), xlab=TeX("$N$"),
     col= "tomato",lwd=2, ylab="",xaxt="n",cex.lab = 1.2, cex.main=1.2, frame.plot = FALSE, axes=F)

box(bty="l")

axis(1, at=seq(300, 5000, by=1000),las=2, labels=seq(300, 5000, by=1000),
     las = 2, cex.lab = 1.2,cex.axis = 1)

lines(bias_ols2[1,]~nvals, col="skyblue3",lwd=2, lty=4)
lines(bias_ols3[1,]~nvals, col="palegreen4",lwd=2, lty=5)
lines(bias_ols4[1,]~nvals, col="hotpink3",lwd=2, lty=6)


abline(a=0,b=0,lty=3)
axis(2, cex.axis=1)

### log(Var)
plot(log(var_ols1[1,])~log(nvals),
     ylim=c(-8,-2), type = "l", main =paste0("log(Variance)"), xlab=TeX("$log(N)$"),
     col= "tomato",lwd=2, ylab="", cex.lab = 1.2, cex.main=1.2, frame.plot = FALSE, axes=F)

box(bty="l")
axis(1, at=c(6,6.5,7,7.5,8,8.5,9),las=2, labels=c(6, 6.5,7,7.5,8,8.5,9),
     las = 2, cex.lab = 1.2,cex.axis = 1)

lines(log(var_ols2[1,])~log(nvals), col="skyblue3",lwd=2, lty=4)
lines(log(var_ols3[1,])~log(nvals), col="palegreen4",lwd=2, lty=5)
lines(log(var_ols4[1,])~log(nvals), col="hotpink3",lwd=2, lty=6)


lines(log(var.est_ols1[1,])~log(nvals), col="tomato")
lines(log(var.est_ols2[1,])~log(nvals), col="skyblue3")
lines(log(var.est_ols3[1,])~log(nvals), col="palegreen4")
lines(log(var.est_ols4[1,])~log(nvals), col="hotpink3")

abline(a=0,b=0,lty=3)
axis(2, cex.axis=1)

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


# setwd("..")

