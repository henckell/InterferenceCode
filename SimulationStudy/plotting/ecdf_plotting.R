################################### 

library(latex2exp)
library(scales)
library(ggplot2)
library(ggpubr)

nvals <- c(300, 600, 1200,2400, 4800)#,9600)

setwd("~/GitHub/InvarianceCode/")

typeofgraph <-  "2dlatt"

# c("rand_pfix","rand_npfix","rand_npfix_growing","family","2dlatt")



filename.save <- paste0("plot_",typeofgraph)

tau <- list()
Res <- list()

pval.normal_ols1 <- list()
pval.normal_ols2 <- list()
pval.normal_ols3 <- list()
pval.normal_ols4 <- list()
  

for(i in 1:length(nvals)){
  nval <- nvals[i]
  
  filename <- paste0("results_data/model_", typeofgraph,"/estimation_",typeofgraph,"_nval",nval,".Rda")
  
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
ecdf.ols1 <- sapply(1:ncol(pval.normal_ols1), function(i){ecdf(pval.normal_ols1[, i])})
names(pval.normal_ols1)<- paste0(nvals)

pval.normal_ols2 <- do.call(cbind, pval.normal_ols2)
pval.normal_ols2 <- data.frame(pval.normal_ols2)
ecdf.ols2 <- sapply(1:ncol(pval.normal_ols2), function(i){ecdf(pval.normal_ols2[, i])})
names(pval.normal_ols2)<- paste0(nvals)

pval.normal_ols3 <- do.call(cbind, pval.normal_ols3)
pval.normal_ols3 <- data.frame(pval.normal_ols3)
ecdf.ols3 <- sapply(1:ncol(pval.normal_ols3), function(i){ecdf(pval.normal_ols3[, i])})
names(pval.normal_ols3)<- paste0(nvals)

pval.normal_ols4 <- do.call(cbind, pval.normal_ols4)
pval.normal_ols4 <- data.frame(pval.normal_ols4)
ecdf.ols4 <- sapply(1:ncol(pval.normal_ols4), function(i){ecdf(pval.normal_ols4[, i])})
names(pval.normal_ols4)<- paste0(nvals)

# boxplot(tau)

pdf(file = paste0("ecdf",typeofgraph , ".pdf"), width = 10, height = 4.5)

par(mfrow=c(1,2))


if(typeofgraph=="2dlatt"){
  nvals <- c(289, 576, 1225, 2401, 4761)
}

g_full <- lapply(seq_len(ncol(pval.normal_ols4)), FUN = function(i) {
  
  g <- ggplot(data = pval.normal_ols4,aes(pval.normal_ols4[,i]))+stat_ecdf(geom = "step",col="black") +
    ggtitle(paste0("N = ",nvals[i],": ecdf of p-values")) +  xlab("p-value") + ylab(TeX("$F_n(x)$"))
  
  for(k in 1:100){
    rand.samp <- runif(50, min=0, max=1)

    g <- g + stat_ecdf(data=data.frame(rand.samp),aes(rand.samp),
                       geom = "step",col="lightgray")

  }

  g <- g + stat_ecdf(geom = "step",col="black",size=2)
  
  return(g)
})

ggarrange(g_full[[1]],g_full[[5]])

dev.off()
