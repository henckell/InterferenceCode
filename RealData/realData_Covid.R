
library(dplyr)
library(Matrix)
library(ggplot2)
library(fastDummies)

# Comments: 
# data are smoothed weekly
# have 24 weeks (X.oneweek) for each of the 26 cantons (X.Canton_3)

source("~/GitHub/InterferenceCode/Estimator/estimator.R")
source("~/GitHub/InterferenceCode/helpers/helpers_realdata.R")

data <- read.csv("RealData/Data/data.casegrowth.csv")
  
# adjacency matrix cantons

A <- create.A.init()

A_large <- kronecker(A,diag(24))
  
treat <- data$W
  
C <- data[, -which(colnames(data) %in% c("Y","W"))]
Y <- data$Y
  
adj <-  c("X.casegrowthlag_we"
          , "X.percentage_age_we","X.Density_we","X.population_we"
           ,"X.restGatherings_we","X.cancEvents_we"
            ,"X.workClosing2a_we","X.sre000d0_we"
            ,"X.tre200d0_we","X.ure200d0_we","X.ferien_we")

factor <- c()  

############################# OLS1  ###########################################################
  
res_OLS1 <- estimator(1,0,treat,Y,C,c(),list(),A_large)
est_OLS1 <- res_OLS1$hat_tau
var_OLS1.est <-  res_OLS1$hat_var_tau
CI_OLS1 <- c(est_OLS1 - qnorm(1- 0.05/2)*sqrt(var_OLS1.est/(24*26)), 
                   est_OLS1 + qnorm(1- 0.05/2)*sqrt(var_OLS1.est/(24*26)))
      
      
  ############################# OLS2 ###########################################################
  
res_OLS2 <- estimator(1,0,treat,Y,C,adj,list(),A_large,factor=factor)
est_OLS2 <- res_OLS2$hat_tau
var_OLS2.est <- res_OLS2$hat_var_tau
CI_OLS2 <- c(est_OLS2 - qnorm(1- 0.05/2)*sqrt(var_OLS2.est/(24*26)), 
                     est_OLS2 + qnorm(1- 0.05/2)*sqrt(var_OLS2.est/(24*26)))
  
  
  ############################# OLS3 ###########################################################
  
res_OLS3 <- estimator(1,0,treat,Y,C,c(),list(feat_X1),A_large)
est_OLS3 <- res_OLS3$hat_tau
var_OLS3.est <- res_OLS3$hat_var_tau
CI_OLS3 <- c(est_OLS3 - qnorm(1- 0.05/2)*sqrt(var_OLS3.est/(24*26)), 
                     est_OLS3 + qnorm(1- 0.05/2)*sqrt(var_OLS3.est/(24*26)))
  
  ############################# OLS4  ###########################################################
  
res_OLS4 <- estimator(pai=1,eta=0,W=treat,Y=Y,C=C,adj=adj,feat_functions=list(feat_X1),A=A_large,factor=factor)
est_OLS4 <- res_OLS4$hat_tau
var_OLS4.est <- res_OLS4$hat_var_tau
CI_OLS4 <- c(est_OLS4 - qnorm(1- 0.05/2)*sqrt(var_OLS4.est/(24*26)), 
                     est_OLS4 + qnorm(1- 0.05/2)*sqrt(var_OLS4.est/(24*26)))
      
      CIs <- rbind(CI_OLS1, CI_OLS2,CI_OLS3, CI_OLS4)
      ests <-  rbind(est_OLS1, est_OLS2, est_OLS3, est_OLS4)
      res <- cbind(ests, CIs)



      Estimator <- c("1","2","3","4")
      
      res <- as.data.frame(res)
      
      res$Estimator <- Estimator
      
      par(mfrow=c(1,1))
      
      pdf("RealData.pdf", width = 5, height = 5)
      
      g <- ggplot(data=res,aes(x=Estimator,y=V1,ymin=V2,ymax=V3)) + geom_point() + geom_pointrange() + scale_y_continuous(name="Estimate",limits=c(-1,0))
      g

      dev.off()

      