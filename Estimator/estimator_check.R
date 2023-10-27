## file to verify property of estimator

setwd("~/GitHub/InterferenceCode")

library("igraph")
library("Matrix")
library("parallel")
library("expm")

source("helpers/helpers_general.R")
source("helpers/helpers_graph.R")
source("helpers/helpers_features.R")
source("helpers/helpers_main.R")
source("helpers/helpers_datageneration.R")
source("helpers/helpers_estimator.R")


typeofgraph <- "rand_npfix"
pai <- 0.7
eta <- 0.2
features <- c(1)

error.type_C <- "rnorm"
error.type_Y <- "runif"  #, "rnorm", "rt", "runif", "chisq"
growth.rate <- 2/3
do.intervene <- c()

prob <- 0.2 
prob.rewiring <- 0.05
delta_Y <- 0
sigma_Y <- 1
const <- 10
growth.rate.WS <- 1/4
B_C <- t(matrix(c(0,0,0,
                  2,0,0,
                  0,0,0),nrow=3))
delta_C <- c(-2,0,0.5)
gamma_C <- c(1.5,0,0)
sigma_C <- c(1,1,1)
eta_C <-  c(0,1,5)

beta_0 = c(2, 1, 0.5)
beta_1 = c(2.4, 2.1, 1.0)

adj <- c(2)

beta <- cbind(beta_0,beta_1)

nval <- nval.eff <- 1200

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




check_hat_tau <- c()
check_var_hat_tau <- c()

for(i in 1:50){

dat <- data_gen(A,beta, 
                B_C,eta_C,gamma_C,delta_C,sigma_C,
                delta_Y,sigma_Y,
                typeofgraph,
                features, 
                nval, 
                error.type, 
                do.intervene)

res <- estimator(pai, eta, W=dat$treat,Y=dat$Y,C=dat$C,adj=c(2),feat_functions=list(feat_X1),A,B=50)

check_hat_tau[i] <- res$hat_tau
check_var_hat_tau[i] <- res$hat_var_tau
}

sqrt(mean((check_hat_tau-dat$tau)^2))

sqrt(mean((var(check_hat_tau)*nval - check_var_hat_tau)^2))
