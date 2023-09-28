# Main: Simulation Study for the Paper
setwd("~/GitHub/InvarianceCode")

library("igraph")
library("Matrix")
library("parallel")
library("expm")

source("helpers/helpers_general.R")
source("helpers/helpers_graph.R")
source("helpers/helpers_features.R")
source("helpers/helpers_main.R")
source("helpers/helpers_datageneration.R")


nval_list <- c(300,600,1200,2400,4800)
graph.type_list <-  c("rand_pfix","rand_npfix","rand_npfix_growing","family","2dlatt")


pai_list <- c(0.7,0.7,0.7,1,0.5)
eta_list<- c(0.2,0.2,0.2,0,0.1)
features_list <- list(c(1),c(1),c(1),c(1),c(1,3))

for(i in 1:length(graph.type_list)){
  
  typeofgraph <- graph.type_list[i]
  pai <- pai_list[i]
  eta <- eta_list[i]
  features <- features_list[[i]]

  # nval <- 600
  n.cores.graph = 5 #12
  n.cores.data = 5
  useseed <- 1
  
  nrep_graph <- 50
  nrep_data <- 100
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
  # beta_0 = c(2, 1, 0.5)
  # beta_1 = c(2.4, 2.1, 1.0)
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
  
  for(j in 1:length(nval_list)){
    nval <- nval_list[j]
  
  foldername <- paste0("results_data/model_", typeofgraph) 

if (file.exists(foldername)) {
  setwd(foldername)
  
} else {
  dir.create(foldername)
  setwd(foldername)
  
}

# save setting
parameter <- list(
  beta = cbind(beta_0,beta_1), 
  B_C=B_C,
  eta_C = eta_C,
  gamma_C = gamma_C,
  delta_C = delta_C,
  sigma_C = sigma_C,
  delta_Y = delta_Y,
  sigma_Y = sigma_Y,
  nrep_graph = nrep_graph,
  nrep_data = nrep_data,
  nval=nval,
  typeofgraph =typeofgraph,
  error.types = c(error.type_C,error.type_Y),
  pai=pai,
  eta=eta,
  prob =prob, 
  const=const,
  prob.rewiring=prob.rewiring,
  growth.rate=growth.rate,
  growth.rate.WS=growth.rate.WS,
  features=features,
  useseed=useseed,
  n.cores.data=n.cores.data,#11
  n.cores.graph=n.cores.graph)

# capture.output(parameter, file=paste0("est_params_", foldername, ".txt"))


simulation(   beta = cbind(beta_0,beta_1), 
              B_C=B_C,
              eta_C = eta_C,
              gamma_C = gamma_C,
              delta_C = delta_C,
              sigma_C = sigma_C,
              delta_Y = delta_Y,
              sigma_Y = sigma_Y,
              nrep_graph = nrep_graph,
              nrep_data = nrep_data,
              nval=nval,
              typeofgraph =typeofgraph,
              error.types = c(error.type_C,error.type_Y),
              pai=pai,
              eta=eta,
              prob =prob, ## double check what this and const do
              const=const,
              prob.rewiring=prob.rewiring,
              growth.rate=growth.rate,
              growth.rate.WS=growth.rate.WS,
              features=features,
              useseed=useseed,
              n.cores.data=n.cores.data,#11
              n.cores.graph=n.cores.graph)

setwd("..")
}}

