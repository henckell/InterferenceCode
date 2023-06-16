# Main: Simulation Study for the Paper
# 27.06.2022
setwd("~/Leonard/Papers/Projects/Interference/Rcode/ML-code")

library("igraph")
library("Matrix")
library("parallel")
library("expm")


#setwd("~/GitHub/Confounded-GATE/SimulationStudy_Paper/Dec22")

source("helpers_paper.R")


nval_list <- c(300,600,1200,2400,4800)
graph.type_list <- c("rand_pfix","rand_npfix","rand_npfix_growing","family","2dlatt")


pai_list <- c(0.7,0.7,0.7,1,0.5)
eta_list<- c(0.2,0.2,0.2,0,0.1)
use.only.feat1_list <- c(1,1,1,1,0)
growth.rate_list <- c(0,2/3,)

# nval <- 600
# model <- 2
# n.cores.graph = 1 #12
# n.cores.data = 1
# useseed <- 1
# 
# eff <- "global" #"EATE"
# 
# # 300  600 1200 2400 4800 9600
# #nvals <- 2^(0:5)*300 #2^(2:4)*100 #2^(2:8)*100
# 
# nrep_graph <- 50
# nrep_data <- 100
# # typeofgraph <- "rand_npfix" #"rand_npfix_growing" #"WS_growing", "rand_npfix_growing", "family"
# error.type_C <- "rnorm" #, "rnorm", "rt", "runif", "chisq"
# error.type_Y <- "runif"
# #B <- 1000 # nr of bootstrap repetitions
# # pai <- 0.7
# # eta <- 0.2
# # growth.rate <- 2/3
# # use.only.feat1 <- TRUE
# 
# 
# prob <- 0.2 #runif(1, 0.001, 0.9) # probability of an edge in random graphs
# prob.rewiring <- 0.05#runif(1, 0.0001, 0.5) # rewiring probability of an edge in WS
# sigma_Y <- 1#runif(1,0.1, 4) #sigma_Y
# const <- 10#2
# eta_C1 <- eta_C2 <- 1#runif(1, 0.5, 5) # effects of C's on treatment
# sigma_C <- 1
# sigma_Z <- 1.5
# sigma_P <- 1.5 
# eta_Z <-  2.3
# delta_P <- 1
# # growth.rate.WS <- 1/4 # need multiple options?
# beta_0 = c(2, 1, 0.5)
# beta_1 = c(2.4, 2.1, 1.0)
# delta_C1 = 2
# delta_C2 = 1.5
# delta_C1C2 = 2






# interference, little C
# setting1 <- list(
#   beta_0 = c(0, 0.3, 0.2),
#   beta_1 = c(1.5, 1, 0.5),
#   delta_C1 = 0.1,
#   delta_C2 = 0.1
# )

# interference, heavier C
# setting2 <- list(
#   beta_0 = c(2, 1, 0.5),
#   beta_1 = c(2.4, 2.1, 1.0),
#   delta_C1 = 2,
#   delta_C2 = 1.5,
#   delta_C1C2 = 2
# )

# approx SUTVA, little C
# setting3 <- list(
#   beta_0 = c(0, 0.1, 0.05),
#   beta_1 = c(1, 0.05, 0.05),
#   delta_C1 = 0.1,
#   delta_C2 = 0.1
# )

# approx SUTVA, heavier C
# setting4 <- list(
#   beta_0 = c(0, 0.1, 0.05),
#   beta_1 = c(1, 0.05, 0.05),
#   delta_C1 = 2,
#   delta_C2 = 2
# )

for(i in 1:length(graph.type_list)){
  
  graph.type <- graph.type_list[i]
  pai <- pai_list[i]
  eta <- eta_list[i]
  use.only.feat1 <- use.only.feat1_list[i]
  
  
  # nval <- 600
  model <- 2
  n.cores.graph = 1 #12
  n.cores.data = 1
  useseed <- 1
  
  eff <- "global" #"EATE"
  
  # 300  600 1200 2400 4800 9600
  #nvals <- 2^(0:5)*300 #2^(2:4)*100 #2^(2:8)*100
  
  nrep_graph <- 50
  nrep_data <- 100
  # typeofgraph <- "rand_npfix" #"rand_npfix_growing" #"WS_growing", "rand_npfix_growing", "family"
  error.type <- "runif"
  # error.type_C <- "rnorm" #, "rnorm", "rt", "runif", "chisq"
  # error.type_Y <- "runif"
  #B <- 1000 # nr of bootstrap repetitions
  # pai <- 0.7
  # eta <- 0.2
  # growth.rate <- 2/3
  # use.only.feat1 <- TRUE
  
  
  prob <- 0.2 #runif(1, 0.001, 0.9) # probability of an edge in random graphs
  prob.rewiring <- 0.05#runif(1, 0.0001, 0.5) # rewiring probability of an edge in WS
  sigma_Y <- 1#runif(1,0.1, 4) #sigma_Y
  const <- 10#2
  eta_C1 <- eta_C2 <- 1#runif(1, 0.5, 5) # effects of C's on treatment
  sigma_C <- 1
  sigma_Z <- 1.5
  sigma_P <- 1.5 
  eta_Z <-  2.3
  delta_P <- 1
  growth.rate.WS <- 1/4
  beta_0 = c(2, 1, 0.5)
  beta_1 = c(2.4, 2.1, 1.0)
  delta_C1 = 2
  delta_C2 = 1.5
  delta_C1C2 = 2
  
  for(j in 1:length(nval_list)){
    nval <- nval_list[i]

if(typeofgraph == "rand_npfix_growing"){
  
  foldername <- paste0(Sys.Date(),eff, "_model", model,typeofgraph,
                       "_growth.rate",round(growth.rate,3),"use.only.feat1",
                       use.only.feat1, typeofgraph,"_pi",pai,"_eta",eta
  ) 
}else{
  
  foldername <- paste0(Sys.Date(),eff,"_model", model,typeofgraph,"_growth.rate.WS",round(growth.rate.WS,3),"use.only.feat1",
                       use.only.feat1, typeofgraph,"_pi",pai,"_eta",eta
  )
}


if (file.exists(foldername)) {
  setwd(foldername)
  
} else {
  dir.create(foldername)
  setwd(foldername)
  
}

# save setting
parameter <- list(
  beta_0 = beta_0, 
  beta_1 = beta_1, 
  delta_C1 = delta_C1, 
  delta_C2= delta_C2, 
  delta_C1C2=delta_C1C2,
  sigma_Y = sigma_Y,
  eta_Z = eta_Z,
  nrep_graph = nrep_graph,
  nrep_data = nrep_data,
  typeofgraph =typeofgraph,
  error.type = error.type,
  pai=pai,
  eta=eta,
  prob =prob, 
  const=const,
  eta_C1 = eta_C1,
  eta_C2 = eta_C2,
  eta_Z  = eta_Z,
  sigma_C =sigma_C,
  sigma_Z = sigma_Z,
  sigma_P = sigma_P,
  prob.rewiring=prob.rewiring,
  #B=B,
  growth.rate=growth.rate,
  growth.rate.WS=growth.rate.WS,
  use.only.feat1 =use.only.feat1,
  useseed=useseed)

# capture.output(parameter, file=paste0("est_params_", foldername, ".txt"))


simulation(  beta_0 = beta_0, 
             beta_1 = beta_1, 
             delta_C1 = delta_C1, 
             delta_C2= delta_C2, 
             delta_C1C2=delta_C1C2,
             nrep_graph=nrep_graph, 
             nrep_data=nrep_data,
             n.cores.data,
             n.cores.graph,
             typeofgraph=typeofgraph,
             error.type=error.type,
             pai=pai,
             eta=eta,
             #B=B,
             prob=prob,
             prob.rewiring=prob.rewiring,
             sigma_Y=sigma_Y, 
             const=const, 
             eta_C1=eta_C1, 
             eta_C2=eta_C2, 
             sigma_C=sigma_C,
             sigma_Z= sigma_Z,
             eta_Z= eta_Z, 
             nval=nval,
             growth.rate=growth.rate,
             growth.rate.WS=growth.rate.WS,
             use.only.feat1=use.only.feat1,
             model=model,
             eff=eff,
             useseed=useseed)

setwd("..")
}}

