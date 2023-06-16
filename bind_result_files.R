
Res <- list()


for(i in 1:2){
  load(paste0(i,"estimation_rand_npfix_pi07_eta02_model3_nval4800.Rda"))
  Res <- c(Res,Results)
  
  rm(Results)
}

Results <- Res

filename <- "estimation_rand_npfix_pi07_eta02_model3_nval4800"
assign(filename, Results)
save(Results, file=paste(filename, ".Rda",sep=""))

for(i in 1:200){
  dat <- data_model2_multiple_C(A,beta_0, beta_1, delta_C1,delta_C2, delta_C1C2,
                         sigma_Y,pai, eta,
                         typeofgraph, use.only.feat1, 
                         nval.eff, 
                         error.type,
                         do.intervene=NULL)
  print(dat$tau)
  print(dat$omega_0)
  print(dat$omega_1)
  
}

