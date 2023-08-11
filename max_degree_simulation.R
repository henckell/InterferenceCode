#script to generate the plot with the slopes of the expected max. degrees
# of the dependency graph
library(igraph)
library(parallel)
library(Matrix)
source("helpers/helpers_paper.R")
source("helpers/helpers_graph.R")

nvals <- c(300,600,1200,2400,4800)
prob <- 0.1
n.cores <- 1
K <- 100


results <- lapply(1:length(nvals), FUN = function(i){
  
  nval <- nvals[i]
  
  res <- mclapply(1:K, function(k){
    
    A <- A_random_graph(nval,0.2)
    A.dep <- Matrix(apply((A+A%*%A),1, function(r) ifelse(r>=1, 1,0)), sparse=T)
    diag(A.dep) <- 0
    max.deg.constp <- max(apply(A.dep, 1, sum))
    rm(A)
    rm(A.dep)
    
    

    A <- A_random_graph(nval, nval^(-2/3))
    A.dep <- Matrix(apply((A+A%*%A),1, function(r) ifelse(r>=1, 1,0)), sparse=T)
    diag(A.dep) <- 0
    max.deg.12 <- max(apply(A.dep, 1, sum))
    rm(A)
    rm(A.dep)

    
    A <- A_random_graph(nval, 10/nval)
    A.dep <- Matrix(apply((A+A%*%A),1, function(r) ifelse(r>=1, 1,0)), sparse=T)
    diag(A.dep) <- 0
    max.deg.non <- max(apply(A.dep, 1, sum))
    rm(A)
    rm(A.dep)
    

    print(k)
    
    return(list(
      max.deg.constp=max.deg.constp,
                max.deg.12=max.deg.12, 
                max.deg.non=max.deg.non 
                ))
  
    }, mc.cores=n.cores)
  
  print(paste0("nval: ", nval))
  return(res)
  
})



filename <- paste0("depgraphs")

filename <- gsub("[.]","",filename)
assign(filename, results)
save(results, file=paste(filename, ".Rda",sep=""))

################################# Auswertung ###############################

nvals <- c(300,600,1200,2400,4800)

load("depgraphs.Rda")

res <- lapply(1:length(results), function(k){

  as.data.frame(do.call(rbind, results[[k]]))

})


max.degs.12 <- sapply(1:length(res), function(k){
  mean(unlist(res[[k]]$max.deg.12))
})

max.degs.non <- sapply(1:length(res), function(k){
  mean(unlist(res[[k]]$max.deg.non))
})


rm(res)

load("prob02depgraphs.Rda")



res <- lapply(1:length(results), function(k){

  as.data.frame(do.call(rbind, results[[k]]))

})


max.degs.constp <- sapply(1:length(res), function(k){
  mean(unlist(res[[k]]$max.deg.constp))
})


fit.constp <- unname(lm(log(max.degs.constp)~log(nvals))$coef[2])
fit.12 <- unname(lm(log(max.degs.12)~log(nvals))$coef[2])
fit.non <- unname(lm(log(max.degs.non)~log(nvals))$coef[2])



par(mfrow=c(1,1))
filename.save <- "max_degree"

pdf(file = paste0(filename.save , ".pdf"), width = 11, height = 4.5)

df <- data.frame(y=c(max.degs.non,max.degs.12,max.degs.constp),x=rep(nvals,3),group=as.factor(c(rep("1",5),rep(2,5),rep(3,5))))

legend_title <- "Graph type"

g <- ggplot(data=df,aes(x=log(x),y=log(y)))+ geom_point(aes(shape = group)) + geom_line(aes(linetype=group)) + ylab("") + xlab("log(N)") +
    scale_shape_manual(legend_title,values=c(0, 1, 2),labels = c(TeX("$I(N,10/N)$: slope=1"), TeX("$I(N,N^{-2/3})$: slope=0.63"), TeX("$I(N,0.2)$: slope=0.17")))+
    scale_linetype_manual(legend_title,values=c(1,2,3),labels = c(TeX("$I(N,10/N)$: slope=1"), TeX("$I(N,N^{-2/3})$: slope=0.63"), TeX("$I(N,0.2)$: slope=0.17")))+
    theme(legend.position="right",legend.key.size = unit(1.5, 'cm'),legend.title = element_text(size=16),legend.text = element_text(size=12))

g

dev.off()
