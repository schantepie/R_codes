library(snow)
library(MCMCglmm)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~ CONDITIONAL EVOLVABILITY ~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



#### internal function Cond_evol

Cond_evol<-function(x) {
  G=matrix(x,numtr,numtr)
  Ginv=solve(G)
  mean(apply(beta_positive,1,function(x) {1/(t(x)%*%Ginv%*%x)}))
}

##### Function FOR A MCMC OUTPUT

Cond_evol_mcmc<-function(x,cpus=1) {
  require(MCMCglmm)
  numtr=sqrt(dim(x)[2])
  trial=10000
  beta_positive=matrix(runif(numtr*trial,-1,1),ncol=numtr)
  clus <- makeCluster(cpus)
  clusterExport(clus, list("beta_positive","numtr"), envir=environment())
  CE=parApply(clus,x,1,Cond_evol)
  stopCluster(clus)
  CE=mcmc(CE)
  sum_CE=cbind(posterior.mode(CE),HPDinterval(CE))
  return(sum_CE)
}


###exemple
### CE1_pos=Cond_evol_mcmc(G_Age1,6) with G_Age1 a mcmc output and 6 the number of cpus

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~ RESPONDABILITY ~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#### internal function respondability

respond<-function(x) {
  G=matrix(x,numtr,numtr)
  G2=G^2
  mean(apply(beta_positive,1,function(x) {sqrt((t(x)%*%G2%*%x))}))
}

##### Function FOR A MCMC OUTPUT

respond_mcmc<-function(x,) {
  require(MCMCglmm)
  numtr=sqrt(dim(x)[2])
  trial=10000
  beta_positive=matrix(runif(numtr*trial,-1,1),ncol=numtr)
  clus <- makeCluster(6)
  clusterExport(clus, list("beta_positive","numtr"), envir=environment())
  CE=parApply(clus,x,1,respond)
  stopCluster(clus)
  CE=mcmc(CE)
  sum_CE=cbind(posterior.mode(CE),HPDinterval(CE))
  return(sum_CE)
}
