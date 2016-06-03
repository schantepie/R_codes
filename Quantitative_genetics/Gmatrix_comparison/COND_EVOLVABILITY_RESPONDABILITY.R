###############################################################
##############" FUNCTION CONDITIONAL EVOLVBILITY #############
################   AND RESPONDABILITY #########################
###############################################################

Cond_Evolvability_Respondability <- function(names,nb_pop,names_pop,cpus=1){
  
  require(MCMCglmm)
  require(snow)
  
  Cond_Evolvability=array(NA,c(nb_pop,3))
  Respondability=array(NA,c(nb_pop,3))
  
for(i in 1:nb_pop){
  
gmat=get(paste("Gmat_",i,sep=""))
    
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

Cond_evol_mcmc<-function(gmat,cpus) {
  require(MCMCglmm)
  numtr=sqrt(dim(gmat)[2])
  trial=10000
  beta_positive=matrix(runif(numtr*trial,-1,1),ncol=numtr)
  clus <- makeCluster(cpus)
  clusterExport(clus, list("beta_positive","numtr"), envir=environment())
  CE=parApply(clus,gmat,1,Cond_evol)
  stopCluster(clus)
  CE=mcmc(CE)
  sum_CE=cbind(posterior.mode(CE),HPDinterval(CE))
  return(sum_CE)
}

Cond_Evolvability[i,]=Cond_evol_mcmc(gmat,cpus)
dimnames(Cond_Evolvability)=list(names_pop,c("mode","lower","upper"))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~ RESPONDABILITY ~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#### internal function respondability

respond<-function(x) {
  G=matrix(x,numtr,numtr)
  G2=G%*%G
  mean(apply(beta_positive,1,function(x) {sqrt((t(x)%*%G2%*%x))}))
}

##### Function FOR A MCMC OUTPUT

respond_mcmc<-function(gmat,cpus) {
  require(MCMCglmm)
  numtr=sqrt(dim(gmat)[2])
  trial=10000
  beta_positive=matrix(runif(numtr*trial,-1,1),ncol=numtr)
  clus <- makeCluster(6)
  clusterExport(clus, list("beta_positive","numtr"), envir=environment())
  CE=parApply(clus,gmat,1,respond)
  stopCluster(clus)
  CE=mcmc(CE)
  sum_CE=cbind(posterior.mode(CE),HPDinterval(CE))
  return(sum_CE)
}

Respondability[i,]=respond_mcmc(gmat,cpus)
dimnames(Respondability)=list(names_pop,c("mode","lower","upper"))
  }
  
  return (list(
    Cond_Evolvability=Cond_Evolvability,
    Respondability=Respondability ))
}



###################  EXEMPLE for 8 population
# 
# 
# 
# load("/media/mnhn/Leca/Gmatrix_project/5_Intra_species_variations/blue_tits/morph/blue_bosh_morph_SI_all_meanstand_OK.Rdata")
# assign(paste("Gmat_",1,sep=""),res)
# load("/media/mnhn/Leca/Gmatrix_project/5_Intra_species_variations/blue_tits/morph/blue_calix_morph_SI_all_meanstand_OK.Rdata")
# assign(paste("Gmat_",2,sep=""),res)
# load("/media/mnhn/Leca/Gmatrix_project/5_Intra_species_variations/blue_tits/morph/blue_peerd_morph_SI_all_meanstand_OK.Rdata")
# assign(paste("Gmat_",3,sep=""),res)
# load("/media/mnhn/Leca/Gmatrix_project/5_Intra_species_variations/blue_tits/morph/Blue_korst_morph_SI_all_meanstand_OK.Rdata")
# assign(paste("Gmat_",4,sep=""),res)
# load("/media/mnhn/Leca/Gmatrix_project/5_Intra_species_variations/blue_tits/morph/blue_montrouv_morph_SI_all_meanstand_OK.Rdata")
# assign(paste("Gmat_",5,sep=""),res)
# load("/media/mnhn/Leca/Gmatrix_project/5_Intra_species_variations/blue_tits/morph/blue_corsmuro_CB_morph_SI_all_meanstand_OK.Rdata")
# assign(paste("Gmat_",6,sep=""),res)
# load("/media/mnhn/Leca/Gmatrix_project/5_Intra_species_variations/blue_tits/morph/blue_corsmuro_CV_morph_SI_all_meanstand_OK.Rdata")
# assign(paste("Gmat_",7,sep=""),res)
# load("/media/mnhn/Leca/Gmatrix_project/5_Intra_species_variations/blue_tits/morph/blue_pirio_morph_SI_all_meanstand_OK.Rdata")
# assign(paste("Gmat_",8,sep=""),res)
# 
# 
# 
# names_pop=c("bosh","calix","peerd","korsten","rouviere","muroCB","muroCV","pirio")
# names="Gmat_"
# nb_pop=8
# 
# CER=Cond_Evolvability_Respondability(names,nb_pop,names_pop,cpus=6)
# 
# 
# 

