
library(MCMCglmm)
library(matrixcalc)
#number sample
nb_sample=1000
nb_trait=3
nb_pop=7 #ou nb_ageclass ou nb_ de time period
nb_tensor=min((nb_trait*(nb_trait+1))/2,nb_pop-1) #

#########################################################################
############# MAtrice S par itération puis estimation du mode #########
#########################################################################

S_mat=array(,c((nb_trait*(nb_trait+1))/2,(nb_trait*(nb_trait+1))/2,nb_sample))
pG_Ani=array(,c(nb_trait,nb_trait,nb_pop,nb_sample))
mat_val=matrix(c(1:(nb_trait^2)),nrow=nb_trait,byrow=T)
varval= diag(mat_val)
covval=sort(mat_val[upper.tri(mat_val)])

for (s in 1:nb_sample){
  #cration matrice  G par ligne de sample MCMC (n sample = smax)
  for (i in 1:nb_pop){#5 time période
    pG_Ani[,,i,s]=matrix(get(paste("G_Age",i,sep=''))[s,],ncol=nb_trait) #remplacer par nom de sorti adéquate
  }
  #variance de variance
  A=cov(t(apply(pG_Ani[,,,s],3, function(pG) pG[varval])),t(apply(pG_Ani[,,,s],3, function(pG) pG[varval])))
  # covariance de covariance
  D=2*cov(t(apply(pG_Ani[,,,s],3, function(pG) pG[covval])),t(apply(pG_Ani[,,,s],3, function(pG) pG[covval])))
  # covariance entre  covariance et variance
  B=sqrt(2)*cov(t(apply(pG_Ani[,,,s],3, function(pG) pG[varval])),t(apply(pG_Ani[,,,s],3, function(pG) pG[covval])))
  #matrice S
  S_mat[,,s]=rbind(cbind(A,B),cbind(t(B),D))
}

postmod=function(x){posterior.mode(x)}

Smean <- apply(S_mat, 1:2, postmod)
#Smean <- apply(S_mat, 1:2, mean) #utilisation de mean a la place du mode
Smean_eig_values=eigen(Smean)$values[1:nb_tensor]
Smean_eig_vector=eigen(Smean)$vectors[,1:nb_tensor]
Smean_eig_values_mcmc=mcmc(t(apply(S_mat,3,function(x) diag(t(Smean_eig_vector[,1:nb_tensor])%*%x%*%Smean_eig_vector[,1:nb_tensor]))))


#############################################################
############# Construction TENSEUR puis eigenanalysis    ###
#############################################################
Etens=array(,c(nb_trait,nb_trait,nb_tensor))
Etens_valvec=array(,c(nb_trait+1,nb_trait,nb_tensor))

for (i in 1:nb_tensor){  
  diag(Etens[,,i])=Smean_eig_vector[1:nb_trait,i]
  Etens[,,i][lower.tri(Etens[,,i])]=(1/sqrt(2))*Smean_eig_vector[(nb_trait+1):((nb_trait*(nb_trait+1))/2),i]
  Etens[,,i][upper.tri(Etens[,,i])]= t(Etens[,,i])[upper.tri(Etens[,,i])]
  Etens_valvec[,,i]=rbind(c(eigen(Etens[,,i])$values),matrix(c(eigen(Etens[,,i])$vectors),ncol=nb_trait))
  Etens_valvec[,,i]=Etens_valvec[,,i][,order(abs(Etens_valvec[1,,i]),decreasing=T)]
}

# Etens_valvec[,,1]
# Etens_valvec[,,2]
######################################################
############# Projection Coordonnée puis  Variance#########
######################################################

Coord1=array(,c(nb_sample,nb_pop,nb_tensor))
Coord2=array(,c(nb_sample,nb_pop,nb_tensor))
for (i in 1:nb_tensor){  
  Coord1[,,i]=t(apply(pG_Ani, 3:4, function(x,y) {frobenius.prod(y,x)},y=Etens[,,i]))
  Coord2[,,i]=t(apply(pG_Ani, 3:4, function(x,y) {I(frobenius.prod(y,x)^2)/(I(norm(x,type="F")^2))*100},y=Etens[,,i]))
}

##############"Summary des coordonnées sur les deux premiers tenseurs

coord1Tensor1=cbind(posterior.mode(mcmc(Coord1[,,1])),HPDinterval(mcmc(Coord1[,,1])))
coord1Tensor1=round(coord1Tensor1,2)
as.data.frame(paste(coord1Tensor1[,1],"[",coord1Tensor1[,2],":",coord1Tensor1[,3],"]",sep=""))

coord1Tensor2=cbind(posterior.mode(mcmc(Coord1[,,2])),HPDinterval(mcmc(Coord1[,,2])))
coord1Tensor2=matrix(as.character(round(coord1Tensor2,2)),8,3)
as.data.frame(paste(coord1Tensor2[,1],"[",coord1Tensor2[,2],":",coord1Tensor2[,3],"]",sep=""))

##############"Summary des % d'explication sur les deux premiers tenseurs

coord2Tensor1=cbind(posterior.mode(mcmc(Coord2[,,1])),HPDinterval(mcmc(Coord2[,,1])))
coord2Tensor1=round(coord2Tensor1,2)
as.data.frame(paste(coord2Tensor1[,1],"[",coord2Tensor1[,2],":",coord2Tensor1[,3],"]",sep=""))

coord2Tensor2=cbind(posterior.mode(mcmc(Coord2[,,2])),HPDinterval(mcmc(Coord2[,,2])))
coord2Tensor2=round(coord2Tensor2,2)
as.data.frame(paste(coord2Tensor2[,1],"[",coord2Tensor2[,2],":",coord2Tensor2[,3],"]",sep=""))


###########
#Prop_explain_by_tensor
####################################################

###############
results_prop_explain_tensor=Smean_eig_values_mcmc/apply(Smean_eig_values_mcmc,1,sum)
Variance_explain_tensor=cbind(posterior.mode(mcmc(results_prop_explain_tensor)),HPDinterval(mcmc(results_prop_explain_tensor)))

Variance_explain_tensor

###########
#Prop_explain_by_eigen_on_significant_tensor
#############

Prop_explain_by_eigen_on_tensor_1=Etens_valvec[1,,1]^2
Prop_explain_by_eigen_on_tensor_2=Etens_valvec[1,,2]^2
Prop_explain_by_eigen_on_tensor_3=Etens_valvec[1,,3]^2


###########
# eigenvector_on_significant_tensor
#############

e11.vect=Etens_valvec[2:(nb_trait+1),1,1]  #eigenvector of the first eigenvalue of the first eigentensor
e12.vect=Etens_valvec[2:(nb_trait+1),2,1]  #eigenvector of the second eigenvalue of the first eigentensor

e21.vect=Etens_valvec[2:(nb_trait+1),1,2]  #eigenvector of the first eigenvalue of the second eigentensor

e31.proj=Etens_valvec[2:(nb_trait+1),1,3]  #eigenvector of the first eigenvalue of the third eigentensor


#Function to do the projection of Va
proj<- function(G, b) t(b) %*% G %*% (b)

e11.proj <- mcmc(t(apply(pG_Ani, 3:4, proj , b = Etens_valvec[2:4,1,1]))) # Va projected on the first eigenvector of  the first eigentensor
e11.proj <- cbind(posterior.mode(mcmc(e11.proj)),HPDinterval(mcmc(e11.proj)))

e12.proj <- mcmc(t(apply(pG_Ani, 3:4, proj , b = Etens_valvec[2:4,2,1]))) # Va projected on the second eigenvector of  the first eigentensor
e12.proj <-cbind(posterior.mode(mcmc(e12.proj)),HPDinterval(mcmc(e12.proj))) 



e21.proj <- mcmc(t(apply(pG_Ani, 3:4, proj , b = Etens_valvec[2:4,1,2])))
e21.proj <-cbind(posterior.mode(mcmc(e21.proj)),HPDinterval(mcmc(e21.proj)))

e31.proj <- mcmc(t(apply(pG_Ani, 3:4, proj , b = Etens_valvec[2:4,1,3])))
e31.proj <-cbind(posterior.mode(mcmc(e31.proj)),HPDinterval(mcmc(e31.proj)))
