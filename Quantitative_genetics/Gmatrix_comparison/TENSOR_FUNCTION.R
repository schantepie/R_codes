

tensor_analyse<-function(names,nb_sample,names_traits,names_pop,allresults=FALSE){
require(MCMCglmm)
require(matrixcalc)
nb_trait=length(names_traits)
nb_pop=length(names_pop)

nb_tensor=min((nb_trait*(nb_trait+1))/2,nb_pop-1) 


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
  if(allresults==TRUE){
    model=get(paste(names,i,sep=''))
    model=model$VCV
    model=model[,grep (".animal",colnames(model))]
    pG_Ani[,,i,s]=matrix(model[s,],ncol=nb_trait)
    }else{
    pG_Ani[,,i,s]=matrix(get(paste(names,i,sep=''))[s,],ncol=nb_trait) #remplacer par nom de sorti adéquate
    }
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

postmod=function(x){posterior.mode(mcmc(x))}

Smean <- apply(S_mat, 1:2, postmod)
#Smean <- apply(S_mat, 1:2, mean) #utilisation de mean a la place du mode
Smean_eig_values=eigen(Smean)$values[1:nb_tensor]
Smean_eig_vector=eigen(Smean)$vectors[,1:nb_tensor,drop=FALSE]
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

#############################################################
############# Projection sur Coordonnée puis  Variance  #########
#############################################################

Coordinate=array(,c(nb_sample,nb_pop,nb_tensor))
Coordinate_percent=array(,c(nb_sample,nb_pop,nb_tensor))
for (i in 1:nb_tensor){  
  Coordinate[,,i]=t(apply(pG_Ani, 3:4, function(x,y) {frobenius.prod(y,x)},y=Etens[,,i]))
  Coordinate_percent[,,i]=t(apply(pG_Ani, 3:4, function(x,y) {I(frobenius.prod(y,x)^2)/(I(norm(x,type="F")^2))*100},y=Etens[,,i]))
}

##############"Summary des coordonnées sur les deux premiers tenseurs

Summa_coordinate=simplify2array(lapply(1:nb_tensor,function(x) cbind(posterior.mode(mcmc(Coordinate[,,x])),HPDinterval(mcmc(Coordinate[,,x])))))
dimnames(Summa_coordinate)=list(names_pop,c("Coord","lower","upper"),paste("tensor",1:nb_tensor,sep=""))

##############"Summary des % d'explication sur les deux premiers tenseurs

Summa_coordinate_percent=simplify2array(lapply(1:nb_tensor,function(x) cbind(posterior.mode(mcmc(Coordinate_percent[,,x])),HPDinterval(mcmc(Coordinate_percent[,,x])))))
dimnames(Summa_coordinate_percent)=list(names_pop,c("Coord","lower","upper"),paste("tensor",1:nb_tensor,sep=""))

####################################################
#Prop_explain_by_tensor
####################################################

results_prop_explain_tensor=Smean_eig_values_mcmc/apply(Smean_eig_values_mcmc,1,sum)
Variance_explain_tensor=cbind(posterior.mode(mcmc(results_prop_explain_tensor)),HPDinterval(mcmc(results_prop_explain_tensor)))
dimnames(Variance_explain_tensor)=list(paste("tensor",1:nb_tensor,sep=""),c("%explained","lower","upper"))

###########
#Prop_explain_by_eigen_on_significant_tensor
#############

Prop_explain_by_eigen_on_tensor=lapply(1:nb_tensor,function(x) Etens_valvec[1,,x]^2)
Prop_explain_by_eigen_on_tensor=t(simplify2array(Prop_explain_by_eigen_on_tensor))
dimnames(Prop_explain_by_eigen_on_tensor)=list(paste("tensor",1:nb_tensor,sep=""),paste("%_eigen_value_", 1:nb_trait,sep=""))

###########
# eigenvector_on_significant_tensor
#############

eigenvect.vect=simplify2array(lapply(1:nb_tensor, function(x) Etens_valvec[2:(nb_trait+1),,x]))
dimnames(eigenvect.vect)=list(names_traits, paste("eigen_vector_", 1:nb_trait,sep=""),paste("tensor",1:nb_tensor,sep=""))

###########
# projection of Va
#############

#Function to do the projection of Va
proj<- function(G, b) t(b) %*% G %*% (b)
eigenvect.proj=array(NA,c(nb_pop,3,nb_trait,nb_tensor))

for (i in 1:nb_tensor){
  for (j in 1:nb_trait){
  project=mcmc(t(apply(pG_Ani, 3:4, proj , b = Etens_valvec[2:4,j,i])))
  project=cbind(posterior.mode(mcmc(project)),HPDinterval(mcmc(project)))
  eigenvect.proj[,,j,i]=project
  }
}

dimnames(eigenvect.proj)=list(names_pop,c("Va","lower","upper"),paste("eigen_vector_", 1:nb_trait,sep=""),paste("tensor",1:nb_tensor,sep=""))

return(list(Variance_explain_tensor=Variance_explain_tensor,
            Summa_coordinate=Summa_coordinate,
            Summa_coordinate_percent=Summa_coordinate_percent,
            Prop_explain_by_eigen_on_tensor=Prop_explain_by_eigen_on_tensor,
            eigenvect.vect=eigenvect.vect,
            eigenvect.proj=eigenvect.proj))
}

# ####### Example
# 
# ### load the data with "names" followed by increasing suite of number
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
# ### set value for function
# 
# nb_sample=1000
# names="Gmat_"
# names_pop=c("bosh","calix","peerd","korsten","rouviere","muroCB","muroCV","pirio")
# names_traits=c("tars","wing","mass")
# 
# ### run function
# bluemorph_tensor=tensor_analyse(names,nb_sample,name_traits,names_pop)
# 
# 
# summary(bluemorph_tensor)
# bluemorph_tensor$Variance_explain_tensor
# bluemorph_tensor$Prop_explain_by_eigen_on_tensor
# bluemorph_tensor$Summa_coordinate
# bluemorph_tensor$eigenvect.vect 
