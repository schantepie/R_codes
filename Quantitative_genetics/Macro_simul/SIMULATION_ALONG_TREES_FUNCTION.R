######## FUNCTION AND EXEMPLE OF A NEUTRAL EVOLUTION ALONG A PHYLOGENETIC TREE



VCVevol_brownian<-function(nb_ind,Vm,prob_mut,loci,nb_repeat=1,cpus=1,length.burn,tips,times){
 
  require(ape)
  require(picante)
  require(phytools)
  require(snowfall) # to comment for cluster
  require(cluster) # to comment for cluster
  require(mvMORPH)

  
  ##########################################################################
  ######################" PART 1 : FUNCTIONS TO LOAD #######################
  ##########################################################################"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#--------------------- Function individual based models : drift without epistasis------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  ind_based_pleio<-function(length,nb_ind,Vm,prob_mut,loci,traits,phase,m0_start) {
    
    require(mvtnorm)
    allele=2
    
    ## dissociate burning and poursuit of ind base model
    
    if (phase=="burning") {
      m0=array(0, c(traits,loci,nb_ind, allele))
    }  
    if (phase=="pursuit") {
      m0=m0_start
    }
    
    ## init random mates and random allele (no linkage) for a given length
   
    union_matrix=t(combn(1:nb_ind,2))
    lengh_union_matrix=dim(union_matrix)[1]
    drift=simplify2array(lapply(1:length,function(x,y,z) union_matrix[sample(x=(1:lengh_union_matrix),replace = F, size=z),],y=lengh_union_matrix,z=nb_ind))
    sample_allele_father_mother=array(sample(x=1:2,replace = TRUE,size=nb_ind*loci*2*length),c(nb_ind,loci,2,length))
    
    Gmat<-array(NA, c(length,traits*traits))
    
    ## core drift function where reproduction are occuring and mutation added 
   
    for(i in 1:length){
      
      #reproduction
      m_fath= m0[,,drift[,1,i],]
      m_moth= m0[,,drift[,2,i],]
      
      # for each individuals random samplind of one of the allele
      for (ind in 1:nb_ind){
        m_fath[,sample_allele_father_mother[ind,,1,i]>1,ind,1]=m_fath[,sample_allele_father_mother[ind,,1,i]>1,ind,2]
        m0[,,ind,1]=m_fath[,,ind,1]
        m_moth[,sample_allele_father_mother[ind,,2,i]>1,ind,1]=m_moth[,sample_allele_father_mother[ind,,2,i]>1,ind,2]
        m0[,,ind,2]=m_moth[,,ind,1]
      }
      
      #defining where a mutation occur using a given mutation probability per allele per generation 
      mut_occur=cbind(rbinom(loci*nb_ind*allele, size=1, prob=prob_mut),matrix(0,ncol=(traits-1),nrow=(loci*nb_ind*allele)))
      
      #adding mutation
      if (!all(mut_occur[,1]==0)){
        where_occur=which(mut_occur[,1]==1)
        mut_occur[where_occur,]=rmvnorm(length(where_occur), mean = rep(0, traits), Vm)
        add_occur=array(c(t(mut_occur)), c(traits,loci,nb_ind, allele))
        m0=m0+add_occur
      }else{
        m0=m0
      }
      Gen_val=t(apply(m0,c(1,3),sum)) # calculate genetic value
      Gmat[i,]=c(cov(Gen_val))    # calculate genetic VCV matrix
    }
    return(list(Gmat=Gmat,m0=m0))
  }
  

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#------------------------------ Core function to drift the VCV matrix ----------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


  evol_vcv_tree<-function(rep,phy,nb_ind,Vm,prob_mut,loci,traits){
  
    ##burning phase
    burned=ind_based_pleio(length=phy$lengthburned,nb_ind=nb_ind,Vm=Vm,prob_mut=prob_mut,loci=loci,traits=traits,phase="burning",m0=NULL)
    phy$node.m0[phy$edge[1,1],,,,] <- burned$m0
    
    ## propagation along tree
    edges = list(NULL)
  
    for (i in 1:phy$Nedge) {
      edges[[i]] = matrix(NA,phy$edge.length[i]+1,(1+(traits*traits))) 
      if (phy$edge[i,1]==phy$Nterm+1) {
        anc.age = 0
        anc.m0= phy$node.m0[phy$edge[1,1],,,,]
        anc.vcv =c(cov(t(apply(anc.m0,c(1,3),sum))))
        phy$node.vcv[phy$edge[1,1],] = anc.vcv
      } else {
        anc.edge = match(phy$edge[i,1],phy$edge[,2])
        anc.age = phy$ages[match(phy$edge[i,1],phy$edge[,2])]
        anc.vcv = phy$node.vcv[phy$edge[i,1],]  
        anc.m0= phy$node.m0[phy$edge[i,1],,,,]
      } 
      edges[[i]][,1]<- anc.age:phy$ages[i]
      edges[[i]][1,2:(1+(traits*traits))] <- anc.vcv 
      
      #--call for "ind_based_pleio" function   ---
      model_evo=ind_based_pleio(length=phy$edge.length[i],nb_ind,Vm,prob_mut,loci,traits,phase="pursuit",m0_start=anc.m0)
      edges[[i]][-1,2:((traits*traits)+1)]<- model_evo$Gmat
      phy$node.vcv[phy$edge[i,2],] <- edges[[i]][nrow(edges[[i]]),2:(1+(traits*traits))]
      phy$node.m0[phy$edge[i,2],,,,] <- model_evo$m0
    }
    results=list(vcv_mat=edges,node.m0=phy$node.m0)
    return(results)
  }
  
  
  ##########################################################################
  ######################" PART 2 : parallelise the simulation ##############
  #########################################################################"
  traits=dim(Vm)[1] 

  #------set random seed

  tseed=as.numeric(Sys.time())
  set.seed((tseed - floor(tseed)) * 1e8 )
  
  #------phy construction following Revell et al 2007
  
  b<-exp((log(tips)-log(2))/times)-1
  phy<-pbtree(b=b,n=tips,t=times,type="discrete")
  
  #----- reorder phy 
  
  phy <- reorder(phy, 'cladewise')
  phy$Nterm = length(phy$tip.label)
  phy$Nedge = nrow(phy$edge)
  phy$node.vcv = matrix(NA,nrow=(phy$Nterm+phy$Nnode),ncol=(traits*traits))
  phy$node.m0 =  array(NA, c((phy$Nterm+phy$Nnode),traits,loci,nb_ind, 2)) #2 pour le nombre d'alleles
  phy$ages = node.age(phy)$ages
  phy$lengthburned = length.burn

  
  #----- Parallel evaluation of the evol_vcv_tree of the function
  
    sfInit(parallel=TRUE, cpus=cpus)
    sfExportAll()
    simu<- sfLapply(1:nb_repeat,function(rep,phy,nb_ind,Vm,prob_mut,loci,traits)  evol_vcv_tree(rep,phy,nb_ind,Vm,prob_mut,loci,traits),
                    phy=phy,
                    nb_ind=nb_ind,
                    Vm=Vm,
                    prob_mut=prob_mut,
                    loci=loci,
                    traits=traits)
    sfStop()

  
##### version for cluster without parallelisation
#   nb_repeat=1
#   simu<-lapply(1:nb_repeat,function(rep,phy,nb_ind,Vm,prob_mut,loci,traits)   evol_vcv_tree(rep,phy,nb_ind,Vm,prob_mut,loci,traits),
#                phy=phy,
#                nb_ind=nb_ind,
#                Vm=Vm,
#                prob_mut=prob_mut,
#                loci=loci,
#                traits=traits)
################
  
  Simu_arraylist=simplify2array(simu, higher = TRUE)
  Simu_array=apply(Simu_arraylist,1,function(x) simplify2array(x, higher = TRUE))

  
  # ~~~ estimate a Gmatrix as mean of Gmatrices estimated along tree
  
  sim_vcv_mat=apply(Simu_array$vcv_mat,1,function(y) array(unlist(y), dim = c(dim(y[[1]]), length(y))))
  rearrange<-function(x){
    x=do.call(rbind,x)
    x=x[-duplicated(x),,drop=FALSE] #remove dupication intrduced in the propagation process
    x=apply(x,2,mean)
  }
  Gmean_simul=as.matrix(do.call(rbind, lapply(apply(Simu_array$vcv_mat,2,function(x) simplify2array(x,higher=TRUE)),rearrange))[,-1])
  if(nb_repeat==1) Gmean_simul=t(as.matrix(Gmean_simul))
  
  
  #  ~~~ estimate a theoretic Gmatrix from Lynch et al. 1986
  Gtheo=matrix(rep(c(2*nb_ind*(2*loci*prob_mut*Vm)),nb_repeat),ncol=traits*traits,byrow=T)
  
  #  ~~~ estimate phenotypic values at the tips
  
  terms=which(phy$ages==times) #get the edge number corresponding to the tips ; or terms=which(phy$edge[,2] <= Ntip(phy)) 
  sim_node.m0=Simu_array$node.m0[terms,,,,,,drop=FALSE] # to keep only info from tips
  sim_node.m0=apply(sim_node.m0,c(1,2,4,6),sum)
  Pheno_mean=apply(sim_node.m0,c(1,2,4),mean)
  row.names(Pheno_mean)=phy$tip.label

  ############## Uncomment to have an unbiased estimate of sigma (biais is due to ML estimate)
#    multBM<-function(x,y) mvBM(data=x,tree=y, model="BM1",method="pic",diagnostic = FALSE,echo = FALSE)$sigma
#     correct=length(phy$tip.label)/(length(phy$tip.label)-1)
#     Rate=t(apply(Pheno_mean,3,multBM ,y=phy))
#     Rate_correct=Rate*correct
#     R=t(apply(Pheno_mean,3,multBM,phy)) #to get sigma estimate
#     print(rbind(R,(Gmean_simul/nb_ind),(Gtheo/nb_ind))) # to check difference in sigma rate)
  ###########################"
  
  likmultBM<-function(x,y) mvBM(data=x,tree=phy ,method="pic", scale.height = TRUE, model="BM1",diagnostic = FALSE,echo = FALSE)$LogLik

  multLL<-function(rep,phy,sigma,Pheno_mean,traits){
    sigma_matrix=matrix(sigma[rep,],traits,traits)
    Pheno_mean_byrep=Pheno_mean[,,rep]
    logl=mvLL(data=Pheno_mean_byrep,tree=phy,method="pic",param=list(estim=FALSE, sigma= sigma_matrix))$logl
    return(logl)
  } 
  
  #estimate of different LogLikelihood
  
  LogLik_R=apply(Pheno_mean,3,likmultBM,phy) #estimation of Rate and LogLik (only LogLik is keet)
  LogLik_R0_Gmean_simul=sapply(1:nb_repeat,function(rep,phy,sigma,Pheno_mean,traits) multLL(rep,phy,sigma,Pheno_mean,traits),phy=phy,sigma=(Gmean_simul/nb_ind),Pheno_mean=Pheno_mean,traits=traits) 
  LogLik_R0_Gtheo=sapply(1:nb_repeat,function(rep,phy,sigma,Pheno_mean,traits) multLL(rep,phy,sigma,Pheno_mean,traits),phy=phy,sigma=(Gtheo/nb_ind),Pheno_mean=Pheno_mean,traits=traits) 
  
  ddl=((traits^2-traits)/2)+traits 
  Signi_chi=cbind(LogLik_R,LogLik_R0_Gmean_simul,LogLik_R0_Gtheo,ddl,pchisq(2*(LogLik_R-LogLik_R0_Gmean_simul),ddl ,lower.tail =FALSE),pchisq(2*(LogLik_R-LogLik_R0_Gtheo),ddl ,lower.tail =FALSE))
  colnames(Signi_chi)=c("LogLik_estim_rate","LogLik_Gmean_simul_rate","LogLik_Gtheo_rate","ddl","CHI2_simul","CHI2_theo")

  #simu_results=list(Signi_chi=Signi_chi,pheno.termbranch=sim_node.m0,Gtheo=Gtheo,Gmean_simul=Gmean_simul,VCVmat=sim_vcv_mat,phy=phy)
  simu_results=list(Signi_chi=Signi_chi,phy=phy,Gtheo=Gtheo,Gmean_simul=Gmean_simul)
  return(simu_results)
}


################### Example


# Parameter used in Revell et al 2007
Vm_var=c(0.05,0.1,0.15,0.20)
Vm_cov=c(0.05303301,0.04330127,0,0.09185587,0.07071068,0.12990381)
# Vm_var=c(0.05,0.05)
# Vm_cov=c(0.0)

nb_ind=100
prob_mut=0.0025
loci=20
nb_repeat=1
length.burn=3000
cpus=1
tips=20
times=10000

#construct the VCV mutational matrix
Vm=diag(Vm_var)
Vm[lower.tri(Vm)]=c(Vm_cov)
Vm=t(Vm)
Vm[lower.tri(Vm)]=c(Vm_cov)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#note all replications use the same tree when nb_repeat>1

results_simu=VCVevol_brownian(length.burn=length.burn,nb_ind=nb_ind,Vm=Vm,prob_mut=prob_mut,loci=loci,nb_repeat=5,cpus=5,tips=tips,times=times)#,length.burn=length.burn

results_simu$Signi_chi
