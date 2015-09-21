
library(ape)
library(picante)
library(truncnorm)
library(phytools)

set.seed(99)
tree<-pbtree(n=10,scale=100)
tr=rtree(7)
plot(tree)


G0=matrix(c(2,0.1,1,
            0.1,2,1.8,
            1,1.8,2.2),3)

a=VCVevol_brownian(G0,tree, 0.01,0.01)

tree=phy_pass

VCVevol_brownian<-function(VCVmatrix,tree,sigma_var,sigma_cor,nb_repeat){
phy=tree
x.root=0 #root value
sigma = 0.3 #brownian motion st. dev.
trend = 0 #brownian motion trend
trait.mode = 'Brownian'
set.minbl = 1 #length of shortest branch


CorG0=cov2cor(G0)
cor= CorG0[lower.tri(CorG0)]
V=diag(G0)

x.root=c(V,cor) 

nb_correl=(dim(CorG0)[1]^2-dim(CorG0)[1])/2
nb_traits=dim(CorG0)[2]

######################" FUNCTION TO LOAD


drift_correl<-function(cormat,len,sd){
  options(digits=8)
  dime=dim(cormat)[1]
  x=0
  store=array(NA,c(dime,dime,len))
  correl=array(NA,c(len,sqrt(length(cormat))))
  while(x<len){
    x=x+1
    if(x==1) {
      cho2=t(chol(cormat))
      cho4=matrix(0,dime,dime)
      colnames(cho4)=1:dime
      row.names(cho4)=1:dime
      row.names(cormat)=1:3
      colnames(cormat)=1:3
    }
    
    ang2=matrix(0,dime,dime)
    cho4[1,1]=1
    cho4[2:dime,1]=apply(matrix(cho2[2:dime,1]),2,FUN=function(x) {rtruncnorm(1, a=-0.999999, b=0.999999, mean = x, sd = sd)} )
    ang2[2:dime,1]=acos(cho4[2:dime,1])
    cho4[2,2]=sin(ang2[2,1])
    
    cho4[3,2]=apply(matrix(cho2[3,2]/sin(ang2[3,1])), 2,FUN=function(x) {rtruncnorm(1, a=-0.999999, b=0.999999, mean = x, sd = sd)})*sin(ang2[3,1])
    ang2[3,2]=acos(cho4[3,2]/sin(ang2[3,1]))
    
    cho4[3,3]=sin(ang2[3,2])*sin(ang2[3,1])
    reconst=cho4%*%t(cho4)
    
    # store[,,x]<-reconst
    store[,,x]<-reconst[order(row.names(reconst)),order(colnames(reconst))]
    cor=store[,,x]
    correl[x,]= cor[lower.tri(cor)]
    rand=sample(1:3)
    # cormat=cormat[rand,rand]
    reconst=reconst[rand,rand]
    cho2=t(chol(reconst)) 
    # cho2=t(chol(cormat))
    # eigen(reconst)
    cho4=matrix(0,3,3)
    colnames(cho4)=colnames(reconst)
    row.names(cho4)=row.names(reconst)
  }
  return(correl)
}

###
cor2cov<-function(cor=NULL,va=NULL){
  d <- sqrt(va)
  COR=diag(length(cor))
  COR[lower.tri(COR)]=cor
  COR[upper.tri(COR)]= cor
  V <- outer(d, d) * COR
  return(V) 
}

### line version
cor2cov2<-function(corva=NULL, nbtraits=3){
  d <- sqrt(corva[1:nbtraits])
  cor=corva[(nbtraits+1):length(corva)]
  COR=diag(length(cor))
  COR[lower.tri(COR)]=cor
  COR[upper.tri(COR)]= cor
  V <- outer(d, d) * COR
  
  
  return(V) 
}

#########################################

    # I think this code requires the tree to be in cladewise order -- pdc
    phy <- reorder(phy, 'cladewise')
    #print(mu)
    phy$edge.length.original = phy$edge.length
    # minbl = min(phy$edge.length)
    phy$edge.length = round(set.minbl*phy$edge.length/minbl)
    phy$Nterm = length(phy$tip.label)
    phy$Nedge = nrow(phy$edge)
    phy$node.traits = matrix(NA,nrow=(phy$Nterm+phy$Nnode),ncol=nb_traits+nb_correl)
    phy$node.traits[phy$edge[1,2],] = rep(0,nb_traits+nb_correl)#x.root
    phy$ages = node.age(phy)$ages
    min.trait = x.root
    max.trait = x.root
    edges = list(NULL)
    simu=list()
    
    nbsim=1
    i=2
    for(nbsim in 1:nb_repeat){
    for (i in 1:phy$Nedge) {
      
    
      edges[[i]] = matrix(NA,phy$edge.length[i]+1,(1+nb_traits+nb_correl))
      if (phy$edge[i,1]==phy$Nterm+1) {
        anc.age = 0
        anc.trait =x.root
      } else {
        anc.edge = match(phy$edge[i,1],phy$edge[,2])
        anc.age = phy$ages[match(phy$edge[i,1],phy$edge[,2])]
        anc.trait = phy$node.traits[phy$edge[i,1],]  
      }
      edges[[i]][,1]=anc.age:phy$ages[i]
      edges[[i]][1,2:(1+nb_traits+nb_correl)] = anc.trait
      
      if (length(grep(trait.mode,'Brownian'))==1) {
        
        #--------------variance bounded ZERO ----------------
        for (v in 2:(1+nb_traits)){
          for (j in 2:nrow(edges[[i]])) {
            repeat {
              edges[[i]][j,v] = edges[[i]][(j-1),v] + rnorm(1,mean=0,sd=sigma)
              if (edges[[i]][j,v] > 0) break
            }
          }
        }
        
        #--------------correlation function for drift ----------------
        mat_cor=diag(3)
        mat_cor[lower.tri(mat_cor)]= edges[[i]][1,((2+nb_traits):(1+nb_traits+nb_correl))]
        mat_cor[upper.tri(mat_cor)]= edges[[i]][1,((2+nb_traits):(1+nb_traits+nb_correl))]
         edges[[i]][-1,(2+nb_traits):(1+nb_traits+nb_correl)]= drift_correl(cor= mat_cor,len=phy$edge.length[i],sd=0.001) 
       } 

      mimi=apply(edges[[i]][,2:(1+nb_traits+nb_correl)],2,"min")
      min.trait[which(mimi<min.trait)]=mimi[which(mimi<min.trait)]
      mama=apply(edges[[i]][,2:(1+nb_traits+nb_correl)],2,"max")
      max.trait[which(mama>max.trait)]=mama[which(mama>max.trait)]
      phy$node.traits[phy$edge[i,2],] = edges[[i]][nrow(edges[[i]]),2:(1+nb_traits+nb_correl)]
    
    }
    cor_edge=edges
    edge2=rapply(edges,function(x) x[,-1],how = "list")
    cov_edge=rapply(edge2,function(x) t(apply(x,1,cor2cov2)),how = "list")
    results=list(cor_edge=cor_edge,cov_edge=cov_edge,min.trait=min.trait,max.trait=max.trait)
    simu[[nbsim]]=results
    }
  
    Simu_arraylist=simplify2array(simu, higher = TRUE)
    Simu_array=apply(Simu_arraylist,1,function(x) simplify2array(x, higher = TRUE))
    sim_cor_edge=apply(Simu_array$cor_edge,1,function(y) array(unlist(y), dim = c(dim(y[[1]]), length(y))))
    sim_cov_edge=apply(Simu_array$cov_edge,1,function(y) array(unlist(y), dim = c(dim(y[[1]]), length(y))))
    
    simu_results=list(sim_cor_edge=sim_cor_edge,sim_cov_edge=sim_cov_edge,min.trait=Simu_array$min.trait,max.trait=Simu_array$max.trait,phy=phy)
    return(simu_results)
}






a$sim_cov_edge[]min(phy$edge.length)

edges=a$cov_edge
phy=a$phy
max.ht = max(vcv.phylo(a$phy))
min.trait=a$min.trait
max.trait=a$max.trait

    par(mfcol=c(4,2),cex.lab=1,cex.axis=1,mar=c(5.1,5.1,1.1,1.1))
    plot.phylo(phy,show.tip.label=TRUE)
    
    i=1
    for (u in 1:6){
      plot(c(0,max.ht),c(min.trait[u],max.trait[u]),type='n',xlab='time',ylab='trait value')
      for (i in 1:phy$Nedge) lines(edges[[i]][,c(1,(u+1))])
    }
    
    
    
    
    
    
    a=VCVevol_brownian(G0,tree, 0.01,0.01)
    
    
    
    
    
    
    #####################################"" TEST OISEAU
    
    setwd("/media/mnhn/Leca/Gmatrix_project/3_Simulations")
    # phy_pass=read.nexus.data("Passeriformes.nex")
    phy_pass=read.nexus("Furnariidae.nex")
    phy_pass=compute.brlen(phy_pass, method="Grafen")
    
    phy_pass$edge.length*100
    
    
    plot(phy_pass)
    phy_pass
    
    plot(phy_pass)
    
    G0=matrix(c(2,0.1,1,
                0.1,2,1.8,
                1,1.8,2.2),3)
    
    a=VCVevol_brownian(G0,phy_pass,0.01,0.01,nb_repeat=100)
    
    a$sim_cov_edge
    