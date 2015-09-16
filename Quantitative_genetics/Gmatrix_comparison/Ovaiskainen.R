Ddivergence_ova<-function(names,nb_Gmatrix){
  nb_Gmatrix=nb_Gmatrix
  names=as.character(names)
 
  ############## Packages requirement
  if (require(mvtnorm) == FALSE)
    stop("mvtnorm not loaded")
  if (require(coda) == FALSE)
    stop("coda not loaded")
  if (require(MCMCglmm) == FALSE)
    stop("coda not loaded")
  
  ############## load first Gmatrix to assess the number of traits and the number of iterations
  model=get(paste(names,1,sep=''))
 # model=model$VCV
 # model=model[,agrep (".animal",colnames(model))]
  nb_trait=sqrt(dim(model)[2])
  nb_iter=dim(model)[1]
  
  ##############Load Gmatrices into an array
  mod=array(,c(nb_iter,nb_trait*nb_trait,nb_Gmatrix)) 
  mod[,,1]=model
  for (i in 2:nb_Gmatrix){ 
    model=get(paste(names,i,sep=''))
   # model=model$VCV
   # mod[,,i]=model[,agrep (".animal",colnames(model))]
    mod[,,i]=model
  }
  sample_mod=mod[sample(nb_iter),,]  #randomize array
  
  ################# Function to get distribtion of trait and distance between densities
  Ddivergence<-function (z = NULL, u = NULL, n = 10000) {
    xi <- rmvnorm(n, rep(0, dim(z)[1]),z)
    fx <- dmvnorm(xi, rep(0, dim(z)[1]),z)
    gx <- dmvnorm(xi, rep(0, dim(z)[1]),u)
    mean(sqrt(0.5 * ((fx - gx)^2)/(fx + gx)))
  }
  
  ################ Applying function between the array of Gmatrix and the randomized array
  dist=array(,c(nb_Gmatrix,nb_Gmatrix,nb_iter)) 

  
    for (i in 1:nb_iter){
      mo=array(mod[i,,],c(nb_trait,nb_trait,nb_Gmatrix))
      mo_sa=array(sample_mod[i,,],c(nb_trait,nb_trait,nb_Gmatrix))
      for (j in 1:nb_Gmatrix){
        dist[,j,i] =apply(mo,3,function(x,y) {Ddivergence(x,y)},y=mo_sa[,,j]) 
      }    
    }

  
  ############## Comparison between intra and intra variation of D
  if(nb_Gmatrix==2){
  Ddiv=mcmc(apply(dist,3, function(x) combn(diag(x), 2,sum,simplify=T)-((x[lower.tri(x)]+t(x)[lower.tri(x)]))))
  Ddist=mcmc(apply(dist,3,function(x) x[lower.tri(x)]))
  difference_summa=cbind(posterior.mode(Ddiv),HPDinterval(Ddiv))
  Pvalue=length(Ddiv[Ddiv>0])/nb_iter#}
  }else{
  Ddiv=mcmc(t(apply(dist,3, function(x) combn(diag(x), 2,sum,simplify=T)-((x[lower.tri(x)]+t(x)[lower.tri(x)])))))
  Ddist=mcmc(t(apply(dist,3,function(x) x[lower.tri(x)])))
  difference_summa=cbind(posterior.mode(Ddiv),HPDinterval(Ddiv)) 
  }
  Ddistance_summa=cbind(posterior.mode(Ddist),HPDinterval(Ddist))
  names_2by2=apply(do.call(rbind,combn(1:(nb_Gmatrix), 2,simplify=F)),1, paste, collapse="_")
  rownames(difference_summa)=names_2by2
  rownames(Ddistance_summa)=names_2by2
  names(Pvalue)=names_2by2
  results<-list(distance=Ddistance_summa,difference=difference_summa,pvalue=Pvalue)
  return(results)
}
 
