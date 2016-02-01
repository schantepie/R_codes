###exemple

# name="gmat" 
#   nb_gmat=10
#   div_thin=2 #increase thining by div_thin factor 
#   traits=c("LD","CS","FL") #trait to use in the raw data to perform standardisation, the file should be names data
#   scale=c(10,1,10)#for example here, data LD and FL were multiplicated x10 in the row data
#   scale=NULL
#   mean_std=FALSE #standardize to the mean
#   var_std=TRUE

thining_scale<-function(name=NULL,nb_gmat=NULL,scale=NULL,traits=NULL,div_thin=NULL,mean_std=FALSE,var_std=FALSE){
  require(MCMCglmm)
  load(paste(1,name,".Rdata" , sep=""))
  data_reduced=data[,traits]
  GMAT_VCV_list=list()
  
  for (i in 1:nb_gmat){
    nomFichier=paste(i,name,".Rdata" , sep="")
    load(nomFichier)
    GMAT_VCV_list[[i]]= assign(paste("Gmat_",i,sep=""),Gmat$VCV)
  }
  GMAT_VCV_list_concat<- do.call(rbind,GMAT_VCV_list)
  
  if (!is.null(scale)){
    data_reduced=t(apply(data_reduced,1,function(x) {as.numeric(x)/(scale)}))
    GMAT_VCV_list_concat=mcmc(t(apply(GMAT_VCV_list_concat[,agrep (".animal",colnames(GMAT_VCV_list_concat))],1,function(x) {x/(scale^2)})))
  }
  
  gmat_thin=GMAT_VCV_list_concat[,agrep (".animal",colnames(GMAT_VCV_list_concat))]
  
  if (!is.null(div_thin))  {
    gmat_thining=dim(gmat_thin)[1]/div_thin
    gmat_thin= mcmc(gmat_thin[seq(1,gmat_thining*div_thin,div_thin),])
  }
  
  gmat_thin=mcmc(gmat_thin)
  autoc=autocorr(gmat_thin)
  
  if (mean_std==TRUE)  {
    mean_data=apply(data_reduced,2,function(x)mean(as.numeric(x[!is.na(x)])))
    tmean_data=t(mean_data)
    std_mean_data=mean_data%*%tmean_data
    colnames(std_mean_data)=colnames(data_reduced)
    rownames(std_mean_data)=colnames(data_reduced)
    std_mean_data=c(std_mean_data)
    gmat_thin=mcmc(t(apply(gmat_thin,1,function(x) {x/std_mean_data})))
  }
  
  if (var_std==TRUE)   {
    var_data=apply(data_reduced,2,function(x)var(as.numeric(x[!is.na(x)])))
    tvar_data=t(var_data)
    std_var_data=var_data%*%tvar_data
    colnames(std_var_data)=colnames(data_reduced)
    rownames(std_var_data)=colnames(data_reduced)
    std_var_data=c(std_var_data)
    gmat_thin=mcmc(t(apply(gmat_thin,1,function(x) {x/std_var_data})))
  }
  
  return(list(gmat_thin=gmat_thin,autocorr=autoc))
}
