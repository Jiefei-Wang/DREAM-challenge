#install.packages("Rfast")
#install.packages("rpgm")
num_parm=11

normalize_scale<-function(mydata){
  mydata$N_insitu=t(scale(t(scale(mydata$insitu))))
  mydata$N_drop=scale(t(scale(t(mydata$drop))))
  return(mydata)
}

normalize_GIWC <- function(mydata){
  mydata$N_drop <- sweep(mydata$drop,2,colSums(mydata$drop),"/")
  mydata$N_insitu <- mydata$insitu
  return(mydata)
}


normalize_LogGIAC <- function(mydata){
  mydata$N_drop <- log2(1+sweep(mydata$drop,1,apply(mydata$drop,1,max),"/"))
  mydata$N_insitu <- log2(1+mydata$insitu)
  return(mydata)
}

normalize_GIAC <- function(mydata){
  mydata$N_drop <- sweep(mydata$drop,1,apply(mydata$drop,1,max),"/")
  mydata$N_insitu <- mydata$insitu
  #sweep(mydata$insitu,2,apply(mydata$insitu,2,max),"/")
  return(mydata)
}


normalize_TMM <- function(mydata){
  normfac_drop <- calcNormFactors(mydata$drop)
  mydata$N_drop<- sweep(mydata$drop,2,normfac_drop,"/")
  normfac_insitu <- calcNormFactors(t(mydata$insitu))
  mydata$N_insitu <-sweep(mydata$insitu,1,normfac_insitu,"/")
  return(mydata)
}



normalize_upqu <- function(mydata){
  parm=mydata$N_parm
  #insitu
  quantileExpressed <- apply(mydata$insitu, 1, function(x){quantile(x[x>0], parm)})
  mydata$N_insitu <- sweep(mydata$insitu,1,quantileExpressed,"/") 
  #dropseq
  quantileExpressed <- apply(mydata$drop, 2, function(x){quantile(x[x>0], parm)})
  mydata$N_drop <- sweep(mydata$drop,2,quantileExpressed,"/")
  return(mydata)  
}

get_parm_upqu <- function(param){
  n=num_parm
  return(as.list(seq(1,n-1)/n))
}

#binarized output
normalize_quantile<-function(mydata){
  parm=mydata$N_parm
  mydata$N_insitu<-apply(mydata$insitu,2,function(x){as.integer(x>quantile(x,parm))})
  mydata$N_drop<-apply(mydata$drop,1,function(x){as.integer(x>quantile(x,parm))})
  
  if(sum(abs(dim(mydata$N_drop)-dim(mydata$drop)))!=0)
    mydata$N_drop=t(mydata$N_drop)
  if(sum(abs(dim(mydata$N_insitu)-dim(mydata$insitu)))!=0)
    mydata$N_insitu=t(mydata$N_insitu)
  mydata
}

get_parm_quantile<-function(){
  n=num_parm
  return(as.list(seq(1,n-1)/n))
}

normalize_author_norm <- function(mydata){
  factor1 <- colSums(mydata$drop)/max(colSums(mydata$drop))
  factor2 <- rowSums(mydata$insitu)/max(rowSums(mydata$insitu))
  mydata$N_drop <- log(1+sweep(mydata$drop, 2, factor1, "/"))
  mydata$N_insitu <- log(1+sweep(mydata$insitu, 1, factor2, "/"))
  return(mydata)
}
#binarized output
normalize_author_quantile<-function(mydata){
  tmp=normalize_author_norm(mydata)
  tmp$drop=tmp$N_drop
  tmp$insitu=tmp$N_insitu
  mydata=normalize_quantile(tmp)
  return(mydata)
}

get_parm_author_quantile<-function(){
  get_parm_quantile()
}


normalize_double_quantile<-function(mydata){
  parm1=mydata$N_parm[1]
  parm2=mydata$N_parm[2]
  mydata$N_insitu<-apply(mydata$insitu,2,function(x){as.integer(x>quantile(x,parm1))})
  mydata$N_drop<-apply(mydata$drop,1,function(x){as.integer(x>quantile(x,parm2))})
  
  if(sum(abs(dim(mydata$N_drop)-dim(mydata$drop)))!=0)
    mydata$N_drop=t(mydata$N_drop)
  if(sum(abs(dim(mydata$N_insitu)-dim(mydata$insitu)))!=0)
    mydata$N_insitu=t(mydata$N_insitu)
  mydata
}

get_parm_double_quantile<-function(){
  n=num_parm
  #parm1=seq(1,n-1)/n
  #parm2=seq(1,n-1)/n
  
  parm1=seq(0.5,0.7,length.out = n)
  parm2=seq(0.8,1,length.out = n)
  
  parm=as.matrix(expand.grid(parm1,parm2))
  parm1=list()
  for(i in 1:nrow(parm)){
    parm1[[i]]=parm[i,]
  }
  return(parm1)
}



normalize_double_upqu<-function(mydata){
  parm1=mydata$N_parm[1]
  parm2=mydata$N_parm[2]
  #insitu
  quantileExpressed <- apply(mydata$insitu, 1, function(x){quantile(x[x>0], parm1)})
  mydata$N_insitu <- sweep(mydata$insitu,1,quantileExpressed,"/") 
  #dropseq
  quantileExpressed <- apply(mydata$drop, 2, function(x){quantile(x[x>0], parm2)})
  mydata$N_drop <- sweep(mydata$drop,2,quantileExpressed,"/")
  return(mydata) 
}

normalize_autoCluster<-function(mydata){
  mydata$N_insitu=apply(mydata$insitu,2,function(x){
    fit=kmeans(x,2)
    gmean=aggregate(x,by=list(fit$cluster),FUN=mean)
    group=fit$cluster-1
    if(gmean[1,2]>gmean[2,2]) 
    {
      group=1-group
    }
    group
  })
  
  mydata$N_drop=t(
    sapply(1:mydata$refNum,function(i,insitu,drop){
    sigNum=mean(insitu[,i])
    drop_data=drop[i,]
    drop_data>quantile(drop_data,sigNum)
  },insitu=mydata$N_insitu,drop=mydata$drop)
  )
  
  mydata
}


normalize_doubleAutoCluster<-function(mydata){
  mydata$N_insitu=apply(mydata$insitu,2,function(x){
    fit=kmeans(x,2)
    gmean=aggregate(x,by=list(fit$cluster),FUN=mean)
    group=fit$cluster-1
    if(gmean[1,2]>gmean[2,2]) 
    {
      group=1-group
    }
    group
  })
  mydata$N_drop=apply(mydata$drop[1:mydata$refNum,],1,function(x){
    fit=kmeans(x,2)
    gmean=aggregate(x,by=list(fit$cluster),FUN=mean)
    group=fit$cluster-1
    if(gmean[1,2]>gmean[2,2]) 
    {
      group=1-group
    }
    group
  })
  
  
  mydata
}


get_parm_double_upqu<-function(){
  get_parm_double_quantile()
}

