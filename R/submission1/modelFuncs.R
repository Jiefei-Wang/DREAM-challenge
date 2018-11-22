buildModel<-function(){
  myModel=list()
  myModel$normalize=normalize_dropClusterOnly_parm
  myModel$computeDist=computeDist_Weighted_mcc
  myModel$compute_pattern=computePattern_author
  
  N_parm=get_parm_dropClusterOnly_parm()
  #The parameter to compute the pattern
  #The number is from the author's code
  myModel$p_parm=0.6
  modelList=list()
  for(i in 1:length(N_parm)){
    myModel$N_parm=N_parm[[i]]
    modelList=c(modelList,list(myModel))
  }
  modelList
}



normalize_dropClusterOnly_parm<-function(mydata){
  parm=mydata$N_parm
  N_driver=84
  mydata$N_insitu=mydata$insitu
  tmp=mydata$drop[1:N_driver,]
  
  #quantileExpressed <- apply(tmp, 2, function(x){quantile(x[x>0], 0.95)})
  #tmp_norm <- sweep(tmp,2,quantileExpressed,"/")
  tmp_norm=tmp
  
  quantileCutoff=apply(tmp_norm,1,function(x){quantile(x[x>0],parm)})
  
  mydata$N_drop<-sapply(1:N_driver,function(ind,gene,cutoff){as.integer(gene[ind,]>cutoff[ind])},gene=tmp_norm,cutoff=quantileCutoff)
  
  if(sum(abs(dim(mydata$N_drop)-dim(mydata$drop)))!=0)
    mydata$N_drop=t(mydata$N_drop)
  mydata
}

get_parm_dropClusterOnly_parm<-function(){
  seq(0.15, 0.5, 0.01)
}


computeDist_Weighted_mcc <-function(mydata){
  all_drop=t(mydata$N_drop[1:84,])
  dropName=rownames(mydata$drop)[1:84]
  driverNum=mydata$refNum
  
  drop_cor=abs(cor(all_drop))
  
  weight=colsums(drop_cor[(driverNum+1):84,1:driverNum])+1
  weight=weight/sum(weight)
  
  drop=mydata$N_drop[1:driverNum,]
  #special treatment of 0 variance cells
  drop_var=apply(drop,2,var)
  ind=which(drop_var==0)
  for(i in ind){
    drop[sample(1:driverNum,2),i]=c(0,1)
  }
  
  insitu=t(mydata$N_insitu)
  mcc = matrix(nrow = nrow(mydata$insitu),ncol = ncol(drop))
  for(i in 1:nrow(mydata$insitu)){
    mcc[i,] <- weighted_cor(insitu[,i],drop,weight)
  }
  
  mydata$distance=mcc
  return (mydata)
}


computePattern_author<-function(mydata, gene) {
  threshold=mydata$p_parm
  gene.expr <- mydata$drop[gene,]
  similarity=t(mydata$distance-min(mydata$distance))
  b1 <- sweep(similarity, 1, gene.expr, '*')
  ind=(gene.expr > 0)
  b2=b1
  b2[ind,]=similarity[ind,]
  #gene.expr=gene.expr[gene.expr > 0] <- 1
  #b2 <- sweep(similarity, 1, gene.expr, '*')
  q <- colsums(b1)/colsums(b2)
  q[is.na(q)] <- 0
  q <- q/(1+q)
  q[q < quantile(q, threshold)] <- 0
  mydata$pattern=cbind(mydata$pattern,matrix(q,ncol=1))
  mydata
}

get_parm_author_mcc<-function(){
  return(list(0.6))
}
