computeDist_cov<-function(mydata){
  mydata$distance=1-abs(cor(t(mydata$N_insitu),as.matrix(mydata$N_drop[1:ncol(mydata$N_insitu),])))
  mydata$distance[is.na(mydata$distance)]=runif(sum(is.na(mydata$distance)))+1
  return(mydata)
}

computeDist_cov_rank<-function(mydata){
  mydata$distance=1-abs(cor(t(mydata$N_insitu),as.matrix(mydata$N_drop[1:ncol(mydata$N_insitu),]),method="kendall"))
  mydata$distance[is.na(mydata$distance)]=runif(sum(is.na(mydata$distance)))+1
  return(mydata)
}
#mean squared error scheme
computeDist_mse <- function(mydata){
  drop=mydata$N_drop[1:mydata$refNum,]
  insitu=t(mydata$N_insitu)
  mydata$distance =apply(drop,2,function(x,insitu){
    colMeans(sweep(insitu,1,x,"-")^2)
  },insitu=insitu)
  return(mydata)
}

computeDist_mcc <-function(mydata){
  drop=mydata$N_drop[1:mydata$refNum,]
  insitu=t(mydata$N_insitu)
  mcc = matrix(nrow = nrow(mydata$insitu),ncol = ncol(drop))
  for(i in 1:nrow(mydata$insitu)){
    mccmtx<- sweep(2*drop, 1,insitu[,i], '-')
    FP = colSums(mccmtx == 2)
    TP = colSums(mccmtx == 1)
    TN = colSums(mccmtx == 0)
    FN = colSums(mccmtx == -1)
    mcc[i,] <- ((TP*TN - FP*FN)/max(0.0001,sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))))
  }
  mydata$distance=max(mcc)-mcc
  return (mydata)
  
}