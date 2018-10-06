computeDist_cov<-function(mydata){
  mydata$distance=1-abs(cor(t(mydata$N_insitu),as.matrix(mydata$N_drop[1:ncol(mydata$N_insitu),])))
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
