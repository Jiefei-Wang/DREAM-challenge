
normalize<-function(model,geneData){
  mydata=c(geneData,model)
  mydata$normalize(mydata)
}
compute_dist<-function(mydata){
  mydata$computeDist(mydata)
}

pred_loc<-function(mydata){
  predict_num=10
  #Compute the rank
  rank_distance=apply(-mydata$distance,2,rank,ties.method="random")
  ind=which(rank_distance<=predict_num,arr.ind=T)
  sort_order=apply(matrix(rank_distance[ind],nrow=predict_num),2,order)
  loc=matrix(ind[,1],nrow=predict_num)
  mydata$loc=sapply(1:ncol(loc),function(i,loc,ind){
    loc[ind[,i],i]
  },loc=loc,ind=sort_order)
  return(mydata)
}

predict_pattern<-function(mydata,gene.start=1,gene.end=mydata$geneNum){
  for(i in gene.start:gene.end){
    mydata=mydata$compute_pattern(mydata,i)
  }
  mydata
}
