#Possibly the location variance information can be used?
computePattern_simple<-function(mydata,gene){
  geneExp=mydata$N_drop[gene,]
  tmp=mydata$distance
  index=(tmp>quantile(tmp,0.1))
  pattern=sweep(1/(tmp+0.0001),2,geneExp,FUN = '*')
  pattern[index]=0
  rowSums(pattern)
}


computePattern_simple1<-function(mydata,gene){
  n=2
  geneExp=mydata$N_drop[gene,]
  tmp=mydata$distance
  indice=apply(tmp,2,function(x,n){rank(x,ties.method="random")<=n},n=n)
  ind=which(indice,arr.ind = T)
  dist=matrix(9999999,nrow(tmp),ncol(tmp))
  dist[ind]=tmp[ind]
  pattern=sweep(1/(dist+0.0001),2,geneExp,FUN = '*')
  rowSums(pattern)
}



