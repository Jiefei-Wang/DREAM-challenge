#Possibly the location variance information can be used?
computePattern_simple<-function(mydata,gene){
  geneExp=mydata$drop[gene,]
  tmp=mydata$distance
  index=(tmp>quantile(tmp,0.1))
  pattern=sweep(1/(tmp+0.0001),2,geneExp,FUN = '*')
  pattern[index]=0
  rowSums(pattern)
}


computePattern_simple1<-function(mydata,gene){
  n=10
  geneExp=mydata$drop[gene,]
  tmp=mydata$distance
  indice=apply(tmp,2,function(x,n){rank(x,ties.method="random")<=n},n=n)
  ind=which(indice,arr.ind = T)
  dist=matrix(9999999,nrow(tmp),ncol(tmp))
  dist[ind]=tmp[ind]
  pattern=sweep(1/(dist+0.0001),2,geneExp,FUN = '*')
  rowSums(pattern)
}


computePattern_author<-function(mydata, gene) {
  threshold=mydata$p_parm
  gene.expr <- mydata$drop[gene,]
  m=max(mydata$distance)
  b1 <- sweep(m-mydata$distance, 2, gene.expr, '*')
  gene.expr[gene.expr > 0] <- 1
  b2 <- sweep(m-mydata$distance, 2, gene.expr, '*')
  q <- rowSums(b1)/rowSums(b2)
  q[is.na(q)] <- 0
  q <- q/(1+q)
  q[q < quantile(q, threshold)] <- 0
  return (q)
}

get_parm_author<-function(){
  n=10
  #return(as.list(seq(1,n-1)/n))
  return(list(0.5,0.6,0.7,0.8))
}
