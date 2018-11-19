#Possibly the location variance information can be used?
computePattern_simple<-function(mydata,gene){
  if(is.null(mydata$tmp_pattern)){
  tmp=mydata$distance
  index=(tmp>quantile(tmp,0.1))
  mydata$tmp_pattern=1/(tmp+0.0001)
  }
  geneExp=mydata$drop[gene,]
  pattern=sweep(mydata$tmp_pattern,2,geneExp,FUN = '*')
  pattern[index]=0
  mydata$pattern=cbind(mydata$pattern,matrix(rowSums(pattern),ncol=1))
}

computePattern_simple1<-function(mydata,gene){
  if(is.null(mydata$tmp_pattern)){
    n=mydata$p_parm
    p=n/nrow(mydata$distance)
    tmp=mydata$distance
    #indice=apply(tmp,2,function(x,n){rank(x,ties.method="random")<=n},n=n)
    indice=apply(tmp,2,function(x,n,p){x<=quantile(x,p)},n=n,p=p)
    ind=which(indice,arr.ind = T)
    dist=matrix(9999999,nrow(tmp),ncol(tmp))
    dist[ind]=tmp[ind]
    dist=1/(dist+0.001)
    dist=sweep(dist,2,colSums(dist),"/")
    mydata$tmp_pattern=dist
  }
  geneExp=mydata$drop[gene,]
  pattern=sweep(mydata$tmp_pattern,2,geneExp,FUN = '*')
  #geneExp2=geneExp
  #geneExp2[geneExp2>0]=1
  #pattern2=sweep(mydata$tmp_pattern,2,geneExp2,FUN = '*')
  
  
  mydata$pattern=cbind(mydata[["pattern"]],matrix(rowSums(pattern),ncol=1))
  mydata
}
get_parm_simple1<-function(){

  return(list(10))
}


computePattern_author<-function(mydata, gene) {
  threshold=mydata$p_parm
  gene.expr <- mydata$drop[gene,]
  m=max(mydata$distance)
  similarity=t(m-mydata$distance)
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

get_parm_author<-function(){
  n=10
  #return(as.list(seq(1,n-1)/n))
  return(list(0.6))
}
