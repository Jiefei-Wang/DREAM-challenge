
nameConvert<-function(geneName,target){
  if(target=="drop"){
    geneName[geneName=="Blimp.1"]="Blimp-1"
    geneName[geneName=="E.spl.m5.HLH"]="E(spl)m5-HLH"
  }
  if(target=="insitu"){
    geneName[geneName=="Blimp-1"]="Blimp.1"
    geneName[geneName=="E(spl)m5-HLH"]="E.spl.m5.HLH"
  }
  return(geneName)
}

weighted_cor<-function(x,y,w){
  x_m=sum(w*x)/sum(w)
  y_m=colsums(sweep(y,1,w,"*"))/sum(w)
  x_tmp=x-x_m
  y_tmp=sweep(y,2,y_m,"-")
  num=colsums(sweep(y_tmp,1,w*x_tmp,"*"))
  den=sqrt(sum(w*x_tmp^2)*colsums(sweep(y_tmp^2,1,w,"*")))
  num/den
}

pickGene<-function(geneData,geneNum=NULL,ind=NULL){
  if(is.null(ind))
    ind=sample(1:84,geneNum)
  ind=sort(ind)
  if(is.numeric(ind[1])){
    geneName=rownames(geneData$drop)[1:84]
    geneName=geneName[ind]
  }else{
    geneName=ind
  }
  
  geneData$insitu=geneData$insitu[,geneName]
  
  index=which(rownames(geneData$drop)%in%geneName)
  drop_ref=geneData$drop[index,]
  drop_nonref=geneData$drop[-index,]
  geneData$drop=rbind(drop_ref,drop_nonref)
  
  geneNum=length(ind)
  geneData$refNum=geneNum
  geneData$geneNum=nrow(geneData$drop)
  geneData$cellNum=ncol(geneData$cellNum)
  geneData$geneInd=ind
  geneData$geneName=geneName
  geneData
}
