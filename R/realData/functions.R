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

author_method<-function(chanllege.data,geneData){
  raw.data=geneData$drop
  gene.name=rownames(raw.data)
  data=as.matrix(chanllege.data$dropSeq.normalized[gene.name,])
  
  
  dm = new("DistMap",
           raw.data=raw.data,
           data=data,
           insitu.matrix=as.matrix(chanllege.data$insitu.binary[,colnames(geneData$insitu)]),
           geometry=as.matrix(chanllege.data$geometry))
  dm <- binarizeSingleCellData(dm, seq(0.15, 0.5, 0.01))
  dm <- mapCells(dm)
  dm
}

author_method_score<-function(dm,cell_loc){
  predLoc=predLocation(list(distance=-dm@mcc.scores))
  mean(predictionScore(predLoc$loc,cell_loc))
}

author_method_pattern<-function(dm){
  gene.name.list=colnames(dm@insitu.matrix)
  pattern_author=foreach(i=1:length(gene.name.list),.combine = cbind,.export ="author_pattern_method")%dopar%{
    gene.level.pred = author_pattern_method(dm, i, threshold=0.75)
    gene.level.pred
  }
  pattern_author
}

author_pattern_method<-function(dm,gene,threshold){
  gene.expr <- as.numeric(dm@data[gene, ])
  b1 <- sweep(dm@mcc.scores, 2, gene.expr, '*')
  gene.expr[gene.expr > 0] <- 1
  b2 <- sweep(dm@mcc.scores, 2, gene.expr, '*')
  q <- rowSums(b1,na.rm = T)/rowSums(b2,na.rm = T)
  q[is.na(q)] <- 0
  q <- q/(1+q)
  q[q < quantile(q, threshold)] <- 0
  q
}


findgeneInd<-function(geneName,geneList){
  for(i in 1:length(geneName)){
    if(!geneName[i]%in% geneList)
      message("This gene is not in the gene list: ",geneList[i])
  }
  (1:length(geneList))[geneList%in% geneName]
}


