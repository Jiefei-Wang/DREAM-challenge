
library(DistMap)
intensityPlot<-function(gene_level,geometry,title="",side="xz",col.range=c(3000,10000)){
  if(title=="")
    title=as.character(substitute(gene_level))
  col.lowest=col.range[1]
  col.highest=col.range[2]
  mydata=data.frame(x=geometry$x,y=geometry$y,z=geometry$z,gene=gene_level)
  mydata$gene=rank(mydata$gene)
  mydata$gene=round((mydata$gene-min(mydata$gene))*(col.highest-col.lowest-1)/(max(mydata$gene)-min(mydata$gene))+1)
  data.color=rev(heat.colors(col.highest))[col.lowest:col.highest]
  mydata$Col <- data.color[mydata$gene]
  if(side=="xz")
    plot(mydata$x,mydata$z,pch = 20,col = mydata$Col,xlab="x",ylab="z",main=title)
  if(side=="zy")
    plot(mydata$x,mydata$y,pch = 20,col = mydata$Col,xlab="x",ylab="y",main=title)
}
#gene_level=simulation$patternData$dropTable[,gene]
intensityPlot2<-function(gene_level,geometry,cut_quantile=0.9,title="",side="xz",col.range=c(80,256)){
  if(title=="")
    title=as.character(substitute(gene_level))
  gene_level[gene_level<quantile(gene_level,cut_quantile)]=quantile(gene_level,cut_quantile)
  col.lowest=col.range[1]
  col.lowest05=col.range[1]+(col.range[2]-col.range[1])*0.05
  col.highest=col.range[2]
  if(var(gene_level)!=0){
    gene_level=scale(gene_level)
    gene_level_min=min(gene_level)
    gene_level_max=max(gene_level)
    gene.col=col.lowest05-col.lowest+(gene_level-gene_level_min)/(gene_level_max-gene_level_min)*(col.highest-col.lowest05)
    gene.col=round(gene.col)
    gene.col[gene.col<1]=1
    gene.col[gene.col>col.highest-col.lowest]=col.highest-col.lowest
  }else
    gene.col=rep(1,length(gene_level))
  
  data.color=rev(heat.colors(col.highest))[col.lowest:col.highest]
  if(side=="xz")
    plot(geometry$x,geometry$z,pch = 20,col = data.color[gene.col],xlab="x",ylab="z",main=title)
  if(side=="zy")
    plot(geometry$x,geometry$y,pch = 20,col = data.color[gene.col],xlab="x",ylab="y",main=title)
}
intensityPlot3<-function(gene_level,geometry,title="",side="xz",col.range=c(3000,10000)){
  if(title=="")
    title=as.character(substitute(gene_level))
  col.lowest=col.range[1]
  col.highest=col.range[2]
  mydata=data.frame(x=geometry$x,y=geometry$y,z=geometry$z,gene=gene_level)
  mydata$gene=round((mydata$gene-min(mydata$gene))*(col.highest-col.lowest-1)/(max(mydata$gene)-min(mydata$gene))+1)
  data.color=rev(heat.colors(col.highest))[col.lowest:col.highest]
  mydata$Col <- data.color[mydata$gene]
  if(side=="xz")
    plot(mydata$x,mydata$z,pch = 20,col = mydata$Col,xlab="x",ylab="z",main=title)
  if(side=="zy")
    plot(mydata$x,mydata$y,pch = 20,col = mydata$Col,xlab="x",ylab="y",main=title)
}

originalMethod<-function(mydata,simulation){
  sim_dropSeq_normalized=t(apply(apply(mydata$drop, 2, scale,center=F), 1, scale,center=F))
  row.names(sim_dropSeq_normalized)=row.names(mydata$drop)
  colnames(sim_dropSeq_normalized)=colnames(mydata$drop)
  dm = new("DistMap",
           raw.data=as.matrix(mydata$drop),
           data=as.matrix(sim_dropSeq_normalized),
           insitu.matrix=as.matrix(mydata$insitu>0.1),
           geometry=as.matrix(geometry))
  dm <- binarizeSingleCellData(dm, seq(0.15, 0.5, 0.01))
  dm <- mapCells(dm)
  dm_res=sapply(1:geneNum,computeVISH,object=dm,threshold=0.75)
  dm_pattern_score=patternScore(dm_res,simulation$patternData$dropTable,mydata$refNum+1,mydata$geneNum)

  dm_loc=predLocation(list(distance=-dm@mcc.scores))
  dm_pred_score=predictionScore(dm_loc$loc,simulation$cell_loc)
  dm_pred_score
  list(pattern=dm_res,score=data.frame(normalization="scale",
             distance="Authors",
             pattern="Authors",
             pattern_score=mean(dm_pattern_score),prediction_score=dm_pred_score))
}

# 
# sweep2<-function(X,MARGIN,STATS,FUN){
#   if(FUN=="*"||FUN=="/"||FUN=="+"||FUN=="-"){
#     if(MARGIN==1){
#       Exp=paste0("X",FUN,"STATS")
#       return(eval(parse(text=Exp)))
#     }
#     if(MARGIN==2){
#       Exp=paste0("t(t(X)",FUN,"STATS)")
#       return(eval(parse(text=Exp)))
#     }
#   }
#   sweep(X,MARGIN,STATS,FUN)
# }


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

