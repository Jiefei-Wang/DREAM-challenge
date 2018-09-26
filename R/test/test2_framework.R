
source("R\\commonFunc\\readData.R")
source("R\\commonFunc\\functions.R")
source("R\\framework.R")
source("R\\generate_fake.R")
geometry=chanllege.data$geometry
refNum=50
geneNum=100
cellNum=1000
simulation=generate_fake1(geometry,refNum,geneNum,cellNum)
mydata=simulation$mydata

normalize_simple<-function(mydata){
  mydata$N_insitu=t(scale(t(scale(mydata$insitu))))
  mydata$N_drop=t(scale(t(scale(mydata$drop))))
  return(mydata)
}
dist_simple<-function(mydata){
  mydata$distance=cor(t(mydata$N_insitu),as.matrix(mydata$N_drop[1:ncol(mydata$N_insitu),]))
  mydata$distance[is.na(mydata$distance)]=runif(sum(is.na(mydata$distance)))+999999
  return(mydata)
}

predLocation_simple<-function(mydata){
  predict_num=1000
  #Compute the rank
  rank_distance=apply(mydata$distance,2,rank,ties.method="random")
  ind=which(rank_distance<=predict_num,arr.ind=T)
  sort_order=apply(matrix(rank_distance[ind],nrow=predict_num),2,order)
  loc=matrix(ind[,1],nrow=predict_num)
  mydata$loc=sapply(1:ncol(loc),function(i,loc,ind){
    loc[ind[,i],i]
  },loc=loc,ind=sort_order)
  return(mydata)
}
#Possibly the location variance information can be used?
computePattern_simple<-function(mydata,gene){
  geneExp=mydata$drop[gene,]
  loc=1:nrow(mydata$insitu)
  similarity=t(sapply(loc,function(x,loc){
    ind=which(loc==x,arr.ind=T)
    ind[,1]=(nrow(loc)-ind[,1]+1)/10
    score=rep(0,ncol(loc))
    score[ind[,2]]=ind[,1]
    score
    },mydata$loc))
  pattern=similarity%*%geneExp
  pattern
}




gene=1
mydata=normalize_simple(mydata)
mydata=dist_simple(mydata)
mydata=predLocation_simple(mydata)
res_pred=computePattern_simple(mydata,gene)

intensityPlot2(simulation$patternData$dropTable[,gene],geometry)
intensityPlot2(res_pred,geometry)

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

mydata$N_insitu=dm@insitu.matrix
mydata$N_drop=dm@data

mydata$distance=dm@mcc.scores


mydata=attachFunc(mydata,normalize,computeDist)

mydata=normalize(mydata)
mydata=compute_dist(mydata)
mydata=pred_loc(mydata)
mydata=compute_pattern(mydata)

a
rank(-a)
