
getParms<-function(dispersion_x=0.8,dispersion_z=0.3,signal=1){
  parms=list(dispersion_x=dispersion_x,dispersion_z=dispersion_z,signal=signal,refGeneSignal=0.01,
             transf=function(x)log((x+1)),maxCount=100,
             speciesDiff=0.01,
             corInflation=1,varInflation_insitu=0.01,varInflation_drop=2,blockNum=40)
  return(parms)
}

#Create the true expression pattern for a gene, each row of geometry has one expression number
createTruePattern_circle<-function(geometry,parms){
#large dispersion means large range
dispersion_x=parms$dispersion_x
dispersion_z=parms$dispersion_z
signal=parms$signal
x=geometry$x
z=geometry$z
anchor_point=c(sample(x,1),sample(z,1))
scale_distance=1/(sqrt(c(var(x),var(z)))*c(dispersion_x,dispersion_z))
if(is.nan(sum(scale_distance))) return(rep(0,length(x)))
intensity=rep(0,length(x))
for(i in 1:length(x)){
  distance=(anchor_point-c(x[i],z[i]))*scale_distance
  distance=sum(distance^2)
  weight=pnorm(distance,sd=signal)
  intensity[i]=min(weight,1-weight)
}
intensity=(intensity-min(intensity))/(max(intensity)-min(intensity))
if(is.nan(sum(intensity))) return(rep(0,length(x)))
return(intensity)
}

#Create the true gene count using the true gene expression pattern
#using a tranfermation to mimic the systematic bias of the sequencing technique
createdropSeqPattern<-function(truePattern,parms){
  count=round(parms$transf(truePattern)*100)
  return(count)
}

#Create the true expression pattern for all genes
#geneNum represents the number of genes in drop-seq
#refNum represents the number of genes in distmap
#The first refNum column is the reference gene
#For the reference gene, adding more variation to make it unique in every location
createPattern<-function(geometry,geneNum,refNum,parms){
  trueTable=matrix(0,length(geometry$x),geneNum)
  dropTable=matrix(0,length(geometry$x),geneNum)
  parmsList=c()
  for(i in 1:geneNum){
  if(i>refNum){
    parms$dispersion_x=parms$dispersion_x/4
    parms$dispersion_y=parms$dispersion_y/4
  } 
  truePattern=createTruePattern_circle(geometry,parms)
  #if it is reference gene, add variation
  if(i<=refNum){
    truePattern=truePattern*2+runif(nrow(geometry),-parms$refGeneSignal,parms$refGeneSignal)
    truePattern[truePattern<0]=0
    truePattern[truePattern>1]=1
  }
  dropPattern=createdropSeqPattern(truePattern,parms)
  trueTable[,i]=truePattern
  dropTable[,i]=dropPattern
  parmsList[i]=0
  }
  gene.name=paste("Gene",1:ncol(trueTable))
  colnames(trueTable)=gene.name
  colnames(dropTable)=gene.name
pattern=list(trueTable=trueTable,dropTable=dropTable,parmsList=parmsList)
return(pattern)
}
#Adding the species effect
#(the species measured in insitu data is not necessarily the species measured in drop-seq)
addSpeciesDiff<-function(patternData,parms){
  refTable=patternData$trueTable
  speciesDiff=matrix(rnorm(length(refTable),sd=parms$speciesDiff),nrow(refTable),ncol(refTable))
  patternData$insituTable=patternData$trueTable+speciesDiff
  patternData$insituTable[patternData$insituTable<0]=0
  return(patternData)
}
#Compute the covariance matrix(Drop-seq only)
computeCov<-function(meanPattern,parms){
  set.seed(sum(meanPattern))
  variance=runif(length(meanPattern),0.5,1.5)
  correlation=matrix(rnorm(length(meanPattern)^2,0,1*parms$corInflation),length(meanPattern),length(meanPattern))
  for(i in 1:length(meanPattern)){
    correlation[i,]=correlation[i,]/sqrt(sum(correlation[i,]*correlation[i,]))
  }
  correlation <- correlation%*%t(correlation)
  if(length(variance)!=1)
    covariance=sqrt(diag(variance))%*%correlation%*%sqrt(diag(variance))
  else
    covariance=correlation*correlation
  set.seed(as.numeric(Sys.time()))
  return(correlation)
}
#simulate the distMap data
sampleInsituData<-function(patternData,refNum,parms){
  insituMean=patternData$insituTable[,1:refNum]
  insituData=insituMean+matrix(rnorm(length(insituMean),0,parms$varInflation_insitu),nrow(insituMean),ncol(insituMean))
  return(insituData)
}
#simulate drop-seq data(With some degree of correlation between genes)
sampleDropData<-function(patternData,cellNum,parms){
  dropMean=patternData$dropTable
  ind=sample(1:nrow(dropMean),cellNum,replace = T)
  dropMean=dropMean[ind,]
  #dropData=matrix(0,nrow(dropMean),ncol(dropMean))
  geneNum=ncol(dropMean)
  packageList=c("mvtnorm")
  #for(i in 1:cellNum){
  sim=foreach(i=1:cellNum,.combine= cbind,.multicombine=TRUE,.inorder=FALSE,.packages =packageList,.export=c("computeCov"))%dopar%{
    meanPattern=dropMean[i,]
    set.seed(sum(meanPattern))
    block=rmultinom(1,geneNum,runif(parms$blockNum,0,1))
    measure=rep(0,geneNum)
    k=1
    for(j in 1:parms$blockNum){
      if(block[j]==0)
        next
      meanBlock=meanPattern[k:(k+block[j]-1)]
      covPattern=computeCov(meanBlock,parms)*parms$varInflation_drop
      measure_block=meanBlock+rmvnorm(n=1,rep(0,block[j]),covPattern)
      measure[k:(k+block[j]-1)]=measure_block
      k=k+block[j]
    }
    measure[measure<0]=0
    measure
    #dropData[,i]=measure
  }
  dropData=sim
  set.seed(as.numeric(Sys.time()))
  dropData=round(dropData)
  row.names(dropData)=colnames(patternData$dropTable)
  colnames(dropData)=paste("Cell",1:ncol(dropData))
  result=list(dropData=dropData,location=ind)
  return(result)
}



