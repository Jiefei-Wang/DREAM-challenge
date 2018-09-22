library(mvtnorm)
getParms<-function(dispersion_x=0.8,dispersion_z=0.3,signal=1){
  parms=list(dispersion_x=dispersion_x,dispersion_z=dispersion_z,signal=signal,
             transf=function(x)log((x+1)),maxCount=100,
             speciesDiff=0.1,
             corInflation=1,varInflation_insitu=0.01,varInflation_drop=1,blockNum=40)
  return(parms)
}

createRefPattern_circle<-function(geometry,parms){
#large dispersion means large range
dispersion_x=parms$dispersion_x
dispersion_z=parms$dispersion_z
signal=parms$signal
x=geometry$x
z=geometry$z
anchor_point=c(sample(x,1),sample(z,1))
scale_distance=1/(sqrt(c(var(x),var(z)))*c(dispersion_x,dispersion_z))

intensity=rep(0,length(x))
for(i in 1:length(x)){
  distance=(anchor_point-c(x[i],z[i]))*scale_distance
  distance=sum(distance^2)
  weight=pnorm(distance,sd=signal)
  intensity[i]=min(weight,1-weight)
}
intensity=(intensity-min(intensity))/(max(intensity)-min(intensity))

return(intensity)
}


createdropSeqPattern<-function(refPattern,parms){
  count=round(parms$transf(refPattern)*100)
  return(count)
}
createPattern<-function(geometry,geneNum,refNum){
  refTable=matrix(0,length(geometry$x),geneNum)
  dropTable=matrix(0,length(geometry$x),geneNum)
  parmsList=list()
  for(i in 1:geneNum){
  parms=getParms(signal = 1)
  refPattern=createRefPattern_circle(geometry,parms)
  if(i<=refNum){
    refPattern=refPattern+runif(nrow(geometry),-0.1,0.5)
    refPattern[refPattern<0]=0
    refPattern[refPattern>1]=1
  }
  dropPattern=createdropSeqPattern(refPattern,parms)
  refTable[,i]=refPattern
  dropTable[,i]=dropPattern
  parmsList[[i]]=0
  }
pattern=list(refTable_true=refTable,dropTable=dropTable,parmsList=parmsList)
return(pattern)
}
addSpeciesDiff<-function(pattern,parms){
  refTable=patternData$refTable_true
  speciesDiff=matrix(rnorm(length(refTable),sd=parms$speciesDiff),nrow(refTable),ncol(refTable))
  patternData$refTable_species=patternData$refTable_true+speciesDiff
  patternData$refTable_species[patternData$refTable_species<0]=0
  return(patternData)
}
computeCov<-function(meanPattern,parms){
  set.seed(sum(meanPattern))
  variance=runif(length(meanPattern),0.5,1.5)
  correlation=matrix(rnorm(length(meanPattern)^2,0,1*parms$corInflation),length(meanPattern),length(meanPattern))
  for(i in 1:length(meanPattern)){
    correlation[i,]=correlation[i,]/sqrt(sum(correlation[i,]*correlation[i,]))
  }
  correlation <- correlation%*%t(correlation)
  #correlation[correlation>1]=1
  #correlation[correlation<-1]=-1
  covariance=sqrt(diag(variance))%*%correlation%*%sqrt(diag(variance))
  set.seed(as.numeric(Sys.time()))
  return(correlation)
}
sampleInsituData<-function(patternData,refNum,parms){
  insituMean=patternData$refTable_species[,1:refNum]
  insituData=insituMean+matrix(rnorm(length(insituMean),0,parms$varInflation_insitu),nrow(insituMean),ncol(insituMean))
  return(insituData)
}

sampleDropData<-function(patternData,cellNum,parms){
  dropMean=patternData$dropTable
  ind=sample(1:ncol(dropMean),cellNum,replace = T)
  dropMean=dropMean[ind,]
  dropData=matrix(0,nrow(dropMean),ncol(dropMean))
  geneNum=ncol(dropMean)
  for(i in 1:cellNum){
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
    dropData[i,]=measure
  }
  set.seed(as.numeric(Sys.time()))
  return(dropData)
}


sum(abs(measure-meanPattern))

source("R\\commonFunc\\readData.R")
source("R\\commonFunc\\functions.R")
refNum=10
geneNum=1000
cellNum=100
parms=getParms()
patternData=createPattern(geometry,geneNum,refNum)
patternData=addSpeciesDiff(patternData,parms)

sim_insitu=sampleInsituData(patternData,refNum,parms)
sim_drop=sampleDropData(patternData,cellNum,parms)



k=2
intensityPlot(patternData$refTable[,k],geometry)


a




