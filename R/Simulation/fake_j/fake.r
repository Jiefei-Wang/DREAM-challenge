getParms<-function(dispersion_x=0.8,dispersion_z=0.3,signal=1){
  parms=list(dispersion_x=dispersion_x,dispersion_z=dispersion_z,signal=signal,
             transf=function(x)log((x+1)),maxCount=100)
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
createPattern<-function(geometry,geneNum){
  refTable=matrix(0,length(geometry$x),geneNum)
  dropTable=matrix(0,length(geometry$x),geneNum)
  for(i in 1:geneNum){
  parms=getParms(signal = 1)
  refPattern=createRefPattern_circle(geometry,parms)
  dropPattern=createdropSeqPattern(refPattern,parms)
  refTable[,i]=refPattern
  dropTable[,i]=dropPattern
  }
pattern=list(refTable=refTable,dropTable=dropTable)
return(pattern)
}

patternData=createPattern(geometry,10)

k=2
intensityPlot(patternData$refTable[,k],geometry)


library(GenomicRanges)





