
generate_fake<-function(geometry,refNum,geneNum,cellNum){
  # 
  # cl_func=NULL
  # if(exists("cl",parent.frame()))
  #   clusterExport(cl,"computeCov")
  # else{
  #   cl_func <- makePSOCKcluster(2)
  #   registerDoParallel(cl_func)
  #   clusterExport(cl_func,"computeCov")
  # }
  #Obtain the simulation parameters
  parms=getParms()
  parms$refGeneSignal=0.01
  #Create true mean expression and dropSeq mean expression
  #First refNum column is alway the reference gene in the distMap dataset
  patternData=createPattern(geometry,geneNum,refNum,parms)
  #Create insitu mean expression
  #Two datasets do not have to measure the exactly same species, may have some gene variation.
  #species effect is added in the insitu data. It is controled by parms$speciesDiff
  patternData=addSpeciesDiff(patternData,parms)
  #===========================simulation=============================
  #Simulate insitu data
  sim_insitu=sampleInsituData(patternData,refNum,parms)
  #simulate dropSeq data as well as the cell's location
  sim_dropSeq_list=sampleDropData(patternData,cellNum,parms)
  sim_dropSeq=sim_dropSeq_list$dropData
  sim_cellLoc=sim_dropSeq_list$location
  
  sim_dropSeq=addInflationError(sim_dropSeq,parms)
  
  mydata=list(refNum=refNum,geneNum=geneNum,cellNum=cellNum,insitu=sim_insitu,drop=sim_dropSeq)
  metaData=list(cell_loc=sim_cellLoc,patternData=patternData,simData=mydata)
  
  return(metaData)
}



