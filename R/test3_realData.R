library(doParallel)
library(ggplot2) 
library(tictoc)
#Set the number of the clusters and the packages that will be export to the clusters.
clusterNum=detectCores()-1
clusterPkg=c("np","rpgm","Rfast","edgeR")
source("R\\commonFunc\\createCluster.R")
source("R\\commonFunc\\readData.R")
source("R\\commonFunc\\functions.R")
source("R\\frameworkFuncs\\framework.R")
source("R\\commonFunc\\authorsData.R")


simulation$patternData$dropTable=as.matrix(simulation$patternData$trueTable)
simulation$simData$insitu=as.matrix(simulation$simData$insitu)
simulation$simData$drop=as.matrix(simulation$simData$drop)
geneData=simulation$simData
refNum=geneData$refNum
geneNum=geneData$geneNum
cellNum=geneData$cellNum


getNormFuncs()
getDistFuncs()
getPatternFuncs()


#attach function to the data and compute the performance
#The parameter can be the function name or the index of the function obtained from the above three functions
modelList=attachFunc_list(normalization=c("GIWC","GIAC","LogGIAC","author_norm","author_quantile","quantile","scale","TMM","upqu","double_upqu"),
                          distance=c("cov","cov_rank","mse"),
                          pattern=c("simple1"))
modelList=c(modelList,attachFunc_list(normalization=c("quantile","author_quantile"),distance=c("mcc"),pattern=c("simple1")))



modelList=attachFunc_list(normalization=c("double_upqu"),distance=c("cov"),pattern=c("simple1"))


modelList=attachFunc_list(normalization=c("GIWC","GIAC"),distance=c("cov"),pattern=c("simple1"))


#Compute the performance
tic()
result=computePerformance(modelList,geneData,simulation,parallel=F)
toc()
result[order(result$prediction_score,decreasing=T),1:6]


mydata=modelList[[2]]
mymodel=predict_all(mydata,geneData,pattern=F)
mymodel=predict_pattern(mymodel,patternInd=1,mymodel$refNum)
gene=2
intensityPlot2(simulation$patternData$dropTable[,gene],geometry,title="true pattern")
intensityPlot2(mymodel$pattern[,gene],geometry,title="predicted pattern")






dm = new("DistMap",
         raw.data=as.matrix(chanllege.data$dropSeq.raw),
         data=as.matrix(chanllege.data$dropSeq.normalized),
         insitu.matrix=as.matrix(chanllege.data$insitu.binary),
         geometry=as.matrix(geometry))
dm <- binarizeSingleCellData(dm, seq(0.15, 0.5, 0.01))
dm <- mapCells(dm)


gene.name.list=rownames(chanllege.data$dropSeq.raw)
reference_gene=colnames(chanllege.data$insitu.raw)
reference_gene[reference_gene=="Blimp.1"]="Blimp-1"
reference_gene[reference_gene=="E.spl.m5.HLH"]="E(spl)m5-HLH"
reference_ind=unlist(sapply(reference_gene,function(x,name.list){which(name.list==x)},name.list=gene.name.list))
gene.level.pred = computeVISH(dm, reference_ind[gene], threshold=0.75)

intensityPlot2(gene.level.pred,geometry,title="predicted pattern")




