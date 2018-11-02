library(doParallel)
library(ggplot2) 
library(tictoc)
#Set the number of the clusters and the packages that will be export to the clusters.
clusterNum=10
clusterPkg=c("np","rpgm","Rfast","edgeR")
#source("R\\commonFunc\\createCluster.R")
source("R\\commonFunc\\readData.R")
source("R\\commonFunc\\functions.R")
source("R\\frameworkFuncs\\framework.R")
source("R\\commonFunc\\authorsData.R")


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


modelList=attachFunc_list(normalization=c("double_upqu"),distance=c("cov"),pattern=c("simple1"))


#Compute the performance
tic()
result=computePerformance(modelList,geneData,simulation,parallel=T)
toc()
result[order(result$prediction_score,decreasing=T),1:6]





