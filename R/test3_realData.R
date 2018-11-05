library(doParallel)
library(ggplot2) 
library(tictoc)
#Set the number of the clusters and the packages that will be export to the clusters.
clusterNum=detectCores()-1
clusterPkg=c("np","rpgm","Rfast","edgeR")
#source("R\\commonFunc\\createCluster.R")
source("R\\commonFunc\\readData.R")
source("R\\commonFunc\\functions.R")
source("R\\frameworkFuncs\\framework.R")
source("R\\commonFunc\\authorsData.R")
source("R\\realData\\scoring_real.R")
source("R\\realData\\functions.R")


refNum=20
set.seed(1)
#geneData=pickGene(simulation$simData,refNum)
ind=c(39,52,34,8,26,76,12,51,3,54,5,4,72,6,28,36,7,60,69,77)
geneData=pickGene(simulation$simData,refNum,ind)
geneNum=geneData$geneNum
cellNum=geneData$cellNum



getNormFuncs()
getDistFuncs()
getPatternFuncs()


#attach function to the data and compute the performance
#The parameter can be the function name or the index of the function obtained from the above three functions
modelList=attachFunc_list(normalization=c("GIWC","GIAC","LogGIAC","author_norm","author_quantile","quantile","scale","TMM","upqu","double_upqu"),
                          distance=c("cov","cov_rank","mse"),
                          pattern=c("author"))
modelList=c(modelList,attachFunc_list(normalization=c("quantile","author_quantile"),distance=c("mcc"),pattern=c("author")))



modelList=attachFunc_list(normalization=c("double_upqu"),distance=c("cov"),pattern=c("simple1"))


modelList=attachFunc_list(normalization=c("double_upqu"),distance=c("mse"),pattern=c("author"))


#Compute the performance
tic()
result=computePerformance_real(modelList,geneData,simulation$cell_loc,parallel=T)
toc()
result[order(result$prediction_score,decreasing=T),1:6]


mydata=modelList[[108]]
mymodel=predict_all(mydata,geneData,pattern=F)
mymodel=predict_pattern(mymodel,patternInd=1,mymodel$refNum)
gene=2
intensityPlot3(geneData$insitu[,gene],geometry,title="true pattern")
intensityPlot3(mymodel$pattern[,gene],geometry,title="predicted pattern")
intensityPlot3(author_pattern[,gene],geometry,title="predicted pattern")
intensityPlot3(pattern_author_res[,geneData$geneInd[gene]],geometry,title="predicted pattern")

dm=author_method(chanllege.data,geneData)
author_method_score(dm,simulation$cell_loc)
author_pattern=author_method_pattern(dm)







pattern_score=patternScore_SS(pattern_author_res,as.matrix(chanllege.data$insitu.raw),1,84)

n=20
ind=which(rank(pattern_score,ties.method="random")<=20)


result1=result[result$normalization=="double_upqu"&result$distance=="cov_rank",]


result1[which.max(result1$pattern_score_train),]
