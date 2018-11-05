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



#=================================20=====================================
refNum=20
set.seed(1)
#geneData=pickGene(simulation$simData,refNum)
ind=c(39,52,34,8,26,76,12,51,3,54,5,4,72,6,28,36,7,60,69,77)
geneData=pickGene(simulation$simData,refNum,ind)
geneNum=geneData$geneNum
cellNum=geneData$cellNum



modelList=attachFunc_list(normalization=c("TMM"),distance=c("mse"),pattern=c("author"))



mydata=modelList[[1]]
mymodel=predict_all(mydata,geneData,pattern=F)
mymodel=predict_pattern(mymodel,patternInd=1,mymodel$refNum)


mean(patternScore_SS(mymodel$pattern,mymodel$insitu,1,refNum))
mean(predictionScore(mymodel$loc,simulation$cell_loc))


dm=author_method(chanllege.data,geneData)
author_method_score(dm,simulation$cell_loc)
author_pattern=author_method_pattern(dm)

colnames(mymodel$insitu)


#=================================40=====================================




