library(doParallel)
library(ggplot2) 
library(tictoc)
library(psych)
#Set the number of the clusters and the packages that will be export to the clusters.
clusterNum=detectCores()-1
clusterPkg=c("np","rpgm","Rfast","edgeR")
source("R\\commonFunc\\createCluster.R")
source("R\\commonFunc\\readData.R")
source("R\\commonFunc\\functions.R")
source("R\\frameworkFuncs\\framework.R")
source("R\\realData\\scoring_real.R")
source("R\\realData\\functions.R")

load("R\\commonFunc\\authorBinaryData.RData")

geneData=simulation$simData
ind=c(4,6,9,12,13,16,29,32,38,40,43,53,59,63,65,66,69,73,79,83)
ind=c(2,3,6,7,8,11,12,13,14,15,17,18,22,24,25,26,32,33,34,37,38,42,43,46,50,52,54,58,59,60,63,64,66,67,68,69,70,76,79,81)
sim1=pickGene(geneData,ind=ind)

colnames(sim1$insitu)
rownames(sim1$drop)[1:length(ind)]

getNormFuncs()
getDistFuncs()
getPatternFuncs()

modelList=attachFunc_list(normalization=c("dropClusterOnly_parm"),distance=c("mcc"),pattern=c("simple1"))




#Compute the performance
tic()
result=computePerformance_real(modelList,sim1,simulation$cell_loc,parallel=T)
toc()
result[which.max(result$pattern_score_train),]

dm=author_method(chanllege.data,sim1)
author_method_score(dm,simulation$cell_loc)

author_pattern=author_method_pattern(dm)




mydata=modelList[[which.max(result$pattern_score_train)]]
mymodel=normalize(mydata,geneData)
mymodel=compute_dist(mymodel)
mymodel=pred_loc(mymodel)
mymodel$p_parm[[1]]=50
mymodel$tmp_pattern=NULL
mymodel=predict_pattern(mymodel,patternInd=1,84)

gene=2

intensityPlot3(geneData$insitu[,gene],geometry,title="true pattern")
intensityPlot3(mymodel$pattern[,gene],geometry,title="predicted pattern")
intensityPlot3(author_pattern[,gene],geometry,title="predicted pattern")


cell=317
loc=simulation$cell_loc[cell]
intensityPlot3(max(mymodel$distance[,cell])-mymodel$distance[,cell],geometry,title="predicted pattern")
points(geometry$x[loc],geometry$z[loc],col="blue",pch=20,cex=1)

index=c()
for(i in 1:ncol(mymodel$drop)){
  index[i]=rank(mymodel$distance[,i])[simulation$cell_loc[i]]
}
median(index)

which(rank(-index)<5)
