#install.packages("Rfast")
#install.packages("rpgm")
#BiocManager::install("edgeR")
library(doParallel)
library(ggplot2) 
library(tictoc)
#Set the number of the clusters and the packages that will be export to the clusters.
clusterNum=10
clusterPkg=c("np","rpgm","Rfast","edgeR")
source("R\\commonFunc\\createCluster.R")
source("R\\commonFunc\\readData.R")
source("R\\commonFunc\\functions.R")
source("R\\frameworkFuncs\\framework.R")


#Create data
refNum=20
geneNum=100
cellNum=1000
simulation=generate_fake(geometry,refNum,geneNum,cellNum)
geneData=simulation$simData
#Check the functions, will be used in the next step
getNormFuncs()
getDistFuncs()
getPatternFuncs()


#attach function to the data and compute the performance
#The parameter can be the function name or the index of the function obtained from the above three functions
modelList=attachFunc_list(normalization=c("GIWC","GIAC","LogGIAC","author_norm","author_quantile","quantile","scale","TMM","upqu","double_upqu"),
                           distance=c("cov","cov_rank","mse"),
                           pattern=c("simple1"))
modelList=c(modelList,attachFunc_list(normalization=c("quantile","author_quantile"),distance=c("mcc"),pattern=c("simple1")))



modelList=attachFunc_list(normalization=c("GIWC"),distance=c("cov"),pattern=c("author"))
modelList=c(modelList,attachFunc_list(normalization=c("upqu"),distance=c("cov"),pattern=c("author")))
modelList=attachFunc_list(normalization=c("GIWC"),distance=c("cov"),pattern=c("simple1"))

modelList=attachFunc_list(normalization=c("double_upqu"),distance=c("cov"),pattern=c("simple1"))



#Compute the performance
tic()
result=computePerformance(modelList,geneData,simulation,parallel=T)
toc()
result[order(result$prediction_score,decreasing=T),1:6]

#The relationship between scores
result$model=factor(paste(result$normalization,result$distance,sep="+"))
ggplot(result, aes(x=prediction_score, y=pattern_score_test,color=model)) + geom_point()+facet_grid(. ~ pattern)





ggplot(result, aes(x=pattern_score_train, y=pattern_score_test))+ geom_point()


#The author's method
dm_list=originalMethod(mydata,simulation)
dm_list$score

#Let's see one pattern
modelList=attachFunc_list(normalization=c("GIAC"),distance=c("mse"),pattern=c("simple1"))
mydata=modelList[[1]]
mymodel=predict_all(mydata,geneData)
mean(patternScore(mymodel$pattern,simulation$patternData$dropTable,refNum+1,geneNum))
mean(predictionScore(mymodel$loc,simulation$cell_loc))

ind=which(rank(mymodel$distance[,gene])<=10)
mymodel$distance[ind,gene]


#The pattern of a gene
gene=68
intensityPlot2(simulation$patternData$dropTable[,gene],geometry,title="true pattern")
intensityPlot2(mymodel$pattern[,gene],geometry,title="predicted pattern")
#intensityPlot2(dm_list$pattern[,gene],geometry,title="predicted pattern, author's method")


#take a look at prediction accuracy
ind=c()
for(i in 1:cellNum){
  ind[i]=rank(mymodel$distance[,i])[simulation$cell_loc[i]]
}
median(ind)

k=20
which(k+5>ind&ind>k)


cell=69
loc=simulation$cell_loc[cell]
intensityPlot2(max(mymodel$distance[,cell])-mymodel$distance[,cell],geometry,title="Dist")
points(geometry$x[loc],geometry$z[loc],col="blue",pch=20,cex=1)
