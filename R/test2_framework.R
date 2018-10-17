#install.packages("Rfast")
#install.packages("rpgm")
#BiocManager::install("edgeR")
library(doParallel)
library(ggplot2) 
#Set the number of the clusters and the packages that will be export to the clusters.
clusterNum=4
clusterPkg=c("rpgm","Rfast","edgeR")
source("R\\commonFunc\\readData.R")
source("R\\commonFunc\\functions.R")
source("R\\frameworkFuncs\\framework.R")


#Create data
geometry=chanllege.data$geometry
refNum=50
geneNum=100
cellNum=1000
simulation=generate_fake(geometry,refNum,geneNum,cellNum)
simData=simulation$simData
#Check the functions, will be used in the next step
getNormFuncs()
getDistFuncs()
getPatternFuncs()


#attach function to the data and compute the performance
#The parameter can be the function name or the index of the function obtained from the above three functions
#modelList=attachFunc_list(simData,normalization=c("quantile"),distance=c("mcc"),pattern=c("simple1","author"))
#modelList=c(modelList,attachFunc_list(simData,normalization=c("upqu"),distance=c("cov"),pattern=c("simple1","author")))

#modelList=attachFunc_list(simData,normalization=c("upqu"),distance=c("cov"),pattern=c("author"))

modelList=attachFunc_list(simData,normalization=c("rowMax","Mingmei"),distance=c("cov","mse"),pattern=c("author"))



#Compute the performance
result=computePerformance(modelList,simulation,parallel=T)
result[order(result$prediction_score,decreasing=T),]

#The relationship between scores
result$model=factor(paste(result$normalization,result$distance,sep="+"))
ggplot(result, aes(x=prediction_score, y=pattern_score_test,color=model)) + geom_point()+facet_grid(. ~ pattern)


#The author's method
dm_list=originalMethod(simData,simulation)
dm_list$score

#Let's see one pattern
myModel=modelList[[2]]
myModel=predict_all(myModel)
mean(patternScore(myModel$pattern,simulation$patternData$dropTable,refNum+1,geneNum))
mean(predictionScore(myModel$loc,simulation$cell_loc))

#The pattern of a gene
gene=62
intensityPlot2(simulation$patternData$dropTable[,gene],geometry,title="true pattern")
intensityPlot2(myModel$pattern[,gene],geometry,title="predicted pattern")
intensityPlot2(dm_list$pattern[,gene],geometry,title="predicted pattern, author's method")


#take a look at prediction accuracy
ind=c()
for(i in 1:cellNum){
  ind[i]=rank(myModel$distance[,i])[simulation$cell_loc[i]]
}
mean(ind)

