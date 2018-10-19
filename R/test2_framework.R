#install.packages("Rfast")
#install.packages("rpgm")
#BiocManager::install("edgeR")
library(doParallel)
library(ggplot2) 
#Set the number of the clusters and the packages that will be export to the clusters.
clusterNum=10
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
mydata=simulation$simData
#Check the functions, will be used in the next step
getNormFuncs()
getDistFuncs()
getPatternFuncs()


#attach function to the data and compute the performance
#The parameter can be the function name or the index of the function obtained from the above three functions
modelList=attachFunc_list(mydata,normalization=c("IndividualCellIntensity","IndividualGeneIntensity","LogIndividualGeneIntensity","Mingmei","quantile","scale","TMM","upqu"),
                           distance=c("cov","mcc","mse"),
                           pattern=c("simple1","author"))
#modelList=c(modelList,attachFunc_list(mydata,normalization=c("upqu"),distance=c("cov"),pattern=c("simple1","author")))

modelList=attachFunc_list(mydata,normalization=c("IndividualCellIntensity"),distance=c("cov"),pattern=c("author"))
modelList=c(modelList,attachFunc_list(mydata,normalization=c("upqu"),distance=c("cov"),pattern=c("author")))

#Compute the performance
result=computePerformance(modelList,simulation,parallel=T)
result[order(result$prediction_score,decreasing=T),]

#The relationship between scores
result$model=factor(paste(result$normalization,result$distance,sep="+"))
ggplot(result, aes(x=prediction_score, y=pattern_score_test,color=model)) + geom_point()+facet_grid(. ~ pattern)



ggplot(result, aes(x=pattern_score_train, y=pattern_score_test))+ geom_point()


#The author's method
dm_list=originalMethod(mydata,simulation)
dm_list$score

#Let's see one pattern
mydata=modelList[[17]]
mydata=predict_all(mydata)
mean(patternScore(mydata$pattern,simulation$patternData$dropTable,refNum+1,geneNum))
mean(predictionScore(mydata$loc,simulation$cell_loc))

#The pattern of a gene
gene=100
intensityPlot(simulation$patternData$dropTable[,gene],geometry,title="true pattern")
intensityPlot(mydata$pattern[,gene],geometry,title="predicted pattern")
intensityPlot2(dm_list$pattern[,gene],geometry,title="predicted pattern, author's method")


#take a look at prediction accuracy
ind=c()
for(i in 1:cellNum){
  ind[i]=rank(mydata$distance[,i])[simulation$cell_loc[i]]
}
mean(ind)

