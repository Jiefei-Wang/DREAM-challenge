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
mydataList=attachFunc_list(mydata,normalization=c("quantile"),distance=c("mcc"),pattern=c("simple1","author"))
mydataList=c(mydataList,attachFunc_list(mydata,normalization=c("upqu"),distance=c("cov"),pattern=c("simple1","author")))


#Compute the performance
result=computePerformance(mydataList,simulation,parallel=T)
result[order(result$prediction_score,decreasing=T),]

#The relationship between scores
ggplot(result, aes(x=prediction_score, y=pattern_score,color=normalization,shape=distance)) + geom_point()+facet_grid(. ~ pattern)


#The author's method
dm_list=originalMethod(mydata,simulation)
dm_list$score

#Let's see one pattern
mydata=mydataList[[6]]
mydata=predict_all(mydata)
mean(patternScore(mydata$pattern,simulation$patternData$dropTable,refNum+1,geneNum))
mean(predictionScore(mydata$loc,simulation$cell_loc))

#The pattern of a gene
gene=62
intensityPlot2(simulation$patternData$dropTable[,gene],geometry,title="true pattern")
intensityPlot2(mydata$pattern[,gene],geometry,title="predicted pattern")
intensityPlot2(dm_list$pattern[,gene],geometry,title="predicted pattern, author's method")


#take a look at prediction accuracy
ind=c()
for(i in 1:cellNum){
  ind[i]=rank(mydata$distance[,i])[simulation$cell_loc[i]]
}
mean(ind)

