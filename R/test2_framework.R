#install.packages("Rfast")
#install.packages("rpgm")
#BiocManager::install("edgeR")
library(doParallel)

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
mydataList=attachFunc_list(mydata,normalization=c("columnSum","rowMax","rowMaxlog","TMM","scale"),distance=c("mse","cov"),pattern=c("simple"))

#Compute the performance
result=computePerformance(mydataList,simulation,parallel=T)
result[order(result$prediction_score,decreasing=T),]
#The relationship between scores
plot(result$prediction_score,result$pattern_score,xlab="Prediction score",ylab="Pattern score")

#The author's method
dm_list=originalMethod(mydata,simulation)
dm_list$score

#Let's see one pattern
mydata=mydataList[[9]]
mydata=predict_all(mydata)
mean(patternScore(mydata$pattern,simulation$patternData$dropTable,refNum+1,geneNum))
mean(predictionScore(mydata$loc,simulation$cell_loc))


gene=61
intensityPlot2(simulation$patternData$dropTable[,gene],geometry,title="true pattern")
intensityPlot2(mydata$pattern[,gene],geometry,title="predicted pattern")
intensityPlot2(dm_list$pattern[,gene],geometry,title="predicted pattern, author's method")





