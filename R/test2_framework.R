

clusterNum=4
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


mydataList=attachFunc_list(mydata,normalization=c("rowMax","simple"),distance=c("mse","simple"),pattern=c("simple"))

#attach function to the data and compute the performance
#The number is the index of the function obtained from the above three functions
mydataList=attachFunc_list(mydata,normalization=1,distance=2,pattern=1)
mydataList=c(mydataList,attachFunc_list(mydata,normalization=3,distance=2,pattern=1))
mydataList=c(mydataList,attachFunc_list(mydata,normalization=1,distance=1,pattern=1))
mydataList=c(mydataList,attachFunc_list(mydata,normalization=3,distance=1,pattern=1))
#Compute the performance
result=computePerformance(mydataList,simulation,parallel=T)
result


#The author's method
dm_list=originalMethod(mydata,simulation)
dm_list$score

#Let's see one pattern
mydata=mydataList[[2]]
mydata=predict_all(mydata)
mean(patternScore(mydata$pattern,simulation$patternData$dropTable,1,geneNum))

gene=51
intensityPlot2(simulation$patternData$dropTable[,gene],geometry,title="true pattern")
intensityPlot2(mydata$pattern[,gene],geometry,title="predicted pattern")
intensityPlot2(dm_list$pattern[,gene],geometry,title="predicted pattern, author's method")





