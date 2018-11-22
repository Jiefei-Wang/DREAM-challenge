library(doParallel)
library(tictoc)
library(caret)
#Set the number of the clusters and the packages that will be export to the clusters.
clusterNum=detectCores()-1
clusterPkg=c("Rfast","caret")
folder="R\\submission1\\"
source(paste0(folder,"createCluster.R"))
source(paste0(folder,"crossValidation.R"))
source(paste0(folder,"framework.R"))
source(paste0(folder,"lossFunc.R"))
source(paste0(folder,"modelFuncs.R"))
source(paste0(folder,"tools.R"))
load(paste0(folder,"authorBinaryData.RData"))

for(pkg in clusterPkg)
  library(pkg,character.only=T)
clusterExport(cl,"weighted_cor")



#20 genes
ind=c(75,14,81,19,10,5,4,20,72,23,78,44,80,50,62,79,11,25,48,2)
#40 genes
ind=c(16,71,84,9,77,69,6,68,1,73,15,24,32,57,64,26,55,54,41,66,75,14,81,19,10,5,4,20,72,23,78,44,80,50,62,79,11,25,48,2)
#80 genes
ind=c(52,29,43,40,47,21,65,37,58,56,28,3,30,74,31,63,36,39,15,49,64,69,16,26,1,68,55,77,9,6,75,10,84,14,53,32,66,5,19,71,57,54,73,78,20,72,24,4,81,41,23,62,80,44,11,50,25,79,2,48)

sim1=pickGene(geneData,ind=ind)

modelList=buildModel()

result1=CV(modelList,sim1,foldNum=10)

mydata=modelList[[which.max(result1)]]
mymodel=normalize(mydata,sim1)
mymodel=compute_dist(mymodel)
mymodel=pred_loc(mymodel)








