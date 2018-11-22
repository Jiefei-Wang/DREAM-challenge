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



sim1=pickGene(geneData,ind=ind)

modelList=buildModel()

result1=CV(modelList,sim1,foldNum=10)

mydata=modelList[[which.max(result1)]]
mymodel=normalize(mydata,sim1)
mymodel=compute_dist(mymodel)
mymodel=pred_loc(mymodel)








