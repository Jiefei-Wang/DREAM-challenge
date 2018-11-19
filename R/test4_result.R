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

set.seed(1)
#ind=c(39,52,34,8,26,76,12,51,3,54,5,4,72,6,28,36,7,60,69,77)
ind=c(39,52,34,8,26,76,12,51,3,54,5,4,72,6,28,36,7,60,69,77)
#ind_char=c("cad","Kr","hb","ftz","odd","h","kni","hkb","tll","eve","sna","oc","zen","disco","twi","fkh","brk","rho","run","gt")
#ind=findgeneInd(ind_char,colnames(simulation$simData$insitu))

refNum=length(ind)
filename=paste0("select",refNum)
geneData=pickGene(simulation$simData,refNum,ind)
geneNum=geneData$geneNum
cellNum=geneData$cellNum



modelList=attachFunc_list(normalization=c("TMM"),distance=c("mse"),pattern=c("author"))
#modelList=attachFunc_list(normalization=c("double_upqu"),distance=c("mse"),pattern=c("author"))
mydata=modelList[[1]]
mymodel=predict_all(mydata,geneData,pattern=F)
mymodel=predict_pattern(mymodel,patternInd=1,mymodel$refNum)


pattern_score=patternScore_SS(mymodel$pattern,mymodel$insitu,1,mymodel$refNum)
mean(pattern_score)
barplot(pattern_score)
mean(predictionScore(mymodel$loc,simulation$cell_loc))


dm=author_method(chanllege.data,geneData)
author_method_score(dm,simulation$cell_loc)

#Refit
ind1=which(pattern_score<0.1)
geneData1=geneData
geneData1$insitu=geneData1$insitu[,-ind1]
geneData1$refNum=geneData$refNum-length(ind1)
geneData1$drop=geneData1$drop[-ind1,]


modelList=attachFunc_list(normalization=c("TMM"),distance=c("mse"),pattern=c("author"))
mydata=modelList[[1]]
mymodel=predict_all(mydata,geneData1,pattern=F)
mymodel=predict_pattern(mymodel,patternInd=1,mymodel$refNum)

pattern_score1=patternScore_SS(mymodel$pattern,mymodel$insitu,1,mymodel$refNum)
mean(pattern_score1)
barplot(pattern_score1)
mean(predictionScore(mymodel$loc,simulation$cell_loc))


reference_gene=colnames(geneData$insitu)
reference_gene[reference_gene=="Blimp.1"]="Blimp-1"
reference_gene[reference_gene=="E.spl.m5.HLH"]="E(spl)m5-HLH"

final_res=rbind(matrix(reference_gene,ncol=10),t(mymodel$loc))
write.csv(final_res,file=paste0(filename,".csv"),quote=F)

#=================================40=====================================

set.seed(1)
ind=c(63,76,59,74,17,15,7,51,46,9,16,21,50,12,2,80,49,67,13,29,60,28,64,40,4,34,18,37,68,31,43,3,66,33,38,1,11,6,30,82)
#ind_char=c("cad","Kr","hb","ftz","odd","h","kni","hkb","tll","eve","sna","oc","zen","disco")
#ind=findgeneInd(ind_char,colnames(simulation$simData$insitu))
refNum=length(ind)
filename=paste0("select",refNum)
geneData=pickGene(simulation$simData,refNum,ind)
geneNum=geneData$geneNum
cellNum=geneData$cellNum



modelList=attachFunc_list(normalization=c("double_upqu"),distance=c("mse"),pattern=c("author"))
mydata=modelList[[1]]
mydata$N_parm=c(0.909090909090909,0.909090909090909)
mymodel=predict_all(mydata,geneData,pattern=F)
mymodel=predict_pattern(mymodel,patternInd=1,mymodel$refNum)


pattern_score=patternScore_SS(mymodel$pattern,mymodel$insitu,1,refNum)
mean(pattern_score)
barplot(pattern_score)
mean(predictionScore(mymodel$loc,simulation$cell_loc))


dm=author_method(chanllege.data,geneData)
author_method_score(dm,simulation$cell_loc)

#Refit
ind1=ind[pattern_score>0.1]
refNum=length(ind1)
geneData1=pickGene(simulation$simData,refNum,ind1)


modelList=attachFunc_list(normalization=c("TMM"),distance=c("mse"),pattern=c("author"))
mydata=modelList[[1]]
mymodel=predict_all(mydata,geneData1,pattern=F)
mymodel=predict_pattern(mymodel,patternInd=1,mymodel$refNum)

pattern_score1=patternScore_SS(mymodel$pattern,mymodel$insitu,1,refNum)
mean(pattern_score1)
barplot(pattern_score1)
mean(predictionScore(mymodel$loc,simulation$cell_loc))





reference_gene=colnames(geneData$insitu)
reference_gene[reference_gene=="Blimp.1"]="Blimp-1"
reference_gene[reference_gene=="E.spl.m5.HLH"]="E(spl)m5-HLH"

final_res=rbind(matrix(reference_gene,ncol=10),t(mymodel$loc))
write.csv(final_res,file=paste0(filename,".csv"),quote=F)




#=================================60=====================================

set.seed(1)
ind_char=c("cad","Kr","hb","ftz","odd","h","kni","hkb","tll","eve","sna","oc","zen","disco","twi","fkh","brk","rho","run","gt","Dfd","tsh","prd","srp","ems","Antp","nub","tkv","bowl","knrl","dpn","htl","Doc3","Doc2","numb","peb","danr","zen2","D","Ance","exex","cnc","E.spl.m5.HLH","NetA","noc","Traf4","apt","fj","Ilp4","CG8147","Ama","ImpL2","Btk29A","Cyp310a1","trn","croc","bmm","dan","mfas","Mdr49")
ind=findgeneInd(ind_char,colnames(simulation$simData$insitu))
ind=c(69,26,71,45,43,39,28,76,65,16,36,79,2,18,5,7,22,77,44,27,59,19,29,21,53,8,1,35,70,49,72,34,83,33,24,15,60,54,37,25,20,31,17,58,67,63,38,3,11,6,46,51,4,66,50,30,13,12,68,64)
refNum=length(ind)
filename=paste0("select",refNum)
geneData=pickGene(simulation$simData,refNum,ind)
geneNum=geneData$geneNum
cellNum=geneData$cellNum



modelList=attachFunc_list(normalization=c("double_upqu"),distance=c("mse"),pattern=c("author"))
res=computePerformance_real(modelList,geneData,simulation$cell_loc,parallel=T)
res[which.max(res$pattern_score_train),]
mydata=modelList[[which.max(res$pattern_score_train)]]
mymodel=predict_all(mydata,geneData,pattern=F)
mymodel=predict_pattern(mymodel,patternInd=1,mymodel$refNum)


pattern_score=patternScore_SS(mymodel$pattern,mymodel$insitu,1,refNum)
mean(pattern_score)
barplot(pattern_score)
mean(predictionScore(mymodel$loc,simulation$cell_loc))


dm=author_method(chanllege.data,geneData)
author_method_score(dm,simulation$cell_loc)

#Refit
ind1=ind[pattern_score>0.15]
refNum=length(ind1)
geneData1=pickGene(simulation$simData,refNum,ind1)


modelList=attachFunc_list(normalization=c("TMM"),distance=c("mse"),pattern=c("author"))
mydata=modelList[[1]]
mymodel=predict_all(mydata,geneData1,pattern=F)
mymodel=predict_pattern(mymodel,patternInd=1,mymodel$refNum)

pattern_score1=patternScore_SS(mymodel$pattern,mymodel$insitu,1,refNum)
mean(pattern_score1)
barplot(pattern_score1)
mean(predictionScore(mymodel$loc,simulation$cell_loc))





reference_gene=colnames(geneData$insitu)
reference_gene[reference_gene=="Blimp.1"]="Blimp-1"
reference_gene[reference_gene=="E.spl.m5.HLH"]="E(spl)m5-HLH"

final_res=rbind(matrix(reference_gene,ncol=10),t(mymodel$loc))
write.csv(final_res,file=paste0(filename,".csv"),quote=F)


