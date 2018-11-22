library(doParallel)
library(tictoc)
library(caret)
#Set the number of the clusters and the packages that will be export to the clusters.
clusterNum=detectCores()-1
clusterPkg=c("np","rpgm","Rfast","caret")
source("R\\commonFunc\\createCluster.R")
source("R\\commonFunc\\readData.R")
source("R\\commonFunc\\functions.R")
source("R\\frameworkFuncs\\framework.R")
source("R\\realData\\scoring_real.R")
source("R\\realData\\functions.R")
source("R\\realData\\crossValidation.R")

load("R\\commonFunc\\authorBinaryData.RData")


geneName=c("Ilp4","zen2","ImpE2","Proc","Poxm","betaTub56D","E(spl)m8-HLH","beta-Spec","CG8147","Chro","Pburs","l(1)sc","CG14427","Su(H)","Hr4","E(spl)mbeta-HLH","Egfr","danr","RluA-1","Taf1","par-1","neur","pnr","CrebB","srp","disco","hth","Thor","Mdr49","CG1136","pros","pdm2","E(spl)m5-HLH","PsGEF","shg","Doc2","Act5C","Traf4","Stat92E","slp1","numb","Ack-like","Ance","ImpL2","zfh1","CenG1A","Cyp310a1","caup","Chi","bowl","Antp","htl","slp2","rho","lack","velo","croc","NetA","salm","side","Btk29A","Mes2","abd-A","CtBP","Abd-B","RhoGEF2","exex","sala","mfas","Doc3","pan")
geneName=nameConvert(geneName,target="drop")
ind=which(geneName%in% colnames(geneData$insitu))
ind=ind[1:20]


geneData=simulation$simData
ind=c(75,14,81,19,10,5,4,20,72,23,78,44,80,50,62,79,11,25,48,2)
ind_1=c(57,54,73,78,20,72,24,4,81,41,23,62,80,44,11,50,25,79,2,48)
ind_jiefei=c(4,6,9,12,13,16,29,32,38,40,43,53,59,63,65,66,69,73,79,83)
ind=c(16,71,84,9,77,69,6,68,1,73,15,24,32,57,64,26,55,54,41,66,75,14,81,19,10,5,4,20,72,23,78,44,80,50,62,79,11,25,48,2)
ind_1=c(64,69,16,26,1,68,55,77,9,6,75,10,84,14,53,32,66,5,19,71,57,54,73,78,20,72,24,4,81,41,23,62,80,44,11,50,25,79,2,48)
ind_jiefei=c(4,6,7,12,13,15,17,18,22,24,26,30,31,32,33,34,36,37,38,40,41,42,43,45,52,55,58,59,60,61,62,63,64,65,66,67,68,69,72,81)
ind=c(52,29,43,40,47,21,65,37,58,56,28,3,30,74,31,63,36,39,15,49,64,69,16,26,1,68,55,77,9,6,75,10,84,14,53,32,66,5,19,71,57,54,73,78,20,72,24,4,81,41,23,62,80,44,11,50,25,79,2,48)
ind_jiefei=c(1,3,5,7,8,10,11,12,13,14,15,16,17,18,19,20,21,22,24,26,27,28,32,33,34,35,36,37,38,39,41,42,44,45,49,50,52,54,55,58,59,60,61,62,63,64,66,67,68,69,70,71,72,75,76,77,78,80,81,83)

sim1=pickGene(geneData,ind=ind)

colnames(sim1$insitu)
rownames(sim1$drop)[1:length(ind)]

getNormFuncs()
getDistFuncs()
getPatternFuncs()

modelList=attachFunc_list(normalization=c("dropClusterOnly_parm"),distance=c("Weighted_mcc","mcc1"),pattern=c("simple1"))

modelList=attachFunc_list(normalization=c("dropClusterOnly_parm"),distance=c("Weighted_mcc"),pattern=c("author"))



#Compute the performance
tic()
result=computePerformance_real(modelList,sim1,simulation$cell_loc,parallel=T)
toc()
result[which.max(result$pattern_score_train),]
result[which.max(result$prediction_score),]

result[which.max(result$pattern_score_train[result$pattern=="simple1(10)"]),]
result[which.max(result$pattern_score_train[result$pattern=="author_mcc(0.6)"]),]

result[result$pattern=="simple1(10)",]
result[result$pattern=="author_mcc(0.6)",]


result1=computePerformance_crossValidation(modelList,sim1,simulation$cell_loc,foldNum=10)
plot(result1$prediction,result1$pattern)
which.max(result1$pattern)
modelList[[which.max(result1$pattern)]]$N_parm


drop=t(geneData$drop[1:84,])
quantileList=seq(0.1,0.5,0.01)
var_record=c()
for(i in 1:length(quantileList)){
  drop_tmp=apply(drop,2,function(x,q){x>quantile(x[x!=0],q)},q=quantileList[i])
  drop_var=diag(var(drop_tmp,drop_tmp))
  var_record=c(var_record,median(drop_var))
}
which.max(var_record)
quantileList[which.max(var_record)]
plot(var_record)


dm=author_method(chanllege.data,geneData)

dm=author_method(chanllege.data,sim1)
author_method_score(dm,simulation$cell_loc)

author_pattern=author_method_pattern(dm)




mydata=modelList[[23]]
mymodel=normalize(mydata,sim1)
mymodel=compute_dist(mymodel)
mymodel=pred_loc(mymodel)
mymodel=predict_pattern(mymodel,patternInd=1,84)



index=sort(ind)

mymodel$insitu=cbind(geneData$insitu[,index],geneData$insitu[,-index])
score=patternScore_SS(mymodel$pattern,mymodel$insitu,1,84)
plot(score)
mean(score[21:84])


mydata=modelList[[which.max(result$prediction_score)]]
mymodel=normalize(mydata,sim1)
mymodel=compute_dist(mymodel)
mymodel=pred_loc(mymodel)
mymodel=predict_pattern(mymodel,patternInd=1,84)



index=sort(ind)

mymodel$insitu=cbind(geneData$insitu[,index],geneData$insitu[,-index])
score=patternScore_SS(mymodel$pattern,mymodel$insitu,1,84)
mean(score[21:84])








geneNameList=rownames(mymodel$drop)[1:84]


gene=8
geneName=geneNameList[gene]
par(mfrow=c(2,2))
intensityPlot3(geneData$insitu[,geneName],geometry,title="true pattern")
intensityPlot3(mymodel$pattern[,gene],geometry,title="predicted pattern")
intensityPlot3(author_pattern[,gene]>quantile(author_pattern[,gene],1-mean(geneData$insitu[,geneName])),geometry,title="predicted pattern")
intensityPlot3(mymodel$pattern[,gene]>quantile(mymodel$pattern[,gene],1-mean(geneData$insitu[,geneName])),geometry,title="quantile")
#intensityPlot3(mymodel1$pattern[,gene],geometry,title="predicted pattern")
par(mfrow=c(1,1))



mydata1=modelList[[which.max(result$prediction_score)]]
mymodel1=predict_all(mydata,sim1,pattern = F)
mymodel1=predict_pattern(mymodel1,patternInd=1,20)



par(mfrow=c(2,2))
intensityPlot3(geneData$insitu[,geneName],geometry,title="true pattern")
intensityPlot3(mymodel$pattern[,gene],geometry,title="predicted pattern")
intensityPlot3(author_pattern[,gene],geometry,title="predicted pattern")
par(mfrow=c(1,1))











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
