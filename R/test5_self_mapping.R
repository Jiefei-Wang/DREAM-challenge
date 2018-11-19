library(doParallel)
library(ggplot2) 
library(tictoc) 
#Set the number of the clusters and the packages that will be export to the clusters.
clusterNum=detectCores()-1
#BiocManager::install("ExpressionNormalizationWorkflow")
library("ExpressionNormalizationWorkflow")

source("R\\commonFunc\\readData.R")
source("R\\commonFunc\\authorsData.R")

geneData=simulation$simData

data_normal=apply(geneData$drop,1,function(x,q){
  x[x>quantile(x,q)]=quantile(x,q)
  x
},q=0.9)
data_normal1=as.data.frame(t(data_normal))
tmp=data.frame(matrix(0,ncol(data_normal1),1))
rownames(tmp)=colnames(data_normal1)
#====================================================
inpData=expSetobj(data_normal1,tmp)

driver_gene=rownames(data_normal1)[1:84]

pct_thrsh <- 0.75 
## Perform the PVCA analysis
res=pvcAnaly(inpData, pct_thrsh, driver_gene)










data_normal=geneData$drop
data_normal[data_normal>5]=5
#data_normal=geneData$drop

empInd=c()
for(i in 1:nrow(data_normal)){
  if(var(data_normal[i,])<=0.1)
    empInd=c(empInd,i)
}

data_normal1=data_normal[-empInd,]


data_reduct=prcomp(t(data_normal1),
                   center = TRUE,
                   scale. = TRUE) 
pca_data=data_reduct$x

final_data=t(pca_data)[1:10,]

#heatmap(data_cov,symm=T)

range(data_cov)

var(as.vector(data_cov))




data_normal=geneData$drop[1:84,]
data_normal[data_normal>5]=1
hist(as.vector(data_normal))
data_cov=cor(data_normal,method="pearson")

#====================================================
cell_loc=simulation$cell_loc[1,]
col=rainbow(max(cell_loc))
cell_col=col[cell_loc]
cell_coor=as.matrix(geometry[cell_loc,c(1,3)])
record=c()
for(i in 1:geneData$cellNum){
  tmp=data_cov[i,]
  tmp_coor=sweep(cell_coor,2,cell_coor[i,],"-")
  tmp_dist=sqrt(tmp_coor[,1]^2+tmp_coor[,2]^2)
  record=rbind(record,cbind(tmp_dist,tmp))
}

plot(tmp_dist,tmp)




ind=sample(1:nrow(record),5000)
plot(record[ind,1],record[ind,2],type="p",pch=20)


#=============================================





plot(cell_coor[,1],cell_coor[,2],col=cell_col)

newdata=list()
newdata$loc=cell_loc[order(cell_coor[,1])]
newdata$drop=geneData$drop[,order(cell_coor[,1])]
newdata$col=rainbow(length(newdata$loc))
newdata$cell_cor=geometry[newdata$loc,c(1,3)]
plot(newdata$cell_cor$x,newdata$cell_cor$z,col=newdata$col)



data_normal=apply(newdata$drop,2,function(x,q){
  x[x>quantile(x,q)]=quantile(x,q)
  x
},q=0.9)

data_cov=cor(data_normal,method="pearson")

n=geneData$cellNum
node=matrix(0,n*(n-1)/2,3)
ind=1
for(i in 1:(n-1)){
  for(j in (i+1):n){
    node[ind,]=c(i,j,data_cov[i,j])
    ind=ind+1
  }
}

library(igraph)
g <- graph.data.frame(node,directed = F)

# Edge weights, will be recycled
E(g)$weight <- 10*node[,3]^2
E(g)$color="white"
V(g)$color=newdata$col
coords <- layout.fruchterman.reingold(g, weights=E(g)$weight)

# Eliminate the margin
plot(g, layout=coords, vertex.size=5,vertex.label=NA)



