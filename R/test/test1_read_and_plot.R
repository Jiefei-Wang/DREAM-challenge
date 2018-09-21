
# install.packages("devtools")
# library(devtools)
# install_github("rajewsky-lab/DistMap")
library(DistMap)

source("R\\commonFunc\\readData.R")
source("R\\commonFunc\\functions.R")


#exclude one reference gene from 84 genes in insitu database
gene.name.list=colnames(insitu.binary)
gene.name=gene.name.list[3]
insitu.binary.test=insitu.binary[,colnames(insitu.binary)!=gene.name]

#fit the model
dm = new("DistMap",
         raw.data=as.matrix(dropSeq.raw),
         data=as.matrix(dropSeq.normalized),
         insitu.matrix=as.matrix(insitu.binary.test),
         geometry=as.matrix(geometry))
dm <- binarizeSingleCellData(dm, seq(0.15, 0.5, 0.01))
dm <- mapCells(dm)

par(mfrow=c(1,2))
#plot the predicted gene level expression patterns
gene.level.pred = computeVISH(dm, gene.name, threshold=0.75)
intensityPlot(gene.level.pred,geometry,title=paste(gene.name,": predicted"))

#plot the true gene level expression patterns
gene.level.true=insitu.raw[,colnames(insitu.raw)==gene.name]
gene.level.true[gene.level.true<quantile(gene.level.true,0.75)]=0
intensityPlot(gene.level.true,geometry,title=paste(gene.name,": true"))

par(mfrow=c(1,1))






