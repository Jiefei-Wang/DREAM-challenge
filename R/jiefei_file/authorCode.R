
source("R\\commonFunc\\readData.R")

dm = new("DistMap",
         raw.data=as.matrix(chanllege.data$dropSeq.raw),
         data=as.matrix(chanllege.data$dropSeq.normalized),
         insitu.matrix=as.matrix(chanllege.data$insitu.binary),
         geometry=as.matrix(geometry))
dm <- binarizeSingleCellData(dm, seq(0.15, 0.5, 0.01))
dm <- mapCells(dm)


gene.name.list=rownames(chanllege.data$dropSeq.raw)
ref.name.list=colnames(chanllege.data$insitu.raw)
gene.name="odd"

if(!gene.name%in%gene.name.list)
  stop("The gene is not in the drop sequence data")
if(gene.name%in% ref.name.list)
  stop("The gene is in the reference data")



#plot the predicted gene level expression patterns
gene.level.pred = computeVISH(dm, gene.name, threshold=0.75)
intensityPlot(gene.level.pred,geometry,title=paste(gene.name,": predicted"))

res=data.frame(x=geometry$x,z=geometry$z,gene.level.pred)
colnames(res)[3]=gene.name

write.csv(res,file="R\\jiefei_file\\plot.csv",row.names = F)

