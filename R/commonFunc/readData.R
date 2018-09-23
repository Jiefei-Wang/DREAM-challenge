
if(FALSE){
dropSeq.raw=read.csv("data\\dge_raw.txt", sep = "\t",header = F)
rownames(dropSeq.raw) = dropSeq.raw$V1
dropSeq.raw$V1 = NULL

dropSeq.normalized=read.csv("data\\dge_normalized.txt", sep = "\t")
gene.name=rownames(dropSeq.normalized)
gene.name[gene.name=="Blimp-1"]="Blimp.1"
gene.name[gene.name=="E(spl)m5-HLH"]="E.spl.m5.HLH"
rownames(dropSeq.normalized)=gene.name

insitu.binary = read.csv("data\\binarized_bdtnp.csv",check.names=T)
insitu.raw=read.csv("data\\bdtnp.txt", sep = "\t")
geometry=read.csv("data\\geometry.txt",sep=" ")
colnames(geometry)=c("x","y","z")
}

load("R\\commonFunc\\GeneData.RData")

