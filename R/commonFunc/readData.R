
if(FALSE){
dropSeq.raw=read.csv("data\\dge_raw.txt", sep = "\t",header = F)
rownames(dropSeq.raw) = dropSeq.raw$V1
dropSeq.raw$V1 = NULL

dropSeq.normalized=read.csv("data\\dge_normalized.txt", sep = "\t")
data.name=rownames(dropSeq.normalized)
data.name[data.name=="Blimp-1"]="Blimp.1"
data.name[data.name=="E(spl)m5-HLH"]="E.spl.m5.HLH"
rownames(dropSeq.normalized)=data.name

insitu.binary = read.csv("data\\binarized_bdtnp.csv",check.names=T)
insitu.raw=read.csv("data\\bdtnp.txt", sep = "\t")
geometry=read.csv("data\\geometry.txt",sep=" ")
colnames(geometry)=c("x","y","z")
}

load("R\\commonFunc\\GeneData.RData")
