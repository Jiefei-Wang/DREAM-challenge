
if(FALSE){
dm = new("DistMap",
         raw.data=as.matrix(chanllege.data$dropSeq.raw),
         data=as.matrix(chanllege.data$dropSeq.normalized),
         insitu.matrix=as.matrix(chanllege.data$insitu.binary),
         geometry=as.matrix(geometry))
dm <- binarizeSingleCellData(dm, seq(0.15, 0.5, 0.01))
dm <- mapCells(dm)


gene.name.list=rownames(chanllege.data$dropSeq.raw)
cell.loc=apply(dm@mcc.scores,2,which.max)
pattern=c()


pattern=foreach(i=1:length(gene.name.list),.combine = cbind,.packages="DistMap")%dopar%{
  gene.level.pred = computeVISH(dm, i, threshold=0.75)
  gene.level.pred
}

insitu=chanllege.data$insitu.raw
reference_gene=colnames(chanllege.data$insitu.raw)

reference_gene[reference_gene=="Blimp.1"]="Blimp-1"
reference_gene[reference_gene=="E.spl.m5.HLH"]="E(spl)m5-HLH"
reference_ind=unlist(sapply(reference_gene,function(x,name.list){which(name.list==x)},name.list=gene.name.list))
drop_ref=chanllege.data$dropSeq.raw[reference_ind,]
drop_nonref=chanllege.data$dropSeq.raw[-reference_ind,]
pattern_ref=insitu
pattern_nonref=pattern[,-reference_ind]


pattern_true=cbind(pattern_ref,pattern_nonref)
#insitu
drop=rbind(drop_ref,drop_nonref)
cell_loc=cell.loc


simulation=list()
simulation$cell_loc=cell_loc
simulation$patternData=list()
simulation$patternData$trueTable=pattern_true
simulation$simData=list()
geneData=list()
geneData$refNum=length(reference_ind)
geneData$geneNum=length(gene.name.list)
geneData$cellNum=ncol(drop)
geneData$insitu=insitu
geneData$drop=drop
simulation$simData=geneData

save(simulation,file="authorData.RData")
}

load("R\\commonFunc\\authorData.RData")
