insitu=chanllege.data$insitu.binary
drop=chanllege.data$dropSeq.normalized

reference_gene_insitu=colnames(insitu)
#reference_gene_drop=nameConvert(reference_gene_insitu,"drop")
#gene_drop=rownames(drop)
reference_gene_drop_ind=which(gene_drop%in% reference_gene_insitu)

drop_ref=drop[reference_gene_drop_ind,]
drop_nonref=drop[-reference_gene_drop_ind,]
drop_new=rbind(drop_ref,drop_nonref)


insitu=as.matrix(insitu)
drop_new=as.matrix(drop_new)
rownames(drop_new)=nameConvert(rownames(drop_new),"insitu")
cell_loc=simulation$cell_loc

simulation=list()
simulation$cell_loc=cell_loc
simulation$simData=list()
simulation$simData$refNum=84
simulation$simData$geneNum=nrow(drop_new)
simulation$simData$cellNum=ncol(drop_new)
simulation$simData$insitu=insitu
simulation$simData$drop=drop_new

save(simulation,file="authorBinaryData.RData")
