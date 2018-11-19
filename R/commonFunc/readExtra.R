library("openxlsx")
extra.data=read.xlsx("R\\commonFunc\\gene result.xlsx")
extra.name=colnames(extra.data)
extra.data=apply(extra.data,2,function(x)x/max(x))
colnames(extra.data)=extra.name

colnames(pattern_author_res)=tolower(colnames(pattern_author_res))

score=patternScore_cor(extra.data,pattern_author_res[,extra.name],1,length(extra.name))
score
#select.name=extra.name[score>0.3]
select.name=c("hsp83")
if(sum(select.name%in% colnames(geneData$insitu))==0){
  
  gene.name=colnames(geneData$insitu)
  gene_name=rownames(chanllege.data$dropSeq.raw)
  gene_name[gene_name=="Blimp-1"]="Blimp.1"
  gene_name[gene_name=="E(spl)m5-HLH"]="E.spl.m5.HLH"
  reference_ind1=unlist(sapply(gene.name,function(x,name.list){which(name.list==x)},name.list=gene_name))
  reference_ind2=unlist(sapply(select.name,function(x,name.list){which(tolower(name.list)==tolower(x))},name.list=gene_name))
  reference_ind=c(reference_ind1,reference_ind2)
  if(length(gene.name)!=length(reference_ind1)){
    stop("subset unsuccessful")
  }
  
  geneData$insitu=cbind(geneData$insitu,extra.data[,select.name])
  drop_ref=chanllege.data$dropSeq.raw[reference_ind,]
  drop_nonref=chanllege.data$dropSeq.raw[-reference_ind,]
  geneData$drop=as.matrix(rbind(drop_ref,drop_nonref))
  geneData$refNum=length(reference_ind)
}else{
  message("The extra gene is already in the insitu data")
}

# 
if(FALSE){
gene=7
par(mfrow=c(1,2))
intensityPlot3(extra.data[,gene],geometry,title="true pattern")
intensityPlot3(pattern_author_res[,extra.name][,gene],geometry,title="true pattern")
par(mfrow=c(1,1))


extra.name %in% colnames(simulation$simData$insitu)
}
