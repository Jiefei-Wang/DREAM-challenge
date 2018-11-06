library("openxlsx")
extra.data=read.xlsx("R\\commonFunc\\gene result.xlsx")
extra.name=colnames(extra.data)
extra.data=apply(extra.data,2,function(x)x/max(x))
colnames(extra.data)=extra.name
geneData$insitu=cbind(geneData$insitu,extra.data)

colnames(pattern_author_res)=tolower(colnames(pattern_author_res))

score=patternScore_cor(extra.data,pattern_author_res[,extra.name],1,length(extra.name))

select.name=extra.name[score>0.3]

if(sum(select.name%in% colnames(geneData$insitu))==0){
  geneData$insitu=cbind(geneData$insitu,extra.data[,select.name])
  gene.name=colnames(geneData$insitu)
  reference_ind=unlist(sapply(gene.name,function(x,name.list){which(name.list==x)},name.list=rownames(chanllege.data$dropSeq.raw)))
  drop_ref=chanllege.data$dropSeq.raw[reference_ind,]
  drop_nonref=chanllege.data$dropSeq.raw[-reference_ind,]
  geneData$drop=rbind(drop_ref,drop_nonref)
}else{
  message("The extra gene is already in the insitu data")
}

# 
# gene=8
# par(mfrow=c(1,2))
# intensityPlot3(extra.data[,gene],geometry,title="true pattern")
# intensityPlot3(pattern_author_res[,extra.name][,gene],geometry,title="true pattern")
# par(mfrow=c(1,1))




