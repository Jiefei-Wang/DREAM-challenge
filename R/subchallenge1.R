#20 genes
ind=c(75,14,81,19,10,5,4,20,72,23,78,44,80,50,62,79,11,25,48,2)
folder="R\\"
source(paste0(folder,"finalSub.R"))

reference_gene=nameConvert(colnames(sim1$insitu),"insitu")
filename=paste0(length(reference_gene),"gene")
final_res=rbind(matrix(reference_gene,ncol=10),t(mymodel$loc))
write.csv(final_res,file=paste0(filename,".csv"),quote=F)
