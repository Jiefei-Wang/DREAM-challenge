#80 genes
ind=c(52,29,43,40,47,21,65,37,58,56,28,3,30,74,31,63,36,39,15,49,64,69,16,26,1,68,55,77,9,6,75,10,84,14,53,32,66,5,19,71,57,54,73,78,20,72,24,4,81,41,23,62,80,44,11,50,25,79,2,48)
folder="R\\"
source(paste0(folder,"finalSub.R"))

reference_gene=nameConvert(colnames(sim1$insitu),"insitu")
filename=paste0(length(reference_gene),"gene")
final_res=rbind(matrix(reference_gene,ncol=10),t(mymodel$loc))
write.csv(final_res,file=paste0(filename,".csv"),quote=F)