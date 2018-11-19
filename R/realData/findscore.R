
refGeneInd=1:84
refNum=40
#geneData=pickGene(simulation$simData,refNum)
index=sample(refGeneInd,refNum)
nchange=1
simData=simulation$simData
cell_loc=simulation$cell_loc
stopRule=0.001
record_ind=c()
record_score=c()
smallest_score=0
repeat{
  clusterExport(cl,"index")
  res=foreach(i=1:(length(cl)*5),.combine = c,.packages = "DistMap")%dopar%{
    ind=sample(1:refNum,nchange)
    ind_sub=sample(refGeneInd[-index],nchange)
    index[ind]=ind_sub
    
    geneData=pickGene(simData,refNum,index)
    geneNum=geneData$geneNum
    cellNum=geneData$cellNum
    dm=author_method(chanllege.data,geneData)
    score=author_method_score(dm,cell_loc)
    result=list(score=score,index=index)
    list(result)
  }
  smallest_score_previous=smallest_score
  smallest_score=0
  smallest_ind=NULL
  for(j in 1:length(res)){
    tmp=res[[j]]
    if(tmp$score>smallest_score){
      smallest_ind=j
      smallest_score=tmp$score
    }
    record_score=c(record_score,tmp$score)
    record_ind=rbind(record_ind,tmp$index)
  }
  if(-(smallest_score_previous-smallest_score)<stopRule)
    break
  index=res[[smallest_ind]]$index
  message("one run")
}

plot(record_score)


record_ind[which.max(record_score),]
ind=record_ind[which.max(record_score),]
