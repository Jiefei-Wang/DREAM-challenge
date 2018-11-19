
drop_data=geneData$drop[1:84,]
cutoff
drop_data[drop_data>cutoff]=cutoff
drop_data1=t(drop_data)


geneNum=20
ind=1:geneNum
score_previous=0
largest_score=0
record_score=c()
record_ind=c()
stopRule=0.001
set.seed(1)
for(i in 1:40){
  res=foreach(i=1:(length(cl)*10),.combine = c,.packages = "DistMap")%dopar%{
    nchange=sample(1:10,1)
    index=ind
    candidate=1:84
    candidate=candidate[-index]
    ind_sub=sample(1:length(index),nchange)
    index[ind_sub]=sample(candidate,nchange)
    
    fit_target=drop_data1[,-index]
    predictor=drop_data1[,index]
    r2=c()
    for(i in 1:ncol(fit_target)){
      out=fit_target[,i]
      fit=summary(lm(out~predictor))
      r2=c(r2,fit$r.squared)
    }
    score=mean(r2)
    list(list(score=score,index=index))
  }
  
  score_previous=largest_score
  largest_score=0
  largest_ind=NULL
  for(j in 1:length(res)){
    tmp=res[[j]]
    if(tmp$score>largest_score){
      largest_ind=j
      largest_score=tmp$score
    }
    record_score=c(record_score,tmp$score)
    record_ind=rbind(record_ind,tmp$index)
  }
  #ind=res[[largest_ind]]$index
  message(largest_score)
  maxInd=which.max(record_score)
  ind=sort(record_ind[maxInd,])
}
plot(record_score)
maxInd=which.max(record_score)
record_score[maxInd]
ind=sort(record_ind[maxInd,])
