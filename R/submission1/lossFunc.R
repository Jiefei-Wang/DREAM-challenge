
patternScore_SS<-function(predPattern,truePattern,gene_start,gene_end){
  turning=0.01
  weight=1
  predPattern=matrix(predPattern[,gene_start:gene_end],nrow(predPattern))
  truePattern=matrix(truePattern[,gene_start:gene_end],nrow(predPattern))
  truePattern1=apply(truePattern,2,function(x){
    fit=kmeans(x,2)
    gmean=aggregate(x,by=list(fit$cluster),FUN=mean)
    group=fit$cluster-1
    if(gmean[1,2]>gmean[2,2]) 
    {
      group=1-group
    }
    group
  })
  weight_sig=1
  n_loc=nrow(truePattern1)
  true_cl=(truePattern1!=0)
  n_sig=colSums(true_cl)
  score=rep(NA,ncol(predPattern))
  for(i in 1:ncol(predPattern)){
    pred_cl=(predPattern[,i]>quantile(predPattern[,i],1-n_sig[i]/n_loc))
    mytable=table(factor(true_cl[,i],levels = c(F,T)),factor(pred_cl,levels = c(F,T)))
    margin1=rowSums(mytable)/n_loc
    margin2=colSums(mytable)/n_loc
    T_specificity=(mytable[1,1]+turning)/(n_loc-n_sig[i]+turning)
    T_sensityvity=(mytable[2,2]+turning)/(n_sig[i]+turning)
    score[i]=weight*min(c(T_specificity,T_sensityvity))+(1-weight)*max(c(T_specificity,T_sensityvity))
  }
  return(score)
}
