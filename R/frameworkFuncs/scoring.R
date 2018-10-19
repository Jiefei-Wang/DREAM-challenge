computePerformance<-function(dataList,simulation,parallel){
  performance=data.frame()
  if(parallel){
    #Add the global function into export list
    g=ls(globalenv())
    global_func=c()
    for(i in 1:length(g))
      if(is.function(get(g[i])))
        global_func=c(global_func,g[i])
    export=c(ls(as.environment("Scoring")),global_func)
    #compute the score
    insitu=dataList[[1]]$insitu
    drop=dataList[[1]]$drop
    for(i in 1:length(dataList)){
      dataList[[i]]$insitu=c()
      dataList[[i]]$drop=c()
    }
    performance=foreach(i=1:length(dataList),.combine= rbind,.multicombine=TRUE,.inorder=FALSE,.export=export,.packages=clusterPkg)%dopar%{
      curData=dataList[[i]]
      curData$insitu=insitu
      curData$drop=drop
      curData=predict_all(curData)
      pattern_score_test=patternScore(curData$pattern,simulation$patternData$dropTable,curData$refNum+1,curData$geneNum)
      pattern_score_train=patternScore(curData$pattern,simulation$patternData$dropTable,1,curData$refNum)
      pattern_score_test=mean(pattern_score_test,na.rm = T)
      pattern_score_train=mean(pattern_score_train,na.rm = T)
      
      pred_score=predictionScore(curData$loc,simulation$cell_loc)
      pred_score=pred_score
      data.frame(normalization=curData$funcName[1],
                 distance=curData$funcName[2],
                 pattern=curData$funcName[3],
                 pattern_score_train=pattern_score_train,
                 pattern_score_test=pattern_score_test
                 ,prediction_score=pred_score)
    }
  }else
    for(i in 1:length(dataList)){
      curData=dataList[[i]]
      curData=predict_all(curData)
      pattern_score_test=patternScore(curData$pattern,simulation$patternData$dropTable,curData$refNum+1,curData$geneNum)
      pattern_score_train=patternScore(curData$pattern,simulation$patternData$dropTable,1,curData$refNum)
      pattern_score_test=mean(pattern_score_test,na.rm = T)
      pattern_score_train=mean(pattern_score_train,na.rm = T)
      pred_score=predictionScore(curData$loc,simulation$cell_loc)
      pred_score=pred_score
      performance=rbind(performance,data.frame(normalization=curData$funcName[1],
                                               distance=curData$funcName[2],
                                               pattern=curData$funcName[3],
                                               pattern_score_train=pattern_score_train,
                                               pattern_score_test=pattern_score_test,
                                               prediction_score=pred_score))
    }
  performance$pattern_score_train=round(performance$pattern_score_train,3)
  performance$pattern_score_test=round(performance$pattern_score_test,3)
  performance$prediction_score=round(performance$prediction_score,3)
  return(performance)
}


#predPattern=mydata$pattern
#truePattern=simulation$patternData$dropTable
#gene_start=51
#gene_end=100

patternScore<-function(predPattern,truePattern,gene_start,gene_end){
  turning=0.01
  weight=0.8
  predPattern=matrix(predPattern[,gene_start:gene_end],nrow(predPattern))
  truePattern=matrix(truePattern[,gene_start:gene_end],nrow(predPattern))
  weight_sig=1
  n_loc=nrow(truePattern)
  true_cl=(truePattern!=0)
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


patternScore1<-function(predPattern,truePattern,gene_start,gene_end){
  turning=0.01
  predPattern=matrix(predPattern[,gene_start:gene_end],nrow(predPattern))
  truePattern=matrix(truePattern[,gene_start:gene_end],nrow(predPattern))
  score=rep(NA,ncol(predPattern))
  for(i in 1:ncol(predPattern)){
    score[i]=cor(predPattern[,i],truePattern[,i],method="spearman")
  }
  return(score)
}

patternScore1<-function(predPattern,truePattern,gene_start,gene_end){
  turning=0.01
  predPattern=matrix(predPattern[,gene_start:gene_end],nrow(predPattern))
  truePattern=matrix(truePattern[,gene_start:gene_end],nrow(predPattern))
  weight_sig=1
  n_loc=nrow(truePattern)
  true_cl=(truePattern!=0)
  n_sig=colSums(true_cl)
  score=rep(NA,ncol(predPattern))
  for(i in 1:ncol(predPattern)){
    pred_cl=(predPattern[,i]>quantile(predPattern[,i],1-n_sig[i]/n_loc))
    mytable=table(factor(true_cl[,i],levels = c(F,T)),factor(pred_cl,levels = c(F,T)))
    margin1=rowSums(mytable)/n_loc
    margin2=colSums(mytable)/n_loc
    E_specificity=margin1[1]*margin2[1]
    E_sensityvity=margin1[2]*margin2[2]
    T_specificity=(mytable[1,1]+turning)/(n_loc-n_sig[i]+turning)
    T_sensityvity=(mytable[2,2]+turning)/(n_sig[i]+turning)
    scale1=min(E_specificity,1-E_specificity)
    scale2=min(E_sensityvity,1-E_sensityvity)
    if(scale1<0.01)scale1=0.01
    if(scale2<0.01)scale2=0.01
    df_specificity=(T_specificity-E_specificity)/scale1
    df_sensityvity=(T_sensityvity-E_sensityvity)/scale2
    score[i]=(df_specificity+weight_sig*df_sensityvity)/(1+weight_sig)
  }
  return(score)
}

predLocation<-function(mydata){
  predict_num=10
  #Compute the rank
  rank_distance=apply(mydata$distance,2,rank,ties.method="random")
  ind=which(rank_distance<=predict_num,arr.ind=T)
  sort_order=apply(matrix(rank_distance[ind],nrow=predict_num),2,order)
  loc=matrix(ind[,1],nrow=predict_num)
  mydata$loc=sapply(1:ncol(loc),function(i,loc,ind){
    loc[ind[,i],i]
  },loc=loc,ind=sort_order)
  return(mydata)
}



#trueLoc=simulation$cell_loc
#predLoc=mydata$loc
predictionScore<-function(predLoc,trueLoc){
  matchPos=sweep(predLoc,2,trueLoc,"==")
  ind=which(matchPos,arr.ind = T)
  sum(1/ind[,1])/length(trueLoc)
}