computePerformance<-function(dataList,simulation,parallel){
  performance=data.frame()
  if(parallel){
    g=ls(globalenv())
    global_func=c()
    for(i in 1:length(g))
      if(is.function(get(g[i])))
        global_func=c(global_func,g[i])
    export=c(ls(as.environment("Scoring")),global_func)
    performance=foreach(i=1:length(dataList),.combine= rbind,.multicombine=TRUE,.inorder=FALSE,.export=export)%dopar%{
      curData=dataList[[i]]
      curData=predict_all(curData)
      pattern_score=patternScore(curData$pattern,simulation$patternData$dropTable,curData$refNum+1,curData$geneNum)
      pattern_score=mean(pattern_score)
      pred_score=predictionScore(curData$loc,simulation$cell_loc)
      pred_score=pred_score
      data.frame(normalization=curData$funcName[1],
                                               distance=curData$funcName[2],
                                               pattern=curData$funcName[3],
                                               pattern_score=pattern_score,prediction_score=pred_score)
    }
  }else
  for(i in 1:length(dataList)){
    curData=dataList[[i]]
    curData=predict_all(curData)
    pattern_score=patternScore(curData$pattern,simulation$patternData$dropTable,curData$refNum+1,curData$geneNum)
    pattern_score=mean(pattern_score)
    pred_score=predictionScore(curData$loc,simulation$cell_loc)
    pred_score=pred_score
    performance=rbind(performance,data.frame(normalization=curData$funcName[1],
                                             distance=curData$funcName[2],
                                             pattern=curData$funcName[3],
                                             pattern_score=pattern_score,prediction_score=pred_score))
  }
  return(performance)
}

#predPattern=mydata$pattern
#truePattern=simulation$patternData$dropTable
#gene_start=92
#gene_end=92
patternScore<-function(predPattern,truePattern,gene_start,gene_end){
  predPattern=matrix(predPattern[,gene_start:gene_end],nrow(predPattern))
  truePattern=matrix(truePattern[,gene_start:gene_end],nrow(predPattern))
  weight_sig=2
  n_loc=nrow(truePattern)
  true_cl=(truePattern!=0)
  n_sig=colSums(true_cl)
  score=rep(NA,ncol(predPattern))
  for(i in 1:ncol(predPattern)){
    pred_cl=(predPattern[,i]>quantile(predPattern[,i],1-n_sig[i]/n_loc))
    mytable=table(true_cl[,i],pred_cl)
    if(nrow(mytable)==1&ncol(mytable)==1){
      score[i]=1
      next
    }else{
      if(nrow(mytable)==1)
        mytable=t(mytable)
      if(ncol(mytable)==1){
        if(colnames(mytable)=="FALSE")
        {score[i]=mytable[1,1]/(n_loc-n_sig[i])/(1+weight_sig)
        next}
        else
        {score[i]=weight_sig*mytable[2,1]/n_sig[i]/(1+weight_sig)
        next}
      }
    }
    
    score[i]=(mytable[1,1]/(n_loc-n_sig[i])+weight_sig*mytable[2,2]/n_sig[i])/(1+weight_sig)
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