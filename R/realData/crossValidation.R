

computePerformance_crossValidation<-function(modelList,geneData,cell_loc,foldNum){
  performance=c()
    #Add the global function into export list
    g=ls(globalenv())
    global_func=c()
    for(i in 1:length(g))
      if(is.function(get(g[i])))
        global_func=c(global_func,g[i])
    export=c(global_func)
    geneData$drop=geneData$drop[1:84,]
    performance=foreach(i=1:length(modelList),.combine= rbind,.multicombine=TRUE,.inorder=FALSE,.export=export,.packages=clusterPkg)%dopar%{
      geneData1=geneData
      performance=c()
      insitu=geneData$insitu
      drop_ref=geneData$drop[1:geneData$refNum,]
      drop_nonref=geneData$drop[-(1:geneData$refNum),]
      
      refNum=geneData$refNum
      
      curModel=modelList[[i]]
      curModel=predict_all(curModel,geneData,pattern=F)
      pred_score=predictionScore(curModel$loc,cell_loc)
      pred_score=mean(pred_score)
      set.seed(1)
      flds <- createFolds(1:refNum, k = foldNum, list = TRUE, returnTrain = FALSE)
      
      score=c()
      for(k in 1:foldNum){
        index=flds[[k]]
        curModel=modelList[[i]]
        geneData1$refNum=refNum-length(index)
        geneData1$insitu=insitu[,-index]
        geneData1$drop=rbind(drop_ref[-index,],drop_ref[index,,drop=F],drop_nonref)
        curModel=normalize(curModel,geneData1)
        curModel=compute_dist(curModel)
        curModel=predict_pattern(curModel,patternInd=1,gene.end=refNum,gene.start = refNum-length(index)+1)
        score=c(score,patternScore_SS(curModel$pattern,insitu[,index,drop=F],length(index),1))
        message(k)
      }
      score=mean(score)
      
      performance=data.frame(pattern=score,prediction=pred_score)
      
      performance
    }
  performance$pattern=round(performance$pattern,3)
  #performance$pattern_score_test=round(performance$pattern_score_test,3)
  performance$prediction=round(performance$prediction,3)
  return(performance)
}

