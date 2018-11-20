computePerformance_real<-function(modelList,geneData,cell_loc,parallel){
  performance=c()
  if(parallel){
    #Add the global function into export list
    g=ls(globalenv())
    global_func=c()
    for(i in 1:length(g))
      if(is.function(get(g[i])))
        global_func=c(global_func,g[i])
    export=c(global_func)
    performance=foreach(i=1:length(modelList),.combine= rbind,.multicombine=TRUE,.inorder=FALSE,.export=export,.packages=clusterPkg)%dopar%{
      curModel=modelList[[i]]
      curModel=predict_all(curModel,geneData,pattern=F)
      performance=c()
      patternFuncList=curModel$compute_pattern
      for(j in 1:length(patternFuncList)){
        tmpModel=curModel
        tmpModel=predict_pattern(tmpModel,patternInd=j,tmpModel$refNum)
        tmpModel$p_parm=tmpModel$p_parm[[j]]
        performance=rbind(performance,getSummaryScore_real(tmpModel,cell_loc))
      }
      performance
    }
  }else{
    for(i in 1:length(modelList)){
      curModel=modelList[[i]]
      curModel=predict_all(curModel,geneData,pattern=F)
      patternFuncList=curModel$compute_pattern
      
      for(j in 1:length(patternFuncList)){
        tmpModel=curModel
        tmpModel=predict_pattern(tmpModel,patternInd=j,tmpModel$refNum)
        tmpModel$p_parm=tmpModel$p_parm[[j]]
        performance=rbind(performance,getSummaryScore_real(tmpModel,cell_loc))
      }
    }
  }
  performance$pattern_score_train=round(performance$pattern_score_train,3)
  #performance$pattern_score_test=round(performance$pattern_score_test,3)
  performance$prediction_score=round(performance$prediction_score,3)
  return(performance)
}

getSummaryScore_real<-function(mydata,cell_loc){
  #pattern_score_test=patternScore(mydata$pattern,simulation$patternData$dropTable,mydata$refNum+1,mydata$geneNum)
  geneNum=mydata$pattern_train_num
  if(is.null(geneNum))
    geneNum=mydata$refNum
  pattern_score_train=patternScore_SS(mydata$pattern,mydata$insitu,1,geneNum)
  #pattern_score_test=mean(pattern_score_test,na.rm = T)
  pattern_score_train=mean(pattern_score_train,na.rm = T)
  
  pred_score=predictionScore(mydata$loc,cell_loc)
  pred_score=mean(pred_score,na.rm = T)
  res=data.frame(normalization=mydata$funcName[1],
                 distance=mydata$funcName[2],
                 pattern=mydata$patternModel,
                 pattern_score_train=pattern_score_train,
                 #pattern_score_test=pattern_score_test,
                 prediction_score=pred_score,
                 normalization_parm=paste0(mydata$N_parm,collapse = "+"),
                 distance_parm=paste0(mydata$d_parm,collapse = "+"),
                 pattern_parm=paste0(mydata$p_parm,collapse = "+")
  )
  res
}





#predPattern=mymodel$pattern
#truePattern=mymodel$insitu
#gene_start=2
#gene_end=2

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


patternScore_cor<-function(predPattern,truePattern,gene_start,gene_end){
  turning=0.01
  predPattern=matrix(predPattern[,gene_start:gene_end],nrow(predPattern))
  truePattern=matrix(truePattern[,gene_start:gene_end],nrow(predPattern))
  score=rep(NA,ncol(predPattern))
  for(i in 1:ncol(predPattern)){
    score[i]=cor(predPattern[,i],truePattern[,i],method="spearman")
  }
  return(score)
}

#top10=simulation$cell_loc
#predLoc=apply(top10,2,function(x)sample(x,10))
#trueLoc=simulation$cell_loc
#predLoc=mydata$loc
predictionScore1<-function(predLoc,top10){
 
  ave_x=diff(range(geometry$x))/2
  ave_z=diff(range(geometry$z))/2
  z_inflate=ave_x/ave_z
  weight=1/(1:nrow(predLoc))
  weight=weight/sum(weight)
  geometry_tmp=geometry
  geometry_tmp$z=geometry_tmp$z*z_inflate
  coordinate_x=matrix(geometry_tmp[predLoc,"x"],nrow=10)
  coordinate_z=matrix(geometry_tmp[predLoc,"z"],nrow=10)
  ave_x=diff(range(geometry_tmp$x))
  ave_z=diff(range(geometry_tmp$z))
  ave_dist=sqrt(ave_x^2+ave_z^2)
  score_final=c()
  for(i in 1:10){
    trueLoc=top10[i,]
    
    true_x=geometry_tmp[trueLoc,"x"]
    true_z=geometry_tmp[trueLoc,"z"]
    
    dist_x=sweep(coordinate_x,2,true_x,"-")
    dist_z=sweep(coordinate_z,2,true_z,"-")
    
    dist=sqrt(dist_x^2+dist_z^2)
    
    score=1-dist/ave_dist
    score=colSums(sweep(score,1,weight,"*"))
    score_final=rbind(score_final,score)
  }
  weight=1/(1:nrow(predLoc))+1
  weight=weight/sum(weight)
  weight%*%score_final
  # matchPos=sweep(predLoc,2,trueLoc,"==")
  # ind=which(matchPos,arr.ind = T)
  # sum(1/ind[,1])/length(trueLoc)
}
predictionScore<-function(predLoc,trueLoc){
  ave_x=diff(range(geometry$x))/2
  ave_z=diff(range(geometry$z))/2
  z_inflate=ave_x/ave_z
  
  geometry_tmp=geometry
  geometry_tmp$z=geometry_tmp$z*z_inflate
  coordinate_x=matrix(geometry_tmp[predLoc,"x"],nrow=10)
  coordinate_z=matrix(geometry_tmp[predLoc,"z"],nrow=10)
  true_x=geometry_tmp[trueLoc,"x"]
  true_z=geometry_tmp[trueLoc,"z"]
  
  
  ave_x=diff(range(geometry_tmp$x))
  ave_z=diff(range(geometry_tmp$z))
  ave_dist=sqrt(ave_x^2+ave_z^2)
  
  weight=1/(1:nrow(predLoc))
  weight=weight/sum(weight)
  
  dist_x=sweep(coordinate_x,2,true_x,"-")
  dist_z=sweep(coordinate_z,2,true_z,"-")
  
  dist=sqrt(dist_x^2+dist_z^2)
  
  score=1-dist/ave_dist
  score=colsums(sweep(score,1,weight,"*"))
  
  
  
  # matchPos=sweep(predLoc,2,trueLoc,"==")
  # ind=which(matchPos,arr.ind = T)
  # sum(1/ind[,1])/length(trueLoc)
}
