computePerformance<-function(modelList,geneData,simulation,parallel){
  performance=data.frame()
  if(parallel){
    #Add the global function into export list
    g=ls(globalenv())
    global_func=c()
    for(i in 1:length(g))
      if(is.function(get(g[i])))
        global_func=c(global_func,g[i])
    export=c(ls(as.environment("Scoring")),global_func)
    performance=foreach(i=1:length(modelList),.combine= rbind,.multicombine=TRUE,.inorder=FALSE,.export=export,.packages=clusterPkg)%dopar%{
      curModel=modelList[[i]]
      curModel=predict_all(curModel,geneData,pattern=F)
      res=c()
      patternFuncList=curModel$compute_pattern
      for(j in 1:length(patternFuncList)){
        curModel$tmp_pattern=NULL
        curModel=predict_pattern(curModel,patternInd=j,curModel$refNum)
        res=rbind(res,getSummaryScore(curModel,simulation))
      }
      res
    }
  }else{
    performance=c()
    for(i in 1:length(modelList)){
      curModel=modelList[[i]]
      curModel=predict_all(curModel,geneData,pattern=F)
      patternFuncList=curModel$compute_pattern
      
      for(j in 1:length(patternFuncList)){
        curModel$tmp_pattern=NULL
        curModel=predict_pattern(curModel,patternInd=j,curModel$refNum)
        performance=rbind(performance,getSummaryScore(curModel,simulation))
      }
    }
  }
  performance$pattern_score_train=round(performance$pattern_score_train,3)
  #performance$pattern_score_test=round(performance$pattern_score_test,3)
  performance$prediction_score=round(performance$prediction_score,3)
  return(performance)
}

getSummaryScore<-function(mydata,simulation){
  #pattern_score_test=patternScore(mydata$pattern,simulation$patternData$dropTable,mydata$refNum+1,mydata$geneNum)
  pattern_score_train=patternScore(mydata$pattern,mydata$insitu,1,mydata$refNum)
  #pattern_score_test=mean(pattern_score_test,na.rm = T)
  pattern_score_train=mean(pattern_score_train,na.rm = T)
  
  pred_score=predictionScore(mydata$loc,simulation$cell_loc)
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

#predPattern=mydata$pattern
#truePattern=simulation$patternData$dropTable
#gene_start=1
#gene_end=84

patternScore1<-function(predPattern,truePattern,gene_start,gene_end){
  turning=0.01
  weight=1
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


patternScore<-function(predPattern,truePattern,gene_start,gene_end){
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




#trueLoc=simulation$cell_loc
#predLoc=mydata$loc
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