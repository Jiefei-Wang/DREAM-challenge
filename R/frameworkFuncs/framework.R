library(doParallel)
for(pkg in clusterPkg)
  library(pkg,character.only=T)
clusterExport(cl,"geometry")
clusterExport(cl,"weighted_cor")

#========================Read functions into a namespace=======================
for(env in search()){
  if(env=="Normalization") detach(env,character.only=T)
  if(env=="Distance") detach(env,character.only=T)
  if(env=="Pattern") detach(env,character.only=T)
  if(env=="Scoring") detach(env,character.only=T)
  if(env=="DataSource") detach(env,character.only=T)
}
rm("env")
Normalization <- new.env()
Distance<-new.env()
Pattern<-new.env()
Scoring<-new.env()
DataSource<-new.env()
sys.source("R\\frameworkFuncs\\normalization.R",envir=Normalization)
sys.source("R\\frameworkFuncs\\distance.R",envir=Distance)
sys.source("R\\frameworkFuncs\\pattern.R",envir=Pattern)
sys.source("R\\frameworkFuncs\\scoring.R",envir=Scoring)
sys.source("R\\simulation\\generate_fake.R",envir=DataSource)
sys.source("R\\Simulation\\fake\\fake.R",envir=DataSource)
attach(Normalization)
attach(Distance)
attach(Pattern)
attach(Scoring)
attach(DataSource)
rm("Normalization")
rm("Distance")
rm("Pattern")
rm("Scoring")
rm("DataSource")



#========================Framework=======================
#Create the dataset associated with the given functions, and function parameters
attachFunc<-function(normalize,normalize.parm,computeDist,dist.parm,predict.pattern,pattern.parm){
  model=list()
  model$normalize=normalize
  model$N_parm=normalize.parm
  model$computeDist=computeDist
  model$d_parm=dist.parm
  #This function has been implemented, do not need to be specified by the inputs.
  model$pred_loc=predLocation
  model$compute_pattern=predict.pattern
  model$p_parm=pattern.parm
  model
}
#Same as attachFunc, but the parameter can be an integer
#The attached function's paramter will be automatically find
attachFunc_list<-function(normalization,distance,pattern){

  model_list=list()
  for(i in 1:length(normalization))
    for(j in 1: length(distance)){
        model_list=c(model_list,attachFunc_hide(normalization[i],distance[j],pattern))}
  model_list
}
attachFunc_hide<-function(normalization,distance,pattern){
  #Check the input
  if(is.character(normalization)){
    normalization=paste0("normalize_",normalization,collapse = "")
  }else
    normalization=paste0("normalize_",getNormFuncs()[normalization],collapse = "")
  
  if(is.character(distance)){
    distance=paste0("computeDist_",distance,collapse = "")
  }else
    distance=paste0("computeDist_",getDistFuncs()[distance],collapse = "")
  
  if(is.character(pattern)){
    pattern=paste0("computePattern_",pattern)
  }else
    pattern=paste0("computePattern_",getPatternFuncs()[pattern])

  #Obtain the function parameters
  normalization_parm_fun=getFunParms(normalization)
  if(is.null(normalization_parm_fun)){
    normalization_parm=list(0)
  }else
    normalization_parm=do.call(normalization_parm_fun,args = list())
  distance_parm_fun=getFunParms(distance)
  if(is.null(distance_parm_fun)){
    distance_parm=list(0)
  }else
    distance_parm=do.call(distance_parm_fun,args = list())
  pattern_parm_list=list()
  pattern_list=c()
  pattern_name_list=c()
  for(i in 1:length(pattern)){
    pattern_parm_fun=getFunParms(pattern[i])
    if(is.null(pattern_parm_fun)){
      pattern_parm_list=c(pattern_parm_list,list(0))
      pattern_list=c(pattern_list,get(pattern[i]))
      pattern_name_list=c(pattern_name_list,pattern[i])
    }else{
      pattern_parm=do.call(pattern_parm_fun,args = list())
      pattern_parm_list=c(pattern_parm_list,pattern_parm)
      for(j in 1:length(pattern_parm)){
        pattern_list=c(pattern_list,get(pattern[i]))
        pattern_name_list=c(pattern_name_list,pattern[i])
      }
    }
  }
  
  
  #Make a data list
  model_list=list()
  k=1
  for (i in 1:length(normalization_parm)) {
    for (j in 1:length(distance_parm)) {
      single_model = attachFunc(get(normalization),normalization_parm[[i]],
                                 get(distance),distance_parm[[j]],pattern_list,pattern_parm_list)
      single_model$funcName=c(gsub("normalize_","",normalization,fixed = T),
                               gsub("computeDist_","",distance,fixed = T),
                               gsub("computePattern_","",pattern_name_list,fixed = T))
      model_list[[k]] = single_model
      k = k + 1
    }
  }
  return(model_list)
}



normalize<-function(model,geneData){
  mydata=c(geneData,model)
  mydata$normalize(mydata)
}
compute_dist<-function(mydata){
  mydata$computeDist(mydata)
}
pred_loc<-function(mydata){
  mydata$pred_loc(mydata)
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

predict_pattern<-function(mydata,patternInd,gene.end=mydata$geneNum,gene.start=1){
  patternFuncList=mydata$compute_pattern
  parmList=mydata$p_parm
  mydata$compute_pattern=patternFuncList[[patternInd]]
  mydata$p_parm=parmList[[patternInd]]
  for(i in gene.start:gene.end){
    mydata=mydata$compute_pattern(mydata,i)
  }
  mydata$patternModel=paste0(mydata$funcName[3+patternInd-1],"(",mydata$p_parm,")")
  mydata$compute_pattern=patternFuncList
  mydata$p_parm=parmList
  mydata
}
predict_all<-function(model,geneData,pattern=T,patternInd=1){
mydata=normalize(model,geneData)
mydata=compute_dist(mydata)
mydata=pred_loc(mydata)
if(!pattern)
  return(mydata)
  
mydata=predict_pattern(mydata,patternInd)
return(mydata)
}


#========================Show all the available functions=======================
getNormFuncs<-function(){
  gsub("normalize_","",getFuncList("Normalization"),fixed = T)
}
getDistFuncs<-function(){
  gsub("computeDist_","",getFuncList("Distance"),fixed = T)
}

getPatternFuncs<-function(){
  gsub("computePattern_","",getFuncList("Pattern"),fixed = T)

}
getFuncList<-function(spaceName){
  funcs=ls(envir=as.environment(spaceName))
  for(i in 1:length(funcs)){
    ind=grep("get_parm_",funcs[i],fixed=T)
    if(length(ind)!=0)
      funcs[i]=""
  }
  funcs[funcs!=""]
}
getFunParms<-function(FuncName){
  spaceNames=c("Normalization","Distance","Pattern")
  funcBegining=c("normalize_","computeDist_","computePattern_")
  for(i in 1:length(spaceNames)){
    funcs=ls(envir=as.environment(spaceNames[i]))
    if(FuncName %in% funcs){
      specialName=gsub(funcBegining[i],"get_parm_",FuncName,fixed=T)
      if(specialName %in% funcs)
        return(specialName)
      return(NULL)
    }
  }
  return(NULL)
}



