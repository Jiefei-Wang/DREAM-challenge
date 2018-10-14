library(doParallel)
for(pkg in clusterPkg)
  library(pkg,character.only=T)

if(!exists("cl",parent.frame())){
  #try(stopCluster(cl))
  cl <- makePSOCKcluster(clusterNum)
  registerDoParallel(cl)
}

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
attachFunc<-function(mydata,normalize,normalize.parm,computeDist,dist.parm,predict.pattern){
  mydata$normalize=normalize
  mydata$N_parm=normalize.parm
  mydata$computeDist=computeDist
  mydata$d_parm=dist.parm
  #This function has been implemented, do not need to be specified by the inputs.
  mydata$pred_loc=predLocation
  mydata$compute_pattern=predict.pattern
  mydata
}
#Same as attachFunc, but the parameter can be an integer
#The attached function's paramter will be automatically find
attachFunc_list<-function(mydata,normalization,distance,pattern){

  model_list=list()
  for(i in 1:length(normalization))
    for(j in 1: length(distance))
      for(k in 1: length(pattern)){
        model_list=c(model_list,attachFunc_hide(mydata,normalization[i],distance[j],pattern[k]))}
  model_list
}
attachFunc_hide<-function(mydata,normalization,distance,pattern){
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
    pattern=paste0("computePattern_",pattern,collapse = "")
  }else
    pattern=paste0("computePattern_",getPatternFuncs()[pattern],collapse = "")

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
  #Make a data list
  mydata_list=list()
  k=1
  for (i in 1:length(normalization_parm)) {
    for (j in 1:length(distance_parm)) {
      mydata_single = attachFunc(mydata,get(normalization),normalization_parm[[i]],
                                 get(distance),distance_parm[[j]],get(pattern))
      mydata_single$funcName=c(gsub("normalize_","",normalization,fixed = T),
                               gsub("computeDist_","",distance,fixed = T),
                               gsub("computePattern_","",pattern,fixed = T))
      mydata_list[[k]] = mydata_single
      k = k + 1
    }
  }
  return(mydata_list)
}



normalize<-function(mydata){
  mydata$normalize(mydata)
}
compute_dist<-function(mydata){
  mydata$computeDist(mydata)
}
pred_loc<-function(mydata){
  mydata$pred_loc(mydata)
}
predict_pattern<-function(mydata,gene){
  mydata$compute_pattern(mydata,gene)
}
predict_all<-function(mydata){
mydata=normalize(mydata)
mydata=compute_dist(mydata)
mydata=pred_loc(mydata)
mydata$pattern=sapply(1:mydata$geneNum,predict_pattern,mydata=mydata)
mydata
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



