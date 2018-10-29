if(exists("cl",parent.frame())){
  try(stopCluster(cl))
}
cl <- makePSOCKcluster(clusterNum)

registerDoParallel(cl)