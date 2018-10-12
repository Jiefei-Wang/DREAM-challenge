#install.packages("Rfast")
#install.packages("rpgm")



normalize_scale<-function(mydata){
  mydata$N_insitu=t(scale(t(scale(mydata$insitu))))
  mydata$N_drop=t(scale(t(scale(mydata$drop))))
  return(mydata)
}

normalize_columnSum <- function(mydata){
  mydata$N_drop <- mydata$drop/colSums(mydata$drop)
  mydata$N_insitu <- mydata$insitu
  return(mydata)
}


normalize_rowMaxlog <- function(mydata){
  
  
  mydata$N_drop <- log2(1+mydata$drop/apply(mydata$drop,1,max))
  
  mydata$N_insitu <- log2(1+mydata$insitu)
  
  
  return(mydata)
}
normalize_rowMax <- function(mydata){
  
  
  mydata$N_drop <- mydata$drop/apply(mydata$drop,1,max)
  
  mydata$N_insitu <- mydata$insitu
  
  
  return(mydata)
}
normalize_TMM <- function(mydata){
  normfac_drop <- calcNormFactors(mydata$drop)
  
  mydata$N_drop<- sweep(mydata$drop,2,normfac_drop,"/")
  
  normfac_insitu <- calcNormFactors(mydata$insitu)
  
  mydata$N_insitu <-sweep(mydata$insitu,2,normfac_insitu,"/")
  
                        
 return(mydata)
}



normalize_upqu <- function(mydata){
  parm=mydata$N_parm
  #insitu
  quantileExpressed <- apply(mydata$insitu, 1, function(x){quantile(x[x>0], parm)})
  mydata$N_insitu <- sweep(mydata$insitu,1,quantileExpressed,"/") 
  #dropseq
  quantileExpressed <- apply(mydata$drop, 1, function(x){quantile(x[x>0], parm)})
  mydata$N_drop <- sweep(mydata$drop,1,quantileExpressed,"/")
  return(mydata)  
}
get_parm_upqu <- function(param){
  n=40
  return(as.list(seq(1,n-1,by=2)/n))
}

