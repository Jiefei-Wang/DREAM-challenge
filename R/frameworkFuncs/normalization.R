normalize_simple<-function(mydata){
  mydata$N_insitu=t(scale(t(scale(mydata$insitu))))
  mydata$N_drop=t(scale(t(scale(mydata$drop))))
  return(mydata)
}

normalize_columnSum <- function(mydata){
  #divide the column by the colum sums
  mydata$N_drop <- mydata$drop/colSums(mydata$drop)
  mydata$N_insitu <- mydata$insitu
  return(mydata)
}

