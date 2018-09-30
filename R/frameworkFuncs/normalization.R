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

normalize_rowMax <- function(mydata){
  rowmax <- as.integer(apply(mydata$drop, 1, function(x) max(x)))
  drop <- mydata$drop
  insitu <- mydata$insitu
  
  #take the maximum of each row
  rowmax <- structure(list(V1 = c(rowmax)), .Names = "V1", class = "data.frame")
  
  #divide each row by the row maximum and get normalized
  mydata$N_drop <- mapply(function(x, y) x/y, drop, rowmax)
  mydata$N_insitu <- mydata$insitu
  
  
  return(mydata)
}

