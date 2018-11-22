
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
}

load("R\\commonFunc\\GeneData.RData")
load("R\\commonFunc\\authorBinaryData.RData")

geometry=chanllege.data$geometry
pred_score=predictionScore(mymodel$loc,simulation$cell_loc)
mean(pred_score)
