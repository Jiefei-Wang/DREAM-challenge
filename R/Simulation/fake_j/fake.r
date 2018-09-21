
sim_circle<-function(geometry,parms){
#large dispersion means large range
dispersion_x=parms$dispersion_x
dispersion_z=parms$dispersion_z
signal=parms$signal
x=geometry$x
z=geometry$z
anchor_point=c(sample(x,1),sample(z,1))
scale_distance=1/(sqrt(c(var(x),var(z)))*c(dispersion_x,dispersion_z))

geometry$intensity=0
for(i in 1:length(x)){
  distance=(anchor_point-c(x[i],z[i]))*scale_distance
  distance=sum(distance^2)
  weight=pnorm(distance)
  geometry$intensity[i]=2*min(weight,1-weight)*signal
}
return(geometry)
}

getParms<-function(dispersion_x=0.8,dispersion_z=0.3,signal=1){
  parms=list(dispersion_x=dispersion_x,dispersion_z=dispersion_z,signal=signal)
  return(parms)
}

parms=getParms()
result=sim_circle(geometry,parms)
intensityPlot(result$intensity,result)
