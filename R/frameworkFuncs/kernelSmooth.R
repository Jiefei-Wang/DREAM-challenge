
install.packages("np")
library("np")
clusterEvalQ(cl,library("np"))
data=data.frame(dist=dist[,1],x=geometry$x,z=geometry$z)
bw <- npregbw(dist~x+z,data)
y=dist[,1]
ghat <- npreg(bws=c(8,2),txdat=x,tydat=y)
yhat=predict(ghat,newdata=x)

dist=mymodel$distance
x=data.frame(x=geometry$x,z=geometry$z)
dist1=parApply(cl,dist,2,function(y,x_geo){
  ghat <- npreg(bws=c(7.6,1.7),txdat=x_geo,tydat=y)
  predict(ghat,newdata=x_geo)
},x_geo=x)





#take a look at prediction accuracy
ind=c()
for(i in 1:cellNum){
  ind[i]=rank(mymodel$distance[,i])[simulation$cell_loc[i]]
}
median(ind)
mean(ind)


#take a look at prediction accuracy
ind1=c()
for(i in 1:cellNum){
  ind1[i]=rank(dist1[,i])[simulation$cell_loc[i]]
}
median(ind1)
mean(ind1)

y=dist[,2]
ghat <- npreg(bws=c(10,5),txdat=x_geo,tydat=y)
pred=predict(ghat,newdata=x_geo)
rank(pred)[simulation$cell_loc[2]]

cell=2
loc=simulation$cell_loc[cell]
intensityPlot2(max(mymodel$distance[,cell])-mymodel$distance[,cell],geometry,title="Dist")
points(geometry$x[loc],geometry$z[loc],col="blue",pch=20,cex=1)


intensityPlot2(max(dist1[,cell])-dist1[,cell],geometry,title="Dist")
points(geometry$x[loc],geometry$z[loc],col="blue",pch=20,cex=1)
