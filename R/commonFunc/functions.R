
intensityPlot<-function(gene_level,geometry,title="",side="xz",col.range=c(3000,10000)){
  if(title=="")
    title=as.character(substitute(gene_level))
  col.lowest=col.range[1]
  col.highest=col.range[2]
  mydata=data.frame(x=geometry$x,y=geometry$y,z=geometry$z,gene=gene_level)
  mydata$gene=rank(mydata$gene)
  mydata$gene=round((mydata$gene-min(mydata$gene))*(col.highest-col.lowest-1)/(max(mydata$gene)-min(mydata$gene))+1)
  data.color=rev(heat.colors(col.highest))[col.lowest:col.highest]
  mydata$Col <- data.color[mydata$gene]
  if(side=="xz")
    plot(mydata$x,mydata$z,pch = 20,col = mydata$Col,xlab="x",ylab="z",main=title)
  if(side=="zy")
    plot(mydata$x,mydata$y,pch = 20,col = mydata$Col,xlab="x",ylab="y",main=title)
}
