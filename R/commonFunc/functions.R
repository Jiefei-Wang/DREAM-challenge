# 
# intensityPlot<-function(gene_level,geometry,plot_quantile=0.5,title="",side="xz",col.range=c(3000,10000)){
#   if(title=="")
#     title=as.character(substitute(gene_level))
#   col.lowest=col.range[1]
#   col.highest=col.range[2]
#   mydata=data.frame(x=geometry$x,y=geometry$y,z=geometry$z,gene=gene_level)
#   mydata$gene=rank(mydata$gene)
#   mydata$gene=round((mydata$gene-min(mydata$gene))*(col.highest-col.lowest-1)/(max(mydata$gene)-min(mydata$gene))+1)
#   data.color=rev(heat.colors(col.highest))[col.lowest:col.highest]
#   mydata$Col <- data.color[mydata$gene]
#   if(side=="xz")
#     plot(mydata$x,mydata$z,pch = 20,col = mydata$Col,xlab="x",ylab="z",main=title)
#   if(side=="zy")
#     plot(mydata$x,mydata$y,pch = 20,col = mydata$Col,xlab="x",ylab="y",main=title)
# }

intensityPlot2<-function(gene_level,geometry,plot_quantile=0.5,title="",side="xz",col.range=c(80,256)){
  if(title=="")
    title=as.character(substitute(gene_level))
  col.lowest=col.range[1]
  col.lowest05=col.range[1]+(col.range[2]-col.range[1])*0.2
  col.highest=col.range[2]
  cutoff=0.0001
  gene_level_nonzero=gene_level[gene_level>=cutoff]
  gene_level_quantile_5=quantile(gene_level_nonzero,0.05)
  gene_level_quantile_95=quantile(gene_level_nonzero,0.95)
  gene.col=rep(1,length(gene_level))
  gene.col[gene_level>=cutoff]=col.lowest05-col.lowest+(gene_level_nonzero-gene_level_quantile_5)/(gene_level_quantile_95-gene_level_quantile_5)*(col.highest-col.lowest05)
  gene.col=round(gene.col)
  gene.col[gene.col<1]=1
  gene.col[gene.col>col.highest-col.lowest]=col.highest-col.lowest
  
  data.color=rev(heat.colors(col.highest))[col.lowest:col.highest]
  if(side=="xz")
    plot(geometry$x,geometry$z,pch = 20,col = data.color[gene.col],xlab="x",ylab="z",main=title)
  if(side=="zy")
    plot(geometry$x,geometry$y,pch = 20,col = data.color[gene.col],xlab="x",ylab="y",main=title)
}
