library(ggplot2)

result$model=gsub("_"," ",factor(paste(result$normalization,result$distance,sep="+")),fixed=T)
ggplot(result, aes(x=prediction_score, y=pattern_score_train,color=model)) + 
  geom_point()+guides(color=guide_legend(ncol=2))+labs(x="Prediction score", y="Pattern score")



model=as.character(unique(result$model))
record=c()
for(i in 1:length(model)){
  curModel=model[i]
  tmp=result[result$model==curModel,]
  #if(nrow(tmp)==1)
  #  next
  ind=which.max(tmp$pattern_score_train)
  #
 res=data.frame(model=curModel,predict_fit=tmp[ind,]$prediction_score,predict_best=max(tmp$prediction_score))
  #res=data.frame(model=curModel,predict_best=max(tmp$prediction_score))
  #res=data.frame(model=curModel,predict_ave=mean(tmp$prediction_score),predict_fit=tmp[ind,]$prediction_score,predict_best=max(tmp$prediction_score))
  record=rbind(record,res)
}

record