
result$model=factor(paste(result$normalization,result$distance,sep="+"))
result1=result[result$pattern=="simple1(0)",]
ggplot(result1, aes(x=prediction_score, y=pattern_score_train,color=model)) + geom_point()+labs(x="Prediction score",y="Pattern score")


ind=which(result1$prediction_score<0.05&result1$pattern_score_test>0.55)
result1[3,]


library(Hmisc)
summarize(as.matrix(result1[,c("prediction_score","pattern_score_test","pattern_score_train")]),result1$model,max)

model=unique(result1$model)
record1=c()
for(i in 1:length(model)){
  tmp=result1[result1$model==model[i],]
  if(nrow(tmp)==1)
    next
  tmp1=tmp$prediction_score
  Ind=which(tmp1==max(tmp1))
  record1=rbind(record1,tmp[Ind[1],])
}



model=unique(result1$model)
record2=c()
for(i in 1:length(model)){
  tmp=result1[result1$model==model[i],]
  if(nrow(tmp)==1)
    next
  tmp1=tmp$pattern_score_train
  Ind=which(tmp1==max(tmp1))
  record2=rbind(record2,tmp[Ind[1],])
}

record3=data.frame(optimal=record2$prediction_score,Best=record1$prediction_score,model=record1$model)



tmp=result[result$model==model[22],]
tmp=tmp[order(rownames(tmp)),]
n=80
tmp$parm=seq(1,n-1)/n

ggplot(tmp,aes(x=parm, y=prediction_score))+ geom_point()+labs(x="Parameters",y="Prediction score")






