# Author: Leonna Szangolies, Potsdam Institute for Climate Impact Research, 2020
# This function is part of the R package CoinCalc3.0().
# CC.for.visroc performes part of a Receiver Operating Characteristics analysis and returns the Output needed for the CC.visroc analysis

require(MESS)
CC.for.visroc<-function(data1,data2, thres1, thres2=seq(from=0.05,to=0.95, by=0.05)){
  visrocvector<-list()
  tpmat<-vector(length=length(thres2))
  fpmat<-vector(length=length(thres2))
  stop=c()
  first<-c()
  for(i in 1:length(data1)){
    if(is.na(data1[i])){
      stop=c(stop,i)
      first[i]=0
      next
    }
    q<-quantile(data1,probs=c(thres1),na.rm=TRUE)
    if(data1[i]>=q[1]){
      first[i]=1
    }else{
      first[i]=0
    }
  }
  k=0
  visrocvector[["p"]]<-sum(first)
  visrocvector[["q"]]<-length(first)-sum(first)-length(stop)
  visrocvector[["N"]]<-length(first)
  q2<-quantile(data2,probs= thres2,na.rm=TRUE)
  for(thresi in c(1:length(thres2))){
    thres<-q2[thresi]
    k=k+1
    tp=0
    fp=0
    for(i in c(1:length(data1))[! c(1:length(data1)) %in% stop]){
      if(is.na(data2[i])){ 
          stop=c(stop,i)
          first[i]=0
          next
        }
      if(data2[i]< thres){
        if(first[i]==1){
          tp=tp+1
        }
        else{
          fp=fp+1
        }
      }
    }
    tpmat[k]<-tp/sum(first)
    fpmat[k]<-fp/(length(first)-sum(first)-length(stop))
  }
  m<-matrix(c(fpmat,tpmat),ncol=2)
  
  d<-sqrt((1-tpmat)^2+(0-fpmat)^2)
  posopt<-which.min(d)
  
  auc<-try(auc(fpmat[!is.na(fpmat)],tpmat[!is.na(tpmat)],type="spline"))
  visrocvector[["uauc"]]<-auc
  visrocvector[["f1"]]<-tpmat[posopt]
  visrocvector[["h1"]]<-fpmat[posopt]
  return(visrocvector)
}
