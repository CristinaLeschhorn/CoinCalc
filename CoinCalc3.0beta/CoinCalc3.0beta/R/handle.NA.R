# Author: Leonna Szangolies, Potsdam Institute for Climate Impact Research, 2020
# This function is part of the R package CoinCalc3.0().
# This function transforms timeseries which contain NAs so that shuffle and surrogate test can be applied in ECA.

handle.NA<-function(seriesA,seriesB,delT,tau){
if(any(is.na(seriesB))){
  if(delT==0){
    seriesB[which(is.na(seriesB))] = 0
    seriesA[which(is.na(seriesB))-tau] = 0
  }
  else{
    bnew<-c()
    for(i in 1:length(seriesB)){
      btest<-c()
      for(j in 0:delT){
        btest<-c(btest,seriesB[i+j])
      }
      bsum<-sum(btest,na.rm=TRUE)
      if(is.numeric(bsum) && bsum>0){
        bnew[i]<-1
      } else{
        if(any(is.na(btest))){
          bnew[i]<-NA
        }
        else{
          bnew[i]<-0
        }
      }
    }
    delT<-0
    seriesB<-bnew
    seriesB[which(is.na(seriesB))] = 0
    seriesA[which(is.na(seriesB))-tau] = 0
  }
}
  if(any(is.na(seriesA))){
    if(delT==0){
      seriesA[which(is.na(seriesA))] = 0
      seriesB[which(is.na(seriesA))-tau] = 0
    }
    else{
      anew<-c()
      for(i in 1:length(seriesA)){
        atest<-c()
        for(j in 0:delT){
          atest<-c(atest,seriesA[i+j])
        }
        asum<-sum(atest,na.rm=TRUE)
        if(is.numeric(asum) && asum>0){
          anew[i]<-1
        } else{
          if(any(is.na(atest))){
            anew[i]<-NA
          }
          else{
            anew[i]<-0
          }
        }
      }
      delT<-0
      seriesA<-anew
      seriesA[which(is.na(seriesA))] = 0
      seriesB[which(is.na(seriesA))-tau] = 0
    }
  }
  CA_out=list(seriesA,seriesB)
  return(CA_out)}