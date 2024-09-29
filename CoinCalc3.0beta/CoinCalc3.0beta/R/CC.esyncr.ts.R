# Autor: Leonna Szangolies, Potsdam Institute for Climate Impact Research, 2020
# This function is part of the R package CoinCalc3.0().
# esyncr.ts performes event synchronisation vor varying tau using time series and some parameters. The significance level using the shuffle tests can be computed.

CC.esyncr.ts<-function(a,b=a,alpha=0.05,max_delT=min(length(a),length(b)),sym=FALSE,offset=30,siglevel=TRUE,reps=1000){
  ecrprec<-vector(length=offset)
  ecrtrig<-vector(length=offset)
  ecrprecsig<-vector(length=offset)
  ecrtrigsig<-vector(length=offset)
  for(k in 0:(offset-1)){
    tau=k
    if (length(table(a)) > 2) {
      print("|-------------    ERROR #1    -------------|")
      print("|       Time seriesA is not binary         |")
      print("|      (or the number of events=0)!        |")
      print("| Use CC.binarize() to preprocess seriesA  |")
      print("|------------------------------------------|")
      return()
    }
    if (length(table(b)) > 2) {
      print("|-------------    ERROR #1    -------------|")
      print("|       Time seriesB is not binary         |")
      print("|      (or the number of events=0)!        |")
      print("| Use CC.binarize() to preprocess seriesB  |")
      print("|------------------------------------------|")
      return()
    }
    if(length(a)!=length(b)){
      print("|-------------    ERROR #2    -------------|")
      print("|    lengths of series A and B differ      |")
      print("|------------------------------------------|")
      return()
    } 
    ta<-which(a==1)
    tb<-which(b==1)
    
    if(length(ta)<3 || length(tb)<3){
      print("|-------------    ERROR #3    -------------|")
      print("| At least 3 events needed in each series. |")
      print("|------------------------------------------|")
      return()
    } 
    
    delT=matrix(nrow=length(ta),ncol=length(tb))
    
    sum_ab=0
    for(i in 2:(length(ta)-1)){
      for(j in 2:(length(tb)-1)){
        delT[i,j]=min(min(ta[i+1]-ta[i],ta[i]-ta[i-1],tb[j+1]-tb[j],tb[j]-tb[j-1],na.rm=TRUE)/2,max_delT)
        if(tau>=tb[j]){next}
        if(0<(ta[i]-tau-tb[j]) && (ta[i]-tau-tb[j])<=delT[i,j]){
          sum_ab=sum_ab+1
        }else{
          if((ta[i]-tau)==tb[j]){
            sum_ab=sum_ab+0.5
          }
        }
      }
    }
    sum_ba=0
    for(i in 2:(length(tb)-1)){
      for(j in 2:(length(ta)-1)){
        if(tau>=ta[j]){next}
        if(0<(tb[i]-tau-ta[j]) && (tb[i]-tau-ta[j])<=delT[j,i]){
          sum_ba=sum_ba+1
        }else{
          if((tb[i]-tau)==ta[j]){
            sum_ba=sum_ba+0.5
          }
        }
      }
    }
    Q=(sum_ba+sum_ab)/sqrt((length(ta)-2)*(length(tb)-2)) #Q=1 -> fully synchronized
    q=(sum_ba-sum_ab)/sqrt((length(ta)-2)*(length(tb)-2)) #q=1 -> x events always precede y events
    
    ecrprec[k+1]<-Q
    ecrtrig[k+1]<-q
    
    if (siglevel == TRUE) {#shuffle test
      surdist = matrix(NA, 2, reps)
      span = seq(1:length(a))
      for (surno in 1:reps) {
        surA = sample(span, size = length(ta))
        surB = sample(span, size = length(tb))
        K_Q_sur = 0
        K_q_sur = 0
        for (step_b in 2:(length(tb)-1)) {
          if(tau>=surB[step_b]){next}
          new1=(0.5*sum((surB[step_b] - surA[2:(length(ta)-1)]-tau)>0 && (surB[step_b] - surA[2:(length(ta)-1)]-tau)<delT[2:(length(ta)-1),step_b])+sum((surB[step_b]-surA[2:(length(ta)-1)]-tau)==0))
          K_Q_sur=K_Q_sur+new1
          K_q_sur=K_q_sur+new1
        }
        for (step_a in 2:(length(ta)-1)) {
          new2=(0.5*sum((surA[step_a] - surB[2:(length(tb)-1)]-tau)>0 && (surA[step_a] - surB[2:(length(tb)-1)]-tau)<delT[step_a,2:(length(tb)-1)])+sum((surA[step_a]-surB[2:(length(tb)-1)]-tau)==0))
          K_Q_sur=K_Q_sur+new2
          K_q_sur=K_q_sur-new2
        }
        surdist[1, surno] = (K_Q_sur/sqrt((length(ta)-2)*(length(tb)-2)))
        surdist[2, surno] = (K_q_sur/sqrt((length(ta)-2)*(length(tb)-2)))
      }
      Pprecsig = quantile(surdist[1, ],1-alpha)
      Ptriggsig = quantile(surdist[2, ],1-alpha) 
      
      ecrprecsig[k+1]<-Pprecsig
      ecrtrigsig[k+1]<-Ptriggsig
    }
  }
  CR_out=list(c(0:(offset-1)),ecrprec,ecrtrig,ecrprecsig,ecrtrigsig)
  names(CR_out)=c("tau","Q","q","Q significance level","q significance level")
  return(CR_out)
}