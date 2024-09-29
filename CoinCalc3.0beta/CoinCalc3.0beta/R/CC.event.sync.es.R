# Autor: Leonna Szangolies, Potsdam Institute for Climate Impact Research, 2020
# This function is part of the R package CoinCalc3.0().
# CC.event.sync.es performes event synchronisation with local delt including a shuffle significance tests for event series.

CC.event.sync.es<-function(a,b,spana,spanb,max_delT=min(spana[2],spanb[2]),tau=0,sigtest=TRUE,reps=1000,alpha=0.05){
  if(!is.vector(a)){
    print("|-------------    ERROR #2    -------------|")
    print("|         series A is not a vector         |")
    print("|------------------------------------------|")
    return()}
  # ---- #3: seriesB is not a vector
  if(!is.vector(b)){
    print("|-------------    ERROR #3    -------------|")
    print("|         series B is not a vector         |")
    print("|------------------------------------------|")
    return()}
  # ---- #4: time spans do not match
  if(spana[1]>spanb[2] || spanb[1]>spana[2]){
    print("|-------------    ERROR #4    -------------|")
    print("|         no common time span found        |")
    print("|          for spanA and spanB             |")       
    print("|------------------------------------------|")
    return()} 
  # ---- #6: list format
  if(is.list(a)|| is.list(b)){
    print("|-------------    ERROR #6    -------------|")
    print("|  seriesA or seriesB seems to be a list   |")
    print("|          but has to be a vector          |")       
    print("|------------------------------------------|")
    return()} 
  ta<-a
  tb<-b
  
  if(length(ta)<3 || length(tb)<3){
    print("|-------------    ERROR #3    -------------|")
    print("| At least 3 events needed in each series. |")
    print("|------------------------------------------|")
    return()
  } 
  
  span = (max(spana[1],spanb[1]):min(spana[2],spanb[2]))
  
  delT=matrix(nrow=length(ta),ncol=length(tb))
  
  sum_ab=0
  for(i in 2:(length(ta)-1)){
    if(!is.element(ta[i],span-tau)){next}
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
    if(!is.element(tb[i],span-tau)){next}
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
  
  
  if (sigtest == TRUE) {#shuffle test
    surdist = matrix(NA, 2, reps)
    for (surno in 1:reps) {
      surA = sample(c(spana[1]:spana[2]), size = length(ta))
      surB = sample(c(spanb[1]:spanb[2]), size = length(tb))
      K_Q_sur = 0
      K_q_sur = 0
      for (step_b in 2:(length(tb)-1)) {
        surA_small=surA[2:(length(ta)-1)]
        if(tau>=surB[step_b]){next}
        new1=(0.5*sum((surB[step_b]  - surA_small[is.element(surA_small,span-tau)] - tau)>0 && (surB[step_b]  - surA_small[is.element(surA_small,span-tau)]-tau)<delT[c(2:(length(ta)-1))[is.element(surA_small,span-tau)],step_b])+sum((surB[step_b]-surA_small[is.element(surA_small,span-tau)]-tau)==0))
        K_Q_sur=K_Q_sur+new1
        K_q_sur=K_q_sur+new1
      }
      for (step_a in 2:(length(ta)-1)) {
        surB_small=surB[2:(length(tb)-1)]
        if(tau>=surA[step_a]){next}
        new2=(0.5*sum((surA[step_a]  - surB_small[is.element(surB_small,span-tau)] - tau)>0 && (surA[step_a]  - surB_small[is.element(surB_small,span-tau)]-tau)<delT[step_a,c(2:(length(tb)-1))[is.element(surB_small,span-tau)]])+sum((surA[step_a]-surB_small[is.element(surB_small,span-tau)]-tau)==0))
        K_Q_sur=K_Q_sur+new2
        K_q_sur=K_q_sur-new2
      }
      surdist[1, surno] = (K_Q_sur/sqrt((length(ta)-2)*(length(tb)-2)))
      surdist[2, surno] = (K_q_sur/sqrt((length(ta)-2)*(length(tb)-2)))
    }
    PQ = 1 - ecdf(surdist[1, ])(Q)
    Pq = 1 - ecdf(surdist[2, ])(q)
    
    sig_testQ = logical
    if (PQ >= alpha) {
      sig_testQ = TRUE
    }
    else {
      sig_testQ = FALSE
    }
    sig_testq = logical
    if (Pq >= alpha) {
      sig_testq = TRUE
    }
    else {
      sig_testq = FALSE
    }
    CA_out = list(sig_testQ, sig_testq, PQ, Pq, 
                  Q, q)
    names(CA_out) = c("NH Q", "NH q", "p-value Q", 
                      "p-value q", "Q", "q")
  }else{
    CA_out = list(Q, q)
    names(CA_out) = c("Q", "q")
  }
  return(CA_out)
}