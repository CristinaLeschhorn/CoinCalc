# Author: Leonna Szangolies, Potsdam Institute for Climate Impact Research, 2020
# This function is part of the R package CoinCalc3.0().
# CC.eca.aggregated.es performes aggregated event coincidence analysis including two different significance tests using several pairs of time series and some parameters.

CC.eca.aggregated.es<-function (seriesListA,seriesListB,spanListA,spanListB, alpha = 0.05, delT = 0, sym = FALSE, tau = 0, sigtest = "poisson", reps = 1000) 
{
  if (length(seriesListA) != length(seriesListB)) {
    print("|-------------    ERROR #1    -------------|")
    print("|  Number of event series A and B differ   |")
    print("|------------------------------------------|")
    return()
  }
  if (length(seriesListA) != length(spanListA) || length(seriesListB) != length(spanListB)) {
    print("|-------------    ERROR #1    -------------|")
    print("| Number of event series and spans differ  |")
    print("|------------------------------------------|")
    return()
  }
  if(tau<0 || delT<0){
    print("|-------------    ERROR #2    -------------|")
    print("|    The offset (tau) or delta T (delT)    |")
    print("|              is negative.                |") 
    print("|------------------------------------------|")
    return()
  }
  CRprec=0
  CRtrigg=0
  N_A=vector(length=length(seriesListA))
  N_B=vector(length=length(seriesListB))
  N_A_edge=vector(length=length(seriesListA))
  N_B_edge=vector(length=length(seriesListB))
  Tlen=vector(length=length(seriesListA))
  K_prec=vector(length=length(seriesListA))
  K_trigg=vector(length=length(seriesListA))
  for(currentSeries in 1:length(seriesListA)){
    seriesA=seriesListA[[currentSeries]]
    seriesB=seriesListB[[currentSeries]]
    spanA=spanListA[[currentSeries]]
    spanB=spanListB[[currentSeries]]
    if(!is.vector(seriesA)){
      print("|-------------    ERROR #3    -------------|")
      print("|         series A is not a vector         |")
      print("|------------------------------------------|")
      return()}
    # ---- #3: seriesB is not a vector
    if(!is.vector(seriesB)){
      print("|-------------    ERROR #3    -------------|")
      print("|         series B is not a vector         |")
      print("|------------------------------------------|")
      return()}
    # ---- #4: time spans do not match
    if(spanA[1]>spanB[2] || spanB[1]>spanA[2]){
      print("|-------------    ERROR #4    -------------|")
      print("|         no common time span found        |")
      print("|          for spanA and spanB             |")       
      print("|------------------------------------------|")
      return()} 
    # ---- #6: list format
    if(is.list(seriesA)|| is.list(seriesB)){
      print("|-------------    ERROR #6    -------------|")
      print("|  seriesA or seriesB seems to be a list   |")
      print("|          but has to be a vector          |")       
      print("|------------------------------------------|")
      return()} 
    
    N_A[currentSeries] = length(seriesA)
    N_B[currentSeries] = length(seriesA)
    span=c(max(spanA[1],spanB[1]):min(spanA[2],spanB[2])) # getting common time span of seriesA and B. Will be used for coincidence analysis.
    Tlen[currentSeries] = length(span)
    N_A_edge[currentSeries]=sum(is.element(seriesA,c((spanA[1]-1):(spanA[1]+tau+delT-1))))
    N_B_edge[currentSeries]=sum(is.element(seriesB,c((spanB[2]-(tau+delT)+1):(spanB[2]+1))))
    K_prec[currentSeries] = 0
    for (step_a in (1+N_A_edge[currentSeries]):N_A[currentSeries]) {
      if( !is.element(seriesA[step_a],span) ){next} 
      if (sym == FALSE) {
        K_prec[currentSeries]=K_prec[currentSeries]+as.numeric(any((0<=((seriesA[step_a] - tau) - seriesB[is.element(seriesB,span)])) & (((seriesA[step_a] - tau) - seriesB[is.element(seriesB,span)])<=delT)))
      }
      if (sym == TRUE) {
        K_prec[currentSeries]=K_prec[currentSeries]+as.numeric(any((-delT<=((seriesA[step_a] - tau) - seriesB[is.element(seriesB,span)])) & (((seriesA[step_a] - tau) - seriesB[is.element(seriesB,span)])<=delT)))
      } 
    }
    CRprec = CRprec + K_prec[currentSeries]
    K_trigg[currentSeries] = 0
    for (step_b in 1:(N_B[currentSeries]-N_B_edge[currentSeries])) {
      if( !is.element(seriesB[step_b],span) ){next}
      if (sym == FALSE) {
        K_trigg[currentSeries]=K_trigg[currentSeries]+as.numeric(any((0<=((seriesA[is.element(seriesA,span)] - tau) - seriesB[step_b])) & (((seriesA[is.element(seriesA,span)] - tau) - seriesB[step_b])<=delT)))
      }
      if (sym == TRUE) {
        K_trigg[currentSeries]=K_trigg[currentSeries]+as.numeric(any((-delT<=((seriesA[is.element(seriesA,span)] - tau) - seriesB[step_b])) & (((seriesA[is.element(seriesA,span)] - tau) - seriesB[step_b])<=delT)))
      }
    }
    CRtrigg = CRtrigg + K_trigg[currentSeries]
  }
  CRprec=CRprec/(sum(N_A)-sum(N_A_edge))
  CRtrigg=CRtrigg/(sum(N_B)-sum(N_B_edge))
  
  if (sigtest == "poisson") {
    Pprec_v = rep(0,length=length(seriesListA))
    Ptrigg_v = rep(0,length=length(seriesListA))
    for(currentSeries in 1:length(seriesListA)){
      if (sym == FALSE) {
        for (Ktmp in K_prec[currentSeries]:(N_A[currentSeries]-N_A_edge[currentSeries])) {
          Ptmp = choose((N_A[currentSeries]-N_A_edge[currentSeries]), Ktmp) * ((1 - ((1 - ((delT + 1)/(Tlen[currentSeries] - tau)))^N_B[currentSeries]))^Ktmp) * (((1 - ((delT + 1)/(Tlen[currentSeries] - tau)))^N_B[currentSeries])^(N_A[currentSeries]-N_A_edge[currentSeries] - Ktmp))
          Pprec_v[currentSeries] = Pprec_v[currentSeries] + Ptmp
        }
        for (Ktmp in K_trigg[currentSeries]:(N_B[currentSeries]-N_B_edge[currentSeries])) {
          Ptmp = choose((N_B[currentSeries]-N_B_edge[currentSeries]), Ktmp) * ((1 - ((1 - ((delT + 1)/(Tlen[currentSeries] - tau)))^N_A[currentSeries]))^Ktmp) * (((1 - ((delT + 1)/(Tlen[currentSeries] - tau)))^N_A[currentSeries])^(N_B[currentSeries]-N_B_edge[currentSeries] - Ktmp))
          Ptrigg_v[currentSeries] = Ptrigg_v[currentSeries] + Ptmp
        }
      }
      if (sym == TRUE) {
        delTsym = delT + delT + 1
        for (Ktmp in K_prec[currentSeries]:(N_A[currentSeries]-N_A_edge[currentSeries])) {
          Ptmp = choose((N_A[currentSeries]-N_A_edge[currentSeries]), Ktmp) * ((1 - ((1 - (delTsym/(Tlen[currentSeries])))^N_B[currentSeries]))^Ktmp) * (((1 - (delTsym/(Tlen[currentSeries])))^N_B[currentSeries])^(N_A[currentSeries]-N_A_edge[currentSeries] - Ktmp))
          Pprec_v[currentSeries] = Pprec_v[currentSeries] + Ptmp
        }
        for (Ktmp in K_trigg[currentSeries]:(N_B[currentSeries]-N_B_edge[currentSeries])) {
          Ptmp = choose((N_B[currentSeries]-N_B_edge[currentSeries]), Ktmp) * ((1 - ((1 - (delTsym/(Tlen[currentSeries])))^N_A[currentSeries]))^Ktmp) * (((1 - (delTsym/(Tlen[currentSeries])))^N_A[currentSeries])^(N_B[currentSeries]-N_B_edge[currentSeries] - Ktmp))
          Ptrigg_v[currentSeries] = Ptrigg_v[currentSeries] + Ptmp
        }
      }
    }
    Pprec=sum(Pprec_v*(N_A-N_A_edge)/sum(N_A-N_A_edge))
    Ptrigg=sum(Ptrigg_v*(N_B-N_B_edge)/sum(N_B-N_B_edge))
    #Pprec=1
    #Ptrigg=1
    #for(currentSeries in c(1:length(seriesListA))){
    #  Pprec=Pprec*(1-Pprec_v[currentSeries])
    #  Ptrigg=Ptrigg*(1-Ptrigg_v[currentSeries])
    #}
    #Pprec=1-Pprec
    #Ptrigg=1-Ptrigg
  }
  
  if (sigtest == "shuffle.surrogate") {
    CRprec_sur=rep(0,reps)
    CRtrigg_sur=rep(0,reps)
    for(surrogates in 1:reps){
      for(currentSeries in 1:length(seriesListA)){
        surA=rep(0,length=Tlen[currentSeries])
        surA[sample(Tlen[currentSeries],size=N_A[currentSeries])]=1
        surB=rep(0,length=Tlen[currentSeries])
        surB[sample(Tlen[currentSeries],size=N_B[currentSeries])]=1
        
        K_prec_sur=0
        for (step_a in (1+N_A_edge[currentSeries]):N_A[currentSeries]) {
          if (sym == FALSE) {
            K_prec_sur=K_prec_sur+as.numeric(any((0<=((surA[step_a] - tau) - surB)) & (((surA[step_a] - tau) - surB)<=delT)))
          }
          if (sym == TRUE) {
            K_prec_sur=K_prec_sur+as.numeric(any((-delT<=((surA[step_a] - tau) - surB)) & (((surA[step_a] - tau) - surB)<=delT)))
          }
        }
        CRprec_sur[surrogates] = CRprec_sur[surrogates] + K_prec_sur
        
        K_trigg_sur = 0
        for (step_b in 1:(N_B[currentSeries]-N_B_edge[currentSeries])) {
          if (sym == FALSE) {
            K_trigg_sur=K_trigg_sur+as.numeric(any((0<=((surA - tau) - surB[step_b])) & (((surA - tau) - surB[step_b])<=delT)))
          }
          if (sym == TRUE) {
            K_trigg_sur=K_trigg_sur+as.numeric(any((-delT<=((surA - tau) - surB[step_b])) & (((surA - tau) - surB[step_b])<=delT)))
          }
        } 
        CRtrigg_sur[surrogates] = CRtrigg_sur[surrogates] + K_trigg_sur
      }
      CRprec_sur[surrogates]=CRprec_sur[surrogates]/(sum(N_A)-sum(N_A_edge))
      CRtrigg_sur[surrogates]=CRtrigg_sur[surrogates]/(sum(N_B)-sum(N_B_edge))
    }
    Pprec=sum(CRprec<=CRprec_sur, na.rm = TRUE)/reps
    Ptrigg=sum(CRtrigg<=CRtrigg_sur, na.rm = TRUE)/reps
  }
  sig_testprec = logical
  if (Pprec >= alpha) {
    sig_testprec = TRUE
  }
  else {
    sig_testprec = FALSE
  }
  sig_testtrigg = logical
  if (Ptrigg >= alpha) {
    sig_testtrigg = TRUE
  }
  else {
    sig_testtrigg = FALSE
  }
  CA_out = list(sig_testprec, sig_testtrigg, Pprec, Ptrigg, 
                CRprec, CRtrigg)
  names(CA_out) = c("NH precursor", "NH trigger", "p-value precursor", 
                    "p-value trigger", "precursor coincidence rate", "trigger coincidence rate")
  return(CA_out)
}