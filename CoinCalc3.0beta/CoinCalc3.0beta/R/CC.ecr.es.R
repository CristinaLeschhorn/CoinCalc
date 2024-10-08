# Autor: Leonna Szangolies, Potsdam Institute for Climate Impact Research, 2020
# This function is part of the R package CoinCalc3.0().
# ecr.es performes event coincidence analysis vor varying tau using event sequences and some parameters. The significance level using two different significance tests can be computed.

CC.ecr.es<-function(seriesA,spanA,seriesB=seriesA,spanB=spanA,alpha=0.05,delT=0,sym=FALSE,offset=30,siglevel=TRUE,sigtest="shuffle.surrogate",reps=1000){
  ecrprec<-vector(length=offset)
  ecrtrig<-vector(length=offset)
  ecrprecsig<-vector(length=offset)
  ecrtrigsig<-vector(length=offset)
  for(i in 0:(offset-1)){
    tau=i
    # ----------------------------------------------- TEST FOR ERRORS IN FUNCTION CALL  ------------------------------------------------------------------------------  
    # ---- #1: delT or tau negative
    if(tau<0 || delT<0){
      print("|-------------    ERROR #1    -------------|")
      print("|    The offset (tau) or delta T (delT)    |")
      print("|              is negative.                |") 
      print("|------------------------------------------|")
      return()}
    # ---- #2: seriesA is not a vector
    if(!is.vector(seriesA)){
      print("|-------------    ERROR #2    -------------|")
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
    
    # ----------------------------------------------- EVENT COINCIDENCE ANALYSIS  ------------------------------------------------------------------------------------  
    # Get Coincidences 
    N_A=length(seriesA) # number of events in seriesA
    N_B=length(seriesB) # number of ecents in seriesB
    span=(max(spanA[1],spanB[1]):min(spanA[2],spanB[2])) # getting common time span of seriesA and B. Will be used for coincidence analysis.
    N_A_edge=sum(is.element(seriesA,c((spanA[1]-1):(spanA[1]+tau+delT-1))))
    N_B_edge=sum(is.element(seriesB,c((spanB[2]-(tau+delT)+1):(spanB[2]+1))))
    
    # Get Precursor Coincidence Rate
    K_prec=0
    for(step_a in (1+N_A_edge):N_A){
      if( !is.element(seriesA[step_a],span) ){next} 
      if (sym == FALSE) {
        K_prec=K_prec+as.numeric(any((0<=((seriesA[step_a] - tau) - seriesB[is.element(seriesB,span)])) & (((seriesA[step_a] - tau) - seriesB[is.element(seriesB,span)])<=delT)))
      }
      if (sym == TRUE) {
        K_prec=K_prec+as.numeric(any((-delT<=((seriesA[step_a] - tau) - seriesB[is.element(seriesB,span)])) & (((seriesA[step_a] - tau) - seriesB[is.element(seriesB,span)])<=delT)))
      } 
    } # end for step_a
    CRprec=K_prec/(N_A-N_A_edge)
    
    # Get Trigger Coincidence Rate 
    K_trigg=0
    for(step_b in 1:(N_B-N_B_edge)){
      if( !is.element(seriesB[step_b],span) ){next}
      if (sym == FALSE) {
        K_trigg=K_trigg+as.numeric(any((0<=((seriesA[is.element(seriesA,span)] - tau) - seriesB[step_b])) & (((seriesA[is.element(seriesA,span)] - tau) - seriesB[step_b])<=delT)))
      }
      if (sym == TRUE) {
        K_trigg=K_trigg+as.numeric(any((-delT<=((seriesA[is.element(seriesA,span)] - tau) - seriesB[step_b])) & (((seriesA[is.element(seriesA,span)] - tau) - seriesB[step_b])<=delT)))
      }
    } # end for step_b
    CRtrigg=K_trigg/(N_B-N_B_edge)      
    ecrprec[i+1]<-CRprec
    ecrtrig[i+1]<-CRtrigg
    
    # ----------------------------------------------- Significance Test ------------------------------------------------------------------------------------------------
    
    if(siglevel){
    if(sigtest=="wt.surrogate"){
      # create 'reps' (standard = 1000) surrogate sime series for seriesA and seriesB, having the same average waiting times between two events. Then performe
      # coincidence analysis as given above for all 1000 pairs of event series. Finally get the percentile of the ampirically found coincidence rate using the 
      # ecdf of the 1000 surrogatecoincidence rates.
      Tlen=length(span)
      # calculate waiting times of seriesA
      seriesA_sort=sort(seriesA)
      wtA=rep(0,N_A)
      wtA[1]=seriesA_sort[1]-spanA[1]
      for(i in 2:N_A){wtA[i]=seriesA_sort[i]-seriesA_sort[i-1]}
      # calculate waiting times of seriesB
      seriesB_sort=sort(seriesB)
      wtB=rep(0,N_B)
      wtB[1]=seriesB_sort[1]-spanB[1]
      for(i in 2:N_B){wtB[i]=seriesB_sort[i]-seriesB_sort[i-1]}
      # create reps times two surrogate event series and perform coincidence analysis
      surdist=matrix(NA,2,reps)
      for(surno in 1:reps){
        # create surrogate for seriesA
        surA=sample(wtA,size=1)
        for(i in 2:length(seriesA)){
          tmp=surA[i-1]+sample(wtA,size=1)
          if(tmp>span[Tlen]){break}else{surA[i]=tmp}}
        # create surrogate for seriesB
        surB=sample(wtB,size=1)
        for(i in 2:length(seriesB)){
          tmp=surB[i-1]+sample(wtB,size=1)
          if(tmp>span[Tlen]){break}else{surB[i]=tmp}}
        # perform event coincidence analysis between the two surrogate event time series as given above:
        # Get surrogate Precursor Coincidence Rate
        K_prec_sur=0
        for(step_a in (1+N_A_edge):N_A){
          if( !is.element(surA[step_a],span) ){next}
          if (sym == FALSE) {
            K_prec_sur=K_prec_sur+as.numeric(any((0<=((surA[step_a] - tau) - surB[is.element(surB,span)])) & (((surA[step_a] - tau) - surB[is.element(surB,span)])<=delT)))
          }
          if (sym == TRUE) {
            K_prec_sur=K_prec_sur+as.numeric(any((-delT<=((surA[step_a] - tau) - surB[is.element(surB,span)])) & (((surA[step_a] - tau) - surB[is.element(surB,span)])<=delT)))
          }
        } # end for step_a
        surdist[1,surno]=K_prec_sur/(N_A-N_A_edge)
        # Get surrogate Trigger Coincidence Rate 
        K_trigg_sur=0
        for(step_b in 1:(N_B-N_B_edge)){
          if( !is.element(surB[step_b],span) ){next}
          if (sym == FALSE) {
            K_trigg_sur=K_trigg_sur+as.numeric(any((0<=((surA[is.element(surA,span)] - tau) - surB[step_b])) & (((surA[is.element(surA,span)] - tau) - surB[step_b])<=delT)))
          }
          if (sym == TRUE) {
            K_trigg_sur=K_trigg_sur+as.numeric(any((-delT<=((surA[is.element(surA,span)] - tau) - surB[step_b])) & (((surA[is.element(surA,span)] - tau) - surB[step_b])<=delT)))
          }
        } # end for step_b
        surdist[2,surno]=K_trigg_sur/(N_B-N_B_edge)  
      }
      Pprecsig = quantile(surdist[1, ],0.95)
      Ptriggsig = quantile(surdist[2, ],0.95) 
    } # end if sigtest surrogate
    
    
    if(sigtest=="shuffle.surrogate"){
      # create 'reps' (standard = 1000) shuffled time series for seriesA and seriesB, having the same number of events. Then performe
      # coincidence analysis as given above for all 1000 pairs of event series. Finally get the percentile of the empirically found coincidence rate using the 
      # ecdf of the 1000 surrogatecoincidence rates.
      Tlen=length(span)
      span=seq(1:Tlen)
      # create reps times two surrogate event series and perform coincidence analysis
      surdist=matrix(NA,2,reps)
      for(surno in 1:reps){
        # create shuffle for seriesA
        surA=sample(span,size=N_A)
        # create shuffle for seriesB
        surB=sample(span,size=N_B)
        # perform event coincidence analysis between the two surrogate event time series as given above:
        # Get surrogate Precursor Coincidence Rate
        K_prec_sur=0
        for (step_a in (1+N_A_edge):N_A) {
          if (sym == FALSE) {
            K_prec_sur=K_prec_sur+as.numeric(any((0<=((surA[step_a] - tau) - surB)) & (((surA[step_a] - tau) - surB)<=delT)))
          }
          if (sym == TRUE) {
            K_prec_sur=K_prec_sur+as.numeric(any((-delT<=((surA[step_a] - tau) - surB)) & (((surA[step_a] - tau) - surB)<=delT)))
          }
        }
        surdist[1, surno] = K_prec_sur/(N_A-N_A_edge)
        K_trigg_sur = 0
        for (step_b in 1:(N_B-N_B_edge)) {
          if (sym == FALSE) {
            K_trigg_sur=K_trigg_sur+as.numeric(any((0<=((surA - tau) - surB[step_b])) & (((surA - tau) - surB[step_b])<=delT)))
          }
          if (sym == TRUE) {
            K_trigg_sur=K_trigg_sur+as.numeric(any((-delT<=((surA - tau) - surB[step_b])) & (((surA - tau) - surB[step_b])<=delT)))
          }
        } 
        surdist[2,surno]=K_trigg_sur/(N_B-N_B_edge)  
      }
      Pprecsig = quantile(surdist[1, ],1-alpha)
      Ptriggsig = quantile(surdist[2, ],1-alpha)       
    } # end if sigtest shuffle 
    ecrprecsig[i+1]<-Pprecsig
    ecrtrigsig[i+1]<-Ptriggsig
  }
    
  }
  CR_out=list(c(0:(offset-1)),ecrprec,ecrtrig,ecrprecsig,ecrtrigsig)
  names(CR_out)=c("tau","precursor rates","trigger rates","precursor significance level","trigger significance level")
  return(CR_out)
}