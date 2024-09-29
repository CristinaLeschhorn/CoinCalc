# Autor: Leonna Szangolies, Potsdam Institute for Climate Impact Research, 2020
# This function is part of the R package CoinCalc3.0().
# ecr.ts performes event coincidence analysis vor varying tau using time series and some parameters. The significance level using two different significance tests can be computed.

CC.ecr.ts<-function(seriesA,seriesB=seriesA,alpha=0.05,delT=0,sym=FALSE,offset=30,siglevel=TRUE,sigtest="shuffle.surrogate",reps=1000){
  ecrprec<-vector(length=offset)
  ecrtrig<-vector(length=offset)
  ecrprecsig<-vector(length=offset)
  ecrtrigsig<-vector(length=offset)
  for(i in 0:(offset-1)){
    tau=i
    # ----------------------------------------------- TEST FOR ERRORS IN FUNCTION CALL  ------------------------------------------------------------------------------  
    # ---- #1: Time series not binary 
    if(offset>length(seriesA)){
      print("|-------------    ERROR #1    -------------|")
      print("|    Offset cant be longer than Series     |")
      print("|------------------------------------------|")
      return()}
    if(length(table(seriesA))>2){
      print("|-------------    ERROR #1    -------------|")
      print("|       Time seriesA is not binary         |")
      print("|      (or the number of events=0)!        |") 
      print("| Use CC.binarize() to preprocess seriesA  |")
      print("|------------------------------------------|")
      return()}
    if(length(table(seriesB))>2){  
      print("|-------------    ERROR #1    -------------|")
      print("|       Time seriesB is not binary         |")
      print("|      (or the number of events=0)!        |") 
      print("| Use CC.binarize() to preprocess seriesA  |")
      print("|------------------------------------------|")
      return()}
    # ---- #2: delT or tau negative
    if(delT<0){
      print("|-------------    ERROR #2    -------------|")
      print("|            The delta T (delT)            |")
      print("|              is negative.                |") 
      print("|------------------------------------------|")
      return()}
    # ---- #4: lengths of series A and B differ
    if(length(seriesA)!=length(seriesB)){
      print("|-------------    ERROR #4    -------------|")
      print("|    lengths of series A and B differ      |")
      print("|------------------------------------------|")
      return()} 
    # ---- #5: seriesA or seriesB is not a vector
    if(!is.vector(seriesA) && !is.matrix(seriesA)){
      print("|-------------    ERROR #5    -------------|")
      print("| series A is neither vector nor matrix    |")
      print("|------------------------------------------|")
      return()} 
    # ---- #5: seriesA or seriesB is not a vector
    if(!is.vector(seriesB) && !is.matrix(seriesB)){
      print("|-------------    ERROR #5    -------------|")
      print("| series B is neither vector nor matrix    |")
      print("|------------------------------------------|")
      return()} 
    
    # Break function if one of the two series only holds NAs
    if(sum(is.na(seriesA))==length(seriesA) || sum(is.na(seriesB))==length(seriesB)){
      CA_out=list(NA,NA,NA,NA,NA,NA)
      names(CA_out)=c("NH precursor","NH trigger","p-value precursor","p-value trigger","precursor coincidence rate","trigger coincidence rate")
      return(CA_out)
    }
    
    # seriesA or seriesB contain NAs, needs special handling
    if ((any(is.na(seriesA)) || any(is.na(seriesB))) && (sigtest == "wt.surrogate" || sigtest =="shuffle.surrogate" || sigtest=="periodic.surrogate")) {
      resultNA=handle.NA(seriesA,seriesB,delT,tau)
      seriesA=resultNA[[1]]
      seriesB=resultNA[[2]]
    }
    
    # write the two binarized time series to bindata
    seriesA[is.na(seriesB)]=NA
    seriesB[is.na(seriesA)]=NA
    bindata=matrix(NA,2,length(seriesA))
    bindata[1,]=seriesA
    bindata[2,]=seriesB
    rownames(bindata)=c("seriesA","seriesB")
    
    
    # ----------------------------------------------- COINCIDENCE ANALYSIS part I ------------------------------------------------------------------  
    # Get Coincidences 
    #  Count number of simultaneous events (Coincicences, K), number of events per series (N_A, N_N) and the maximum number of coincidences (Kmax)
    Tlen=length(bindata[1,!is.na(bindata[1,])])
    steplen=length(seriesA)
    N_A=as.numeric(Tlen-table(bindata[1,]==1)[1])
    N_B=as.numeric(Tlen-table(bindata[2,]==1)[1])
    N_A_edge=sum(bindata[1,][0:tau+delT],na.rm=TRUE)
    N_B_edge=sum(bindata[2,][(length(bindata[2,])-(tau+delT)+1):(length(bindata[2,])+1)],na.rm=TRUE)
    
    # Get Precursor Coincidence Rate
    K_prec=0
    for(step in (1+N_A_edge):steplen){    
      if(is.na(bindata[1,step])){next}
      if(bindata[1,step]==1){
        start=step-tau-(delT)
        if(sym==TRUE){
          end=step-tau+(delT)}
        if(sym==FALSE){end=step-tau}
        if(start<1 & end>=1){start=1};if(end<1){next}
        if(start>steplen){next};if(end>steplen){end=steplen}
        if( is.element(1,bindata[2,start:end])){K_prec=K_prec+1}
      } # end if bindata[1,]==1
    } # end for step
    CRprec=K_prec/(N_A-N_A_edge)
    
    # Get Trigger Coincidence Rate 
    K_trigg=0
    for(step in 1:(steplen-N_B_edge)){    if(is.na(bindata[2,step])){next}
      if(bindata[2,step]==1){
        if(sym==TRUE){
          start=step+tau-(delT)}
        if(sym==FALSE){
          start=step+tau}
        end=step+tau+(delT)
        if(start<1 & end>=1){start=1};if(end<1){next}
        if(start>steplen){next};if(end>steplen){end=steplen}
        if( is.element(1,bindata[1,start:end] )){K_trigg=K_trigg+1}
      } # end if bindata[2,]==1
    } # end for step 
    CRtrigg=K_trigg/(N_B-N_B_edge)
    
    ecrprec[i+1]<-CRprec
    ecrtrig[i+1]<-CRtrigg
    
    # ----------------------------------------------- Significance Level  -------------------------------------------------------------------------------------------    
    if(siglevel){
    if(sigtest=="wt.surrogate"){
      # create 'reps' (standard = 1000) surrogate sime series for seriesA and seriesB, having the same average waiting times between two events. Then performe
      # coincidence analysis as given above for all 1000 pairs of event series. Finally get the percentile of the ampirically found coincidence rate using the 
      # ecdf of the 1000 surrogatecoincidence rates.
      if(delT==0){delT=1}
      seriesA=CC.ts2es(seriesA) # transform time series to event series
      seriesB=CC.ts2es(seriesB)
      span=(max(seriesA$span[1],seriesB$span[1]):min(seriesA$span[2],seriesB$span[2]))
      Tlen=length(span)
      # calculate waiting times of seriesA
      seriesA_sort=sort(seriesA$es)
      wtA=rep(0,N_A)
      wtA[1]=seriesA_sort[1]-seriesA$span[1]
      for(i in 2:N_A){wtA[i]=seriesA_sort[i]-seriesA_sort[i-1]}
      # calculate waiting times of seriesB
      seriesB_sort=sort(seriesB$es)
      wtB=rep(0,N_B)
      wtB[1]=seriesB_sort[1]-seriesB$span[1]
      for(i in 2:N_B){wtB[i]=seriesB_sort[i]-seriesB_sort[i-1]}
      # create reps times two surrogate event series and perform coincidence analysis
      surdist=matrix(NA,2,reps)
      for(surno in 1:reps){
        # create surrogate for seriesA
        surA=sample(wtA,size=1)
        for(i in 2:length(seriesA$es)){
          tmp=surA[i-1]+sample(wtA,size=1)
          if(tmp>span[Tlen]){break}else{surA[i]=tmp}}
        # create surrogate for seriesB
        surB=sample(wtB,size=1)
        for(i in 2:length(seriesB$es)){
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
      # create reps times two surrogate event series and perform coincidence analysis
      surdist = matrix(NA, 2, reps)
      span = seq(1:steplen)
      for (surno in 1:reps) {
        surA = sample(span, size = N_A)
        surB = sample(span, size = N_B)
        K_prec_sur = 0
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
        surdist[2, surno] = K_trigg_sur/(N_B-N_B_edge)
      }
      Pprecsig = quantile(surdist[1, ],0.95)
      Ptriggsig = quantile(surdist[2, ],0.95)      
    } # end if sigtest shuffle 
    
    }
  }
  CR_out=list(c(0:(offset-1)),ecrprec,ecrtrig,ecrprecsig,ecrtrigsig)
  names(CR_out)=c("tau","precursor rates","trigger rates","precursor significance level","trigger significance level")
  return(CR_out)
}