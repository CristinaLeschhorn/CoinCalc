# Autor: Leonna Szangolies, Potsdam Institute for Climate Impact Research, 2020
# This function is part of the R package CoinCalc3.0().
# CC.ceca.es performes conditioned event coincidence analysis for event sequences

CC.ceca.es <- function(seriesA,seriesB,seriesC,spanA,spanB,spanC,alpha=0.05,delT=0,delT.cond=0,sym=FALSE,sym.cond=FALSE,tau=0,tau.cond=0,sigtest="poisson",reps=1000){ 
  
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
  # ---- #3: seriesB is not a vector
  if(!is.vector(seriesC)){
    print("|-------------    ERROR #3    -------------|")
    print("|         series C is not a vector         |")
    print("|------------------------------------------|")
    return()}
  # ---- #4: time spans do not match
  if(spanA[1]>spanB[2] || spanB[1]>spanA[2] || spanA[1]>spanC[2] || spanC[1]>spanA[2] || spanB[1]>spanC[2] || spanC[1]>spanB[2]){
    print("|-------------    ERROR #4    -------------|")
    print("|         no common time span found        |")
    print("|       for spanA and spanB and spanC      |")       
    print("|------------------------------------------|")
    return()} 
  # ---- #6: list format
  if(is.list(seriesA)|| is.list(seriesB) || is.list(seriesC)){
    print("|-------------    ERROR #6    -------------|")
    print("|   seriesA or seriesB or seriesC seems    |")
    print("|   to be a list but has to be a vector    |")       
    print("|------------------------------------------|")
    return()} 
  
  N_B=length(seriesB) # number of events in seriesA
  N_C=length(seriesC) # number of ecents in seriesB
  seriesBC=c()
  span=(max(spanB[1],spanC[1]):min(spanB[2],spanC[2])) # getting common time span of seriesA and B. Will be used for coincidence analysis.
  spanBC=c(max(spanB[1],spanC[1]),min(spanB[2],spanC[2]))
  N_B_edge=sum(is.element(seriesB,c((spanB[1]-1):(spanB[1]+tau.cond+delT.cond-1))))
  N_C_edge=sum(is.element(seriesC,c((spanC[2]-(tau.cond+delT.cond)+1):(spanC[2]+1))))
  
  # Get Precursor Coincidence Rate
  for(step_a in (1+N_B_edge):N_B){
    if( !is.element(seriesB[step_a],span) ){next} 
    if (sym == FALSE) {
      if(any((0<=((seriesB[step_a] - tau.cond) - seriesB[is.element(seriesC,span)])) & (((seriesB[step_a] - tau.cond) - seriesB[is.element(seriesC,span)])<=delT.cond))){
        seriesBC<-c(seriesBC,seriesB[step_a])
      }
    }
    if (sym == TRUE) {
      if(any((-delT.cond<=((seriesB[step_a] - tau.cond) - seriesC[is.element(seriesC,span)])) & (((seriesB[step_a] - tau.cond) - seriesC[is.element(seriesC,span)])<=delT.cond))){
        seriesBC<-c(seriesBC,seriesB[step_a])
      }
    } 
  } 
  
  # Are there events in B found, that have been conditioned by C? If not, end function!
  if(length(seriesBC)<1){
    sig_testprec=TRUE
    sig_testtrigg=TRUE
    Pprec=1
    Ptrigg=1
    CRprec=0
    CRtrigg=0
    CA_out=list(sig_testprec,sig_testtrigg,Pprec,Ptrigg,CRprec,CRtrigg)
    print(paste0("No C-conditioned B-events found using delT.cond = ",delT.cond,"and tau.cond = ",tau.cond))
    names(CA_out)=c("NH precursor","NH trigger","p-value precursor","p-value trigger","conditioned precursor coincidence rate","conditioned trigger coincidence rate")
    return(CA_out)}
  
  CA_out=CC.eca.es(seriesA,seriesBC,spanA,spanBC,alpha=alpha,delT=delT,sym=sym,tau=tau,sigtest=sigtest,reps=reps)
  #sigtest if(delT==0){delT=1}
  names(CA_out)=c("NH precursor","NH trigger","p-value precursor","p-value trigger","conditioned precursor coincidence rate","conditioned trigger coincidence rate")
  return(CA_out)
  
} # end function CC.ceca.es
