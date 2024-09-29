# Autor: Leonna Szangolies, Jonatan Siegmund, Potsdam Institute for Climate Impact Research, 2020
# This function is part of the R package CoinCalc().
# CC.ceca.ts performes conditioned event coincidence analysis like described in 
# Siegmund et al. 2016: Meteorological Drivers of Extremes in Daily Stem Radius 
# Variations of Beech, Oak, and Pine in Northeastern Germany: An Event Coincidence Analysis.
# Frontiers in Plant Sciences, 7, 733.

CC.ceca.ts <- function(seriesA,seriesB,seriesC,alpha=0.05,delT=0,delT.cond=0,sym=FALSE,sym.cond=FALSE,tau=0,tau.cond=0,sigtest="poisson",reps=1000){ 
        
    # ----------------------------------------------- TEST FOR ERRORS IN FUNCTION CALL  ------------------------------------------------------------------------------  
    # ---- #1: Time series not binary 
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
        print("| Use CC.binarize() to preprocess seriesB  |")
        print("|------------------------------------------|")
        return()}
    if(length(table(seriesC))>2){  
        print("|-------------    ERROR #1    -------------|")
        print("|       Time seriesC is not binary         |")
        print("|      (or the number of events=0)!        |") 
        print("| Use CC.binarize() to preprocess seriesC  |")
        print("|------------------------------------------|")
        return()}
    # ---- #2: delT or tau negative
    if(tau<0 || delT<0){
        print("|-------------    ERROR #2    -------------|")
        print("|    The offset (tau) or delta T (delT)    |")
        print("|              is negative.                |") 
        print("|------------------------------------------|")
        return()}
    # ---- #4: lengths of series A and B differ
    if(length(seriesA)!=length(seriesB) || length(seriesA)!=length(seriesC) || length(seriesB)!=length(seriesC)){
        print("|-------------    ERROR #4    -------------|")
        print("|    lengths of series A, B and C differ   |")
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
    if(!is.vector(seriesC) && !is.matrix(seriesC)){
        print("|-------------    ERROR #5    -------------|")
        print("| series C is neither vector nor matrix    |")
        print("|------------------------------------------|")
        return()} 
    # ---- #6: seriesA or seriesB contain NAs, surrogate doesnt work
    if(any(is.na(seriesA)) && sigtest=="surrogate"){
        print("|-------------    ERROR #6    -------------|")
        print("|          seriesA contains NAs.           |")
        print("| surrogate test not allowed, use poisson  |")       
        print("|------------------------------------------|")
        return()} 
    if(any(is.na(seriesB)) && sigtest=="surrogate"){
        print("|-------------    ERROR #6    -------------|")
        print("|          seriesB contains NAs.           |")
        print("| surrogate test not allowed, use poisson  |")       
        print("|------------------------------------------|")
        return()} 
    if(any(is.na(seriesB)) && sigtest=="shuffle"){
        print("|-------------    ERROR #6    -------------|")
        print("|    seriesB or seriesA contains NAs.      |")
        print("| shuffle test not allowed, use poisson    |")       
        print("|------------------------------------------|")
        return()} 
    if(any(is.na(seriesC)) && sigtest=="surrogate"){
        print("|-------------    ERROR #6    -------------|")
        print("|          seriesC contains NAs.           |")
        print("| surrogate test not allowed, use poisson  |")       
        print("|------------------------------------------|")
        return()} 
    
    # Break function if one of the two series only holds NAs
    if(sum(is.na(seriesA))==length(seriesA) || sum(is.na(seriesB))==length(seriesB) || sum(is.na(seriesC))==length(seriesC) ){
        CA_out=list(NA,NA,NA,NA,NA,NA)
        names(CA_out)=c("NH precursor","NH trigger","p-value precursor","p-value trigger","precursor coincidence rate","trigger coincidence rate")
    return(CA_out)}
    
    
    bindata=rbind(seriesB,seriesC)
    #bindata[is.na(seriesC),1]=NAv                        !!!
    #bindata[is.na(seriesB),2]=NA                         !!!
    rownames(bindata)=c("seriesB","seriesC")

    Tlen=length(bindata[1,!is.na(bindata[1,])])
    steplen=length(seriesB)
    seriesBC=matrix(0,1,steplen)
    seriesBC[1,is.na(bindata[1,])]=NA
    N_B=as.numeric(Tlen-table(bindata[1,]==1)[1])
    N_C=as.numeric(Tlen-table(bindata[2,]==1)[1])
    N_B_edge=sum(bindata[1,][0:tau.cond+delT.cond],na.rm=TRUE)
    N_C_edge=sum(bindata[2,][(length(bindata[2,])-(tau.cond+delT.cond)+1):(length(bindata[2,])+1)],na.rm=TRUE)
    
    
    # Get simoultaneous events for Precursor Coincidence
    
    for(step in (1+N_B_edge):steplen){    
      if(is.na(bindata[1,step])){next}
        if(bindata[1,step]==1){
            start=step-tau.cond-(delT.cond)
            if(sym.cond==TRUE){
                end=step-tau.cond+(delT.cond)}
            if(sym.cond==FALSE){end=step-tau.cond}
            if(start<1 & end>=1){start=1};if(end<1){next}
            if(start>steplen){next};if(end>steplen){end=steplen}
            if( is.element(1,bindata[2,start:end])){seriesBC[1,step]=1}
        } # end if bindata[1,]==1
    } # end for step
    
    # Are there events in B found, that have been conditioned by C? If not, end function!
    if(length(table(seriesBC))<2){
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
    
    CA_out=CC.eca.ts(seriesA,seriesBC,alpha=alpha,delT=delT,sym=sym,tau=tau,sigtest=sigtest,reps=reps)
    #sigtest if(delT==0){delT=1}
    names(CA_out)=c("NH precursor","NH trigger","p-value precursor","p-value trigger","conditioned precursor coincidence rate","conditioned trigger coincidence rate")
    return(CA_out)
        
    } # end function CC.ceca.ts
