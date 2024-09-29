# Author: Leonna Szangolies, Potsdam Institute for Climate Impact Research, 2020
# This function is part of the R package CoinCalc3.0().
# CC.eca.aggregated.ts performes aggregated event coincidence analysis including two different significance tests using several pairs of time series and some parameters.

CC.eca.aggregated.ts<-function (seriesList, alpha = 0.05, delT = 0, sym = FALSE, tau = 0, sigtest = "poisson", reps = 1000) 
{
CRprec=0
CRtrigg=0
N_A=vector(length=length(seriesList))
N_B=vector(length=length(seriesList))
N_A_edge=vector(length=length(seriesList))
N_B_edge=vector(length=length(seriesList))
Tlen=vector(length=length(seriesList))
K_prec=vector(length=length(seriesList))
K_trigg=vector(length=length(seriesList))
for(currentSeries in 1:length(seriesList)){
  seriesA=seriesList[[currentSeries]][,1]
  seriesB=seriesList[[currentSeries]][,2]
  if (length(table(seriesA)) > 2) {
    print("|-------------    ERROR #1    -------------|")
    print("|  One of the Time seriesA is not binary   |")
    print("|      (or the number of events=0)!        |")
    print("| Use CC.binarize() to preprocess seriesA  |")
    print("|------------------------------------------|")
    return()
  }
  if (length(table(seriesB)) > 2) {
    print("|-------------    ERROR #1    -------------|")
    print("|  One of the Time seriesB is not binary   |")
    print("|      (or the number of events=0)!        |")
    print("| Use CC.binarize() to preprocess seriesB  |")
    print("|------------------------------------------|")
    return()
  }
  if (tau < 0 || delT < 0) {
    print("|-------------    ERROR #2    -------------|")
    print("|    The offset (tau) or delta T (delT)    |")
    print("|              is negative.                |")
    print("|------------------------------------------|")
    return()
  }
  if (length(seriesA) != length(seriesB)) {
    print("|-------------    ERROR #4    -------------|")
    print("|    lengths of series A and B differ      |")
    print("|------------------------------------------|")
    return()
  }
  if (!is.vector(seriesA) && !is.matrix(seriesA)) {
    print("|-------------    ERROR #5    -------------|")
    print("| series A is neither vector nor matrix    |")
    print("|------------------------------------------|")
    return()
  }
  if (!is.vector(seriesB) && !is.matrix(seriesB)) {
    print("|-------------    ERROR #5    -------------|")
    print("| series B is neither vector nor matrix    |")
    print("|------------------------------------------|")
    return()
  }
  if (any(is.na(seriesA)) && sigtest == "surrogate") {
    print("|-------------    ERROR #6    -------------|")
    print("|          seriesA contains NAs.           |")
    print("| surrogate test not allowed, use poisson  |")
    print("|------------------------------------------|")
    return()
  }
  if (any(is.na(seriesB)) && sigtest == "surrogate") {
    print("|-------------    ERROR #6    -------------|")
    print("|          seriesB contains NAs.           |")
    print("| surrogate test not allowed, use poisson  |")
    print("|------------------------------------------|")
    return()
  }
  if (any(is.na(seriesB)) && sigtest == "shuffle") {
    print("|-------------    ERROR #6    -------------|")
    print("|    seriesB or seriesA contains NAs.      |")
    print("| shuffle test not allowed, use poisson    |")
    print("|------------------------------------------|")
    return()
  }
  if (sum(is.na(seriesA)) == length(seriesA) || sum(is.na(seriesB)) == 
      length(seriesB)) {
    CA_out = list(NA, NA, NA, NA, NA, NA)
    names(CA_out) = c("NH precursor", "NH trigger", "p-value precursor", 
                      "p-value trigger", "precursor coincidence rate", 
                      "trigger coincidence rate")
    return(CA_out)
  }
  seriesA[is.na(seriesB)] = NA
  seriesB[is.na(seriesA)] = NA
  bindata = matrix(NA, 2, length(seriesA))
  bindata[1, ] = seriesA
  bindata[2, ] = seriesB
  rownames(bindata) = c("seriesA", "seriesB")
  Tlen[currentSeries] = length(bindata[1, !is.na(bindata[1, ])])
  steplen = length(seriesA)
  N_A[currentSeries] = as.numeric(Tlen[currentSeries] - table(bindata[1, ] == 1)[1])
  N_B[currentSeries] = as.numeric(Tlen[currentSeries] - table(bindata[2, ] == 1)[1])
  N_A_edge[currentSeries]=sum(bindata[1,][0:tau+delT],na.rm=TRUE)
  N_B_edge[currentSeries]=sum(bindata[2,][(length(bindata[2,])-(tau+delT)+1):(length(bindata[2,])+1)],na.rm=TRUE)
  K_prec[currentSeries] = 0
  for (step in (1+N_A_edge[currentSeries]):steplen) {
    if (is.na(bindata[1, step])) {
      next
    }
    if (bindata[1, step] == 1) {
      start = step - tau - (delT)
      if (sym == TRUE) {
        end = step - tau + (delT)
      }
      if (sym == FALSE) {
        end = step - tau
      }
      if (start < 1 & end >= 1) {
        start = 1
      }
      if (start < 1 & end < 1) {
        next
      }
      if (start > steplen) {
        next
      }
      if (end > steplen) {
        end = steplen
      }
      if (is.element(1, bindata[2, start:end])) {
        K_prec[currentSeries] = K_prec[currentSeries] + 1
      }
    }
  }
  CRprec = CRprec + K_prec[currentSeries]
  K_trigg[currentSeries] = 0
  for (step in 1:(steplen-N_B_edge[currentSeries])) {
    if (is.na(bindata[2, step])) {
      next
    }
    if (bindata[2, step] == 1) {
      if (sym == TRUE) {
        start = step + tau - (delT)
      }
      if (sym == FALSE) {
        start = step + tau
      }
      end = step + tau + (delT)
      if (start < 1 & end >= 1) {
        start = 1
      }
      if (start < 1 & end < 1) {
        next
      }
      if (start > steplen) {
        next
      }
      if (end > steplen) {
        end = steplen
      }
      if (is.element(1, bindata[1, start:end])) {
        K_trigg[currentSeries] = K_trigg[currentSeries] + 1
      }
    }
  }
  CRtrigg = CRtrigg + K_trigg[currentSeries]
}
CRprec=CRprec/(sum(N_A)-sum(N_A_edge))
CRtrigg=CRtrigg/(sum(N_B)-sum(N_B_edge))

if (sigtest == "poisson") {
  Pprec_v = rep(0,length=length(seriesList))
  Ptrigg_v = rep(0,length=length(seriesList))
  for(currentSeries in 1:length(seriesList)){
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
}

if (sigtest == "shuffle.surrogate") {
CRprec_sur=rep(0,reps)
CRtrigg_sur=rep(0,reps)
for(surrogates in 1:reps){
for(currentSeries in 1:length(seriesList)){
    surA=sample(Tlen[currentSeries],size=N_A[currentSeries])
    surB=sample(Tlen[currentSeries],size=N_B[currentSeries])
    K_prec_sur = 0
    if(N_A[currentSeries]>0){
    for (step_a in (1+N_A_edge[currentSeries]):N_A[currentSeries]) {
      if (sym == FALSE) {
        K_prec_sur=K_prec_sur+as.numeric(any((0<=((surA[step_a] - tau) - surB)) & (((surA[step_a] - tau) - surB)<=delT)))
      }
      if (sym == TRUE) {
        K_prec_sur=K_prec_sur+as.numeric(any((-delT<=((surA[step_a] - tau) - surB)) & (((surA[step_a] - tau) - surB)<=delT)))
      }
    }
    }
    if(is.na(K_prec_sur)){print(currentSeries)}
    CRprec_sur[surrogates] = CRprec_sur[surrogates] + K_prec_sur
    K_trigg_sur = 0
    if(N_B[currentSeries]>0){
    for (step_b in 1:(N_B[currentSeries]-N_B_edge[currentSeries])) {
      if (sym == FALSE) {
        K_trigg_sur=K_trigg_sur+as.numeric(any((0<=((surA - tau) - surB[step_b])) & (((surA - tau) - surB[step_b])<=delT)))
      }
      if (sym == TRUE) {
        K_trigg_sur=K_trigg_sur+as.numeric(any((-delT<=((surA - tau) - surB[step_b])) & (((surA - tau) - surB[step_b])<=delT)))
      }
    }
    }
    CRtrigg_sur[surrogates] = CRtrigg_sur[surrogates] + K_trigg_sur
  }
  
  CRprec_sur[surrogates]=CRprec_sur[surrogates]/(sum(N_A)-sum(N_A_edge))
  CRtrigg_sur[surrogates]=CRtrigg_sur[surrogates]/(sum(N_B)-sum(N_B_edge))
}
Pprec = 1 - ecdf(CRprec_sur)(CRprec)
Ptrigg = 1 - ecdf(CRtrigg_sur)(CRtrigg) 
#Pprec=sum(CRprec<=CRprec_sur, na.rm = TRUE)/reps
#Ptrigg=sum(CRtrigg<=CRtrigg_sur, na.rm = TRUE)/reps
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