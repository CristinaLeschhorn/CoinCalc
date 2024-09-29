# Author: Leonna Szangolies, Potsdam Institute for Climate Impact Research, 2020
# This function is part of the R package CoinCalc3.0().
# CC.double.aaft performs coupled Amplitude Adjusted Fourier Transformation for two timeseries and returns surrogates

require("lomb")
CC.double.aaft<-function (data1,data2,steps=c(1:length(data1)), nsur=100) 
{
  if(length(data1)!=length(data2)){
    print("|-------------    ERROR #1    -------------|")
    print("|     Time seriesA and Time seriesB        |")
    print("|      need to have equal length           |")
    print("|------------------------------------------|")
    return()
  }
  if(!identical(steps,c(1:length(data1)))){
    data1=Re(fft(Re(lsp(data1,steps,from=0,to=(length(data1)/10),plot=FALSE)$power),inverse=TRUE))
    data2=Re(fft(Re(lsp(data2,steps,from=0,to=(length(data2)/10),plot=FALSE)$power),inverse=TRUE))
    turn=FALSE
  }else{
    turn=TRUE
  }
  n = length(data1)
  ixV1 = order(data1)
  rxV1 = rank(data1)
  oxV1 = data1[ixV1]
  surrogates1 = matrix(data = 0, n, nsur)
  ixV2 = order(data2)
  rxV2 = rank(data2)
  oxV2 = data2[ixV2]
  surrogates2 = matrix(data = 0, n, nsur)
  for (count in 1:nsur) {
    rV = rnorm(n)
    orV = rV[order(rV)]
    yV1 = orV[rxV1]
    if (n%%2 == 0) {
      n2 = n/2
    }else {
      n2 = (n - 1)/2
    }
    tmpV1 = fft(yV1, 2 * n2)
    magnV1 = abs(tmpV1)
    fiV1 = Arg(tmpV1)
    rfiV = runif(n2 - 1) * 2 * pi
    nfiV1 = c(0, rfiV, fiV1[n2 + 1], -rev(rfiV))
    tmpV1 = c(magnV1[1:(n2 + 1)], rev(magnV1[2:n2]))
    c.exp1 = cos(nfiV1) + (0+1i) * sin(nfiV1)
    tmpV1 = tmpV1 * c.exp1
    yftV1 = Re(fft(tmpV1, inverse = TRUE))
    iyftV1 = rank(yftV1)
    surrogates1[, count] = oxV1[iyftV1]
    
    yV2 = orV[rxV2]
    tmpV2 = fft(yV2, 2 * n2)
    magnV2 = abs(tmpV2)
    fiV2 = Arg(tmpV2)
    nfiV2 = c(0, rfiV, fiV2[n2 + 1], -rev(rfiV))
    tmpV2 = c(magnV2[1:(n2 + 1)], rev(magnV2[2:n2]))
    c.exp2 = cos(nfiV2) + (0+1i) * sin(nfiV2)
    tmpV2 = tmpV2 * c.exp2
    yftV2 = Re(fft(tmpV2, inverse = TRUE))
    if(turn==FALSE){
      iyftV2 = rank(yftV2)
    }
    else{
      iyftV2 = rank(-yftV2)
    }
    surrogates2[, count] = oxV2[iyftV2]
  }
  return(list("data1.series"=surrogates1,"data2.series"=surrogates2))
}