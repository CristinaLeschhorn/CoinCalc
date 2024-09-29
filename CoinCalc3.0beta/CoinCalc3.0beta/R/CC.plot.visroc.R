# Author: Leonna Szangolies, Potsdam Institute for Climate Impact Research, 2020
# This function is part of the R package CoinCalc3.0().
# CC.plot plots the Output of the VISROC anlysis.

require(plot3D)
CC.plot.visroc<-function(CL=c(),field=c(),p=c(),F1H1=c()){
  if (length(CL[1,]) != 4) {
    print("|-------------    ERROR #1    -------------|")
    print("|      Input CL has the wrong format       |")
    print("|------------------------------------------|")
    return()
  }
  if (length(field[1,]) != 3) {
    print("|-------------    ERROR #2    -------------|")
    print("|    Input field has the wrong format      |")
    print("|------------------------------------------|")
    return()
  }
  if (length(p) != 5) {
    print("|-------------    ERROR #3    -------------|")
    print("|      Input p has the wrong format        |")
    print("|------------------------------------------|")
    return()
  }
  if (length(F1H1) != 7) {
    print("|-------------    ERROR #4    -------------|")
    print("|    Input F1H1 has the wrong format       |")
    print("|------------------------------------------|")
    return()
  }
par(mar=c(4,4,4,4),mfrow=c(1,1),las=1)

if(length(field)!=0){
colfunc <- colorRampPalette(c("white","green", "blue"))
colvec<-colfunc(1000000)[round(field[,3]*2000000)+1] #colfunc(max([round(field[,3]*2000000)+1]))[round(field[,3]*2000000)+1]

plot(field[,1],field[,2],col=colvec,pch=16,ylim=c(0,1),xlim=c(0,1),main=paste("ROC Analysis"),xlab="False positive rate",ylab="True positive rate")  #"\n p at rectangle point (d18O>50%) = ",round(pval[1],3),
colkey(colfunc(1000000),c(0,0.5),add=TRUE)
}

if(length(CL)!=0){
lines(CL[,1],CL[,2],col="orange",cex=0.5)
lines(CL[,1],CL[,3],col="green",cex=0.5)
lines(CL[,1],CL[,4],col="blue",cex=0.5)
}

if(length(p)!=0){
points(p[2],p[3],pch=2)
}

legend("bottomright",legend=c("p=10%","p=5%","p=1%"),col=c("orange","green","blue"),lty=1,cex=0.7)
}