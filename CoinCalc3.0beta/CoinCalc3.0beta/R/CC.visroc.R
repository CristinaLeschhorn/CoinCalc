# Author: Leonna Szangolies, Potsdam Institute for Climate Impact Research, 2020
# This function is part of the R package CoinCalc3.0().
# CC.visroc computation p-values for a Receiver Operating Characteristics (ROC) analysis

require(pracma)
CC.visroc<-function(p,q,N,uauc,f1,h1, outCL=TRUE, outp=TRUE, F1H1=TRUE, outfield=TRUE){
work<-vector(length=900)
freq<-vector(length=900)
auccr<-vector(length=3)
auc1<-vector(length=900)
pdf<-vector(length=900)
cdf<-vector(length=900)
kcr<-vector(length=3)
val<-vector(length=3)
#Check the validity of input parameters
if (p<=0) p=1
if (q<=0) q=1
if (N<1) N=1
if ((uauc<0)||(uauc>1)) uauc=0.5
if ((p==1)||(q==1)||(N==1)||(uauc<=0.5)||(h1<=f1)) ifault=1
if (h1<=f1){
  mauc=h1
  h1=f1
  f1=mauc
}
#Check on the use of the Gaussian approximation
pp=as.double(p)
qq=as.double(q)
case=0
if ((max(p,q)>=30)&&(p+q>=40)) case=1
if (case==0){
#Use of AS62 to estimate the distribution of AUC
piqo=min(p,q)
pq1=p*q+1
paxo=max(p,q)
q1=paxo+1
for(i in 1:q1){
  freq[i]=1
}
q1=q1+1
for(i in q1:pq1){
  freq[i]=0
}
work[1]=0
iq=paxo
for(i in 2:piqo){
  work[i]=0
  iq=iq+paxo
  q1=iq+2
  l=1+iq/2
  w=i
  for(j in 1:l){
    w=w+1
    q1=q1-1
    sumAS62=freq[j]+work[j]
    freq[j]=sumAS62
    work[w]=sumAS62-freq[q1]
    freq[q1]=sumAS62
  }
}
sumAS62=0
for(i in 1:p*q+1){
  sumAS62=sumAS62+freq[i]
  auc1[i]=1-(i-1)/(p*q)
}
cdf[1]=1
for(i in 1:(p*q+1)){
  pdf[i]=freq[i]/sumAS62
  auc1[i]=1-(i-1)/(p*q)
  cdf[i+1]=cdf[i]-pdf[i]
  if (cdf[i]>=0.9) auccr[1]=auc1[i]
  if (cdf[i]>=0.95) auccr[2]=auc1[i]
  if (cdf[i]>=0.99) auccr[3]=auc1[i]
}
}
if (case==1){
#Use of Gaussian approximation
s=sqrt(1/qq+1/pp+1/(pp*qq))/sqrt(12)
auccr[1]=0.5+1.2816*s
auccr[2]=0.5+1.6449*s
auccr[3]=0.5+2.3264*s
}
final<-list()
#Calculate k-values for three AUC confidence levels
if(outCL){
for(l in 1:3){
  y0=auccr[l]
  know=0
  kold=1
  while(abs(know-kold)>0.001*kold){
  k=know
  x1=(pp*qq-sqrt(p^2*q^2+qq*(pp+k)*(k*k+qq*(k-pp))))/(2*qq*(pp+k))+0.5
  r=sqrt(k+4*qq*(x1-x1*x1))
  r0=sqrt(k)
  auc=1-x1/2+(qq/(qq+k))*(x1-1)*x1/2+sqrt(k*(k+qq+pp)/pp)/(2*(k+qq))*((k+qq)*atan(sqrt(qq)*(2*x1-1)/r)/(4*sqrt(qq))+(2*x1-1)*r/4+(k+qq)*atan(sqrt(qq)/r0)/(4*sqrt(qq))+r0/4)
  kold=know
  x1old=x1
  aucold=auc
  k=kold+0.001
  x1=(pp*qq-sqrt(pp**2*qq**2+qq*(pp+k)*(k*k+qq*(k-pp))))/(2*qq*(pp+k))+0.5
  r=sqrt(k+4*qq*(x1-x1*x1))
  r0=sqrt(k)
  auc=1-x1/2+(qq/(qq+k))*(x1-1)*x1/2+sqrt(k*(k+qq+pp)/pp)/(2*(k+qq))*((k+qq)*atan(sqrt(qq)*(2*x1-1)/r)/(4*sqrt(qq))+(2*x1-1)*r/4+(k+qq)*atan(sqrt(qq)/r0)/(4*sqrt(qq))+r0/4)
  dauc=(auc-aucold)/0.001
  know=kold-(aucold-y0)/dauc
  }
  kcr[l]=know
}

CL<-matrix(ncol=4)
for(ix in 0:1000){
  x=ix/1000
  for(l in 1:3){
    k=kcr[l]
    valy=0.5+qq/(qq+k)*(x-0.5)+0.5/(k+qq)*sqrt(k*(k+qq+pp)*(k+4*qq*(x-x*x))/pp)
    val[l]=valy
  }
CL<-rbind(CL,round(c(x,val[1],val[2],val[3]),8))
}
CL<-CL[-1,]
final[["CL"]]<-CL
#write.table(CL,paste(outputdirectory,'/outCL.dat',sep=""),quote=FALSE,sep=" ",row.names=FALSE,col.names=FALSE)
}

#Calculation of the p field on the ROC diagram
if(outfield){
field<-matrix(ncol=3)
x=0
y=0
for(ix in 0:N){
  for(iy in 0:N){
    x=ix/N
    y=iy/N
    xv=x-0.5
    yv=y-0.5
    k=2*(pp*yv**2+qq*xv**2)-0.5*(pp+qq)+sqrt((2*(pp*yv**2+qq*xv**2)-0.5*(pp+qq))**2+4*pp*qq*(xv-yv)**2)
    x1=(pp*qq-sqrt(pp**2*qq**2+qq*(pp+k)*(k*k+qq*(k-pp))))/(2*qq*(pp+k))+0.5
    r=sqrt(k+4*qq*(x1-x1*x1))
    r0=sqrt(k)
    auc=1-x1/2+(qq/(qq+k))*(x1-1)*x1/2+sqrt(k*(k+qq+pp)/pp)/(2*(k+qq))*((k+qq)*atan(sqrt(qq)*(2*x1-1)/r)/(4*sqrt(qq))+(2*x1-1)*r/4+(k+qq)*atan(sqrt(qq)/r0)/(4*sqrt(qq))+r0/4)
    if (case==1) valy=(erf((auc-0.5)/sqrt(2)/s)+1)/2
    if (case==0){
      valy=cdf[1]
      for(i in 2:p*q+1){
        if (auc1[i]>=auc) valy=cdf[i]
      }
    }
    field<-rbind(field,round(c(x,y,(1-valy)),8))
  }
}
final[["field"]]<-field[-1,]
#write(field,paste(outputdirectory,'/outfield.dat',sep=""),quote=FALSE,sep=" ",row.names=FALSE,col.names=FALSE)
}

#Calculation of the p value for the user defined uauc
if(outp){
auc=uauc
if (case==1) valy=(erf((auc-0.5)/sqrt(2)/s)+1)/2
if (case==0){
  valy=cdf[1]
  for(i in 2:p*q+1){
    if (auc1[i]==auc) valy=cdf[i]
  }
}
final[["p"]]<-c(1-valy,uauc,p,q,ifault)
#write(c(1-valy,uauc,p,q,ifault),paste(outputdirectory,'/outp.dat',sep=""))
}

#Calculations based on the k-ellipse passing through h1,f1
if(F1H1){
x=f1
y=h1
xv=x-0.5
yv=y-0.5
k=2*(pp*yv**2+qq*xv**2)-0.5*(pp+qq)+sqrt((2*(pp*yv**2+qq*xv**2)-0.5*(pp+qq))**2+4*pp*qq*(xv-yv)**2)
x1=(pp*qq-sqrt(pp**2*qq**2+qq*(pp+k)*(k*k+qq*(k-pp))))/(2*qq*(pp+k))+0.5
r=sqrt(k+4*qq*(x1-x1*x1))
r0=sqrt(k)
auc=1-x1/2+(qq/(qq+k))*(x1-1)*x1/2+sqrt(k*(k+qq+pp)/pp)/(2*(k+qq))*((k+qq)*atan(sqrt(qq)*(2*x1-1)/r)/(4*sqrt(qq))+(2*x1-1)*r/4+(k+qq)*atan(sqrt(qq)/r0)/(4*sqrt(qq))+r0/4)
mauc=auc
auc=mauc
if (case==1) valy=(erf((auc-0.5)/sqrt(2)/s)+1)/2
if (case==0){
  valy=cdf[1]
  for(i in 2:p*q+1){
    if (auc1[i]==auc) valy=cdf[i]
  }
}
mval=valy
final[["F1H1"]]<-c(1-mval,f1,h1,ifault,p,q,mauc)
#write(c(1-mval,f1,h1,ifault,p,q,mauc),paste(outputdirectory,'/out-F1H1.dat',sep=""))
}
return(final)
#END OF CALCULATIONS
}
