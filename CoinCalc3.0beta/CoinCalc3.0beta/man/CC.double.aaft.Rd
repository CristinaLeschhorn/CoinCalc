\name{CC.double.aaft}
\alias{CC.double.aaft}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function: Coupled Amplitude Adjusted Fourier Transformation for two timeseries
}
\description{
This function performs a bivariate Amplitude Adjusted Fourier Transformation for two timeseries of the same length.
}
\usage{
CC.double.aaft(data1,data2,steps=c(1:length(data1)), nsur=100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data1}{
Time series A. Numeric vector containing a time series.
SeriesA and seriesB need to be of same length. Both may contain missing values.
}
  \item{data2}{
Time series B. Numeric vector containing a time series.
SeriesA and seriesB need to be of same length. Both may contain missing values.
}
 \item{steps}{
 Steps to analyse. Numeric vector containing the vector indices of interest.
 }
 \item{nsur}{
 Positive value indicating the number of surrogates to be returned.
 }
}
\details{
The function is designed for evenly-spaced timeseries with the same second-order properties. The coupling of the two timeseries in the AAFT is achieved by performing the same phase-scrambling on both timeseries.
}
\value{
##Output: list, containing two elements##
 \item{data1.series}{List, nsur bivariate surrogates of data1}
 \item{data2.series}{List, nsur bivariate surrogates of data2}
}
\references{
Barnett, Adrian, Peter Baker, and Maintainer Adrian Barnett (2012): Package 'season', RDocumentation
}
\author{
Leonna Szangolies, Potsdam Institute for Climate Impact Research
}
\examples{
# Example:
# First we load example data set 1: April.Tmean.Niederrimbach supplied by CoinCalc
data(CC.Example.Data1)
#Apply coupled AAFT
output = CC.double.aaft(CC.Example.Data1[,2],CC.Example.Data1[,3])
#Perform ECA on the surrogates
prec<-vector(length=length(output$data1.series[1,]))
trig<-vector(length=length(output$data1.series[1,]))
for(i in 1:length(output$data1.series[1,])){
  a<-CC.binarize(output$data1.series[,i])
  b<-CC.binarize(output$data2.series[,i],thres=0.1,event="lower")
  eca<-CC.eca.ts(a,b)
  prec[i]<-eca[[1]]
  trig[i]<-eca[[2]]
}
result<-list("True NH Precursor"=sum(prec),"False NH Precursor"=length(prec)-sum(prec),"True NH Trigger"=sum(trig),"False NH Trigger"=length(trig)-sum(trig))
result
#Result: For the majority of the 100 surrogate pairs the null hypothesis for precurser and trigger rates are rejected and therefor the coincidence rates are significant
}
