\name{CC.for.visroc}
\alias{CC.for.visroc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function: Performs part of a Receiver Operating Characteristics analysis and returns the Output needed for the CC.visroc analysis
}
\description{
This function performs the ROC analysis as in the R package "MASS" and generates an Output as needed for the VISROC analysis which allows to compute p-values for the ROC analysis using a k-ellipse-family.
}
\usage{
CC.for.visroc(data1, data2, thres1, thres2=seq(from=0.05,to=0.95, by=0.05))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data1}{
Time series A. Numeric vector containing a time series.
}
  \item{data2}{
Time series B. Numeric vector containing a time series.
}
  \item{thres1}{
  Numeric value. Threshold percentile defining events in data1.
  }
  \item{thres2}{
  Numeric vector. Threshold percentiles for data2 to use for the analysis. If not defined this will be 19 steps between 0.05 and 0.95.
  }
}
\details{
 This function should be used as preparation for CC.visroc.
}
\value{
  \item{p}{ Number of positive events. }
\item{q}{Number of negative events.}
\item{N}{ Number of segments in which the interval [0,1] is divided for the calculation of the p-values. }
\item{uauc}{ Value of the Area Under the Curve (AUC) of the Receiver Operating Characteristics (ROC) diagram.}
\item{f1 and h1}{Values of the point (f_1,h_1) which is the optimum point of the diagramm and for which the p-value should be computed in the VISROC analysis.}
}
\author{
Leonna Szangolies, Potsdam Institute for Climate Impact Research
}
\examples{
# Example:
# First we load example data set 1: April.Tmean.Niederrimbach supplied by CoinCalc
data(CC.Example.Data1)
# Compute input for the VISROC function
output<-CC.for.visroc(CC.Example.Data1[,2],CC.Example.Data1[,3],thres1=0.9)
output
# Result: 6 positive events, 53 negative events, 60 segments, AUC= 0.86, optimal point (1,0.12)
# Now the VISROC function can be applied and the results can be plotted
visroc<-do.call(CC.visroc,output)
do.call(CC.plot.visroc,visroc)
}
