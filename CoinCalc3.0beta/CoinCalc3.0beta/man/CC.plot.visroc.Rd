\name{CC.plot.visroc}
\alias{CC.plot.visroc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function: Visualizes the significance of ROC results using the Output of CC.visroc
}
\description{
This function plots the Output of the VISROC anlysis, needs this as Input. 
}
\usage{
CC.plot.visroc(CL=c(),field=c(),p=c(),F1H1=c())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{CL}{Matrix, the F and H coordinates of the points with p-values = 10\% / 5\% / 1\%}
  \item{field}{Matrix, the coordinates of the lattice points in the ROC diagram for which calculation has been made together with the corresponding p-values}
  \item{p}{Vector: p-value of 'AUC0', input value of the user defined 'AUC0', 'P', 'Q', error code 'ifault'}
 \item{F1H1}{Vector: p-value of the k-ellipse passing through the point (F_1,H_1), 'F_1', 'H_1', value of 'ifault', 'P', 'Q', AUC of the k-ellipse}
}
\author{
Leonna Szangolies, Potsdam Institute for Climate Impact Research
}
\examples{
# Example:
# First we load example data set 1: April.Tmean.Niederrimbach supplied by CoinCalc
data(CC.Example.Data1)
# Compute input for the VISROC function and apply the VISROC function
output<-CC.for.visroc(CC.Example.Data1[,2],CC.Example.Data1[,3],thres1=0.9)
visroc<-do.call(CC.visroc,output)
# Plot the VISROC results to visualize the significance of the ROC analysis
do.call(CC.plot.visroc,visroc)
}
