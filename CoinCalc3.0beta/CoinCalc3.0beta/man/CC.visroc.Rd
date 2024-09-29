\name{CC.visroc}
\alias{CC.visroc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function: Computation of p-values for ROC analysis
}
\description{
This function computes p-values for the ROC analysis using a k-ellipse-family. Useful to visualize the significance of ROC results.
}
\usage{
CC.visroc(p,q,N,uauc,f1,h1, outCL=TRUE, outp=TRUE, F1H1=TRUE, outfield=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{p}{ 
Number of positive events. 
}
\item{q}{
Number of negative events.
}
\item{N}{ 
Number of segments in which the interval [0,1] is divided for the calculation of the p-values. 
}
\item{uauc}{ 
Value of the Area Under the Curve (AUC) of a Receiver Operating Characteristics (ROC) diagram.
}
\item{f1, h1}{
Values of the point (f_1,h_1) for which the p-value of the k-ellipse should be computed.
}
\item{outCL}{ 
TRUE or FALSE, defining whether the F and H coordinates of the points with p-values = 10\% / 5\% / 1\% should be computed and returned.
}
\item{outp}{ 
TRUE or FALSE, defining whether a vector with: p-value of 'AUC0', input value of the user defined 'AUC0', 'P', 'Q' and error code 'ifault' should be computed and returned.
}
\item{F1H1}{ 
TRUE or FALSE, defining whether a vector with: p-value of the k-ellipse passing through the point (F_1,H_1), 'F_1', 'H_1', value of 'ifault', 'P', 'Q' and AUC of the k-ellipse should be computed and returned.
}
\item{outfield}{ 
TRUE or FALSE, defining whether the coordinates of the lattice points in the ROC diagram for which calculation has been made and the corresponding p-values should becomputed and returned.
}
}
\details{
This function evaluates the significance of a binary prediction for a family of ellipses which are based on confidence and cover the whole Receiver Operating Characteristics space. It follows after a ROC analysis, which can be done with other R packages eg. "MESS".
}
\value{
##Output: list, containing the desired specified elements##
  \item{CL}{Matrix, the F and H coordinates of the points with p-values = 10\% / 5\% / 1\%}
  \item{field}{Matrix, the coordinates of the lattice points in the ROC diagram for which calculation has been made together with the corresponding p-values}
  \item{p}{Vector: p-value of 'AUC0', input value of the user defined 'AUC0', 'P', 'Q', error code 'ifault'}
 \item{F1H1}{Vector: p-value of the k-ellipse passing through the point (F_1,H_1), 'F_1', 'H_1', value of 'ifault', 'P', 'Q', AUC of the k-ellipse}
}
\references{
Sarlis, Nicholas V. (2018): Statistical significance of Earth's electric ans magnetic field variationspreceding earthquakes in Grece and Japan revisited, In:Entropy.
}
\author{
Leonna Szangolies, Potsdam Institute for Climate Impact Research
}
\examples{
# Example:
# First we load example data set 1: April.Tmean.Niederrimbach supplied by CoinCalc
data(CC.Example.Data1)
# Compute input for the VISROC function
input<-CC.for.visroc(CC.Example.Data1[,2],CC.Example.Data1[,3],thres1=0.9)
# Now the VISROC function can be applied
visroc<-do.call(CC.visroc,input)
visroc$p
#Result: The p-value is significant.
# The results can be plotted with the function CC.plot.visroc
do.call(CC.plot.visroc,visroc)
}
