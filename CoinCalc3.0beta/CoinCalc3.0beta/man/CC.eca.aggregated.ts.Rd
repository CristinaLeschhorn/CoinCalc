\name{CC.eca.aggregated.ts}
\alias{CC.eca.aggregated.ts}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function: Perform aggregated ECA for event timeseries
}
\description{
This functions performs an aggregated Event Coincidence Analysis for a sequence of pairs of binary timeseries.
}
\usage{
CC.eca.aggregated.ts(seriesList, alpha = 0.05, delT = 0, sym = FALSE, tau = 0, sigtest = "poisson", reps = 1000) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{seriesList}{ 
  A list where each entry is a numeric matrix with a pair of binary timeseries (0,1). A '1' denotes a time step with, and a '0' a time step without event. 
}
  \item{alpha}{
Value between [0 1], default=0.05. 
The significance level 'alpha' of the significance test. null hypothesis: The found coincidence rates result from randomly, independently distributed event time series.
}
  \item{delT}{
Positive value [1 inf], default = 1. 
Tolerance window delta T. delT describes the temporal tolerance window in which coincidences may occur. 
If e.g. delT=3, an event in seriesA occurring simultaneously or one, two or three time steps after an event in seriesB
also counts for a coincidence (K). With delT=0 only simultaneous events are counted for coincidences, for e.g. delT=1 the tolerance window is extended by one time step. 
}
  \item{sym}{
  TRUE or FALSE, default = FALSE.
  FALSE: non-symmetrically tolerance window. In case of the precursor coincidence, the window starts at t-tau and ends at t-tau-delT, for trigger coincidence, the window starts at t+tau and ends at t+tau+delT.
  TRUE: symmetrical tolerance window. The tolerance window is symmetrically placed around t-tau and t+tau, respectively. For precursor coincidence, the window starts at t-tau+delT and ends at t-tau-delT, for trigger coincidence the window starts at t+tau-delT and ends at t+tau+delT
 }
  \item{tau}{
Positive value [0 inf], default = 0. 
Time lag parameter tau allows for shifted events to be counted for coincidences (K). If e.g. tau=5, an event
in seriesA occurring exactly 5 time steps after an event in seriesB also counts for coincidence (K).
}
  \item{sigtest}{
    Character specifying the type of significance test, default = "poisson". 
    "poisson": significance test based on the calculation of the probability that the empirical coincidence rate would occur for two random (poisson-process) time series.
    "shuffle.surrogate": significance test based on the calculation of event coincidence analysis for a large number of randomly shuffled time series having the same number of events like the original time series.
}
  \item{reps}{
    Positive whole number [100, inf], default = "1000". 
    Number of surrogate/shuffled time series to be produced for the "surrogate" or "shuffle" significance test. 
}
}
\details{
For the coupled ECA the coincidences are counted separately and are summed up afterwards.
}
\value{
## Output: list, containing six elements ##
\item{NH precursor}{                 Logical, null hypothesis for precursor coincidence rate accepted}
\item{NH trigger}{                Logical, null hypothesis for trigger coincidence rate accepted}
\item{p-value precursor}{                    Numeric, p-value used for the precursor significance test}
\item{p-value trigger}{                   Numeric, p-value used for the trigger significance test}
\item{precursor coincidence rate}{Numeric, precursor coincidence rate}
\item{trigger coincidence rate}{  Numeric, trigger coincidence rate}
}
\references{
Schleussner, Carl-Friedrich, Jonathan F. Donges, Reik V. Donner, and Hans Joachim Schellnhuber (2016): Armed-conflict risks enhanced by climate-related disasters in ethnically fractionalized countries, In:PNAS(113 (33)), pp. 9216-9221.
}
\author{
Leonna Szangolies, Potsdam Institute for Climate Impact Research
}
\examples{
# Example:
# First we load example data set 3 with Land surface temperature and NDVI of several spots in italy, supplied by CoinCalc.
data(CC.Example.Data3)
# Binarize data
x1<-CC.bin.array(CC.Example.Data3$LST,dims=c("x","y","t"),ev.def = "percentile",thres = 0.95,event = "higher")
x2<-CC.bin.array(CC.Example.Data3$NDVI,dims=c("x","y","t"),ev.def = "percentile",thres = 0.05,event = "lower")
# Transform data to usable format
newdata=list()
ix=1
for(i in 1:131){
  for(j in 1:131){
  one<-matrix(ncol=2,nrow=96)
  one[,1]<-x1[i,j,]
  one[,2]<-x2[i,j,]
  if(!all(is.na(one[,1])) && !all(is.na(one[,2]))){
  newdata[[ix]]<-one
  ix=ix+1
  }
  }
}
# Apply aggregated ECA
CC.eca.aggregated.ts(newdata,sigtest="shuffle.surrogate")
# Result: The aggregated coincidance rates are significant even though the precursor coincidence rate is only 0.18 and the trigger coincidence rate is 0.17. However the sample size is large.
}
