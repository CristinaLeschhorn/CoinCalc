\name{CC.ecr.es}
%-------------------------------------------------------------------------------------------------
\alias{CC.ecr.es}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function: Event Coincidence Rates for varying tau for two event sequences
}
\description{
This function calculates the event coincidence rates using two data sets of event sequence format (es) for different taus and the corresponding sigificance level.
}
\usage{
CC.ecr.es(seriesA, spanA, seriesB=seriesA, spanB=spanA, alpha=0.05, delT=0, sym=FALSE, offset=30, siglevel="TRUE",sigtest="shuffle.surrogate", reps=1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{seriesA}{
Time series A. Numeric vector containing a set of values specifying only the (e.g. time-) steps that denote an event.
SeriesA and seriesB do not need to be of same length. Both may contain missing values.
}
  \item{seriesB}{
Time series B. Numeric vector containing a set of values specifying only the (e.g. time-) steps that denote an event.
SeriesA and seriesB do not need to be of same length. Both may contain missing values.
}
  \item{spanA}{
Vector containing two integers. spanA specifies the start and end point of the data set given in seriesA.
}
  \item{spanB}{
Vector containing two integers. spanB specifies the start and end point of the data set given in seriesB.
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
  \item{offset}{
Integer defining the maximum tau to be used. Event rates for tau in [0 offset] are calculated.
}
  \item{siglevel}{
Logical [TRUE FALSE]. Defining wether a sigificance level should be calculated using shuffled time series.
}
  \item{sigtest}{
    Character specifying the type of significance test, default = "shuffle.surrogate". 
    "wt.surrogate": significance test based on the calculation of event coincidence analysis for a large number of surrogate time series having the same average waiting time between two events.
    "shuffle.surrogate": significance test based on the calculation of event coincidence analysis for a large number of randomly shuffled time series having the same number of events like the original time series.
}
  \item{reps}{
    Positive whole number [100, inf], default = "1000". 
    Number of surrogate/shuffled time series to be produced for the "surrogate" or "shuffle" significance test. 
}
}

\value{
## Output: list, containing five elements ##
\item{tau}{                 Series of taus}
\item{precursor coincidence rates}{Numeric vector, series of precursor coincidence rates for different taus}
\item{trigger coincidence rates}{  Numeric vector, series of trigger coincidence rates for different taus}
\item{precursor significance levels}{                    Numeric vector, series of precursor coincidence rates which are significance levels}
\item{trigger significance levels}{                   Numeric vector, series of trigger coincidence rates which are significance levels}


}
\references{
Siegmund, J., Seigmund, N. and Donner, R. (2015): CoinCalc - An R package for quantifying simultaneities 
of events in two event (time) series. Submitted to Computers and Geosciences.
Donges, J. F., Schleussner, C. F., Siegmund, J. F., and Donner, R. V. (2015): Coincidence analysis
for quantifying statistical interrelationships between event time series, arXiv:1508.03534.
}
\author{
Leonna Szangolies, Jonatan Siegmund, Potsdam Institute for Climate Impact Research
}


\examples{

library(CoinCalc3.0)
data(CC.Example.Data1)
# Then, we binarize both time series:
Flow_bin=CC.binarize(CC.Example.Data1$DOY_Lilac_Flowering,ev.def="percentile",thres=0.1,event="lower")
Tmean_bin=CC.binarize(CC.Example.Data1$April_Tmean,ev.def="percentile",thres=0.9,event="higher")
out <- CC.ecr.ts(Flow_bin,Tmean_bin)
plot(out$tau,out$'precursor rates','l',main='precursor coincidence rates')
lines(out$tau,out$'precursor significance level','l',col='red')
legend("topright",lwd=1,col="red",legend="significance level")

# Result: Some tau definitions have coincidence rates over the significance level.

}
