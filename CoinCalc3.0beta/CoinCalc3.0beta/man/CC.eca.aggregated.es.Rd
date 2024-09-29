\name{CC.eca.aggregated.es}
%- Also NEED an '\alias' for EACH other topic documented here.
\alias{CC.eca.aggregated.es}
\title{
Function: Perform aggregated ECA for event series
}
\description{
This functions performs an aggregated Event Coincidence Analysis for a sequence of pairs of event series.
}
\usage{
CC.eca.aggregated.es(seriesListA,seriesListB,spanListA,spanListB, alpha = 0.05, delT = 0, sym = FALSE, tau = 0, sigtest = "poisson", reps = 1000) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{seriesListA}{ 
  A list containing the different series of event timings.
}
\item{seriesListB}{ 
  A list containing the to the series A corresponding series B of event timings in the same order as in seriesListA 
}
\item{spanListA}{
A list containing the spans (start and end points) of the series A's. Spans are vectors of two integers. Has to have the same order as seriesListA. 
}
\item{spanListB}{
A list containing the spans (start and end points) of the series B's. Spans are vectors of two integers. Has to have the same order as seriesListB. 
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
# We use a simple artificial example here. 
artA=list(c(5,10,15,20,25,30,35,40,45,50,55,60,65,70,75),c(4,6,20,23,28,38,39,45,50,53,71,79))
artB=list(c(12,22,32,42,52,62,72),c(6,19,23,45,47,53,80))
span_artA=list(c(1,100),c(1,99))
span_artB=list(c(1,99),c(1,100))
# Apply aggregated ECA
CC.eca.aggregated.es(artA,artB,span_artA,span_artB)
# Result: The aggregated coincidance rates are not significant.
}
