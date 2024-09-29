\name{CC.esyncr.ts}
%-------------------------------------------------------------------------------------------------
\alias{CC.esyncr.ts}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function: Event Synchronisation Rates for varying tau for two time series 
}
\description{
This function calculates the event synchronisation rates using two data sets of event time series format (ts) for different taus and the corresponding sigificance level.
}
\usage{
CC.esyncr.ts(a, b=a, alpha=0.05, max_delT=min(length(a),length(b)), sym=FALSE, offset=30, siglevel="TRUE", reps=1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{a}{
Time series A. Numeric vector or one-row-matrix containing a time series of binary data (0,1) with all time steps of interest, where a '1' denotes a time step with, and a '0' a time step without event.
SeriesA and seriesB need to be of same length. Both may contain missing values.
}
  \item{b}{
Time series B. Numeric vector or one-row-matrix containing a time series of binary data (0,1) with all time steps of interest, where a '1' denotes a time step with, and a '0' a time step without event.
SeriesA and seriesB need to be of same length. Both may contain missing values.
}
  \item{alpha}{
Value between [0 1], default=0.05. 
The significance level 'alpha' of the significance test. null hypothesis: The found coincidence rates result from randomly, independently distributed event time series.
}
\item{max_delT}{
Positive value [1 inf], default = min(length(a),length(b)). 
Maximum tolerance window delta T. delT is adapted for every timestep within the function, here a maximum value max_delT describing the maximum temporal tolerance window in which coincidences may occur can be defined.
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
Logical. Defining wether a sigificance level should be calculated using shuffled time series.
}
  \item{reps}{
    Positive whole number [100, inf], default = "1000". 
    Number of surrogate/shuffled time series to be produced for the "surrogate" or "shuffle" significance test. 
}
}

\value{
## Output: list, containing five elements ##
\item{tau}{                 Series of taus}
\item{Q synchronisation rates}{Numeric vector, series of Q synchronisation rates for different taus}
\item{q synchronisation rates}{  Numeric vector, series of q synchronisation rates for different taus}
\item{Q significance levels}{                    Numeric vector, series of Q synchronisation rates which are significance levels}
\item{q significance levels}{                   Numeric vector, series of q synchronisation rates which are significance levels}


}
\references{
Siegmund, J., Seigmund, N. and Donner, R. (2015): CoinCalc - An R package for quantifying simultaneities 
of events in two event (time) series. Submitted to Computers and Geosciences.
Donges, J. F., Schleussner, C. F., Siegmund, J. F., and Donner, R. V. (2015): Coincidence analysis
for quantifying statistical interrelationships between event time series, arXiv:1508.03534.
}
\author{
Leonna Szangolies, Potsdam Institute for Climate Impact Research
}
