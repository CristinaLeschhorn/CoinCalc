\name{CC.event.sync.ts}
%- Also NEED an '\alias' for EACH other topic documented here.
\alias{CC.event.sync.ts}
\title{
Function: Computation of Event Synchronization with significance test
}
\description{
This function performs Event Synchronisation with self-adapting delT.
}
\usage{
CC.event.sync.ts(a,b,max_delT=min(length(a),length(b)),tau=0,sigtest=TRUE,reps=1000,alpha=0.05)
}
\arguments{
  \item{a}{Time series A. Numeric vector containing a time series of binary data (0,1) with all time steps of interest, where a '1' denotes a time step with, and a '0' a time step without event.
}
  \item{b}{
Time series B. Numeric vector containing a time series of binary data (0,1) with all time steps of interest, where a '1' denotes a time step with, and a '0' a time step without event.
}
\item{max_delT}{
Positive value [1 inf], default = min(length(a),length(b)). 
Maximum tolerance window delta T. delT is adapted for every timestep within the function, here a maximum value max_delT describing the maximum temporal tolerance window in which coincidences may occur can be defined.
}
\item{tau}{
Positive value [0 inf], default = 0. 
Time lag parameter tau allows for shifted events to be counted for coincidences (K). If e.g. tau=5, an event in seriesA occurring exactly 5 time steps after an event in seriesB also counts for coincidence (K).
}
\item{sigtest}{
Logical. Defines whether significance test should be performed. Currently only a shuffle test is available.  
}
\item{reps}{
Positive whole number [100, inf], default = "1000". 
    Number of surrogate/shuffled time series to be produced for the "shuffle" significance test.
}
\item{alpha}{
Value between [0 1], default=0.05. 
The significance level 'alpha' of the significance test. null hypothesis: The found synchronisation rates result from randomly, independently distributed event time series.
}
}
\details{
CC.event.sync performs Event Synchronisation with self-adapting tolerance window.
}
\value{
## Output: list, containing six elements ##
\item{NH Q}{                 Logical, null hypothesis for synchronisation rate accepted}
\item{NH q}{                Logical, null hypothesis for preceding rate accepted}
\item{p-value Q}{                    Numeric, p-value used for the synchronisation significance test}
\item{p-value q}{                   Numeric, p-value used for the preceding significance test}
\item{Synchronisation rate Q}{ Numeric [0,1], synchronisation rate, Q=1 if the events are fully synchronized}
\item{Preceding rate q}{  Numeric [-1,1], preceding rate, q=1 if events in a always precede those in b}
}
\references{
Quian Quiroga, R., Kreuz, T., Grassberger, P. (2002): Event synchronization: a simple and fast method to measure synchronicity and time delay patterns, In: Physical Review
}
\author{
Leonna Szangolies, Potsdam Institute for Climate Impact Research
}
\examples{
# Example:
# First we load example data set 1 supplied by CoinCalc
data(CC.Example.Data1)
# Binarize the timeseries
Flow_bin=CC.binarize(CC.Example.Data1$DOY_Lilac_Flowering,ev.def="percentile",thres=0.1,event="lower")
Tmean_bin=CC.binarize(CC.Example.Data1$April_Tmean,ev.def="percentile",thres=0.9,event="higher")
# Compute the Synchronisation of events
sync<-CC.event.sync.ts(Tmean_bin,Flow_bin)
sync
# Result: There is a significant synchronization of the events and the events of April.Tmean are significantly always preceding the Lilac Flowering
}
