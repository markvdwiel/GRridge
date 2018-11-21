\name{dataSimlin}
\alias{dataSimlin}
\docType{data}
\title{
Small simulated data set with linear response
}

\description{The three objects are:

Y: vector with response values (linear)

simdata: matrix with simulated variables

part5: list with 5 components representing partition of variables
}

\value{
Three R objects (see Description)
}

\usage{data(dataSimlin)}



\details{
This is simulated data allowing a quick test for grridge
}

\examples{
data(dataSimlin)

#apply grridge
grsim <- grridge(highdimdata=simdata, response=Y, partitions=part5, unpenal = ~1, innfold=10)

#apply CV
grsimcv <- grridgeCV(grsim,simdata,Y,outerfold=5)

}

\keyword{datasets}
