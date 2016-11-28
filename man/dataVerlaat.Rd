\name{dataVerlaat}
\alias{dataVerlaat}
\alias{datcenVerlaat}
\alias{CpGann}
\alias{diffmeanFarkas}
\alias{respVerlaat}
\alias{pvalFarkas}
\docType{data}
\title{
Contains 5 R-objects, including the data and the binary response
}
\description{
The five objects are:
datcenVerlaat: methylation data for cervix samples (arcsine-transformed beta values); respVerlaat: binary response;
diffmeanFarkas: effect size external study, Cases minus Controls; pvalFarkas: p-values from external study;
and CpGann: annotation of probes according to location 
}
\usage{data(dataVerlaat)}
\format{
The formats are:

datcenVerlaat: data frame [1:9691,1:44];

respVerlaat: numeric [1:44], 0 = Normal, 1 = Precursor;

diffmeanFarkas: numeric [1:44] ;

pvalFarkas: numeric [1:44];

CpGann: Factor w/ 6 levels "CpG-Island", "North-Shelf", "South-Shelf", "North-Shore", "South-Shore", "Distant"
}

\details{
This data is used for illustration in the statistical paper below.
}

\value{
Five R objects (see description)
}

\references{
Mark van de Wiel, Tonje Lien, Wina Verlaat, Wessel van Wieringen, Saskia Wilting. (2016). 
Better prediction by use of co-data: adaptive group-regularized ridge regression.
Statistics in Medicine, 35(3), 368-81. 
}
\examples{
data(dataVerlaat)
}
\keyword{datasets}

 
