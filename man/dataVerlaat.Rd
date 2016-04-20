\name{dataVerlaat}
\alias{dataVerlaat}
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
\usage{data(datcenVerlaat)}
\format{
The formats are:
datcenVerlaat: data frame [1:9691,1:44]\\
respVerlaat: numeric [1:44], 0 = Normal, 1 = Precursor\\
diffmeanFarkas: numeric [1:44] \\
pvalFarkas: numeric [1:44]\\
CpGann: Factor w/ 6 levels "CpG-Island", "North-Shelf", "South-Shelf", "North-Shore", "South-Shore", "Distant"
}
\details{
This data is used for illustration in the statistical paper below.
}

\references{
Mark A. van de Wiel, Tonje G. Lien, Wina Verlaat, Wessel N. van Wieringen, Saskia M. Wilting (2015).  
Better prediction by use of co-data: Adaptive group-regularized ridge regression. 
Statistics in Medicine. Preliminary version: http://arxiv.org/abs/1411.3496.
}
\examples{
data(dataVerlaat)
## maybe str(datcenVerlaat) ; plot(datcenVerlaat) ...
}
\keyword{datasets}

 
