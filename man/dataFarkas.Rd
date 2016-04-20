\name{dataFarkas}
\alias{dataFarkas}
\docType{data}
\title{
Contains three R-objects, including the data and the binary response

}
\description{The three objects are:
datcenFarkas: methylation data for cervix samples (arcsine-transformed beta values); respFarkas: binary response; 
and CpGannFarkas: annotation of probes according to location 
}
\usage{data(respFarkas)}
\format{
The formats are:
datcenFarkas: data frame [1:40000,1:37]\\
respfarkas: Factor w/ 2 levels "Normal","Precursor"\\
CpGannFarkas: Factor w/ 6 levels "CpG-Island", "North-Shelf", "South-Shelf", "North-Shore", "South-Shore", "Distant"
}
\details{
This data is used for illustration in the statistical paper below.
}
\source{
Farkas, S. et al. (2013). Genome-wide DNA methylation assay reveals novel candidate biomarker
genes in cervical cancer. Epigenetics, 8, 1213.
}
\references{
Mark A. van de Wiel, Tonje G. Lien, Wina Verlaat, Wessel N. van Wieringen, Saskia M. Wilting (2015).  
Better prediction by use of co-data: Adaptive group-regularized ridge regression. 
Statistics in Medicine. Preliminary version: http://arxiv.org/abs/1411.3496.
}
\examples{
data(dataFarkas)
## maybe str(respFarkas) ; plot(respFarkas) ...
}
\keyword{datasets}
