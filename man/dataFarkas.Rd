\name{dataFarkas}
\alias{dataFarkas}
\alias{CpGannFarkas}
\alias{datcenFarkas}
\alias{respFarkas}
\docType{data}
\title{
Contains three R-objects, including the data and the binary response

}
\description{The three objects are:
datcenFarkas: methylation data for cervix samples (arcsine-transformed beta values); respFarkas: binary response;and CpGannFarkas: annotation of probes according to location 
}

\format{
\itemize{
\item{datcenFarkas: data frame [1:40000,1:37]}
\item{respfarkas: Factor w/ 2 levels "Normal","Precursor"}
\item{CpGannFarkas: Factor w/ 6 levels "CpG-Island", "North-Shelf", "South-Shelf", "North-Shore", "South-Shore", "Distant"}
}
}

\details{
This data is used for illustration in the statistical paper below.
}
\value{
Three R objects, i.e. a matrix contains methylation data, a vector contains binary responses and an annotation matrix. 
}
\source{
Farkas, S. et al. (2013). Genome-wide DNA methylation assay reveals novel candidate biomarker
genes in cervical cancer. Epigenetics, 8, 1213.
}
\references{
Mark van de Wiel, Tonje Lien, Wina Verlaat, Wessel van Wieringen, Saskia Wilting. (2016). 
Better prediction by use of co-data: adaptive group-regularized ridge regression.
Statistics in Medicine, 35(3), 368-81.
}
\examples{
data(dataFarkas)
}
\keyword{datasets}
