\name{dataWurdinger}
\alias{dataWurdinger}
\alias{datWurdinger_BC}
\alias{respWurdinger}
\alias{annotationWurdinger}
\alias{coDataWurdinger}
\docType{data}
\title{
R-objects related to the mRNAseq data 
}

\description{The four objects are:

datWurdinger_BC: A matrix containing preprocessed mRNA sequencing data (quasi-gausian scale, normalized). Columns are samples (81 samples with Breast Cancer and Colorectal Cancer) and rows are features (18410 features).

respWurdinger: A factor containing responses for samples with Breast cancer (n=40) and  colorectal cancer (n=41)

annotationWurdinger: A list containing ensembleID, geneSymbol, entrezID and chromosome location.; 

coDataWurdinger: A list containing co-data sets from external sources, namely (i) a list of genes that are expressed in platelets; (ii) immunologic signature pathway and; (iii) transcription factor based pathway and a list of cancer somatic genes.
}

\value{
Four R objects (see Description)
}

\usage{data(dataWurdinger)}

\source{
Wurdinger,T., et al (2015). RNA-seq of tumor-educated platelets enables blood-based pan-cancer, multiclass, and molecular pathway cancer diagnostics. Cancer Cell, 28(5), 666-676.
}


\details{
This data is used for illustration in the application paper below.
}

\examples{
data(dataWurdinger)

# Transform the data set to the square root scale
dataSqrtWurdinger <- sqrt(datWurdinger_BC)
#
#Standardize the transformed data
datStdWurdinger <- t(apply(dataSqrtWurdinger,1,
                           function(x){(x-mean(x))/sd(x)}))
#
# A list of gene names in the primary RNAseq data
genesWurdinger <- annotationWurdinger$geneSymbol

}

\keyword{datasets}
