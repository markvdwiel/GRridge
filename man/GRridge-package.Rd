\name{GRridge-package}
\alias{GRridge-package}
\alias{GRridge}
\docType{package}

\title{
Implements adaptive group-regularized (logistic) ridge regression by use of co-data. 
}

\description{ 
This package implements adaptive group-regularized (logistic) ridge regression by use of co-data. It uses co-data to improve 
predictions of binary and continuous response from high-dimension (e.g. genomics) data. Here, co-data is 
auxiliary information on variables (e.g. genes), such as annotation or p-values from other studies. The package 
includes convenience functions to convert such co-data to the correct input format. In additiion, in includes
functions for evaluating the predictive performance.
}

\details{
\tabular{ll}{
Package: \tab GRridge\cr
Type: \tab Package\cr
Version: \tab 1.5.4\cr
Date: \tab 2018-10-31\cr
License: \tab GPL \cr
}

Main functions in the GRridge package are:
\describe{
 \code{auc}: Computes Area-under-the-ROC-curve \cr
 \code{CreatePartition}: Creates a partition (groups) of variables \cr
 \code{dataFarkas}: Large data set plus external information \cr
 \code{dataSimlin}: Small simulated data set with linear response \cr
 \code{dataVerlaat}: Methylation data plus external information \cr
 \code{dataWurdinger}: RNAseq data plus external information \cr
 \code{simlinsmall}: Simulated data for linear regression \cr
 \code{grridge}: Group-regularized (logistic, survival) ridge regression \cr
 \code{grridgeCV}: Cross-validated predictions for a \code{grridge} (logistic, survival) regression. \cr
 \code{matchGeneSets}: Creates a grouping of variables (genes) from gene sets \cr
 \code{mergeGroups}: Merge groups in a partition \cr
 \code{PartitionsSelection}: Co-data selection in a Group-regularized ridge regression model \cr
 \code{predict.grridge}: Predictions for new samples from a \code{grridge} object \cr
 \code{roc}: Computes an ROC-curve for probabilistic classifiers \cr
}
}

\author{
Mark A. van de Wiel (mark.vdwiel@vumc.nl), Putri Novianti (p.novianti@vumc.nl)
}
\references{
Mark van de Wiel, Tonje Lien, Wina Verlaat, Wessel van Wieringen, Saskia Wilting. (2016). 
Better prediction by use of co-data: adaptive group-regularized ridge regression.
Statistics in Medicine, 35(3), 368-81.
 
Novianti PW, Snoek B, Wilting SM, van de Wiel MA (2017). 
Better diagnostic signatures from RNAseq data through use of auxiliary co-data. 
Bioinformatics, 33, 1572-1574.
}

\keyword{ package }
\seealso{\code{GRridge} depends on: \code{penalized}. Examples: \code{grridge}}


