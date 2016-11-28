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
Version: \tab 0.99.1\cr
Date: \tab 2015-05-12\cr
License: \tab GPL \cr
}

Main functions in the GRridge package are:
\describe{
 \code{auc}: Computes Area-under-the-ROC-curve \cr
 \code{CreatePartition}: Creates a partition (groups) of variables \cr
 \code{dataFarkas}: Large data set plus external information \cr
 \code{dataVerlaat}: Data plus external information \cr
 \code{grridge}: Group-regularized (logistic) ridge regression \cr
 \code{grridge.cv}: Cross-validated predictions for a \code{grridge} (logistic) regression. \cr
 \code{matchGeneSets}: Creates a grouping of variables (genes) from gene sets \cr
 \code{mergeGroups}: Merge groups in a partition \cr
 \code{PartitionsSelection}: Co-data selection in a Group-regularized ridge regression model \cr
 \code{predict.grridge}: Predictions for new samples from a \code{grridge} object \cr
 \code{roc}: Computes an ROC-curve for probabilistic classifiers \cr
}
}

\author{
Mark A. van de Wiel (mark.vdwiel@vumc.nl)
}
\references{
Mark van de Wiel, Tonje Lien, Wina Verlaat, Wessel van Wieringen, Saskia Wilting. (2016). 
Better prediction by use of co-data: adaptive group-regularized ridge regression.
Statistics in Medicine, 35(3), 368-81.
}

\keyword{ package }
\seealso{\code{GRridge} depends on: \code{penalized}. Examples: \code{grridge}}


