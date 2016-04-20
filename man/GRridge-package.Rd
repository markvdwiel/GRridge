\name{GRridge-package}
\alias{GRridge-package}
\alias{GRridge}
\docType{package}
\title{
Implements adaptive group-regularized (logistic) ridge regression by use of co-data. 
}
\description{
This package implements adaptive group-regularized (logistic, survival, linear) ridge regression by use of co-data. It uses co-data to improve 
predictions of binary, survival and continuous response from high-dimension (e.g. genomics) data. Here, co-data is 
auxiliary information on variables (e.g. genes), such as annotation or p-values from other studies. The package 
includes convenience functions to convert such co-data to the correct input format. In additiion, in includes
functions for evaluating the predictive performance.
}
\details{
\tabular{ll}{
Package: \tab GRridge\cr
Type: \tab Package\cr
Version: \tab 2.2\cr
Date: \tab 2015-12-09\cr
License: \tab GPL \cr
}

\cr
\code{\link{auc}}: Computes Area-under-the-ROC-curve \cr
\code{\link{CreatePartition}}: Creates a partition (groups) of variables \cr
\code{\link{dataFarkas}}: Large data set plus external information \cr
\code{\link{dataVerlaat}}: Data plus external information \cr
\code{\link{grridge}}: Group-regularized (logistic, survival, linear) ridge regression \cr
\code{\link{grridge.cv}}: Cross-validated predictions for a \code{\link[GRridge]{grridge}} (logistic) regression. \cr
\code{\link{matchGeneSets}}: Creates a grouping from gene sets \cr
\code{\link{predict.grridge}}: Predictions for new samples from a \code{\link[GRridge]{grridge}} object \cr
\code{\link{roc}}: Computes an ROC-curve for probabilistic classifiers \cr
}

\author{
Mark A. van de Wiel (mark.vdwiel@vumc.nl)
}
\references{
Mark A. van de Wiel, Tonje G. Lien, Wina Verlaat, Wessel N. van Wieringen, Saskia M. Wilting (2015).  
Better prediction by use of co-data: Adaptive group-regularized ridge regression. 
Statistics in Medicine. Preliminary version: http://arxiv.org/abs/1411.3496.
}

\keyword{ package }
\seealso{\code{GRridge} depends on: \code{\link[penalized]{penalized}. Examples: \code{\link[GRridge]{grridge}}}
}

