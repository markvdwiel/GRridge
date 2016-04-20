\name{grridge.cv}
\alias{grridge.cv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Returns the cross-validated predictions 
}
\description{
Returns the cross-validated predictions for a \code{\link[GRridge]{grridge}} (logistic) regression. 
}
\usage{
grridge.cv(grr, highdimdata, response, outerfold = length(response), fixedfolds = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{grr}{
List. Output of GRridge function.
}
  \item{highdimdata}{
Matrix or numerical data frame. Contains the primary data of the study. Columns are samples, rows are variables (features).
}
  \item{response}{
Factor, numeric or binary. Response values. The \code{length(response)} should equal \code{ncol(highdimdata)}.
}
  \item{outerfold}{
Integer. Fold used for cross-validation loop.
}
  \item{fixedfolds}{
Boolean. Use fixed folds for cross-validation?
}
}
\details{
This convenience function returns cross-validated predictions from \code{\link[GRridge]{grridge}}, including those from 
ordinary (logistic) ridge regression. It can be used to compute ROC-curves. 
}
\value{
A matrix containing the predictions. The first column contains the sample indices, 
the second the prediction by ordinary ridge, the third the predictions by group-regularized ridge, the 
fourth (optional) the predictions by group-regularized ridge plus selection. Finally, it may contain
predictions by lasso and/or a regression model with unpenalized covariates only. 
}
\references{
Mark A. van de Wiel, Tonje G. Lien, Wina Verlaat, Wessel N. van Wieringen, Saskia M. Wilting (2015).  
Better prediction by use of co-data: Adaptive group-regularized ridge regression. 
Statistics in Medicine. Preliminary version: http://arxiv.org/abs/1411.3496.
}
\author{
Mark A. van de Wiel
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
ROC-curves: \code{\link[GRridge]{roc}}. AUC: \code{\link[GRridge]{auc}}. Example: \code{\link[GRridge]{grridge}}.
}

