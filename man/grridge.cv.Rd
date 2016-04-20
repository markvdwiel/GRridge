\name{grridge.cv}
\alias{grridge.cv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Returns the cross-validated predictions
}
\description{
Returns the cross-validated predictions for a \code{\link[GRridge]{grridge}} logistic, linear or Cox regression. 
}
\usage{
grridge.cv(grr, highdimdata, response, outerfold = length(response), fixedfolds = TRUE, recalibrate=FALSE)
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
Factor, numeric, binary or survival. Response values. The number of response values should equal \code{ncol(highdimdata)}.
}
  \item{outerfold}{
Integer. Fold used for cross-validation loop.
}
  \item{fixedfolds}{
Boolean. Use fixed folds for cross-validation?
}

\item{recalibrate}{
Boolean. Should the prediction model be recalibrated on the test samples? Only implemented for 
logistic and linear regression with only penalized covariates.
}
}
\details{
This convenience function returns cross-validated predictions from \code{\link[GRridge]{grridge}}, including those from 
ordinary (logistic) ridge regression. It can be used to compute ROC-curves. About \code{recalibrate}: this option allows to compare recalibrated models, but only if the test sample size is large enough. See 
\code{\link[GRridge]{predict.grridge}}
}
\value{
For linear and logistic regression: A matrix containing the predictions. The first column contains the sample indices, 
the second the prediction by ordinary ridge, the third the predictions by group-regularized ridge, the 
fourth (optional) the predictions by group-regularized ridge plus selection. Finally, it may contain
predictions by lasso and/or a regression model with unpenalized covariates only. 
For Cox regression: a list with three components. First, a matrix with the linear predictions, in the same format as above.
Second, a list of survival probability matrices (one for each prediction model), the rows of which indicate the individuals and the columns indicate the observed times. Third, a vector with the observed times. 
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
For logistic regression: ROC-curves: \code{\link[GRridge]{roc}}. AUC: \code{\link[GRridge]{auc}}. 

Example: \code{\link[GRridge]{grridge}}.

For survival: additional R-scripts for performance evaluation ([integrated] Brier scores, [integrated] survival AUC) are available from the author's web-site. 
}

