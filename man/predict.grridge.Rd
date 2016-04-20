\name{predict.grridge}
\alias{predict.grridge}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Predictions for new samples
}
\description{
Returns predictions for new samples for a \code{\link[GRridge]{grridge}} object
}
\usage{
predict.grridge(grr, datanew, dataunpennew=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{grr}{
List. Output of GRridge function.
}

  \item{datanew}{
Vector or data frame. Contains the new data. For a data frame: columns are samples, rows are variables (features).
}

\item{dataunpennew}{
Vector or data frame. Optional new data for unpenalized variables. NOTE: columns are covariates, rows are samples.
}
}
\details{
This function returns predictions of the response using the \code{\link[GRridge]{grridge}} output. It should be applied 
to samples NOT used for fitting of the models.
}
\value{
A matrix containing the predictions from all models available in \code{grr}. 
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Mark A. van de Wiel
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
Cross-validated predictions: \code{\link[GRridge]{grridge.cv}}. Examples: \code{\link[GRridge]{grridge}}.
}

