\name{roc}
\alias{roc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Produces ROC curve for probabilistic classifiers (e.g. logistic regression)
}
\description{
Computes an ROC-curve for probabilistic classifiers.}
\usage{
roc(probs, true, cutoffs)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{probs}{
Numeric vector, with values between 0 and 1
}
  \item{true}{
Binary vector.
}
  \item{cutoffs}{
Numeric vector, with increasing values between 0 and 1.
}
}
\details{
The vector probs contains predicted probabilities for the response to equal 1, as produced by 
a probabilistic classifier like logistic regression. The cutoffs can simply represent a grid of values between 0 and 1. 
}
\value{
A matrix with two rows which contain corresponding False Positive and True Positive Rates for all cutoffs.
}

\author{
Mark A. van de Wiel
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
For area-under-the ROC-curve: \code{\link[GRridge]{auc}}. Examples: \code{\link[GRridge]{grridge}}.
}

