\name{auc}
\alias{auc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Area-under-the-ROC-curve
}
\description{
Computes Area-under-the-ROC-curve}
\usage{
auc(rocout)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rocout}{
Matrix with two rows
}
}
\details{
Operates on the output of \code{\link[GRridge]{roc}}. The rows of the input matrix represent the False Positive Rates (FPR) and 
corresponding True Positive Rates (TPR) at fixed thresholds. 
}
\value{
Numeric. Interpolated Area-under-the-ROC-curve.
}

\author{
Mark A. van de Wiel
}

\seealso{
ROC-curves: \code{\link[GRridge]{roc}}. Examples: \code{\link[GRridge]{grridge}}.
}

