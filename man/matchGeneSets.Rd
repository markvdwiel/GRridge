\name{matchGeneSets}
\alias{matchGeneSets}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creates a grouping of variables (genes) from gene sets
}
\description{
Creates a grouping of variables (genes) from gene sets by matching the IDs of the genes with 
the IDs of the members of the gene sets
}
\usage{
matchGeneSets(GeneIds, GeneSets, minlen = 25, remain = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{GeneIds}{
Character vector. Vector of gene IDs. Can be any ID (gene symbol, entrezID, etc) as long as it matches with
IDs used in \code{GeneSets}.
}
  \item{GeneSets}{
Named list of character vectors. Each component of the list represents a named gene set. Each vector a list of member genes.
}
  \item{minlen}{
Integer. Minimum number of members of a gene set.
}
  \item{remain}{
Boolean. If \code{remain=TRUE}, all genes that are not in any list will be grouped in one remainder group. 
}
}
\details{
About \code{minlen}: to avoid overfitting in the \code{\link[GRridge]{grridge}} function, we recommend to not use groups with less than 25 members, unless \code{monotone=TRUE} is used in the \code{\link[GRridge]{grridge}} function, in which case 10 members may suffice as a lower bound. About \code{remain}:
it is often beneficial to down-weight genes that are not part of any gene set, so we recommend to use 
\code{remain=TRUE}
}
\value{
A list the components of which contain the indices of the variables belonging to each of the groups.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Mark A. van de Wiel
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[GRridge]{grridge}}, \code{\link[GRridge]{CreatePartition}}
}
