\name{matchGeneSets}
\alias{matchGeneSets}
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
About \code{minlen}: to avoid overfitting in the \code{\link{grridge}} function, we recommend to not use groups with less than 25 members, unless \code{monotone=TRUE} is used in the \code{\link{grridge}} function, in which case 10 members may suffice as a lower bound. About \code{remain}:
it is often beneficial to down-weight genes that are not part of any gene set, so we recommend to use 
\code{remain=TRUE}
}
\value{
A list the components of which contain the indices of the variables belonging to each of the groups.
}

\author{
Mark A. van de Wiel
}

\seealso{
\code{\link{grridge}}, \code{\link{CreatePartition}}
}

\examples{

# Load data objects
data(dataWurdinger)

# Transform the data set to the square root scale
dataSqrtWurdinger <- sqrt(datWurdinger_BC)

#Standardize the transformed data
datStdWurdinger <- t(apply(dataSqrtWurdinger,1,function(x){(x-mean(x))/sd(x)}))

# A list of gene names in the primary RNAseq data
genesWurdinger <- annotationWurdinger$geneSymbol

# We  show an example of GRridge classification model by using overlapping groups,
# i.e. pathway-based grouping. Transcription factor based pathway was extracted from 
# the MSigDB (Section C3: motif gene sets; subsection: transcription factor targets;
# file's name: "c3.tft.v5.0.symbols.gmt"). The gene sets are based on 
# TRANSFAC version 7.5 database (http://www.gene-regulation.com/). 

# Some features may belong to more than one group. The argument minlen=25 implies 
# the minimum number of features in a gene set. If remain=TRUE, gene sets with less 
# than 25 members are grouped to the "remainder" group. "genesWurdinger" is an object
# containing gene names from the mRNA sequencing data set. 
# See help(matchGeneSets) for more detail information.
# Also see Vignette for more detail examples

## The "TFsym" is available on https://github.com/markvdwiel/GRridgeCodata
# gseTF <- matchGeneSets(genesWurdinger,TFsym,minlen=25,remain=TRUE)

}