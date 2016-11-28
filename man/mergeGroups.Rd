\name{mergeGroups}
\alias{mergeGroups}

\title{
Merge groups in a partition 
}

\description{
Pathway-based partition often contains a considerable number of gene sets (or groups). 
This function merges groups resulted by \code{\link{matchGeneSets}}.
The first principal component in each group is calculated.
Hierarchical clustering analysis is then performed on the first principal components from all groups.
Important Note: re-grouping is only done in the non-reminder group.
}

\usage{mergeGroups(highdimdata, initGroups=initGroups,maxGroups=maxGroups,
                    methodDistance="manhattan", methodClust="complete")
}
\arguments{

  \item{highdimdata}{
  Matrix or numerical data frame. Contains the primary data of the study. Columns are samples, rows are features.
}

  \item{initGroups}{
A list of initial given groups resulted from \code{\link{matchGeneSets}}.
}

  \item{maxGroups}{
Numeric. The new desired number of groups.
}

  \item{methodDistance}{
The distance method used for clustering. See \code{\link{dist}} for futher options. Default: manhattan distance.
}

\item{methodClust}{
The agglomeration method used for Grouping. See \code{\link{hclust}} for further options. Default: complete.
}
}


\value{
A list object containing:
 \item{newGroups}{A list the components of which contain the indices of the features belonging to each of the group.
This object is the same as the object created by \code{\link{matchGeneSets}} and \code{\link{CreatePartition}}}
  \item{newGroupMembers}{A list of members in the new merged groups.}
}

\author{
Putri W. Novianti
}

\seealso{
Creating partitions: \code{\link{CreatePartition}}.
Creating partitions based on overlaping groups (gene sets/pathways): \code{\link{matchGeneSets}}.

}

\examples{
# Source: http://software.broadinstitute.org/gsea/msigdb/collections.jsp
# A GMT file containing information about transcription factor targets should 
# be downloaded first from the aforementioned source.
# Section C3: motif gene sets; subsection: transcription factor targets;
# file: "c3.tft.v5.0.symbols.gmt"
# Details of the gene sets:
# Gene sets contain genes that share a transcription factor binding site 
# defined in the TRANSFAC (version 7.4, http://www.gene-regulation.com/) database.
# Each of these gene sets is annotated by a TRANSFAC record.

# Load data objects
data(dataWurdinger)

# Transform the data set to the square root scale
dataSqrtWurdinger <- sqrt(datWurdinger_BC)

#Standardize the transformed data
datStdWurdinger <- t(apply(dataSqrtWurdinger,1,function(x){(x-mean(x))/sd(x)}))

# A list of gene names in the primary RNAseq data
genesWurdinger <- annotationWurdinger$geneSymbol

## Creating partitions based on pathways information (e.g. GSEA object)
## Some variables may belong to more than one groups (gene sets).
## The argument minlen=25 implies the minimum number of members in a gene set
## If remain=TRUE, gene sets with less than 25 members are grouped to the
## "remainder" group.
## The "TFsym" is available on https://github.com/markvdwiel/GRridgeCodata
# gseTF <- matchGeneSets(genesWurdinger,TFsym,minlen=25,remain=TRUE)

## Regrouping gene sets by hierarchical clustering analysis.
## The number of gene sets from the GSEA database is relatively too high to be used 
## in the GRridge model. Here, the initial gene sets are re-grouped into maxGroups=5, using
## information from the primary data set.
# gseTF_newGroups <- mergeGroups(highdimdata=datStdWurdinger, initGroups =gseTF, maxGroups=5);

## Extracting indices of new groups
## This following object (gseTF2) can be used further as an input
## in the "partitions" argument in the "grridge" function
# gseTF2 <- gseTF_newGroups$newGroups

## Members of the new groups
# newGroupMembers <- gseTF_newGroups$newGroupMembers
}
