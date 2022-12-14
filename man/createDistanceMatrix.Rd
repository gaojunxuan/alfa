% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/distMatrix.R
\name{createDistanceMatrix}
\alias{createDistanceMatrix}
\title{Create a pairwise distance matrix for the given set of strings}
\usage{
createDistanceMatrix(dnaSet, metric = "euclidean", k = 10)
}
\arguments{
\item{dnaSet}{an object of type DNAStringSet or AAStringSet}

\item{metric}{metric for computing the distance, for now only 'euclidean'
and 'standardizedEuclidean' are supported}

\item{k}{size of k-mers (10 by default for DNA strings)}
}
\value{
an length(dnaSet) by length(dnaSet) matrix containing the pairwise
distance of every pair of strings in the set with the diagonal entries set
to NA.
}
\description{
Although the parameter is named dnaSet, it can be any set of strings of
the type DNAStringSet and AAStringSet over any well-defined alphabet
}
\examples{
library("ape")
data("woodmouse")
woodmouse <- as.list(woodmouse)
woodmouse <- lapply(woodmouse, ape::as.character.DNAbin)
woodmouse <- lapply(woodmouse, paste0, collapse = "")
woodmouseSet <- Biostrings::DNAStringSet(unlist(woodmouse))
alfa::createDistanceMatrix(woodmouseSet, metric = "euclidean")

}
