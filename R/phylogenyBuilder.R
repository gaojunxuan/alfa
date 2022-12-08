#' Construct phylogenetic tree using the neighbor-joining algorithm.
#'
#' This function currently uses the function `nj` of the package `ape` to
#' build the tree internally.
#'
#' @param distMat distance matrix for constructing the phylogenetic tree
#' @param ... additional parameters for the plotting function `plot`
#'
#' @return a phylogenetic tree object
#'
#' @examples
#' library("ape")
#' data("woodmouse")
#' woodmouse <- as.list(woodmouse)
#' woodmouse <- lapply(woodmouse, ape::as.character.DNAbin)
#' woodmouse <- lapply(woodmouse, paste0, collapse = "")
#' woodmouseSet <- Biostrings::DNAStringSet(unlist(woodmouse))
#' distMat <- alfa::createDistanceMatrix(woodmouseSet, metric = "euclidean")
#' alfa::neighborJoiningTree(distMat)
#'
#' @references
#' Paradis E, Schliep K (2019). “ape 5.0: an environment for modern
#' phylogenetics and evolutionary analyses in R.” Bioinformatics, 35, 526-528.
#'
#' @import ape
#' @export
neighborJoiningTree <- function(distMat, ...) {
  njTree <- ape::bionj(distMat)
  plot(njTree, ...)
  return(njTree)
}

#' Construct phylogenetic tree using the neighbor-joining algorithm.
#'
#' @param distMat distance matrix for constructing the phylogenetic tree
#' @param ... additional parameters for the plotting function `plot`
#'
#' @return a phylogenetic tree object
#'
#' @examples
#' library("ape")
#' data("woodmouse")
#' woodmouse <- as.list(woodmouse)
#' woodmouse <- lapply(woodmouse, ape::as.character.DNAbin)
#' woodmouse <- lapply(woodmouse, paste0, collapse = "")
#' woodmouseSet <- Biostrings::DNAStringSet(unlist(woodmouse))
#' distMat <- alfa::createDistanceMatrix(woodmouseSet, metric = "euclidean")
#' alfa::upgmaTree(distMat)
#'
#' @references
#' Paradis E, Schliep K (2019). “ape 5.0: an environment for modern
#' phylogenetics and evolutionary analyses in R.” Bioinformatics, 35, 526-528.
#'
#' @import ape
#' @export
upgmaTree <- function(distMat, ...) {
  cluster <- ape::as.phylo(hclust(as.dist(distMat),
                                  method = "average",
                                  members = NULL))
  plot(cluster, ...)
  return(cluster)
}

#' Plot the pairwise distance in the distance matrix against
#' the pairwise distance on the given phylogenetic tree.
#'
#' This can be used to compare the effectiveness of a given alignment-free method
#' when `tree` is a phylogenetic tree created using alternative methods
#' such as MSA.
#'
#' @param distMat pairwise distance matrix
#' @param tree phylogenetic tree to compare against
#'
#' @examples
#' library("ape")
#' data("woodmouse")
#' woodmouse <- as.list(woodmouse)
#' woodmouse <- lapply(woodmouse, ape::as.character.DNAbin)
#' woodmouse <- lapply(woodmouse, paste0, collapse = "")
#' woodmouseSet <- Biostrings::DNAStringSet(unlist(woodmouse))
#' distMat <- alfa::createDistanceMatrix(woodmouseSet, metric = "euclidean")
#' njtree <- alfa::neighborJoiningTree(distMat)
#' plotDistanceCorrelation(distMat, njtree)
#'
#' @export
plotDistanceCorrelation <- function(distMat, tree) {
  x <- as.vector(as.dist(distMat))
  y <- as.vector(as.dist(cophenetic(tree)))
  plot(x, y)
  abline(lm(y ~ x), col = "red")
}
