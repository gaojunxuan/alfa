neighborJoiningTree <- function(distMat, ...) {
  njTree <- ape::bionj(distMat)
  plot(njTree, ...)
  return(njTree)
}

upgmaTree <- function(distMat, ...) {
  cluster <- ape::as.phylo(hclust(as.dist(distMat),
                                  method = "average",
                                  members = NULL))
  plot(cluster, ...)
  return(cluster)
}

#' Plot the pairwise distance in the distance matrix against
#' the pairwise distance on the given phylogenetic tree. This can
#' be used to compare the effectiveness of a given alignment-free method
#' when `tree` is a phylogenetic tree created using alternative methods
#' such as MSA.
#'
#' @export
plotDistanceCorrelation <- function(distMat, tree) {
  x <- as.vector(as.dist(distMat))
  y <- as.vector(as.dist(cophenetic(tree)))
  plot(x, y)
  abline(lm(y ~ x), col = "red")
}
