
#' Launch Shiny App for alfa
#'
#' A function that launches the Shiny app for alfa
#' The purpose of this app is only to illustrate how a Shiny
#' app works. The code has been placed in \code{./inst/shiny-scripts}.
#' (Codes adapted from the runTestingPackage.R file for the TestingPackage package
#' (Silva, A. (2022).)
#'
#' @return No return value but open up a Shiny page.
#'
#' @examples
#' \dontrun{
#'
#' alfa::runAlfa()
#' }
#'
#' @references
#' Silva, A. (2022). Anjalisilva/TestingPackage: A Simple R Package Illustrating Components of an R Package: 2019-2022 BCB410H - Applied Bioinformatics, University of Toronto, Canada. GitHub. \href{https://github.com/anjalisilva/TestingPackage}{Link}
#'
#' @export
#'
#' @importFrom shiny runApp

runAlfa <- function() {
  appDir <- system.file("shiny-scripts",
                        package = "alfa")
  actionShiny <- shiny::runApp(appDir, display.mode = "normal")
  return(actionShiny)
}
# [END]
