

#' Run shiny application in DEP2
#'
#' @export
#' @import shiny
#' @examples
run_app <- function() {
  # Launch the app
  appDir <- system.file("DEP2shiny221013", package = "DEP2")
  suppressWarnings(shiny::runApp(appDir, display.mode = "normal"))
}




