

#' Run shiny application in DEP2
#'
#' @export
#'
#' @examples
run_app <- function() {
  # Launch the app
  appDir <- system.file("DEP2shiny20221013", package = "DEP2")
  suppressWarnings(runApp(appDir, display.mode = "normal"))
}




