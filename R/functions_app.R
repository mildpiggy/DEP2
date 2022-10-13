run_app <- function() {
  # Launch the app
  appDir <- system.file("DEP2_shiny20220926", package = "DEP2")
  suppressWarnings(runApp(appDir, display.mode = "normal"))
}




