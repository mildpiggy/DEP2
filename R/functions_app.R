
# The DEP2 shiny app version, to verify the log file
app_version = 0503

#' Run shiny application in DEP2
#'
#' Run the shiny application in DEP2
#' @export
#' @examples
#' library(DEP2)
#' DEP2::run_app()
run_app <- function() {
  check_depends <- check_shiny_depends()

  if(isTRUE(check_depends)){
    # Launch the app
    appDir <- system.file("DEP2shiny230202", package = "DEP2")
    suppressWarnings(shiny::runApp(appDir, display.mode = "normal"))
  }else{
    stop("Packages ",paste0(check_depends,collapse = ", ")," are required but not found")
  }
}




