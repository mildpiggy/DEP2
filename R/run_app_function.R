#' DEP shiny apps
#'
#' \code{run_app} launches an interactive shiny app
#' for interactive differential enrichment/expression analysis
#' of proteomics data.
#'
#' @param app 'LFQ' or 'TMT', The name of the app.
#' @return Launches a browser with the shiny app
#' @examples
#' \dontrun{
#' # Run the app
#' run_app('LFQ')
#'
#' run_app('TMT')
#'
#' }
#' @export
run_app <- function(app) {
  assertthat::assert_that(is.character(app),
                          length(app) == 1)

  # Locate all the shiny apps that exist
  valid_apps <- list.files(system.file("shiny_apps", package = "DEP2"))

  valid_apps_msg <- paste0("Valid apps are: '",
                           paste(valid_apps, collapse = "', '"), "'")

  # Show error if an unvalid app-name is given
  if (!app %in% valid_apps) {
    stop("Please run `run_app()` with a valid app as argument\n",
         valid_apps_msg, call. = FALSE)
  }

  # Launch the app
  appDir <- system.file("shiny_apps", app, package = "DEP2")
  suppressWarnings(runApp(appDir, display.mode = "normal"))
}
