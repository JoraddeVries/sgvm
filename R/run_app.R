#' Run SGVM Shiny App
#' @export

run_sgvm <- function() {
  app_dir <- system.file("shiny", package = "SGVM")
  if (app_dir == "") stop("Could not find app directory. Try reinstalling the package.")
  shiny::runApp(app_dir, launch.browser = TRUE)
}