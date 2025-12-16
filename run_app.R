run_app <- function() {
  shiny::runApp(
    system.file("app", package = "GVMR"),
    launch.browser = TRUE
  )
}