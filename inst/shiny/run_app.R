run_sgvm <- function() {
  shiny::runApp(
    system.file("app", package = "sgvm"),
    launch.browser = TRUE
  )
}