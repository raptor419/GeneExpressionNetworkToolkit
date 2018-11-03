#' @title Launch Application
#' @export
launch_application <- function()
{
  shiny::runApp(system.file('GENT', package = "GeneExpressionNetworkToolkit"))
}
