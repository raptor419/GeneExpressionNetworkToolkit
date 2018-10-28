.launch_application <- function()
{
  shiny::runApp(appDir = system.file(package = GeneExpressionNetworkToolkit))
}
