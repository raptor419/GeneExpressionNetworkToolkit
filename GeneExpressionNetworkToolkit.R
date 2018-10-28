#' @title ShinyApp
#'
#' @description Launch the ShinyApp
#'
#' @param 
#'
#' @return NULL
#'
#' @examples launch_application()
#'
#' @export

launch_application <- function()
{
  shiny::runApp(appDir = system.file(package = [GeneExpressionNetworkToolkit]))
}
