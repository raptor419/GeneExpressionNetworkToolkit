#' @title Launch Application
#' @export
launch_application <- function()
{
  library(GEOmetadb)
  getSQLiteFile(destdir = getwd(), destfile = "GEOmetadb.sqlite.gz")
  shiny::runApp(system.file('GENT', package = "GeneExpressionNetworkToolkit"))
}
