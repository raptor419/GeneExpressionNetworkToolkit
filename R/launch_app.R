#' @title Launch Application
#' @export
launch_app <- function()
{
  input.good = paste(getwd(), "GEOmetadb.sqlite", sep="/")
  library(GEOmetadb)
  getSQLiteFile(destdir = getwd(), destfile = "GEOmetadb.sqlite.gz")
  shiny::runApp(system.file('GENT', package = "GeneExpressionNetworkToolkit"))
}

