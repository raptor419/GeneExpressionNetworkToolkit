#' @title Launch Application
#' @export
launch_app <- function()
{
  input.good = paste(getwd(), "inst/GENT/GEOmetadb.sqlite", sep="/")
  if(!file_test("-f", input.good)){
    library(GEOmetadb)
    getSQLiteFile(destdir = getwd(), destfile = "GEOmetadb.sqlite.gz")
  }
  shiny::runApp(system.file('GENT', package = "GeneExpressionNetworkToolkit"))
}

