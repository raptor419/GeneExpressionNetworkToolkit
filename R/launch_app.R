#' @title Launch Application
#' @export
launch_app <- function()
{
  input.good = paste(getwd(), "inst/GENT/GEOmetadb.sqlite", sep="/")
  m = paste(getwd(), "inst/GENT", sep="/")
  if(!file_test("-f", input.good)){
    library(GEOmetadb)
    dirs = paste(getwd(), "inst/GENT", sep="/")
    getSQLiteFile(destdir = dirs, destfile = "GEOmetadb.sqlite.gz")
  }
  setwd(m)
  shiny::runApp(system.file('GENT', package = "GeneExpressionNetworkToolkit"))
}

