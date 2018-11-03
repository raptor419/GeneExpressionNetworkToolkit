#' @title Launch Application
#' @export
launch_application <- function()
{
  input.good = "~/inst/GENT/GEOmetadb.sqlite"
  if(!file_test("-f", input.good)){
    library(GEOmetadb)
    dirs = paste(getwd(), "inst/GENT", sep="/")
    getSQLiteFile(destdir = dirs, destfile = "GEOmetadb.sqlite.gz")
  }
  shiny::runApp(system.file('GENT', package = "GeneExpressionNetworkToolkit"))
}

