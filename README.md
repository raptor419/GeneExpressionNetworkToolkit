# GENT: An R Package working as a Gene expression network toolkit

Authors: Tavpritesh Sethi and Devishi Kesar

Please email all comments/questions to devishi15024 [AT] iiitd.ac.in

### Summary

This repository will host the development version of the package. It is also available on CRAN. It acts as a toolkit to help find the most important genes for a given disease after gene expression analysis.

The package currently includes functionality to:
* extract data from GDB database
* employ various statistical tools to manipulate data
* predict most important genes using algorithms like bayesian neural networks
* generate plots to compare models using different algorithms

### Installation Instructions
You can install the most recent development version using the devtools package.  First you have
to install devtools using the following code.  Note that you only have to do this once
```
if(!require(devtools)) install.packages("devtools")
```
Then you can load the package and use the function `install_github`

```
library(devtools)
install_github("SAFE-ICU/GeneExpressionNetworkToolkit",dependencies=TRUE)
```

Note that this will install all the packages suggested and required to run our package.  It may take a few minutes the first time, but this only needs to be done on the first use.  In the future you can update to the most recent development version using the same code.

### Getting Started
How to run code

Run server.R or ui.R to view the app on shiny server <br />
Make sure to keep all modueles in the same working directory as server/ui.R <br />

Instructions to run the code for first time is given in comments of server.R <br />
For running the code for the first time :- <br />
- See server.R - inital lines commented out which include - Downloading GEOsqlite.db in the system <br />
- Make sure GEOmetadb.sqlite is present in the working directory. Also update the GEOsqlite  <br />
- Install required libraries <br />
