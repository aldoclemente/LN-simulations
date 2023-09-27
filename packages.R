if(system.file(package = "diffusionMap") == ""){
  cat("diffusionMap package is not installed.\n
       Download the tar.gz file at https://github.com/RonBarry/diffusionMaps\n
        run from terminal: R CMD INSTALL diffusionMaps_2.0.0.tar.gz")
}

if(system.file(package = "KrigLinCaution")==""){
  devtools::install_github("jayverhoef/KrigLinCaution")
}

library(fdaPDE)
library(spatstat)
library(maptools)
library(shp2graph)
library(igraph)
library(rgeos)
library(spam)
library(GWmodel)
library(diffusionMaps)
library(KrigLinCaution)