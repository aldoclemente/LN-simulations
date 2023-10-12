if(!require(diffusionMaps)){
  cat("diffusionMap package is not installed.\n
       Download the tar.gz file at https://github.com/RonBarry/diffusionMaps\n
        run from terminal: R CMD INSTALL diffusionMaps_2.0.0.tar.gz")
}

if(!require(KrigLinCaution)){
  devtools::install_github("jayverhoef/KrigLinCaution")
}

if(!require(pacman)) install.packages("pacman")
pacman::p_load("fdaPDE", "spatstat", "maptools", "shp2graph", 
               "igraph", "rgeos", "spam", "GWmodel")

library(diffusionMaps)
library(KrigLinCaution)
