if(!require(pacman, quietly = TRUE)) install.packages("pacman")
pacman::p_load("devtools")

if(!require(diffusionMaps, quietly = TRUE)){
  cat("diffusionMap package is not installed.\n
       Downloaded the tar.gz file at https://github.com/RonBarry/diffusionMaps\n")
  
  # install diffusionMaps dependencies 
  pacman::p_load("RANN", "colorspace", "R.utils", "sqldf", "geometry", "plot3Drgl")
  
  system("git clone git@github.com:RonBarry/diffusionMaps.git")
  system("R CMD INSTALL diffusionMaps/diffusionMaps_2.0.0.tar.gz")
  
  unlink("diffusionMaps/", recursive=TRUE)
  cat("\ndiffusionMaps R package installed.\n")
}

if(!require(KrigLinCaution, quietly = TRUE)){
  devtools::install_github("jayverhoef/KrigLinCaution")
}

# removed from CRAN
if(!require(maptools, quietly = TRUE)){
  devtools::install_version("maptools", version="1.1-8", repos="https://cran.stat.unipd.it/")
}

# depends on maptools, removed from CRAN (London data are available with shp2graph)
if(!require(shp2graph, quietly = TRUE)){
  devtools::install_version("shp2graph", version="0-5", repos="https://cran.stat.unipd.it/" )
}

pacman::p_load("fdaPDE", "spatstat", "shp2graph", 
               "igraph", "spam", "GWmodel", "MetricGraph", "rSPDE")

library(diffusionMaps, quietly = TRUE)
library(KrigLinCaution, quietly = TRUE)
