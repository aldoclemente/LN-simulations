# London House Pricing Case Study  ---------------------------------------------

setwd("SpatialRegression/")

# loading packages and auxiliary functions
source("../packages.R")
source("../utils.R")
source("../settings.R")
source("../plot.R")
source("../Simulation.R")

data(LNNT) # LN.nt   # SpatialLinesDataFrame <- linear network
data(LNHP) # LN.prop # data

# convert LN.nt into mesh.1.5D
london_RN <- as.mesh.1.5D(LN.nt)
london_RN <- 
  nrow(london_RN$nodes)
nrow(london_RN$edges)

FEMbasis  = create.FEM.basis(london_RN)

locs = cbind(LN.prop$X, LN.prop$Y)
