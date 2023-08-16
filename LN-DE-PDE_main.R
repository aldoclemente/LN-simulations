graphics.off()
rm(list=ls())

# Density Estimation over Linear Networks --------------------------------------

# loading packages and auxiliary functions
source("packages.R")
source("utils.R")
source("settings.R")

tests.names = c("test_1", "test_2")
ntest = 1

# available domains
domains = c("simplenet", "ontario") 

# Competing methods (available in spatstat package) ----------------------------
# methods[1] -> DE-PDE
# methods[2] -> KDE-PDE
# methods[3] -> KDE-ES   (very slow !)
# methods[4] -> KDE-2D
# methods[5] -> VORONOI  (slow !)

methods.names = c("DE-PDE", "KDE-PDE", "KDE-ES", "KDE-2D", "VORONOI")
methods = c(T,T,F,T,F)

# Fixing domain ----------------------------------------------------------------
sett = setting(domains[ntest]) 
mesh = sett$mesh
FEMbasis = sett$FEMbasis
nnodes = sett$nnodes
spatstat.linnet = sett$spatstat.linnet

# plot(mesh,pch="."); plot(mesh)
# plot(spatstat.linnet)

# test locations
locs.test = runiflpp(1000, spatstat.linnet)
locs.test = cbind(locs.test$data$x, locs.test$data$y)

# Test Hyperparameters ---------------------------------------------------------

# number of simulations
nsim =  20 

if(ntest==1){
  n = c(50, 100) #, 150, 250) # numbers of occurences
  sources = c(6,8)         
}
if(ntest==2){ 
  n = c(100, 250, 500, 1000)
  sources = c(32, 185, 400)
}
set.seed(1234)

# True Density Function --------------------------------------------------------

# aux_test is defined in settings.R
auxiliary_test = aux_test[[ntest]]

DENSITY = linfun(auxiliary_test, L=spatstat.linnet) 

# Nb. CPP_get.FEM.Mass.Matrix has to be exported !
Mass = CPP_get.FEM.Mass.Matrix(FEMbasis)
true.density = DENSITY(x=mesh$nodes[,1], y=mesh$nodes[,2])
true.density.FEM = true.density / sum( Mass %*% true.density)

# density plot 
plot(FEM(true.density.FEM,FEMbasis))

# Storing density evaluated at test locations
true.density = DENSITY(x=locs.test[,1], y=locs.test[,2]) / sum( Mass %*% true.density)

# Building folders -------------------------------------------------------------
date_ = gsub(":","_",gsub(" ","-",Sys.time()))
if(!dir.exists("data/")) {
  dir.create("data/")
}

if( !dir.exists(paste("data/", tests.names[ntest],"/",sep=""))){
  dir.create(paste("data/", tests.names[ntest],"/",sep=""))
}

folder.name = paste("data/", tests.names[ntest],"/",date_,"/",sep="")

if(!dir.exists(folder.name)) {
  dir.create(folder.name)
}

# Performing simulation ----------------------------------------------------
if(ntest==1)lambda = 10^seq(from=-4, to=-3,length.out = 20)
if(ntest==2)lambda = 10^seq(from=-6, to=-3,length.out = 20)

source("LN-DE-PDE_loop.R")

# Post-processing
source("LN-DE-PDE_post_processing.R")
