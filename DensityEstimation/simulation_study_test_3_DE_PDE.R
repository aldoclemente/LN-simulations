graphics.off()
rm(list=ls())

# Density Estimation over Linear Networks --------------------------------------
setwd("DensityEstimation/")

# loading packages and auxiliary functions
source("../packages.R")
source("../utils.R")
source("../settings.R")
source("../plot.R")
source("../Simulation.R")

tests.names = c("test_1", "test_2", "test_3")
ntest = 3

# available domains
domains = c("estevan", "ontario",  "simplenet") 

# Competing methods (available in spatstat package) ----------------------------
# methods[1] -> DE-PDE
# methods[2] -> KDE-PDE
# methods[3] -> KDE-ES   (very slow !)
# methods[4] -> KDE-2D
# methods[5] -> VORONOI  (slow !)

method_names = c("DE-PDE", "KDE-HEAT", "KDE-ES", "KDE-2D", "VORONOI")
methods = c(T,T,F,T,F)
method_names = method_names[methods]

# Fixing domain ----------------------------------------------------------------
sett = setting(domains[ntest]) 
mesh = sett$mesh
FEMbasis = sett$FEMbasis
nnodes = sett$nnodes
spatstat.linnet = sett$spatstat.linnet

plot(mesh, linewidth=0.75)

# Test Hyperparameters ---------------------------------------------------------

if(ntest==1){
  n_obs = as.integer(c(50, 75, 100, 150))
  lambda = 10^seq(from=-5,to=0.,length.out=20)
  n_sim = 30L
  
}
if(ntest==2){
  FAMILY="poisson"
  
  l<-make.link("log")
  link<-l$linkfun
  inv.link<-l$linkinv
  n_obs = c(100, 250, 500, 1000)
  lambda = 10^seq(from=-5,to=0.,length.out=20)
  n_sim=30L
}
if(ntest==3){
  n_obs = c(50, 100, 150, 250) # numbers of occurences
  lambda = 10^seq(from=-4, to=-3,length.out = 20)
  sources = c(6,8)         
  n_sim = 1L
  # test locations
  test_locs = runiflpp(1000, spatstat.linnet)
  test_locs = cbind(test_locs$data$x, test_locs$data$y)
}
set.seed(1234)

# True Density Function --------------------------------------------------------

# aux_test defined in settings.R
auxiliary_test = aux_density

DENSITY = linfun(auxiliary_test, L=spatstat.linnet) 

Mass = fdaPDE:::CPP_get.FEM.Mass.Matrix(FEMbasis)
true_density = DENSITY(x=mesh$nodes[,1], y=mesh$nodes[,2])
true_density.FEM = true_density / sum( Mass %*% true_density)

# density plot 
plot(FEM(true_density.FEM,FEMbasis), linewidth=0.75) + 
  viridis::scale_color_viridis()

# Storing density evaluated at test locations
true_density = DENSITY(x=test_locs[,1], y=test_locs[,2]) / sum( Mass %*% true_density)

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

# DE-PDE -----------------------------------------------------------------------
DE_PDE <- DensityEstimationSimulation(method_name=method_names[1],
                                      n_obs = n_obs, n_sim = n_sim,
                                      FEMbasis = FEMbasis)  
# KDE-HEAT ---------------------------------------------------------------------
KDE_HEAT <- DensityEstimationSimulation(method_name=method_names[2],
                                   n_obs = n_obs, n_sim = n_sim,
                                   FEMbasis = FEMbasis)  
# KDE-2D -----------------------------------------------------------------------
KDE_2D <- DensityEstimationSimulation(method_name=method_names[3],
                                       n_obs = n_obs, n_sim = n_sim,
                                       FEMbasis = FEMbasis) 
# VORONOI ----------------------------------------------------------------------
VORONOI <- DensityEstimationSimulation(method_name=method_names[4],
                                       n_obs = n_obs, n_sim = n_sim,
                                       FEMbasis = FEMbasis)  
# Loop -------------------------------------------------------------------------
for(j in 1:length(n_obs)){ 
  cat(paste("-------------------  n = ", n_obs[j], "  -------------------\n", sep="") ) 
  for(i in 1:n_sim){
    cat(paste("-------------------  ", i, " / ", nsim,"  -------------------\n", sep="") )
    PP = rlpp(n=n_obs[j], f = DENSITY)  
    data = cbind(PP$data$x, PP$data$y)
    
    ### DE-PDE ### -------------------------------------------------------------
    start = Sys.time()
    invisible(capture.output(output_CPP <- DE.FEM(data = data, 
                              FEMbasis = FEMbasis,
                              lambda = lambda[1])))
    cat(paste("- DE-PDE DONE, time elapsed = ", 
              difftime(Sys.time(),start, units = "mins")," mins \n"))
    
    DE_PDE$update_estimate(estimate=FEM(exp(output_CPP$g), FEMbasis),i,j)
    DE_PDE$update_error(true_density, test_locs=test_locs,i,j) 
    DE_PDE$compute_mean_field(j)

    ### KDE-HEAT ### -----------------------------------------------------------
    start = Sys.time()
    invisible(capture.output(bw <- bw.lppl(X = PP) ))
    invisible(capture.output(output_KDE_HEAT <- densityHeat(x = as.lpp(PP), 
                                                               sigma = as.numeric(bw), 
                                                               iterMax = 1e+9) )) 
    
    coef_ <- as.linfun(output_KDE_HEAT/n_obs[j])(mesh$nodes[,1], mesh$nodes[,2])
    KDE_HEAT$update_estimate(estimate = FEM(coef_, FEMbasis),i,j)
    KDE_HEAT$update_error(true_density, test_locs=test_locs,i,j) 
    KDE_HEAT$compute_mean_field(j)
    
    cat(paste0("- KDE-HEAT DONE, time elapsed = ", 
               difftime(Sys.time(),start, units = "mins")," mins \n"))
    
    # KDE-2D ------------- -----------------------------------------------------
    start = Sys.time()
    invisible(capture.output(bw <- bw.scott(X = PP) ))
    invisible(capture.output(output_KDE_2D <-  densityQuick.lpp(X = PP, sigma = bw))) 
    
    cat(paste0("- KDE-2D DONE, time elapsed = ", 
              difftime(Sys.time(),start, units = "mins")," mins \n"))
    coef_ <- as.linfun(output_KDE_2D/n_obs[j])(mesh$nodes[,1], mesh$nodes[,2])
    KDE_2D$update_estimate(estimate = FEM(coef_, FEMbasis),i,j)
    KDE_2D$update_error(true_density, test_locs=test_locs,i,j=j) 
    KDE_2D$compute_mean_field(j)
    
    ### VORONOI--------------------------- -------------------------------------
    start = Sys.time()
    invisible(capture.output(bw <- bw.voronoi(X = PP) ))
    invisible(capture.output(output_VORONOI <- densityVoronoi(X = PP, sigma = bw) ))
    
    cat(paste0("- VORONOI DONE, time elapsed = ", 
              difftime(Sys.time(),start, units = "mins")," mins \n"))
    coef_ <- as.linfun(output_VORONOI/n_obs[j])(mesh$nodes[,1], mesh$nodes[,2])
    VORONOI$update_estimate(estimate = FEM(coef_, FEMbasis),i,j)
    VORONOI$update_error(true_density, test_locs=test_locs,i,j=j)  
    VORONOI$compute_mean_field(j)
    
  }
  DE_PDE$compute_mean_field(j)
  KDE_HEAT$compute_mean_field(j)
  KDE_2D$compute_mean_field(j)
  VORONOI$compute_mean_field(j)
}                                     

save(DE_PDE, KDE_HEAT, KDE_2D, VORONOI,
     file = paste0(folder.name,date_,".RData"))

SimulationBlock <- BlockSimulation(list(SR_PDE, GWR, Lattice, RR_Krig))
SimulationBlock$method_names
ORDER = c(1,3,2,4)
boxplot(SimulationBlock, ORDER)

SimulationBlock <- BlockSimulation(list(SR_PDE, GWR))
SimulationBlock$method_names
ORDER = c(1,2)
boxplot(SimulationBlock,ORDER)
