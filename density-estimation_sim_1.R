# nonparametric density estimation - simulation 1
graphics.off()
rm(list=ls())

# Density Estimation over Linear Networks --------------------------------------

if(!require(pacman)) install.packages("pacman")

# loading packages and auxiliary functions
source("utils/packages.R")
source("utils/utils.R")
source("utils/plot.R")
source("utils/Simulation.R")
source("utils/create_knots.R")
source("utils/rank_reduced_kriging.R")
source("utils/aux.R")

data("simplenet")
delta = 0.025
mesh = as.mesh.1.5D(simplenet)
mesh = refine.mesh.1.5D(mesh,delta)

FEMbasis = create.FEM.basis(mesh)
nnodes = nrow(mesh$nodes)
spatstat.linnet = simplenet #as.linnet(mesh)

# Competing methods ------------------------------------------------------------
# methods[1] -> DE-PDE
# methods[2] -> KDE-HEAT  (available in spatstat package))
# methods[3] -> KDE-ES    (available in spatstat package, very slow !)
# methods[4] -> KDE-2D    (available in spatstat package)
# methods[5] -> VORONOI   (available in spatstat package, slow      !)
method_names = c("DE-PDE", "KDE-HEAT", "KDE-ES", "KDE-2D", "VORONOI")
methods = c(T,T,F,T,T)
method_names = method_names[methods]

# Test Hyperparameters ---------------------------------------------------------
n_obs = as.integer(c(50, 100, 150, 250)) # numbers of occurences
lambda = 10^seq(from=-4, to=-3,length.out = 20)
sources = c(6,8)         
n_sim = 30L
# test locations
test_locs = runiflpp(1000, spatstat.linnet)
test_locs = cbind(test_locs$data$x, test_locs$data$y)
set.seed(1234)

# True Density Function --------------------------------------------------------

DENSITY = linfun(aux, L=spatstat.linnet) 

Mass = fdaPDE:::CPP_get.FEM.Mass.Matrix(FEMbasis)
true_density = DENSITY(x=mesh$nodes[,1], y=mesh$nodes[,2])
true_density.FEM = true_density / sum( Mass %*% true_density)

# Storing density evaluated at test locations
true_density = DENSITY(x=test_locs[,1], y=test_locs[,2]) / sum( Mass %*% true_density)

# Building folders -------------------------------------------------------------
folder.name = "density-estimation/"

if(!dir.exists(folder.name)) {
  dir.create(folder.name)
}

folder.name = paste0(folder.name, "simulation_1/")
if(!dir.exists(folder.name)) {
  dir.create(folder.name)
}

# plot point pattern plot ------------------------------------------------------
PP = rlpp(n=n_obs[2], f = DENSITY)  
plot(mesh, linewidth=3, color="gray") + geom_point(data=data.frame(x=PP$data$x,y=PP$data$y),
                             aes(x=x, y=y), color="red3", size=4)
ggsave(paste0(folder.name, "point_pattern.pdf"),width = 7, height = 7)

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
    cat(paste("-------------------  ", i, " / ", n_sim,"  -------------------\n", sep="") )
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

    # ### KDE-HEAT ### -----------------------------------------------------------
    start = Sys.time()
    invisible(capture.output(bw <- bw.lppl(X = PP) ))
    invisible(capture.output(output_KDE_HEAT <- densityHeat(x = as.lpp(PP),
                                                               sigma = as.numeric(bw),
                                                               iterMax = 1e+8) ))

    coef_ <- as.linfun(output_KDE_HEAT/n_obs[j])(mesh$nodes[,1], mesh$nodes[,2])
    KDE_HEAT$update_estimate(estimate = FEM(coef_, FEMbasis),i,j)
    KDE_HEAT$update_error(true_density, test_locs=test_locs,i,j)

    cat(paste0("- KDE-HEAT DONE, time elapsed = ",
               difftime(Sys.time(),start, units = "mins")," mins \n"))

    # KDE-2D ------------- -----------------------------------------------------
    start = Sys.time()
    invisible(capture.output(bw <- bw.scott(X = PP) ))
    invisible(capture.output(output_KDE_2D <-  densityQuick.lpp(x = PP, sigma = bw)))

    cat(paste0("- KDE-2D DONE, time elapsed = ",
              difftime(Sys.time(),start, units = "mins")," mins \n"))
    coef_ <- as.linfun(output_KDE_2D/n_obs[j])(mesh$nodes[,1], mesh$nodes[,2])
    KDE_2D$update_estimate(estimate = FEM(coef_, FEMbasis),i,j)
    KDE_2D$update_error(true_density, test_locs=test_locs,i,j=j)

    ### VORONOI--------------------------- -------------------------------------
    start = Sys.time()
    invisible(capture.output(bw <- bw.voronoi(X = PP) ))
    invisible(capture.output(output_VORONOI <- densityVoronoi(X = PP, sigma = bw) ))

    cat(paste0("- VORONOI DONE, time elapsed = ",
              difftime(Sys.time(),start, units = "mins")," mins \n"))
    coef_ <- as.linfun(output_VORONOI/n_obs[j])(mesh$nodes[,1], mesh$nodes[,2])
    VORONOI$update_estimate(estimate = FEM(coef_, FEMbasis),i,j)
    VORONOI$update_error(true_density, test_locs=test_locs,i,j=j)
  }
  DE_PDE$compute_mean_field(j)
  KDE_HEAT$compute_mean_field(j)
  KDE_2D$compute_mean_field(j)
  VORONOI$compute_mean_field(j)
}                                     

save(folder.name, DE_PDE, KDE_HEAT, KDE_2D, VORONOI, folder.name,
     file = paste0(folder.name,"data",".RData"))

# Post processing --------------------------------------------------------------
SimulationBlock <- BlockSimulation(list(DE_PDE, KDE_HEAT, KDE_2D, VORONOI))
SimulationBlock$method_names

title.size <- 26
MyTheme <- theme(
  axis.text = element_text(size=title.size-5),
  axis.title = element_text(size=title.size),
  title = element_text(size=title.size),
  plot.title = element_text(hjust = 0.5),
  legend.text = element_text(size=title.size-5),
  legend.key.size = unit(1,"cm"),
  legend.key.height = unit(1,"cm"),
  legend.title = element_blank(),
  legend.background = element_rect(fill="white", color="black",
                                   linewidth =c(1,0.5))
)

ORDER = c(1,2,3,4)
plt <- boxplot(SimulationBlock, ORDER) +
  labs(title="RMSE", x="observations") +
  theme(legend.position.inside = c(0.90,0.80)) +
  MyTheme
ggsave(paste0(folder.name, "RMSE.pdf"), plot=plt, width = 7, height = 7) 

# setting same color scale
color.min <- rep(0., times = length(SimulationBlock$n_obs))
color.max <- rep(max(true_density), times = length(SimulationBlock$n_obs))

# No nel caso dei LN devi dare i coeffs.
#smooth_lim(FEM(true_density.FEM, FEMbasis))
colors = matrix(nrow=4, ncol=2)

colors[,1] = rep(0, times=4)
colors[,2] = rep(max(true_density.FEM), times = length(SimulationBlock$n_obs))

for(j in 1:length(SimulationBlock$n_obs)){
  for(i in 1:SimulationBlock$num_methods){
    colors[j,1] = min(min(SimulationBlock$Simulations[[i]]$meanField[[j]]$coeff), 
                      color.min[j])
    colors[j,2] = max(max(SimulationBlock$Simulations[[i]]$meanField[[j]]$coeff), 
                      color.max[j])
    
  }
}

for(i in 1:SimulationBlock$num_methods){
  for(j in 1:length(SimulationBlock$n_obs)){

    SimulationBlock$Simulations[[i]]$plot_mean_field(j,linewidth=3)+
            viridis::scale_color_viridis(limits=colors[i,]) + 
            theme( legend.position = "none")
    ggsave(paste0(folder.name,"estimate_", 
                  SimulationBlock$Simulations[[i]]$method_name,"_",
                  SimulationBlock$n_obs[j],".pdf"),width = 7, height = 7)
    
  }
}

for(i in 1:length(SimulationBlock$n_obs)){
plot(FEM(true_density.FEM, FEMbasis), linewidth=3) +
        scale_color_viridis(limits=colors[i,]) +
        theme( legend.position = "none")
ggsave(paste0(folder.name, "true_field_", n_obs[i], ".pdf"),width = 7, height = 7)

plot.colorbar(FEM(true_density.FEM, FEMbasis), 
              colorscale =  viridis, limits = colors[i,],
              width = 2,
              file = paste0(folder.name,"colorbar_", n_obs[i]))

}

# table ------------------------------------------------------------------------
# DE-PDE
cat("--- DE-PDEs ---\n")
rmse_table_de_pde <- matrix(DE_PDE$errors, nrow=SimulationBlock$num_sim, 
                     ncol=length(SimulationBlock$n_obs))
colMeans(rmse_table_de_pde)
apply(rmse_table_de_pde, MARGIN = 2, sd)

# KDE HEAT
cat("--- KDE HEAT ---\n")
rmse_table_heat <- matrix(KDE_HEAT$errors, nrow=SimulationBlock$num_sim, 
                     ncol=length(SimulationBlock$n_obs))
colMeans(rmse_table_heat)
apply(rmse_table_heat, MARGIN = 2, sd)

# KDE 2D
cat("--- KDE 2D ---\n")
rmse_table_2d <- matrix(KDE_2D$errors, nrow=SimulationBlock$num_sim, 
                     ncol=length(SimulationBlock$n_obs))
colMeans(rmse_table_2d)
apply(rmse_table_2d, MARGIN = 2, sd)

# VORONOI
cat("--- VORONOI ---\n")
rmse_table_voro <- matrix(VORONOI$errors, nrow=SimulationBlock$num_sim, 
                     ncol=length(SimulationBlock$n_obs))
colMeans(rmse_table_voro)
apply(rmse_table_voro, MARGIN = 2, sd)

tmp <- cbind(colMeans(rmse_table_de_pde), apply(rmse_table_de_pde, MARGIN = 2, sd),
      colMeans(rmse_table_heat), apply(rmse_table_heat, MARGIN = 2, sd),
      colMeans(rmse_table_2d), apply(rmse_table_2d, MARGIN = 2, sd),
      colMeans(rmse_table_voro), apply(rmse_table_voro, MARGIN = 2, sd))

colnames(tmp) <- c("DEG-mean", "DEG-sd", "HEAT-mean", "HEAT-sd", 
                   "2d-mean", "2d-sd", "voronoi-mean", "voronoi-sd")

rownames(tmp) <- n_obs
write.table( round(tmp, digits = 4), paste0(folder.name, "table.txt"))





