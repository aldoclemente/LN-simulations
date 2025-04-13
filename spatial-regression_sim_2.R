# nonparametric regression - simulation 2 - cluster data + matern error
graphics.off()
rm(list=ls())
set.seed(0)

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

inla.edges = list()
for(i in 1:length(simplenet$from)){
  inla.edges[[i]] = rbind(c(simplenet$vertices$x[simplenet$from[i]],
                            simplenet$vertices$y[simplenet$from[i]]),
                          c(simplenet$vertices$x[simplenet$to[i]],
                            simplenet$vertices$y[simplenet$to[i]]))
}

inla.graph <- metric_graph$new(edges = inla.edges)
inla.graph$build_mesh(h = delta)

mesh = as.mesh.1.5D(inla.graph)
FEMbasis = create.FEM.basis(mesh)
nnodes = nrow(mesh$nodes)
spatstat.linnet = as.linnet(mesh)

# NO RR-Krig 
method_names = c("NRG", "GWR", "Lattice", "WMG")#, "IsoExp")
# Test Hyperparameters ---------------------------------------------------------

n_obs = as.integer(c(100, 150, 250, 500)) 
lambda = 10^seq(from=-8, to=-6,length.out = 250)
sources = c(6,7)         
n_sim = 30L

## Building true signal --------------------------------------------------------
source("utils/aux.R")
aux2 = linfun(cluster, L=spatstat.linnet) 
plot(FEM(aux2(mesh$nodes[,1], mesh$nodes[,2]), FEMbasis), linewidth=3) + scale_color_viridis()

# {
# centers = matrix(nrow=0,ncol=2)
# for(i in 1:length(simplenet$from)){
#   centers = rbind(centers, c(0,0))
#   centers[i,1] = mean(c(simplenet$lines$ends[i,1], simplenet$lines$ends[i,3]))
#   centers[i,2] = mean(c(simplenet$lines$ends[i,2], simplenet$lines$ends[i,4]))
# }
# 
# lengths = matrix(0,nrow=nrow(simplenet$lines$ends), ncol=1)
# for( i in 1:nrow(lengths)){
#   lengths[i] = crossdist(simplenet$lines$ends[i,1], simplenet$lines$ends[i,2],
#                          simplenet$lines$ends[i,3], simplenet$lines$ends[i,4])
# }
# }

PP = rlpp(n=150, f = aux2)
plot(PP)
# in utils.R 
ND = compute_dist_matrix(points1= mesh$nodes, 
                         points2 = mesh$nodes, 
                         L =spatstat.linnet)

# aux is defined in settings.R
true_signal <- aux(mesh$nodes[,1], mesh$nodes[,2])
#plot(FEM(true_signal, FEMbasis), linewidth=0.75)  + scale_color_viridis()

observations <- true_signal + rnorm(nnodes, mean=0, sd = 0.05*diff(range(true_signal)))

test_locations = mesh$nodes
test_true = aux(test_locations[,1], test_locations[,2])

# ---

sigma=0.1*abs(diff(range(true_signal))); range=0.3;
nu=1.5; rspde.order=2;
op = matern.operators(nu=nu, range=range, sigma=sigma, parameterization = "matern", m=rspde.order,
                      graph = inla.graph)
#u = simulate(op)
#range(u)

# ------------------------------------------------------------------------------
folder.name = "spatial-regression/"
if(!dir.exists(folder.name)) {
  dir.create(folder.name)
}

folder.name = paste0(folder.name, "simulation_2/")
if(!dir.exists(folder.name)) {
  dir.create(folder.name)
}

# NRG -----------------------------------------------------------------------
SR_PDE <- SpatialRegressionSimulation(method_name=method_names[1],
                                      n_obs = n_obs, n_sim = n_sim,
                                      FEMbasis = FEMbasis)  
# GWR -- -----------------------------------------------------------------------
GWR <- SpatialRegressionSimulation(method_name=method_names[2],
                                   n_obs = n_obs, n_sim = n_sim,
                                   FEMbasis = FEMbasis)  
# Lattice method ---------------------------------------------------------------
tmp = as.lattice(mesh)
nodes.lattice = tmp$nodes.lattice
adj_matrix = tmp$adj_matrix
T_matrix = makeTranMatrix(adj_matrix, M = 0.5)

Lattice <- SpatialRegressionSimulation(method_name=method_names[3],
                                       n_obs = n_obs, n_sim = n_sim,
                                       FEMbasis = FEMbasis) 

# Whittle Matern fields on graphs (INLA) ---------------------------------------
WMG <- SpatialRegressionSimulation(method_name=method_names[4],
                                   n_obs = n_obs, n_sim = n_sim,
                                   FEMbasis = FEMbasis)

# Isotropic Exponential Cov Function (Andares 2020 AoS)  -----------------------
# IsoExp <- SpatialRegressionSimulation(method_name=method_names[5],
#                                       n_obs = n_obs, n_sim = n_sim,
#                                       FEMbasis = FEMbasis)
#inla.graph$check_euclidean()
inla.predict <- data.frame(edge_number = inla.graph$mesh$VtE[,1],
                            distance_on_edge = inla.graph$mesh$VtE[,2])


locations = rep(list(), times=n_sim*length(n_obs))
OBSERVATIONS = rep(list(), times=n_sim * length(n_obs))
# Loop -------------------------------------------------------------------------
pdf(paste0(folder.name, "gcv.pdf"))
for(j in 1:length(n_obs)){
  cat(paste("-------------------  n = ", n_obs[j], "  -------------------\n", sep="") )
  for(i in 1:n_sim){
    cat(paste("-------------------  ", i, " / ", n_sim,"  -------------------\n", sep="") )
    
    PP = rlpp(n=n_obs[j], f = aux2) #runiflpp(n=n_obs[j], L=spatstat.linnet)
    locs = cbind(PP$data$x, PP$data$y)
    net_dist = compute_dist_matrix(points1 = locs, points2 = locs, L = spatstat.linnet)
    
    # obs = aux(locs[,1], locs[,2]) + rnorm(n=n_obs[j], mean=0, sd=0.05*abs(diff(range(true_signal))))
    
    u = simulate(op)
    err = eval.FEM(FEM(u, FEMbasis), locations = locs)
    obs = aux(locs[,1], locs[,2]) + err
    
    OBSERVATIONS[[(j-1)*n_sim + i]] = obs
    locations[[(j-1)*n_sim + i]] = locs
    
    ### NRG ### ------------------------------------------------------------- 
    invisible(capture.output(output_CPP <- smooth.FEM(observations = obs, 
                                                      locations = locs,
                                                      FEMbasis = FEMbasis,
                                                      lambda = lambda,
                                                      lambda.selection.criterion = "grid",
                                                      lambda.selection.lossfunction = "GCV",
                                                      DOF.evaluation = "stochastic"))) # "stochastic"
    
    print(plot(lambda, output_CPP$optimization$GCV_vector, type="l", lwd=2,
               xlab=expression(lambda), ylab="", main=paste0("obs: ", n_obs[j], ", sim: ", i)))
    
    y_hat = eval.FEM(output_CPP$fit.FEM, locations = locs)
    SR_PDE$update_estimate(estimate = output_CPP$fit.FEM,i = i, j=j)
    SR_PDE$update_y_hat(vec =y_hat, i = i, j = j)
    
    SR_PDE$update_error(y_hat=eval.FEM(output_CPP$fit.FEM, test_locations), y_true=test_true,i,j)
    
    ### GWR ### ----------------------------------------------------------------
    
    data_ = data.frame(y = obs)
    Sp.data = SpatialPointsDataFrame(coords = locs,
                                     data = data_)
    
    # ND
    invisible(capture.output(bw.ND <- bw.gwr(y ~ 1,
                                             data = Sp.data, 
                                             approach="AIC", 
                                             kernel="gaussian",
                                             dMat = net_dist)))
    
    predict_net_dist = compute_dist_matrix(locs, test_locations, spatstat.linnet)
    Sp.predict = SpatialPointsDataFrame(coords=test_locations, 
                                        data = data.frame(matrix(NA, nrow=nrow(test_locations),ncol=1)))
    invisible(capture.output(GWR.ND <- gwr.predict(y ~ 1, 
                                                   data = Sp.data,
                                                   predictdata = Sp.predict ,
                                                   kernel = "gaussian",
                                                   bw = bw.ND,
                                                   dMat1 = predict_net_dist, 
                                                   dMat2 = net_dist)))
    
    GWR$update_error(y_hat=GWR.ND$SDF$prediction, y_true=test_true,i,j)
    GWR$update_estimate(estimate = FEM(GWR.ND$SDF$prediction, FEMbasis), i = i, j=j)
    
    
    # lattice based method -----------------------------------------------------
    # dummy z-coordinates !
    locations.lattice = cbind(locs, rep(0, times=n_obs[j]))
    output_press = crossvalNparSmoother(T=T_matrix,
                                        nodelocs = nodes.lattice,
                                        locs = locations.lattice,
                                        Z = obs,
                                        k_max=1500)
    
    output_lattice = nparSmoother(T=T_matrix,
                                  nodelocs = nodes.lattice,
                                  locs = locations.lattice,
                                  Z = obs,
                                  k = output_press$k)
    prediction.latt = eval.FEM(FEM(output_lattice[,4], FEMbasis), locations=locs)
    Lattice$update_estimate(estimate = FEM(output_lattice[,4], FEMbasis), i = i, j=j)
    
    Lattice$update_error(y_hat=eval.FEM(FEM(output_lattice[,4], FEMbasis), test_locations), 
                         y_true=test_true,i,j)
    
    Lattice$update_y_hat(vec = prediction.latt, i=i, j=j)

    # Whittle Matern fields (INLA) ---------------------------------------------
    inla.graph$add_observations(Sp.data)
    
    fit <- graph_lme(y ~ 1, graph = inla.graph, model = "WM")
  
    inla.coeffs <- predict(fit, newdata= inla.predict, normalized = TRUE)
    
    y_hat = eval.FEM(FEM(inla.coeffs$mean, FEMbasis), locations = locs)
    WMG$update_estimate(estimate = FEM(inla.coeffs$mean, FEMbasis),i = i, j=j)
    WMG$update_y_hat(vec =y_hat, i = i, j = j)
    WMG$update_error(y_hat = eval.FEM(FEM(inla.coeffs$mean, FEMbasis), locations = test_locations),
                     y_true=test_true,i,j)
    inla.graph$clear_observations()
    
    # Iso Exp (Andares (2020) AoS ) --------------------------------------------
    #inla.graph$add_observations(Sp.data)
    # fit <- graph_lme(y ~ 1, graph = inla.graph, model = "isoExp")
    # 
    # inla.coeffs <- predict(fit, newdata= inla.predict, normalized = TRUE)
    # 
    # y_hat = eval.FEM(FEM(inla.coeffs$mean, FEMbasis), locations = locs)
    # IsoExp$update_estimate(estimate = FEM(inla.coeffs$mean, FEMbasis),i = i, j=j)
    # IsoExp$update_y_hat(vec =y_hat, i = i, j = j)
    # IsoExp$update_error(y_hat = eval.FEM(FEM(inla.coeffs$mean, FEMbasis), locations = test_locations),
    #                     y_true=test_true,i,j)
    # inla.graph$clear_observations()
  }
  SR_PDE$compute_mean_field(j)
  GWR$compute_mean_field(j)
  Lattice$compute_mean_field(j)
  WMG$compute_mean_field(j)
  #IsoExp$compute_mean_field(j)
}                                     
dev.off()

save(SR_PDE, GWR, Lattice, WMG, locations, OBSERVATIONS,
     folder.name,
     file = paste0(folder.name,"data",".RData"))

# Post processing --------------------------------------------------------------
SimulationBlock <- BlockSimulation(list(SR_PDE, GWR, Lattice, WMG))#, IsoExp))

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
SimulationBlock$method_names
ORDER = c(1,4,2,3) #, 5)

{
  plt <- boxplot(SimulationBlock, ORDER) +
    labs(title="RMSE", x="observations") +
    theme(legend.position = "inside",
          legend.position.inside = c(0.80,0.80)) +
    MyTheme
  ggsave(paste0(folder.name, "RMSE.pdf"), plot=plt, width = 8, height = 7) 
}

{
  plt <- boxplot(SimulationBlock, ORDER) +
    labs(title="RMSE", x="observations") + ylim(c(0,0.02)) +
    theme(legend.position = "none")  +
    MyTheme
  ggsave(paste0(folder.name, "RMSE-no-legend.pdf"), plot=plt, width = 8, height = 7) 
}

# ---
SimulationBlock <- BlockSimulation(list(SR_PDE, GWR, Lattice, WMG))
{
  plt <- plot(mesh, linewidth=0.75)
  ggsave(paste0(folder.name, "domain.pdf"), plot=plt,width = 7, height = 7)
}

# setting same color scale
color.min <- rep(min(true_signal), times = length(SimulationBlock$n_obs))
color.max <- rep(max(true_signal), times = length(SimulationBlock$n_obs))

for(j in 1:length(SimulationBlock$n_obs)){
  for(i in 1:SimulationBlock$num_methods){
    color.min[j] <- min(min(SimulationBlock$Simulations[[i]]$meanField[[j]]$coeff), 
                        color.min[j])
    color.max[j] <- max(max(SimulationBlock$Simulations[[i]]$meanField[[j]]$coeff), 
                        color.max[j])
  }
}
#cbind(color.min, color.max)

for(i in 1:SimulationBlock$num_methods){
  for(j in 1:length(SimulationBlock$n_obs)){
    plt <- SimulationBlock$Simulations[[i]]$plot_mean_field(j,linewidth=3)+
      viridis::scale_color_viridis(limits=c(color.min[j],color.max[j])) + 
      theme( legend.position = "none")
    ggsave(paste0(folder.name,"estimate_", SimulationBlock$Simulations[[i]]$method_name,"_",
                  SimulationBlock$n_obs[j],".pdf"), plot=plt,width = 7, height = 7)
  }
}

for(j in 1:length(SimulationBlock$n_obs)){
  print(plot.colorbar(FEM(aux(mesh$nodes[,1], mesh$nodes[,2]), FEMbasis), 
                      colorscale =  viridis, limits = c(color.min[j], color.max[j]),
                      width = 2,
                      file = paste0(folder.name,paste0("colorbar_", n_obs[j]))))
}


for(j in 1:length(SimulationBlock$n_obs)){
  plot(FEM(aux(mesh$nodes[,1], mesh$nodes[,2]), FEMbasis), linewidth=3) +
    scale_color_viridis(limits=c(color.min[j],color.max[j])) +
    theme( legend.position = "none")
  ggsave(paste0(folder.name, "true_field_", n_obs[j], ".pdf"),width = 7, height = 7)
}

#n_obs = SimulationBlock$n_obs
for(j in 1:length(n_obs)){
  locs = locations[[(j-1)*n_sim + 1]]
  obs = OBSERVATIONS[[(j-1)*n_sim + 1]]
  
  plot(mesh, linewidth=3, color="gray") + geom_point(data=data.frame(x=locs[,1],y=locs[,2]),
                                                     aes(x=x, y=y, color=obs), size=5) + 
    scale_color_viridis(limits=c(min(color.min[j], min(obs)), max(color.max[j], max(obs)))) + 
    theme( legend.position = "none")
  ggsave(paste0(folder.name, "observations_", n_obs[j], ".pdf"),width = 7, height = 7)
}

# table ------------------------------------------------------------------------
# NRG
cat("--- NRG ---\n")
rmse_table_sr_pde <- matrix(SR_PDE$errors, nrow=SimulationBlock$num_sim, 
                            ncol=length(SimulationBlock$n_obs))
colMeans(rmse_table_sr_pde)
apply(rmse_table_sr_pde, MARGIN = 2, sd)

# GWR
cat("--- GWR ---\n")
rmse_table_gwr <- matrix(GWR$errors, nrow=SimulationBlock$num_sim, 
                         ncol=length(SimulationBlock$n_obs))
colMeans(rmse_table_gwr)
apply(rmse_table_gwr, MARGIN = 2, sd)

# Lattice
cat("--- Lattice ---\n")
rmse_table_lattice <- matrix(Lattice$errors, nrow=SimulationBlock$num_sim, 
                             ncol=length(SimulationBlock$n_obs))
colMeans(rmse_table_lattice)
apply(rmse_table_lattice, MARGIN = 2, sd)

# RR Krig
# cat("--- RR Krig ---\n")
# rmse_table_krig <- matrix(RR_Krig$errors, nrow=SimulationBlock$num_sim, 
#                           ncol=length(SimulationBlock$n_obs))
# colMeans(rmse_table_krig)
# apply(rmse_table_krig, MARGIN = 2, sd)

# WMG
cat("--- WMG ---\n")
rmse_table_wmg <- matrix(WMG$errors, nrow=SimulationBlock$num_sim, 
                         ncol=length(SimulationBlock$n_obs))
colMeans(rmse_table_wmg)
apply(rmse_table_wmg, MARGIN = 2, sd)

tmp <- cbind(colMeans(rmse_table_sr_pde), apply(rmse_table_sr_pde, MARGIN = 2, sd),
             colMeans(rmse_table_gwr), apply(rmse_table_gwr, MARGIN = 2, sd),
             colMeans(rmse_table_lattice), apply(rmse_table_lattice, MARGIN = 2, sd),
             colMeans(rmse_table_wmg), apply(rmse_table_wmg, MARGIN = 2, sd))

colnames(tmp) <- c("NRG-mean", "NRG-sd", "GWR-mean", "GWR-sd", 
                   "Lattice-mean", "Lattice-sd",
                   "WMG-mean", "WMG-sd")
rownames(tmp) <- n_obs

write.table( round(tmp, digits = 4), paste0(folder.name, "table-mean-sd.txt"))

# ---

tmp <- cbind(apply(rmse_table_sr_pde, MARGIN=2, median), apply(rmse_table_sr_pde, MARGIN = 2, IQR),
             apply(rmse_table_gwr, MARGIN = 2, median), apply(rmse_table_gwr, MARGIN = 2, IQR),
             apply(rmse_table_lattice, MARGIN = 2, median), apply(rmse_table_lattice, MARGIN = 2, IQR),
             apply(rmse_table_wmg, MARGIN = 2, median), apply(rmse_table_wmg, MARGIN = 2, IQR))

colnames(tmp) <- c("NRG-Q2", "NRG-IQR", "GWR-Q2", "GWR-IQR", 
                   "Lattice-Q2", "Lattice-IQR",
                   "WMG-Q2", "WMG-IQR")
rownames(tmp) <- n_obs

write.table( round(tmp, digits = 4), paste0(folder.name, "table-median-IQR.txt"))
