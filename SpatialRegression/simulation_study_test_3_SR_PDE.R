graphics.off()
rm(list=ls())

# Spatial Regression over Linear Networks --------------------------------------

if(!require(pacman)) install.packages("pacman")
pacman::p_load("rstudioapi")

# setting working directory 
setwd(dirname(getActiveDocumentContext()$path))

# loading packages and auxiliary functions
source("../packages.R")
source("../settings.R")
source("../utils/utils.R")
source("../utils/plot.R")
source("../utils/Simulation.R")
source("../utils/CaseStudy.R")
source("create_knots.R")
source("rank_reduced_kriging.R")

tests.names = c("test_1", "test_2", "test_3")
ntest = 3

# available domains
domains = c("estevan", "ontario",  "simplenet") 

# Competing methods (available in spatstat package) ----------------------------

method_names = c("SR-PDE", "GWR", "Lattice", "RR-Krig")
methods = c(T,T,T,T)

# Fixing domain ----------------------------------------------------------------
sett = setting(domains[ntest]) 
mesh = sett$mesh

mesh = fdaPDE::refine.mesh.1.5D(mesh, delta=0.0025)
FEMbasis = create.FEM.basis(mesh)
nnodes = nrow(mesh$nodes)
spatstat.linnet = as.linnet(mesh) # sett$spatstat.linnet

plot(mesh, linewidth=0.75)
plot(spatstat.linnet)

# Test Hyperparameters ---------------------------------------------------------

if(ntest==1){
  n_obs = as.integer(c(25, 50, 75, 100))
  lambda = 10^seq(from=-5,to=0.,length.out=20)
  n_sim = 30L
  
}
if(ntest==2){
  FAMILY="poisson"
  l<-make.link("log")
  link<-l$linkfun
  inv.link<-l$linkinv
  
  n_obs = as.integer(c(100, 250, 500, 1000))
  lambda = 10^seq(from=-5,to=0.,length.out=20)
  n_sim=30L
}
if(ntest==3){
  n_obs = as.integer(c(25, 50, 75, 100)) 
  lambda = 10^seq(from=-4, to=-3,length.out = 20)
  sources = c(6,8)         
  n_sim = 30L
  # test locations
  #locs.test = runiflpp(1000, spatstat.linnet)
  #locs.test = cbind(locs.test$data$x, locs.test$data$y)
}
set.seed(1234)

## Building true signal --------------------------------------------------------

# in ../utils.R 
ND = compute_dist_matrix(points1= mesh$nodes, 
                         points2 = mesh$nodes, 
                         L =spatstat.linnet)

# aux_density is defined in settings.R
true_signal <- aux_density(mesh$nodes[,1], mesh$nodes[,2])
plot(FEM(true_signal, FEMbasis), linewidth=0.75)  + scale_color_viridis()

observations <- true_signal + rnorm(nnodes, mean=0, sd = 0.05*diff(range(true_signal)))

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

# SR-PDE -----------------------------------------------------------------------
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
# Reduced Rank Kriging ---------------------------------------------------------
RR_Krig <- SpatialRegressionSimulation(method_name=method_names[4],
                                       n_obs = n_obs, n_sim = n_sim,
                                       FEMbasis = FEMbasis)  


locations = rep(list(), times=n_sim*length(n_obs))
OBSERVATIONS = rep(list(), times=n_sim * length(n_obs))
# Loop -------------------------------------------------------------------------
for(j in 1:length(n_obs)){
  cat(paste("-------------------  n = ", n_obs[j], "  -------------------\n", sep="") )
  for(i in 1:n_sim){
    cat(paste("-------------------  ", i, " / ", n_sim,"  -------------------\n", sep="") )
    
    sample_ = sample(1:nnodes, size=n_obs[j])
    
    locs = mesh$nodes[sample_,]
    obs = observations[sample_]
    OBSERVATIONS[[(j-1)*n_sim + i]] = obs
    locations[[(j-1)*n_sim + i]] = locs
    
    net_dist = ND[sample_, sample_]
    
    obs = true_signal[sample_] + rnorm(n_obs[j], mean=0, sd = 0.05*diff(range(true_signal)))
    OBSERVATIONS[[j]] = OBSERVATIONS[[j]] + obs/n_sim
    ### SR-PDE ### ------------------------------------------------------------- 
    invisible(capture.output(output_CPP <- smooth.FEM(observations = obs, 
                                                      locations = locs,
                                                      FEMbasis = FEMbasis,
                                                      lambda = lambda,
                                                      lambda.selection.criterion = "grid",
                                                      lambda.selection.lossfunction = "GCV",
                                                      DOF.evaluation = "stochastic"))) # "stochastic"
    
    y_hat = eval.FEM(output_CPP$fit.FEM, locations = locs)
    SR_PDE$update_estimate(estimate = output_CPP$fit.FEM,i = i, j=j)
    SR_PDE$update_error(y_hat=y_hat , y_true=true_signal[sample_] , i=i, j=j)
    SR_PDE$update_y_hat(vec =y_hat, i = i, j = j)
    ### GWR ### ----------------------------------------------------------------
    
    data_ = data.frame(observations = obs)
    Sp.data = SpatialPointsDataFrame(coords = locs,
                                     data = data_)
    
    # ND
    invisible(capture.output(bw.ND <- bw.gwr(observations ~ 1,
                                             data = Sp.data, 
                                             approach="AIC", 
                                             kernel="gaussian",
                                             dMat = net_dist)))
    
    invisible(capture.output(GWR.ND <- gwr.basic(observations ~ 1, 
                                                 data = Sp.data, 
                                                 kernel = "gaussian",
                                                 bw = bw.ND,
                                                 dMat = net_dist)))
    
    GWR$update_error(y_hat = GWR.ND$SDF$yhat,
                     y_true = true_signal[sample_],i,j)
    GWR$update_y_hat(vec = GWR.ND$SDF$yhat, i=i, j=j)
    
    predict_dist = ND[sample_, ]
    
    Sp.predict = SpatialPointsDataFrame(coords=mesh$nodes, 
                                        data = data.frame(matrix(NA, nrow=nrow(mesh$nodes),ncol=1)))
    invisible(capture.output(GWR.ND <- gwr.predict(observations ~ 1, 
                          data = Sp.data,
                          predictdata = Sp.predict ,
                          kernel = "gaussian",
                          bw = bw.ND,
                          dMat1 = predict_dist, 
                          dMat2 = net_dist)))
    
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
    Lattice$update_error(y_hat = prediction.latt,
                         y_true = true_signal[sample_],i,j)
    Lattice$update_y_hat(vec = prediction.latt, i=i, j=j)
    ### Reduced Rank Kriging (Ver Hoef) -------------------------------------
    
    matNames = as.character(1:nrow(locs))
    rownames(net_dist) = matNames
    colnames(net_dist) = matNames
    
    knots = create_knots(locations=locs)
    
    predict_net_dist = ND[, sample_]
    rownames(predict_net_dist) =as.character(1:nrow(mesh$nodes))
    colnames(predict_net_dist) = matNames
    
    RR.krig = CovMat_RRKrig(obs, net_dist, knots, 
                            predict_net_dist = predict_net_dist,   
                            cov_model="sph")
    RR.krig$prediction
    RR.krig
    
    
    RR_Krig$update_error(y_hat = LOOCV(RR.krig$CovMat, obs)$cv.pred,
                         y_true = true_signal[sample_],i,j)
    RR_Krig$update_y_hat(RR.krig$cv.pred, i=i, j=j)
    RR_Krig$update_estimate(FEM(RR.krig$prediction, FEMbasis), i, j)
  }
  SR_PDE$compute_mean_field(j)
  GWR$compute_mean_field(j)
  Lattice$compute_mean_field(j)
  RR_Krig$compute_mean_field(j)
}                                     

save(SR_PDE, GWR, Lattice, RR_Krig, locations, OBSERVATIONS,
     folder.name,
     file = paste0(folder.name,"data",".RData"))

# Post processing --------------------------------------------------------------

SimulationBlock <- BlockSimulation(list(SR_PDE, GWR, Lattice, RR_Krig))

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
ORDER = c(1,3,2,4)
pdf(paste0(folder.name,"RMSE.pdf"),width = 12)
boxplot(SimulationBlock, ORDER) +
  labs(title="RMSE", x="observations") +
  theme(legend.position = c(0.90,0.80)) +
  MyTheme
dev.off()

pdf(paste0(folder.name, "domain.pdf"))
plot(mesh, linewidth=0.75)
dev.off()


# setting same color scale
color.min <- rep(0., times = length(SimulationBlock$n_obs))
color.max <- rep(max(true_signal), times = length(SimulationBlock$n_obs))

for(j in 1:length(SimulationBlock$n_obs)){
  for(i in 1:SimulationBlock$num_methods){
    color.min[j] <- min(min(SimulationBlock$Simulations[[i]]$meanField[[j]]$coeff), 
                        color.min[j])
    color.max[j] <- max(max(SimulationBlock$Simulations[[i]]$meanField[[j]]$coeff), 
                        color.max[j])
  }
}
cbind(color.min, color.max)

for(i in 1:SimulationBlock$num_methods){
  for(j in 1:length(SimulationBlock$n_obs)){
    pdf(paste0(folder.name,"estimated_field_same_scale_", 
               SimulationBlock$Simulations[[i]]$method_name,"_",
               SimulationBlock$n_obs[j],".pdf"))
    print(SimulationBlock$Simulations[[i]]$plot_mean_field(j,linewidth=3)+
            viridis::scale_color_viridis(limits=c(color.min[i],color.max[i])))
    print(SimulationBlock$Simulations[[i]]$plot_mean_field(j,linewidth=3)+
            viridis::scale_color_viridis(limits=c(color.min[i],color.max[i])) + 
            theme( legend.position = "none"))
    dev.off()
  }
}

plot.colorbar(FEM(aux_density(mesh$nodes[,1], mesh$nodes[,2]), FEMbasis), 
              colorscale =  viridis, limits = c(color.min[2], color.max[2]),
              width = 2,
              file = paste0(folder.name,"colorbar"))

pdf(paste0(folder.name, "true_field.pdf"), family = "serif", width = 10, height = 10)
plot(FEM(aux_density(mesh$nodes[,1], mesh$nodes[,2]), FEMbasis), linewidth=3) +
  scale_color_viridis(limits=c(color.min[i],color.max[i])) +
  theme( legend.position = "none")

plot(FEM(aux_density(mesh$nodes[,1], mesh$nodes[,2]), FEMbasis), linewidth=4) +
  scale_color_viridis(limits=c(color.min[i],color.max[i])) +
  theme( legend.position = "none")

dev.off()

n_obs = SimulationBlock$n_obs
for(j in 1:length(n_obs)){
  locs = locations[[(j-1)*n_sim + 1]]
  obs = OBSERVATIONS[[(j-1)*n_sim + 1]] 
  pdf(paste0(folder.name,"observations_" , n_obs[j],".pdf"), family = "serif", width = 10, height = 10)
  print(plot(mesh, linewidth=3, color="gray") + geom_point(data=data.frame(x=locs[,1],y=locations[[i]][,2]),
                                          aes(x=x, y=y, color=obs), size=4) + 
    scale_color_viridis())
  
  print(plot(mesh, linewidth=3, color="gray") + geom_point(data=data.frame(x=locs[,1],y=locs[,2]),
                                                aes(x=x, y=y, color=obs), size=5) + 
          scale_color_viridis())
  
  print(plot(mesh, linewidth=3, color="gray") + geom_point(data=data.frame(x=locs[,1],y=locs[,2]),
                                                aes(x=x, y=y, color=obs), size=6) + 
          scale_color_viridis())
    
  print(plot(mesh, linewidth=3, color="gray") + geom_point(data=data.frame(x=locs[,1],y=locs[,2]),
                                          aes(x=x, y=y, color=obs), size=4) + 
    scale_color_viridis() + theme( legend.position = "none"))
  
  
  print(plot(mesh, linewidth=3, color="gray") + geom_point(data=data.frame(x=locs[,1],y=locs[,2]),
                                                aes(x=x, y=y, color=obs), size=5) + 
          scale_color_viridis() + theme( legend.position = "none"))
  
  print(plot(mesh, linewidth=3, color="gray") + geom_point(data=data.frame(x=locs[,1],y=locs[,2]),
                                                aes(x=x, y=y, color=obs), size=6) + 
          scale_color_viridis() + theme( legend.position = "none"))
  dev.off()
}


tmp = list(SR_PDE, GWR, Lattice, RR_Krig)
names(tmp) <- SimulationBlock$method_names

for(j in 1:length(n_obs)){
  locs = locations[[j]]
  for(method in SimulationBlock$method_names){
    #tmp[[method]]$compute_mean_y_hat(j)
    obs = tmp[[method]]$y_hat[[(j-1)*n_sim + 1]]
    pdf(paste0(folder.name,"y_hat_", method , n_obs[j],".pdf"), family = "serif", width = 10, height = 10)
    print(plot(mesh, linewidth=3, color="gray") + geom_point(data=data.frame(x=locs[,1],y=locs[,2]),
                                                  aes(x=x, y=y, color=obs), size=4) + 
            scale_color_viridis())
    
    print(plot(mesh, linewidth=0.75) + geom_point(data=data.frame(x=locs[,1],y=locs[,2]),
                                                  aes(x=x, y=y, color=obs), size=5) + 
            scale_color_viridis())
    
    print(plot(mesh, linewidth=0.75) + geom_point(data=data.frame(x=locs[,1],y=locs[,2]),
                                                  aes(x=x, y=y, color=obs), size=6) + 
            scale_color_viridis())
    
    print(plot(mesh, linewidth=0.75) + geom_point(data=data.frame(x=locs[,1],y=locs[,2]),
                                                  aes(x=x, y=y, color=obs), size=4) + 
            scale_color_viridis() + theme( legend.position = "none"))
    
    
    print(plot(mesh, linewidth=0.75) + geom_point(data=data.frame(x=locs[,1],y=locs[,2]),
                                                  aes(x=x, y=y, color=obs), size=5) + 
            scale_color_viridis() + theme( legend.position = "none"))
    
    print(plot(mesh, linewidth=0.75) + geom_point(data=data.frame(x=locs[,1],y=locs[,2]),
                                                  aes(x=x, y=y, color=obs), size=6) + 
            scale_color_viridis() + theme( legend.position = "none"))
    dev.off()
  }
}

# tabelle
# SR-PDE
cat("--- SR-PDE ---\n")
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
cat("--- KDE 2D ---\n")
rmse_table_lattice <- matrix(Lattice$errors, nrow=SimulationBlock$num_sim, 
                        ncol=length(SimulationBlock$n_obs))
colMeans(rmse_table_lattice)
apply(rmse_table_lattice, MARGIN = 2, sd)

# RR Krig
cat("--- RR Krig ---\n")
rmse_table_krig <- matrix(RR_Krig$errors, nrow=SimulationBlock$num_sim, 
                          ncol=length(SimulationBlock$n_obs))
colMeans(rmse_table_krig)
apply(rmse_table_krig, MARGIN = 2, sd)

tmp <- cbind(colMeans(rmse_table_sr_pde), apply(rmse_table_sr_pde, MARGIN = 2, sd),
      colMeans(rmse_table_gwr), apply(rmse_table_gwr, MARGIN = 2, sd),
      colMeans(rmse_table_lattice), apply(rmse_table_lattice, MARGIN = 2, sd),
      colMeans(rmse_table_krig), apply(rmse_table_krig, MARGIN = 2, sd))

tmp[1,]
