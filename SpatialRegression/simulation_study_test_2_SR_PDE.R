graphics.off()
rm(list=ls())

# Generalize LinearSpatial Regression over Linear Networks ---------------------
setwd("SpatialRegression/")

# loading packages and auxiliary functions
source("../packages.R")
source("../utils.R")
source("../settings.R")
source("../plot.R")
source("../Simulation.R")
source("create_knots.R")
source("rank_reduced_kriging.R")

tests.names = c("test_1", "test_2", "test_3")
ntest = 2

# available domains
domains = c("estevan", "ontario",  "simplenet") 

# Competing methods (available in spatstat package) ----------------------------

method_names = c("SR-PDE", "GWR", "Lattice", "RR-Krig")
methods = c(T,T,F,F)
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
  
  n_obs = as.integer(c(100, 250, 500, 1000))
  lambda = 10^seq(from=-5,to=0.,length.out=20)
  n_sim=30L
}
if(ntest==3){
  n_obs = as.integer(c(50, 100, 150, 250)) # numbers of occurences
  lambda = 10^seq(from=-4, to=-3,length.out = 20)
  sources = c(6,8)         
  n_sim = 30L
  # test locations
  locs.test = runiflpp(1000, spatstat.linnet)
  locs.test = cbind(locs.test$data$x, locs.test$data$y)
}
set.seed(1234)

## Building true signal --------------------------------------------------------

# in ../utils.R 
ND = compute_dist_matrix(points1= mesh$nodes, 
                         points2 = mesh$nodes, 
                         L =spatstat.linnet)

# aux_test is defined in settings.R
if(ntest==1){
  sample.idx =c(2378, 1271, 1802, 63)
  field = aux_test_regression(ND,source = sample.idx[1], sigma=0.1)
  field.2 = aux_test_regression(ND,source = sample.idx[2], sigma=0.125)
  field.3 = aux_test_regression(ND,source = sample.idx[3], sigma=0.1)
  field.4 = aux_test_regression(ND,source = sample.idx[4], sigma=0.25)
  true_signal <- field + field.2 + field.3 + field.4 
  plot(FEM(true_signal, FEMbasis), linewidth=0.75)  + scale_color_viridis()
  
  observations <- true_signal + rnorm(nnodes, mean=0, sd = 0.05*diff(range(true_signal)))
}
if(ntest==2){
  # true spatial field
  sample.idx = c(250, 1000, 950, 63)
  field = aux_test_regression(ND,source = sample.idx[1],sigma = 0.175)
  field.2 = aux_test_regression(ND,source = sample.idx[2],sigma = 0.175)
  field.3 = aux_test_regression(ND,source = sample.idx[3],sigma = 0.35)
  field.4 = aux_test_regression(ND,source = sample.idx[4],sigma = 0.35)
  true_signal <- field + field.2 + field.3 + field.4 
  plot(FEM(true_signal, FEMbasis), linewidth=0.75)  + scale_color_viridis()
  
  # covariates
  X = matrix(nrow=nnodes, ncol=2)
  X[,1] = rnorm(nnodes, mean=0.5,sd=0.25)
  X[,2] = 1/4 * aux_test_regression_cov(ND, source=63, sigma=1.5)
  plot(FEM(X[,2],FEMbasis), linewidth=0.75) + viridis::scale_color_viridis()
  betas = c(1.,1.)
  
  true_signal = field + X%*%betas # g(mu) = true_signal
  
  mu<-inv.link(true_signal)
  observations <- as.numeric(rpois(nnodes, lambda = mu))
}


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
# Loop -------------------------------------------------------------------------
for(j in 1:length(n_obs)){  
  for(i in 1:n_sim){
    
    sample_ = sample(1:nnodes, size=n_obs[j])
    
    locs = mesh$nodes[sample_,]
    obs = observations[sample_]
    
    net_dist = ND[sample_, sample_]
    
    ### SR-PDE ### ------------------------------------------------------------- 
    output_CPP = smooth.FEM(observations = obs, 
                            locations = locs, covariates = X[sample_,],
                            FEMbasis = FEMbasis,
                            lambda = lambda,
                            lambda.selection.criterion = "grid",
                            lambda.selection.lossfunction = "GCV",
                            DOF.evaluation = "stochastic",
                            family = FAMILY) 
    
    lambda_pos <- output_CPP$optimization$lambda_position
    y_hat = inv.link(eval.FEM(FEM(output_CPP$solution$f[,lambda_pos], FEMbasis), locations = locs) + 
            X[sample_,]%*% output_CPP$solution$beta[,lambda_pos])
    SR_PDE$update_estimate(estimate = output_CPP$fit.FEM,i = i, j=j)
    SR_PDE$update_error(y_hat=y_hat , y_true=inv.link(true_signal[sample_]), i=i, j=j)
    
    ### GWR ### ----------------------------------------------------------------
    data_ = data.frame(observations = obs,
                       X1 = X[sample_,1],
                       X2 = X[sample_,2])
    Sp.data = SpatialPointsDataFrame(coords = locs,
                                     data = data_)
    
    bw.ND = bw.ggwr(observations ~ X1 + X2,
                    family=FAMILY,
                    data = Sp.data, 
                    approach="AIC", 
                    kernel="gaussian",
                    dMat = net_dist)
    
    GWR.ND = ggwr.basic(observations ~ X1 + X2,
                        family = FAMILY,
                        data = Sp.data, 
                        kernel = "gaussian",
                        bw = bw.ND,
                        dMat = net_dist)
    
    GWR$update_error(y_hat = GWR.ND$SDF$yhat,
                     y_true = inv.link(true_signal[sample_]),i,j) 
  }
  SR_PDE$compute_mean_field(j)
}                                     

save(SR_PDE, GWR, folder.name,
     file = paste0(folder.name,"data",".RData"))

# Post processing --------------------------------------------------------------

SimulationBlock <- BlockSimulation(list(SR_PDE, GWR))

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
pdf(paste0(folder.name,"test_2_RMSE.pdf"))
boxplot(SimulationBlock, ORDER=c(1,2)) + 
          labs(title="RMSE", x="observations") +
          theme(legend.position = c(0.90,0.90)) +
          MyTheme
dev.off()

pdf(paste0(folder.name, "test_2_domain.pdf"))
plot(mesh, linewidth=0.75)
dev.off()

for(i in 1:length(n_obs)){
  pdf(paste0(folder.name,"test_2_estimated_field_",n_obs[i],".pdf"))
  print(SR_PDE$plot_mean_field(i,linewidth=0.75))
  dev.off()
}
