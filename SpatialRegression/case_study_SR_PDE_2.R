# London House Pricing ---------------------------------------------------------

if(!require(pacman)) install.packages("pacman")
pacman::p_load("rstudioapi", "shp2graph", "sfnetworks", "tidygraph", "GWmodel", 
               "fdaPDE", "spatstat")

# setting working directory 
setwd(dirname(getActiveDocumentContext()$path))

source("../utils/utils.R")
source("../utils/plot.R")
source("../utils/Simulation.R")
source("../utils/CaseStudy.R")
# preprocessing ----------------------------------------------------------------
data(LNNT)
network <- st_as_sf(LN.nt) 

# setting the coordinate reference system
st_crs(network) = 27700

# lat/lon reference system 
network <- st_transform(network, 4326)

# building sfnetwork
sfnetwork <- as_sfnetwork(network, directed = FALSE, edges_as_lines = TRUE)

# simplifying 
sfnetwork <- sfnetwork %>%
  activate("edges") %>%
  filter(!edge_is_multiple()) %>%
  filter(!edge_is_loop())

# cleaning
sfnetwork <- sfnetwork %>% 
  convert(to_spatial_subdivision, .clean = TRUE)

# selecting full connected graph
sfnetwork <- sfnetwork %>% 
  convert(to_components, .clean = TRUE, .select = 1L)

mesh <- as.mesh.1.5D(sfnetwork)
FEMbasis <- create.FEM.basis(mesh)

# data 
data(LNHP)
data <- st_as_sf(LN.prop)
data <- unique(data)
# setting the coordinate reference system
st_crs(data) = 27700

# lat/lon reference system 
data <- st_transform(data, 4326)

locs <- st_coordinates(data)

locs <- projection.points.1.5D(mesh, locs)
data$X <- locs[,1]
data$Y <- locs[,2]

linnet <- as.linnet(mesh)
LPP = lpp(locs, linnet)

# Network distance matrix
ND <- pairdist.lpp(LPP) 

# data analysis ----------------------------------------------------------------

data$DATA.IDX = 1:nrow(data)
data$response = log(data$PURCHASE) 
data <- as.data.frame(data)

# Building folders -------------------------------------------------------------
date_ = gsub(":","_",gsub(" ","-",Sys.time()))
if(!dir.exists("data/")) {
  dir.create("data/")
}

if( !dir.exists("data/case-study/")){
  dir.create("data/case-study/")
}

folder.name = paste("data/case-study/",date_,"/",sep="")

if(!dir.exists(folder.name)) {
  dir.create(folder.name)
}

# 10-folds Cross Validation ---------------------------------------------------- 
K <- 10L
KfoldObj <- KfoldsCrossValidation(data, seed=3145L, K=K)
# SR-PDE -----------------------------------------------------------------------
SR_PDE <- SpatialRegressionCaseStudy(method_name="SR-PDE", 
                                     n_obs=KfoldObj$num_obs_kFold,
                                     FEMbasis = FEMbasis)
# GWR -- -----------------------------------------------------------------------
GWR <- SpatialRegressionCaseStudy(method_name="GWR",
                                  n_obs=KfoldObj$num_obs_kFold,
                                  FEMbasis = FEMbasis)
for(j in 1:K){
  tmp = KfoldObj$get_data(j)
  train_data = tmp$train_data
  test_data = tmp$test_data
  
  # GWR ------------------------------------------------------------------------
  train_ND = ND[train_data$DATA.IDX, train_data$DATA.IDX] 
  cross_ND = ND[train_data$DATA.IDX,-train_data$DATA.IDX]
  
  Sp.data.train = SpatialPointsDataFrame(coords = cbind(train_data$X, train_data$Y),
                                         data = train_data)
  
  Sp.data.test = SpatialPointsDataFrame(coords = cbind(test_data$X, test_data$Y),
                                        data = test_data)
  
  bw.ND = bw.gwr(response ~ FLOORSZ + PROF + BATH2, 
                 data = Sp.data.train, 
                 approach="AIC", 
                 kernel="gaussian",
                 dMat = train_ND)
  
  GWR.ND = gwr.predict(response ~ FLOORSZ + PROF + BATH2, 
                       data = Sp.data.train, 
                       predictdata = Sp.data.test,
                       kernel = "gaussian",
                       bw = bw.ND,
                       dMat1 = cross_ND,
                       dMat2 = train_ND)
  
  GWR$update_error(GWR.ND$SDF$prediction, test_data$response, j)
  
  # SR-PDE ---------------------------------------------------------------------
  X = cbind( train_data$FLOORSZ, 
             train_data$PROF,   #, 
             train_data$BATH2) #, 
  lambda = 10^seq(from=-3,to=-1.5,length.out=20) 
  output_CPP = smooth.FEM(observations = train_data$response, 
                          locations = cbind(train_data$X, train_data$Y),
                          FEMbasis = FEMbasis,
                          covariates = X,
                          lambda = lambda,
                          lambda.selection.criterion = "grid",
                          lambda.selection.lossfunction = "GCV",
                          DOF.evaluation = "stochastic")
  
  beta1 = output_CPP$solution$beta[1]
  beta2 = output_CPP$solution$beta[2]
  beta3 = output_CPP$solution$beta[3]
  
  prediction = beta1*test_data$FLOORSZ + 
    beta2*test_data$PROF + 
    beta3*test_data$BATH2 +
    eval.FEM(output_CPP$fit.FEM, cbind(test_data$X, test_data$Y))
  
  SR_PDE$update_estimate(output_CPP$fit.FEM, j)
  SR_PDE$compute_mean_field(j)
  SR_PDE$update_error(prediction, test_data$response,j)
}

save(SR_PDE, GWR, folder.name,
     file = paste0(folder.name,"data",".RData"))

# Post processing --------------------------------------------------------------

SimulationBlock <- BlockCaseStudy(list(SR_PDE, GWR))

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
# SimulationBlock$method_names <- c("SR-PDE", "GWR")
# SimulationBlock$Simulations[[1]]$method_name <- "SR-PDE"
# 
# SimulationBlock$method_names
pdf(paste0(folder.name,"case_study_CV_error.pdf"))
boxplot(SimulationBlock, ORDER=c(1,2)) + 
  labs(title="CV error", x="observations") +
  MyTheme
dev.off()

pdf(paste0(folder.name, "case_study_domain.pdf"))
plot(mesh, linewidth=0.25)
dev.off()

pdf(paste0(folder.name, "case_study_estimate.pdf"))
SR_PDE$plot_mean_field(1,linewidth=0.25)
dev.off()
