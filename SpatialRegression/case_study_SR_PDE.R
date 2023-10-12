
# London House Pricing Case Study  ---------------------------------------------

if(!require(pacman)) install.packages("pacman")
pacman::p_load("rstudioapi")

# setting working directory 
setwd(dirname(getActiveDocumentContext()$path))

# loading packages and auxiliary functions
source("../packages.R")
source("../utils.R")
source("../settings.R")
source("../plot.R")
source("../Simulation.R")
source("../CaseStudy.R")

data(LNNT) # LN.nt   # SpatialLinesDataFrame <- linear network
data(LNHP) # LN.prop # data

spat.stat.linnet = maptools::as.linnet.SpatialLines(LN.nt)

x.min = min(spat.stat.linnet$vertices$x)
y.min = min(spat.stat.linnet$vertices$y)

x.max = max(spat.stat.linnet$vertices$x)
y.max = max(spat.stat.linnet$vertices$y)

x.norm = (spat.stat.linnet$vertices$x - x.min)/(x.max-x.min)
y.norm = (spat.stat.linnet$vertices$y - y.min)/(y.max-y.min)
coords_ = cbind(x.norm, y.norm)

Windows_ = owin(xrange=c(min(x.norm), max(x.norm)), 
                yrange=c(min(y.norm), max(y.norm)) )

spat.stat.linnet = linnet(vertices=as.ppp(coords_, W = Windows_), 
                          edges = cbind(spat.stat.linnet$from, spat.stat.linnet$to),
                          sparse = T)

locs = LN.prop@coords
#which.duplicated = which(duplicated(locs))
which.duplicated = which(duplicated(LN.prop@data))
locs = cbind( (locs[,1]-x.min)/(x.max - x.min), (locs[,2]-y.min)/(y.max-y.min) )

x11()
plot(spat.stat.linnet)
points(locs[,1],locs[,2], pch=16, col="red4")
LPP = lpp(locs, spat.stat.linnet)

# Network distance matrix
ND <- pairdist.lpp(LPP) 

# fdaPDE mesh 
mesh = as.mesh.1.5D(spat.stat.linnet)
FEMbasis = create.FEM.basis(mesh)

# data
dataFrame = LN.prop
dataFrame@coords = cbind(LPP$data$x, LPP$data$y)
dataFrame$X <- dataFrame@coords[,1]; dataFrame$Y <- dataFrame@coords[,2]; 
dataFrame$DATA.IDX = 1:nrow(dataFrame)
dataFrame$PURCHASE = log(dataFrame$PURCHASE) #/10^3 # k pounds

data <- as.data.frame(dataFrame)

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
SR_PDE <- SpatialRegressionCaseStudy(method_name="SR_PDE", 
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
  
  bw.ND = bw.gwr(PURCHASE ~ FLOORSZ + PROF + BATH2, 
                 data = Sp.data.train, 
                 approach="AIC", 
                 kernel="gaussian",
                 dMat = train_ND)
  
  GWR.ND = gwr.predict(PURCHASE ~ FLOORSZ + PROF + BATH2, 
                       data = Sp.data.train, 
                       predictdata = Sp.data.test,
                       kernel = "gaussian",
                       bw = bw.ND,
                       dMat1 = cross_ND,
                       dMat2 = train_ND)
  
  GWR$update_error(GWR.ND$SDF$prediction, test_data$PURCHASE, j)
  
  # SR-PDE ---------------------------------------------------------------------
  X = cbind( train_data$FLOORSZ, 
             train_data$PROF,   #, 
             train_data$BATH2) #, 
  lambda = 10^seq(from=-3,to=-1.5,length.out=20) 
  output_CPP = smooth.FEM(observations = train_data$PURCHASE, 
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
  SR_PDE$update_error(prediction, test_data$PURCHASE,j)
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
SimulationBlock$method_names
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
