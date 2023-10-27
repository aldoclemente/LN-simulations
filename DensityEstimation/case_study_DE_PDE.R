# Chicago crimes case study ----------------------------------------------------

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
data("chicago")

# Competing methods ------------------------------------------------------------
# methods[1] -> DE-PDE
# methods[2] -> KDE-HEAT  (available in spatstat package))
# methods[3] -> KDE-ES    (available in spatstat package, very slow !)
# methods[4] -> KDE-2D    (available in spatstat package)
# methods[5] -> VORONOI   (available in spatstat package, slow      !)

method_names = c("DE-PDE", "KDE-HEAT", "KDE-ES", "KDE-2D", "VORONOI")
methods = c(T,T,F,T,T)
method_names = method_names[methods]

# Fixing domain ----------------------------------------------------------------

  vertices = cbind(chicago$domain$vertices$x, chicago$domain$vertices$y)
  edges = cbind( chicago$domain$from, chicago$domain$to)
  
  mesh = create.mesh.1.5D(nodes = vertices, edges = edges)
  
  #mesh = normalize_mesh(mesh)
  res <- normalize_mesh_unit(mesh)
  mesh <- res$mesh
  LN = as.linnet(mesh)
  
  
  x.norm = (chicago$data$x - res$x.min)/(res$x.max-res$x.min)
  y.norm = (chicago$data$y - res$y.min)/(res$y.max-res$y.min)
  
  chicago.norm = spatstat.linnet::lpp(X= ppp(x.norm, y=y.norm, 
                                             window= LN$window), 
                                      L= LN)
  
spat.stat.linnet = chicago.norm
mesh = as.mesh.1.5D(spat.stat.linnet$domain)
FEMbasis = create.FEM.basis(mesh)

lambda = 10^seq(from=-6, to=-3,length.out = 20)
n = nrow(spat.stat.linnet$data)

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
data <- as.data.frame(cbind(spat.stat.linnet$data$x, spat.stat.linnet$data$y))
colnames(data) <- c("x","y")

K <- 10L
KfoldObj <- KfoldsCrossValidation(data, seed=3145L, K=K)
# DE-PDE -----------------------------------------------------------------------
DE_PDE <- DensityEstimationCaseStudy(method_name=method_names[1],
                                      n_obs = KfoldObj$num_obs_kFold,
                                      FEMbasis = FEMbasis)  
# KDE-HEAT ---------------------------------------------------------------------
KDE_HEAT <- DensityEstimationCaseStudy(method_name=method_names[2],
                                        n_obs = KfoldObj$num_obs_kFold,
                                        FEMbasis = FEMbasis)  
# KDE-2D -----------------------------------------------------------------------
KDE_2D <- DensityEstimationCaseStudy(method_name=method_names[3],
                                      n_obs = KfoldObj$num_obs_kFold,
                                      FEMbasis = FEMbasis) 
# VORONOI ----------------------------------------------------------------------
VORONOI <- DensityEstimationCaseStudy(method_name=method_names[4],
                                       n_obs = KfoldObj$num_obs_kFold,
                                       FEMbasis = FEMbasis)  
for(j in 1:K){
  cat(paste("-------------------  ", j, " / ", K,"  -------------------\n", sep="") )
  tmp = KfoldObj$get_data(j)
  train_data = tmp$train_data
  test_data = tmp$test_data
  
  # DE-PDE ---------------------------------------------------------------------
  start = Sys.time()
  invisible(capture.output(output_CPP <- DE.FEM(data = cbind(train_data$x, train_data$y), FEMbasis = FEMbasis,
                          lambda = lambda,
                          preprocess_method ="RightCV",
                          nfolds = 10)))
  cat(paste("- DE-PDE DONE, time elapsed = ", 
            difftime(Sys.time(),start, units = "mins")," mins \n"))
  
  DE_PDE$update_estimate(estimate=FEM(exp(output_CPP$g), FEMbasis),j)
  DE_PDE$update_error(test_data,j) 
  
  # Training Point Pattern over spat.stat.linnet Road Network
  PP_train = lpp(X = ppp(x = train_data$x, y = train_data$y, window = spat.stat.linnet$domain$window),
                 L = spat.stat.linnet$domain)
  # KDE-HEAT  ------------------------------------------------------------------
  start = Sys.time()
  invisible(capture.output(bw <- bw.lppl(X = PP_train)))
  invisible(capture.output(output_KDE_HEAT <- densityHeat(x = PP_train, sigma = as.numeric(bw)))) 
  cat(paste0("- KDE-HEAT DONE, time elapsed = ", 
             difftime(Sys.time(),start, units = "mins")," mins \n"))
  
  coef_ <- as.linfun(output_KDE_HEAT/nrow(train_data))(mesh$nodes[,1], mesh$nodes[,2])
  KDE_HEAT$update_estimate(estimate = FEM(coef_, FEMbasis),j)
  KDE_HEAT$update_error(test_data,j) 
  
  # KDE-2D ---------------------------------------------------------------------
  start = Sys.time()
  invisible(capture.output(bw <- bw.scott(X = PP_train)))
  invisible(capture.output(output_KDE_2D <- densityQuick.lpp(X = PP_train, sigma = bw))) #, at = points)
  cat(paste0("- KDE-2D DONE, time elapsed = ", 
             difftime(Sys.time(),start, units = "mins")," mins \n"))
  
  coef_ <- as.linfun(output_KDE_2D/nrow(train_data))(mesh$nodes[,1], mesh$nodes[,2])
  KDE_2D$update_estimate(estimate = FEM(coef_, FEMbasis),j)
  KDE_2D$update_error(test_data,j=j) 
  
  # KDE-VORONOI ----------------------------------------------------------------
  start = Sys.time()
  invisible(capture.output(bw = bw.voronoi(X = PP_train) ))
  invisible(capture.output(output_VORONOI <- densityVoronoi(X = PP_train, sigma = bw)))
  cat(paste0("- VORONOI DONE, time elapsed = ", 
             difftime(Sys.time(),start, units = "mins")," mins \n"))
  
  coef_ <- as.linfun(output_VORONOI/nrow(train_data))(mesh$nodes[,1], mesh$nodes[,2])
  VORONOI$update_estimate(estimate = FEM(coef_, FEMbasis),j)
  VORONOI$update_error(test_data,j)  
  
  DE_PDE$compute_mean_field(j)
  KDE_HEAT$compute_mean_field(j)
  KDE_2D$compute_mean_field(j)
  VORONOI$compute_mean_field(j)
}

save(DE_PDE, KDE_HEAT, KDE_2D, VORONOI, folder.name,
     file = paste0(folder.name,"data",".RData"))

# Post processing --------------------------------------------------------------
SimulationBlock <- BlockCaseStudy(list(DE_PDE, KDE_HEAT, KDE_2D, VORONOI))

title.size <- 26
MyTheme <- theme(
  axis.ticks = element_blank(),
  axis.title = element_text(size=title.size),
  title = element_text(size=title.size),
  plot.title = element_text(hjust = 0.5),
  legend.text = element_text(size=title.size-5),
  legend.key.size = unit(1,"cm"),
  legend.key.height = unit(1,"cm"),
  legend.title = element_blank(),
  legend.background = element_rect(fill="white", color="white",
                                   linewidth =c(1,0.5))
)
SimulationBlock$method_names
pdf(paste0(folder.name,"case_study_CV_error.pdf"), width = 7)
boxplot(SimulationBlock, ORDER=c(1,2,3,4)) + 
  labs(title="CV error", x="observations") +
  MyTheme
dev.off()

pdf(paste0(folder.name, "case_study_domain.pdf"))
plot(mesh, linewidth=1)
dev.off()

pdf(paste0(folder.name, "case_study_domain_obs.pdf"))
plot(mesh, linewidth=1) + geom_point(data=data.frame(x=x.norm,y=y.norm),
                                     aes(x=x, y=y), color="red3", size=3)
dev.off()

png(paste0(folder.name, "case_study_domain_obs.png"), bg="transparent")
plot(mesh, linewidth=1) + geom_point(data=data.frame(x=x.norm,y=y.norm),
                                     aes(x=x, y=y), color="red3", size=3)
dev.off()

folder.estimates <- paste0(folder.name,"estimates/") 
if(!dir.exists(folder.estimates))
  dir.create(folder.estimates)

# setting same color scale
color.min <- rep(1e5, times = length(SimulationBlock$n_obs))
color.max <- rep(-1e5, times = length(SimulationBlock$n_obs))

for(j in 1:length(SimulationBlock$n_obs)){
  for(i in 1:SimulationBlock$num_methods){
    color.min[j] <- min(min(SimulationBlock$Simulations[[i]]$meanField[[j]]$coeff), 
                        color.min[j])
    color.max[j] <- max(max(SimulationBlock$Simulations[[i]]$meanField[[j]]$coeff), 
                        color.max[j])
  }
}
colors <- cbind(color.min, color.max)

for(i in 1:SimulationBlock$num_methods){
  for(j in 1:length(SimulationBlock$n_obs)){
    pdf(paste0(folder.estimates,"test_3_estimated_field_", 
               SimulationBlock$Simulations[[i]]$method_name,"_",j,".pdf"))
    print(SimulationBlock$Simulations[[i]]$plot_mean_field(j,linewidth=1.) +
            viridis::scale_color_viridis(limits=c(colors[j,1],colors[j,2])) +  # option = "turbo"
            labs(title=SimulationBlock$Simulations[[i]]$method_name) + MyTheme )
    dev.off()
  }
}

# Chicago map 
library(mapview)
library(leaflet)
library(leaflet.providers)

its_bbox = st_bbox(c(xmin = -87.608582 , ymin = 41.785805 , 
                     xmax = -87.5745 , ymax = 41.81), crs = 4326) %>% 
  st_as_sfc()

map.type <- "Esri.WorldTopoMap" # Esri.WorldGrayCanvas "Stadia.AlidadeSmooth" bella ma
map <- mapview(st_as_sf(its_bbox), color="black", alpha=0.7, lwd=0.45,
               layer.name="road-network", map.type=map.type)
map

rect_coords <- data.frame(lng1 =-87.608582 , lat1= 41.785805 , 
                          lng2=-87.5745, lat2=41.81)


map@map <- map@map %>% fitBounds(lng1= (rect_coords$lng1),
                                 lat1= (rect_coords$lat1),
                                 lng2= (rect_coords$lng2),
                                 lat2= (rect_coords$lat2)) %>% hideGroup("road-network")

map

html_fl = paste0(folder.name, "case-study-map.html")
png_fl = paste0(folder.name, "case-study-map.png")

mapshot2(map, url = html_fl, file = png_fl, 
         cliprect = c(260,140, 450,450), zoom=2)#"viewport")
