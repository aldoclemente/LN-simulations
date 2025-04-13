# London House Pricing ---------------------------------------------------------
set.seed(0)
if(!require(pacman)) install.packages("pacman")
pacman::p_load("shp2graph", "sfnetworks", "tidygraph", "GWmodel", 
               "fdaPDE", "spatstat", "sf", "MetricGraph")

source("utils/utils.R")
source("utils/plot.R")
source("utils/Simulation.R")
source("utils/CaseStudy.R")
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
FEMbasis <- create.FEM.basis(mesh) # 75234 nodes, 96677 edeges

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

# Building folders -------------------------------------------------------------

folder.name = "spatial-regression/"

if(!dir.exists(folder.name)) {
  dir.create(folder.name)
}

folder.name = paste0(folder.name, "london_case_study/")
if(!dir.exists(folder.name)) {
  dir.create(folder.name)
}

# plot -------------------------------------------------------------------------
# pacman::p_load("leaflet", "leaflet.providers", "mapview", "webshot2")
# 
# map.type <- "Esri.WorldTopoMap" # Esri.WorldGrayCanvas "Stadia.AlidadeSmooth" bella ma
# map_1 <- mapview(st_as_sf(sfnetwork, "edges"), color="black", alpha=0.7, lwd=0.45,
#                  layer.name="road-network", map.type=map.type)
# map_2 <- mapview(st_as_sf(data, crs=4326), zcol="response", legend=F, alpha=1, cex=3, 
#                  layer.name="data", map.type=map.type)
# map <- map_1 + map_2
# map
# 
# nodes <- sfnetwork %>%
#   activate("nodes") %>%
#   st_as_sf()
# nodes <- as_Spatial(nodes)
# 
# cntr_crds <- c(mean(coordinates(nodes)[, 1]),
#                mean(coordinates(nodes)[, 2]))
# 
# map@map <- map@map %>% setView(cntr_crds[1], cntr_crds[2], zoom = 10) 
# 
# html_fl = paste0(folder.name, "case-study-map.html")
# png_fl = paste0(folder.name, "case-study-map.png")
# 
# mapshot2(map, url = html_fl, file = png_fl, delay=10, cliprect = c(190,135, 600,500), zoom=1.5  )#"viewport")
# 
# # pdf_fl = paste0(folder.name, "case-study-map.pdf")
# # mapshot2(map, url = html_fl, file = pdf_fl, delay=15, cliprect = c(300,200, 600,500), zoom=1  )#"viewport")
# 
# # add rect
# rect_coords <- data.frame(lng1 =-0.13900 , lat1=51.52954 , lng2=-0.08737, lat2=51.51076)
# map@map <- map@map %>% 
#   addRectangles(lng1=rect_coords$lng1, lat1=rect_coords$lat1, 
#                 lng2=rect_coords$lng2, lat2=rect_coords$lat2, 
#                 fillColor = "transparent",color = "red", 
#                 opacity = 1, weight = 2, group = "rect")
# 
# html_fl = paste0(folder.name, "case-study-map_rect.html")
# png_fl = paste0(folder.name, "case-study-map_rect.png")
# 
# mapshot2(map, url = html_fl, file = png_fl, delay=10, cliprect = c(190,135, 600,500), zoom=1.5  )#"viewport")
# 
# # zoom 
# map1 <- mapview(st_as_sf(sfnetwork, "edges"), color="black", alpha=0.7, lwd=1,
#         layer.name="road-network", map.type=map.type) # Stadia.AlidadeSmooth OpenStreetMap
# 
# map1@map <- map1@map %>% fitBounds(lng1= (rect_coords$lng1-0.00025*rect_coords$lng1),
#                       lat1= (rect_coords$lat1-0.00025*rect_coords$lat1),
#                       lng2= (rect_coords$lng2+0.00025*rect_coords$lng2),
#                       lat2= (rect_coords$lat2+0.00025*rect_coords$lat2))%>% 
#   addRectangles(lng1=rect_coords$lng1, lat1=rect_coords$lat1, 
#                           lng2=rect_coords$lng2, lat2=rect_coords$lat2, 
#                           fillColor = "transparent",color = "red", 
#                           opacity = 1, weight = 2, group = "rect") %>%
#             addLayersControl(
#                   baseGroups ="road-network",
#                   overlayGroups =c("rect", "data"),
#                   options = layersControlOptions(collapsed = FALSE)) %>% 
#             hideGroup(group="rect")
# map1
# 
# html_fl = paste0(folder.name, "case-study-map_zoom.html")
# png_fl = paste0(folder.name, "case-study-map_zoom.png")
# 
# mapshot2(map1, url = html_fl, file = png_fl, delay=10, cliprect = c(190,135, 600,500), zoom=1.5  )#"viewport")

# 10-folds Cross Validation ---------------------------------------------------- 
data <- as.data.frame(data)
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

# WMG <- SpatialRegressionCaseStudy(method_name="WMG", 
#                                   n_obs=KfoldObj$num_obs_kFold,
#                                   FEMbasis = FEMbasis)
# inla.graph = metric_graph$new(V = mesh$nodes, E = mesh$edges) # Too memory demanding


for(j in 1:K){
  cat(paste("-------------------  ", j, " / ", K,"  -------------------\n", sep="") )
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
pdf(paste0(folder.name,"CV_error.pdf"))
boxplot(SimulationBlock, ORDER=c(1,2)) + 
  labs(title="CV error", x="") +
  MyTheme
dev.off()

head(SimulationBlock$results)
col_names <- names(summary(c(1)))
rmse_table <- matrix(0, nrow=SimulationBlock$num_methods, 
                     ncol= (length(col_names) + 1))

for(j in 1:length(SimulationBlock$method_names)){
  rmse_table_method <- SimulationBlock$results[SimulationBlock$results$method == SimulationBlock$method_names[j],1]
  
  rmse_table[j,1:length(col_names)] = summary(rmse_table_method)
  rmse_table[j,length(col_names)+1] = sd(rmse_table_method)
}

colnames(rmse_table) <- c(col_names, "Sd")
rownames(rmse_table) <- SimulationBlock$method_names


write.table(round(rmse_table, digits=4),
            file=paste0(folder.name,"CV_error.txt"))


# pdf(paste0(folder.name, "case_study_domain.pdf"))
# plot(mesh, linewidth=0.25)
# dev.off()
# 
# pdf(paste0(folder.name, "case_study_estimate.pdf"))
# SR_PDE$plot_mean_field(1,linewidth=0.25)
# dev.off()

# library(plotrix)
# png(paste0(folder.name,"legend.png"), bg = "transparent", width=1100, height=300, family = "serif")
# #x11(width=11,height=3)
# par(mai=c(1,0.75,0,0))
# plot(c(0, 112.5), c(0, 30), type= "n", 
#      xlab = "", ylab = "", xaxt = 'n', yaxt = 'n', 
#      frame.plot = F, bg='transparent')
# gradient.rect(0, 0, 100, 5, col = viridis(1000), border = NA)
# axis(1,at=c(0, 200/(675-40)*100,400/(675-40)*100,600/(675-40)*100, 100), 
#      labels=c('','5.3', '6', '6.4', ''), lwd.ticks = 0, cex.axis=4, lwd=0)
# dev.off()


# tmp <- as_Spatial(st_as_sf(sfnetwork, "edges"))
# 
# num_edges = nrow(tmp)
# FEM <- SR_PDE$estimates[[1]]
# coeff <- matrix(nrow=num_edges, ncol=1)
# for(e in 1:num_edges){
#   coeff[e]= (FEM$coeff[mesh$edges[e,1]] + FEM$coeff[mesh$edges[e,2]])/2  
# }
# 
# 
# tmp <-  st_as_sf(tmp) %>% mutate(coeff = as.vector(coeff)) 
# 
# 
# map.type <- "Esri.WorldTopoMap" # Esri.WorldGrayCanvas "Stadia.AlidadeSmooth" bella ma
# map <- mapview(tmp, zcol="coeff", alpha=0.7, lwd=0.45,
#                   layer.name="road-network", map.type = map.type, legend=F)
#  
# cntr_crds <- c(mean(coordinates(nodes)[, 1]),
#                mean(coordinates(nodes)[, 2]))
# 
# nodes <- sfnetwork %>%
#   activate("nodes") %>%
#   st_as_sf()
# nodes <- as_Spatial(nodes)
# 
# map@map <- map@map %>% setView(cntr_crds[1], cntr_crds[2], zoom = 10) 
# 
# html_fl = paste0(folder.name, "case-study-map_estimate.html")
# png_fl = paste0(folder.name, "case-study-map_estimate.png")
# 
# mapshot2(map, url = html_fl, file = png_fl, delay=10, cliprect = c(190,135, 600,500), zoom=1.5  )
