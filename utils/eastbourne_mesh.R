source("plot.R")
load("eastbourne_mesh.RData")

destdir = "Space-Time-LN/"
if(!dir.exists(destdir)) dir.create(destdir)

linewidth = 1.5
size = 2
# network
{
  pdf(paste0(destdir, "eastbourne_raw.pdf"))
  for(i in 1:length(linewidth)){
    print(plot(mesh, linewidth = linewidth[i], color="gray50"))
    print(plot(mesh, linewidth = linewidth[i], color="gray50") +
            geom_point(data = data.frame(x = mesh$nodes[,1],
                                         y = mesh$nodes[,2]),
                       aes(x = x, y = y), color = "black", size=size))
  }
  dev.off()
}


plot(mesh, linewidth = linewidth[1], color="gray50", linetype=2) + 
  geom_point(data = data.frame(x = mesh$nodes[,1],
                               y = mesh$nodes[,2]),
             aes(x = x, y = y), color = "black", size=4) +
  geom_point(data = data.frame(x = 0.5,
                               y = 0.5),
             aes(x = x, y = y), color = "green3", size=4) +
  geom_point(data = data.frame(x = 0.375,
                               y = 0.375),
             aes(x = x, y = y), color = "red3", size=4) +
  geom_point(data = data.frame(x = (0.5+0.125/2), # 0.125
                               y = 0.5- 0.125/2*7),
             aes(x = x, y = y), color = "blue3", size=4) +
  theme(panel.grid = element_line(color = "#8ccde3",
                                  size = 0.75,
                                  linetype = 2))

library(fdaPDE)
delta = 0.025
mesh_ref <- refine.mesh.1.5D(mesh, delta)
source("utils.R")
source("../packages.R")
pacman::p_load("shp2graph", "sfnetworks", "tidygraph", 
               "fdaPDE", "spatstat", "sf")

spatstat.linnet <- as.linnet(mesh_ref)
plot(spatstat.linnet)

plot(mesh_ref, linewidth = linewidth[1], color="gray50") + 
  geom_point(data = data.frame(x = mesh_ref$nodes[,1],
                               y = mesh_ref$nodes[,2]),
             aes(x = x, y = y), color = "black", size=4)

#spatstat.linnet <- as.linnet(mesh)
sf_network <- as_sfnetwork(spatstat.linnet)

x0 = (0.5+0.125/2) #0.35 # caso 1
y0 = 0.5- 0.125/2*7#0.345 # caso 1
Lx = 0.115 # 0.125 # caso 1
Ly = 0.1#0.08 # caso 1
p1 = st_point(c(x0-Lx, y0-Ly))
p2 = st_point(c(x0+Lx, y0-Ly))
p3 = st_point(c(x0+Lx, y0+Ly))
p4 = st_point(c(x0-Lx, y0+Ly))

poly = st_multipoint(c(p1, p2, p3, p4)) %>%
  st_cast("POLYGON") 

{
  x11()
  plot(sf_network, col="gray", cex=0, lwd=linewidth[3])
  plot(poly, border = "red", lty = 4, lwd=linewidth[1] , add = TRUE)
}

sf_filtered = sf_network %>%
  activate("edges") %>%
  st_intersection(poly) %>%
  activate("nodes") %>%
  filter(!node_is_isolated())

#sf_filtered = st_filter(sf_network, poly, .pred = st_intersects)

{ 
  x11()
  par(mfrow=c(1,2))
  plot(sf_network, col = "grey", cex=0, lwd=linewidth)
  plot(poly, border = "red", lty = 4, lwd = 4, add = TRUE)
  plot(sf_network, col = "grey", cex=0, lwd=linewidth)
  plot(sf_filtered, add = TRUE, cex = 0, lwd = linewidth+1)
}

filtered <- as.mesh.1.5D(sf_filtered)

plot(mesh, linewidth = linewidth, color="gray50") +  
  geom_rect(aes(xmin = x0-Lx, xmax = x0 + Lx, 
                ymin = y0-Ly, ymax = y0 + Ly),
            fill = "transparent", color = "red", linewidth = linewidth/2) #, linetype=2)

plot(filtered, linewidth = linewidth) +
  geom_rect(aes(xmin = x0-Lx, xmax = x0 + Lx, 
                ymin = y0-Ly, ymax = y0 + Ly),
            fill = "transparent", color = "red", linewidth = linewidth/2)


{
  pdf(paste0(destdir,"eastbourne_filtered.pdf"))
  for(i in 1:length(linewidth)){
    print(plot(mesh, linewidth = linewidth[i], color="gray50") +  
            geom_rect(aes(xmin = x0-Lx, xmax = x0 + Lx, 
                          ymin = y0-Ly, ymax = y0 + Ly),
                      fill = "transparent", color = "red", linewidth = linewidth/2))
    
    print(plot(mesh, linewidth = linewidth[i], color="gray50") +  
            geom_rect(aes(xmin = x0-Lx, xmax = x0 + Lx, 
                          ymin = y0-Ly, ymax = y0 + Ly),
                      fill = "transparent", color = "red", linewidth = linewidth/2, linetype=2))
  }
  dev.off()
}


{
  pdf(paste0(destdir,"eastbourne_LN.pdf"))
  for(i in 1:length(linewidth)){
    
    print(plot(mesh, linewidth = linewidth[i], color="gray50"))
    print(plot(mesh, linewidth = linewidth[i], color="gray50") +
            geom_point(data = data.frame(x = mesh$nodes[,1],
                                         y = mesh$nodes[,2]),
                       aes(x = x, y = y), color = "black", size=size))
    
    print(plot(mesh_ref, linewidth = linewidth[i], color="gray50") +
            geom_point(data = data.frame(x = mesh_ref$nodes[,1],
                                         y = mesh_ref$nodes[,2]),
                       aes(x = x, y = y), color = "black", size=size))
    
    print(plot(mesh, linewidth = linewidth[i], color="gray50") +  
            geom_rect(aes(xmin = x0-Lx, xmax = x0 + Lx, 
                          ymin = y0-Ly, ymax = y0 + Ly),
                      fill = "transparent", color = "red3", linewidth = linewidth/2))
    
    print(plot(mesh, linewidth = linewidth[i], color="gray50") +  
            geom_rect(aes(xmin = x0-Lx, xmax = x0 + Lx, 
                          ymin = y0-Ly, ymax = y0 + Ly),
                      fill = "transparent", color = "red3", linewidth = linewidth/2, linetype=2))
    
    print(plot(mesh, linewidth = linewidth[i], color="gray50") +
            geom_point(data = data.frame(x = mesh$nodes[,1],
                                         y = mesh$nodes[,2]),
                       aes(x = x, y = y), color = "black", size=size) +   
            geom_rect(aes(xmin = x0-Lx, xmax = x0 + Lx, 
                          ymin = y0-Ly, ymax = y0 + Ly),
                      fill = "transparent", color = "red3", linewidth = linewidth/2))
    
    print(plot(mesh, linewidth = linewidth[i], color="gray50") +
            geom_point(data = data.frame(x = mesh$nodes[,1],
                                         y = mesh$nodes[,2]),
                       aes(x = x, y = y), color = "black", size=size) +  
            geom_rect(aes(xmin = x0-Lx, xmax = x0 + Lx, 
                          ymin = y0-Ly, ymax = y0 + Ly),
                      fill = "transparent", color = "red3", linewidth = linewidth/2, linetype=2))
    
    print(plot(mesh_ref, linewidth = linewidth[i], color="gray50") +
            geom_point(data = data.frame(x = mesh_ref$nodes[,1],
                                         y = mesh_ref$nodes[,2]),
                       aes(x = x, y = y), color = "black", size=size) +   
            geom_rect(aes(xmin = x0-Lx, xmax = x0 + Lx, 
                          ymin = y0-Ly, ymax = y0 + Ly),
                      fill = "transparent", color = "red3", linewidth = linewidth/2))
    
    print(plot(mesh_ref, linewidth = linewidth[i], color="gray50") +
                 geom_point(data = data.frame(x = mesh_ref$nodes[,1],
                                              y = mesh_ref$nodes[,2]),
                            aes(x = x, y = y), color = "black", size=size) +  
            geom_rect(aes(xmin = x0-Lx, xmax = x0 + Lx, 
                          ymin = y0-Ly, ymax = y0 + Ly),
                      fill = "transparent", color = "red3", linewidth = linewidth/2, linetype=2))
  }
  dev.off()
}


{
  pdf(paste0(destdir,"eastbourne_filtered_nodes.pdf"))
  for(i in 1:nrow(filtered$nodes)){
    print(plot(filtered, linewidth = linewidth, color="gray50") +  
            geom_point(data = data.frame(x = filtered$nodes[i,1],
                                         y = filtered$nodes[i,2]),
                       aes(x = x, y = y), color = "red3", size=4) )
  }
dev.off()
}

hat_function <- function(x, mesh, center, thresold, dist_mat){
  res <- matrix(0, nrow=nrow(x), ncol=1)
  for(i in 1:nrow(x)){
    #if(dist[center,i])
    #l <- sqrt(sum((x[i,]-mesh$nodes[center,])^2))
    l <- dist_mat[center,i]
    res[i] <- ifelse( l < thresold, 1-l/thresold,0) 
  }
  
  return(res)
}

filtered_ref = refine.mesh.1.5D(filtered, delta = delta/10)
spatstat_filtered = as.linnet(filtered_ref)
dist_mat = pairdist( lpp(X=filtered_ref$nodes, L=spatstat.linnet) )
FEMbasis <- create.FEM.basis(filtered_ref)
w <- 7 # width
h <- 7 # height

node_idxs = 1:nrow(filtered$nodes)
#node_idxs = node_idxs[-c(30,31,33)] # caso 1
node_idxs = node_idxs[-c(28,29,30)] # caso 2 (in basso)
{
pdf(paste0(destdir,"eastbourne_zoom_hat_functions.pdf"),family = "serif", width = h, height = w)
for(i in node_idxs){
  center = i
  coeff <- hat_function(filtered_ref$nodes, filtered_ref, center, delta, dist_mat)
  coeff[center] <- 1

  print(plot(FEM(coeff, FEMbasis), linewidth = linewidth) + 
    scale_color_gradientn(colors = jet.col(100), limits = c(0, 1)) + theme( legend.position = "none"))
}
dev.off()
}

{
  pdf(paste0(destdir,"eastbourne_zoom_hat_functions_border.pdf"),family = "serif", width = h, height = w)
  for(i in node_idxs){
    center = i
    coeff <- hat_function(filtered_ref$nodes, filtered_ref, center, delta, dist_mat)
    coeff[center] <- 1
    
    print(plot(FEM(coeff, FEMbasis), linewidth = linewidth) + 
            scale_color_gradientn(colors = jet.col(100), limits = c(0, 1)) +
            geom_rect(aes(xmin = x0-Lx, xmax = x0 + Lx, 
                          ymin = y0-Ly, ymax = y0 + Ly),
                      fill = "transparent", color = "red3", linewidth = linewidth/2, linetype=2) + 
            theme( legend.position = "none"))
  }
  dev.off()
}

plot.colorbar(FEM(coeff, FEMbasis), 
              colorscale =  jet.col, width = 2, limits = c(0,1),
              file = paste0(destdir,"eastbourne_hat_colorbar"))


# labels ---
spatstat.linnet <- as.linnet(mesh)
sf_network <- as_sfnetwork(spatstat.linnet)

sf_filtered = sf_network %>%
  activate("edges") %>%
  st_intersection(poly) %>%
  activate("nodes") %>%
  filter(!node_is_isolated())
filtered <- as.mesh.1.5D(sf_filtered)

lab_nodes  = list()
idxs_nodes = 1:nrow(filtered$nodes)
idxs_nodes = idxs_nodes[-c(7,8,9,11)]
for(i in idxs_nodes) lab_nodes[[i]] = substitute(v[i], list(i = as.numeric(i) ))

lab_edges  = list()
idxs_edges = 1:nrow(filtered$edges)
idxs_edges = idxs_edges[-c(3, 5,6,7)]
for(i in idxs_edges) lab_edges[[i]] = substitute(e[i], list(i = as.numeric(i) ))

plot(filtered, linewidth = linewidth, color="gray50") +
  geom_point(data = data.frame(x = filtered$nodes[11,1],
                               y = filtered$nodes[11,2]),
             aes(x = x, y = y), color = "black", size=size)


sizes = c(5,6,7,8)
{
  pdf(paste0(destdir, "eastbourne_zoom.pdf"), family="serif",  width = h, height = w)
  for(i in 1:length(linewidth)){
    print(plot(filtered, linewidth = linewidth[i], color="gray50"))
    RAW = plot(filtered, linewidth = linewidth[i], color="gray50")
    
    for(j in idxs_nodes){
      RAW = RAW +
        geom_point(data = data.frame(x = filtered$nodes[j,1],
                                     y = filtered$nodes[j,2]),
                   aes(x = x, y = y), color = "black", size=size)
    }
    # RAW = RAW + geom_point(data = data.frame(x = filtered$nodes[7,1],
    #                                          y = filtered$nodes[7,2]),
    #                        aes(x = x, y = y), color = "black", size=size)
    # 
    print(RAW)
    #TMP = plot(filtered, linewidth = linewidth[i], color="gray50")
    for(k in sizes){
      TMP =RAW
    for(j in idxs_nodes){
      TMP = TMP + annotate("text", x = filtered$nodes[j,1]+delta/6, y = (filtered$nodes[j,2]-delta/4.5), 
                           label = lab_nodes[[j]],
                           size=k)
    }
    for(j in idxs_edges){
         TMP = TMP + annotate("text", x = (filtered$nodes[filtered$edges[j,1],1]+filtered$nodes[filtered$edges[j,2],1])/2, 
                           y = (filtered$nodes[filtered$edges[j,1],2]+filtered$nodes[filtered$edges[j,2],2])/2 + delta/6, 
                           label = lab_edges[[j]],
                           size=k)
      
    }
    print(TMP)
    }
    
  }
  dev.off()
}

# ---

hat_function <- function(x, mesh, center, thresold, dist_mat){
  res <- matrix(0, nrow=nrow(x), ncol=1)
  for(i in 1:nrow(x)){
    #if(dist[center,i])
    #l <- sqrt(sum((x[i,]-mesh$nodes[center,])^2))
    l <- dist_mat[center,i]
    res[i] <- ifelse( l < thresold, 1-l/thresold,0) 
  }
  
  return(res)
}


destdir = "Space-Time-LN/"
x = c(0., 0.5, 0.5+0.5*sqrt(2)/2, 0.5+0.5*sqrt(2)/2)
y = c(0., 0, 0.5*sqrt(2)/2, -0.5*sqrt(2)/2)
vertices = cbind(x, y)
linewidth = 2.5
size = 5
edges = matrix(c(1,2,2,3,2,4), nrow=3,ncol=2, byrow=T)
source("plot.R")
source("utils.R")
source("../packages.R")
pacman::p_load("shp2graph", "sfnetworks", "tidygraph", 
               "fdaPDE", "spatstat", "sf")

mesh_coarse = create.mesh.1.5D(vertices, edges, order = 1)
delta=0.125
mesh = refine.mesh.1.5D(mesh_coarse, delta = delta)
mesh_ref = refine.mesh.1.5D(mesh, delta/100)
spatstat_filtered = as.linnet(mesh_ref)
dist_mat = pairdist( lpp(X=mesh_ref$nodes, L=spatstat_filtered) )
FEMbasis <- create.FEM.basis(mesh_ref)
w <- 7 # width
h <- 7 # height

node_idxs = 1:nrow(mesh$nodes)
#node_idxs = node_idxs[-c(30,31,33)] # caso 1
#node_idxs = node_idxs[-c(28,29,30)] # caso 2 (in basso)

color.palette = c("jet.col", "inferno", "magma", "plasma")
palette = list(jet.col, inferno, magma, plasma)
for(COL in 1:length(color.palette)){

  {
  pdf(paste0(destdir,"hat_functions_", color.palette[COL] ,".pdf"),family = "serif", width = h, height = w)
  
  for(i in node_idxs){
    center = i
    coeff <- hat_function(mesh_ref$nodes, mesh, center, delta, dist_mat)
    coeff[center] <- 1
    
    print(plot(FEM(coeff, FEMbasis), linewidth = linewidth) + 
            scale_color_gradientn(colors = palette[[COL]](100), limits = c(0, 1)) + theme( legend.position = "none"))
  }
  
  print(plot(mesh, linewidth = linewidth, color="gray50") +
          geom_point(data = data.frame(x = mesh$nodes[,1],
                                       y = mesh$nodes[,2]),
                     aes(x = x, y = y), color = "black", size=size))
  
  for(i in node_idxs){
    center = i
    coeff <- hat_function(mesh_ref$nodes, mesh, center, delta, dist_mat)
    coeff[center] <- 1
    
    print(plot(FEM(coeff, FEMbasis), linewidth = linewidth) + 
            geom_point(data = data.frame(x = mesh$nodes[,1],y = mesh$nodes[,2]),
                      aes(x = x, y = y), color = "gray50", size=size) + 
            scale_color_gradientn(colors = palette[[COL]](100), limits = c(0, 1)) + theme( legend.position = "none"))
  }
  dev.off()
  }
  
  
  plot.colorbar(FEM(coeff, FEMbasis), 
                colorscale =  palette[[COL]], width = 2, limits = c(0,1),
                file = paste0(destdir,"colorbar_", color.palette[COL]))
  
}


{
  pdf(paste0(destdir,"hat_functions_viridis.pdf"),family = "serif", width = h, height = w)
  
  for(i in node_idxs){
    center = i
    coeff <- hat_function(mesh_ref$nodes, mesh, center, delta, dist_mat)
    coeff[center] <- 1
    
    print(plot(FEM(coeff, FEMbasis), linewidth = linewidth) + 
            scale_color_gradientn(colors = viridis(100), limits = c(0, 1)) + theme( legend.position = "none"))
  }
  
  print(plot(mesh, linewidth = linewidth, color="gray50") +
          geom_point(data = data.frame(x = mesh$nodes[,1],
                                       y = mesh$nodes[,2]),
                     aes(x = x, y = y), color = "black", size=size))
  
  dev.off()
  
  pdf(paste0(destdir,"hat_functions2_viridis.pdf"),family = "serif", width = h, height = w)
  for(i in node_idxs){
    center = i
    coeff <- hat_function(mesh_ref$nodes, mesh, center, delta, dist_mat)
    coeff[center] <- 1
    
    print(plot(FEM(coeff, FEMbasis), linewidth = linewidth) + 
            geom_point(data = data.frame(x = mesh$nodes[,1],y = mesh$nodes[,2]),
                       aes(x = x, y = y), color = "black", size=size) +
            scale_color_gradientn(colors = viridis(100), limits = c(0, 1)) + theme( legend.position = "none"))
  }
  dev.off()
}

{
  pdf(paste0(destdir,"hat_functions_magma.pdf"),family = "serif", width = h, height = w)
  
  for(i in node_idxs){
    center = i
    coeff <- hat_function(mesh_ref$nodes, mesh, center, delta, dist_mat)
    coeff[center] <- 1
    
    print(plot(FEM(coeff, FEMbasis), linewidth = linewidth) + 
            scale_color_gradientn(colors = magma(100), limits = c(0, 1)) + theme( legend.position = "none"))
  }
  
  print(plot(mesh, linewidth = linewidth, color="gray50") +
          geom_point(data = data.frame(x = mesh$nodes[,1],
                                       y = mesh$nodes[,2]),
                     aes(x = x, y = y), color = "black", size=size))
  
  dev.off()
  
  pdf(paste0(destdir,"hat_functions2_magma.pdf"),family = "serif", width = h, height = w)
  for(i in node_idxs){
    center = i
    coeff <- hat_function(mesh_ref$nodes, mesh, center, delta, dist_mat)
    coeff[center] <- 1
    
    print(plot(FEM(coeff, FEMbasis), linewidth = linewidth) + 
            geom_point(data = data.frame(x = mesh$nodes[,1],y = mesh$nodes[,2]),
                       aes(x = x, y = y), color = "black", size=size) +
            scale_color_gradientn(colors = magma(100), limits = c(0, 1)) + theme( legend.position = "none"))
  }
  dev.off()
}

{
  pdf(paste0(destdir,"hat_functions_inferno.pdf"),family = "serif", width = h, height = w)
  
  for(i in node_idxs){
    center = i
    coeff <- hat_function(mesh_ref$nodes, mesh, center, delta, dist_mat)
    coeff[center] <- 1
    
    print(plot(FEM(coeff, FEMbasis), linewidth = linewidth) + 
            scale_color_gradientn(colors = inferno(100), limits = c(0, 1)) + theme( legend.position = "none"))
  }
  
  print(plot(mesh, linewidth = linewidth, color="gray50") +
          geom_point(data = data.frame(x = mesh$nodes[,1],
                                       y = mesh$nodes[,2]),
                     aes(x = x, y = y), color = "black", size=size))
  
  dev.off()
  
  pdf(paste0(destdir,"hat_functions2_inferno.pdf"),family = "serif", width = h, height = w)
  for(i in node_idxs){
    center = i
    coeff <- hat_function(mesh_ref$nodes, mesh, center, delta, dist_mat)
    coeff[center] <- 1
    
    print(plot(FEM(coeff, FEMbasis), linewidth = linewidth) + 
            geom_point(data = data.frame(x = mesh$nodes[,1],y = mesh$nodes[,2]),
                       aes(x = x, y = y), color = "black", size=size) +
            scale_color_gradientn(colors = inferno(100), limits = c(0, 1)) + theme( legend.position = "none"))
  }
  dev.off()
}


