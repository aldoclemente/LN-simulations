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

sf_network <- as_sfnetwork(spatstat.linnet)

x0 = 0.35 
y0 = 0.345
Lx = 0.125
Ly = 0.08
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

sf_filtered = st_filter(sf_network, poly, .pred = st_intersects)

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

lab_nodes  = list()
idxs = 1:nrow(filtered$nodes)
for(i in idxs) lab_nodes[[i]] = substitute(xi[i], list(i = as.numeric(i) ))

lab_edges  = list()
idxs = 1:nrow(filtered$edges)
for(i in idxs) lab_edges[[i]] = substitute(e[i], list(i = as.numeric(i) ))

{
  pdf(paste0(destdir, "eastbourne_zoom.pdf"))
  for(i in 1:length(linewidth)){
    print(plot(filtered, linewidth = linewidth[i], color="gray50"))
    print(plot(filtered, linewidth = linewidth[i], color="gray50") +
            geom_point(data = data.frame(x = filtered$nodes[,1],
                                         y = filtered$nodes[,2]),
                       aes(x = x, y = y), color = "black", size=size))
    
    TMP = plot(filtered, linewidth = linewidth[i], color="gray50") +
            geom_point(data = data.frame(x = filtered$nodes[,1],
                                         y = filtered$nodes[,2]),
                       aes(x = x, y = y), color = "black", size=size)
    for(j in 1:nrow(filtered$nodes)){
      TMP = TMP + annotate("text", x = filtered$nodes[j,1], y = (filtered$nodes[j,2]-delta/4.5), 
                     label = lab_nodes[[j]],
                     size=4)
    }
    for(j in 1:nrow(filtered$edges)){
      TMP = TMP + annotate("text", x = (filtered$nodes[filtered$edges[j,1],1]+filtered$nodes[filtered$edges[j,2],1])/2, 
                                   y = (filtered$nodes[filtered$edges[j,1],2]+filtered$nodes[filtered$edges[j,2],2])/2 + delta/6, 
                           label = lab_edges[[j]],
                           size=4)
    }
    print(TMP)
            
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

{
pdf(paste0(destdir,"hat_function.pdf"),family = "serif", width = h, height = w)
for(i in 1:nrow(filtered$nodes)){
  center = i
  coeff <- hat_function(filtered_ref$nodes, filtered_ref, center, delta, dist_mat)
  coeff[center] <- 1

  print(plot(FEM(coeff, FEMbasis), linewidth = linewidth) + 
    scale_color_gradientn(colors = jet.col(100), limits = c(0, 1)) + theme( legend.position = "none"))
}
dev.off()
}
plot.colorbar(FEM(coeff, FEMbasis), 
              colorscale =  viridis, width = 2, limits = c(0,1),
              file = paste0(destdir,"colorbar"))






