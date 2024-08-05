# plot simplenet 

source("utils.R")
source("plot.R")
source("../settings.R")
library(spatstat)
library(fdaPDE)

hat_function <- function(x, center, dist){
  res <- matrix(0, nrow=nrow(x), ncol=1)
  for(i in 1:nrow(x)){
    l <- sqrt(sum((x[i,]-mesh$nodes[center,])^2))
    res[i] <- ifelse( l < dist, 1-l/dist,0) 
  }

  return(res)
}


load(simplenet)
mesh = as.mesh.1.5D(simplenet)
#fdaPDE:::plot.mesh.1.5D(mesh)

pdf("simplenet.pdf", family = "serif", width = 10, height = 10)
plot(mesh, linewidth=3, color="gray") + 
  geom_point(data=data.frame(x=mesh$nodes[,1],y=mesh$nodes[,2]),
             aes(x=x, y=y), size=4, color="black") + theme( legend.position = "none")
mesh = refine.mesh.1.5D(mesh,0.125)
plot(mesh, linewidth=3, color="gray") + 
  geom_point(data=data.frame(x=mesh$nodes[,1],y=mesh$nodes[,2]),
             aes(x=x, y=y), size=4, color="black") + theme( legend.position = "none")


mesh_ref = refine.mesh.1.5D(mesh,0.015)
FEMbasis = create.FEM.basis(mesh_ref)
coeff = hat_function(mesh_ref$nodes, center = 4L, dist = 0.125)
coeff[4] = 1
plot(FEM(coeff, FEMbasis), linewidth = 3) + scale_color_viridis(limits = c(0,1)) + theme( legend.position = "none")
plot(FEM(coeff, FEMbasis), linewidth = 4) + scale_color_viridis(limits = c(0,1)) + theme( legend.position = "none")
dev.off()

plot.colorbar(FEM(coeff, FEMbasis), 
              colorscale =  viridis, width = 2, limits = c(0,1),
              file = "simplenet_colorbar")
