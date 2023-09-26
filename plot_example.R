library(fdaPDE)
library(ggplot2)
library(ggforce)
library(viridis)

eps = 1 / 2
x = c(0., 1)
y = c(0.,eps)
vertices = expand.grid(x,y)
vertices = cbind(vertices[,1], vertices[,2])
edges = matrix(c(1,2,1,3,3,4), nrow=3,ncol=2, byrow=T)

mesh = create.mesh.1.5D(vertices, edges, order = 1)
mesh = refine.mesh.1.5D(mesh,delta=0.0125)

nnodes=dim(mesh$nodes)[1]
FEMbasis=create.FEM.basis(mesh)
plot(mesh,size=2)

# Exact solution (pointwise at nodes)
aux.4 = function(x,y){
  h = 1
  source = 4 
  points_ = cbind(x,y)
  idx.ok = which( abs(points_[,2]-mesh$nodes[source,2]) < 10 * .Machine$double.eps )
  coef = vector(mode="numeric", length=length(x))
  
  for(i in idx.ok){
    delta = abs(points_[i,1] - mesh$nodes[source,1])
    if(delta < h ){
      coef[i] = 3*exp(1/(delta^2 / h^2 -1)+1 ) - 1
      
    }
  }
  
  return(coef)
}
aux.3 = function(x,y){
  
  h = eps
  source = 3 
  points_ = cbind(x,y)
  idx.ok = which( abs(points_[,1]-mesh$nodes[source,1]) < 10 * .Machine$double.eps )
  coef = vector(mode="numeric", length=length(x))
  
  for(i in idx.ok){
    delta = abs(points_[i,2] - mesh$nodes[source,2])
    if(delta < h ){
      coef[i] = -1 - 1/h*delta
    }
    
    
  }
  return(coef)
}
aux.1 = function(x,y){
  
  h = 1
  source = 1 
  points_ = cbind(x,y)
  idx.ok = which( abs(points_[,2]-mesh$nodes[source,2]) < 10 * .Machine$double.eps )
  coef = vector(mode="numeric", length=length(x))
  
  for(i in idx.ok){
    delta = abs(points_[i,1] - mesh$nodes[source,1])
    if(delta <= h ){
      coef[i] = -2 - 1/h*delta
    }
  }
  return(coef)
  
}  
AUX = function(x,y){
  
  res = aux.1(x,y) + aux.3(x,y) + aux.4(x,y)
  return(res)
}

sol_exact=AUX(mesh$nodes[,1],mesh$nodes[,2])
FEM = FEM(sol_exact, FEMbasis) 
plot(FEM)
source("plot.R")
plot(FEM)

# Set smoothing parameter
lambda = 10^seq(-4,-2,length.out=10)

source("Simulation.R")
locs = list()
locs[[1]] <-mesh$nodes
locs[[2]] <-refine.by.splitting.mesh.1.5D(mesh)$nodes
n_obs = c(nrow(locs[[1]]), nrow(locs[[2]]))
n_sim = 5L

SR_PDE <- Simulation(method_name = "SR-PDE", 
                     n_obs = n_obs, n_sim =  n_sim,FEMbasis =  FEMbasis)

SR_PDE2 <- Simulation(method_name = "SR-PDE2", 
                      n_obs = n_obs, n_sim =  n_sim,FEMbasis =  FEMbasis)
for(j in 1:length(n_obs)){
for(i in 1:n_sim){
  sol_exact=AUX(locs[[j]][,1],locs[[j]][,2])
  data = sol_exact + rnorm(n_obs[j], mean=0, sd=0.05*abs(diff(range(sol_exact))))
  #### Test 1.2: grid with exact GCV
  output_CPP<-smooth.FEM(observations=data, locations = locs[[j]], FEMbasis=FEMbasis, lambda=lambda,
                         lambda.selection.criterion='grid', DOF.evaluation='exact', lambda.selection.lossfunction='GCV')
  SR_PDE$update_estimate(output_CPP$fit.FEM, i,j)
  SR_PDE$update_error(true_field=sol_exact, test_locs=locs[[j]], i,j)
  
  SR_PDE2$update_estimate(output_CPP$fit.FEM, i,j)
  SR_PDE2$update_error(true_field=sol_exact, test_locs=locs[[j]], i,j)
}
SR_PDE$compute_mean_field(j)
SR_PDE2$compute_mean_field(j)
}

SR_PDE$plot_mean_field(1L)
SR_PDE$plot_mean_field(2L)

SR_PDE2$plot_mean_field(1L)
SR_PDE2$plot_mean_field(2L)

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

boxplot(SR_PDE) + MyTheme
boxplot(SR_PDE2) + MyTheme

plot(FEM)
library(viridis)
plot(FEM, linewidth=2) + scale_color_viridis()

# BlockSimulation

BLOCK <- BlockSimulation(list(SR_PDE, SR_PDE2))
boxplot(BLOCK)

{
library(colorspace)
begin=0.25
end=0.95
border_col = darken(viridis(length(BLOCK$method_names), begin=begin,end=end), amount=0.25)
fill_col = viridis(length(BLOCK$method_names), begin=begin, end=end)
BORDER = c()
FILL = c()
for(i in 1:length(BLOCK$method_names)){
    FILL = append(FILL, fill_col[i])
    BORDER = append(BORDER, border_col[i])
}
ggFILL <-scale_fill_manual(values = FILL) 
ggBORDER <- scale_color_manual(values= BORDER) 
}

boxplot(BLOCK) + ggFILL + ggBORDER
boxplot(BLOCK) + ggFILL + ggBORDER + MyTheme
