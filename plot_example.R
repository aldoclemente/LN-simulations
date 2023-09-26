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

DE_PDE <- DensityEstimationSimulation(method_name = "DE-PDE", n_obs = n_obs, n_sim =  n_sim,FEMbasis =  FEMbasis)

DE_PDE2 <- DensityEstimationSimulation(method_name = "DE-PDE2", n_obs = n_obs, n_sim =  n_sim,FEMbasis =  FEMbasis)

for(j in 1:length(n_obs)){
for(i in 1:n_sim){
  sol_exact=AUX(locs[[j]][,1],locs[[j]][,2])
  data = sol_exact + rnorm(n_obs[j], mean=0, sd=0.05*abs(diff(range(sol_exact))))
  #### Test 1.2: grid with exact GCV
  output_CPP<-smooth.FEM(observations=data, locations = locs[[j]], FEMbasis=FEMbasis, lambda=lambda,
                         lambda.selection.criterion='grid', DOF.evaluation='exact', lambda.selection.lossfunction='GCV')
  DE_PDE$update_estimate(output_CPP$fit.FEM, i,j)
  DE_PDE$update_error(true_field=sol_exact, test_locs=locs[[j]], i,j)
  
  DE_PDE2$update_estimate(output_CPP$fit.FEM, i,j)
  DE_PDE2$update_error(true_field=sol_exact, test_locs=locs[[j]], i,j)
}
  DE_PDE$compute_mean_field(j)
  DE_PDE2$compute_mean_field(j)
}

DE_PDE$plot_mean_field(1L)
DE_PDE$plot_mean_field(2L)

DE_PDE2$plot_mean_field(1L)
DE_PDE2$plot_mean_field(2L)

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

boxplot(DE_PDE) + MyTheme
boxplot(DE_PDE2) + MyTheme

plot(FEM)
library(viridis)
plot(FEM, linewidth=2) + scale_color_viridis()

# BlockSimulation

BLOCK <- BlockSimulation(list(DE_PDE, DE_PDE2))
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
boxplot(BLOCK) + ggFILL + ggBORDER + labs(title="ciao") +  MyTheme

### Case Study ############################
locs <- fdaPDE::refine.by.splitting.mesh.1.5D(mesh)$nodes
colnames(locs) <- c("x","y") 
locs <- as.data.frame(locs)
source("CaseStudy.R")
K <- 3
KfoldObj <- KfoldsCrossValidation(locs,seed=3145L, K=3L)

DE_PDE <- DensityEstimationCaseStudy(method_name = "DE-PDE",n_obs = KfoldObj$num_obs_kFold, FEMbasis = FEMbasis)
DE_PDE2 <- DensityEstimationCaseStudy(method_name = "DE-PDE2",n_obs = KfoldObj$num_obs_kFold, FEMbasis = FEMbasis)

lambda = 10^seq(from=-6, to=-3,length.out = 10)

for(j in 1:K){
  tmp = KfoldObj$get_data(j)
  train_data = tmp$train_data
  test_data = tmp$test_data
  
  output_CPP = fdaPDE::DE.FEM(data = cbind(train_data$x, train_data$y), FEMbasis = FEMbasis,
                          lambda = lambda,
                          preprocess_method ="RightCV",
                          nfolds = 10)
  
  DE_PDE$update_estimate(estimate = FEM(exp(output_CPP$g), FEMbasis),j = j)
  DE_PDE$update_error(test_locs = test_data, j=j) 
  DE_PDE$compute_mean_field(j)
  
  DE_PDE2$update_estimate(estimate = FEM(exp(output_CPP$g), FEMbasis),j = j)
  DE_PDE2$update_error(test_locs = test_data, j=j) 
  DE_PDE2$compute_mean_field(j)
}

boxplot(DE_PDE) + MyTheme

BLOCK_CASE_STUDY <- BlockCaseStudy(c(DE_PDE,DE_PDE2))
boxplot(BLOCK_CASE_STUDY)

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

boxplot(BLOCK_CASE_STUDY) + ggFILL + ggBORDER

### SE VOLESSIMO INVERTIRE ORDINE 
x <- BLOCK_CASE_STUDY$results
x
order_ <- sort(unique(x$method),decreasing = T)
x$method <- factor(x$method, levels=order_)

p<-ggplot(x)+
  geom_boxplot(aes( x=method, y=errors, group=interaction(method),
                    fill=method, color=method))+
  scale_x_discrete(limits=order_)+
  labs(x="", y="") +
  theme(
    axis.ticks.x = element_blank(),
    legend.position = "none")
p  

{
  library(colorspace)
  begin=0.95 # 0.25
  end=0.25   # 0.95
  border_col = darken(viridis(length(order_), begin=begin,end=end), amount=0.25)
  fill_col = viridis(length(order_), begin=begin, end=end)
  BORDER = c()
  FILL = c()
  for(i in 1:length(order_)){
    FILL = append(FILL, fill_col[i])
    BORDER = append(BORDER, border_col[i])
  }
  ggFILL <-scale_fill_manual(values = FILL) 
  ggBORDER <- scale_color_manual(values= BORDER) 
}
p + ggFILL + ggBORDER
