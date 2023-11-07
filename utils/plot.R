if(!require(pacman)) install.packages("pacman")
pacman::p_load("ggplot2", "ggforce", "viridis", "plotrix")

# Overload plot function for class mesh.1.5D
plot.mesh.1.5D <- function(x, ...){
  mesh <- x
  num_edges= dim(mesh$edges)[1]
  
  x=vector(mode="double", length=2*num_edges)
  y=vector(mode="double", length=2*num_edges)
  grp.nodes=vector(mode="integer", length=2*num_edges)
  
  for(e in 1:num_edges){
    x[(2*(e-1)+1):(2*(e-1)+2)] = c(mesh$nodes[mesh$edges[e,1], 1], mesh$nodes[mesh$edges[e,2],1])
    y[(2*(e-1)+1):(2*(e-1)+2)]=  c(mesh$nodes[mesh$edges[e,1], 2], mesh$nodes[ mesh$edges[e,2],2])
    grp.nodes[(2*(e-1)+1):(2*(e-1)+2)] = rep(e,times=2)
  }
  
  data_plot <- data.frame(x=x, y=y, grp.nodes)
  
  ggplot(data_plot) +
    geom_link2(aes(x = x, y = y, group = grp.nodes),
               lineend = 'round', n = 1, ...) +
    labs(x="",y="",color="", title="") +  
    coord_fixed(ratio=1) + theme_void() 
}

# Overloaded plot function for class FEM
plot.FEM <-function(x, ...){
  FEM <- x
  mesh <- FEM$FEMbasis$mesh
  num_edges= dim(mesh$edges)[1]
  
  x=vector(mode="double", length=2*num_edges)
  y=vector(mode="double", length=2*num_edges)
  coeff=vector(mode="double", length=2*num_edges)
  grp.nodes=vector(mode="integer", length=2*num_edges)
  
  for(e in 1:num_edges){
    x[(2*(e-1)+1):(2*(e-1)+2)] = c(mesh$nodes[mesh$edges[e,1], 1], mesh$nodes[mesh$edges[e,2],1])
    y[(2*(e-1)+1):(2*(e-1)+2)]=  c(mesh$nodes[mesh$edges[e,1], 2], mesh$nodes[ mesh$edges[e,2],2])
    coeff[(2*(e-1)+1):(2*(e-1)+2)]= c(FEM$coeff[mesh$edges[e,1]],
                                      FEM$coeff[mesh$edges[e,2]])  
    grp.nodes[(2*(e-1)+1):(2*(e-1)+2)] = rep(e,times=2)
  }
  
  data_plot <- data.frame(x=x, y=y,
                          coeff=coeff, grp.nodes)
  
  ggplot(data_plot) +
    geom_link2(aes(x = x, y = y, colour = coeff, group = grp.nodes),
               lineend = 'round', n = 10, ...) +
    labs(x="",y="",color="", title="") +  
    coord_fixed(ratio=1) + 
    theme(legend.key.width = unit(0.05,"cm"),
          legend.key.height = unit(2, "cm"))+ theme_void() 
}