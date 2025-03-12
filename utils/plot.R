if(!require(pacman)) install.packages("pacman")
pacman::p_load("ggplot2", "ggforce", "viridis", "plotrix", "latex2exp", "colorspace", "grid", "gridExtra")

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

plot.colorbar <- function(x, limits= NULL, colorscale = jet.col, 
                          horizontal = FALSE, cex.axis = 2, width= 5,
                          file = "plot.pdf"){
  
  mesh <- x$FEMbasis$mesh
  
  coeff <- apply(mesh$edges, MARGIN=1, FUN = function(edge){
    mean(x$coeff[edge,])
  })
  
  if(is.null(limits)) limits = c(min(coeff, na.rm = T), max(coeff, na.rm = T))
  cmin = limits[1]; cmax=limits[2]
  
  exps <- -15:15
  range_ <- round( diff(range(limits)), digits = 2)
  cmin_exp <- which(floor(log10(abs(signif(signif(cmin, digits = 2) / 10^exps, digits = 0))))==0)
  cmax_exp <- which(floor(log10(abs(signif(signif(cmax, digits = 2) / 10^exps, digits = 0))))==0)
  k <- exps[max(cmin_exp, cmax_exp)]
  
  at = seq(0, 100, length.out=5)
  labels = as.character(round(seq(cmin*10^(-k), cmax*10^(-k), length.out=5), 2))
  text_ <- ifelse(k != 0, paste0("$\\times 10^{", k,"}$"), " ")
  
  #x11(width = 11, height = 3)
  pdf(paste0(file, "_horiziontal.pdf"), family = "serif", width = 11, height = 3)
  par(mai = c(1,0.75,0,0))
  plot(c(0, 112.5), c(0, 15), type = "n", xlab = "", ylab = "", xaxt = "n", yaxt = "n", frame.plot = F)
  gradient.rect(0, 0, 100, width, col = colorscale(1000), border = "black")
  axis(1, at = at, labels = labels, # at = c(0,33.33,66.66,100)
       cex.axis = cex.axis, lwd.ticks = 0, lwd = 0) # lwd.ticks = 2, 
  text(107, 2, TeX(text_), cex = 2)
  dev.off()
  
  pdf(paste0(file, "_vertical.pdf"), family = "serif", width = 3, height = 11)
  #x11(width = 3, height = 11)
  par(mai = c(1,0.75,0,0))
  plot(c(0, 15), c(0, 112.5), type = "n", xlab = "", ylab = "", xaxt = "n", yaxt = "n", frame.plot = F)
  gradient.rect(0, 0, width, 100, col = colorscale(1000), border = "black", gradient = "y")
  axis(4, at = at, labels = labels, # at = c(0,33.33,66.66,100)
       cex.axis = cex.axis, lwd.ticks = 0, lwd = 0,  line=-11.5+width) # lwd.ticks = 2, 
  text(2.5, 107, TeX(text_), cex = 2)
  dev.off()
  
}
