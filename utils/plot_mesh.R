x = c(0., 0.5, 0.5+0.5*sqrt(2)/2, 0.5+0.5*sqrt(2)/2)
y = c(0., 0, 0.5*sqrt(2)/2, -0.5*sqrt(2)/2)
vertices = cbind(x, y)
edges = matrix(c(1,2,2,3,2,4), nrow=3,ncol=2, byrow=T)

mesh = create.mesh.1.5D(vertices, edges, order = 1)
linewidth = 3
source("plot.R")

nnodes <- nrow(mesh$nodes)
lab <- vector(length = nnodes)
lab    <- c(expression(v[1]),
            expression(v[2]),
            expression(v[3]),
            expression(v[4]))
inner <- c(2)
boundary <- c(1,3,4)

pdf("network_notation.pdf")
plot(mesh, linewidth = linewidth, color="gray50") +
  geom_point(data = data.frame(x = mesh$nodes[boundary,1],
                               y = mesh$nodes[boundary,2]),
             aes(x = x, y = y), color = "red2", size = 5) +
  geom_point(data = data.frame(x = mesh$nodes[inner,1],
                               y = mesh$nodes[inner,2]),
             aes(x = x, y = y), color = "green2", shape=15, size = 5) +
  annotate("text", x = mesh$nodes[,1], y = (mesh$nodes[,2]-0.035), label = lab,
           size=7)
dev.off()

pdf("network_discretization_1.pdf")
nnodes <- nrow(mesh$nodes)
lab    <- c(expression(xi[1]),
            expression(xi[2]),
            expression(xi[3]),
            expression(xi[4]))
plot(mesh, linewidth = linewidth, color="gray50") +
  geom_point(data = data.frame(x = mesh$nodes[,1],
                               y = mesh$nodes[,2]),
             aes(x = x, y = y), color = "black", size = 5) +
  annotate("text", x = mesh$nodes[,1], y = (mesh$nodes[,2]-0.035), label = lab,
           size=7)
dev.off()

mesh = refine.mesh.1.5D(mesh,delta=0.125)
nnodes <- nrow(mesh$nodes)
lab    <- c(expression(xi[1]),
            expression(xi[2]),
            expression(xi[3]),
            expression(xi[4]),
            expression(xi[5]),
            expression(xi[6]),
            expression(xi[7]),
            expression(xi[8]),
            expression(xi[9]),
            expression(xi[10]),
            expression(xi[11]),
            expression(xi[12]),
            expression(xi[13]))

pdf("network_discretization_2.pdf")
plot(mesh, linewidth = linewidth, color="gray50") +
  geom_point(data = data.frame(x = mesh$nodes[,1],
                               y = mesh$nodes[,2]),
             aes(x = x, y = y), color = "black", size = 5)+
  annotate("text", x = mesh$nodes[,1], y = (mesh$nodes[,2]-0.035), label = lab,
           size=7)
dev.off()

pdf("network_discretization_3.pdf")
plot(mesh, linewidth = linewidth, color="gray50") + 
  geom_point(data = data.frame(x = mesh$nodes[,1],
                               y = mesh$nodes[,2]),
             aes(x = x, y = y), color = "black", size = 5)+
  annotate("text", x = mesh$nodes[1:4,1], y = (mesh$nodes[1:4,2]-0.035), 
           label = lab[1:4],
           size=7) +
  annotate("text", x = mesh$nodes[6,1], y = (mesh$nodes[6,2]-0.035), 
           label = "...",
           size=7) +
  annotate("text", x = (mesh$nodes[9,1]+0.01), y = (mesh$nodes[9,2]-0.02), 
           label = "...",
           size=7, angle=45) + 
  annotate("text", x = (mesh$nodes[12,1]-0.025), y = (mesh$nodes[12,2]-0.025), 
           label = "...",
           size=7, angle=135)  
dev.off()
