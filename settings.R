setting <-function(network = "" ){
  
  if(network == "ontario"){
  data("ORN")
  
  mesh = as.mesh.1.5D(ORN.nt)
  mesh = normalize_mesh_unit(mesh)$mesh
  mesh = refine.mesh.1.5D(mesh, delta=0.0125)
 
  FEMbasis = create.FEM.basis(mesh)
  
  nnodes = nrow(mesh$nodes)
  spatstat.linnet = as.linnet(mesh)
  
  res = list(mesh = mesh, FEMbasis=FEMbasis, nnodes = nnodes,
             spatstat.linnet = spatstat.linnet)
  return(res)
  }else if( network == "estevan"){
    
    data("ERN_OSM_correct")
    
    mesh = as.mesh.1.5D(ERN_OSM_cor.nt)
    mesh = normalize_mesh_unit(mesh)$mesh
    mesh = refine.mesh.1.5D(mesh, delta=0.0125)
      
    FEMbasis = create.FEM.basis(mesh)
    
    nnodes = nrow(mesh$nodes)
    spatstat.linnet = as.linnet(mesh)
    
    res = list(mesh = mesh, FEMbasis=FEMbasis, nnodes = nnodes,
               spatstat.linnet = spatstat.linnet)

  }else if(network=="simplenet"){
    data("simplenet")
    mesh = as.mesh.1.5D(simplenet)
    mesh = refine.mesh.1.5D(mesh,0.025)
    
    FEMbasis = create.FEM.basis(mesh)
    nnodes = nrow(mesh$nodes)
    spatstat.linnet = as.linnet(mesh)
    
    res = list(mesh=mesh, FEMbasis = FEMbasis, nnodes=nnodes,
               spatstat.linnet = spatstat.linnet)  
  }
}

auxiliary_test1 = function(x, y, seg, tp, sigma= 0.125, 
                           nodes.lpp = ppp(x = mesh$nodes[,1], y = mesh$nodes[,2], 
                                           window = owin(xrange = c(min(mesh$nodes[,1]),max(mesh$nodes[,1])),
                                                         yrange = c(min(mesh$nodes[,2]),max(mesh$nodes[,2])))),
                           L = spatstat.linnet,
                           source = sources)
{ 
  PP = ppp(x = x, y = y, window = nodes.lpp$window)
  ND = crossdist.lpp(lpp(nodes.lpp, L), lpp(PP, L))
  
  return(   0.25 * 1/sqrt(2*pi*sigma^2) * exp(-ND[source[1],]^2/(2*sigma^2)) + 
              0.25 * 1/sqrt(2*pi*sigma^2) * exp(-ND[source[2],]^2/(2*sigma^2)))
  
  
}

aux_density = auxiliary_test1

f.exp = function(ND,source=63,sigma=0.125){
  
  nnodes = nrow(ND)
  res = vector(mode="numeric", length=nnodes)
  
  source.1 = source
  distances.1 = ND[source.1,]
  other.1 = vector(mode="integer")
  for(i in 1:nnodes ){
    res[i] = exp(-distances.1[i]/sigma)  
  }
  
  return(res)
}

aux_test_regression =f.exp

f.sin = function(ND,source=63,sigma=1e-5){
  
  #sigma = 0.125/24
  nnodes = nrow(ND)
  res = vector(mode="numeric", length=nnodes)
  
  source.1 = source
  distances.1 = ND[source.1,]
  other.1 = vector(mode="integer")
  for(i in 1:nnodes){
    
    res[i] = sin(2*pi*distances.1[i]/sigma)*cos(2*pi*distances.1[i]/sigma)
    
  }
  
  return(res)
}

aux_test_regression_cov = f.sin
