setting <-function(network = "" ){
  
  if(network == "ontario"){
  data("ORN")
  
  mesh = as.fdaPDE.SpatialLinesDataFrame(ORN.nt)
  mesh = normalize_mesh_unit(mesh)$mesh
 
  FEMbasis = create.FEM.basis(mesh)
  
  nnodes = nrow(mesh$nodes)
  spatstat.linnet = as.spatstat.linnet.fdaPDE(mesh)
  
  res = list(mesh = mesh, FEMbasis=FEMbasis, nnodes = nnodes,
             spatstat.linnet = spatstat.linnet)
  return(res)
  }else if( network == "estevan"){
    
    data("ERN_OSM_correct")
    
    mesh = as.fdaPDE.SpatialLinesDataFrame(ERN_OSM_cor.nt)
    mesh = normalize_mesh_unit(mesh)$mesh
    mesh = refine.mesh.1.5D(mesh, delta=0.0125)
      
    FEMbasis = create.FEM.basis(mesh)
    
    nnodes = nrow(mesh$nodes)
    spatstat.linnet = as.spatstat.linnet.fdaPDE(mesh)
    
    res = list(mesh = mesh, FEMbasis=FEMbasis, nnodes = nnodes,
               spatstat.linnet = spatstat.linnet)

  }else if(network=="simplenet"){
    data("simplenet")
    mesh = as.fdaPDE.spatstat.linnet(simplenet)
    mesh = refine.mesh.1.5D(mesh,0.025)
    
    FEMbasis = create.FEM.basis(mesh)
    nnodes = nrow(mesh$nodes)
    spatstat.linnet = as.spatstat.linnet.fdaPDE(mesh)
    
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

auxiliary_test2 = function(x, y, seg, tp, sigma= 0.125, 
                           nodes.lpp = ppp(x = mesh$nodes[,1], y = mesh$nodes[,2], 
                                           window = owin(xrange = c(min(mesh$nodes[,1]),max(mesh$nodes[,1])),
                                                         yrange = c(min(mesh$nodes[,2]),max(mesh$nodes[,2])))),
                           L = spatstat.linnet,
                           source = sources)
{ 
  PP = ppp(x = x, y = y, window = nodes.lpp$window)
  ND = crossdist.lpp(lpp(nodes.lpp, L), lpp(PP, L))
  
  return(     1./3 * 1/sqrt(2*pi*sigma^2) * exp(-ND[source[1],]^2/(2*sigma^2)) + 
                1./3 * 1/sqrt(2*pi*sigma^2) * exp(-ND[source[2],]^2/(2*sigma^2)) + 
                1./3 * 1/sqrt(2*pi*sigma^2) * exp(-ND[source[3],]^2/(2*sigma^2)) )
  
  
}

aux_test = c(auxiliary_test1, auxiliary_test2)
