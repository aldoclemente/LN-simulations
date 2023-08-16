# LN-DE-PDE: simulation loop ---------------------------------------------------

rmse.DE_PDE = rmse.KDE_PDE = rmse.KDE_ES = rmse.KDE_2D = rmse.VORONOI = 
              matrix(nrow=nsim, ncol = length(n))

DE_PDE.FEM = KDE_PDE.FEM = KDE_ES.FEM = KDE_2D.FEM = VORONOI.FEM = 
             matrix(0,nrow=nnodes,ncol= length(n))

for(j in 1:length(n)){
  cat(paste("###############  n = ", n[j], "  ###############\n", sep="") ) 
  for( i in 1:nsim){
    cat(paste("###############  ", i, " / ", nsim,"  ###############\n", sep="") )
    PP = rlpp(n=n[j], f = DENSITY)  
    data = cbind(PP$data$x, PP$data$y)
    
    # DE-PDE #
    if(methods[1]){
      start = Sys.time()
      invisible(capture.output( DE_PDE <-  fdaPDE::DE.FEM(data = data, 
                                                          FEMbasis = FEMbasis,
                                                          lambda = lambda[1]) ) )
      cat(paste("- DE-PDE DONE, time elapsed = ", difftime(Sys.time(),start, units = "mins")," mins \n", sep="") )
      #plot(DE_PDE$CV_err, main=paste("n = ", n[j], sep=""))
      rmse.DE_PDE[i,j] = sqrt(mean((true.density - eval.FEM(FEM(coeff=exp(DE_PDE$g),FEMbasis),
                                                            locations = locs.test))^2 ))
      
      DE_PDE.FEM[,j] = DE_PDE.FEM[,j] + exp(DE_PDE$g) / nsim
    }
    # spatstat returns the INTENSITY function. 
    # The estimation of DENSITY or INTENSITY is equivalent, if n is fixed,
    # INTENSITY(p) = n DENSITY(p) for all p. 
    
    # See 
    # McSwiggan, Greg, Adrian Baddeley, and Gopalan Nair. 
    # "Kernel density estimation on a linear network." 
    # Scandinavian Journal of Statistics 44.2 (2017): 324-345.
    
    # KDE-PDE
    if(methods[2]){  
      start = Sys.time()  
      invisible(capture.output(bw <- bw.lppl(X = PP) ))
      invisible(capture.output(KDE_PDE <- densityHeat(x = as.lpp(PP), sigma = as.numeric(bw), iterMax = 1e+9) )) 
      cat(paste("- KDE-PDE DONE, time elapsed = ", difftime(Sys.time(),start, units = "mins")," mins \n", sep="") )
      rmse.KDE_PDE[i,j] = sqrt(mean( (true.density - eval.FEM(FEM(coeff=as.linfun(KDE_PDE/n[j])(mesh$nodes[,1], mesh$nodes[,2]),
                                                                  FEMbasis),
                                                              locations = locs.test))^2 ))
    
      KDE_PDE.FEM[,j] = KDE_PDE.FEM[,j] + as.linfun(KDE_PDE/n[j])(mesh$nodes[,1], mesh$nodes[,2])/ nsim
    }
    
    # KDE-ES 
    if(methods[3]){
      start = Sys.time()
      invisible(capture.output(bw <- bw.lppl(X = PP) ))
      invisible(capture.output(KDE_ES <- densityEqualSplit(x = PP, sigma = bw) ))
      cat(paste("- KDE-ES DONE, time elapsed = ", difftime(Sys.time(),start, units = "mins")," mins \n", sep="") )
      
      rmse.KDE_ES[i,j] = sqrt(mean( (true.density - eval.FEM(FEM(coeff=as.linfun(KDE_ES/n[j])(mesh$nodes[,1], mesh$nodes[,2]),
                                                                 FEMbasis),
                                                             locations = locs.test))^2 ))
      KDE_ES.FEM[,j] = KDE_ES.FEM[,j] + as.linfun(KDE_ES.FEM/n[j])(mesh$nodes[,1], mesh$nodes[,2])/ nsim 
    }
    
    # KDE-2D
    if(methods[4]){
      start = Sys.time()
      invisible(capture.output(bw <- bw.scott(X = PP) ))
      invisible(capture.output(KDE_2D <-  densityQuick.lpp(X = PP, sigma = bw) )) #, at = points)
      cat(paste("- KDE-2D DONE, time elapsed = ", difftime(Sys.time(),start, units = "mins")," mins \n", sep="") )
      
      rmse.KDE_2D[i,j] = sqrt(mean( (true.density - eval.FEM(FEM(coeff=as.linfun(KDE_2D/n[j])(mesh$nodes[,1], mesh$nodes[,2]),
                                                                 FEMbasis),
                                                             locations = locs.test))^2 ))
      KDE_2D.FEM[,j] = KDE_2D.FEM[,j] + as.linfun(KDE_2D/n[j])(mesh$nodes[,1], mesh$nodes[,2])/ nsim
    }
    
    # KDE-VORONOI #
    if(methods[5]){
      start = Sys.time()
      invisible(capture.output(bw <- bw.voronoi(X = PP) ))
      invisible(capture.output(VORONOI <- densityVoronoi(X = PP, sigma = bw) ))
      cat(paste("- VORONOI DONE, time elapsed = ", difftime(Sys.time(),start, units = "mins")," mins \n", sep="") )
      
      rmse.VORONOI[i,j] = sqrt(mean( (true.density - eval.FEM(FEM(coeff=as.linfun(VORONOI/n[j])(mesh$nodes[,1], mesh$nodes[,2]),
                                                                  FEMbasis),
                                                              locations = locs.test))^2 ))
      VORONOI.FEM[,j] = KDE_2D.FEM[,j] + as.linfun(VORONOI/n[j])(mesh$nodes[,1], mesh$nodes[,2])/ nsim
    }
    cat(paste("###############  ############### ###############\n", sep="") )
    
    save(rmse.DE_PDE, rmse.KDE_PDE, rmse.KDE_ES,
         rmse.KDE_2D, rmse.VORONOI, i ,j, methods.names, methods, n,
         folder.name, date_, file = paste(folder.name,"RMSE",".RData", sep=""))
    
    save(DE_PDE.FEM, KDE_PDE.FEM, KDE_2D.FEM, VORONOI.FEM, KDE_ES.FEM, true.density.FEM,
         FEMbasis,
         file = paste(folder.name,"estimates",".RData", sep=""))
    
  }
}

