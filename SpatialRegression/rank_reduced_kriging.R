
rank_reduced_kriging <- function(obs, 
                                    net_dist, knots,
                                    model=c(T,T,T,T)){
  
  rrDist = net_dist[,rownames(knots)]
  knDist = as.matrix(dist(knots, diag = TRUE, upper = TRUE))
  
  RRexp = NULL
  RRsph = NULL
  RRgau = NULL
  RRcau = NULL
  if(model[1]){
    # -------------------------------------------------------------------------
    #               Reduce Rank Exponential Model
    # -------------------------------------------------------------------------
    theta = c(log(2), log(15000), log(10000), log(.7))
    RRexpEst = optim(theta, m2LLrr, z = obs,  
                     rrDist = rrDist, knDist = knDist)
    sigmapRRexp = exp(RRexpEst$par)[1]
    alphaRRexp = exp(RRexpEst$par)[2]
    rhoRRexp = exp(RRexpEst$par)[3]
    sigma0RRexp = exp(RRexpEst$par)[4]
    SigRRexp = sigmapRRexp^2*exp(-rrDist/alphaRRexp) %*% 
      solve(exp(-knDist/rhoRRexp)) %*% t(exp(-rrDist/alphaRRexp)) + 
      diag(rep(sigma0RRexp^2, times = length(rrDist[,1])))
    RRexp = LOOCV(SigRRexp, obs)
  }
  if(model[2]){
    # -------------------------------------------------------------------------
    #               Reduce Rank Spherical Model
    # -------------------------------------------------------------------------
    theta = c(log(2), log(15000), log(10000), log(.7))
    RRsphEst = optim(theta, m2LLrr, z = obs,  
                     rrDist = rrDist, knDist = knDist, corMod = 'sph')
    sigmapRRsph = exp(RRsphEst$par)[1]
    alphaRRsph = exp(RRsphEst$par)[2]
    rhoRRsph = exp(RRsphEst$par)[3]
    sigma0RRsph = exp(RRsphEst$par)[4]
    SigRRsph = sigmapRRsph^2*acor.sph(rrDist,alphaRRsph) %*% 
      solve(acor.sph(knDist,rhoRRsph)) %*% 
      t(acor.sph(rrDist,alphaRRsph)) + 
      diag(rep(sigma0RRsph^2, times = length(rrDist[,1])))
    RRsph = LOOCV(SigRRsph,obs)
  }
  
  if(model[3]){
    # -------------------------------------------------------------------------
    #               Reduce Rank Gaussian Model
    # -------------------------------------------------------------------------
    theta = c(log(2), log(15000), log(10000), log(.7))
    RRgauEst = optim(theta, m2LLrr, z = obs,  
                     rrDist = rrDist, knDist = knDist, corMod = 'gau')
    sigmapRRgau = exp(RRgauEst$par)[1]
    alphaRRgau = exp(RRgauEst$par)[2]
    rhoRRgau = exp(RRgauEst$par)[3]
    sigma0RRgau = exp(RRgauEst$par)[4]
    SigRRgau = sigmapRRgau^2*acor.gau(rrDist,alphaRRgau) %*% 
      solve(acor.gau(knDist,rhoRRgau)) %*% 
      t(acor.gau(rrDist,alphaRRgau)) + 
      diag(rep(sigma0RRgau^2, times = length(rrDist[,1])))
    RRgau = LOOCV(SigRRgau,obs)
  }
  if(model[4]){
    # -------------------------------------------------------------------------
    #               Reduce Rank Cauchy Model
    # -------------------------------------------------------------------------
    theta = c(log(2), log(15000), log(10000), log(.7))
    RRcauEst = optim(theta, m2LLrr, z = obs,  
                     rrDist = rrDist, knDist = knDist, corMod = 'cau')
    sigmapRRcau = exp(RRcauEst$par)[1]
    alphaRRcau = exp(RRcauEst$par)[2]
    rhoRRcau = exp(RRcauEst$par)[3]
    sigma0RRcau = exp(RRcauEst$par)[4]
    SigRRcau = sigmapRRcau^2*acor.cau(rrDist,alphaRRcau) %*% 
      solve(acor.cau(knDist,rhoRRcau)) %*% 
      t(acor.cau(rrDist,alphaRRcau)) + 
      diag(rep(sigma0RRcau^2, times = length(rrDist[,1])))
    RRcau = LOOCV(SigRRcau,obs)
  }
  
  ret_ = list(RRexp = RRexp,
              RRsph = RRsph,
              RRgau = RRgau,
              RRcau = RRcau)
  
  return(ret_)
}

CovMat_RRKrig <- function(obs, net_dist, knots, cov_model="sph",
                          predict_net_dist=NULL){ # "exp" "sph" "gau" "cau" 
  
  
  i = match(cov_model, c("exp", "sph", "gau", "cau"))
  
  rrDist = net_dist[,rownames(knots)]
  knDist = as.matrix(dist(knots, diag = TRUE, upper = TRUE))
  
  if(!is.null(predict_net_dist)){
    predict_rrDist = predict_net_dist[,rownames(knots)]
  }
  
  CovMat = NULL  
  if(i == 1){
    # -------------------------------------------------------------------------
    #               Reduce Rank Exponential Model
    # -------------------------------------------------------------------------
    theta = c(log(2), log(15000), log(10000), log(.7))
    RRexpEst = optim(theta, m2LLrr, z = obs,  
                     rrDist = rrDist, knDist = knDist)
    param_estimates=list(sigmap = exp(RRexpEst$par)[1], 
                         alpha = exp(RRexpEst$par)[2], 
                         rho = exp(RRexpEst$par)[3],
                         sigma0 = exp(RRexpEst$par)[4])
    
    CovMat = param_estimates$sigmap^2*exp(-rrDist/param_estimates$alpha) %*% 
      solve(exp(-knDist/param_estimates$rho)) %*% t(exp(-rrDist/param_estimates$alpha)) + 
      diag(rep(param_estimates$sigma0^2, times = length(rrDist[,1])))
    if(!is.null(predict_net_dist)){
      predict_CovMat = param_estimates$sigmap^2*exp(-predict_rrDist/param_estimates$alpha) %*% 
        solve(exp(-knDist/param_estimates$rho)) %*% t(exp(-rrDist/param_estimates$alpha)) 
    #  + diag(rep(param_estimates$sigma0^2, times = length(predict_rrDist[,1]))) 
    }
  }
  if(i==2){
    # -------------------------------------------------------------------------
    #               Reduce Rank Spherical Model
    # -------------------------------------------------------------------------
    theta = c(log(2), log(15000), log(10000), log(.7))
    RRsphEst = optim(theta, m2LLrr, z = obs,  
                     rrDist = rrDist, knDist = knDist, corMod = 'sph')
    param_estimates=list(
          sigmap = exp(RRsphEst$par)[1],
          alpha = exp(RRsphEst$par)[2],
          rho = exp(RRsphEst$par)[3],
          sigma0 = exp(RRsphEst$par)[4])
    
    CovMat = param_estimates$sigmap^2*acor.sph(rrDist,param_estimates$alpha) %*% 
      solve(acor.sph(knDist,param_estimates$rho)) %*% 
      t(acor.sph(rrDist,param_estimates$alpha)) + 
      diag(rep(param_estimates$sigma0^2, times = length(rrDist[,1])))
    if(!is.null(predict_net_dist)){
      predict_CovMat = param_estimates$sigmap^2*acor.sph(predict_rrDist,param_estimates$alpha) %*% 
        solve(acor.sph(knDist,param_estimates$rho)) %*% 
        t(acor.sph(rrDist,param_estimates$alpha)) #+ 
        #diag(rep(param_estimates$sigma0^2, times = length(predict_rrDist[,1]))) 
    }
  }
  
  if(i==3){
    # -------------------------------------------------------------------------
    #               Reduce Rank Gaussian Model
    # -------------------------------------------------------------------------
    theta = c(log(2), log(15000), log(10000), log(.7))
    RRgauEst = optim(theta, m2LLrr, z = obs,  
                     rrDist = rrDist, knDist = knDist, corMod = 'gau')
    param_estimates=list(
                  sigmap = exp(RRgauEst$par)[1],
                  alpha = exp(RRgauEst$par)[2],
                  rho = exp(RRgauEst$par)[3],
                  sigma0 = exp(RRgauEst$par)[4])
    CovMat = param_estimates$sigmap^2*acor.gau(rrDist,param_estimates$alpha) %*% 
      solve(acor.gau(knDist,param_estimates$rho)) %*% 
      t(acor.gau(rrDist,param_estimates$alpha)) + 
      diag(rep(param_estimates$sigma0^2, times = length(rrDist[,1])))
    if(!is.null(predict_net_dist)){
      predict_CovMat = param_estimates$sigmap^2*acor.gau(predict_rrDist,param_estimates$alpha) %*% 
        solve(acor.gau(knDist,param_estimates$rho)) %*% 
        t(acor.gau(rrDist,param_estimates$alpha)) 
      #+ 
      #  diag(rep(param_estimates$sigma0^2, times = length(predict_rrDist[,1]))) 
    }
  }
  if(i==4){
    # -------------------------------------------------------------------------
    #               Reduce Rank Cauchy Model
    # -------------------------------------------------------------------------
    theta = c(log(2), log(15000), log(10000), log(.7))
    RRcauEst = optim(theta, m2LLrr, z = obs,  
                     rrDist = rrDist, knDist = knDist, corMod = 'cau')
    param_estimates=list(
            sigmap = exp(RRcauEst$par)[1],
            alpha = exp(RRcauEst$par)[2],
            rho = exp(RRcauEst$par)[3],
            sigma0 = exp(RRcauEst$par)[4])
    CovMat = param_estimates$sigmap^2*acor.cau(rrDist,param_estimates$alpha) %*% 
      solve(acor.cau(knDist,param_estimates$rho)) %*% 
      t(acor.cau(rrDist,param_estimates$alpha)) + 
      diag(rep(param_estimates$sigma0^2, times = length(rrDist[,1])))
    if(!is.null(predict_net_dist)){
      predict_CovMat = param_estimates$sigmap^2*acor.cau(predict_rrDist,param_estimates$alpha) %*% 
        solve(acor.cau(knDist,param_estimates$rho)) %*% 
        t(acor.cau(predict_rrDist,param_estimates$alpha)) + 
        diag(rep(param_estimates$sigma0^2, times = length(predict_rrDist[,1]))) 
    }
  }
  
  
  
  if(is.null(predict_net_dist)){
    ret_ = list(CovMat = CovMat, param_estimates=param_estimates, 
                predict_CovMat = NULL, prediction = NULL,
                cov_model=c("exp", "sph", "gau", "cau")[i])
  }else{
    
    prediction = vector(mode="double", length = nrow(predict_net_dist))
    for(i in 1:nrow(predict_CovMat)){
      k = predict_CovMat[i,]
      w = solve(CovMat, k)
      prediction[i] = t(w) %*% obs 
    }
  
    ret_ = list(CovMat = CovMat, param_estimates=param_estimates,
                predict_CovMat = predict_CovMat, prediction = prediction,
                cov_model=c("exp", "sph", "gau", "cau")[i])
  }
  return(ret_)
}
