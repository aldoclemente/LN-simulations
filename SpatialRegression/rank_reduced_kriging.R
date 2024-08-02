
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

CovMat_RRKrig <- function(obs, net_dist, knots, cov_model="sph", # "exp" "sph" "gau" "cau"
                          predict_net_dist=NULL,
                          theta = c(log(2), log(15000), log(10000), log(.7))){ 
  
  
  i = match(cov_model, c("exp", "sph", "gau", "cau"))
  
  rrDist = net_dist[,rownames(knots)]
  knDist = as.matrix(dist(knots, diag = TRUE, upper = TRUE))
  
  CovMat = NULL  
  if(i == 1){
    # -------------------------------------------------------------------------
    #               Reduce Rank Exponential Model
    # -------------------------------------------------------------------------
    RRexpEst = optim(theta, m2LLrr, z = obs,  
                     rrDist = rrDist, knDist = knDist,
                     control = list(maxit = 650))
    param_estimates=list(sigmap = exp(RRexpEst$par)[1], 
                         alpha = exp(RRexpEst$par)[2], 
                         rho = exp(RRexpEst$par)[3],
                         sigma0 = exp(RRexpEst$par)[4])
    
    CovMat = param_estimates$sigmap^2*exp(-rrDist/param_estimates$alpha) %*% 
      solve(exp(-knDist/param_estimates$rho)) %*% t(exp(-rrDist/param_estimates$alpha)) + 
      diag(rep(param_estimates$sigma0^2, times = length(rrDist[,1])))
  }
  if(i==2){
    # -------------------------------------------------------------------------
    #               Reduce Rank Spherical Model
    # -------------------------------------------------------------------------
    RRsphEst = optim(theta, m2LLrr, z = obs,  
                     rrDist = rrDist, knDist = knDist, corMod = 'sph',
                     control = list(maxit = 1000))
    param_estimates=list(
          sigmap = exp(RRsphEst$par)[1],
          alpha = exp(RRsphEst$par)[2],
          rho = exp(RRsphEst$par)[3],
          sigma0 = exp(RRsphEst$par)[4])
        
    CovMat = param_estimates$sigmap^2*acor.sph(rrDist,param_estimates$alpha) %*% 
      solve(acor.sph(knDist,param_estimates$rho)) %*% 
      t(acor.sph(rrDist,param_estimates$alpha)) + 
      diag(rep(param_estimates$sigma0^2, times = length(rrDist[,1])))
  }
  if(i==3){
    # -------------------------------------------------------------------------
    #               Reduce Rank Gaussian Model
    # -------------------------------------------------------------------------
    RRgauEst = optim(theta, m2LLrr, z = obs,  
                     rrDist = rrDist, knDist = knDist, corMod = 'gau',
                     control = list(maxit = 1000))
    param_estimates=list(
                  sigmap = exp(RRgauEst$par)[1],
                  alpha = exp(RRgauEst$par)[2],
                  rho = exp(RRgauEst$par)[3],
                  sigma0 = exp(RRgauEst$par)[4])
    CovMat = param_estimates$sigmap^2*acor.gau(rrDist,param_estimates$alpha) %*% 
      solve(acor.gau(knDist,param_estimates$rho)) %*% 
      t(acor.gau(rrDist,param_estimates$alpha)) + 
      diag(rep(param_estimates$sigma0^2, times = length(rrDist[,1])))
  }
  if(i==4){
    # -------------------------------------------------------------------------
    #               Reduce Rank Cauchy Model
    # -------------------------------------------------------------------------
    RRcauEst = optim(theta, m2LLrr, z = obs,  
                     rrDist = rrDist, knDist = knDist, corMod = 'cau',
                     control = list(maxit = 1000))
    param_estimates=list(
            sigmap = exp(RRcauEst$par)[1],
            alpha = exp(RRcauEst$par)[2],
            rho = exp(RRcauEst$par)[3],
            sigma0 = exp(RRcauEst$par)[4])
    CovMat = param_estimates$sigmap^2*acor.cau(rrDist,param_estimates$alpha) %*% 
      solve(acor.cau(knDist,param_estimates$rho)) %*% 
      t(acor.cau(rrDist,param_estimates$alpha)) + 
      diag(rep(param_estimates$sigma0^2, times = length(rrDist[,1])))
  }

  if(is.null(predict_net_dist)){
    ret_ = list(CovMat = CovMat, param_estimates=param_estimates, 
                predict_CovMat = NULL, prediction = NULL,
                cov_model=c("exp", "sph", "gau", "cau")[i])
  }else{
    
    prediction = vector(mode="double", length = nrow(predict_net_dist))
    
    Rho = list(exp=ExpMod, sph=SphMod , gau=GauMod, cau=CauMod)
  
    #invC = solve(CovMat)
    # da modificare
    
    One_ = matrix(1, nrow=nrow(CovMat))
    c_ =  param_estimates$sigmap^2*Rho[[cov_model]](predict_net_dist, param_estimates$alpha)
    z_ = solve(CovMat,c_)
    w_ = solve(CovMat, One_)
    v_ = One_ %*% ( 1. - t(One_)%*%z_)/as.numeric((t(One_)%*%w_))
    b_ = ( c_ + v_ )

    lambda_ = solve(t(CovMat), b_)
    prediction = t(lambda_) %*% obs
    
    ret_ = list(CovMat = CovMat, param_estimates=param_estimates, 
                prediction = prediction,
                cov_model=c("exp", "sph", "gau", "cau")[i])
  }
  return(ret_)
}


# A = rbind(CovMat, rep(1, times=nrow(predict_net_dist)))
# A = cbind(A, c( rep(1, times=nrow(predict_net_dist)), 0))
# tmp = vector(mode="double", length = ncol(predict_net_dist))
# for(k in 1:nrow(predict_net_dist)){
#   
#   Gamma0 = param_estimates$sigmap^2*Rho[[cov_model]](predict_net_dist[,k], param_estimates$alpha) 
#   b = c(Gamma0, 1)
#   lam_ = solve(A,b)  
#   tmp[k] = t(lam_[1:length(obs)]) %*% obs
# }
# 
# range(tmp)
