
rank_reduced_kriging <- function(observations, 
                                    ND, knots,
                                    model=c(T,T,T,T)){
  
  rrDist = ND[,rownames(knots)]
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
    RRexpEst = optim(theta, m2LLrr, z = observations,  
                     rrDist = rrDist, knDist = knDist)
    sigmapRRexp = exp(RRexpEst$par)[1]
    alphaRRexp = exp(RRexpEst$par)[2]
    rhoRRexp = exp(RRexpEst$par)[3]
    sigma0RRexp = exp(RRexpEst$par)[4]
    SigRRexp = sigmapRRexp^2*exp(-rrDist/alphaRRexp) %*% 
      solve(exp(-knDist/rhoRRexp)) %*% t(exp(-rrDist/alphaRRexp)) + 
      diag(rep(sigma0RRexp^2, times = length(rrDist[,1])))
    RRexp = LOOCV(SigRRexp, observations)
  }
  if(model[2]){
    # -------------------------------------------------------------------------
    #               Reduce Rank Spherical Model
    # -------------------------------------------------------------------------
    theta = c(log(2), log(15000), log(10000), log(.7))
    RRsphEst = optim(theta, m2LLrr, z = observations,  
                     rrDist = rrDist, knDist = knDist, corMod = 'sph')
    sigmapRRsph = exp(RRsphEst$par)[1]
    alphaRRsph = exp(RRsphEst$par)[2]
    rhoRRsph = exp(RRsphEst$par)[3]
    sigma0RRsph = exp(RRsphEst$par)[4]
    SigRRsph = sigmapRRsph^2*acor.sph(rrDist,alphaRRsph) %*% 
      solve(acor.sph(knDist,rhoRRsph)) %*% 
      t(acor.sph(rrDist,alphaRRsph)) + 
      diag(rep(sigma0RRsph^2, times = length(rrDist[,1])))
    RRsph = LOOCV(SigRRsph,observations)
  }
  
  if(model[3]){
    # -------------------------------------------------------------------------
    #               Reduce Rank Gaussian Model
    # -------------------------------------------------------------------------
    theta = c(log(2), log(15000), log(10000), log(.7))
    RRgauEst = optim(theta, m2LLrr, z = observations,  
                     rrDist = rrDist, knDist = knDist, corMod = 'gau')
    sigmapRRgau = exp(RRgauEst$par)[1]
    alphaRRgau = exp(RRgauEst$par)[2]
    rhoRRgau = exp(RRgauEst$par)[3]
    sigma0RRgau = exp(RRgauEst$par)[4]
    SigRRgau = sigmapRRgau^2*acor.gau(rrDist,alphaRRgau) %*% 
      solve(acor.gau(knDist,rhoRRgau)) %*% 
      t(acor.gau(rrDist,alphaRRgau)) + 
      diag(rep(sigma0RRgau^2, times = length(rrDist[,1])))
    RRgau = LOOCV(SigRRgau,observations)
  }
  if(model[4]){
    # -------------------------------------------------------------------------
    #               Reduce Rank Cauchy Model
    # -------------------------------------------------------------------------
    theta = c(log(2), log(15000), log(10000), log(.7))
    RRcauEst = optim(theta, m2LLrr, z = observations,  
                     rrDist = rrDist, knDist = knDist, corMod = 'cau')
    sigmapRRcau = exp(RRcauEst$par)[1]
    alphaRRcau = exp(RRcauEst$par)[2]
    rhoRRcau = exp(RRcauEst$par)[3]
    sigma0RRcau = exp(RRcauEst$par)[4]
    SigRRcau = sigmapRRcau^2*acor.cau(rrDist,alphaRRcau) %*% 
      solve(acor.cau(knDist,rhoRRcau)) %*% 
      t(acor.cau(rrDist,alphaRRcau)) + 
      diag(rep(sigma0RRcau^2, times = length(rrDist[,1])))
    RRcau = LOOCV(SigRRcau,observations)
  }
  
  ret_ = list(RRexp = RRexp,
              RRsph = RRsph,
              RRgau = RRgau,
              RRcau = RRcau)
  
  return(ret_)
}