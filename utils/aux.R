aux = function(x, y, seg, tp, sigma= 0.125, 
                           nodes.lpp = ppp(x = mesh$nodes[,1], y = mesh$nodes[,2], 
                                           window = owin(xrange = c(min(mesh$nodes[,1]),max(mesh$nodes[,1])),
                                                         yrange = c(min(mesh$nodes[,2]),max(mesh$nodes[,2])))),
                           L = spatstat.linnet,
                           source = sources)
{ 
  PP = ppp(x = x, y = y, window = nodes.lpp$window)
  ND = crossdist.lpp(lpp(nodes.lpp, L), lpp(PP, L))
  
  res = 0.
  for( k in 1:length(source)){
    res = res + 1./(2*length(source)) * 1./sqrt(2*pi*sigma^2) * exp(-ND[source[k],]^2/(2*sigma^2))
  }
  
  return(res)
  #return(   0.25 * 1/sqrt(2*pi*sigma^2) * exp(-ND[source[1],]^2/(2*sigma^2)) + 
  #            0.25 * 1/sqrt(2*pi*sigma^2) * exp(-ND[source[2],]^2/(2*sigma^2)))
  
  
}

cluster = function(x, y, seg, tp, sigma= 0.065, # 0.075
                   nodes.lpp = ppp(x = mesh$nodes[,1], y = mesh$nodes[,2], 
                                   window = owin(xrange = c(min(mesh$nodes[,1]),max(mesh$nodes[,1])),
                                                 yrange = c(min(mesh$nodes[,2]),max(mesh$nodes[,2])))),
                   L = spatstat.linnet)
{ 
  #set.seed(314156)
  vertices = cbind(simplenet$vertices$x, simplenet$vertices$y) 
  centers = matrix(nrow=0,ncol=2)
  for(i in 1:length(simplenet$from)){
    centers = rbind(centers, c(0,0))
    centers[i,1] = mean(c(simplenet$lines$ends[i,1], simplenet$lines$ends[i,3]))
    centers[i,2] = mean(c(simplenet$lines$ends[i,2], simplenet$lines$ends[i,4]))
  }
    
  lengths = matrix(0,nrow=nrow(simplenet$lines$ends), ncol=1)
  for( i in 1:nrow(lengths)){
    lengths[i] = crossdist(simplenet$lines$ends[i,1], simplenet$lines$ends[i,2],
                           simplenet$lines$ends[i,3], simplenet$lines$ends[i,4])
  }
  #runiflpp(30, spatstat.linnet)
  #source = cbind(centers$data$x, centers$data$y)
  PP = ppp(x = x, y = y, window = nodes.lpp$window)
  ND = crossdist.lpp(lpp(centers, L), lpp(PP, L))
  
  ND2 = crossdist.lpp(lpp(vertices, L), lpp(PP,L))
  
  res = 0
  for(i in 1:nrow(centers)){
    res = res + lengths[i] * 1/sqrt(2*pi*sigma^2) * exp(-ND[i,]^2/(2*sigma^2))
  }
  res = res + 0.125 * 1/sqrt(2*pi*sigma^2) * exp(-ND2[8,]^2/(2*sigma^2))
  
  return( res )
}
