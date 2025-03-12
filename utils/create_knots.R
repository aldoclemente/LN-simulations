
create_knots <- function(locations, seed=NULL, nknots=round(nrow(locations)/2) ){
  if(!is.null(seed))
    set.seed(seed)
  
  nlocs = nrow(locations)
  locs = cbind(1:nlocs, locations)
  km = kmeans(locs[,2:3],nknots,40)
  dall = as.matrix(dist(rbind(locs[,2:3],km$centers), diag = TRUE, upper = TRUE))
  mins = apply(dall[1:nlocs,(nlocs+1):(nlocs+nknots)],2,min)
  knots = NULL
  for(i in 1:nknots)
    knots = rbind(knots, locs[which(dall[1:nlocs,nlocs + i] == mins[i])[1],])
  
  rownames(knots) = as.character(knots[,1])
  
  ret_knots = knots[,2:3] 
  return(ret_knots)
}

