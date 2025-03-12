if(!require(pacman)) install.packages("pacman")
pacman::p_load("ggplot2", "ggforce", "viridis")

# performing n_sim time length(n_obs) simulations!!
.SimulationObjectCtr <- setRefClass("SimulationObject", 
                                    fields = c(method_name= "character", 
                                               n_obs="vector",      # number of observations
                                               n_sim="integer",      # number of repetitions
                                               estimates = "list",   # list of FEMs
                                               meanField = "ANY",          # mean estimated field
                                               y_hat = "list",
                                               mean_Y_HAT = "ANY",
                                               errors = "vector",        # error vector (RMSE)
                                               FEMbasis = "ANY"     
                                    ), 
                                    methods = list(
                                      compute_mean_field = function(j){
                                        coef <- rep(0, times=nrow(FEMbasis$mesh$nodes))
                                        for(i in 1:n_sim){
                                          coef <- coef + estimates[[(((j-1)*n_sim+i))]]$coeff / n_sim
                                        }
                                        meanField[[j]] <<- fdaPDE::FEM(coef, FEMbasis)
                                      },
                                      update_estimate = function(estimate,i,j){
                                        estimates[[(((j-1)*n_sim+i))]] <<- estimate
                                      },
                                      update_y_hat = function(vec, i, j){
                                        y_hat[[(((j-1)*n_sim+i))]] <<- vec
                                      },
                                      compute_mean_y_hat = function(j){
                                        coef <- rep(0, times=n_obs[j])
                                        for(i in 1:n_sim){
                                          coef <- coef + y_hat[[(((j-1)*n_sim+i))]] / n_sim
                                        }
                                        mean_Y_HAT[[j]] <<- coef
                                      },
                                      plot_mean_field = function(j, ...){
                                        
                                        plot(meanField[[j]], ...) + scale_color_viridis()
                                      })
                                    
)

.SimulationObjectCtrHelper <- function(method_name, n_obs, n_sim, FEMbasis) {
  num_nodes <- nrow(FEMbasis$mesh$nodes)
  nullFEM <- fdaPDE::FEM(rep(NA,times=num_nodes), FEMbasis)
  estimates <- rep(list(), length=n_sim*length(n_obs))
  y_hat <- rep(list(), length=n_sim*length(n_obs))
  errors <- rep(NA, times=n_sim*length(n_obs))
  meanField <- rep(list(), length=length(n_obs))
  mean_Y_HAT <- rep(list(), length=length(n_obs))
  return(list(
    method_name=method_name, n_obs=n_obs, n_sim=n_sim,
    estimates = estimates, meanField=meanField, errors=errors,FEMbasis=FEMbasis)
  )
}
setGeneric("Simulation", function(method_name, n_obs, n_sim, FEMbasis) standardGeneric("Simulation"))
setMethod("Simulation", signature=c(method_name="character",n_obs="vector",n_sim="integer", FEMbasis="ANY"),
          function(method_name,n_obs,n_sim, FEMbasis){
            aux_list = .SimulationObjectCtrHelper(method_name,n_obs,n_sim, FEMbasis)
            return(.SimulationObjectCtr(
                       method_name=aux_list$method_name, n_obs=aux_list$n_obs, n_sim=aux_list$n_sim,
                       estimates = aux_list$estimates, meanField=aux_list$meanField, errors=aux_list$errors,
                       FEMbasis=aux_list$FEMbasis))
          })

setMethod("boxplot", "SimulationObject", function(x,...){
  
  data_plot <- data.frame(errors = x$errors, 
                          n_obs=as.character(rep(x$n_obs, each=n_sim)))
  
  p<-ggplot(data_plot)+
    geom_boxplot(aes(x=n_obs, y=errors, group=n_obs, ...))+
    scale_x_discrete(limits=as.character(n_obs))+
    labs(x="", y="")+
    theme(axis.ticks.x = element_blank())
  p
  
})

### Density Estimation Simulation specialization ###
.DensityEstimationSimulationObjectCtr <- setRefClass("DensityEstimationSimulationObject", contains = "SimulationObject",
                                                     methods=list(
                                                       update_error = function(true_field, test_locs, i, j){
                                                         estimated <- fdaPDE::eval.FEM(FEM=estimates[[(((j-1)*n_sim+i))]], locations = test_locs)
                                                         errors[(((j-1)*n_sim+i))] <<- mean( (true_field - estimated)^2)
                                                       }
                                                     ))

setGeneric("DensityEstimationSimulation", function(method_name, n_obs, n_sim, FEMbasis) standardGeneric("DensityEstimationSimulation"))
setMethod("DensityEstimationSimulation", signature=c(method_name="character",n_obs="vector",n_sim="integer", FEMbasis="ANY"),
          function(method_name, n_obs, n_sim, FEMbasis){
            aux_list = .SimulationObjectCtrHelper(method_name,n_obs,n_sim, FEMbasis)
            .DensityEstimationSimulationObjectCtr(method_name=aux_list$method_name, n_obs=aux_list$n_obs, n_sim=aux_list$n_sim,
                                                  estimates = aux_list$estimates, meanField=aux_list$meanField, errors=aux_list$errors,
                                                  FEMbasis=aux_list$FEMbasis)
  
})

### SpatialRegression Simulation specialization ###
.SpatialRegressionSimulationObjectCtr <- setRefClass("SpatialRegressionSimulationObject", contains = "SimulationObject",
                                                     methods=list(
                                                       update_error = function(y_hat, y_true,i, j){
                                                         errors[(((j-1)*n_sim+i))] <<- mean((y_hat - y_true)^2)
                                                       }
                                                     ))

setGeneric("SpatialRegressionSimulation", function(method_name, n_obs, n_sim, FEMbasis) standardGeneric("SpatialRegressionSimulation"))
setMethod("SpatialRegressionSimulation", signature=c(method_name="character",n_obs="vector",n_sim="integer", FEMbasis="ANY"),
          function(method_name, n_obs, n_sim, FEMbasis){
            aux_list = .SimulationObjectCtrHelper(method_name,n_obs,n_sim, FEMbasis)
            .SpatialRegressionSimulationObjectCtr(method_name=aux_list$method_name, n_obs=aux_list$n_obs, n_sim=aux_list$n_sim,
                                                  estimates = aux_list$estimates, meanField=aux_list$meanField, errors=aux_list$errors,
                                                  FEMbasis=aux_list$FEMbasis)
            
          })


.BlockSimulationCtr <- setRefClass("BlockSimulation", 
                                   fields = c(Simulations = "list",
                                              num_methods = "integer",
                                              num_sim = "integer",
                                              n_obs ="vector",
                                              length_obs = "integer",
                                              method_names ="character",
                                              results ="data.frame"
                                              ))

.BlockSimulationCtrHelper <- function(x){
  num_methods = length(x)
  num_sim = x[[1]]$n_sim
  n_obs = x[[1]]$n_obs
  length_obs = length(n_obs)
  method_names = vector(mode="character", length=num_methods)
  # for each method, num_sim * length_obs errors
  colNames = c("errors","n_obs","method")
  results = matrix(NA,nrow=num_methods*num_sim*length_obs,ncol=length(colNames))
  
  for(meth in 1:num_methods){
    tmp<- cbind(x[[meth]]$errors, 
                as.character(rep(x[[meth]]$n_obs, 
                                 each=num_sim)),
                as.character(rep(x[[meth]]$method_name,
                                 times=num_sim*length_obs)))
    
    results[(num_sim*length_obs*(meth-1) + 1):(num_sim*length_obs*(meth)),] = tmp
    
    method_names[meth] <- x[[meth]]$method_name
    
  }
  
  colnames(results) <- colNames
  results <- as.data.frame(results)
  results[,1] <- as.numeric(results[,1])
  return(list(Simulations=x, num_methods=num_methods,
              num_sim=num_sim,n_obs=n_obs,length_obs=length_obs,
              method_names=method_names, results = results))
}

setGeneric("BlockSimulation", function(x) standardGeneric("BlockSimulation"))
setMethod("BlockSimulation",signature=c(x="list"),
          function(x){
            aux_list <- .BlockSimulationCtrHelper(x)
            .BlockSimulationCtr(Simulations=aux_list$Simulations, num_methods=aux_list$num_methods,
                                num_sim=aux_list$num_sim,n_obs=aux_list$n_obs,length_obs=aux_list$length_obs,
                                method_names=aux_list$method_names, results = aux_list$results)
          }
)

setMethod("boxplot", "BlockSimulation", function(x, ORDER=NULL ,...){
  
  if(!is.null(ORDER)){
  ORDER <- x$method_names[ORDER]
  x$results$method <- factor(x$results$method, levels= ORDER)
  
  begin=0.25
  end=0.95
  border_col = darken(viridis(length(x$method_names), begin=begin,end=end), amount=0.25)
  fill_col = viridis(length(x$method_names), begin=begin, end=end)
  BORDER = c()
  FILL = c()
  for(i in 1:length(x$method_names)){
    FILL = append(FILL, fill_col[i])
    BORDER = append(BORDER, border_col[i])
  }
  ggFILL <-scale_fill_manual(values = FILL) 
  ggBORDER <- scale_color_manual(values= BORDER) 
  }
  
  p<-ggplot(x$results)+
    geom_boxplot(aes(x=n_obs,
                     y=errors, group=interaction(method,n_obs),
                     fill=method, color = method))+
    scale_x_discrete(limits=as.character(x$n_obs))+
    labs(x="", y="") +
    theme(
      axis.ticks.x = element_blank(),
      legend.position = c(0.95,0.95), 
      legend.background = element_rect(fill="white", color="black",
                                       linewidth =c(1,0.5)),
      legend.title = element_blank())
  if(!is.null(ORDER)){
  p <- p + ggFILL + ggBORDER
  }
  p
})

