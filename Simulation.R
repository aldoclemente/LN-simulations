library(ggplot2)
library(ggforce)
library(viridis)

# performing n_sim time length(n_obs) simulations!!
.SimulationObjectCtr <- setRefClass("SimulationObject", 
                                    fields = c(method_name= "character", 
                                               n_obs="vector",      # number of observations
                                               n_sim="integer",      # number of repetitions
                                               estimates = "list",   # list of FEMs
                                               meanField = "ANY",          # mean estimated field        
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
                                      update_error = function(true_field, test_locs, i, j){
                                        estimated <- fdaPDE::eval.FEM(FEM=estimates[[(((j-1)*n_sim+i))]], locations = test_locs)
                                        errors[(((j-1)*n_sim+i))] <<- mean( (true_field - estimated)^2)
                                      },
                                      plot_mean_field = function(j, ...){
                                        
                                        plot(meanField[[j]], ...) + scale_color_viridis()
                                      })
)

setGeneric("Simulation", function(method_name, n_obs, n_sim, FEMbasis) standardGeneric("Simulation"))
setMethod("Simulation", signature=c(method_name="character",n_obs="vector",n_sim="integer", FEMbasis="ANY"),
          function(method_name,n_obs,n_sim, FEMbasis){
            num_nodes <- nrow(FEMbasis$mesh$nodes)
            nullFEM <- fdaPDE::FEM(rep(NA,times=num_nodes), FEMbasis)
            estimates <- rep(list(), length=n_sim*length(n_obs))
            errors <- rep(NA, times=n_sim*length(n_obs))
            meanField <- rep(list(), length=length(n_obs))
            return(.SimulationObjectCtr(
                       method_name=method_name, n_obs=n_obs, n_sim=n_sim,
                       estimates = estimates, meanField=meanField, errors=errors,FEMbasis=FEMbasis))
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

.BlockSimulationCtr <- setRefClass("BlockSimulation", 
                                   fields = c(Simulations = "list",
                                              num_methods = "integer",
                                              num_sim = "integer",
                                              n_obs ="integer",
                                              length_obs = "integer",
                                              method_names ="character",
                                              results ="data.frame"
                                              ))

setGeneric("BlockSimulation", function(x) standardGeneric("BlockSimulation"))
setMethod("BlockSimulation",signature=c(x="list"),
           function(x){
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
             return(.BlockSimulationCtr(Simulations=x, num_methods=num_methods,
                                        num_sim=num_sim,n_obs=n_obs,length_obs=length_obs,
                                        method_names=method_names, results = results))
           }
)

setMethod("boxplot", "BlockSimulation", function(x,...){
  
  # data
  p<-ggplot(x$results)+
    geom_boxplot(aes(x=n_obs,
                     y=errors, group=interaction(method,n_obs),
                     fill=method, color = method))+
    scale_x_discrete(limits=as.character(n_obs))+
    labs(x="", y="") +
    theme(
      axis.ticks.x = element_blank(),
      legend.position = c(0.85,0.85), 
      legend.background = element_rect(fill="white", color="black",
                                       linewidth =c(1,0.5)),
      legend.title = element_blank())
  p  
  
})

