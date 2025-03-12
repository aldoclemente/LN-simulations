if(!require(pacman)) install.packages("pacman")
pacman::p_load("ggplot2", "ggforce", "viridis")

.KfoldsCrossValidationCtr <- setRefClass("KfoldsCrossValidationObject",
                                         fields = c(K="integer",       # number of folds
                                                      data="data.frame", # raw data
                                                      folds="list",      # k Folds
                                                      num_obs_kFold="integer" # num_obs_kFold[i] = nrow(folds[[i]])
                                                      ),      
                                         methods = list(
                                           get_Kfold = function(j){
                                             folds[[j]]
                                           },
                                           get_data = function(j){
                                             train_data = list()
                                             for(i in 1:K){
                                               if( i == j){
                                                 test_data = folds[[i]]
                                               }else{
                                                 train_data = rbind(train_data, folds[[i]])
                                                 
                                               }
                                             }
                                             return(list(train_data = train_data, test_data = test_data))
                                           }
                                         ))

setGeneric("KfoldsCrossValidation", function(data, seed, K) standardGeneric("KfoldsCrossValidation"))
setMethod("KfoldsCrossValidation", signature=c(data="data.frame", seed="integer", K="integer"),
          function(data,seed=27182L,K=10L){
            K <- as.integer(K)
            seed <- as.integer(seed)
            folds <- list()
            num_obs_kFold <- vector(mode="integer", length=K)
            data = data[sample(1:nrow(data)), ]
            num_data = round(nrow(data)/K)
            num_obs_kFold[1:(K-1)] <- rep(num_data, times=(K-1))
            
            for(i in 1:(K-1)){
              folds[[i]] = data[(1 + num_data*(i-1)):(num_data*i),]
              
            }
            folds[[K]] = data[(num_data*(K-1) + 1):nrow(data), ]
            num_obs_kFold[K] <- nrow(folds[[K]])
            storage.mode(num_obs_kFold) <- "integer" 
            .KfoldsCrossValidationCtr(K=K, data=data, folds=folds, num_obs_kFold=num_obs_kFold)
          })

.CaseStudyObjectCtr <- setRefClass("CaseStudyObject", contains = "SimulationObject",
                             methods = list(
                               update_estimate = function(estimate,j){
                                 i=1
                                 estimates[[(((j-1)*n_sim+i))]] <<- estimate
                               }
                             ))

setGeneric("CaseStudy", function(method_name, n_obs, FEMbasis) standardGeneric("CaseStudy"))
setMethod("CaseStudy", signature=c(method_name="character",n_obs="vector", FEMbasis="ANY"),
          function(method_name, n_obs, FEMbasis){
            aux_list = .SimulationObjectCtrHelper(method_name,n_obs,1L, FEMbasis)
            return(.CaseStudyObjectCtr(
              method_name=aux_list$method_name, n_obs=aux_list$n_obs, n_sim=aux_list$n_sim,
              estimates = aux_list$estimates, meanField=aux_list$meanField, errors=aux_list$errors,
              FEMbasis=aux_list$FEMbasis))
})

setMethod("boxplot", "CaseStudyObject", function(x,...){
  
  data_plot <- data.frame(errors = x$errors, method_name = rep(x$method_name, times=length(x$n_obs)))
                          
  p<-ggplot(data_plot)+
    geom_boxplot(aes(x=method_name, y=errors, group=method_name, ...))+
    scale_x_discrete(limits=x$method_name)+
    labs(x="", y="")+
    theme(axis.ticks.x = element_blank())
  p
})

# Density Estimation specialization (CROSS VALIDATION ERROR!).
.DensityEstimationCaseStudyObjectCtr <-setRefClass("DensityEstimationCaseStudyObject",
                                                   contains = "CaseStudyObject",
                                                   methods=list(
                                                     update_error = function(test_locs, j){
                                                       i=1
                                                       eval_locs <- fdaPDE::refine.by.splitting.mesh.1.5D(FEMbasis$mesh)$nodes
                                                       coef <- estimates[[(((j-1)*n_sim+i))]]$coeff^2
                                                       errors[(((j-1)*n_sim+i))] <<- mean(fdaPDE::eval.FEM(FEM(coef,FEMbasis), locations=eval_locs)) - 2*mean(fdaPDE::eval.FEM(estimates[[(((j-1)*n_sim+i))]],locations=as.matrix(test_locs)))
                                                     }
                                                   )
)

setGeneric("DensityEstimationCaseStudy", function(method_name, n_obs, FEMbasis) standardGeneric("DensityEstimationCaseStudy"))
setMethod("DensityEstimationCaseStudy", signature=c(method_name="character",n_obs="vector", FEMbasis="ANY"),
          function(method_name, n_obs, FEMbasis){
            aux_list = .SimulationObjectCtrHelper(method_name,n_obs,1L, FEMbasis)
            return(.DensityEstimationCaseStudyObjectCtr(
              method_name=aux_list$method_name, n_obs=aux_list$n_obs, n_sim=aux_list$n_sim,
              estimates = aux_list$estimates, meanField=aux_list$meanField, errors=aux_list$errors,
              FEMbasis=aux_list$FEMbasis))
})

# Spatial Regression Specialization (ERROR ...)
.SpatialRegressionCaseStudyObjectCtr <-setRefClass("SpatialRegressionCaseStudyObject",
                                                   contains = "CaseStudyObject",
                                                   methods=list(
                                                     update_error = function(y_hat, y_true, j){
                                                       i=1
                                                       errors[(((j-1)*n_sim+i))] <<- mean((y_hat - y_true)^2)
                                                     }
                                                   )
)

setGeneric("SpatialRegressionCaseStudy", function(method_name, n_obs, FEMbasis) standardGeneric("SpatialRegressionCaseStudy"))
setMethod("SpatialRegressionCaseStudy", signature=c(method_name="character",n_obs="vector", FEMbasis="ANY"),
          function(method_name, n_obs, FEMbasis){
            aux_list = .SimulationObjectCtrHelper(method_name,n_obs,1L, FEMbasis)
            return(.SpatialRegressionCaseStudyObjectCtr(
              method_name=aux_list$method_name, n_obs=aux_list$n_obs, n_sim=aux_list$n_sim,
              estimates = aux_list$estimates, meanField=aux_list$meanField, errors=aux_list$errors,
              FEMbasis=aux_list$FEMbasis))
          })

.BlockCaseStudyObjectCtr <- setRefClass("BlockCaseStudyObject", contains="BlockSimulation")

setGeneric("BlockCaseStudy", function(x) standardGeneric("BlockCaseStudy"))
setMethod("BlockCaseStudy",signature=c(x="list"),
          function(x){
             aux_list <- .BlockSimulationCtrHelper(x)
            .BlockCaseStudyObjectCtr(Simulations=aux_list$Simulations, num_methods=aux_list$num_methods,
                                     num_sim=aux_list$num_sim,n_obs=aux_list$n_obs,length_obs=aux_list$length_obs,
                                     method_names=aux_list$method_names, results = aux_list$results)
          }
)

setMethod("boxplot", "BlockCaseStudyObject", function(x, ORDER=NULL, ...){
  
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
  
  # data
  p<-ggplot(x$results)+
    geom_boxplot(aes(x=method,
                     y=errors, group=interaction(method),
                     fill=method, color=method))+
    scale_x_discrete(limits=x$method_names)+
    labs(x="", y="") +
    theme(
        axis.ticks.x = element_blank(),
        legend.position = "none")
  
  if(!is.null(ORDER)){
    p <- p + ggFILL + ggBORDER
  }
  p  
})
