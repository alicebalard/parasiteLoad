## Get it somehow from the formula
##' glm.hybrid is used to  fit a generalised linear model, using a negative binomial distribution 
##' for which the parameter mu is defined by the function "MeanLoad".
##'
##' @title Fit a Negative Binomial Generalized Linear Model along a gradient between 0 and 1
##' @param formula an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted.
##' @param data an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. If not found in data, the variables are taken from environment(formula), typically the environment from which glm.hybrid is called.
##' @param alpha.along an object of class character defining along which gradient the response is estimated
##' @param alpha.start optional starting value for alpha. By default 0.1.
##' @param start.mod optional function used to calculate the starting values before optimisation. By default, glm.nb() from the package MASS.
##' @param start.values optional values for the parameters. By default, NA.
##' @return An object of class "glm.hybrid", a list (see \code{details})
##' @details The list returned contains:
##' \itemize{
##' \item twologlik. Twice the maximum likelihood calculated
##' \item start.mod. The start.mod that were given as argument
##' \item start.param. blabla
##' \item override.start.values. The equivalent start.values, when given as argument (otherwise, returns NA)
##' \item opt.param. The optimised parameters
##' \item opt.lower. Lower bounds calculated thanks to the hessian matrix genereated by 
##' \code{optim} function. So far, extremely approximative...
##' \item opt.upper. Upper bounds calculated thanks to the hessian matrix genereated by 
##' \code{optim} function. So far, extremely approximative...
##' \item df.residual. nb$df.residual-1 [...]
##' \item converged. ogical. \code{TRUE} if converged
##' }
##' @references Baird, S. J. E., Ribas A., Macholán M., Albrecht T., Pialek J. and Goüy de Bellocq J. (2012) 

##' \emph{Where are the wormy mice? Reexamination of hybrid parasitism in the European House Mouse Hybrid Zone}
##' @author Emanuel Heitlinger, Alice Balard
##' @export

glm.hybrid <- function(formula, data, 
                       alpha.along,
                       alpha.start = 0.1, 
                       start.mod = MASS::glm.nb, 
                       start.values = NA){
  ## create the formula in the environment of our function
  formula <- formula(substitute(formula))
  response <- all.vars(formula)[1]
  if(!class(start.mod)%in%"function"){
    stop("supply a function to estimate starting parameters, even if you supply parameters verbatim (via start.values) use this for structure")
  }
  nb <- start.mod(formula, data=data)
  nb.e  <- effects::allEffects(nb, xlevels=2)
  if(length(nb.e) > 1){
    stop("glm.hybrid not only currently only implements models with all potential interactions defined")
  } else {
    nb.e <- nb.e[[1]]
    var <- nb.e[["variables"]]
    is.factor.var <- sapply(var, function (x) x[["is.factor"]])
    factor.var <- sort(names(is.factor.var[is.factor.var]))
    alpha.var <- names(is.factor.var[!is.factor.var])
    
    if (length(alpha.var) != 1 ||
        !all(alpha.var %in% alpha.along) ||
        max(data[, alpha.var]) > 1 ||
        min(data[, alpha.var]) < 0){
      stop("glm.hybrid is currently only implemented for one continuous variable scaled between 0 and 1, along which a non-linar effect (of intensity alpha) is tested")
    } else {
      original.inverse <- nb.e$transformation$inverse
      nb.t <- as.data.frame(nb.e, transform = original.inverse)
      start.param <- nb.t[, "fit"]
 
      ## if > 1 group, then nb.t[ , factor.var] is a data.frame
      if (length(factor.var) == 1){ 
        names(start.param) <- nb.t[, factor.var] 
      } else { 
        names(start.param) <- apply(nb.t[, factor.var], 1, paste, collapse=":")
      }
      names(start.param) <- paste0(names(start.param),
                                   rep(c(".inter", ".growth"),
                                       times=length(start.param)/2))
      ## we could reduce this by simply using ever even number parameter as intercept uneven as growth
      ## same in the ML loglik function
      start.param[grepl("growth$", names(start.param))] <-
        start.param[grepl("growth$", names(start.param))] -
        start.param[grepl("inter$", names(start.param))]
      param <- c(k=nb$theta, alpha=alpha.start, start.param)
      param[names(start.values)] <- start.values
      
      opt <- hybrid.maxim(param = param, data = data,
                          group.name = factor.var,
                          response = response,
                          alpha.along = alpha.along)
      ## add proxy of 95%CI
      # bounds <- ML_bounds_Wald(param = param, data = data,
      #                         group.name = factor.var,
      #                        response = response,
      #                       alpha.along = alpha.along)
      out <- list(twologlik = opt$value*2,
                  start.mod = substitute(start.mod),
                  start.param = param[names(opt$par)],
                  override.start.values = start.values[names(opt$par)],
                  opt.param = opt$par,
                  #           opt.lower = bounds[,"LowerBounds"],
                  #           opt.upper = bounds[,"UpperBounds"],
                  df.residual = nb$df.residual-1,
                  converged = as.logical(opt$convergence))
      class(out) <- append(class(out),"hybrid.glm")
      return(out)
    }
  }
}