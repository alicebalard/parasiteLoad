## Get it somehow from the formula
##' Here we describe our main formula interface
##'
##' And here we put all the details
##' @title glm.hybrid
##' @param formula 
##' @param data 
##' @param alpha.along 
##' @param alpha.start 
##' @param start.mod 
##' @param start.values 
##' @return an object of class glm.hybrid
##' @author Emanuel Heitlinger
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
        max(data[, alpha.var])!=1 ||
        min(data[, alpha.var])!=0){
      stop("glm.hybrid is currently only implemented for one continuous variable scaled between 0 and 1, along which a non-linar effect (of intensity alpha) is tested")
    } else {
      original.inverse <- nb.e$transformation$inverse
      nb.t <- as.data.frame(nb.e, transform = original.inverse)
      start.param <- nb.t[, "fit"]
      names(start.param) <-
        apply(nb.t[, c(factor.var, alpha.var)], 1,
              paste, collapse=":")
      ## this works because our maximal value of alpha.along is 1
      ## get the growth relative to zero
      start.param[grepl("\\:1$", names(start.param))] <-
        start.param[grepl("\\:1$", names(start.param))] - 
        start.param[grepl("\\:0$", names(start.param))]
      ## rename
      names(start.param) <- gsub("\\:1$", "\\.growth",
                                 names(start.param))
      names(start.param) <- gsub("\\:0$", "\\.inter",
                                 names(start.param))
      
      param <- c(k=nb$theta, alpha=alpha.start, start.param)
      param[names(start.values)] <- start.values
      opt <- hybrid.maxim(param = param, data = data,
                          group.name = factor.var,
                          response = response,
                          alpha.along = alpha.along)
      out <- list(twologlik = opt$value*2,
                  start.mod = substitute(start.mod),
                  start.param = param[names(opt$par)],
                  override.start.values = start.values[names(opt$par)], 
                  opt.param = opt$par,
                  df.residual = nb$df.residual-1,
                  converged = as.logical(opt$convergence))
      class(out) <- append(class(out),"hybrid.glm")
      return(out)
    }
  }
}