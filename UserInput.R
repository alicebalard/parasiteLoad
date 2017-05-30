library(MASS)
library(effects)

## Get it somehow from the formula
glm.hybrid <- function(formula, data, 
                       alpha.along, alpha.start = 0.1, 
                       start=MASS::glm.nb){
    ## create the formula in the environment of our function
    formula <- formula(substitute(formula))
    nb <- glm.nb(formula, data=data)
    nb.e  <- allEffects(nb, xlevels=2)
    if(length(nb.e) > 1){
        stop("glm.hybrid not only currently only implements models with all potential interactions defined")
    } else {
        nb.e <- nb.e[[1]]
    }
    var <- nb.e[["variables"]]
    conti.var <- sapply(var, function (x) !x[["is.factor"]])
    if (sum(conti.var) != 1 ||
        !names(conti.var)[conti.var] %in% alpha.along){
        stop("glm.hybrid is currently only implemented for one continuous variable, which along which a non-linar effect is tested")
    } else {
        var.names <- sort(names(nb.e[["variables"]]))
        original.inverse <- nb.e$transformation$inverse
        nb.t <- as.data.frame(nb.e, transform = original.inverse)
        nb.param <- nb.t[, "fit"]
        names(nb.param) <- apply(nb.t[,var.names], 1,
                                 paste, collapse=":")
        nb.param[grepl("\\:1", names(nb.param))] <-
            nb.param[grepl("\\:1", names(nb.param))] - 
            nb.param[grepl("\\:0", names(nb.param))]
        names(nb.param) <- gsub("\\:1", "\\.growth", names(nb.param))
        names(nb.param) <- gsub("\\:0", "\\.inter", names(nb.param))
        nb.param <- c(k=nb$theta, alpha=alpha.start, nb.param)
        optim(par = nb.param,
              fn = LogLik, ## function to be maximized
              control = list(fnscale=-1),
              data = simdata,
              group.name=names(conti.var[!conti.var]))
    }
}




