#########################
## MeanLoad model: How does the mean load vary depending on the parameters?
MeanLoad <- function(intercept, growth, alpha, HI){
  (intercept + growth*HI)*(1 - alpha*2*HI*(1 - HI))
}

IndividualLogLik_negbin <- function() {
  function(x) {
    ## by makes sure we get all levels: get the name of the parameter
    ## from the values within the by "loop"
    param.pattern <- unique(interaction(x[, gname], sep=":"))
    ## construct a regex with it
    par.regex <- paste0("^k$|alpha|^", param.pattern)
    ## select from our ugly paramter collection
    sub.param <- param[grepl(par.regex, names(param))]
    l.lik <- stats::dnbinom(x[, response],
                            size=abs(sub.param[names(sub.param) %in% "k"]),
                            mu=abs(MeanLoad(alpha=sub.param[names(sub.param) %in% "alpha"],
                                            intercept=sub.param[grepl("inter",
                                                                      names(sub.param))],
                                            growth=sub.param[grepl("growth", names(sub.param))],
                                            HI=x[, alpha.along])),
                            log = TRUE)
    l.lik
  }
}