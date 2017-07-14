## Aim for such an interface... in input functions in a seperately
## tested script before this!

## INPUT
## alice.ml.nb(bristles ~ latitude + sex + strain, data = drosophila)

## OUTPUT
## alpha                          -   estimate, lower, upper
## k                              -   estimate, lower, upper
## factor1..N * intercept         -   estimate, lower, upper
## factor1..N * growth            -   estimate, lower, upper
##                                    ML estimate


#########################
## MeanLoad model: How does the mean load vary depending on the parameters?
MeanLoad <- function(intercept, growth, alpha, HI){
    (intercept + growth*HI)*(1 - alpha*2*HI*(1 - HI))
}

## By row
LogLik <- function(data, param, group.name, response){
  ## split the name into two
  gname <- sort(group.name)
  split.L<- by(data, data[, gname], function(x)  {
    ## by makes sure we get all levels: get the name of the paramter
    ## from the values within the by "loop"
    param.pattern <- unique(interaction(x[, gname], sep=":"))
    ## construct a regex with it
    par.regex <- paste0("^k$|alpha|^", param.pattern)
    ## select from our ugly paramter collection
    sub.param <- param[grepl(par.regex, names(param))]
    l.lik <- dnbinom(x[, response],
                     size=abs(sub.param[names(sub.param) %in% "k"]),
                     mu=abs(MeanLoad(alpha=sub.param[names(sub.param) %in% "alpha"],
                                     intercept=sub.param[grepl("inter",
                                                               names(sub.param))],
                                     growth=sub.param[grepl("growth", names(sub.param))],
                                     HI=x$HI)),
                     log = TRUE)
    l.lik
  })
  all.l.lik <- unlist(split.L)
  if(length(all.l.lik)!=nrow(data)){
    stop("Not all likelihoods considered, group/paramter matching problem")
  } else{
    sum(all.l.lik)
  }
}

hybrid.maxim <- function (param, data, group.name, response = response,
                          hessian=FALSE){
  optim(par = param, 
        fn = LogLik, ## function to be maximized
        control = list(fnscale=-1),
        method = "L-BFGS-B",
        data = data,
        group.name = group.name,
        response = response,
        hessian = hessian)
}


