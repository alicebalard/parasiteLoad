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
LogLik <- function(data, param, group.name, response, alpha.along, whichsign = 1){
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
    l.lik <- stats::dnbinom(x[, response],
                     size=abs(sub.param[names(sub.param) %in% "k"]),
                     mu=abs(MeanLoad(alpha=sub.param[names(sub.param) %in% "alpha"],
                                     intercept=sub.param[grepl("inter",
                                                               names(sub.param))],
                                     growth=sub.param[grepl("growth", names(sub.param))],
                                     HI=x[, alpha.along])),
                     log = TRUE)
    l.lik
  })
  all.l.lik <- unlist(split.L)
  if(length(all.l.lik)!=nrow(data)){
    stop("Not all likelihoods considered, group/parameter matching problem")
  } else{
    sum(all.l.lik) * whichsign
  }
}

hybrid.maxim <- function (param, data, group.name, response = response,
                          alpha.along, hessian=FALSE, control = list(fnscale=-1),
                          whichsign = 1){
  stats::optim(par = param, 
        fn = LogLik, ## function to be maximized
        control = control, ## maximise by default
        method = "L-BFGS-B",
        data = data,
        group.name = group.name,
        response = response,
        alpha.along = alpha.along, 
        whichsign =  whichsign,
        hessian = hessian)
}

##Approximation of the CI by hessian matrix
# Wald test (cf "Max Lik estimation and Inference book) p46:

ML_bounds_Wald <- function(param, data, group.name,
                           response, alpha.along){
  # use start values inferred from glm.nb:
  fit.include.hessian <- hybrid.maxim(param = param, data = data,
                                      group.name = group.name,
                                      response = response,
                                      alpha.along = alpha.along,
                                      hessian = TRUE, whichsign = -1, 
                                      control = list())
  MLE <- fit.include.hessian$par
  ObsInfo <- fit.include.hessian$hessian # observed Fisher information matrix
  Vhat <- solve(ObsInfo) # inverse of observed Fisher information matrix
  Std.errors <- sqrt(diag(Vhat))
  # obtain the MLEs, estimated std errors, and approx Wald 95% CIs
  Wald.table <- cbind(MLE,
                      Std.errors,
                      LowerBounds = MLE - stats::qnorm(0.975)*Std.errors,
                      UpperBounds = MLE + stats::qnorm(0.975)*Std.errors)
  round(Wald.table, 4)
}

#  Likelihood Ratio Test NB think abou the dDF definition in this case...
anova.hybrid <- function(m1, m2){
  ## Test if the difference between 2 likelihood is significant
  dLL = abs(m1$twologlik/2 - m2$twologlik/2)
  dDF = length(m1$opt.param) - length(m2$opt.param)
  p = 1 - stats::pchisq(2*dLL, df = dDF) # G-test
  print(list(c(dLL = dLL, dDF = dDF, p = p)))
}
