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
    max((intercept + growth*HI)*(1 - alpha*2*HI*(1 - HI)), 0.01)
}

## The likelihood function over a set of inds
LogLik <- function(sub.data.by.group, param.for.sub.data, response, alpha.along){
  # Here data is already splitted by groups
  l.lik <- stats::dnbinom(
    sub.data.by.group[, response],
    size = abs(param.for.sub.data[names(param.for.sub.data) %in% "k"]),
    mu = MeanLoad(
      alpha = param.for.sub.data[names(param.for.sub.data) %in% "alpha"],
      intercept = param.for.sub.data[grepl("inter", names(param.for.sub.data))],
      growth = param.for.sub.data[grepl("growth", names(param.for.sub.data))],
      HI = sub.data.by.group[, alpha.along]
    ),
    log = TRUE
  )
  sum(l.lik)
}

get.param.for.sub.data <- function(param, sub.data.by.group, sorted.group.name) {
  ## by makes sure we get all levels: get the name of the parameter
  ## from the values within the by "loop"
  param.pattern <- unique(
    interaction(
      sub.data.by.group[, sorted.group.name], 
      sep=":"
      )
    )
  ## construct a regex with it
  param.regex <- paste0("^k$|alpha|^", param.pattern)
  ## select from our ugly paramter collection
  param.for.sub.data <- param[grepl(param.regex, names(param))]
  param.for.sub.data
}

# The likelihood analysis
hybrid.maxim <- function (param, data, group.name, response = response,
                          alpha.along, hessian=FALSE, control = list(fnscale=-1)){
  ## split the name into two
  sorted.group.name <- sort(group.name)
  splitted.optimResults <- by(data, data[, sorted.group.name], function(sub.data.by.group)  {
    ## we now have the data splitted by groups
    ## we are still missing the params splitted by groups
    param.for.sub.data <- get.param.for.sub.data(param, sub.data.by.group, sorted.group.name)
    print(param.for.sub.data)
    stats::optim(par = param.for.sub.data,
                 fn = LogLik, ## function to be maximized
                 control = control, ## maximise by default
                 method = "L-BFGS-B",
                 sub.data.by.group = sub.data.by.group,
                 response = response,
                 alpha.along = alpha.along,
                 hessian = hessian)
  })
  print(splitted.optimResults)
}

##Approximation of the CI by hessian matrix
# Wald test (cf "Max Lik estimation and Inference book) p46:

## Give marginal CI, = for this parameters, without considering the correlations with the others!
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

#  Likelihood Ratio Test NB think about the dDF definition in this case...
anova.hybrid <- function(m1, m2){
  ## Test if the difference between 2 likelihood is significant
  dLL = abs(m1$twologlik/2 - m2$twologlik/2)
  dDF = length(m1$opt.param) - length(m2$opt.param)
  p = 1 - stats::pchisq(2*dLL, df = dDF) # G-test
  print(list(c(dLL = dLL, dDF = dDF, p = p)))
}
