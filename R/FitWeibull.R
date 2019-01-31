#' Fit the model for Weibull distribution of the expected load
#'
#' @param data A data frame
#' @param response A character string. Response (e.g. "worm_count")
#' @param hybridIndex A vector of points representing the index used as x axis
#' @param paramBounds A vector of parameters (upper, lower, start) for the optimisation
#' @param config A list containing an optimizer (default: "optimx"), a method (default "bobyqa", "L-BFGS-B") and a control (default list(follow.on = TRUE))
#' @return A fit for binomial distributed data
#' @export

# no difference between subspecies, no hybrid effect
FitBasicNoAlphaWeibull <- function(data, response, hybridIndex, paramBounds, config){
  print("Fitting model basic without alpha")
  data$response <- data[[response]] # little trick
  start <-  list(L1 = paramBounds[["L1start"]],
                 myshape = paramBounds[["myshapeStart"]])
  HI <- hybridIndex
  fit <- bbmle::mle2(
    response ~ dweibull(shape = myshape,
                        scale = MeanLoad(L1, L1, 0, HI) /
                          gamma(1 + (1 / myshape))),
    data = data,
    start = start,
    lower = c(L1 = paramBounds[["L1LB"]],
              myshape = paramBounds[["myshapeLB"]]),
    upper = c(L1 = paramBounds[["L1UB"]],
              myshape = paramBounds[["myshapeUB"]]),
    optimizer = config$optimizer,
    method = config$method,
    control = config$control)
  printConvergence(fit)
  return(fit)
}

# no difference between subspecies, flexible hybrid effect
FitBasicAlphaWeibull <- function(data, response, hybridIndex, paramBounds, config){
  print("Fitting model basic with alpha")
  data$response <- data[[response]] # little trick
  HI <- hybridIndex
  start <-  list(L1 = paramBounds[["L1start"]],
                 alpha = paramBounds[["alphaStart"]],
                 myshape = paramBounds[["myshapeStart"]])
  fit <- bbmle::mle2(
    response ~ dweibull(shape = myshape,
                        scale = MeanLoad(L1, L1, alpha, HI) /
                          gamma(1 + (1 / myshape))),
    data = data,
    start = start,
    lower = c(L1 = paramBounds[["L1LB"]],
              alpha = paramBounds[["alphaLB"]],
              myshape = paramBounds[["myshapeLB"]]),
    upper = c(L1 = paramBounds[["L1UB"]],
              alpha = paramBounds[["alphaUB"]],
              myshape = paramBounds[["myshapeUB"]]),
    optimizer = config$optimizer,
    method = config$method,
    control = config$control)
  printConvergence(fit)
  return(fit)
}

# difference between subspecies, flexible hybrid effect
FitAdvancedNoAlphaWeibull <- function(data, response, hybridIndex, paramBounds, config){
  print("Fitting model advanced without alpha")
  data$response <- data[[response]]
  HI <- hybridIndex
  start <-  list(L1 = paramBounds[["L1start"]],
                 L2 = paramBounds[["L2start"]],
                 myshape = paramBounds[["myshapeStart"]])
  fit <- bbmle::mle2(
    response ~ dweibull(shape = myshape,
                        scale = MeanLoad(L1, L2, 0, HI) /
                          gamma(1 + (1 / myshape))),
    data = data,
    start = start,
    lower = c(L1 = paramBounds[["L1LB"]],
              L2 = paramBounds[["L2LB"]],
              myshape = paramBounds[["myshapeLB"]]),
    upper = c(L1 = paramBounds[["L1UB"]],
              L2 = paramBounds[["L2UB"]],
              myshape = paramBounds[["myshapeUB"]]),
    optimizer = config$optimizer,
    method = config$method,
    control = config$control)
  printConvergence(fit)
  return(fit)
}

# difference between subspecies, flexible hybrid effect
FitAdvancedAlphaWeibull <- function(data, response, hybridIndex, paramBounds, config){
  print("Fitting model advanced with alpha")
  data$response <- data[[response]]
  HI <- hybridIndex
  start <-  list(L1 = paramBounds[["L1start"]],
                 L2 = paramBounds[["L2start"]],
                 alpha = paramBounds[["alphaStart"]],
                 myshape = paramBounds[["myshapeStart"]])
  fit <- bbmle::mle2(
    response ~ dweibull(shape = myshape,
                        scale = MeanLoad(L1, L2, alpha, HI) /
                          gamma(1 + (1 / myshape))),
    data = data,
    start = start,
    lower = c(L1 = paramBounds[["L1LB"]],
              L2 = paramBounds[["L2LB"]],
              alpha = paramBounds[["alphaLB"]],
              myshape = paramBounds[["myshapeLB"]]),
    upper = c(L1 = paramBounds[["L1UB"]],
              L2 = paramBounds[["L2UB"]],
              alpha = paramBounds[["alphaUB"]],
              myshape = paramBounds[["myshapeUB"]]),
    optimizer = config$optimizer,
    method = config$method,
    control = config$control)
  printConvergence(fit)
  return(fit)
}
