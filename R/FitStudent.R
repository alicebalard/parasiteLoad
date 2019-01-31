#' Fit the model for Student distribution of the expected load
#'
#' @param data A data frame
#' @param response A character string. Response (e.g. "worm_count")
#' @param hybridIndex A vector of points representing the index used as x axis
#' @param paramBounds A vector of parameters (upper, lower, start) for the optimisation
#' @param config A list containing an optimizer (default: "optimx"), a method (default "bobyqa", "L-BFGS-B") and a control (default list(follow.on = TRUE))
#' @return A fit for Student distributed data
#' @export

FitBasicNoAlphaStudent <- function(data, response, hybridIndex, paramBounds, config){
  print("Fitting model basic without alpha")
  data$response <- data[[response]] # little trick
  HI <- hybridIndex
  start <-  list(L1 = paramBounds[["L1start"]],
                 mydf = paramBounds[["mydfStart"]])
  fit <- bbmle::mle2(
    response ~ dt(ncp = MeanLoad(L1, L1, 0, HI),
                  df = mydf),
    data = data,
    start = start,
    lower = c(L1 = paramBounds[["L1LB"]],
              mydf = paramBounds[["mydfLB"]]),
    upper = c(L1 = paramBounds[["L1UB"]],
              mydf = paramBounds[["mydfUB"]]),
    optimizer = config$optimizer,
    method = config$method,
    control = config$control)
  printConvergence(fit)
  return(fit)
}

FitBasicAlphaStudent <- function(data, response, hybridIndex, paramBounds, config){
  print("Fitting model basic with alpha")
  data$response <- data[[response]] # little trick
  HI <- hybridIndex
  start <-  list(L1 = paramBounds[["L1start"]],
                 alpha = paramBounds[["alphaStart"]],
                 mydf = paramBounds[["mydfStart"]])
  fit <- bbmle::mle2(
    response ~ dt(ncp = MeanLoad(L1, L1, alpha, HI),
                  df = mydf),
    data = data,
    start = start,
    lower = c(L1 = paramBounds[["L1LB"]],
              alpha = paramBounds[["alphaLB"]],
              mydf = paramBounds[["mydfLB"]]),
    upper = c(L1 = paramBounds[["L1UB"]],
              alpha = paramBounds[["alphaUB"]],
              mydf = paramBounds[["mydfUB"]]),
    optimizer = config$optimizer,
    method = config$method,
    control = config$control)
  printConvergence(fit)
  return(fit)
}

FitAdvancedNoAlphaStudent <- function(data, response, hybridIndex, paramBounds, config){
  print("Fitting model advanced without alpha")
  data$response <- data[[response]]
  HI <- hybridIndex
  start <-  list(L1 = paramBounds[["L1start"]],
                 L2 = paramBounds[["L2start"]],
                 mydf = paramBounds[["mydfStart"]])
  fit <- bbmle::mle2(
    response ~ dt(ncp = MeanLoad(L1, L2, 0, HI),
                  df = mydf),
    data = data,
    start = start,
    lower = c(L1 = paramBounds[["L1LB"]],
              L2 = paramBounds[["L2LB"]],
              mydf = paramBounds[["mydfLB"]]),
    upper = c(L1 = paramBounds[["L1UB"]],
              L2 = paramBounds[["L2UB"]],
              mydf = paramBounds[["mydfUB"]]),
    optimizer = config$optimizer,
    method = config$method,
    control = config$control)
  printConvergence(fit)
  return(fit)
}

FitAdvancedAlphaStudent <- function(data, response, hybridIndex, paramBounds, config){
  print("Fitting model advanced with alpha")
  data$response <- data[[response]]
  HI <- hybridIndex
  start <-  list(L1 = paramBounds[["L1start"]],
                 L2 = paramBounds[["L2start"]],
                 alpha = paramBounds[["alphaStart"]],
                 mydf = paramBounds[["mydfStart"]])
  fit <- bbmle::mle2(
    response ~ dt(ncp = MeanLoad(L1, L2, alpha, HI),
                  df = mydf),
    data = data,
    start = start,
    lower = c(L1 = paramBounds[["L1LB"]],
              L2 = paramBounds[["L2LB"]],
              alpha = paramBounds[["alphaLB"]],
              mydf = paramBounds[["mydfLB"]]),
    upper = c(L1 = paramBounds[["L1UB"]],
              L2 = paramBounds[["L2UB"]],
              alpha = paramBounds[["alphaUB"]],
              mydf = paramBounds[["mydfUB"]]),
    optimizer = config$optimizer,
    method = config$method,
    control = config$control)
  printConvergence(fit)
  return(fit)
}

