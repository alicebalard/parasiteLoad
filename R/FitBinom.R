#' Fit the model for binomial distribution
#'
#' @param data A data frame
#' @param response A character string. Response (e.g. "worm_count")
#' @param hybridIndex A vector of points representing the index used as x axis
#' @param paramBounds A vector of parameters (upper, lower, start) for the optimisation
#' @param config A list containing an optimizer (default: "optimx"), a method (default "bobyqa", "L-BFGS-B") and a control (default list(follow.on = TRUE))
#' @return A fit for binomial distributed data
#' @export

FitBasicNoAlphaBinomial <- function(data, response, hybridIndex, paramBounds, config){
  print("Fitting model basic without alpha")
  data$response <- data[[response]] # little trick
  HI <- data[[hybridIndex]]
  start <-  list(L1 = paramBounds[["L1start"]])
  fit <- bbmle::mle2(
    response ~ dbinom(prob = MeanLoad(L1, L1, 0, HI),
                      size = 1),
    data = data,
    start = start,
    lower = c(L1 = paramBounds[["L1LB"]]),
    upper = c(L1 = paramBounds[["L1UB"]]),
    optimizer = config$optimizer,
    method = config$method,
    control = config$control)
  printConvergence(fit)
  return(fit)
}

FitBasicAlphaBinomial <- function(data, response, hybridIndex, paramBounds, config){
  print("Fitting model basic with alpha")
  data$response <- data[[response]] # little trick
  HI <- data[[hybridIndex]]
  start <-  list(L1 = paramBounds[["L1start"]],
                 alpha = paramBounds[["alphaStart"]])
  fit <- bbmle::mle2(
    response ~ dbinom(prob = MeanLoad(L1, L1, alpha, HI),
                      size = 1),
    data = data,
    start = start,
    lower = c(L1 = paramBounds[["L1LB"]],
              alpha = paramBounds[["alphaLB"]]),
    upper = c(L1 = paramBounds[["L1UB"]],
              alpha = paramBounds[["alphaUB"]]),
    optimizer = config$optimizer,
    method = config$method,
    control = config$control)
  printConvergence(fit)
  return(fit)
}

FitAdvancedNoAlphaBinomial <- function(data, response, hybridIndex, paramBounds, config){
  print("Fitting model advanced without alpha")
  data$response <- data[[response]]
  HI <- data[[hybridIndex]]
  start <-  list(L1 = paramBounds[["L1start"]],
                 L2 = paramBounds[["L2start"]])
  fit <- bbmle::mle2(
    response ~ dbinom(prob = MeanLoad(L1, L2, 0, HI),
                      size = 1),
    data = data,
    start = start,
    lower = c(L1 = paramBounds[["L1LB"]],
              L2 = paramBounds[["L2LB"]]),
    upper = c(L1 = paramBounds[["L1UB"]],
              L2 = paramBounds[["L2UB"]]),
    optimizer = config$optimizer,
    method = config$method,
    control = config$control)
  printConvergence(fit)
  return(fit)
}

FitAdvancedAlphaBinomial <- function(data, response, hybridIndex, paramBounds, config){
  print("Fitting model advanced with alpha")
  data$response <- data[[response]]
  HI <- data[[hybridIndex]]
  start <-  list(L1 = paramBounds[["L1start"]],
                 L2 = paramBounds[["L2start"]],
                 alpha = paramBounds[["alphaStart"]])
  fit <- bbmle::mle2(
    response ~ dbinom(prob = MeanLoad(L1, L2, alpha, HI),
                      size = 1),
    data = data,
    start = start,
    lower = c(L1 = paramBounds[["L1LB"]],
              L2 = paramBounds[["L2LB"]],
              alpha = paramBounds[["alphaLB"]]),
    upper = c(L1 = paramBounds[["L1UB"]],
              L2 = paramBounds[["L2UB"]],
              alpha = paramBounds[["alphaUB"]]),
    optimizer = config$optimizer,
    method = config$method,
    control = config$control)
  printConvergence(fit)
  return(fit)
}
