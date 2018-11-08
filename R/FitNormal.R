#' Fit the model for normal distribution
#'
#' @param data A data frame
#' @param response A character string. Response (e.g. "worm_count")
#' @param hybridIndex A vector of points representing the index used as x axis
#' @param paramBounds A vector of parameters (upper, lower, start) for the optimisation
#' @param config A list containing an optimizer (default: "optimx"), a method (default "bobyqa", "L-BFGS-B") and a control (default list(follow.on = TRUE))
#' @return A fit for normally distributed data
#' @export

FitBasicNoAlphaNormal <- function(data, response, hybridIndex, paramBounds, config){
  print("Fitting model basic without alpha")
  data$response <- data[[response]] # little trick
  start <-  list(L1 = paramBounds[["L1start"]],
                 mysd = paramBounds[["mysdStart"]])
  fit <- bbmle::mle2(
    response ~ dnorm(mean = MeanLoad(L1, L1, 0, HI),
                     sd = mysd),
    data = data,
    start = start,
    lower = c(L1 = paramBounds[["L1LB"]],
              mysd = paramBounds[["mysdLB"]]),
    upper = c(L1 = paramBounds[["L1UB"]],
              mysd = paramBounds[["mysdUB"]]),
    optimizer = config$optimizer,
    method = config$method,
    control = config$control)
  printConvergence(fit)
  return(fit)
}

FitBasicAlphaNormal <- function(data, response, hybridIndex, paramBounds, config){
  print("Fitting model basic with alpha")
  data$response <- data[[response]] # little trick
  start <-  list(L1 = paramBounds[["L1start"]],
                 mysd = paramBounds[["mysdStart"]],
                 alpha = paramBounds[["alphaStart"]])
  fit <- bbmle::mle2(
    response ~ dnorm(mean = MeanLoad(L1, L1, alpha, HI),
                     sd = mysd),
    data = data,
    start = start,
    lower = c(L1 = paramBounds[["L1LB"]],
              mysd = paramBounds[["mysdLB"]],
              alpha = paramBounds[["alphaLB"]]),
    upper = c(L1 = paramBounds[["L1UB"]],
              mysd = paramBounds[["mysdUB"]],
              alpha = paramBounds[["alphaUB"]]),
    optimizer = config$optimizer,
    method = config$method,
    control = config$control)
  printConvergence(fit)
  return(fit)
}

FitAdvancedNoAlphaNormal <- function(data, response, hybridIndex, paramBounds, config){
  print("Fitting model advanced without alpha")
  data$response <- data[[response]]
  start <-  list(L1 = paramBounds[["L1start"]],
                 L2 = paramBounds[["L2start"]],
                 mysd = paramBounds[["mysdStart"]])
  fit <- bbmle::mle2(
    response ~ dnorm(mean = MeanLoad(L1, L2, 0, HI),
                     sd = mysd),
    data = data,
    start = start,
    lower = c(L1 = paramBounds[["L1LB"]],
              mysd = paramBounds[["mysdLB"]],
              L2 = paramBounds[["L2LB"]]),
    upper = c(L1 = paramBounds[["L1UB"]],
              mysd = paramBounds[["mysdUB"]],
              L2 = paramBounds[["L2UB"]]),
    optimizer = config$optimizer,
    method = config$method,
    control = config$control)
  printConvergence(fit)
  return(fit)
}

FitAdvancedAlphaNormal <- function(data, response, hybridIndex, paramBounds, config){
  print("Fitting model advanced with alpha")
  data$response <- data[[response]]
  start <-  list(L1 = paramBounds[["L1start"]],
                 L2 = paramBounds[["L2start"]],
                 alpha = paramBounds[["alphaStart"]],
                 mysd = paramBounds[["mysdStart"]])
  fit <- bbmle::mle2(
    response ~ dnorm(mean = MeanLoad(L1, L2, alpha, HI),
                     sd = mysd),
    data = data,
    start = start,
    lower = c(L1 = paramBounds[["L1LB"]],
              mysd = paramBounds[["mysdLB"]],
              L2 = paramBounds[["L2LB"]],
              alpha = paramBounds[["alphaLB"]]),
    upper = c(L1 = paramBounds[["L1UB"]],
              mysd = paramBounds[["mysdUB"]],
              L2 = paramBounds[["L2UB"]],
              alpha = paramBounds[["alphaUB"]]),
    optimizer = config$optimizer,
    method = config$method,
    control = config$control)
  printConvergence(fit)
  return(fit)
}
