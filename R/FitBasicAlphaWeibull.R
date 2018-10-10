#' Fit the model "no difference between subspecies, flexible hybrid effect" for Weibull distribution of the expected load
#'
#' @param data A data frame
#' @param response A character string. Response (e.g. "worm_count")
#' @param hybridIndex A vector of points representing the index used as x axis
#' @param paramBounds A vector of parameters (upper, lower, start) for the optimisation
#' @param config A list containing an optimizer (default: "optimx"), a method (default "bobyqa", "L-BFGS-B") and a control (default list(follow.on = TRUE))
#' @return A fit for Student distributed data for no difference between subspecies, flexible hybrid effect
#' @export

FitBasicAlphaWeibull <- function(data, response, hybridIndex, paramBounds, config){
  print("Fitting model basic with alpha")
  data$response <- data[[response]] # little trick
  start <-  list(L1 = paramBounds[["L1start"]],
                 alpha = paramBounds[["alphaStart"]],
                 myshape = paramBounds[["myshapeStart"]])
  fit <- bbmle::mle2(
    response ~ dweibull(shape = myshape,
                        scale = MeanLoad(L1, L1, alpha, HI)/myshape),
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
