#' Fit the model "no difference between subspecies, no hybrid effect" for binomial distribution of the expected load
#'
#' @param data A data frame
#' @param response A character string. Response (e.g. "worm_count")
#' @param hybridIndex A vector of points representing the index used as x axis
#' @param paramBounds A vector of parameters (upper, lower, start) for the optimisation
#' @param config A list containing an optimizer (default: "optimx"), a method (default "bobyqa", "L-BFGS-B") and a control (default list(follow.on = TRUE))
#' @return A fit for binomial distributed data for no difference between subspecies, no hybrid effect
#' @export

FitBasicNoAlphaBinomial <- function(data, response, hybridIndex, paramBounds, config){
  print("Fitting model basic without alpha")
  data$response <- data[[response]] # little trick
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
