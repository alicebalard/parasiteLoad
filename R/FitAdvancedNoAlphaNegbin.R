#' Fit the model "difference between subspecies, no hybrid effect" for negative binomial distribution of the expected load
#'
#' @param data A data frame
#' @param response A character string. Response (e.g. "worm_count")
#' @param hybridIndex A vector of points representing the index used as x axis
#' @param paramBounds A vector of parameters (upper, lower, start) for the optimisation
#' @param config A list containing an optimizer (default: "optimx"), a method (default "bobyqa", "L-BFGS-B") and a control (default list(follow.on = TRUE))
#' @return A fit for negative binomial distributed data for difference between subspecies, no hybrid effect
#' @export

FitAdvancedNoAlphaNegbin <- function(data, response, hybridIndex, paramBounds, config){
  print("Fitting model advanced without alpha")
  data$response <- data[[response]]
  start <-  list(L1 = paramBounds[["L1start"]],
                 L2 = paramBounds[["L2start"]],
                 A1 = paramBounds[["A1start"]],
                 A2 = paramBounds[["A2start"]],
                 Z = paramBounds[["Zstart"]])
  fit <- mle2(
    response ~ dnbinom(mu = MeanLoad(L1, L2, 0, HI),
                       size = SizeNegBin(A1, A2, Z, HI)),
    data = data,
    start = start,
    lower = c(L1 = paramBounds[["L1LB"]],
              L2 = paramBounds[["L2LB"]],
              A1 = paramBounds[["A1LB"]],
              A2 = paramBounds[["A2LB"]],
              Z = paramBounds[["ZLB"]]),
    upper = c(L1 = paramBounds[["L1UB"]],
              L2 = paramBounds[["L2UB"]],
              A1 = paramBounds[["A1UB"]],
              A2 = paramBounds[["A2UB"]],
              Z = paramBounds[["ZUB"]]),
    optimizer = config$optimizer,
    method = config$method,
    control = config$control)
  printConvergence(fit)
  return(fit)
}
