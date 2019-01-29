#' Fit the four cases for weibull distribution of the expected load
#'
#' Four hypotheses tested : no difference between subspecies, no hybrid effect;
#' no difference between subspecies, flexible hybrid effect
#' difference between subspecies, no hybrid effect
#' difference between subspecies, flexible hybrid effect
#'
#' @param data A data frame
#' @param response A character string. Response (e.g. "worm_count")
#' @param hybridIndex A vector of points representing the index used as x axis
#' @param paramBounds A vector of parameters (upper, lower, start) for the optimisation
#' @param config A list containing an optimizer (default: "optimx"), a method (default "bobyqa", "L-BFGS-B") and a control (default list(follow.on = TRUE))
#' @return A fit for binomial distributed data for all four cases
#' @export

runWeibullShifted <- function (data, response, hybridIndex, paramBounds, config) {
  results = list()
  methods = c(
    fitBasicNoAlpha = FitBasicNoAlphaWeibullShifted,
    fitBasicAlpha = FitBasicAlphaWeibullShifted,
    fitAdvancedNoAlpha = FitAdvancedNoAlphaWeibullShifted,
    fitAdvancedAlpha = FitAdvancedAlphaWeibullShifted
  )
  for (methodName in names(methods)){
    method <- methods[[methodName]]
    results[[methodName]] <- method(
      data,
      response,
      hybridIndex,
      paramBounds,
      config
    )
  }
  return(results)
}
