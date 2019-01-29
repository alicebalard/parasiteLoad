#' Get our parameters at strat, lower, upper, for optimisation
#'
#' @param model Method to be used in fitting the model
#' So far implemented for "binomial", "negbin", "student", "normal", "weibull"
#' @param data A data frame
#' @param response A character string. Response (e.g. "worm_count")
#' @return A vector of parameters (upper, lower, start) for the optimisation
#' @export

getParamBounds <- function(model, data, response){
  if (model == "binomial"){
    paramBounds <- c(L1start = 0,
                     L1LB = 0,
                     L1UB = 1,
                     L2start = 0,
                     L2LB = 0,
                     L2UB = 1,
                     alphaStart = 0, alphaLB = -5, alphaUB = 5)
  } else if (model == "normal"){
    paramBounds <- c(L1start = mean(stats::na.omit(data[[response]])),
                     L1LB = min(stats::na.omit(data[[response]])),
                     L1UB = max(stats::na.omit(data[[response]])),
                     L2start = mean(stats::na.omit(data[[response]])),
                     L2LB = min(stats::na.omit(data[[response]])),
                     L2UB = max(stats::na.omit(data[[response]])),
                     alphaStart = 0, alphaLB = -5, alphaUB = 5,
                     mysdStart = 1, mysdLB = 0, mysdUB = 10)
    } else if (model == "student"){
    paramBounds <- c(L1start = mean(stats::na.omit(data[[response]])),
                     L1LB = min(stats::na.omit(data[[response]])),
                     L1UB = max(stats::na.omit(data[[response]])),
                     L2start = mean(stats::na.omit(data[[response]])),
                     L2LB = min(stats::na.omit(data[[response]])),
                     L2UB = max(stats::na.omit(data[[response]])),
                     alphaStart = 0, alphaLB = -5, alphaUB = 5,
                     mydfStart = 1, mydfLB = 1, mydfUB = 10)
  } else if (model == "negbin"){
    paramBounds <- c(L1start = mean(stats::na.omit(data[[response]])),
                     L1LB = 0,
                     L1UB = max(stats::na.omit(data[[response]])),
                     L2start = mean(stats::na.omit(data[[response]])),
                     L2LB = 0,
                     L2UB = max(stats::na.omit(data[[response]])),
                     alphaStart = 0, alphaLB = -5, alphaUB = 5,
                     A1start = 10, A1LB = 1e-9, A1UB = 1000,
                     A2start = 10, A2LB = 1e-9, A2UB = 1000,
                     Zstart = 0, ZLB = -20, ZUB = 20)
  } else if (model == "weibull"){
    paramBounds <- c(L1start = mean(stats::na.omit(data[[response]])),
                     L1LB = 1e-9,
                     L1UB = max(stats::na.omit(data[[response]])),
                     L2start = mean(stats::na.omit(data[[response]])),
                     L2LB = 1e-9,
                     L2UB = max(stats::na.omit(data[[response]])),
                     alphaStart = 0, alphaLB = -5, alphaUB = 5,
                     myshapeStart = 1, myshapeLB = 1e-9, myshapeUB = 5)
  } else if (model == "weibullshifted"){
    paramBounds <- c(L1start = mean(stats::na.omit(data[[response]])),
                     L1LB = 1e-9,
                     L1UB = max(stats::na.omit(data[[response]])),
                     L2start = mean(stats::na.omit(data[[response]])),
                     L2LB = 1e-9,
                     L2UB = max(stats::na.omit(data[[response]])),
                     alphaStart = 0, alphaLB = -5, alphaUB = 5,
                     myshapeStart = 1, myshapeLB = 1e-9, myshapeUB = 5,
                     SHIFTStart = 1, SHIFTLB = 1e-9, SHIFTUB = 10)
  }
  return(paramBounds)
}
