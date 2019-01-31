#' Fit the model for Weibull distribution of the expected load with shift of response values
#'
#' @param data A data frame
#' @param response A character string. Response (e.g. "worm_count")
#' @param hybridIndex A vector of points representing the index used as x axis
#' @param paramBounds A vector of parameters (upper, lower, start) for the optimisation
#' @param config A list containing an optimizer (default: "optimx"), a method (default "bobyqa", "L-BFGS-B") and a control (default list(follow.on = TRUE))
#' @return A fit for binomial distributed data
#' @export

# no difference between subspecies, no hybrid effect
FitBasicNoAlphaWeibullShifted <-
  function(data, response, hybridIndex, paramBounds, config){
    print("Fitting model basic without alpha")
    data$response <- data[[response]] # little trick
    HI <- data[[hybridIndex]]
    start <-  list(L1 = paramBounds[["L1start"]],
                   myshape = paramBounds[["myshapeStart"]],
                   SHIFT = paramBounds[["SHIFTStart"]])
    fit <- bbmle::mle2(
      response ~ dweibull(shape = myshape,
                          scale = (MeanLoad(L1, L1, 0, HI)+ SHIFT) /
                            gamma(1 + (1 / myshape))),
      data = data,
      start = start,
      lower = c(L1 = paramBounds[["L1LB"]],
                myshape = paramBounds[["myshapeLB"]],
                SHIFT = paramBounds[["SHIFTLB"]]),
      upper = c(L1 = paramBounds[["L1UB"]],
                myshape = paramBounds[["myshapeUB"]],
                SHIFT = paramBounds[["SHIFTUB"]]),
      optimizer = config$optimizer,
      method = config$method,
      control = config$control)
    printConvergence(fit)
    return(fit)
  }

# no difference between subspecies, flexible hybrid effect
FitBasicAlphaWeibullShifted <-
  function(data, response, hybridIndex, paramBounds, config){
    print("Fitting model basic with alpha")
    data$response <- data[[response]] # little trick
    HI <- data[[hybridIndex]]
    start <-  list(L1 = paramBounds[["L1start"]],
                   alpha = paramBounds[["alphaStart"]],
                   myshape = paramBounds[["myshapeStart"]],
                   SHIFT = paramBounds[["SHIFTStart"]])
    fit <- bbmle::mle2(
      response ~ dweibull(shape = myshape,
                          scale = (MeanLoad(L1, L1, alpha, HI)+ SHIFT) /
                            gamma(1 + (1 / myshape))),
      data = data,
      start = start,
      lower = c(L1 = paramBounds[["L1LB"]],
                alpha = paramBounds[["alphaLB"]],
                myshape = paramBounds[["myshapeLB"]],
                SHIFT = paramBounds[["SHIFTLB"]]),
      upper = c(L1 = paramBounds[["L1UB"]],
                alpha = paramBounds[["alphaUB"]],
                myshape = paramBounds[["myshapeUB"]],
                SHIFT = paramBounds[["SHIFTUB"]]),
      optimizer = config$optimizer,
      method = config$method,
      control = config$control)
    printConvergence(fit)
    return(fit)
  }

# difference between subspecies, flexible hybrid effect
FitAdvancedNoAlphaWeibullShifted <-
  function(data, response, hybridIndex, paramBounds, config){
    print("Fitting model advanced without alpha")
    data$response <- data[[response]]
    HI <- data[[hybridIndex]]
    start <-  list(L1 = paramBounds[["L1start"]],
                   L2 = paramBounds[["L2start"]],
                   myshape = paramBounds[["myshapeStart"]],
                   SHIFT = paramBounds[["SHIFTStart"]])
    fit <- bbmle::mle2(
      response ~ dweibull(shape = myshape,
                          scale = (MeanLoad(L1, L2, 0, HI)+ SHIFT) /
                            gamma(1 + (1 / myshape))),
      data = data,
      start = start,
      lower = c(L1 = paramBounds[["L1LB"]],
                L2 = paramBounds[["L2LB"]],
                myshape = paramBounds[["myshapeLB"]],
                SHIFT = paramBounds[["SHIFTLB"]]),
      upper = c(L1 = paramBounds[["L1UB"]],
                L2 = paramBounds[["L2UB"]],
                myshape = paramBounds[["myshapeUB"]],
                SHIFT = paramBounds[["SHIFTUB"]]),
      optimizer = config$optimizer,
      method = config$method,
      control = config$control)
    printConvergence(fit)
    return(fit)
  }

# difference between subspecies, flexible hybrid effect
FitAdvancedAlphaWeibullShifted <-
  function(data, response, hybridIndex, paramBounds, config){
    print("Fitting model advanced with alpha")
    data$response <- data[[response]]
    HI <- data[[hybridIndex]]
    start <-  list(L1 = paramBounds[["L1start"]],
                   L2 = paramBounds[["L2start"]],
                   alpha = paramBounds[["alphaStart"]],
                   myshape = paramBounds[["myshapeStart"]],
                   SHIFT = paramBounds[["SHIFTStart"]])
    fit <- bbmle::mle2(
      response ~ dweibull(shape = myshape,
                          scale = (MeanLoad(L1, L2, alpha, HI) + SHIFT) /
                            gamma(1 + (1 / myshape))),
      data = data,
      start = start,
      lower = c(L1 = paramBounds[["L1LB"]],
                L2 = paramBounds[["L2LB"]],
                alpha = paramBounds[["alphaLB"]],
                myshape = paramBounds[["myshapeLB"]],
                SHIFT = paramBounds[["SHIFTLB"]]),
      upper = c(L1 = paramBounds[["L1UB"]],
                L2 = paramBounds[["L2UB"]],
                alpha = paramBounds[["alphaUB"]],
                myshape = paramBounds[["myshapeUB"]],
                SHIFT = paramBounds[["SHIFTUB"]]),
      optimizer = config$optimizer,
      method = config$method,
      control = config$control)
    printConvergence(fit)
    return(fit)
  }
