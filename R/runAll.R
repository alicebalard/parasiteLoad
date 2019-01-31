#' Run full analysis for the 4 hypotheses
#'
#' @param data A data frame
#' @param response A character string. Response (e.g. "worm_count")
#' @param model Method to be used in fitting the model
#' So far implemented for "binomial", "negbin", "student"
#' @param group A character string. Which group is considered (e.g. "Sex")
#' @param hybridIndex The hybrid index
#' @param paramBounds A named vector giving the start, lower and upper parameters for optimization
#' @return The full result of the fit
#' @export

runAll <- function (data, response, model, group, paramBounds, hybridIndex) {
  print(paste0("Fit for the response: ", response))
  defaultConfig <- list(optimizer = "optimx",
                        method = c("L-BFGS-B", "bobyqa"),
                        control = list(follow.on = TRUE))
  paramBounds <- paramBounds
  marshalledData <- marshallData(data, response, group)
  # Choose model
  if (model == "binomial"){
    run <- runBinomial
  } else if (model == "student"){
    run <- runStudent
  } else if (model == "negbin"){
    run <- runNegbin
  } else if (model == "weibull"){
    run <- runWeibull
  } else if (model == "normal"){
    run <- runNormal
  } else if (model == "weibullshifted"){
    run <- runWeibullShifted
  }
  # then run the analysis
  print("Fitting for all")
  FitAll <- run(
    data = marshalledData[["all"]],
    response = response,
    hybridIndex = hybridIndex,
    paramBounds = paramBounds,
    config = defaultConfig
  )
  print(paste0("Fitting for groupA : ", levels(data[[group]])[1]))
  FitGroupA <- run(
    data = marshalledData[["groupA"]],
    response = response,
    hybridIndex = hybridIndex,
    paramBounds = paramBounds,
    config = defaultConfig
  )
  print(paste0("Fitting for groupB : ", levels(data[[group]])[2]))
  FitGroupB <- run(
    data = marshalledData[["groupB"]],
    response = response,
    hybridIndex = hybridIndex,
    paramBounds = paramBounds,
    config = defaultConfig
  )
  return(list(FitAll = FitAll, FitGroupA = FitGroupA, FitGroupB = FitGroupB))
}
