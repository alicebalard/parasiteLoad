#' Run analysis in user friendly formula interface
#'
#' @param formula An object of class "formula"
#' @param data A data frame containing the variables in the model
#' @param model Method to be used in fitting the model
#' So far implemented for "binomial", "negbin", "student"
#' @param group A character string. Which group is considered (e.g. "Sex")
#' @param hybridIndex The hybrid index
#' @param myparamBounds A named vector giving the start, lower and upper parameters for optimization
#' @param hybridEffect Logical. TRUE if hybrid effect (deviation from additivity between both parents) is tested, FALSE if not.
#' @param config
#' @return mlHyb returns an object of class "mlHyb"
#' The functions summary and anova are used to obtain and print a summary and analysis of variance table of the results. The generic accessor functions coefficients, effects, fitted.values and residuals extract various useful features of the value returned by lm.
#' An object of class "mlHyb" is a list containing at least the following components:
#'   coefficients
#' TO BE UPDATED
#' @export

mlHyb <- function(formula, data, model, hybridIndex = "HI", myparamBounds = "default",
                  hybridEffect = TRUE,
                  config = list(optimizer = "optimx", method = c("L-BFGS-B", "bobyqa"), control = list(follow.on = TRUE))){
  # extract response from formula
  response <- all.vars(formula)[1]
  # remove NAs
  data <- data[!is.na(data[[response]]) & !is.na(data[[hybridIndex]]),]
  # Choose model
  if (model == "negbin"){
    fit1 <- FitAdvancedAlphaNegbin
    fit0 <- FitAdvancedNoAlphaNegbin
  } else if (model == "weibull"){
    fit1 <- FitAdvancedAlphaWeibull
    fit0 <- FitAdvancedNoAlphaWeibull
  } else if (model == "normal"){
    fit1 <- FitAdvancedAlphaNormal
    fit0 <- FitAdvancedNoAlphaNormal
    fit <- FitNormal
  } else if (model == "weibullshifted"){
    fit1 <- FitAdvancedAlphaWeibullShifted
    fit0 <- FitAdvancedNoAlphaWeibullShifted
  }
  # auto or manual parameter bounds
  if (myparamBounds == "default"){
    paramBounds <- getParamBounds(model, data, response)
  } else {
    paramBounds <- myparamBounds
  }
  # is the hybrid effect tested?
  if (hybridEffect == FALSE){
    fit <- fit0(data, response, hybridIndex, paramBounds, config)
  } else {
    fit <- fit1(data, response, hybridIndex, paramBounds, config)
  }
  # create an S3 object in list
  out <- list(fit = fit, call = match.call())
  class(out) <- append(class(out), "mlHyb")
  return(out)
}

print.mlHyb <- function(x){
  coefShort <- list(intercept = x$fit@coef[["L1"]], slope = x$fit@coef[["L2"]] - x$fit@coef[["L1"]])
  if (x$call$hybridEffect == TRUE){
    coefShort <- c(coefShort, bend = x$fit@coef[["alpha"]])}
  coefDetailed <- x$fit@coef
  cat("\nCall:\n",
      paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  cat("Coefficients:\n")
  print.default(format(coefShort, digits = 3L), print.gap = 2L, quote = FALSE)
  cat("\nFull coefficients detailed:\n")
  print.default(format(coefDetailed, digits = 3L), print.gap = 2L, quote = FALSE)
}

