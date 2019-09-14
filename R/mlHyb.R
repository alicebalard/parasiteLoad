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
#' The function print creates a user friendly output
#' The functions lrt.mlHyb is used to obtain a likelihood ratio test (G-test) between models with and without bend (hybridization effect)
#' An object of class "mlHyb" is a list containing at least the following components:
#'   fit : the ML fit considering hybridization effect
#'   fitNull : the ML fit considering perfect additivity between both sides (null hypothesis on hybridization effect)
#'   call : the call
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
  # test with or without hybrid effect
  fitNull <- fit0(data, response, hybridIndex, paramBounds, config)
  fitHE <- fit1(data, response, hybridIndex, paramBounds, config)
  # create an S3 object in list
  out <- list(fit = fitHE, fitNull = fitNull, call = match.call())
  class(out) <- append(class(out), "mlHyb")
  return(out)
}

print.mlHyb <- function(x){
  if (x$call$hybridEffect == TRUE){
    coefShort <- list(intercept = x$fit@coef[["L1"]], slope = x$fit@coef[["L2"]] - x$fit@coef[["L1"]],
                      bend = x$fit@coef[["alpha"]])
    coefDetailed <- x$fit@coef
  } else {
    coefShort <- list(intercept = x$fitNull@coef[["L1"]], slope = x$fitNull@coef[["L2"]] - x$fitNull@coef[["L1"]])
    coefDetailed <- x$fitNull@coef
  }
  cat("\nCall:\n",
      paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  cat("Coefficients:\n")
  print.default(format(coefShort, digits = 3L), print.gap = 2L, quote = FALSE)
  cat("\nFull coefficients detailed:\n")
  print.default(format(coefDetailed, digits = 3L), print.gap = 2L, quote = FALSE)
}

# make it look like summary(lm(Aspiculuris_Syphacia~Sex, data = toyData))
lrt.mlHyb <- function(x){
  cat("Likelihood ratio test (G-test) for significance of bend (hybridization effect):\n\n")
  Gtest(x$fit, x$fitNull)
}
