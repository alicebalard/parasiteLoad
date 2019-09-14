#' Run analysis in user friendly formula interface
#'
#' @param formula An object of class "formula"
#' @param data A data frame containing the variables in the model
#' @param model Method to be used in fitting the model
#' So far implemented for "binomial", "negbin", "student"
#' @param group A character string. Which group is considered (e.g. "Sex")
#' @param hybridIndex The hybrid index
#' @param myparamBounds A named vector giving the start, lower and upper parameters for optimization
#' @param config
#' @return mlHyb returns an object of class "mlHyb"
#' The functions summary and anova are used to obtain and print a summary and analysis of variance table of the results. The generic accessor functions coefficients, effects, fitted.values and residuals extract various useful features of the value returned by lm.
#' An object of class "mlHyb" is a list containing at least the following components:
#'   coefficients
#' TO BE UPDATED
#' @export

mlHyb <- function(formula, data, model, hybridIndex = "HI", myparamBounds = "default",
                  config = list(optimizer = "optimx", method = c("L-BFGS-B", "bobyqa"), control = list(follow.on = TRUE))){
  # extract response from formula
  response <- all.vars(formula)[1]
  # remove NAs
  data <- data[!is.na(data[[response]]) & !is.na(data[[hybridIndex]]),]
  # auto or manual parameter bounds
  if (myparamBounds == "default"){
    paramBounds <- getParamBounds(model, data, response)
  } else {
    paramBounds <- myparamBounds
  }
  # if...alpha, non alpha TODO
  fit <- FitAdvancedAlphaNegbin(data, response, hybridIndex, paramBounds, config)
  # create an S3 object in list
  out <- list(fit = fit, call = match.call())
  class(out) <- append(class(out), "mlHyb")
  return(out)
}

print.mlHyb <- function(x){
  coef <- list(intercept = x$fit@coef[["L1"]], slope = x$fit@coef[["L2"]] - x$fit@coef[["L1"]],
               bend = x$fit@coef[["alpha"]])
  cat("\nCall:\n",
      paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  cat("Coefficients:\n")
  print.default(format(coef, digits = 3L), print.gap = 2L, quote = FALSE)
}





