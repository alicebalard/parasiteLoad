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





