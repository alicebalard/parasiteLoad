#' Get bananaDF for a given model
#'
#' @param mod A model to be ploted
#' @param hybridIndex A vector of points representing the index used as x axis
#' @return A data frame to be ploted
#' @export

getBananaDF_alpha <- function(mod, hybridIndex){
  ## Fitted coefficients
  fittedCoef <- bbmle::coef(mod)

  ## profile investigates behavior of objective function near the MLE
  myProf <- bbmle::profile(mod)
  ## Marginal confidence interval
  myConfInt <- bbmle::confint(myProf)

  ## Get marginal confidence interval for a alpha
  alphaInf <- myConfInt[rownames(myConfInt) == "alpha"][1]
  alphaSup <- myConfInt[rownames(myConfInt) == "alpha"][2]

  # create dataframe to fill up
  bananaDF = data.frame(HI = hybridIndex)

  # if no L2 calculated, set L1
  if("L2" %in% names(fittedCoef) == FALSE){
    fittedCoef <- c(fittedCoef, fittedCoef[names(fittedCoef) %in% "L1"])
    names(fittedCoef)[1] <- "L2"

    # calculate response expected by the model
    bananaDF$fit <- MeanLoad(L1 = fittedCoef[["L1"]],
                             L2 = fittedCoef[["L1"]],
                             alpha = fittedCoef[["alpha"]],
                             hybridIndex = hybridIndex)
  } else {
    # calculate response expected by the model
    bananaDF$fit <- MeanLoad(L1 = fittedCoef[["L1"]],
                             L2 = fittedCoef[["L2"]],
                             alpha = fittedCoef[["alpha"]],
                             hybridIndex = hybridIndex)
  }

  # Run over HI values and optimise max and min for each
  bananaDF2 = data.frame(HI = numeric(), min = numeric(), max = numeric())

  for(i in hybridIndex){
    maxLoad <- optim(par = alphaSup - alphaInf,
                     fn = MeanLoad,
                     lower = alphaInf,
                     upper = alphaSup,
                     L1 = fittedCoef[["L1"]],
                     L2 = fittedCoef[["L2"]],
                     method = "L-BFGS-B",
                     control = list(fnscale=-1), # maximize
                     hybridIndex = i)
    minLoad <- optim(par = alphaSup - alphaInf,
                     fn = MeanLoad,
                     lower = alphaInf,
                     upper = alphaSup,
                     L1 = fittedCoef[["L1"]],
                     L2 = fittedCoef[["L2"]],
                     method = "L-BFGS-B", # minimize
                     hybridIndex = i)
    bananaDF2 = rbind(bananaDF2, data.frame(HI = i, min = minLoad$value, max = maxLoad$value))
  }
  bananaDF = merge(bananaDF, bananaDF2)
  return(bananaDF)
}
