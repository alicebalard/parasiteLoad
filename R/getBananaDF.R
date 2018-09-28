#' Get bananaDF for a given model
#'
#' @param mod A model to be ploted
#' @param hybridIndex A vector of points representing the index used as x axis
#' @return A data frame to be ploted
#' @export

getBananaDF <- function(mod, hybridIndex){
  ## profile investigates behavior of objective function near the MLE
  myProf <- profile(mod)
  ## Marginal confidence interval
  myConfInt <- confint(myProf)
  ## Get marginal confidence interval for a given parameter
  getInf <- function(paramname){
    myConfInt[rownames(myConfInt) == paramname][1]}
  getSup <- function(paramname){
    myConfInt[rownames(myConfInt) == paramname][2]}
  # Optimisation of the expectedResponse (inf and sup) for all parameters varying in their confidence interval
  expectedResponse <- function(v, hybridIndex){
    L1 = v[1]; L2 = v[2]; alpha = v[3]
    heterozygoty <- 2 * hybridIndex * (1 - hybridIndex)
    expectedResponse <- (L1 + (L2 - L1) * hybridIndex) * (1 - alpha * heterozygoty)
    return(expectedResponse)
  }
  # create dataframe to fill up
  bananaDF = data.frame(HI = hybridIndex)
  # if no L2 calculated, set L1
  if("L2" %in% rownames(myConfInt) == FALSE){
    myConfInt <- rbind(myConfInt, myConfInt[rownames(myConfInt) %in% "L1"])
    rownames(myConfInt)[rownames(myConfInt) %in% ""] <- "L2"
    # calculate response expected by the model
    bananaDF$fit <- expectedResponse(c(coef(mod)[names(coef((mod))) %in% "L1"],
                                       coef(mod)[names(coef((mod))) %in% "L1"],
                                       coef(mod)[names(coef((mod))) %in% "alpha"]),
                                     hybridIndex)
  } else {
    # calculate response expected by the model
    bananaDF$fit <- expectedResponse(c(coef(mod)[names(coef((mod))) %in% "L1"],
                                       coef(mod)[names(coef((mod))) %in% "L2"],
                                       coef(mod)[names(coef((mod))) %in% "alpha"]),
                                     hybridIndex)

  }
  # Run over HI values and optimise max and min for each
  bananaDF2 = data.frame(HI = numeric(), min = numeric(), max = numeric())
  for(i in hybridIndex){
    maxLoad <- optim(par = c(L1 = getSup("L1") - getInf("L1"),
                             L2 = getSup("L2") - getInf("L2"),
                             alpha = getSup("alpha") - getInf("alpha")),
                     fn = expectedResponse,
                     lower = c(L1 = getInf("L1"), L2 = getInf("L2"), alpha = getInf("alpha")),
                     upper = c(L1 = getSup("L1"), L2 = getSup("L2"), alpha = getSup("alpha")),
                     method = "L-BFGS-B",
                     control = list(fnscale=-1), # maximize
                     hybridIndex = i)
    minLoad <- optim(par = c(L1 = getSup("L1") - getInf("L1"),
                             L2 = getSup("L2") - getInf("L2"),
                             alpha = getSup("alpha") - getInf("alpha")),
                     fn = expectedResponse,
                     lower = c(L1 = getInf("L1"), L2 = getInf("L2"), alpha = getInf("alpha")),
                     upper = c(L1 = getSup("L1"), L2 = getSup("L2"), alpha = getSup("alpha")),
                     method = "L-BFGS-B",
                     hybridIndex = i)
    bananaDF2 = rbind(bananaDF2, data.frame(HI = i, min = minLoad$value, max = maxLoad$value))
  }
  bananaDF = merge(bananaDF, bananaDF2)
  return(bananaDF)
}
