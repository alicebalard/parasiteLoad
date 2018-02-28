## source the functions defining meanload and aggregation for the negative binomial
source("ModelParasiteLoad_negBin.R")

library(bbmle)
library(ggplot2)
library(optimx)

########## Function to fit all hypotheses
fit <- function(data, response, hybridIndex, paramBounds){
  data$response <- data[[response]] # little trick
  optimizer = "optimx"
  method = c("bobyqa", "L-BFGS-B")
  control = list(follow.on = TRUE)
  ## H0 and H2: considering just one L, just one A (no differences between subspecies)
  MaxLikelihoodFunBasic <- function(start, data, fixed = NULL, parameters = NULL){
    if (!is.null(fixed)) {
      fit <- mle2(
        response ~ dnbinom(mu = MeanLoad(L1, L1, alpha, HI),
                           size = SizeNegBin(A1, A1, Z, HI)),
        data = data, start = start, parameters = parameters,
        lower = c(L1 = paramBounds[["L1LB"]], 
                  A1 = paramBounds[["A1LB"]],
                  Z = paramBounds[["ZLB"]]),
        upper = c(L1 = paramBounds[["L1UB"]],
                  A1 = paramBounds[["A1UB"]],
                  Z = paramBounds[["ZUB"]]),
        fixed = fixed, optimizer = optimizer, method = method, control = control)
    } else {
      fit <- mle2(
        response ~ dnbinom(mu = MeanLoad(L1, L1, alpha, HI),
                           size = SizeNegBin(A1, A1, Z, HI)),
        data = data, start = start, parameters = parameters,
        lower = c(L1 = paramBounds[["L1LB"]],
                  A1 = paramBounds[["A1LB"]],
                  alpha = paramBounds[["alphaLB"]],
                  Z = paramBounds[["ZLB"]]),
        upper = c(L1 = paramBounds[["L1UB"]],
                  A1 = paramBounds[["A1UB"]],
                  alpha = paramBounds[["alphaUB"]],
                  Z = paramBounds[["ZUB"]]),
        optimizer = optimizer, method = method, control = control) 
    }
    convergence <- fit@details$convergence
    print(ifelse(convergence == 0, "Did converge", "Did not converge"))
    return(fit)
  }
  print("alpha basic")
  start <-  list(L1 = paramBounds[["L1start"]],
                 alpha = paramBounds[["alphaStart"]],
                 A1 = paramBounds[["A1start"]],
                 Z = paramBounds[["Zstart"]])
  myFitAlphaBasic <- MaxLikelihoodFunBasic(start, data)
  print("fixed alpha basic")
  myFitNoAlphaBasic <- MaxLikelihoodFunBasic(start, data, fixed = list(alpha = 0))
  ## H1 and H3: considering 2 L, 2 A ( differences between subspecies)
  MaxLikelihoodFunAdvanced <- function(start, data, fixed = NULL, parameters = NULL){
      if (!is.null(fixed)) {
      fit <- mle2(
        response ~ dnbinom(mu = MeanLoad(L1, L2, alpha, HI),
                           size = SizeNegBin(A1, A2, Z, HI)),
        data = data, start = start, parameters = parameters,
        lower = c(L1 = paramBounds[["L1LB"]], L2 = paramBounds[["L2LB"]],
                  A1 = paramBounds[["A1LB"]], A2 = paramBounds[["A2LB"]],
                  Z = paramBounds[["ZLB"]]),
        upper = c(L1 = paramBounds[["L1UB"]], L2 = paramBounds[["L2UB"]],
                  A1 = paramBounds[["A1UB"]], A2 = paramBounds[["A2UB"]],
                  Z = paramBounds[["ZUB"]]),
        fixed = fixed,
        optimizer = optimizer, method = method, control = control)
      } else {
        fit <- mle2(
          response ~ dnbinom(mu = MeanLoad(L1, L2, alpha, HI),
                             size = SizeNegBin(A1, A2, Z, HI)),
          data = data, start = start, parameters = parameters,
          lower = c(L1 = paramBounds[["L1LB"]], L2 = paramBounds[["L2LB"]],
                    A1 = paramBounds[["A1LB"]], A2 = paramBounds[["A2LB"]],
                    alpha = paramBounds[["alphaLB"]],
                    Z = paramBounds[["ZLB"]]),
          upper = c(L1 = paramBounds[["L1UB"]], L2 = paramBounds[["L2UB"]],
                    A1 = paramBounds[["A1UB"]], A2 = paramBounds[["A2UB"]],
                    alpha = paramBounds[["alphaUB"]],
                    Z = paramBounds[["ZUB"]]),
          optimizer = optimizer, method = method, control = control) 
      }
    convergence <- fit@details$convergence
    print(ifelse(convergence == 0, "Did converge", "Did not converge"))
    return(fit)
  }
  print("alpha advanced")
  start <-  list(L1 = paramBounds[["L1start"]], L2 = paramBounds[["L2start"]],
                 alpha = paramBounds[["alphaStart"]],
                 A1 = paramBounds[["A1start"]], A2 = paramBounds[["A2start"]],
                 Z = paramBounds[["Zstart"]])
  myFitAlphaAdvanced <- MaxLikelihoodFunAdvanced(start, data)
  print("fixed alpha advanced")
  myFitNoAlphaAdvanced <- MaxLikelihoodFunAdvanced(start, data, fixed = list(alpha = 0))
  ################## Output ##################
  return(list(fitNoAlphaBasic = myFitNoAlphaBasic,
              fitAlphaBasic = myFitAlphaBasic,
              fitNoAlphaAdvanced = myFitNoAlphaAdvanced,
              fitAlphaAdvanced = myFitAlphaAdvanced))
}

## Test the significance of alpha for each hypothesis
isAlphaSignif <- function(modelAlpha, modelNoAlpha){
  withOrWithoutalpha <- anova(modelAlpha, modelNoAlpha)
  test1 = withOrWithoutalpha[colnames(withOrWithoutalpha) == "Pr(>Chisq)"][2] < 0.05
  test2 = logLik(modelAlpha) > logLik(modelNoAlpha)
  print(paste0("Is the anova between model with or without alpha significant? ", test1))
  print(paste0(" Does the model with alpha have a lower likelihood? ", test2))
}
          
## Compare the hypotheses between each other : G-test
## Test if the difference between 2 likelihood is significant
Gtest <- function(model0, model1){
  LL0 <- Reduce(x = c(model0), f = function(accum, model){
    accum + logLik(model)}, init = 0)
  LL1 <- Reduce(x = c(model1), f = function(accum, model){
    accum + logLik(model)}, init = 0)
  dLL <- LL1 - LL0
  N0 <- Reduce(x = c(model0), f = function(accum, model){
    accum + length(coef(model))}, init = 0)
  N1 <- Reduce(x = c(model1), f = function(accum, model){
    accum + length(coef(model))}, init = 0)
  dDF <- N1 - N0
  print(paste0("Likelihood difference : ", round(dLL, 2)))
  print(paste0("Difference of degrees of freedom : ", dDF))
  pvalue <- 1 - pchisq(2*dLL, df=dDF)
  if (pvalue < 0.01){
    print("p-value < 0.01")
  } else {
    print(paste0("p-value : ", round(pvalue, 3)))
  }
}