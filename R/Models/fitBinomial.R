library(bbmle)
library(optimx)

## Function defining the distribution of the prob argument of the Binomial distribution
MeanLoad <- function(L1, L2, alpha, hybridIndex){
  heterozygoty <- 2 * hybridIndex * (1 - hybridIndex)
  mean <- (L1 + (L2 - L1) * hybridIndex) * (1 - alpha * heterozygoty)
  mean <- sapply(mean, function(x) {
    return(max(x, 0.01))
  })
  return(mean)
}

printConvergence <- function(fit) {
  convergence <- fit@details$convergence
  print(ifelse(convergence == 0, "Did converge", "Did not converge"))
}

## Fit functions for each model parameters set

# no difference between subspecies, no hybrid effect
FitBasicNoAlpha <- function(data, response, hybridIndex, paramBounds, config){
  print("Fitting model basic without alpha")
  data$response <- data[[response]] # little trick
  start <-  list(L1 = paramBounds[["L1start"]])
  fit <- mle2(
    response ~ dbinom(prob = MeanLoad(L1, L1, 0, HI),
                       size = 1),
    data = data, 
    start = start,
    lower = c(L1 = paramBounds[["L1LB"]]),
    upper = c(L1 = paramBounds[["L1UB"]]),
    optimizer = config$optimizer, 
    method = config$method, 
    control = config$control)
  printConvergence(fit)
  return(fit)
}

# no difference between subspecies, flexible hybrid effect
FitBasicAlpha <- function(data, response, hybridIndex, paramBounds, config){
  print("Fitting model basic with alpha")
  data$response <- data[[response]] # little trick
  start <-  list(L1 = paramBounds[["L1start"]],
                 alpha = paramBounds[["alphaStart"]])
  fit <- mle2(
    response ~ dbinom(prob = MeanLoad(L1, L1, alpha, HI),
                      size = 1),
    data = data,
    start = start,
    lower = c(L1 = paramBounds[["L1LB"]],
              alpha = paramBounds[["alphaLB"]]),
    upper = c(L1 = paramBounds[["L1UB"]],
              alpha = paramBounds[["alphaUB"]]),
    optimizer = config$optimizer, 
    method = config$method, 
    control = config$control)
  printConvergence(fit)
  return(fit)
}

# difference between subspecies, no hybrid effect
FitAdvancedNoAlpha <- function(data, response, hybridIndex, paramBounds, config){
  print("Fitting model advanced without alpha")
  data$response <- data[[response]]
  start <-  list(L1 = paramBounds[["L1start"]],
                 L2 = paramBounds[["L2start"]])
  fit <- mle2(
    response ~ dbinom(prob = MeanLoad(L1, L2, 0, HI),
                      size = 1),
    data = data, 
    start = start,
    lower = c(L1 = paramBounds[["L1LB"]], 
              L2 = paramBounds[["L2LB"]]),
    upper = c(L1 = paramBounds[["L1UB"]],
              L2 = paramBounds[["L2UB"]]),
    optimizer = config$optimizer, 
    method = config$method, 
    control = config$control)
  printConvergence(fit)
  return(fit)
}

# difference between subspecies, flexible hybrid effect
FitAdvancedAlpha <- function(data, response, hybridIndex, paramBounds, config){
  print("Fitting model advanced with alpha")
  data$response <- data[[response]]
  start <-  list(L1 = paramBounds[["L1start"]],
                 L2 = paramBounds[["L2start"]],
                 alpha = paramBounds[["alphaStart"]])
  fit <- mle2(
    response ~ dbinom(prob = MeanLoad(L1, L2, alpha, HI),
                      size = 1),
    data = data, 
    start = start,
    lower = c(L1 = paramBounds[["L1LB"]], 
              L2 = paramBounds[["L2LB"]], 
              alpha = paramBounds[["alphaLB"]]),
    upper = c(L1 = paramBounds[["L1UB"]],
              L2 = paramBounds[["L2UB"]],
              alpha = paramBounds[["alphaUB"]]),
    optimizer = config$optimizer, 
    method = config$method, 
    control = config$control)
  printConvergence(fit)
  return(fit)
}

run <- function (data, response, hybridIndex, paramBounds, config) {
  results = list()
  methods = c(
    fitBasicNoAlpha = FitBasicNoAlpha, 
    fitBasicAlpha = FitBasicAlpha,
    fitAdvancedNoAlpha = FitAdvancedNoAlpha,
    fitAdvancedAlpha = FitAdvancedAlpha
  )
  for (methodName in names(methods)){
    method <- methods[[methodName]]
    results[[methodName]] <- method(
      data, 
      response, 
      hybridIndex, 
      paramBounds, 
      config
    )
  }
  return(results)
}
