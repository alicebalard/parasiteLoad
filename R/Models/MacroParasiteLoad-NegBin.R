## Functions defining the distribution of mu and 1/k of the Negative binomial distribution
MeanLoad <- function(L1, L2, alpha, hybridIndex){
  heterozygoty <- 2 * hybridIndex * (1 - hybridIndex)
  mean <- (L1 + (L2 - L1) * hybridIndex) * (1 - alpha * heterozygoty)
  mean <- sapply(mean, function(x) {
    return(max(x, 0.01))
  })
  return(mean)
}

Aggregation <- function(A1, A2, Z, hybridIndex){
  heterozygoty <- 2 * hybridIndex * (1 - hybridIndex)
  aggregation <- (A1 + (A2 - A1) * hybridIndex) + Z * heterozygoty 
  return(aggregation)
} 

SizeNegBin <- function(A1, A2, Z, hybridIndex){
  aggregation <- Aggregation(A1, A2, Z, hybridIndex)
  aggregation <- sapply(aggregation, function(x) {
    return(max(x, 0.01))
  })
  size <- 1/aggregation
  return(size)
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
  start <-  list(L1 = paramBounds[["L1start"]],
                 A1 = paramBounds[["A1start"]],
                 Z = paramBounds[["Zstart"]])
  fit <- mle2(
    response ~ dnbinom(mu = MeanLoad(L1, L1, 0, HI),
                       size = SizeNegBin(A1, A1, Z, HI)),
    data = data, 
    start = start,
    lower = c(L1 = paramBounds[["L1LB"]], 
              A1 = paramBounds[["A1LB"]],
              Z = paramBounds[["ZLB"]]),
    upper = c(L1 = paramBounds[["L1UB"]],
              A1 = paramBounds[["A1UB"]],
              Z = paramBounds[["ZUB"]]),
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
                 alpha = paramBounds[["alphaStart"]],
                 A1 = paramBounds[["A1start"]],
                 Z = paramBounds[["Zstart"]])
  fit <- mle2(
    response ~ dnbinom(mu = MeanLoad(L1, L1, alpha, HI),
                       size = SizeNegBin(A1, A1, Z, HI)),
    data = data,
    start = start,
    lower = c(L1 = paramBounds[["L1LB"]],
              A1 = paramBounds[["A1LB"]],
              alpha = paramBounds[["alphaLB"]],
              Z = paramBounds[["ZLB"]]),
    upper = c(L1 = paramBounds[["L1UB"]],
              A1 = paramBounds[["A1UB"]],
              alpha = paramBounds[["alphaUB"]],
              Z = paramBounds[["ZUB"]]),
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
                 L2 = paramBounds[["L2start"]],
                 A1 = paramBounds[["A1start"]],
                 A2 = paramBounds[["A2start"]],
                 Z = paramBounds[["Zstart"]])
  fit <- mle2(
    response ~ dnbinom(mu = MeanLoad(L1, L2, 0, HI),
                       size = SizeNegBin(A1, A2, Z, HI)),
    data = data, 
    start = start,
    lower = c(L1 = paramBounds[["L1LB"]], 
              L2 = paramBounds[["L2LB"]], 
              A1 = paramBounds[["A1LB"]],
              A2 = paramBounds[["A2LB"]],
              Z = paramBounds[["ZLB"]]),
    upper = c(L1 = paramBounds[["L1UB"]],
              L2 = paramBounds[["L2UB"]],
              A1 = paramBounds[["A1UB"]],
              A2 = paramBounds[["A2UB"]],
              Z = paramBounds[["ZUB"]]),
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
                 A1 = paramBounds[["A1start"]],
                 A2 = paramBounds[["A2start"]],
                 alpha = paramBounds[["alphaStart"]],
                 Z = paramBounds[["Zstart"]])
  fit <- mle2(
    response ~ dnbinom(mu = MeanLoad(L1, L2, alpha, HI),
                       size = SizeNegBin(A1, A2, Z, HI)),
    data = data, 
    start = start,
    lower = c(L1 = paramBounds[["L1LB"]], 
              L2 = paramBounds[["L2LB"]], 
              A1 = paramBounds[["A1LB"]],
              A2 = paramBounds[["A2LB"]],
              alpha = paramBounds[["alphaLB"]],
              Z = paramBounds[["ZLB"]]),
    upper = c(L1 = paramBounds[["L1UB"]],
              L2 = paramBounds[["L2UB"]],
              A1 = paramBounds[["A1UB"]],
              A2 = paramBounds[["A2UB"]],
              alpha = paramBounds[["alphaUB"]],
              Z = paramBounds[["ZUB"]]),
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
