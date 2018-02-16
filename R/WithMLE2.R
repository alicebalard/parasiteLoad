# install.packages("bbmle")
library(bbmle)
library(ggplot2)

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

##### Input end #####
myFun <- function(data, hybridIndex, response, 
                  L1start = 10, L1LB = 0, L1UB = 700, 
                  L2start = 10, L2LB = 0, L2UB = 700, 
                  alphaStart = 0, alphaLB = -5, alphaUB = 5,
                  A1start = 10, A1LB = 0, A1UB = 1000, 
                  A2start = 10, A2LB = 0, A2UB = 1000, 
                  Zstart = 0, ZLB = -5, ZUB = 5){
  data$response <- data[[response]] # little trick
  ################## H1 ##################
  startH1 <- list(L1 = L1start, L2 = L2start, alpha = alphaStart,
                  A1 = A1start, A2 = A2start, Z = Zstart)
  myFitAlphaH1 <- mle2(
    response ~ dnbinom(mu = MeanLoad(L1, L2, alpha, HI),
                       size = 1/abs(Aggregation(A1, A2, Z, HI))),
    data = data, method = "L-BFGS-B",
    start = startH1
  )
  myFitNoAlphaH1 <- mle2(
    response ~ dnbinom(mu = MeanLoad(L1, L2, alpha, HI),
                       size = 1/abs(Aggregation(A1, A2, Z, HI))),
    data = data, method = "L-BFGS-B",
    fixed = list(alpha = 0), start = startH1
  )
  ################## H3 ##################
  startH3 <- list(L1 = L1start, L2 = L2start, alpha = alphaStart,
                  A1 = A1start, A2 = A2start, Z = Zstart)
  myFitAlphaH3 <- mle2(
    response ~ dnbinom(mu = MeanLoad(L1, L2, alpha, HI),
                       size = 1/abs(Aggregation(A1, A2, Z, HI))),
    data = data, method = "L-BFGS-B",
    parameters = list(L1~Sex, L2~Sex, alpha~Sex, A1~Sex, A2~Sex, Z~Sex),
    start = startH3
  )
  myFitNoAlphaH3 <- mle2(
    response ~ dnbinom(mu = MeanLoad(L1, L2, alpha, HI),
                       size = 1/abs(Aggregation(A1, A2, Z, HI))),
    data = data, method = "L-BFGS-B",
    parameters = list(L1~Sex, L2~Sex, A1~Sex, A2~Sex, Z~Sex),
    fixed = list(alpha = 0), start = startH3
  )
  ################## Output ##################
  return(list(
    H1 = list(fitNoAlpha = myFitNoAlphaH1,
              fitAlpha = myFitAlphaH1),
    H3 = list(fitNoAlpha = myFitNoAlphaH3,
              fitAlpha = myFitAlphaH3)))
}

# For a given model, downstream analyses and plot
finalFun <- function(model){
  ## H1
  ## Difference between with and without alpha?
  pValueAlpha <- anova(model$H1$fitNoAlpha, model$H1$fitAlpha)[10]
  LLdropIfNoAlpha <- as.numeric(logLik(model$H1$fitNoAlpha) - logLik(model$H1$fitAlpha))
  print(paste("For H1, the p-value of the anova test when we remove alpha is", round(pValueAlpha, 5), 
              "It corresponds to a drop of likelihood of", round(LLdropIfNoAlpha,2)))
  ## H3
  ## Difference between with and without alpha?
  pValueAlpha <- anova(model$H3$fitNoAlpha, model$H3$fitAlpha)[10]
  LLdropIfNoAlpha <- as.numeric(logLik(model$H3$fitNoAlpha) - logLik(model$H3$fitAlpha))
  print(paste("For H3, the p-value of the anova test when we remove alpha is", round(pValueAlpha, 5), 
        "It corresponds to a drop of likelihood of", round(LLdropIfNoAlpha,2)))
  ## Difference H3-H1
  p <- anova(model$H3$fitAlpha, ModelPinworm$H1$fitAlpha)[10]
  print(paste("The anova between H0 and H3 has a p-value of", round(p,3)))
  if(p >= 0.05){print(paste("Therefore we keep the model without sex"))
  } else {
    print(paste("Therefore we consider sex as a significant variable"))
  }
}
