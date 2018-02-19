# install.packages("bbmle")
library(bbmle)
library(ggplot2)
library(optimx)

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

## Import data WATWM
Joelle_data <- read.csv("../examples/Reproduction_WATWM/EvolutionFinalData.csv")
# to check
Joelle_data[is.na(Joelle_data)] <- 0
# pinworms (A. tetraptera and S. obvelata)
# Trichuris muris (whipworm)
# Taenia taeniaeformis (tapeworm) 

##### Input end #####
myFun <- function(data, hybridIndex, response, 
                  L1start = 10, L1LB = 0, L1UB = 700, 
                  L2start = 10, L2LB = 0, L2UB = 700, 
                  alphaStart = 0, alphaLB = -5, alphaUB = 5,
                  A1start = 10, A1LB = 0, A1UB = 1000, 
                  A2start = 10, A2LB = 0, A2UB = 1000, 
                  Zstart = 0, ZLB = -5, ZUB = 5){
  data$response <- data[[response]] # little trick
  ###
  MaxLikelihoodFun <- function(start, data, fixed = NULL, parameters = NULL){
    if (!is.null(fixed)) {
      fit <- mle2(
        response ~ dnbinom(mu = MeanLoad(L1, L2, alpha, HI),
                                     size = 1/abs(Aggregation(A1, A2, Z, HI))),
        data = data, 
        optimizer = "optimx",
        method = "bobyqa",
        start = start, 
        fixed = fixed,
        parameters = parameters
      )
    }
    else {
      fit <- mle2(
        response ~ dnbinom(mu = MeanLoad(L1, L2, alpha, HI),
                           size = 1/abs(Aggregation(A1, A2, Z, HI))),
        data = data, 
        optimizer = "optimx",
        method = "bobyqa",
        start = start, 
        parameters = parameters
      ) 
    }
    convergence <- fit@details$convergence
    print(ifelse(convergence == 0, "Did converge", "Did not converge"))
    return(fit)
  }
  ################## H1 ##################
  print("H1 alpha")
  startH1 <- list(L1 = L1start, L2 = L2start, alpha = alphaStart,
                  A1 = A1start, A2 = A2start, Z = Zstart)
  myFitAlphaH1 <- MaxLikelihoodFun(startH1, data)
  print("H1 no alpha")
  myFitNoAlphaH1 <- MaxLikelihoodFun(startH1, data, fixed = list(alpha = 0))
  ################## H3 ##################
  print("H3 alpha")
  startH3 <- list(L1 = L1start, L2 = L2start, alpha = alphaStart,
                  A1 = A1start, A2 = A2start, Z = Zstart)
  myFitAlphaH3 <- MaxLikelihoodFun(
    startH3, 
    data, 
    parameters = list(L1~Sex, L2~Sex, alpha~Sex, A1~Sex, A2~Sex, Z~Sex))
  print("H3 no alpha")
  myFitNoAlphaH3 <- MaxLikelihoodFun(
    startH3,
    data,
    parameters = list(L1~Sex, L2~Sex, A1~Sex, A2~Sex, Z~Sex),
    fixed = list(alpha = 0)
  )
  ################## Output ##################
  return(list(
    H1 = list(fitNoAlpha = myFitNoAlphaH1,
              fitAlpha = myFitAlphaH1),
    H3 = list(fitNoAlpha = myFitNoAlphaH3,
              fitAlpha = myFitAlphaH3)))
}

################## Data analysis ################## 
system.time(ModelPinworm <- myFun(data = Joelle_data, 
                                  hybridIndex = HI, 
                                  response = "Aspiculuris.Syphacia"))

system.time(ModelWhipworm <- myFun(data = Joelle_data, 
                                  hybridIndex = HI, 
                                  response = "Trichuris"))

system.time(ModelTapeworm <- myFun(data = Joelle_data, 
                                  hybridIndex = HI, 
                                  response = "Taenia"))

system.time(ModelMasto <- myFun(data = Joelle_data, 
                                  hybridIndex = HI, 
                                  response = "Mastophorus"))

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

finalFun(ModelPinworm)

################## Plotting ################## 

## profile investigates behavior of objective function near the MLE
system.time(myProf <- profile(ModelPinworm$H1$fitAlpha))

## Marginal confidence interval
myConfInt <- confint(myProf)

## Marginal confidence interval for alpha
alphaCILB <- myConfInt[rownames(myConfInt) == "alpha"][1]
alphaCIUB <- myConfInt[rownames(myConfInt) == "alpha"][2]

## Draw the line for the parameters at their MLE, alpha varying
DF <- data.frame(HI = seq(0,1,0.01),
                 loadMLE = MeanLoad(L1 = coef(ModelPinworm$H1$fitAlpha)[names(coef(ModelPinworm$H1$fitAlpha)) == "L1"],
                                    L2 =  coef(ModelPinworm$H1$fitAlpha)[names(coef(ModelPinworm$H1$fitAlpha)) == "L2"],
                                    alpha =  coef(ModelPinworm$H1$fitAlpha)[names(coef(ModelPinworm$H1$fitAlpha)) == "alpha"], 
                                    hybridIndex = seq(0,1,0.01)),
                 loadMLEnoAlpha = MeanLoad(L1 = coef(ModelPinworm$H1$fitNoAlpha)[names(coef(ModelPinworm$H1$fitNoAlpha)) == "L1"],
                                           L2 =  coef(ModelPinworm$H1$fitNoAlpha)[names(coef(ModelPinworm$H1$fitNoAlpha)) == "L2"],
                                           alpha =  coef(ModelPinworm$H1$fitNoAlpha)[names(coef(ModelPinworm$H1$fitNoAlpha)) == "alpha"], 
                                           hybridIndex = seq(0,1,0.01)),
                 loadMLEAlphaLB = MeanLoad(L1 = coef(ModelPinworm$H1$fitAlpha)[names(coef(ModelPinworm$H1$fitAlpha)) == "L1"],
                                           L2 =  coef(ModelPinworm$H1$fitAlpha)[names(coef(ModelPinworm$H1$fitAlpha)) == "L2"],
                                           alpha =  alphaCILB,
                                           hybridIndex = seq(0,1,0.01)),
                 loadMLEAlphaUB = MeanLoad(L1 = coef(ModelPinworm$H1$fitAlpha)[names(coef(ModelPinworm$H1$fitAlpha)) == "L1"],
                                           L2 =  coef(ModelPinworm$H1$fitAlpha)[names(coef(ModelPinworm$H1$fitAlpha)) == "L2"],
                                           alpha =  alphaCIUB,
                                           hybridIndex = seq(0,1,0.01)))
     
ggplot() +
  geom_point(data = Joelle_data, aes(x = HI, y = log10(Aspiculuris.Syphacia + 1), color = Sex)) +
  geom_ribbon(aes(x = DF$HI, 
                  ymin = log10(DF$loadMLEAlphaUB + 1), 
                  ymax = log10(DF$loadMLEAlphaLB + 1)),
              fill = "grey", alpha = .5) +
  geom_line(aes(x = DF$HI, y = log10(DF$loadMLE + 1))) +
  geom_line(aes(x = DF$HI, y = log10(DF$loadMLEnoAlpha + 1)), linetype="dotted") +
  scale_color_manual(values = c("red", "darkblue")) +
  theme_linedraw()
