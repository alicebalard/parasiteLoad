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

## Import data WATWM
Joelle_data <- read.csv("../examples/Reproduction_WATWM/EvolutionFinalData.csv")
# to check
Joelle_data[is.na(Joelle_data)] <- 0
# pinworms (A. tetraptera and S. obvelata)
# Trichuris muris (whipworm)
# Taenia taeniaeformis (tapeworm) 

## Simulate data for more testing
simparams <- c(k = 2,
               "male:old.alpha" = 1.4,
               "male:young.alpha" = 1.2,
               "male:baby.alpha" = 1.0,
               "female:old.alpha" = 2.0,
               "female:young.alpha" = -1.8,
               "female:baby.alpha" = 1.1,
               "male:old.inter" = 14,
               "male:young.inter" = 12,
               "male:baby.inter" = 10,
               "female:old.inter" = 20,
               "female:young.inter" = 18,
               "female:baby.inter" = 11,
               "male:old.growth" = 2,
               "male:young.growth" = 1,
               "male:baby.growth" = -4,
               "female:old.growth" = 2,
               "female:young.growth" = 0,
               "female:baby.growth" = -1)

SimulatedData <- function(param, n){
  gdata <- data.frame(
    group1 = rep(c("male", "female"), each=n/2),
    group2 = sample(c("old", "young", "baby"), n, replace=TRUE)
  )
  gdata$HI<- round(runif(n), 2)
  xloads <- by(
    gdata, 
    gdata$group1:gdata$group2, 
    function (x) {
      pattern <- paste0("^", unique(x$group1), ":", unique(x$group2))
      this.param <- param[grepl(pattern, names(param))]
      loads <- rnbinom(
        n = nrow(x), 
        size = param["k"],
        mu = MeanLoad(
          L1 = this.param[grepl("\\.inter", names(this.param))],
          L2 = this.param[grepl("\\.growth", names(this.param))],
          alpha = this.param[grepl("\\.alpha", names(this.param))],
          hybridIndex = x$HI
        )
      )
      cbind(x, loads)
    }
  )
  as.data.frame(do.call("rbind", xloads))
}

set.seed(5)
simdata <- SimulatedData(simparams, 1000)

##### Input end #####

myFun <- function(data, hybridIndex, response){
  data$response <- response
  ### H0 hypothesis
  # ## WITH hybridization effect on the mean load
  # myFitAlphaH0 <- mle2(
  #   response ~ dnbinom(mu = MeanLoad(L1, L2, alpha, HI),
  #                      size = 1/abs(Aggregation(A1, A2, Z, HI))),
  #   data = data,
  #   fixed = list(L2 = 0, A2 = 0),
  #   start = list(L1 = 10, alpha = 0, A1 = 1, Z = 0),
  #   method="L-BFGS-B",
  #   lower = c(L1 = 0, A1 = 0),
  #   upper = c(L1 = 200, A1 = 1000))
  # ## WITHOUT Hybridization effect on the mean load
  # myFitNoAlphaH0 <- mle2(
  #   response ~ dnbinom(mu = MeanLoad(L1, L2, alpha, HI),
  #                      size = 1/abs(Aggregation(A1, A2, Z, HI))),
  #   data = data,
  #   fixed = list(alpha = 0, L2 = 0, A2 = 0),
  #   start = list(L1 = 10, alpha = 0, A1 = 1, Z = 0),
  #   method="L-BFGS-B",
  #   lower = c(L1 = 0, A1 = 0),
  #   upper = c(L1 = 200, A1 = 1000))
  # ##### TO DO H1,H2
  ### H3 hypothesis
  ## WITH hybridization effect on the mean load
  myFitAlphaH3 <- mle2(
    response ~ dnbinom(mu = MeanLoad(L1, L2, alpha, HI),
                       size = 1/abs(Aggregation(A1, A2, Z, HI))),
    parameters = list(L1~Sex, L2~Sex, alpha~Sex, A1~Sex, A2~Sex, Z~Sex),
    data = data,
    start = list(L1 = 10, L2 = 10, alpha = 0, A1 = 1, A2 = 1, Z = 0),
    method="L-BFGS-B",
    lower = c(L1 = 0, L2 = 0, A1 = 0, A2 = 0),
    upper = c(L1 = 200, L2 = 200, A1 = 1000, A2 = 1000))
  ## WITHOUT Hybridization effect on the mean load
  myFitNoAlphaH3 <- mle2(
    response ~ dnbinom(mu = MeanLoad(L1, L2, alpha, HI),
                       size = 1/abs(Aggregation(A1, A2, Z, HI))),
    data = data,
    fixed = list(alpha = 0),
    parameters = list(L1~Sex, L2~Sex, A1~Sex, A2~Sex, Z~Sex),
    start = list(L1 = 10, L2 = 10, alpha = 0, A1 = 1, A2 = 1, Z = 0),
    method="L-BFGS-B",
    lower = c(L1 = 0, L2 = 0, A1 = 0, A2 = 0),
    upper = c(L1 = 200, L2 = 200, A1 = 1000, A2 = 1000))
  return(list(H0 = list(fitNoAlpha = myFitNoAlphaH0,
                        predictNoAlpha = predict(myFitNoAlphaH0),
                        fitAlpha = myFitAlphaH0,
                        predictAlpha = predict(myFitAlphaH0)),
              H3 = list(fitNoAlpha = myFitNoAlphaH3,
                        predictNoAlpha = predict(myFitNoAlphaH3),
                        fitAlpha = myFitAlphaH3,
                        predictAlpha = predict(myFitAlphaH3))))
}

system.time(ModelPinworm <- myFun(Joelle_data, HI, response = Aspiculuris.Syphacia))

ModelPinworm$fitAlpha

## Difference between with and without alpha?
pValueAlpha <- anova(ModelPinworm$fitNoAlpha, ModelPinworm$fitAlpha)[10]
LLdropIfNoAlpha <- as.numeric(logLik(ModelPinworm$fitNoAlpha) - logLik(ModelPinworm$fitAlpha))

paste("The p-value of the anova test when we remove alpha is", round(pValueAlpha, 5), 
      "It corresponds to a drop of likelihood of", round(LLdropIfNoAlpha,2))

## Plotting time!
ggplot() +
  geom_point(data = Joelle_data, aes(x = HI, y = log10(Aspiculuris.Syphacia + 1), color = Sex)) +
  geom_line(aes(x = HI, y = log10(ModelPinworm$predictAlpha + 1), color = Sex)) +
  geom_line(aes(x = HI, y = log10(ModelPinworm$predictNoAlpha +1), color = Sex), linetype="dotted") +
  scale_color_manual(values = c("red", "darkblue")) +
  scale_y_continuous(trans='log10')+
  theme_linedraw()