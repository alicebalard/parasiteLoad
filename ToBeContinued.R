# install.packages("bbmle")
library(bbmle)
library(ggplot2)
library(optimx)

## Import data WATWM
Joelle_data <- read.csv("examples/Reproduction_WATWM/EvolutionFinalData.csv")
# to check
Joelle_data <- Joelle_data[complete.cases(Joelle_data$HI),]
# pinworms (A. tetraptera and S. obvelata)
# Trichuris muris (whipworm)
dataTrichuris <- Joelle_data[complete.cases(Joelle_data$Trichuris),]
dataTrichuris_F <- dataTrichuris[dataTrichuris$Sex == "F",]
dataTrichuris_F$Sex <- droplevels(dataTrichuris_F$Sex)

# Taenia taeniaeformis (tapeworm) 

##### Input end #####

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
        fixed = fixed,
        optimizer = optimizer, method = method, control = control
      )
    }
    else {
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
        optimizer = optimizer, method = method, control = control
      ) 
    }
    convergence <- fit@details$convergence
    print(ifelse(convergence == 0, "Did converge", "Did not converge"))
    return(fit)
  }
  print("alpha")
  start <-  list(L1 = paramBounds[["L1start"]],
                 alpha = paramBounds[["alphaStart"]],
                 A1 = paramBounds[["A1start"]],
                 Z = paramBounds[["Zstart"]])
  myFitAlphaBasic <- MaxLikelihoodFunBasic(start, data)
  print("fixed alpha")
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
        optimizer = optimizer, method = method, control = control
      )
    }
    else {
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
        optimizer = optimizer, method = method, control = control
      ) 
    }
    convergence <- fit@details$convergence
    print(ifelse(convergence == 0, "Did converge", "Did not converge"))
    return(fit)
  }
  print("alpha")
  start <-  list(L1 = paramBounds[["L1start"]], L2 = paramBounds[["L2start"]],
                 alpha = paramBounds[["alphaStart"]],
                 A1 = paramBounds[["A1start"]], A2 = paramBounds[["A2start"]],
                 Z = paramBounds[["Zstart"]])
  myFitAlphaAdvanced <- MaxLikelihoodFunAdvanced(start, data)
  print("fixed alpha")
  myFitNoAlphaAdvanced <- MaxLikelihoodFunAdvanced(start, data, fixed = list(alpha = 0))
  ################## Output ##################
  return(list(fitNoAlphaBasic = myFitNoAlphaBasic,
              fitAlphaBasic = myFitAlphaBasic,
              fitNoAlphaAdvanced = myFitNoAlphaAdvanced,
              fitAlphaAdvanced = myFitAlphaAdvanced))
}

################## Data analysis ################## 
paramBounds <- c(L1start = 10, L1LB = 0, L1UB = 700, 
                 L2start = 10, L2LB = 0, L2UB = 700, 
                 alphaStart = 0, alphaLB = -5, alphaUB = 5,
                 A1start = 10, A1LB = 0, A1UB = 1000, 
                 A2start = 10, A2LB = 0, A2UB = 1000, 
                 Zstart = 0, ZLB = -5, ZUB = 5)
tTrichurisFit <- fit(data = dataTrichuris, response = "Trichuris",
                     hybridIndex = HI, paramBounds = paramBounds)

tTrichurisFit


anova(ModelTestWhipworm$fitNoAlpha, ModelTestWhipworm$fitAlpha)
logLik(ModelTestWhipworm$fitNoAlpha)
logLik(ModelTestWhipworm$fitAlpha)

## Draw the line for the parameters at their MLE, alpha varying

myPlot <- function(myFitAlpha, myFitNoAlpha, response, HI, Sex){
DF <- data.frame(HI = seq(0,1,0.01),
                 loadMLE = MeanLoad(L1 = coef(myFitAlpha)[names(coef(myFitAlpha)) == "L1"],
                                    L2 =  coef(myFitAlpha)[names(coef(myFitAlpha)) == "L2"],
                                    alpha =  coef(myFitAlpha)[names(coef(myFitAlpha)) == "alpha"], 
                                    hybridIndex = seq(0,1,0.01)),
                 loadMLEnoAlpha = MeanLoad(L1 = coef(myFitNoAlpha)[names(coef(myFitNoAlpha)) == "L1"],
                                           L2 =  coef(myFitNoAlpha)[names(coef(myFitNoAlpha)) == "L2"],
                                           alpha =  coef(myFitNoAlpha)[names(coef(myFitNoAlpha)) == "alpha"], 
                                           hybridIndex = seq(0,1,0.01)),
                 loadWATWM = MeanLoad(L1 = 1.19,
                                      L2 = 8.45,
                                      alpha = 0.54,
                                      hybridIndex = seq(0,1,0.01)))
ggplot() +
  geom_point(aes(x = 1- HI, y = log10(response + 1), color = Sex)) +
  geom_line(aes(x = 1- DF$HI, y = log10(DF$loadMLE + 1))) +
  geom_line(aes(x = 1- DF$HI, y = log10(DF$loadMLEnoAlpha + 1)), linetype="dotted") +
  geom_line(aes(x = 1- DF$HI, y = log10(DF$loadWATWM + 1)), color = "pink") +
  scale_color_manual(values = c("red", "darkblue")) +
  theme_linedraw()
}

myPlot(ModelTestWhipworm$fitAlpha, ModelTestWhipworm$fitNoAlpha, 
       Joelle_data$Trichuris, Joelle_data$HI, Joelle_data$Sex)



myprofle <- profile(ModelTestWhipworm$fitAlpha)








plot(y = log10(dataTrichuris_F$Trichuris +1), 
     x = 1 - dataTrichuris_F$HI)

range(dataTrichuris_F$Trichuris[dataTrichuris_F$HI < 0.3])

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

## Present in table
printModelAsTable <- function(model, modelName, interceptGroupName){
  modelTable <- t(data.frame(model@fullcoef))
  rownames(modelTable) <- modelName
  colNames <- colnames(modelTable)[grep("(Intercept)", colnames(modelTable))]
  parNames <- gsub('\\(Intercept\\)', "", colNames)
  for (parName in parNames){
    oldInterceptColname <- paste0(parName, "(Intercept)")
    newInterceptColname <- paste0(parName, interceptGroupName)
    colnames(modelTable)[colnames(modelTable) == oldInterceptColname]  <- newInterceptColname 
    interceptValue <- modelTable[,newInterceptColname]
    otherColumnsPattern <- paste0(parName, '(?!', interceptGroupName, ')')
    otherColumnsSelector <- grepl(otherColumnsPattern, colnames(modelTable), perl = TRUE)
    otherColumns <- colnames(modelTable)[otherColumnsSelector]
    for (columnToChange in otherColumns) {
      modelTable[,columnToChange] <- modelTable[,columnToChange] + interceptValue
    }
  }
  return(modelTable)
}

df1 <- printModelAsTable(ModelPinworm$H1$fitNoAlpha, "noAlpha", "SexF")
df2 <- printModelAsTable(ModelPinworm$H1$fitAlpha, "alpha", "SexF")
df3 <- printModelAsTable(ModelPinworm$H3$fitNoAlpha, "noAlpha", "SexF")
df4 <- printModelAsTable(ModelPinworm$H3$fitAlpha, "alpha", "SexF")

H1df <- gtools::smartbind(df1, df2)
H3df <- gtools::smartbind(df3, df4)

df1 <- printModelAsTable(ModelWhipworm$H1$fitNoAlpha, "noAlpha", "SexF")
df2 <- printModelAsTable(ModelWhipworm$H1$fitAlpha, "alpha", "SexF")
df3 <- printModelAsTable(ModelWhipworm$H3$fitNoAlpha, "noAlpha", "SexF")
df4 <- printModelAsTable(ModelPinworm$H3$fitAlpha, "alpha", "SexF")

H1df <- gtools::smartbind(df1, df2)
H3df <- gtools::smartbind(df3, df4)

# For a given model, downstream analyses and plot
finalFun <- function(models){
  extractpValue <- function(hypothesis){
    myAnova <- anova(hypothesis$fitNoAlpha, hypothesis$fitAlpha)
    mypValueAlpha <- myAnova[colnames(myAnova) == "Pr(>Chisq)"][2]
  }
  
  LLIncreaseIfAlpha <- function(hypothesis){
    as.numeric(logLik(hypothesis$fitAlpha) - logLik(hypothesis$fitNoAlpha))
  }
  
  isAlphaSignificant <- function(hypothesis, name){
    pValue <- extractpValue(hypothesis)
    LLIncrease <- LLIncreaseIfAlpha(hypothesis)
    print(
      paste0(
        "For ", name, " the p-value of the anova test when we add alpha is ", round(pValue, 3), 
        ". It corresponds to an increase of likelihood of ", round(LLIncrease,2), "."
      )
    )
    return(pValue < 0.05 & LLIncrease > 0)
  }
  
  H1 <- models$H1
  H3 <- models$H3
  
  isAlphaH1Significant <- isAlphaSignificant(H1, "H1")
  isAlphaH3Significant <- isAlphaSignificant(H3, "H3")
  
  selectBestHypothesis <- function(H1, isAlphaH1Significant, H3, isAlphaH3Significant) {
    if(isAlphaH1Significant) {
      print("Using H1 with alpha")
      H1ModelToTest <- H1$fitAlpha
    } else {
      print("Using H1 without alpha")
      H1ModelToTest <- H1$fitNoAlpha
    }
    if(isAlphaH3Significant) {
      print("Using H3 with alpha")
      H3ModelToTest <- H3$fitAlpha
    } else {
      print("Using H3 without alpha")
      H3ModelToTest <- H3$fitNoAlpha
    }
    H1toH3Anova <- anova(H1ModelToTest, H3ModelToTest)
    pValue <- H1toH3Anova[colnames(H1toH3Anova) == "Pr(>Chisq)"][2]
    LLIncrease <- as.numeric(logLik(H3ModelToTest) - logLik(H1ModelToTest))
    
    print(paste("The anova between H1 and H3 has a p-value of", round(pValue, 3)))
    print(paste("The likelihood increase between H1 and H3 is ", round(LLIncrease, 2)))
    isH3Better <- pValue < 0.05 & LLIncrease > 0
    
    if(isH3Better){
      print("Therefore we consider sex as a significant variable")
    } else {
      print("Therefore we keep the model without sex.")
    }
  }
  
  selectBestHypothesis(H1, isAlphaH1Significant, H3, isAlphaH3Significant)
}

finalFun(ModelPinworm)
finalFun(ModelWhipworm)
finalFun(ModelTapeworm)
finalFun(ModelMasto)

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



### Test plots
DF <- data.frame(HI = seq(0,1,0.01),
                 loadMLE = MeanLoad(L1 = coef(ModelWhipworm$H1$fitAlpha)[names(coef(ModelWhipworm$H1$fitAlpha)) == "L1"],
                                    L2 =  coef(ModelWhipworm$H1$fitAlpha)[names(coef(ModelWhipworm$H1$fitAlpha)) == "L2"],
                                    alpha =  coef(ModelWhipworm$H1$fitAlpha)[names(coef(ModelWhipworm$H1$fitAlpha)) == "alpha"], 
                                    hybridIndex = seq(0,1,0.01)),
                 loadMLEnoAlpha = MeanLoad(L1 = coef(ModelWhipworm$H1$fitNoAlpha)[names(coef(ModelWhipworm$H1$fitNoAlpha)) == "L1"],
                                           L2 =  coef(ModelWhipworm$H1$fitNoAlpha)[names(coef(ModelWhipworm$H1$fitNoAlpha)) == "L2"],
                                           alpha =  coef(ModelWhipworm$H1$fitNoAlpha)[names(coef(ModelWhipworm$H1$fitNoAlpha)) == "alpha"], 
                                           hybridIndex = seq(0,1,0.01))
)

ggplot() +
  geom_point(data = Joelle_data, aes(x = HI, y = log10(Trichuris + 1), color = Sex)) +
  # geom_ribbon(aes(x = DF$HI, 
  #                 ymin = log10(DF$loadMLEAlphaUB + 1), 
  #                 ymax = log10(DF$loadMLEAlphaLB + 1)),
  #             fill = "grey", alpha = .5) +
  geom_line(aes(x = DF$HI, y = log10(DF$loadMLE + 1))) +
  geom_line(aes(x = DF$HI, y = log10(DF$loadMLEnoAlpha + 1)), linetype="dotted") +
  scale_color_manual(values = c("red", "darkblue")) +
  theme_linedraw()

DF <- data.frame(HI = seq(0,1,0.01),
                 loadMLE_M = MeanLoad(L1 = coef(ModelWhipworm$H3$fitAlpha)[names(coef(ModelWhipworm$H3$fitAlpha)) == "L1.(Intercept)"] +
                                        coef(ModelWhipworm$H3$fitAlpha)[names(coef(ModelWhipworm$H3$fitAlpha)) == "L1.SexM"],
                                      L2 = coef(ModelWhipworm$H3$fitAlpha)[names(coef(ModelWhipworm$H3$fitAlpha)) == "L2.(Intercept)"] +
                                        coef(ModelWhipworm$H3$fitAlpha)[names(coef(ModelWhipworm$H3$fitAlpha)) == "L2.SexM"],
                                      alpha =   coef(ModelWhipworm$H3$fitAlpha)[names(coef(ModelWhipworm$H3$fitAlpha)) == "alpha.(Intercept)"] +
                                        coef(ModelWhipworm$H3$fitAlpha)[names(coef(ModelWhipworm$H3$fitAlpha)) == "alpha.SexM"],
                                      hybridIndex = seq(0,1,0.01)),
                 loadMLE_F = MeanLoad(L1 = coef(ModelWhipworm$H3$fitAlpha)[names(coef(ModelWhipworm$H3$fitAlpha)) == "L1.(Intercept)"],
                                      L2 = coef(ModelWhipworm$H3$fitAlpha)[names(coef(ModelWhipworm$H3$fitAlpha)) == "L2.(Intercept)"], 
                                      alpha =   coef(ModelWhipworm$H3$fitAlpha)[names(coef(ModelWhipworm$H3$fitAlpha)) == "alpha.(Intercept)"],
                                      hybridIndex = seq(0,1,0.01)),
                 loadMLE_M_noalpha = MeanLoad(L1 = coef(ModelWhipworm$H3$fitNoAlpha)[names(coef(ModelWhipworm$H3$fitNoAlpha)) == "L1.(Intercept)"] +
                                                coef(ModelWhipworm$H3$fitNoAlpha)[names(coef(ModelWhipworm$H3$fitNoAlpha)) == "L1.SexM"],
                                              L2 = coef(ModelWhipworm$H3$fitNoAlpha)[names(coef(ModelWhipworm$H3$fitNoAlpha)) == "L2.(Intercept)"] +
                                                coef(ModelWhipworm$H3$fitNoAlpha)[names(coef(ModelWhipworm$H3$fitNoAlpha)) == "L2.SexM"],
                                              alpha = 0,
                                              hybridIndex = seq(0,1,0.01)),
                 loadMLE_F_noalpha = MeanLoad(L1 = coef(ModelWhipworm$H3$fitNoAlpha)[names(coef(ModelWhipworm$H3$fitNoAlpha)) == "L1.(Intercept)"],
                                              L2 = coef(ModelWhipworm$H3$fitNoAlpha)[names(coef(ModelWhipworm$H3$fitNoAlpha)) == "L2.(Intercept)"], 
                                              alpha = 0,
                                              hybridIndex = seq(0,1,0.01)))

# loadMLEAlphaLB_M = MeanLoad(L1 = coef(ModelWhipworm$H1$fitAlpha)[names(coef(ModelWhipworm$H1$fitAlpha)) == "L1"],
#                             L2 =  coef(ModelWhipworm$H1$fitAlpha)[names(coef(ModelWhipworm$H1$fitAlpha)) == "L2"],
#                             alpha =  alphaCIUB_M,
#                             hybridIndex = seq(0,1,0.01)),
# loadMLEAlphaLB_F = MeanLoad(L1 = coef(ModelWhipworm$H1$fitAlpha)[names(coef(ModelWhipworm$H1$fitAlpha)) == "L1"],
#                             L2 =  coef(ModelWhipworm$H1$fitAlpha)[names(coef(ModelWhipworm$H1$fitAlpha)) == "L2"],
#                             alpha =  alphaCIUB_F,
#                             hybridIndex = seq(0,1,0.01)),
# loadMLEAlphaUB_M = MeanLoad(L1 = coef(ModelWhipworm$H1$fitAlpha)[names(coef(ModelWhipworm$H1$fitAlpha)) == "L1"],
#                             L2 =  coef(ModelWhipworm$H1$fitAlpha)[names(coef(ModelWhipworm$H1$fitAlpha)) == "L2"],
#                             alpha =  alphaCILB_M,
#                             hybridIndex = seq(0,1,0.01)),
# loadMLEAlphaUB_F = MeanLoad(L1 = coef(ModelWhipworm$H1$fitAlpha)[names(coef(ModelWhipworm$H1$fitAlpha)) == "L1"],
#                             L2 =  coef(ModelWhipworm$H1$fitAlpha)[names(coef(ModelWhipworm$H1$fitAlpha)) == "L2"],
#                             alpha =  alphaCILB_F,
#                             hybridIndex = seq(0,1,0.01)))


ggplot() +
  geom_point(data = Joelle_data, aes(x = HI, y = log10(Trichuris + 1), color = Sex)) +
  # geom_ribbon(aes(x = DF$HI, 
  #                 ymin = log10(DF$loadMLEAlphaUB + 1), 
  #                 ymax = log10(DF$loadMLEAlphaLB + 1)),
  #             fill = "grey", alpha = .5) +
  geom_line(aes(x = DF$HI, y = log10(DF$loadMLE_M + 1)), col = "blue") +
  geom_line(aes(x = DF$HI, y = log10(DF$loadMLE_F + 1)), col = "red") +
  #  geom_line(aes(x = DF$HI, y = log10(DF$loadMLEnoAlpha + 1)), linetype="dotted") +
  scale_color_manual(values = c("red", "darkblue")) +
  theme_linedraw()

