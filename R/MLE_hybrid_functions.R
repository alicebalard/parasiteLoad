library(ggplot2)

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

## Plot functions
plotAll <- function(mod, data, response){
  data$response <- data[[response]]
  ## profile investigates behavior of objective function near the MLE
  myProf <- profile(TriJo$H1)

  ## Marginal confidence interval
  myConfInt <- confint(myProf)

  ## Marginal confidence interval for alpha
  alphaCILB <- myConfInt[rownames(myConfInt) == "alpha"][1]
  alphaCIUB <- myConfInt[rownames(myConfInt) == "alpha"][2]

  ## Draw the line for the parameters at their MLE, alpha varying 
  DF <- data.frame(HI = seq(0,1,0.01), 
                   loadMLE = MeanLoad(L1 = coef(mod)[names(coef(mod)) == "L1"], 
                                      L2 =  coef(mod)[names(coef(mod)) == "L2"], 
                                      alpha =  coef(mod)[names(coef(mod)) == "alpha"],  
                                      hybridIndex = seq(0,1,0.01)),
                   loadMLEAlphaLB = MeanLoad(L1 = coef(mod)[names(coef(mod)) == "L1"],
                                             L2 =  coef(mod)[names(coef(mod)) == "L2"],
                                             alpha =  alphaCILB,
                                             hybridIndex = seq(0,1,0.01)),
                   loadMLEAlphaUB = MeanLoad(L1 = coef(mod)[names(coef(mod)) == "L1"],
                                             L2 =  coef(mod)[names(coef(mod)) == "L2"],
                                             alpha =  alphaCIUB,
                                             hybridIndex = seq(0,1,0.01)))
  ggplot() + 
    geom_point(data = data, aes(x = HI, y = log10(response + 1), color = Sex)) + 
    geom_ribbon(aes(x = DF$HI,
                    ymin = log10(DF$loadMLEAlphaUB + 1),
                    ymax = log10(DF$loadMLEAlphaLB + 1)),
                fill = "grey", alpha = .5) +
    geom_line(aes(x = DF$HI, y = log10(DF$loadMLE + 1))) +
    theme_bw()
}

#### Case 2 : 2 sexes
plot2sexes <- function(modF, modM, data, response){
  data$response <- data[[response]]
  ## profile investigates behavior of objective function near the MLE 
  myProfF <- profile(modF)
  myProfM <- profile(modM)
  
  ## Marginal confidence interval 
  myConfIntF <- confint(myProfF) 
  myConfIntM <- confint(myProfM) 
  
  ## Marginal confidence interval for alpha 
  alphaCILBF <- myConfIntF[rownames(myConfIntF) == "alpha"][1] 
  alphaCIUBF <- myConfIntM[rownames(myConfIntM) == "alpha"][2] 
  alphaCILBM <- myConfIntF[rownames(myConfIntF) == "alpha"][1] 
  alphaCIUBM <- myConfIntM[rownames(myConfIntM) == "alpha"][2] 
  
  ## Draw the line for the parameters at their MLE, alpha varying 
  DF <- data.frame(HI = seq(0,1,0.01), 
                   loadMLEF = MeanLoad(L1 = coef(modF)[names(coef(modF)) == "L1"], 
                                       L2 =  coef(modF)[names(coef(modF)) == "L2"], 
                                       alpha =  coef(modF)[names(coef(modF)) == "alpha"],  
                                       hybridIndex = seq(0,1,0.01)), 
                   loadMLEAlphaLBF = MeanLoad(L1 = coef(modF)[names(coef(modF)) == "L1"],
                                              L2 =  coef(modF)[names(coef(modF)) == "L2"], 
                                              alpha =  alphaCILB, 
                                              hybridIndex = seq(0,1,0.01)), 
                   loadMLEAlphaUBF = MeanLoad(L1 = coef(modF)[names(coef(modF)) == "L1"], 
                                              L2 =  coef(modF)[names(coef(modF)) == "L2"], 
                                              alpha =  alphaCIUB, 
                                              hybridIndex = seq(0,1,0.01)), 
                   loadMLEM = MeanLoad(L1 = coef(modM)[names(coef(modM)) == "L1"], 
                                       L2 =  coef(modM)[names(coef(modM)) == "L2"], 
                                       alpha =  coef(modM)[names(coef(modM)) == "alpha"],  
                                       hybridIndex = seq(0,1,0.01)), 
                   loadMLEAlphaLBM = MeanLoad(L1 = coef(modM)[names(coef(modM)) == "L1"],
                                              L2 =  coef(modM)[names(coef(modM)) == "L2"], 
                                              alpha =  alphaCILB, 
                                              hybridIndex = seq(0,1,0.01)), 
                   loadMLEAlphaUBM = MeanLoad(L1 = coef(modM)[names(coef(modM)) == "L1"], 
                                              L2 =  coef(modM)[names(coef(modM)) == "L2"], 
                                              alpha =  alphaCIUB, 
                                              hybridIndex = seq(0,1,0.01))) 
  
  ggplot() + 
    geom_point(data = data, aes(x = HI, y = log10(response + 1), color = Sex)) + 
    geom_ribbon(aes(x = DF$HI,  
                    ymin = log10(DF$loadMLEAlphaUBF + 1),  
                    ymax = log10(DF$loadMLEAlphaLBF + 1)), 
                fill = "pink", alpha = .5) + 
    geom_ribbon(aes(x = DF$HI,  
                    ymin = log10(DF$loadMLEAlphaUBM + 1),  
                    ymax = log10(DF$loadMLEAlphaLBM + 1)), 
                fill = "blue", alpha = .5) + 
    geom_line(aes(x = DF$HI, y = log10(DF$loadMLEF + 1)), col = "pink") + 
    geom_line(aes(x = DF$HI, y = log10(DF$loadMLEM + 1)), col = "blue") + 
    theme_bw()
}