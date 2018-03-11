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
  pvalue <- 1 - pchisq(2*dLL, df=dDF)
  out <- data.frame(dLL = round(dLL, 2),
                    dDF = dDF,
                    pvalue = 
                      if (pvalue < 0.01){
                      "< 0.01"
                    } else {
                      round(pvalue, 3)
                    }
  )
  print(out)
  return(out)
}

## Plot functions
plotAll <- function(mod, data, response, CI, cols = c("red", "blue"),
                    mygroup = "Sex", switchlevels = FALSE, labelfory = "log10 response +1"){
  data$response <- data[[response]]
  data$log10resp <- log10(data$response + 1)  
  if (switchlevels == TRUE){
    data[[mygroup]] <- factor(data[[mygroup]], levels(data[[mygroup]])[c(2,1)])
  }
  
  # Do we plot CI ?
  if (CI == TRUE){
    
    ## profile investigates behavior of objective function near the MLE
    myProf <- profile(mod)
    
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
      geom_point(data = data, aes_string(x = "HI", y = "log10resp", color = mygroup)) + 
      scale_color_manual(values = cols) +
      geom_ribbon(aes(x = DF$HI,
                      ymin = log10(DF$loadMLEAlphaUB + 1),
                      ymax = log10(DF$loadMLEAlphaLB + 1)),
                  fill = "grey", alpha = .5) +
      geom_line(aes(x = DF$HI, y = log10(DF$loadMLE + 1))) +
      theme_bw(base_size = 20)+
      ylab(label = labelfory)
  } else {
    ## Draw the line for the parameters at their MLE, alpha varying 
    DF <- data.frame(HI = seq(0,1,0.01), 
                     loadMLE = MeanLoad(L1 = coef(mod)[names(coef(mod)) == "L1"], 
                                        L2 =  coef(mod)[names(coef(mod)) == "L2"], 
                                        alpha =  coef(mod)[names(coef(mod)) == "alpha"],  
                                        hybridIndex = seq(0,1,0.01)))
    ggplot() + 
      geom_point(data = data, aes_string(x = "HI", y = "log10resp", color = mygroup)) + 
      scale_color_manual(values = cols) +
      geom_line(aes(x = DF$HI, y = log10(DF$loadMLE + 1))) +
      theme_bw(base_size = 20)+
      ylab(label = labelfory)
  }
}

#### Case 2 : 2 sexes
plot2sexes <- function(modF, modM, data, response, CI, cols = c("red", "blue"), 
                       mygroup = "Sex", switchlevels = FALSE, labelfory = "log10 response +1"){
  data$response <- data[[response]]
  data$log10resp <- log10(data$response + 1)
  if (switchlevels == TRUE){
    data[[mygroup]] <- factor(data[[mygroup]], levels(data[[mygroup]])[c(2,1)])
  }
  
  # Do we plot CI ?
  if (CI == TRUE){
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
                     loadMLEM = MeanLoad(L1 = coef(modM)[names(coef(modM)) == "L1"], 
                                         L2 =  coef(modM)[names(coef(modM)) == "L2"], 
                                         alpha =  coef(modM)[names(coef(modM)) == "alpha"],  
                                         hybridIndex = seq(0,1,0.01)), 
                     loadMLEAlphaLBF = MeanLoad(L1 = coef(modF)[names(coef(modF)) == "L1"],
                                                L2 =  coef(modF)[names(coef(modF)) == "L2"], 
                                                alpha =  alphaCILB, 
                                                hybridIndex = seq(0,1,0.01)), 
                     loadMLEAlphaUBF = MeanLoad(L1 = coef(modF)[names(coef(modF)) == "L1"], 
                                                L2 =  coef(modF)[names(coef(modF)) == "L2"], 
                                                alpha =  alphaCIUB, 
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
      geom_point(data = data, aes_string(x = "HI", y = "log10resp", color = mygroup)) + 
      scale_color_manual(values = cols) +
      geom_ribbon(aes(x = DF$HI,  
                      ymin = log10(DF$loadMLEAlphaUBF + 1),  
                      ymax = log10(DF$loadMLEAlphaLBF + 1)), 
                  fill = cols[1], alpha = .5) + 
      geom_ribbon(aes(x = DF$HI,  
                      ymin = log10(DF$loadMLEAlphaUBM + 1),  
                      ymax = log10(DF$loadMLEAlphaLBM + 1)), 
                  fill = cols[2], alpha = .5) + 
      geom_line(aes(x = DF$HI, y = log10(DF$loadMLEF + 1)), col = cols[1]) + 
      geom_line(aes(x = DF$HI, y = log10(DF$loadMLEM + 1)), col = cols[2]) + 
      theme_bw(base_size = 20)+
      ylab(label = labelfory)
  } else {
    ## Draw the line for the parameters at their MLE, alpha varying 
    DF <- data.frame(HI = seq(0,1,0.01), 
                     loadMLEF = MeanLoad(L1 = coef(modF)[names(coef(modF)) == "L1"], 
                                         L2 =  coef(modF)[names(coef(modF)) == "L2"], 
                                         alpha =  coef(modF)[names(coef(modF)) == "alpha"],  
                                         hybridIndex = seq(0,1,0.01)), 
                     loadMLEM = MeanLoad(L1 = coef(modM)[names(coef(modM)) == "L1"], 
                                         L2 =  coef(modM)[names(coef(modM)) == "L2"], 
                                         alpha =  coef(modM)[names(coef(modM)) == "alpha"],  
                                         hybridIndex = seq(0,1,0.01))) 
    ggplot() + 
      geom_point(data = data, aes_string(x = "HI", y = "log10resp", color = mygroup)) + 
      scale_color_manual(values = cols) +
      geom_line(aes(x = DF$HI, y = log10(DF$loadMLEF + 1)), col = cols[1], size = 3) + 
      geom_line(aes(x = DF$HI, y = log10(DF$loadMLEM + 1)), col = cols[2], size = 3) + 
      theme_bw(base_size = 20)+
      ylab(label = labelfory)
  }
}
