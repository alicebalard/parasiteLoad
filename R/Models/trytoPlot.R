scale = 2
shape = 4

data <- rweibull(n = 100, shape, scale)

scale

mean(data) / gamma(1 + 1/shape)

hist(data)
plot(1:100 ~ data)


mod = fit$H0
data = data4stats
response = "delta_ct_MminusE"
CI = T
labelfory = "delta_ct_MminusE"
isLog10 = F

mygroup = "Sex"
cols = c("red", "blue")
# ## Plot functions
# plotAll <- function(mod, data, response, CI, cols = c("red", "blue"),
#                     mygroup = "Sex", switchlevels = FALSE, labelfory = "log10 response +1", isLog10 = TRUE){
  data$response <- data[[response]]
  data$log10resp <- log10(data$response + 1)  
  if (switchlevels == TRUE){
    data[[mygroup]] <- factor(data[[mygroup]], levels(data[[mygroup]])[c(2,1)])
  }
  
  # Do we plot CI ?
  # if (CI == TRUE){
    
    ## profile investigates behavior of objective function near the MLE
    myProf <- profile(mod)
    
    ## Marginal confidence interval
    myConfInt <- confint(myProf)
    
    ## Get marginal confidence interval for a given parameter
    getMin <- function(paramname){
      myConfInt[rownames(myConfInt) == paramname][1]}
    getMax <- function(paramname){
      myConfInt[rownames(myConfInt) == paramname][2]}
    
    # Optimisation of the ExpectedLoad (min and max) for all parameters varying in their confidence interval
    # FitExpectedResponse <- function(data, response, config){
      print("Fitting expected response for H0")
      # data$response <- data[[response]] # little trick
      
      fit <- mle2(MeanLoad(L1, L2, alpha, HI),
                  data = data, 
        start = list(L1 = mean(getMax("L1") - getMin("L1")),
                     L2 = mean(getMax("L1") - getMin("L1")),
                     alpha = mean(getMax("alpha") - getMin("alpha"))),
        lower = c(L1 = getMin("L1"), L2 = getMin("L1"), alpha = getMin("alpha")),
        upper = c(L1 = getMax("L1"), L2 = getMax("L1"), alpha = getMax("alpha")),
        optimizer = config$optimizer, 
        method = config$method, 
        control = config$control)
      printConvergence(fit)
      return(fit)
    }
    
    MeanLoad(L1 = coef(mod)[names(coef(mod)) == "L1"], 
             L2 =  coef(mod)[names(coef(mod)) == "L2"], 
             alpha =  coef(mod)[names(coef(mod)) == "alpha"],  
             hybridIndex = seq(0,1,0.01))
    
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
    # Do we plot on log10 scale?
    if(isLog10 == TRUE){
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
      ggplot() + 
        geom_point(data = data, aes_string(x = "HI", y = "response", color = mygroup)) + 
        scale_color_manual(values = cols) +
        geom_ribbon(aes(x = DF$HI,
                        ymin = DF$loadMLEAlphaUB,
                        ymax = DF$loadMLEAlphaLB),
                    fill = "grey", alpha = .5) +
        geom_line(aes(x = DF$HI, y = DF$loadMLE)) +
        theme_bw(base_size = 20)+
        ylab(label = labelfory)
    }
  } else {
    ## Draw the line for the parameters at their MLE, alpha varying 
    DF <- data.frame(HI = seq(0,1,0.01), 
                     loadMLE = MeanLoad(L1 = coef(mod)[names(coef(mod)) == "L1"], 
                                        L2 =  coef(mod)[names(coef(mod)) == "L2"], 
                                        alpha =  coef(mod)[names(coef(mod)) == "alpha"],  
                                        hybridIndex = seq(0,1,0.01)))
    # Do we plot on log10 scale?
    if(isLog10 == TRUE){
      ggplot() + 
        geom_point(data = data, aes_string(x = "HI", y = "log10resp", color = mygroup)) + 
        scale_color_manual(values = cols) +
        geom_line(aes(x = DF$HI, y = log10(DF$loadMLE + 1))) +
        theme_bw(base_size = 20)+
        ylab(label = labelfory)
    } else {
      ggplot() + 
        geom_point(data = data, aes_string(x = "HI", y = "response", color = mygroup)) + 
        scale_color_manual(values = cols) +
        geom_line(aes(x = DF$HI, y = DF$loadMLE)) +
        theme_bw(base_size = 20)+
        ylab(label = labelfory)
    }
  }
}