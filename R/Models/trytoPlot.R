

mod = fit$H0

data = data4stats
response = "delta_ct_MminusE"
CI = T
labelfory = "delta_ct_MminusE"
isLog10 = F
config = list(optimizer = "optimx",
     method = c("bobyqa", "L-BFGS-B"),
     control = list(follow.on = TRUE))

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
    myProf
    ## Marginal confidence interval
    myConfInt <- confint(myProf)
    myConfInt
    
    ## Get marginal confidence interval for a given parameter
    getInf <- function(paramname){
      myConfInt[rownames(myConfInt) == paramname][1]}
    getSup <- function(paramname){
      myConfInt[rownames(myConfInt) == paramname][2]}
    
    # Optimisation of the ExpectedLoad (inf and sup) for all parameters varying in their confidence interval
    # FitExpectedResponse <- function(data, response, config){
      print("Fitting expected response for H0")
      # data$response <- data[[response]] # little trick
      
      toOptimise <- function(v){
        L1 = v[1]; L2 = v[2]; alpha = v[3]
        heterozygoty <- 2 * hybridIndex * (1 - hybridIndex)
        mean <- (L1 + (L2 - L1) * hybridIndex) * (1 - alpha * heterozygoty)
        return(mean)
      }
      
      fitted <- function(HI){
        hybridIndex = HI
        toOptimise(v = c(coef(fit$H0)["L1"], coef(fit$H0)["L1"], coef(fit$H0)["alpha"]))
      }
      
      fitted(3)
      
      
      # Run over HI values
      bananaDF = data.frame(HI = NA, min = NA, max = NA)
      for(hybridIndex in seq(0,1, 0.1)){
        maxLoad <- optim(par = c(L1 = getMax("L1"),
                                 L2 = getMax("L1"),
                                 alpha = mean(getMax("alpha") - getMin("alpha"))),
                         fn = toOptimise,
                         lower = c(L1 = getInf("L1"), L2 = getInf("L1"), alpha = getInf("alpha")),
                         upper = c(L1 = getSup("L1"), L2 = getSup("L1"), alpha = getSup("alpha")),
                         method = "L-BFGS-B",
                         control = list(fnscale=-1)) # maximize
        minLoad <- optim(par = c(L1 = getMax("L1"),
                                 L2 = getMax("L1"),
                                 alpha = mean(getMax("alpha") - getMin("alpha"))),
                         fn = toOptimise,
                         lower = c(L1 = getInf("L1"), L2 = getInf("L1"), alpha = getInf("alpha")),
                         upper = c(L1 = getSup("L1"), L2 = getSup("L1"), alpha = getSup("alpha")),
                         method = "L-BFGS-B")
        bananaDF = rbind(bananaDF, data.frame(HI = hybridIndex, min = minLoad$value, max = maxLoad$value))
      }
      
      fitDF <- data.frame(HI = seq(0, 1, 0.1), 
                          fit = toOptimise(v = c(coef(fit$H0)["L1"], coef(fit$H0)["L1"], coef(fit$H0)["alpha"])))

      
      fitted(0.2)
      
      
      # Draw the line for the parameters at their MLE, alpha varying 
      ggplot() + 
        geom_point(data = data, aes_string(x = "HI", y = "response", color = mygroup)) + 
        scale_color_manual(values = cols) +
        geom_ribbon(aes(x = bananaDF$HI,
                        ymin = bananaDF$min,
                        ymax = bananaDF$max),
                    fill = "grey", alpha = .5) #+
      
      
      
        geom_line(aes(x = bananaDF$HI, y = log10(DF$loadMLE + 1))) +
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
