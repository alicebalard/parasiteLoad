

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
      
      expectedResponse <- function(v, hybridIndex){
        L1 = v[1]; L2 = v[2]; alpha = v[3]
        heterozygoty <- 2 * hybridIndex * (1 - hybridIndex)
        expectedResponse <- (L1 + (L2 - L1) * hybridIndex) * (1 - alpha * heterozygoty)
        return(expectedResponse)
      }
      #test
      expectedResponse(c(1,2,.3), .4)
      # Run over HI values and optimise max and min for each
      bananaDF = data.frame(HI = NA, min = NA, max = NA)
      for(hybridIndex in seq(0,1, 0.1)){
        maxLoad <- optim(par = c(L1 = getMax("L1"),
                                 L2 = getMax("L1"),
                                 alpha = mean(getMax("alpha") - getMin("alpha"))),
                         fn = expectedResponse,
                         lower = c(L1 = getInf("L1"), L2 = getInf("L1"), alpha = getInf("alpha")),
                         upper = c(L1 = getSup("L1"), L2 = getSup("L1"), alpha = getSup("alpha")),
                         method = "L-BFGS-B",
                         control = list(fnscale=-1),
                         hybridIndex = hybridIndex) # maximize
        minLoad <- optim(par = c(L1 = getMax("L1"),
                                 L2 = getMax("L1"),
                                 alpha = mean(getMax("alpha") - getMin("alpha"))),
                         fn = expectedResponse,
                         lower = c(L1 = getInf("L1"), L2 = getInf("L1"), alpha = getInf("alpha")),
                         upper = c(L1 = getSup("L1"), L2 = getSup("L1"), alpha = getSup("alpha")),
                         method = "L-BFGS-B",
                         hybridIndex = hybridIndex)
        bananaDF = rbind(bananaDF, data.frame(HI = hybridIndex, min = minLoad$value, max = maxLoad$value))
      }
      bananaDF <- bananaDF[-1,]
      bananaDF$fit <- expectedResponse(v = c(coef(fit$H0)["L1"], coef(fit$H0)["L1"], coef(fit$H0)["alpha"]),
                                       hybridIndex = seq(0,1,.1))
      # Draw the line for the parameters at their MLE, alpha varying 
      ggplot() + 
        geom_point(data = data, aes_string(x = "HI", y = "response", fill = mygroup), pch = 21, size = 3) + 
        scale_fill_manual(values = cols) +
        geom_ribbon(aes(x = bananaDF$HI,
                        ymin = bananaDF$min,
                        ymax = bananaDF$max),
                    fill = "grey", alpha = .5) +
        geom_line(aes(x = bananaDF$HI, y = bananaDF$fit)) +
        theme_bw(base_size = 20) +
        ylab(label = labelfory) +
        geom_text(aes(.5, bananaDF$fit[bananaDF$HI == .5] +1,
                      label = paste("alpha = ", round(coef(fit$H0)["alpha"], 2)))) 
    
      
      
      
      
      