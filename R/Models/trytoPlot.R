## Get bananaDF for a given model
getBananaDF <- function(mod){
  ## profile investigates behavior of objective function near the MLE
  myProf <- profile(submod)
  ## Marginal confidence interval
  myConfInt <- confint(myProf)
  # if no L2 calculated, set L1
  if("L2" %in% rownames(myConfInt) == FALSE){
    myConfInt <- rbind(myConfInt, myConfInt[rownames(myConfInt) %in% "L1"])
    rownames(myConfInt)[rownames(myConfInt) %in% ""] <- "L2"
    # coef(mod)[["L2"]] <- coef(mod)[["L1"]]
  }
  ## Get marginal confidence interval for a given parameter
  getInf <- function(paramname){
    myConfInt[rownames(myConfInt) == paramname][1]}
  getSup <- function(paramname){
    myConfInt[rownames(myConfInt) == paramname][2]}
  # Optimisation of the ExpectedLoad (inf and sup) for all parameters varying in their confidence interval
  expectedResponse <- function(v, hybridIndex){
    L1 = v[1]; L2 = v[2]; alpha = v[3]
    heterozygoty <- 2 * hybridIndex * (1 - hybridIndex)
    expectedResponse <- (L1 + (L2 - L1) * hybridIndex) * (1 - alpha * heterozygoty)
    return(expectedResponse)
  }
  # Run over HI values and optimise max and min for each
  bananaDF = data.frame(HI = NA, min = NA, max = NA)
  for(hybridIndex in seq(0,1, 0.1)){
    maxLoad <- optim(par = c(L1 = getMax("L1"), L2 = getMax("L2"),
                             alpha = mean(getMax("alpha") - getMin("alpha"))),
                     fn = expectedResponse,
                     lower = c(L1 = getInf("L1"), L2 = getInf("L2"), alpha = getInf("alpha")),
                     upper = c(L1 = getSup("L1"), L2 = getSup("L2"), alpha = getSup("alpha")),
                     method = "L-BFGS-B",
                     control = list(fnscale=-1),
                     hybridIndex = hybridIndex) # maximize
    minLoad <- optim(par = c(L1 = getMax("L1"),
                             L2 = getMax("L1"),
                             alpha = mean(getMax("alpha") - getMin("alpha"))),
                     fn = expectedResponse,
                     lower = c(L1 = getInf("L1"), L2 = getInf("L2"), alpha = getInf("alpha")),
                     upper = c(L1 = getSup("L1"), L2 = getSup("L2"), alpha = getSup("alpha")),
                     method = "L-BFGS-B",
                     hybridIndex = hybridIndex)
    bananaDF = rbind(bananaDF, data.frame(HI = hybridIndex, min = minLoad$value, max = maxLoad$value))
  }
  bananaDF <- bananaDF[-1,]
  return(bananaDF)
}

## Plot function
bananaPlots <- function(mod, data, response, cols = c("red", "blue"), mygroup = "Sex", isLog10 = F){
  data$response <- data[[response]]
  getBananaDF(mod)
  # Draw the line for the parameters at their MLE, alpha varying 
  ggplot() + 
    geom_point(data = data, aes_string(x = "HI", y = "response", fill = mygroup), pch = 21, size = 3) + 
    scale_fill_manual(values = cols) +
    geom_ribbon(aes(x = bananaDF$HI, ymin = bananaDF$min, ymax = bananaDF$max),
                fill = "grey", alpha = .5) +
    geom_line(aes(x = bananaDF$HI, y = bananaDF$fit)) +
    theme_bw(base_size = 20) +
    ylab(label = response) +
    geom_text(aes(.5, bananaDF$fit[bananaDF$HI == .5] +1,
                  label = paste("alpha = ", round(coef(mod)["alpha"], 2))), cex = 5) 
}  

bananaPlots(mod = fit$H0, data = data4stats, response = "delta_ct_MminusE")


getBananaDF <- function(mod){
  ## profile investigates behavior of objective function near the MLE
  myProf <- profile(submod)
  ## Marginal confidence interval
  myConfInt <- confint(myProf)
  # if no L2 calculated, set L1
  if("L2" %in% rownames(myConfInt) == FALSE){
    myConfInt <- rbind(myConfInt, myConfInt[rownames(myConfInt) %in% "L1"])
    rownames(myConfInt)[rownames(myConfInt) %in% ""] <- "L2"
    coef(mod) <- c(coef(mod), L2 = coef(mod)[["L1"]])
  }
  return(coef(mod))
}
getBananaDF(fit$H0)
bananaPlots(mod = fit$H1, data = data4stats, response = "delta_ct_MminusE")


bananaDF$fit <- expectedResponse(v = c(coef(mod)["L1"], coef(mod)["L2"], coef(mod)["alpha"]),
                                 hybridIndex = seq(0,1,.1))
bananaPlots(mod = fit$H2, data = data4stats, response = "delta_ct_MminusE")

fit$H2$male

## Plot function
bananaPlots <- function(mod, data, response, cols = c("red", "blue"), mygroup = "Sex", isLog10 = F){
  data$response <- data[[response]]
  if(!is.list(mod)){ # case of no group separation
    # what we have UP
  } else {
    for (i in length(names(mod))){
      submod = mod[[names(mod)[i]]]
      
      # Draw the line for the parameters at their MLE, alpha varying 
      ggplot() + 
        geom_point(data = data, aes_string(x = "HI", y = "response", fill = mygroup), pch = 21, size = 3) + 
        scale_fill_manual(values = cols) +
        geom_ribbon(aes(x = bananaDF$HI, ymin = bananaDF$min, ymax = bananaDF$max),
                    fill = "grey", alpha = .5) +
        geom_line(aes(x = bananaDF$HI, y = bananaDF$fit)) +
        theme_bw(base_size = 20) +
        ylab(label = response) +
        geom_text(aes(.5, bananaDF$fit[bananaDF$HI == .5] +1,
                      label = paste("alpha = ", round(coef(mod)["alpha"], 2))), cex = 5) 
      
    }
    
    
    
    
  }
}









get(paste0(deparse(substitute(fit$H2)), "[[\", N, "\"]]"))

var <- "mpg"
#Doesn't work
mtcars$var
#These both work, but note that what they return is different
# the first is a vector, the second is a data.frame
mtcars[[var]]
mtcars[var]

  fit$H2[["female"]]
  
  ## profile investigates behavior of objective function near the MLE
  myProf <- profile(mod)
  ## Marginal confidence interval
  myConfInt <- confint(myProf)
  # if no L2 calculated, set L1
  if("L2" %in% rownames(myConfInt) == FALSE){
    myConfInt <- rbind(myConfInt, myConfInt[rownames(myConfInt) %in% "L1"])
    rownames(myConfInt)[rownames(myConfInt) %in% ""] <- "L2"
  }
  ## Get marginal confidence interval for a given parameter
  getInf <- function(paramname){
    myConfInt[rownames(myConfInt) == paramname][1]}
  getSup <- function(paramname){
    myConfInt[rownames(myConfInt) == paramname][2]}
  # Optimisation of the ExpectedLoad (inf and sup) for all parameters varying in their confidence interval
  expectedResponse <- function(v, hybridIndex){
    L1 = v[1]; L2 = v[2]; alpha = v[3]
    heterozygoty <- 2 * hybridIndex * (1 - hybridIndex)
    expectedResponse <- (L1 + (L2 - L1) * hybridIndex) * (1 - alpha * heterozygoty)
    return(expectedResponse)
  }
  # Run over HI values and optimise max and min for each
  bananaDF = data.frame(HI = NA, min = NA, max = NA)
  for(hybridIndex in seq(0,1, 0.1)){
    maxLoad <- optim(par = c(L1 = getMax("L1"), L2 = getMax("L2"),
                             alpha = mean(getMax("alpha") - getMin("alpha"))),
                     fn = expectedResponse,
                     lower = c(L1 = getInf("L1"), L2 = getInf("L2"), alpha = getInf("alpha")),
                     upper = c(L1 = getSup("L1"), L2 = getSup("L2"), alpha = getSup("alpha")),
                     method = "L-BFGS-B",
                     control = list(fnscale=-1),
                     hybridIndex = hybridIndex) # maximize
    minLoad <- optim(par = c(L1 = getMax("L1"),
                             L2 = getMax("L1"),
                             alpha = mean(getMax("alpha") - getMin("alpha"))),
                     fn = expectedResponse,
                     lower = c(L1 = getInf("L1"), L2 = getInf("L2"), alpha = getInf("alpha")),
                     upper = c(L1 = getSup("L1"), L2 = getSup("L2"), alpha = getSup("alpha")),
                     method = "L-BFGS-B",
                     hybridIndex = hybridIndex)
    bananaDF = rbind(bananaDF, data.frame(HI = hybridIndex, min = minLoad$value, max = maxLoad$value))
  }
  bananaDF <- bananaDF[-1,]
  bananaDF$fit <- expectedResponse(v = c(coef(mod)["L1"], coef(mod)["L2"], coef(mod)["alpha"]),
                                   hybridIndex = seq(0,1,.1))
  # Draw the line for the parameters at their MLE, alpha varying 
  ggplot() + 
    geom_point(data = data, aes_string(x = "HI", y = "response", fill = mygroup), pch = 21, size = 3) + 
    scale_fill_manual(values = cols) +
    geom_ribbon(aes(x = bananaDF$HI, ymin = bananaDF$min, ymax = bananaDF$max),
                fill = "grey", alpha = .5) +
    geom_line(aes(x = bananaDF$HI, y = bananaDF$fit)) +
    theme_bw(base_size = 20) +
    ylab(label = response) +
    geom_text(aes(.5, bananaDF$fit[bananaDF$HI == .5] +1,
                  label = paste("alpha = ", round(coef(mod)["alpha"], 2))), cex = 5) 
}  
