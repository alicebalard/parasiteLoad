#' Plot function : the fit (continuous line), hybrid effect variation (dotted line), and 95% CI of the fit, all parameters varying
#'
#' @param mod A model to be ploted
#' @param data A data frame
#' @param response A character string. Response (e.g. "worm_count")
#' @param hybridIndex A vector of points representing the index used as x axis
#' @param cols A vector of two colors for the two groups
#' @param group A character string. Which group is considered (e.g. "Sex")
#' @param islog10 logical. Is the y axis ploted with log10 transformation
#' @return A plot of our model
#' @export


bananaPlots <- function(mod, data, response, hybridIndex = seq(0,1, 0.05),
                              cols = wes_palette("IsleofDogs1")[c(1,3)], group, islog10 = F){
  # for aes_string
  data$response = data[[response]]
  data$group = data[[group]]

  # Create a dataframe for plotting
  getBananaDF <- function(mod, hybridIndex){

    ## Fitted coefficients
    fittedCoef <- bbmle::coef(mod)
    # if no L2 calculated, set L1
    if("L2" %in% names(fittedCoef) == FALSE){
      fittedCoef <- c(fittedCoef, fittedCoef[names(fittedCoef) %in% "L1"])
      names(fittedCoef)[length(names(fittedCoef))] <- "L2"
    }

    # Confidence interval for all coefficients based on quadratic approximation to the curvature at the maximum likelihood estimate
    myConfInt <- bbmle::confint(mod, method="quad")
    # if no L2 calculated, set L1
    if("L2" %in% rownames(myConfInt) == FALSE){
      myConfInt <- rbind(myConfInt[rownames(myConfInt) %in% "L1"], myConfInt)
      rownames(myConfInt)[1] <- "L2"
    }

    ## Get marginal confidence interval for a given parameter
    getInf <- function(paramname){
      myConfInt[rownames(myConfInt) == paramname][1]}
    getSup <- function(paramname){
      myConfInt[rownames(myConfInt) == paramname][2]}

    # Function to optimise for parameters varying in their confidance interval
    expectedResponse <- function(v, hybridIndex){
      L1 = v[1]; L2 = v[2]; alpha = v[3]
      heterozygoty <- 2 * hybridIndex * (1 - hybridIndex)
      expectedResponse <- (L1 + (L2 - L1) * hybridIndex) * (1 - alpha * heterozygoty)
      return(expectedResponse)
    }

    # create dataframe to fill up
    bananaDF = data.frame(HI = hybridIndex)

    # step 1 : calculate response expected by the model
    bananaDF$fit <- expectedResponse(c(fittedCoef[["L1"]],
                                       fittedCoef[["L2"]],
                                       fittedCoef[["alpha"]]),
                                     hybridIndex = hybridIndex)
    # step 2 : run over HI values and optimise max and min for each variable varying in their 95% CI
    bananaDF2 = data.frame(HI = numeric(), min = numeric(), max = numeric())
    for(i in hybridIndex){
      maxLoad <- optim(par = c(L1 = getSup("L1") - getInf("L1"),
                               L2 = getSup("L2") - getInf("L2"),
                               alpha = getSup("alpha") - getInf("alpha")),
                       fn = expectedResponse,
                       lower = c(L1 = getInf("L1"), L2 = getInf("L2"), alpha = getInf("alpha")),
                       upper = c(L1 = getSup("L1"), L2 = getSup("L2"), alpha = getSup("alpha")),
                       method = "L-BFGS-B",
                       control = list(fnscale=-1), # maximize
                       hybridIndex = i)
      minLoad <- optim(par = c(L1 = getSup("L1") - getInf("L1"),
                               L2 = getSup("L2") - getInf("L2"),
                               alpha = getSup("alpha") - getInf("alpha")),
                       fn = expectedResponse,
                       lower = c(L1 = getInf("L1"), L2 = getInf("L2"), alpha = getInf("alpha")),
                       upper = c(L1 = getSup("L1"), L2 = getSup("L2"), alpha = getSup("alpha")),
                       method = "L-BFGS-B",
                       hybridIndex = i)
      bananaDF2 = rbind(bananaDF2, data.frame(HI = i, min = minLoad$value, max = maxLoad$value))
    }
    bananaDF$minAll <- bananaDF2$min
    bananaDF$maxAll <- bananaDF2$max
    # step 3 : run over HI values and optimise max and min for hybrid effect varying in it's 95% CI and other variable fixed at fit
    bananaDF2 = data.frame(HI = numeric(), min = numeric(), max = numeric())
    for(i in hybridIndex){
      maxLoad <- optim(par = getSup("alpha") - getInf("alpha"),
                       fn = parasiteLoad::MeanLoad,
                       lower = getInf("alpha"),
                       upper = getSup("alpha"),
                       L1 = fittedCoef[["L1"]],
                       L2 = fittedCoef[["L2"]],
                       method = "L-BFGS-B",
                       control = list(fnscale=-1), # maximize
                       hybridIndex = i)
      minLoad <- optim(par = getSup("alpha") - getInf("alpha"),
                       fn = parasiteLoad::MeanLoad,
                       lower = getInf("alpha"),
                       upper = getSup("alpha"),
                       L1 = fittedCoef[["L1"]],
                       L2 = fittedCoef[["L2"]],
                       method = "L-BFGS-B", # minimize
                       hybridIndex = i)
      bananaDF2 = rbind(bananaDF2, data.frame(HI = i, min = minLoad$value, max = maxLoad$value))
    }
    bananaDF$minAlpha <- bananaDF2$min
    bananaDF$maxAlpha <- bananaDF2$max
    return(bananaDF)
  }

  if(is.list(mod) == FALSE){ # we do not have differences between groups
    bananaDFtoplot <- getBananaDF(mod = mod, hybridIndex = hybridIndex)
    bananaDFtoplot$group <- "all"
  } else {
    bananaList <- lapply(mod, FUN = getBananaDF, hybridIndex = hybridIndex)
    bananaDFA <- bananaList$groupA
    bananaDFB <- bananaList$groupB
    bananaDFA$group <- levels(data$group)[1]
    bananaDFB$group <- levels(data$group)[2]
    bananaDFtoplot <- rbind(bananaDFA, bananaDFB)
  }

  ## And now plot!
  p <- ggplot2::ggplot() +
    ggplot2::geom_ribbon(data = bananaDFtoplot,
                         ggplot2::aes_string(x = "HI", ymin = "minAll", ymax = "maxAll",
                                             group = "group"),
                         fill = "grey", alpha = .4) +
    ggplot2::geom_line(data = bananaDFtoplot, size = 2,
                       ggplot2::aes_string(x = "HI", y = "fit", col = "group")) +
    ggplot2::geom_line(data = bananaDFtoplot,
                       ggplot2::aes_string(x = "HI", y = "minAlpha", col = "group"),
                       linetype="dashed", size = 2) +
    ggplot2::geom_line(data = bananaDFtoplot,
                       ggplot2::aes_string(x = "HI", y = "maxAlpha", col = "group"),
                       linetype="dashed", size = 2) +
    ggplot2::geom_point(data = data, ggplot2::aes_string(x = "HI", y = "response", fill = "group"),
                        pch = 21, size = 3, alpha = .7) +
    ggplot2::scale_fill_manual(values = cols) +
    ggplot2::theme_classic(base_size = 20) + {
      if(islog10 == TRUE) ggplot2::scale_y_log10()
    } +
    ggplot2::ylab(label = response)

  if(is.list(mod) == TRUE){
    p <- p + ggplot2::scale_color_manual(values = cols)
  } else {
    p <- p + ggplot2::scale_color_manual(values = "black")
  }
  return(p)
}
