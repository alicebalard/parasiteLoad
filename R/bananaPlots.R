#' Plot function
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

bananaPlots <- function(mod, data, response, hybridIndex = seq(0,1, 0.01),
                        cols = c("green", "orange"), group, islog10 = F){
  data$response = data[[response]]
  if(is.list(mod) == FALSE){ # we do not have differences between groups
    bananaDF = getBananaDF(mod, hybridIndex)
    # Draw the line for the parameters at their MLE, alpha varying
    ggplot() +
      geom_point(data = data, aes_string(x = "HI", y = "response"), fill = "lightgrey", pch = 21, size = 3) +
      scale_fill_manual(values = cols) +
      geom_ribbon(aes(x = bananaDF$HI, ymin = bananaDF$min, ymax = bananaDF$max),
                  fill = "grey", alpha = .5) +
      geom_line(aes(x = bananaDF$HI, y = bananaDF$fit)) +
      theme_classic(base_size = 20) + {
        if(islog10 == TRUE) scale_y_log10()
        } +
      ylab(label = response)
  }else{
    bananaDF = data.frame(HI = numeric(), fit = numeric(), min = numeric(), max = numeric(), group = factor())
    mygroupDFnames = vector()
    for (n in names(mod)){
      thisbanana = getBananaDF(mod[[n]], hybridIndex)
      thisbanana$group <- n
      bananaDF <- rbind(bananaDF, thisbanana)
    }
    # Draw the line for the parameters at their MLE, alpha varying
    ggplot() +
      geom_point(data = data, aes_string(x = "HI", y = "response"), pch = 21, size = 3) +
      geom_ribbon(aes(x = bananaDF$HI, ymin = bananaDF$min, ymax = bananaDF$max, fill = bananaDF$group),
                  alpha = .2) +
      geom_line(aes(x = bananaDF$HI, y = bananaDF$fit, col = bananaDF$group)) +
      scale_fill_manual(values = cols) +
      scale_color_manual(values = cols) +
      theme_classic(base_size = 20) + {
        if(islog10 == TRUE) scale_y_log10()
      } +
      ylab(label = response)
  }
}
