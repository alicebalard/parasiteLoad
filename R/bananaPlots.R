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
    ggplot2::ggplot() +
      ggplot2::geom_point(data = data, ggplot2::aes_string(x = "HI", y = "response"), fill = "lightgrey", pch = 21, size = 3) +
      ggplot2::scale_fill_manual(values = cols) +
      ggplot2::geom_ribbon(ggplot2::aes(x = bananaDF$HI, ymin = bananaDF$min, ymax = bananaDF$max),
                  fill = "grey", alpha = .5) +
      ggplot2::geom_line(ggplot2::aes(x = bananaDF$HI, y = bananaDF$fit)) +
      ggplot2::theme_classic(base_size = 20) + {
        if(islog10 == TRUE) ggplot2::scale_y_log10()
        } +
      ggplot2::ylab(label = response)
  }else{
    bananaDF = data.frame(HI = numeric(), fit = numeric(), min = numeric(), max = numeric(), group = factor())
    mygroupDFnames = vector()
    for (n in names(mod)){
      thisbanana = getBananaDF(mod[[n]], hybridIndex)
      thisbanana$group <- n
      bananaDF <- rbind(bananaDF, thisbanana)
    }
    # Draw the line for the parameters at their MLE, alpha varying
    ggplot2::ggplot() +
      ggplot2::geom_point(data = data, ggplot2::aes_string(x = "HI", y = "response"), pch = 21, size = 3) +
      ggplot2::geom_ribbon(ggplot2::aes(x = bananaDF$HI, ymin = bananaDF$min, ymax = bananaDF$max, fill = bananaDF$group),
                  alpha = .2) +
      ggplot2::geom_line(ggplot2::aes(x = bananaDF$HI, y = bananaDF$fit, col = bananaDF$group)) +
      ggplot2::scale_fill_manual(values = cols) +
      ggplot2::scale_color_manual(values = cols) +
      ggplot2::theme_classic(base_size = 20) + {
        if(islog10 == TRUE) ggplot2::scale_y_log10()
      } +
      ggplot2::ylab(label = response)
  }
}
