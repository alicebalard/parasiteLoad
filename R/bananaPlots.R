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
  # for aes_string
  data$response = data[[response]]
  data$group = data[[group]]
  if(is.list(mod) == FALSE){ # we do not have differences between groups
    bananaDF = getBananaDF(mod, hybridIndex)
    bananaDF$group[bananaDF$group == "groupA"] <- levels(data$group)[1]
    bananaDF$group[bananaDF$group == "groupB"] <- levels(data$group)[2]
    # Draw the line for the parameters at their MLE, alpha varying
    ggplot2::ggplot() +
      ggplot2::geom_point(data = data, ggplot2::aes_string(x = "HI", y = "response", fill = "group"),
                          pch = 21, size = 3, alpha = .5) +
      ggplot2::geom_ribbon(data = bananaDF,
                           ggplot2::aes_string(x = "HI", ymin = "min", ymax = "max"),
                           fill = "grey", alpha = .5) +
      ggplot2::geom_line(data = bananaDF,
                         ggplot2::aes_string(x = "HI", y = "fit")) +
      ggplot2::scale_fill_manual(values = cols) +
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
    bananaDF$group[bananaDF$group == "groupA"] <- levels(data$group)[1]
    bananaDF$group[bananaDF$group == "groupB"] <- levels(data$group)[2]
    # Draw the line for the parameters at their MLE, alpha varying
    ggplot2::ggplot()+
      ggplot2::geom_point(data = data, ggplot2::aes_string(x = "HI", y = "response", fill = "group"),
                          pch = 21, size = 3, alpha = .5) +
      ggplot2::geom_ribbon(data = bananaDF,
                           ggplot2::aes_string(x = "HI", ymin = "min", ymax = "max",
                                               fill = "group"), alpha = .5) +
      ggplot2::geom_line(data = bananaDF,
                         ggplot2::aes_string(x = "HI", y = "fit", col = "group")) +
      ggplot2::scale_fill_manual(values = cols) +
      ggplot2::scale_color_manual(values = cols) +
      ggplot2::theme_classic(base_size = 20) + {
        if(islog10 == TRUE) ggplot2::scale_y_log10()
      } +
      ggplot2::ylab(label = response)
  }
}
