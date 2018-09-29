#' Prepare our data for fitting by group
#'
#' @param data A data frame
#' @param response A character string. Response (e.g. "worm_count")
#' @param group A character string. Which group is considered (e.g. "Sex")
#' @return A list containing the full dataset, a dataset for groupA, a dataset for groupB
#' @export

marshallData <- function (data, response, group) {
  # define group of 2 factors
  # data$group <- data[[group]]
  # it must be a group of 2 levels only
  if (length(levels(data[[group]])) != 2){
    stop("your group has too many levels, you must provide only 2")
  }
  dataForResponse <- data[complete.cases(data[[response]]),]
  dataForResponse_A <- data[data[[group]] == levels(data[[group]])[1],]
  dataForResponse_A <- dataForResponse_A[complete.cases(dataForResponse_A[[response]]),]
  dataForResponse_B <- data[data[[group]] == levels(data[[group]])[2],]
  dataForResponse_B <- dataForResponse_B[complete.cases(dataForResponse_B[[response]]),]
  return(list(
    all = dataForResponse,
    groupA = dataForResponse_A,
    groupB = dataForResponse_B
  ))
}
