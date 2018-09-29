#' Run full analysis for the 4 hypotheses and results of G-test
#'
#' @param data A data frame
#' @param response A character string. Response (e.g. "worm_count")
#' @param model Method to be used in fitting the model
#' So far implemented for "binomial", "negbin", "student"
#' @param group A character string. Which group is considered (e.g. "Sex")
#' @return The full result of the analysis after G-test tests
#' @export

analyse <- function(data, response, model, group) {
  print(paste0("Analysing data for response: ", response))
  FitForResponse <- runAll(data, response, model, group)

  ####### Is alpha significant for each hypothesis?

  # H0: the expected load for the subspecies and between 2 groups is the same
  print("Testing H0 no alpha vs alpha")
  Gtest(model0 = FitForResponse$FitAll$fitBasicNoAlpha,
        model1 = FitForResponse$FitAll$fitBasicAlpha)
  H0 <- FitForResponse$FitAll$fitBasicAlpha

  # H1: the mean load across 2 groups is the same, but can differ across subspecies
  print("Testing H1 no alpha vs alpha")
  Gtest(model0 = FitForResponse$FitAll$fitAdvancedNoAlpha,
        model1 = FitForResponse$FitAll$fitAdvancedAlpha)

  H1 <- FitForResponse$FitAll$fitAdvancedAlpha

  # H2: the mean load across subspecies is the same, but can differ between the 2 groups
  print("Testing H2 groupA no alpha vs alpha")
  Gtest(model0 = FitForResponse$FitGroupA$fitBasicNoAlpha,
        model1 = FitForResponse$FitGroupA$fitBasicAlpha)

  print("Testing H2 groupB no alpha vs alpha")
  Gtest(model0 = FitForResponse$FitGroupB$fitBasicNoAlpha,
        model1 = FitForResponse$FitGroupB$fitBasicAlpha)

  H2 <- list(groupA = FitForResponse$FitGroupA$fitBasicAlpha,
             groupB = FitForResponse$FitGroupB$fitBasicAlpha)

  # H3: the mean load can differ both across subspecies and between 2 groups
  print("Testing H3 groupA no alpha vs alpha")
  Gtest(model0 = FitForResponse$FitGroupA$fitAdvancedNoAlpha,
        model1 = FitForResponse$FitGroupA$fitAdvancedAlpha)

  print("Testing H3 groupB no alpha vs alpha")
  Gtest(model0 = FitForResponse$FitGroupB$fitAdvancedNoAlpha,
        model1 = FitForResponse$FitGroupB$fitAdvancedAlpha)

  H3 <- list(groupA = FitForResponse$FitGroupA$fitAdvancedAlpha,
             groupB = FitForResponse$FitGroupB$fitAdvancedAlpha)

  ####### Compare the hypotheses with G-tests
  # H1 vs H0
  print("Testing H1 vs H0")
  Gtest(model0 = H0, model1 = H1)

  # H2 vs H0
  print("Testing H2 vs H0")
  Gtest(model0 = H0, model1 = H2)

  # H3 vs H1
  print("Testing H3 vs H1")
  Gtest(model0 = H1, model1 = H3)

  # H3 vs H2
  print("Testing H3 vs H2")
  Gtest(model0 = H2, model1 = H3)

  return(list(H0 = H0, H1 = H1, H2 = H2, H3 = H3))
}