#' Run full analysis for the 4 hypotheses and results of G-test
#'
#' @param data A data frame
#' @param response A character string. Response (e.g. "worm_count")
#' @param model Method to be used in fitting the model
#' So far implemented for "binomial", "negbin", "student"
#' @param group A character string. Which group is considered (e.g. "Sex")
#' @param hybridIndex The hybrid index
#' @param myparamBounds A named vector giving the start, lower and upper parameters for optimization
#' @return The full result of the analysis after G-test tests
#' @export

analyse <- function(data, response, model, group,
                    hybridIndex = "HI",
                    myparamBounds = "default"){
  if (myparamBounds == "default"){
    paramBounds <- getParamBounds(model, data, response)
  } else {
    paramBounds <- myparamBounds
  }
  print(paste0("Analysing data for response: ", response))
  FitForResponse <- runAll(data, response, model, group,
                           paramBounds = paramBounds,
                           hybridIndex = hybridIndex)

  ####### Is alpha significant for each hypothesis?

  # H0: the expected load for the subspecies and between 2 groups is the same
  GH0 = Gtest(model0 = FitForResponse$FitAll$fitBasicNoAlpha,
              model1 = FitForResponse$FitAll$fitBasicAlpha)
  print("Testing H0 no alpha vs alpha")
  GH0
  H0 <- list(fitH0 = FitForResponse$FitAll$fitBasicAlpha, Gtest = GH0)

  # H1: the mean load across 2 groups is the same, but can differ across subspecies
  GH1 = Gtest(model0 = FitForResponse$FitAll$fitAdvancedNoAlpha,
              model1 = FitForResponse$FitAll$fitAdvancedAlpha)
  print("Testing H1 no alpha vs alpha")
  GH1

  H1 <- list(fitH1 = FitForResponse$FitAll$fitAdvancedAlpha, Gtest = GH1)

  # H2: the mean load across subspecies is the same, but can differ between the 2 groups
  GH2_GA = Gtest(model0 = FitForResponse$FitGroupA$fitBasicNoAlpha,
                 model1 = FitForResponse$FitGroupA$fitBasicAlpha)
  print("Testing H2 groupA no alpha vs alpha")
  GH2_GA

  GH2_GB = Gtest(model0 = FitForResponse$FitGroupB$fitBasicNoAlpha,
                 model1 = FitForResponse$FitGroupB$fitBasicAlpha)
  print("Testing H2 groupB no alpha vs alpha")
  GH2_GB

  H2 <- list(fitH2 = list(groupA = FitForResponse$FitGroupA$fitBasicAlpha,
                          groupB = FitForResponse$FitGroupB$fitBasicAlpha),
             Gtests = list(groupA = GH2_GA, groupB = GH2_GB))

  # H3: the mean load can differ both across subspecies and between 2 groups
  GH3_GA = Gtest(model0 = FitForResponse$FitGroupA$fitAdvancedNoAlpha,
                 model1 = FitForResponse$FitGroupA$fitAdvancedAlpha)
  print("Testing H3 groupA no alpha vs alpha")
  GH3_GA

  GH3_GB = Gtest(model0 = FitForResponse$FitGroupB$fitAdvancedNoAlpha,
                 model1 = FitForResponse$FitGroupB$fitAdvancedAlpha)
  print("Testing H3 groupB no alpha vs alpha")
  GH3_GB

  H3 <- list(fitH3 = list(groupA = FitForResponse$FitGroupA$fitAdvancedAlpha,
                          groupB = FitForResponse$FitGroupB$fitAdvancedAlpha),
             Gtests = list(groupA = GH3_GA, groupB = GH3_GB))

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
