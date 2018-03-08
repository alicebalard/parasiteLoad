## source the functions defining meanload and aggregation for the negative binomial
source("Models/BCI_qPCR-NormalDistrib.R")
source("MLE_hybrid_functions.R")

# if (!require("devtools")) install.packages("devtools")
# devtools::install_github("sjmgarnier/viridis")
library(rms)
# First statistical model
# 
# model0 <- lm(formula = BCI ~ HI, data = data4stats)
# 
# model1 <- lm(formula = BCI ~ HI * EimeriaDetected, data = data4stats)
# 
# model2 <- lm(formula = BCI ~ HI * EimeriaDetected * Sex, data = data4stats)
# 
# anova(model0, model1) # So we choose model 0
# anova(model1, model2)
# 
# summary(model1)
# 
# # In rms vocabulary
# dd <- datadist(data4stats)
# options(datadist="dd")
# 
# fit <- ols(formula = BCI ~ HI * EimeriaDetected, data = data4stats, x = TRUE, y = TRUE)
# 
# p <- Predict(fit, HI, EimeriaDetected)
# 
# ggplot(p) +
#   geom_point(data = miceTable, aes(x = HI, y = BCI)) +
#   coord_cartesian(ylim = c(min(na.omit(miceTable$BCI)) - 0.005,
#                            max(na.omit(miceTable$BCI)) + 0.005)) +
#   theme_bw()

## Import data 

eimeria_detect <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/raw_data/Eimeria_detection/Summary_eimeria.csv")
miceTable <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/raw_data/MiceTable_2014to2017.csv")

miceTable <- merge(miceTable, eimeria_detect)

data4stats <- miceTable[names(miceTable) %in% c("BCI", "HI", "EimeriaDetected", "Sex")]
data4stats <- na.omit(data4stats)
data4stats$EimeriaDetected <- factor(data4stats$EimeriaDetected)

qplot(data4stats$BCI) + theme_bw()

#### Our model
marshallData <- function (data, response) {
  dataForResponse <- data[complete.cases(data[[response]]),]
  dataForResponse_T <- dataForResponse[dataForResponse$EimeriaDetected == TRUE,]
  dataForResponse_F <- dataForResponse[dataForResponse$EimeriaDetected == FALSE,]
  dataForResponse_T$EimeriaDetected <- droplevels(dataForResponse_T$EimeriaDetected)
  dataForResponse_F$EimeriaDetected <- droplevels(dataForResponse_F$EimeriaDetected)
  return(list(
    all = dataForResponse,
    true = dataForResponse_T,
    false = dataForResponse_F
  ))
}

runAll <- function (data, response) {
  print(paste0("Fit for the response: ", response))
  defaultConfig <- list(optimizer = "optimx",
                        method = c("bobyqa", "L-BFGS-B"),
                        control = list(follow.on = TRUE))
  paramBounds <- c(L1start = 0.5, L1LB = 0, L1UB = 1, 
                   L2start = 0.5, L2LB = 0, L2UB = 1, 
                   alphaStart = 0, alphaLB = -5, alphaUB = 5)
  marshalledData <- marshallData(data, response)
  print("Fitting for all")
  FitAll <- run(
    data = marshalledData$all,
    response = response,
    hybridIndex = HI, 
    paramBounds = paramBounds, 
    config = defaultConfig
  )
  print("Fitting for true")
  FitTrue <- run(
    data = marshalledData$true,
    response = response,
    hybridIndex = HI, 
    paramBounds = paramBounds,
    config = defaultConfig
  )
  print("Fitting for false")
  FitFalse <- run(
    data = marshalledData$false,
    response = response,
    hybridIndex = HI, 
    paramBounds = paramBounds,
    config = defaultConfig
  )
  return(list(FitAll = FitAll, FitTrue = FitTrue, FitFalse = FitFalse))
}

analyse <- function(data, response) {
  print(paste0("Analysing data for response: ", response))
  FitForResponse <- runAll(data, response)
  
  ####### Is alpha significant for each hypothesis?
  
  # H0: the expected load for the subspecies and between sexes is the same
  print("Testing H0 no alpha vs alpha")
  Gtest(model0 = FitForResponse$FitAll$fitBasicNoAlpha, 
        model1 = FitForResponse$FitAll$fitBasicAlpha)
  H0 <- FitForResponse$FitAll$fitBasicAlpha
  
  # H1: the mean load across sexes is the same, but can differ across subspecies
  print("Testing H1 no alpha vs alpha")
  Gtest(model0 = FitForResponse$FitAll$fitAdvancedNoAlpha, 
        model1 = FitForResponse$FitAll$fitAdvancedAlpha)
  
  H1 <- FitForResponse$FitAll$fitAdvancedAlpha
  
  # H2: the mean load across subspecies is the same, but can differ between the sexes
  print("Testing H2 true no alpha vs alpha")
  Gtest(model0 = FitForResponse$FitTrue$fitBasicNoAlpha, 
        model1 = FitForResponse$FitTrue$fitBasicAlpha)
  
  print("Testing H2 false no alpha vs alpha")
  Gtest(model0 = FitForResponse$FitFalse$fitBasicNoAlpha, 
        model1 = FitForResponse$FitFalse$fitBasicAlpha)
  
  H2 <- list(true = FitForResponse$FitTrue$fitBasicAlpha,
             false = FitForResponse$FitFalse$fitBasicAlpha)
  
  # H3: the mean load can differ both across subspecies and between sexes
  print("Testing H3 true no alpha vs alpha")
  Gtest(model0 = FitForResponse$FitTrue$fitAdvancedNoAlpha, 
        model1 = FitForResponse$FitTrue$fitAdvancedAlpha)
  
  print("Testing H3 false no alpha vs alpha")
  Gtest(model0 = FitForResponse$FitFalse$fitAdvancedNoAlpha, 
        model1 = FitForResponse$FitFalse$fitAdvancedAlpha)
  
  H3 <- list(true = FitForResponse$FitTrue$fitAdvancedAlpha,
             false = FitForResponse$FitFalse$fitAdvancedAlpha)
  
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

fit <- analyse(data4stats, "BCI")

plotAll(mod = fit$H1, data = data4stats, response = "BCI", CI = FALSE)
plot2sexes(modF = fit$H3$true, modM = fit$H3$false, data = data4stats, response = "BCI", CI = FALSE)
