source("MLE_hybrid_functions.R")
# source the functions defining meanload and aggregation for the negative binomial
source("Models/MacroParasiteLoad-NegBin.R")


## Import data WATWM
Joelle_data <- read.csv("../data/EvolutionFinalData.csv")
Joelle_data <- Joelle_data[complete.cases(Joelle_data$HI),]

# pinworms (A. tetraptera and S. obvelata (Joelle_data$Aspiculuris.Syphacia))
# Trichuris muris (whipworm (Joelle_data$Trichuris))
# Taenia taeniaeformis (tapeworm (Joelle_data$Taenia))
# Mastophorus muris (Joelle_data$Mastophorus)

marshallData <- function (data, response) {
  dataForResponse <- data[complete.cases(data[[response]]),]
  dataForResponse_F <- dataForResponse[dataForResponse$Sex == "F",]
  dataForResponse_M <- dataForResponse[dataForResponse$Sex == "M",]
  dataForResponse_F$Sex <- droplevels(dataForResponse_F$Sex)
  dataForResponse_M$Sex <- droplevels(dataForResponse_M$Sex)
  return(list(
    all = dataForResponse,
    female = dataForResponse_F,
    male = dataForResponse_M
  ))
}

runAll <- function (data, response) {
  print(paste0("Fit for the response: ", response))
  defaultConfig <- list(optimizer = "optimx",
                        method = c("bobyqa", "L-BFGS-B"),
                        control = list(follow.on = TRUE))
  paramBounds <- c(L1start = 10, L1LB = 0, L1UB = 700, 
                   L2start = 10, L2LB = 0, L2UB = 700, 
                   alphaStart = 0, alphaLB = -5, alphaUB = 5,
                   A1start = 10, A1LB = 0, A1UB = 1000, 
                   A2start = 10, A2LB = 0, A2UB = 1000, 
                   Zstart = 0, ZLB = -5, ZUB = 5)
  marshalledData <- marshallData(data, response)
  print("Fitting for all")
  FitAll <- run(
    data = marshalledData$all,
    response = response,
    hybridIndex = HI, 
    paramBounds = paramBounds, 
    config = defaultConfig
  )
  print("Fitting for female")
  FitFemale <- run(
    data = marshalledData$female,
    response = response,
    hybridIndex = HI, 
    paramBounds = paramBounds,
    config = defaultConfig
  )
  print("Fitting for male")
  FitMale <- run(
    data = marshalledData$male,
    response = response,
    hybridIndex = HI, 
    paramBounds = paramBounds,
    config = defaultConfig
  )
  return(list(FitAll = FitAll, FitFemale = FitFemale, FitMale = FitMale))
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
  print("Testing H2 female no alpha vs alpha")
  Gtest(model0 = FitForResponse$FitFemale$fitBasicNoAlpha, 
        model1 = FitForResponse$FitFemale$fitBasicAlpha)
  
  print("Testing H2 male no alpha vs alpha")
  Gtest(model0 = FitForResponse$FitMale$fitBasicNoAlpha, 
        model1 = FitForResponse$FitMale$fitBasicAlpha)
  
  H2 <- list(female = FitForResponse$FitFemale$fitBasicAlpha,
             male = FitForResponse$FitMale$fitBasicAlpha)
  
  # H3: the mean load can differ both across subspecies and between sexes
  print("Testing H3 female no alpha vs alpha")
  Gtest(model0 = FitForResponse$FitFemale$fitAdvancedNoAlpha, 
        model1 = FitForResponse$FitFemale$fitAdvancedAlpha)
  
  print("Testing H3 male no alpha vs alpha")
  Gtest(model0 = FitForResponse$FitMale$fitAdvancedNoAlpha, 
        model1 = FitForResponse$FitMale$fitAdvancedAlpha)
  
  H3 <- list(female = FitForResponse$FitFemale$fitAdvancedAlpha,
             male = FitForResponse$FitMale$fitAdvancedAlpha)
  
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
}

analyse(Joelle_data, "Trichuris")
