source("Gtest.R")
source("plotBananas.R")

## Import data WATWM
Joelle_data <- read.csv("../data/EvolutionFinalData.csv")
Joelle_data <- Joelle_data[complete.cases(Joelle_data$HI),]
levels(Joelle_data$Sex) <- c(levels(Joelle_data$Sex), "female", "male")
Joelle_data$Sex[Joelle_data$Sex == "F"] <- "female"
Joelle_data$Sex[Joelle_data$Sex == "M"] <- "male"

## Separate in all, male, female the data frames
marshallData <- function (data, response) {
  dataForResponse <- data[complete.cases(data[[response]]),]
  dataForResponse_F <- dataForResponse[dataForResponse$Sex == "female",]
  dataForResponse_M <- dataForResponse[dataForResponse$Sex == "male",]
  dataForResponse_F$Sex <- droplevels(dataForResponse_F$Sex)
  dataForResponse_M$Sex <- droplevels(dataForResponse_M$Sex)
  return(list(
    all = dataForResponse,
    female = dataForResponse_F,
    male = dataForResponse_M
  ))
}

######## Choose a correct distribution for our data : negative binomial ########
# source the functions defining meanload and aggregation for the negative binomial
source("Models/fitNegBin.R")

# Fit per all, female or male our model
runAll <- function (data, response, paramBounds) {
  print(paste0("Fit for the response: ", response))
  defaultConfig <- list(optimizer = "optimx",
                        method = c("bobyqa", "L-BFGS-B"),
                        control = list(follow.on = TRUE))
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

# compare the inner hypotheses (with alpha or without) and the nested hypotheses
analyse <- function(data, response, paramBounds) {
  print(paste0("Analysing data for response: ", response))
  FitForResponse <- runAll(data, response, paramBounds)
  
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
  
  return(list(H0 = H0, H1 = H1, H2 = H2, H3 = H3))
}

## Run the analysis
giveParamBounds <- function(data, response){
  return(c(L1start = mean(na.omit(data[[response]])), 
           L1LB = 0, 
           L1UB = max(na.omit(data[[response]])),
           L2start = mean(na.omit(data[[response]])), 
           L2LB = 0,
           L2UB = max(na.omit(data[[response]])),
           alphaStart = 0, alphaLB = -5, alphaUB = 5, 
           A1start = 10, A1LB = 1e-9, A1UB = 1000,
           A2start = 10, A2LB = 1e-9, A2UB = 1000,
           Zstart = 0, ZLB = -20, ZUB = 20))
}

# pinworms (A. tetraptera and S. obvelata (Joelle_data$Aspiculuris.Syphacia))
AspJo <- analyse(Joelle_data, "Aspiculuris.Syphacia", 
                 giveParamBounds(Joelle_data, "Aspiculuris.Syphacia"))
bananaPlots(mod = AspJo$H1, data = Joelle_data, response = "Aspiculuris.Syphacia", islog10 = TRUE) # error in Hessian
# alpha = 1.39 in WATWM paper
coef(AspJo$H1)[["alpha"]]

# Trichuris muris (whipworm (Joelle_data$Trichuris))
TriJo <- analyse(Joelle_data, "Trichuris", 
                 giveParamBounds(Joelle_data, "Trichuris"))
logLik(TriJo$H0)
# bananaPlots(mod = TriJo$H2, data = Joelle_data, response = "Trichuris") # error in Hessian

# Taenia taeniaeformis (tapeworm (Joelle_data$Taenia))
TaeJo <- analyse(Joelle_data, "Taenia",
                 giveParamBounds(Joelle_data, "Taenia"))

# Mastophorus muris (Joelle_data$Mastophorus)
MasJo <- analyse(Joelle_data, "Mastophorus",
                 giveParamBounds(Joelle_data, "Mastophorus"))

# TriJe <- analyse(Jenny_data, "Trichuris", 
#                  giveParamBounds(Jenny_data, "Trichuris"))
# AspJe <- analyse(Jenny_data, "Aspiculuris_Syphacia",
#                  giveParamBounds(Jenny_data, "Aspiculuris_Syphacia"))
# TaeJe <- analyse(Jenny_data, "Taenia",
#                  giveParamBounds(Jenny_data, "Taenia"))
# MasJe <- analyse(Jenny_data, "Mastophorus",
#                  giveParamBounds(Jenny_data, "Mastophorus"))
