source("Gtest.R")
source("plotBananas.R")

## Import data
HeitlingerFieldData <- read.csv("../../Data_important/FinalFullDF_flotationPcrqPCR.csv")
Flotation_data <- HeitlingerFieldData[!is.na(HeitlingerFieldData$OPG) &
                                        !is.na(HeitlingerFieldData$HI) &
                                        !is.na(HeitlingerFieldData$Sex), ]
# Works if OPG are integers
Flotation_data$OPG <- round(Flotation_data$OPG)

# changes to do for avoiding plotting error
levels(Flotation_data$Sex) <- c(levels(Flotation_data$Sex), "female", "male")
Flotation_data$Sex[Flotation_data$Sex == "F"] <- "female"
Flotation_data$Sex[Flotation_data$Sex == "M"] <- "male"

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

## First question: can we consider the zeros out? If our hypothesis is : all get infected the same way, 
# but hybrids got lower loads?
source("Models/fitBinomial.R")

Flotation_data$flotStatus <- 0
Flotation_data$flotStatus[Flotation_data$OPG > 0] <- 1

ggplot(Flotation_data, aes(x = HI, y = flotStatus)) +
  geom_point() +
  geom_smooth() +
  theme_classic()

# we test the hypothesis 
runAll <- function (data, response) {
  print(paste0("Fit for the response: ", response))
  defaultConfig <- list(optimizer = "optimx",
                        method = c("bobyqa", "L-BFGS-B"),
                        control = list(follow.on = TRUE))
  paramBounds <- c(L1start = mean(na.omit(data[[response]])), 
                   L1LB = min(na.omit(data[[response]])), 
                   L1UB = max(na.omit(data[[response]])), 
                   L2start = mean(na.omit(data[[response]])), 
                   L2LB = min(na.omit(data[[response]])), 
                   L2UB = max(na.omit(data[[response]])), 
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
  
  return(list(H0 = H0, H1 = H1, H2 = H2, H3 = H3))
}

## Run the fit
fit <- analyse(Flotation_data, "flotStatus")
# plot H0
bananaPlots(mod = fit$H0, data = Flotation_data, response = "flotStatus")

######## Choose a correct distribution for our data : negative binomial ########
# source the functions defining meanload and aggregation for the negative binomial
source("Models/fitNegBin.R")
opg <- Flotation_data[Flotation_data$OPG > 0, "OPG"]
hist(opg, 100)

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

fit_flotation <- analyse(Flotation_data, "OPG", 
                         paramBounds = giveParamBounds(Flotation_data, "OPG"))

# bananaPlots(mod = fit_flotation$H0, data = Flotation_data, response = "OPG")
# problems for the profiling...

# Let's focus ONLY on the infected individuals. Hyp : HI makes the LOAD vary, not the infection
fit_flotation_positive <- analyse(Flotation_data[Flotation_data$OPG > 0,], "OPG", 
                                  paramBounds = giveParamBounds(Flotation_data[Flotation_data$OPG > 0,], "OPG"))

# bananaPlots(mod = fit_flotation_positive$H0, data = Flotation_data[Flotation_data$OPG > 0,], response = "OPG")
# problems for the profiling...
