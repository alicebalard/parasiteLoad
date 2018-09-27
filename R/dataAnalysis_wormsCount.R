source("Gtest.R")
source("plotBananas.R")
library(ggplot2)

## Import data WATWM
Joelle_data <- read.csv("../data/EvolutionFinalData.csv")
Joelle_data <- Joelle_data[complete.cases(Joelle_data$HI),]
levels(Joelle_data$Sex) <- c(levels(Joelle_data$Sex), "female", "male")
Joelle_data$Sex[Joelle_data$Sex == "F"] <- "female"
Joelle_data$Sex[Joelle_data$Sex == "M"] <- "male"

## Import our field data Jenny
HeitlingerFieldData <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/FinalFullDF_flotationPcrqPCR_Threshold_5.csv")
HeitlingerFieldData <- HeitlingerFieldData[!is.na(HeitlingerFieldData$Aspiculuris_Syphacia) &
                                             !is.na(HeitlingerFieldData$HI) &
                                             !is.na(HeitlingerFieldData$Sex),]
levels(HeitlingerFieldData$Sex) <- c(levels(HeitlingerFieldData$Sex), "female", "male")
HeitlingerFieldData$Sex[HeitlingerFieldData$Sex == "F"] <- "female"
HeitlingerFieldData$Sex[HeitlingerFieldData$Sex == "M"] <- "male"

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
Joelle_data <- Joelle_data[!is.na(Joelle_data$Aspiculuris.Syphacia),]
Joelle_data$presenceAspiSypha[Joelle_data$Aspiculuris.Syphacia > 0] <- 1
Joelle_data$presenceAspiSypha[Joelle_data$Aspiculuris.Syphacia == 0] <- 0

ggplot(Joelle_data, aes(x = HI, y = presenceAspiSypha)) +
  geom_point() +
  geom_smooth() +
  theme_classic()

## Check the prevalence along HI
source("Models/fitBinomial.R")

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
fit <- analyse(Joelle_data,"presenceAspiSypha")

# plot the different hypothesis
plot <- bananaPlots(mod = fit$H0, data = Joelle_data, response = "presenceAspiSypha") 
plot + coord_cartesian(ylim = c(0.45, 0.9)) # zoom in

HeitlingerFieldData <- HeitlingerFieldData[!is.na(HeitlingerFieldData$Aspiculuris_Syphacia),]
HeitlingerFieldData$presenceAspiSypha[HeitlingerFieldData$Aspiculuris_Syphacia > 0] <- 1
HeitlingerFieldData$presenceAspiSypha[HeitlingerFieldData$Aspiculuris_Syphacia == 0] <- 0

ggplot(HeitlingerFieldData, aes(x = HI, y = presenceAspiSypha)) +
  geom_point() +
  geom_smooth() +
  theme_classic()

fit <- analyse(HeitlingerFieldData, "presenceAspiSypha")
plot <- bananaPlots(mod = fit$H0, data = HeitlingerFieldData, response = "presenceAspiSypha") 
plot# + coord_cartesian(ylim = c(0.45, 0.9)) # zoom in

## And ectoparasites?
EctoData <- HeitlingerFieldData[!is.na(HeitlingerFieldData$Ectoparasites),]
EctoData$presenceEcto[EctoData$Ectoparasites == TRUE] <- 1
EctoData$presenceEcto[EctoData$Ectoparasites == FALSE] <- 0

ggplot(EctoData, aes(x = HI, y = presenceEcto)) +
  geom_point() +
  geom_smooth() +
  theme_classic()

fitEcto <- analyse(EctoData, "presenceEcto")
plot <- bananaPlots(mod = fitEcto$H0, data = EctoData, response = "presenceEcto") 
plot# + coord_cartesian(ylim = c(0.45, 0.9)) # zoom in

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
bananaPlots(mod = AspJo$H1, data = Joelle_data, response = "Aspiculuris.Syphacia", islog10 = TRUE) 
# alpha = 1.39 in WATWM paper
coef(AspJo$H1)[["alpha"]]

# Try only infected!!

# pinworms (A. tetraptera and S. obvelata (Joelle_data$Aspiculuris.Syphacia))
AspJo2 <- analyse(Joelle_data[Joelle_data$Aspiculuris.Syphacia > 0, ], "Aspiculuris.Syphacia", 
                 giveParamBounds(Joelle_data[Joelle_data$Aspiculuris.Syphacia > 0, ], "Aspiculuris.Syphacia"))
bananaPlots(mod = AspJo2$H1, data = Joelle_data[Joelle_data$Aspiculuris.Syphacia > 0, ],
            response = "Aspiculuris.Syphacia", islog10 = TRUE) 
# alpha = 1.216812 now no big difference
coef(AspJo2$H1)[["alpha"]]

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

## And ours?
AspJe <- analyse(HeitlingerFieldData, "Aspiculuris_Syphacia",
                  giveParamBounds(HeitlingerFieldData, "Aspiculuris_Syphacia"))
bananaPlots(mod = AspJe$H1, data = HeitlingerFieldData, 
            response = "Aspiculuris_Syphacia",
            islog10 = T)

# And only positive
AspJePOS <- analyse(HeitlingerFieldData[HeitlingerFieldData$Aspiculuris_Syphacia >0,],
                    "Aspiculuris_Syphacia",
                    giveParamBounds(
                      HeitlingerFieldData[HeitlingerFieldData$Aspiculuris_Syphacia >0,],
                      "Aspiculuris_Syphacia"))
bananaPlots(mod = AspJePOS$H3, data = HeitlingerFieldData[HeitlingerFieldData$Aspiculuris_Syphacia >0,],
            response = "Aspiculuris_Syphacia", islog10 = T)


## qPCR vs worms (is there a correlation?)
ggplot(HeitlingerFieldData, 
       aes(delta_ct_MminusE, Aspiculuris_Syphacia + 1)) +
  geom_point(aes(fill = HI), pch = 21, size = 5) + 
  scale_fill_gradient(low = "blue", high = "red") +
  scale_y_log10() +
  geom_smooth(method = "lm") +
  theme_classic()

summary(lm(Aspiculuris_Syphacia ~ delta_ct_MminusE, data = HeitlingerFieldData))
