library(ggplot2)
library(reshape2)
library(MASS)
source("Gtest.R")
source("plotBananas.R")

## Import data
HeitlingerFieldData <- read.csv("../../Data_important/FinalFullDF_flotationPcrqPCR.csv")
miceTable <- HeitlingerFieldData[!is.na(HeitlingerFieldData$HI) &
                                   !is.na(HeitlingerFieldData$Sex) &
                                   !is.na(HeitlingerFieldData$delta_ct_MminusE), ]

data4stats <- miceTable[names(miceTable) %in% 
                          c("HI", "OPG", "delta_ct_MminusE", "PCRstatus", "Sex", "Status", "Year")]

# To pass positive I add 6 to all
data4stats$delta_ct_MminusE <- data4stats$delta_ct_MminusE + 6

# changes to do for avoiding plotting error
levels(data4stats$Sex) <- c(levels(data4stats$Sex), "female", "male")
data4stats$Sex[data4stats$Sex == "F"] <- "female"
data4stats$Sex[data4stats$Sex == "M"] <- "male"

# First look
ggplot(data4stats, 
       aes(x = HI, y = delta_ct_MminusE, fill = Sex, group = Sex)) +
  geom_point(pch = 21, size = 3, alpha = .5)+
  geom_smooth(aes(col = Sex)) +
  theme_bw()

# prepare our data for fitting by sex
marshallData <- function (data, response) {
  dataForResponse <- data[complete.cases(data[[response]]),]
  dataForResponse_F <- data[data$Sex == "female",]
  dataForResponse_F <- dataForResponse_F[complete.cases(dataForResponse_F[[response]]),]
  dataForResponse_M <- data[data$Sex == "male",]
  dataForResponse_M <- dataForResponse_M[complete.cases(dataForResponse_M[[response]]),]
  return(list(
    all = dataForResponse,
    female = dataForResponse_F,
    male = dataForResponse_M 
  ))
}

## First question: can we consider the zeros out? If our hypothesis is : all get infected the same way, 
# but hybrids got lower loads?

data4stats$qPCRstatus <- 0
data4stats$qPCRstatus[data4stats$delta_ct_MminusE > 0] <- 1

ggplot(data4stats, aes(x = HI, y = qPCRstatus)) +
  geom_point() +
  geom_smooth() +
  theme_classic()

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
fit <- analyse(data4stats, "qPCRstatus")

# plot the different hypothesis
bananaPlots(mod = fit$H0, data = data4stats, response = "qPCRstatus")

######## Choose a correct distribution for our data ########
data4stats <- data4stats[data4stats$delta_ct_MminusE > 0, ]

dat <- data4stats$delta_ct_MminusE

# let's compute some fits...
fits <- list(
  normal = fitdistr(dat,"normal"),
  student = fitdistr(dat, "t", start = list(m = mean(dat), s = sd(dat), df=3), lower=c(-1, 0.001,1))
)

# get the logliks for each model...
sapply(fits, function(i) i$loglik)
# STUDENT is the way to go!

source("Models/fitStudent.R")

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
                   alphaStart = 0, alphaLB = -5, alphaUB = 5,
                   mydfStart = 1, mydfLB = 1, mydfUB = 10)
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
fit <- analyse(data4stats, "delta_ct_MminusE")

fit$H0

# plot the different hypothesis
bananaPlots(mod = fit$H0, data = data4stats, response = "delta_ct_MminusE")
bananaPlots(mod = fit$H1, data = data4stats, response = "delta_ct_MminusE")
bananaPlots(mod = fit$H2, data = data4stats, response = "delta_ct_MminusE")
bananaPlots(mod = fit$H3, data = data4stats, response = "delta_ct_MminusE")

# check by year
fit2016 <- analyse(data4stats[data4stats$Year %in% 2016, ], "delta_ct_MminusE")
bananaPlots(mod = fit2016$H0, data = data4stats[data4stats$Year %in% 2016, ], response = "delta_ct_MminusE")
# not significant

fit2017 <- analyse(data4stats[data4stats$Year %in% 2017, ], "delta_ct_MminusE")
bananaPlots(mod = fit2017$H0, data = data4stats[data4stats$Year %in% 2017, ], response = "delta_ct_MminusE")
# significant
