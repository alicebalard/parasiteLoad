source("Models/BCI_qPCR-NormalDistrib.R")
source("MLE_hybrid_functions.R")
library(ggplot2)
library(reshape2)

## Import data
HeitlingerFieldData <- read.csv("../../Data_important/FinalFullDF_flotationPcrqPCR.csv")
miceTable <- HeitlingerFieldData[!is.na(HeitlingerFieldData$qPCRstatus) &
                                        !is.na(HeitlingerFieldData$HI) &
                                        !is.na(HeitlingerFieldData$Sex), ]

# To pass positive I add 6 to all
miceTable$delta_ct_MminusE <- miceTable$delta_ct_MminusE + 6

## Separate in all, male, female the data frames
qplot(miceTable$delta_ct_MminusE[miceTable$delta_ct_MminusE > 0]) + theme_bw()

# Which distribution to choose?
library(MASS)

dat <- miceTable$delta_ct_MminusE[miceTable$delta_ct_MminusE > 0]

# let's compute some fits...
fits <- list(
  normal = fitdistr(dat,"normal"),
  logistic = fitdistr(dat,"logistic"),
  cauchy = fitdistr(dat,"cauchy"),
  weibull = fitdistr(dat, "weibull"),
  student = fitdistr(dat, "t", start = list(m = mean(dat), s = sd(dat), df=3), lower=c(-1, 0.001,1))
)

# get the logliks for each model...
sapply(fits, function(i) i$loglik)
# WEIBULL is the way to go!

ggplot(miceTable[miceTable$delta_ct_MminusE > 0, ], aes(miceTable$delta_ct_MminusE[miceTable$delta_ct_MminusE > 0])) +
  geom_histogram(aes(y=..density..), bins = 100) +
  stat_function(fun = dnorm, n = 1e3, args = list(mean = fits$normal$estimate[1], sd = fits$normal$estimate[2]),
                aes(color = "normal"), size = 2) +
  stat_function(fun = dweibull, n = 1e3, args = list(shape = fits$weibull$estimate[1], scale = fits$weibull$estimate[2]),
                aes(color = "weibull"), size = 2) +
  theme_bw(base_size = 24)

#### Our model
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

range(miceTable$delta_ct_MminusE)

runAll <- function (data, response, paramBounds = c(L1start = 0, L1LB = 0, L1UB = 10, 
                                                    L2start = 0, L2LB = 0, L2UB = 0, 
                                                    alphaStart = 0, alphaLB = -5, alphaUB = 5)) {
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
  return(list(FitAll = FitAll))
}

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
  
  # # H2: the mean load across subspecies is the same, but can differ between the sexes
  # print("Testing H2 female no alpha vs alpha")
  # Gtest(model0 = FitForResponse$FitFemale$fitBasicNoAlpha, 
  #       model1 = FitForResponse$FitFemale$fitBasicAlpha)
  # 
  # print("Testing H2 male no alpha vs alpha")
  # Gtest(model0 = FitForResponse$FitMale$fitBasicNoAlpha, 
  #       model1 = FitForResponse$FitMale$fitBasicAlpha)
  # 
  # H2 <- list(female = FitForResponse$FitFemale$fitBasicAlpha,
  #            male = FitForResponse$FitMale$fitBasicAlpha)
  # 
  # # H3: the mean load can differ both across subspecies and between sexes
  # print("Testing H3 female no alpha vs alpha")
  # Gtest(model0 = FitForResponse$FitFemale$fitAdvancedNoAlpha, 
  #       model1 = FitForResponse$FitFemale$fitAdvancedAlpha)
  # 
  # print("Testing H3 male no alpha vs alpha")
  # Gtest(model0 = FitForResponse$FitMale$fitAdvancedNoAlpha, 
  #       model1 = FitForResponse$FitMale$fitAdvancedAlpha)
  # 
  # H3 <- list(female = FitForResponse$FitFemale$fitAdvancedAlpha,
  #            male = FitForResponse$FitMale$fitAdvancedAlpha)
  
  ####### Compare the hypotheses with G-tests 
  # H1 vs H0
  print("Testing H1 vs H0")
  Gtest(model0 = H0, model1 = H1)
  
  # # H2 vs H0
  # print("Testing H2 vs H0")
  # Gtest(model0 = H0, model1 = H2)
  # 
  # # H3 vs H1
  # print("Testing H3 vs H1")
  # Gtest(model0 = H1, model1 = H3)
  # 
  # # H3 vs H2
  # print("Testing H3 vs H2")
  # Gtest(model0 = H2, model1 = H3)
  
  return(list(H0 = H0, H1 = H1))
}

# remove zeros
miceTable <- miceTable[miceTable$delta_ct_MminusE > 0,]

## 2016
miceTable2016 <- miceTable[miceTable$Year == 2016,]

## 2017
miceTable2017 <- miceTable[miceTable$Year == 2017,]

perYear <- function(miceTable){
  # not enough values for males vs females
  paramBounds = c(L1start = mean(miceTable$delta_ct_MminusE), L1LB = min(miceTable$delta_ct_MminusE), L1UB = max(miceTable$delta_ct_MminusE), 
                  L2start = mean(miceTable$delta_ct_MminusE), L2LB = min(miceTable$delta_ct_MminusE), L2UB = max(miceTable$delta_ct_MminusE),  
                  alphaStart = 0, alphaLB = -5, alphaUB = 5)
  
  fit <- analyse(miceTable, "delta_ct_MminusE", paramBounds = paramBounds)
  
  plot <- plotAll(mod = fit$H1, data = miceTable, response = "delta_ct_MminusE", CI = TRUE) +
    ylab(label = "delta CT mouse vs eimeria") +
    annotate("text", x = 0.5, y = 0.58, col = "black", cex = 7,
             label = as.character(round(fit$H1@coef[["alpha"]], 2)))
  
  return(list(fit, plot))
}

perYear(miceTable2017)
