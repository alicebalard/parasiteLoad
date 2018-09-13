source("MLE_hybrid_functions.R")
library(ggplot2)
library(reshape2)
library(MASS)

## Import data
HeitlingerFieldData <- read.csv("../data/FinalFullDF_flotationPcrqPCR.csv")
miceTable <- HeitlingerFieldData[!is.na(HeitlingerFieldData$qPCRstatus) &
                                   !is.na(HeitlingerFieldData$HI) &
                                   !is.na(HeitlingerFieldData$Sex), ]

# To pass positive I add 6 to all
miceTable$delta_ct_MminusE <- miceTable$delta_ct_MminusE + 6

######## Choose a correct distribution for our data ########
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

# Normal parameters: mean and sd
fits$normal$estimate["mean"]
# Weibull parameters: shape and scale
# mean = scale * gamma(1 + 1/shape)
fits$weibull$estimate["scale"] * gamma(1 + 1/fits$weibull$estimate["shape"])

# Observed mean:
mean(miceTable$delta_ct_MminusE[miceTable$delta_ct_MminusE > 0])

######## Fit a normal then a weibull distribution to data
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

## 1. Normal distribution
source("Models/fitNormal.R")

runAll <- function (data, response, paramBounds, mysd) {
  data$response <- data[[response]] # little trick
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
    config = defaultConfig,
    mysd = mysd
  )
  return(list(FitAll = FitAll))
}

analyse <- function(data, response, paramBounds, mysd) {
  print(paste0("Analysing data for response: ", response))
  FitForResponse <- runAll(data, response, paramBounds, mysd)
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
  ####### Compare the hypotheses with G-tests 
  # H1 vs H0
  print("Testing H1 vs H0")
  Gtest(model0 = H0, model1 = H1)
  return(list(H0 = H0, H1 = H1))
}

# remove zeros
miceTable <- miceTable[miceTable$delta_ct_MminusE > 0,]

# NB not enough values for males vs females
paramBounds = c(L1start = mean(miceTable$delta_ct_MminusE), 
                L1LB = min(miceTable$delta_ct_MminusE), 
                L1UB = max(miceTable$delta_ct_MminusE), 
                L2start = mean(miceTable$delta_ct_MminusE), 
                L2LB = min(miceTable$delta_ct_MminusE), 
                L2UB = max(miceTable$delta_ct_MminusE),  
                alphaStart = 0, alphaLB = -5, alphaUB = 5)
mysd = fits$normal$estimate[2]

fitNormal <- analyse(miceTable, "delta_ct_MminusE", paramBounds = paramBounds, mysd)

## 2. Weibull distribution
source("Models/fitWeibull.R")

runAll <- function (data, response, paramBounds, shape) {
  data$response <- data[[response]] # little trick
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
    config = defaultConfig,
    shape = shape
  )
  return(list(FitAll = FitAll))
}

analyse <- function(data, response, paramBounds, shape) {
  print(paste0("Analysing data for response: ", response))
  FitForResponse <- runAll(data, response, paramBounds, mysd)
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
  ####### Compare the hypotheses with G-tests 
  # H1 vs H0
  print("Testing H1 vs H0")
  Gtest(model0 = H0, model1 = H1)
  return(list(H0 = H0, H1 = H1))
}

# remove zeros
miceTable <- miceTable[miceTable$delta_ct_MminusE > 0,]

# NB not enough values for males vs females
paramBounds = c(L1start = mean(miceTable$delta_ct_MminusE), 
                L1LB = min(miceTable$delta_ct_MminusE), 
                L1UB = max(miceTable$delta_ct_MminusE), 
                L2start = mean(miceTable$delta_ct_MminusE), 
                L2LB = min(miceTable$delta_ct_MminusE), 
                L2UB = max(miceTable$delta_ct_MminusE),  
                alphaStart = 0, alphaLB = -5, alphaUB = 5)
shape = fits$weibull$estimate[1]

fitWeibull <- analyse(miceTable, "delta_ct_MminusE", paramBounds = paramBounds, shape)

fitNormal
fitWeibull

## 3. Plot

plotAll(mod = fitNormal$H1, data = miceTable, response = "delta_ct_MminusE", CI = TRUE) +
  ylab(label = "delta CT mouse vs eimeria") +
  annotate("text", x = 0.5, y = 0.58, col = "black", cex = 7,
           label = as.character(round(fitNormal$H1@coef[["alpha"]], 2)))

plotAll(mod = fitWeibull$H1, data = miceTable, response = "delta_ct_MminusE", CI = TRUE) +
  ylab(label = "delta CT mouse vs eimeria") +
  annotate("text", x = 0.5, y = 0.58, col = "black", cex = 7,
           label = as.character(round(fitWeibull$H1@coef[["alpha"]], 2)))

# return(list(fit, plot))
