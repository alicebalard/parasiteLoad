source("MLE_hybrid_functions.R")
library(ggplot2)
library(reshape2)
library(MASS)

## Import data
HeitlingerFieldData <- read.csv("../../Data_important/FinalFullDF_flotationPcrqPCR.csv")
miceTable <- HeitlingerFieldData[!is.na(HeitlingerFieldData$HI) &
                                   !is.na(HeitlingerFieldData$Sex) &
                                   !is.na(HeitlingerFieldData$delta_ct_MminusE), ]

data4stats <- miceTable[names(miceTable) %in% 
                          c("HI", "OPG", "delta_ct_MminusE", "PCRstatus", "Sex", "Status")]

# To pass positive I add 6 to all
data4stats$delta_ct_MminusE <- data4stats$delta_ct_MminusE + 6

# First look
ggplot(data4stats, 
       aes(x = HI, y = delta_ct_MminusE, fill = Sex, group = Sex)) +
  geom_point(pch = 21, size = 3, alpha = .5)+
  geom_smooth(aes(col = Sex)) +
  theme_bw()

######## Choose a correct distribution for our data ########
data4stats <- data4stats[data4stats$delta_ct_MminusE > 0, ]
dat <- data4stats$delta_ct_MminusE

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
ggplot(data4stats, aes(delta_ct_MminusE)) +
  geom_histogram(aes(y=..density..), bins = 100) + 
  stat_function(fun = dnorm, n = 1e3, args = list(mean = fits$normal$estimate[1], sd = fits$normal$estimate[2]),
                aes(color = "normal"), size = 2) +
  stat_function(fun = dweibull, n = 1e3, args = list(shape = fits$weibull$estimate[1],
                                                     scale = fits$weibull$estimate[2]),
                aes(color = "weibull"), size = 2) +
  stat_function(fun = dweibull, n = 1e3, args = list(shape = fits$weibull$estimate[1],
                                                     scale = 2),
                aes(color = "weibull2"), size = 2) +
  stat_function(fun = dweibull, n = 1e3, args = list(shape = fits$weibull$estimate[1],
                                                     scale = 5),
                aes(color = "weibull5"), size = 2) +
  stat_function(fun = dweibull, n = 1e3, args = list(shape = fits$weibull$estimate[1],
                                                     scale = 10),
                aes(color = "weibull10"), size = 2) +
  stat_function(fun = dweibull, n = 1e3, args = list(shape = fits$weibull$estimate[1],
                                                     scale = 1),
                aes(color = "weibull1"), size = 2) +
  theme_bw(base_size = 24) 

#### Our model
source("Models/fitWeibull.R")

marshallData <- function (data, response) {
  dataForResponse <- data[complete.cases(data[[response]]),]
  dataForResponse_F <- data[data$Sex == "F",]
  dataForResponse_F <- dataForResponse_F[complete.cases(dataForResponse_F[[response]]),]
  dataForResponse_M <- data[data$Sex == "M",]
  dataForResponse_M <- dataForResponse_M[complete.cases(dataForResponse_M[[response]]),]
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
  paramBounds <- c(L1start = mean(na.omit(data[[response]])), 
                   L1LB = min(na.omit(data[[response]])), 
                   L1UB = max(na.omit(data[[response]])), 
                   L2start = mean(na.omit(data[[response]])), 
                   L2LB = min(na.omit(data[[response]])), 
                   L2UB = max(na.omit(data[[response]])), 
                   alphaStart = 0, alphaLB = -5, alphaUB = 5,
                   myshapestart = 1, myshapeLB = 1, myshapeUB = 10)
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

# plot all
plotAll(mod = fit$H1, data = data4stats, response = "delta_ct_MminusE", CI = F, 
        labelfory = "delta_ct_MminusE", isLog10 = F)

# plot 2 groups
modP = fit$H3$female
modN = fit$H3$male
DF <- data.frame(HI = seq(0,1,0.01), 
                 loadMLEP = MeanLoad(L1 = coef(modP)[names(coef(modP)) == "L1"], 
                                     L2 =  coef(modP)[names(coef(modP)) == "L2"], 
                                     alpha =  coef(modP)[names(coef(modP)) == "alpha"],  
                                     hybridIndex = seq(0,1,0.01)), 
                 loadMLEN = MeanLoad(L1 = coef(modN)[names(coef(modN)) == "L1"], 
                                     L2 =  coef(modN)[names(coef(modN)) == "L2"], 
                                     alpha =  coef(modN)[names(coef(modN)) == "alpha"],  
                                     hybridIndex = seq(0,1,0.01))) 

ggplot() +
  geom_point(data = data4stats, aes_string(x = "HI", y = "delta_ct_MminusE", col = "Sex"), size = 3) +
  scale_color_manual(values = c("grey", "red")) +
  geom_line(aes(x = DF$HI, y = DF$loadMLEN), col = "grey32", size = 2) +
  geom_line(aes(x = DF$HI, y = DF$loadMLEP), col = "red", size = 2) +
  theme_bw(base_size = 20)+
  ylab(label = "delta_ct_MminusE") +
  annotate("text", x = 0.5, y = 4, col = "red", cex = 7,
           label = as.character(round(fit$H3$female@coef[["alpha"]], 2))) +
  annotate("text", x = 0.5, y = 6, col = "grey32", cex = 7,
           label = as.character(round(fit$H3$male@coef[["alpha"]], 2)))

