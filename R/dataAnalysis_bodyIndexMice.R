## source the functions defining meanload and aggregation for the negative binomial
source("Models/fitStudent.R")
source("MLE_hybrid_functions.R")

## Import data
HeitlingerFieldData <- read.csv("../../../Data_important/FinalFullDF_flotationPcrqPCR.csv")
miceTable <- HeitlingerFieldData[!is.na(HeitlingerFieldData$Body_weight) &
                                   !is.na(HeitlingerFieldData$Body_length) &
                                   !is.na(HeitlingerFieldData$HI) &
                                   !is.na(HeitlingerFieldData$Sex) &
                                   (
                                     !is.na(HeitlingerFieldData$OPG) |
                                       !is.na(HeitlingerFieldData$delta_ct_MminusE) |
                                       !is.na(HeitlingerFieldData$PCRstatus)
                                   ), ]

# Works if OPG are integers
# miceTable$OPG <- round(miceTable$OPG)

getDF <- function(rawDF){
  data4stats <- rawDF[names(rawDF) %in% 
                        c("Body_weight", "Body_length", "HI", "OPG", "delta_ct_MminusE", "PCRstatus", "Sex", "Status")]
  data4stats$EimeriaDetected <- NA
  # data4stats$EimeriaDetected[data4stats$OPG == 0] <- "negative"
  # data4stats$EimeriaDetected[data4stats$OPG > 0] <- "positive"
  data4stats$EimeriaDetected[data4stats$delta_ct_MminusE <= -6] <- "negative"
  data4stats$EimeriaDetected[data4stats$delta_ct_MminusE > -6] <- "positive"
  data4stats <- data4stats[!is.na(data4stats$EimeriaDetected),]
  # Use index as in paper https://www.researchgate.net/publication/259551749_Which_body_condition_index_is_best
  # log body mass/log body length
  data4stats$BCI <- log(data4stats$Body_weight, base = 10) / log(data4stats$Body_length, base = 10)
  return(data4stats)
}

data4stats <- getDF(miceTable)

# Remove pregnant/post partum and juveniles
data4stats$status[data4stats$Sex %in% c("F")] <- "non pregnant/lactating female"
data4stats$status[data4stats$Status %in% c("post partum", "post partum (lactating)", "pregnant")] <- "pregnant/lactating female"
data4stats$status[data4stats$Sex %in% c("M")] <- "male"

# Test  our detection of pregnancy in females using BCI
ggplot(data4stats, 
       aes(x = HI, y = BCI, fill = status, group = status)) +
  geom_point(pch = 21, size = 3, alpha = .5)+
  geom_smooth(aes(col = status)) +
  theme_bw()

# Remove pregnant females
d <- data4stats[!data4stats$status %in% c("pregnant/lactating female"), ]

# Regression of BM/BS for males and females (all together, then separate subsp.) 
# Advantage: independant of size!!

# Step 1: fit the model
fit <- lm(Body_weight ~ Body_length * Sex, data = d)

# Step 2: obtain predicted and residual values
d$predicted <- predict(fit)   # Save the predicted values
d$residuals <- residuals(fit) # Save the residual values

# Step 3: plot the actual and predicted values
ggplot(d, aes(x = Body_length, y = Body_weight)) +
  geom_smooth(method = "lm", se = FALSE, color = "lightgrey") +  # Plot regression slope
  geom_segment(aes(xend = Body_length, yend = predicted), alpha = .2) +  # alpha to fade lines
  geom_point(aes(col = delta_ct_MminusE), size = 3) +
  scale_color_gradient(low = "lightgrey", high = "red") +
  geom_point(aes(y = predicted), shape = 1) +
  facet_grid(~ Sex, scales = "free_x") +  # Split panels here by `iv`
  theme_bw()  # Add theme for cleaner look

# Step 4: use residuals as indice
hist(d$residuals[d$Sex =="F"], breaks = 100) # remove outliers, keep [-5,5] interval

d <- d[d$residuals <= 5,]
d$resBMBL <- d$residuals

# give positive values only
d$resBMBL <- d$resBMBL + 5

# Which distribution to choose?
library(MASS)
dat <- d$resBMBL
# let's compute some fits...
fits <- list(
  normal = fitdistr(dat,"normal"),
  # logistic = fitdistr(dat,"logistic"), #tested
  # cauchy = fitdistr(dat,"cauchy"), #tested
  # weibull = fitdistr(dat, "weibull"), #tested (needs positive values)
  student = fitdistr(dat, "t", start = list(m = mean(dat), s = sd(dat), df = 3), lower=c(-1, 0.001,1))
)
# get the logliks for each model...
sapply(fits, function(i) i$loglik)
# STUDENT is the way to go!
ggplot(d, aes(resBMBL)) +
  geom_histogram(aes(y=..density..), bins = 100) + 
  stat_function(fun = dnorm, n = 1e3, args = list(mean = fits$normal$estimate[1], sd = fits$normal$estimate[2]),
                aes(color = "normal"), size = 2) +
  stat_function(fun = dt, n = 1e3, args = list(ncp = fits$student$estimate[1], 
                                               df = fits$student$estimate[3]),
                aes(color = "student"), size = 2) +
  stat_function(fun = dt, n = 1e3, args = list(ncp = fits$student$estimate[1], 
                                               df = 10),
                aes(color = "student10"), size = 2) +
  stat_function(fun = dt, n = 1e3, args = list(ncp = fits$student$estimate[1], 
                                               df = 400),
                aes(color = "student500"), size = 2) +
  stat_function(fun = dt, n = 1e3, args = list(ncp = fits$student$estimate[1], 
                                               df = 1000),
                aes(color = "student1000"), size = 2) +
  stat_function(fun = dt, n = 1e3, args = list(ncp = fits$student$estimate[1], 
                                               df = 5),
                aes(color = "student5"), size = 2) +
  theme_bw(base_size = 24) 

#### Our model
marshallData <- function (data, response) {
  dataForResponse <- data[complete.cases(data[[response]]),]
  dataForResponse_P <- dataForResponse[dataForResponse$EimeriaDetected == "positive",]
  dataForResponse_N <- dataForResponse[dataForResponse$EimeriaDetected == "negative",]
  # dataForResponse_P$EimeriaDetected <- droplevels(dataForResponse_P$EimeriaDetected)
  # dataForResponse_N$EimeriaDetected <- droplevels(dataForResponse_N$EimeriaDetected)
  return(list(
    all = dataForResponse,
    positive = dataForResponse_P,
    negative = dataForResponse_N 
  ))
}

runAll <- function (data, response, mydf) {
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
                   mydfStart = 4, mydfLB = 1, mydfUB = 400)
  marshalledData <- marshallData(data, response)
  print("Fitting for all")
  FitAll <- run(
    data = marshalledData$all,
    response = response,
    hybridIndex = HI, 
    paramBounds = paramBounds, 
    config = defaultConfig
  )
  print("Fitting for positive")
  FitPositive <- run(
    data = marshalledData$positive,
    response = response,
    hybridIndex = HI, 
    paramBounds = paramBounds,
    config = defaultConfig
  )
  print("Fitting for negative")
  FitNegative <- run(
    data = marshalledData$negative,
    response = response,
    hybridIndex = HI, 
    paramBounds = paramBounds,
    config = defaultConfig
  )
  return(list(FitAll = FitAll, FitPositive = FitPositive, FitNegative = FitNegative))
}

analyse <- function(data, response, mydf) {
  print(paste0("Analysing data for response: ", response))
  FitForResponse <- runAll(data, response, mydf)
  
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
  print("Testing H2 positive no alpha vs alpha")
  Gtest(model0 = FitForResponse$FitPositive$fitBasicNoAlpha, 
        model1 = FitForResponse$FitPositive$fitBasicAlpha)
  
  print("Testing H2 negative no alpha vs alpha")
  Gtest(model0 = FitForResponse$FitNegative$fitBasicNoAlpha, 
        model1 = FitForResponse$FitNegative$fitBasicAlpha)
  
  H2 <- list(positive = FitForResponse$FitPositive$fitBasicAlpha,
             negative = FitForResponse$FitNegative$fitBasicAlpha)
  
  # H3: the mean load can differ both across subspecies and between sexes
  print("Testing H3 poitive no alpha vs alpha")
  Gtest(model0 = FitForResponse$FitPositive$fitAdvancedNoAlpha, 
        model1 = FitForResponse$FitPositive$fitAdvancedAlpha)
  
  print("Testing H3 negative no alpha vs alpha")
  Gtest(model0 = FitForResponse$FitNegative$fitAdvancedNoAlpha, 
        model1 = FitForResponse$FitNegative$fitAdvancedAlpha)
  
  H3 <- list(positive = FitForResponse$FitPositive$fitAdvancedAlpha,
             negative = FitForResponse$FitNegative$fitAdvancedAlpha)
  
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

# choose dataset
# data <- d[d$Sex == "F",]
# data <- d[d$Sex == "M",]
data <- d

fit <- analyse(data, "resBMBL")

# plot all
plotAll(mod = fit$H1, data = data, response = "resBMBL", CI = F, 
        labelfory = "resBMBL", isLog10 = F)

# plot 2 groups
modP = fit$H3$positive
modN = fit$H3$negative
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
  geom_point(data = data, aes_string(x = "HI", y = "resBMBL", col = "EimeriaDetected"), size = 3) +
  scale_color_manual(values = c("grey", "red")) +
  geom_line(aes(x = DF$HI, y = DF$loadMLEN), col = "grey32", size = 2) +
  geom_line(aes(x = DF$HI, y = DF$loadMLEP), col = "red", size = 2) +
  theme_bw(base_size = 20)+
  ylab(label = "resBMBL") +
  annotate("text", x = 0.5, y = 4, col = "red", cex = 7,
           label = as.character(round(fit$H3$positive@coef[["alpha"]], 2))) +
  annotate("text", x = 0.5, y = 6, col = "grey32", cex = 7,
           label = as.character(round(fit$H3$negative@coef[["alpha"]], 2)))
