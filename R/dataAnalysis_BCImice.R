## source the functions defining meanload and aggregation for the negative binomial
source("Models/BCI_qPCR-NormalDistrib.R")
source("MLE_hybrid_functions.R")

## Import data
HeitlingerFieldData <- read.csv("../../Data_important/FinalFullDF_flotationPcrqPCR.csv")
miceTable <- HeitlingerFieldData[!is.na(HeitlingerFieldData$Body_weight) &
                                   !is.na(HeitlingerFieldData$Body_length) &
                                   !is.na(HeitlingerFieldData$HI) &
                                   !is.na(HeitlingerFieldData$Sex) &
                                   (
                                     !is.na(HeitlingerFieldData$OPG) |
                                       !is.na(HeitlingerFieldData$delta_ct_MminusE) |
                                       !is.na(HeitlingerFieldData$PCRstatus)
                                   ), ]

# # Works if OPG are integers
# miceTable$OPG <- round(miceTable$OPG)

getDF <- function(rawDF){
  data4stats <- rawDF[names(rawDF) %in% 
                        c("Body_weight", "Body_length", "HI", "OPG", "delta_ct_MminusE", "PCRstatus", "Sex", "Status")]
  
  # data4stats$EimeriaDetected[data4stats$OPG == 0] <- "negative"
  # data4stats$EimeriaDetected[data4stats$OPG > 0] <- "positive"
  data4stats$EimeriaDetected <- NA
  data4stats$EimeriaDetected[data4stats$delta_ct_MminusE <= -6] <- "negative"
  data4stats$EimeriaDetected[data4stats$delta_ct_MminusE > -6] <- "positive"

  data4stats <- data4stats[!is.na(data4stats$EimeriaDetected),]
  
  # Try index as in paper https://www.researchgate.net/publication/259551749_Which_body_condition_index_is_best
  # log body mass/log body length
  data4stats$BCI <- log(data4stats$Body_weight, base = 10) / log(data4stats$Body_length, base = 10)

  # Remove pregnant/post partum and juveniles
  data4stats <- data4stats[!data4stats$Status %in% c("post partum", "post partum (lactating)", "pregnant", "young"), ]
  
  return(data4stats)
}

# Fit distribution Normal vs Student
library(MASS)
# normal fit
fit <- fitdistr(data4stats$BCI, densfun="normal")  # we assume my_data ~ Normal(?,?)fit
hist(data4stats$BCI, pch=20, breaks=50, prob=TRUE, main="")
curve(dnorm(x, fit$estimate[1], fit$estimate[2]), col="red", lwd=2, add=T)
# student fit
fit2 <- fitdistr(data4stats$BCI, "t", start = list(m=mean(my_data),s=sd(my_data), df=3), lower=c(-1, 0.001,1))
mu.std = fit2$estimate[["m"]]
lambda = fit2$estimate[["s"]]
nu = fit2$estimate[["df"]]
curve(dt((x-mu.std)/lambda, nu)/lambda, col="blue", lwd=2, add=TRUE, yaxt="n")


# getDFPregnantVsNotPregnant <- function(rawDF){
#   data4stats <- rawDF[names(rawDF) %in% 
#                         c("Body_weight", "Body_length", "HI", "OPG", "delta_ct_MminusE", "PCRstatus", "Sex", "Status")]
#   
#   # data4stats$EimeriaDetected[data4stats$OPG == 0] <- "negative"
#   # data4stats$EimeriaDetected[data4stats$OPG > 0] <- "positive"
#   data4stats$EimeriaDetected <- NA
#   data4stats$EimeriaDetected[is.na(data4stats$Status) & data4stats$Sex == "F"] <- "negative"
#   data4stats$EimeriaDetected[data4stats$Status %in% c("post partum", "post partum (lactating)", "pregnant") &
#                                       data4stats$Sex == "F"] <- "positive"
#   
#   data4stats <- data4stats[!is.na(data4stats$EimeriaDetected),]
#   
#   # Try index as in paper https://www.researchgate.net/publication/259551749_Which_body_condition_index_is_best
#   # log body mass/log body length
#   data4stats$BCI <- log(data4stats$Body_weight, base = 10) / log(data4stats$Body_length, base = 10)
#   
#   # # Remove pregnant/post partum and juveniles
#   # data4stats <- data4stats[!data4stats$Status %in% c("post partum", "post partum (lactating)", "pregnant", "young"), ]
#   
#   return(data4stats)
# }

ggplot(data4stats, aes(x = HI, y = BCI, fill = Sex, group = Sex)) +
  geom_point(pch = 21, size = 3, alpha = .5)+
  geom_smooth(aes(col = Sex)) +
  theme_bw()

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

runAll <- function (data, response) {
  print(paste0("Fit for the response: ", response))
  defaultConfig <- list(optimizer = "optimx",
                        method = c("bobyqa", "L-BFGS-B"),
                        control = list(follow.on = TRUE))
  paramBounds <- c(L1start = mean(na.omit(data[[response]])), 
                   L1LB = 0, 
                   L1UB = max(na.omit(data[[response]])), 
                   L2start = mean(na.omit(data[[response]])), 
                   L2LB = 0, 
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

fit <- analyse(data4stats, "BCI")

plotAll(mod = fit$H1, data = data4stats, response = "BCI", CI = FALSE)

plot2sexes(modF = fit$H3$positive, modM = fit$H3$negative, data = data4stats, 
           response = "BCI", CI = FALSE, cols = c("grey", "black"), 
           mygroup = "EimeriaDetected", switchlevels = TRUE)

plot2groups <- function(modP, modN, data, response, mygroup = "EimeriaDetected",
                        cols = c("grey", "black")){
  fit <- analyse(data, response)
  data$response <- data[[response]]
  data$log10resp <- log10(data$response + 1)
  ## Draw the line for the parameters at their MLE, alpha varying 
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
   # geom_point(data = data, aes_string(x = "HI", y = "log10resp", color = mygroup)) + 

   geom_point(data = data, aes_string(x = "HI", y = "response", color = mygroup)) +
   scale_color_manual(values = cols) +
   # geom_line(aes(x = DF$HI, y = log10(DF$loadMLEN + 1)), col = "grey32", size = 3) +
   # geom_line(aes(x = DF$HI, y = log10(DF$loadMLEP + 1)), col = "red", size = 3) +
   geom_line(aes(x = DF$HI, y = DF$loadMLEN), col = "grey32", size = 3) +
   geom_line(aes(x = DF$HI, y = DF$loadMLEP), col = "red", size = 3) +
   theme_bw(base_size = 20)+
   ylab(label = "BCI") +
   annotate("text", x = 0.5, y = 0.20, col = "red", cex = 7,
            label = as.character(round(fit$H3$positive@coef[["alpha"]], 2))) +
   annotate("text", x = 0.5, y = 0.2125, col = "grey32", cex = 7,
            label = as.character(round(fit$H3$negative@coef[["alpha"]], 2)))
}

plot2groups(modP = fit$H3$positive, modN = fit$H3$negative, 
           data = data4stats, 
           response = "BCI", cols = c("grey32", "red"))


