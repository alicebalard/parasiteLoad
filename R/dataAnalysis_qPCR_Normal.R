source("Models/BCI_qPCR-NormalDistrib.R")
source("MLE_hybrid_functions.R")
library(ggplot2)
library(reshape2)

## Import data
qpcrData <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/raw_data/Eimeria_detection/qPCR_2016.csv")
names(qpcrData)[1] <- "Mouse_ID"

# Did Enas calculate the other way around? 
qpcrData$delta_ct_cewe[qpcrData$observer_qpcr == "Enas"] <- 
  - qpcrData$delta_ct_cewe[qpcrData$observer_qpcr == "Enas"]

qpcrData$delta_ct_ilwe[qpcrData$observer_qpcr == "Enas"] <- 
  - qpcrData$delta_ct_ilwe[qpcrData$observer_qpcr == "Enas"]

# deltaCT = ct eimeria - ct mouse. If high infection, low deltaCT
# -deltaCT = ct mouse - ct eimeria
qpcrData$status[qpcrData$delta_ct_cewe > 6 & qpcrData$delta_ct_ilwe > 6] <- "non infected"
qpcrData$status[qpcrData$delta_ct_cewe < 6 & qpcrData$delta_ct_ilwe > 6] <- "infected cecum"
qpcrData$status[qpcrData$delta_ct_cewe > 6 & qpcrData$delta_ct_ilwe < 6] <- "infected ileum"

qpcrData$status[
  qpcrData$delta_ct_cewe < 6 & 
    qpcrData$delta_ct_ilwe < 6 & 
    qpcrData$delta_ct_cewe < qpcrData$delta_ct_ilwe] <- "cecum stronger"
qpcrData$status[
  qpcrData$delta_ct_cewe < 6 & 
    qpcrData$delta_ct_ilwe < 6 & 
    qpcrData$delta_ct_cewe > qpcrData$delta_ct_ilwe] <- "ileum stronger"

# and keep the infected segment value OR the higher value 
qpcrData$delta_ct[
  qpcrData$status %in% c("infected cecum", "cecum stronger")] <- 
  qpcrData$delta_ct_cewe[
    qpcrData$status %in% c("infected cecum", "cecum stronger")] 

qpcrData$delta_ct[
  qpcrData$status %in% c("infected ileum", "ileum stronger")] <- 
  qpcrData$delta_ct_ilwe[
    qpcrData$status %in% c("infected ileum", "ileum stronger")] 

# Turn around
qpcrData$delta_ct_MminusE <- - qpcrData$delta_ct

# Set floor values
qpcrData$delta_ct_MminusE[is.na(qpcrData$delta_ct_MminusE)] <- -6

# Add infos
miceTable <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/raw_data/MiceTable_2014to2017.csv")

miceTable <- merge(miceTable, qpcrData)

# To pass positive I add 6 to all
miceTable$delta_ct_MminusE <- miceTable$delta_ct_MminusE + 6

## Separate in all, male, female the data frames
qplot(miceTable$delta_ct_MminusE) + theme_bw()

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
miceTable <- miceTable[miceTable$status != "non infected",]

# not enough values for males vs females
paramBounds = c(L1start = mean(miceTable$delta_ct_MminusE), L1LB = min(miceTable$delta_ct_MminusE), L1UB = max(miceTable$delta_ct_MminusE), 
                L2start = mean(miceTable$delta_ct_MminusE), L2LB = min(miceTable$delta_ct_MminusE), L2UB = max(miceTable$delta_ct_MminusE),  
                alphaStart = 0, alphaLB = -5, alphaUB = 5)

fit <- analyse(miceTable, "delta_ct_MminusE", paramBounds = paramBounds)

plotAll(mod = fit$H1, data = miceTable, response = "delta_ct_MminusE", CI = TRUE) +
  ylab(label = "delta CT mouse vs eimeria") +
  annotate("text", x = 0.5, y = 0.58, col = "black", cex = 7,
           label = as.character(round(fit$H1@coef[["alpha"]], 2))) 

summary(fit$H1)

