## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----install-------------------------------------------------------------
library(parasiteLoad)
library(bbmle)
require(optimx) # for bbmle it needs to be required(?)
library(ggplot2)
library(MASS)

## ----prepData------------------------------------------------------------
# add 1 for worms/oocysts count!!
WATWMdata$`Aspiculuris.Syphacia+1` <- WATWMdata$Aspiculuris.Syphacia + 1
BALdata$`Aspiculuris.Syphacia+1` <- BALdata$Aspiculuris_Syphacia + 1
BALdata$`OPG+1` <- BALdata$OPG + 1

# and select variables for each analysis
qpcr_data <- BALdata[!is.na(BALdata$HI) &
                       !is.na(BALdata$Sex) &
                       !is.na(BALdata$delta_ct_MminusEAllPos), ]

flotation_data <- BALdata[!is.na(BALdata$`OPG+1`) &
                            !is.na(BALdata$HI) &
                            !is.na(BALdata$Sex), ]


pinworms_data <- BALdata[!is.na(BALdata$`Aspiculuris.Syphacia+1`) &
                           !is.na(BALdata$HI) &
                           !is.na(BALdata$Sex),]

body_data <- BALdata[!is.na(BALdata$Body_weight) &
                       !is.na(BALdata$Body_length) &
                       !is.na(BALdata$HI) &
                       !is.na(BALdata$Sex) &
                       !is.na(BALdata$qPCRstatus) ,]

WATWMdata <- WATWMdata[!is.na(WATWMdata$`Aspiculuris.Syphacia+1`) &
                         !is.na(WATWMdata$HI) &
                         !is.na(WATWMdata$Sex), ]

# all worms together for comparison
d1 <- WATWMdata[c("Sex", "Aspiculuris.Syphacia+1", "HI")]
d2 <- pinworms_data[c("Sex", "Aspiculuris.Syphacia+1", "HI")]
d1$batch <- "WATWM"
d2$batch <- "JENNY"
allWorms <- rbind(d1,d2) 
allWorms$batch <- as.factor(allWorms$batch)

## ---- choosefitsqpcr-----------------------------------------------------

qpcr_data[qpcr_data$delta_ct_MminusEAllPos > 0,"delta_ct_MminusEAllPos"]

fits_test_qpcr <- list(
  normal = MASS::fitdistr(qpcr_data[qpcr_data$delta_ct_MminusEAllPos > 0,"delta_ct_MminusEAllPos"],
                          "normal"),
  student = MASS::fitdistr(qpcr_data[qpcr_data$delta_ct_MminusEAllPos > 0, "delta_ct_MminusEAllPos"], "t",
                     start = list(m = mean(qpcr_data[qpcr_data$delta_ct_MminusEAllPos > 0, "delta_ct_MminusEAllPos"]),
                                  s = sd(qpcr_data[qpcr_data$delta_ct_MminusEAllPos > 0, "delta_ct_MminusEAllPos"]),
                                  df=3), lower=c(-1, 0.001,1))
)

#Get the logliks for each model...

sapply(fits_test_qpcr, function(i) i$loglik)

## ----runfitqpcr----------------------------------------------------------
fit_qpcr_student <- analyse(qpcr_data[qpcr_data$delta_ct_MminusEAllPos > 0, ],
                               response = "delta_ct_MminusEAllPos", 
                            model = "student", group = "Sex")
fit_qpcr_student

